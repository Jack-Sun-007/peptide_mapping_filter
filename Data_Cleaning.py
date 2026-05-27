import os
import pandas as pd


# 氨基酸序列数字补全至3位方便排序
def pad_positions(position):
    try:
        start, end = position.split('-')
        return f"{start.zfill(3)}-{end.zfill(3)}"
    except ValueError:
        # 如果无法正常分割成两部分，就返回原始字符串
        return position


# 查找文件夹的csv或excel
def get_csv(filepath):
    try:
        files = os.listdir(filepath)  # 列出目录中的所有文件和文件夹
        files = [f for f in files if os.path.isfile(os.path.join(filepath, f))]  # 过滤出所有文件
        for i, file in enumerate(files, 1):  # 用自然数序号展示
            print(f'{i}. {file}')
        index = int(input("请输入文件的序号: "))
        file = files[index - 1]
        file_path = os.path.abspath(os.path.join(filepath, file))
        return file_path
    except FileNotFoundError:
        print("ERROR")
        return []


# 对csv进行数据清洗
def csv_dataclean(path):
    df = pd.read_csv(path)
    # 指定删除包含列'Identification'空值的行
    df = df.dropna(subset=['Identification'])
    # 优化单个样品的情况:单个样品列数为23，两个样品列数为54
    if len(df.columns) < 30:
        df['Max MS Area'] = df['MS Area'].copy()
        df['MS Area:1'] = df['MS Area'].copy()
        df['∆ ppm:1'] = df['∆ ppm'].copy()
        df['ID Type:1'] = df['ID Type'].copy()
        df['Mono Mass Exp.:1'] = df['Mono Mass Exp.'].copy()
    else:
        pass
    # 过滤mass error ±10 ppm
    df = df.loc[df['∆ ppm'] >= -10]
    df = df.loc[df['∆ ppm'] <= 10]
    # 先按轻重链排序，再按氨基酸序列排序
    df['Positions'] = df['Positions'].apply(pad_positions)
    df = df.sort_values(by=["Protein", "Positions", "Max MS Area"], ascending=[False, True, False])
    # 删掉多余列
    columns_to_drop = list()
    sample_number = 1
    for column in df.columns:
        if 'Ratio' in column:
            sample_number += 1
            columns_to_drop.append(column)
        if '%CV' in column:
            columns_to_drop.append(column)
        if 'Cond.' in column:
            columns_to_drop.append(column)
        if 'Avg' in column:
            columns_to_drop.append(column)
        if 'RT' in column and '.raw' in column:
            columns_to_drop.append(column)
        if 'M/Z' in column and '.raw' in column:
            columns_to_drop.append(column)
    columns_to_drop.extend(['Level', '∆ ppm', 'Conf. Score', 'Integration Type', 'Protein', 'Best ASR'])
    df.drop(columns=columns_to_drop, inplace=True)
    # 重新排序列
    columns = df.columns.tolist()
    msarea, deltappm, idtype, monomass, none = list(), list(), list(), list(), list()
    for i in columns:
        if 'MS Area:' in i:
            msarea.append(i)
        elif '∆ ppm:' in i:
            deltappm.append(i)
        elif 'ID Type:' in i:
            idtype.append(i)
        elif 'Mono Mass Exp.:' in i:
            monomass.append(i)
        else:
            none.append(i)
    monoppm_list = [item for pair in zip(monomass, deltappm) for item in pair]
    new_columns = none + msarea + monoppm_list + idtype
    df = df.reindex(columns=new_columns)
    # 过滤MS Area<1E4的nonspecific
    df = df[~((df['Mod'].str.contains("nonspecific", na=False)) & (df['Max MS Area'] < 10000))]
    # 过滤杂质修饰
    df = df[~df['Mod'].str.contains("Al3\\+|Ca\\+|Ca2\\+|Fe2\\+|Fe3\\+|Na\\+|2x|GasPhase", na=False)]
    # 筛选Missed Cleavages < 4
    df = df.loc[df['Missed Cleavages'] < 4]
    df = df.reset_index(drop=True)
    # 字体字号设置
    styled_df = df.style.set_properties(**{
        'font-family': 'Times New Roman',
        'font-size': '9pt'
    })
    # 导出
    filename = "1.Raw_" + str(path.split("\\")[-1].split(".")[0]) + ".xlsx"
    styled_df.to_excel(filename, engine='openpyxl', index=True)
    # 筛选ID Type=MS2供后续肽段查找用
    df = df.loc[(df['ID Type'] == 'MS2')]
    return df, filename[6:]


# 对肽段模板进行数据合并
def merge_data(Template_path, df, df_filename, trypsin_or_chymo):
    # trypsin_or_chymo: trypsin=1 chymo=0
    a_df = pd.read_excel(Template_path, sheet_name=trypsin_or_chymo)
    b_df = df
    # 依据理论分子量和肽段序列筛选肽段唯一结果，按最大丰度筛选 20251208修改
    # b_unique_df = b_df.drop_duplicates(subset=["Mono Mass Theo."], keep='first')
    b_unique_df = b_df.drop_duplicates(
        subset=["Mono Mass Theo.", "Peptide Sequence"],
        keep='first'
    )
    # 合并数据，依据理论分子量列合并
    merged_df = pd.merge(a_df, b_unique_df, on=["Mono Mass Theo.", "Peptide Sequence"], how="left")
    # 删掉多余列
    columns_to_drop = list()
    for column in merged_df.columns:
        if 'ID Type' in column:
            columns_to_drop.append(column)
    columns_to_drop.extend(['Positions', 'Mono Mass Exp.', 'Max MS Area', 'Comment'])
    merged_df.drop(columns=columns_to_drop, inplace=True)
    # 字体字号设置
    style_df = merged_df.style.set_properties(**{
        'font-family': 'Times New Roman',
        'font-size': '9pt'
    })
    # 将合并后的结果保存到一个新的Excel文件中
    filename = "2.Check_" + df_filename
    style_df.to_excel(filename, engine='openpyxl', index=False)
    # 删掉多余列
    columns_to_drop2 = list()
    for column in merged_df.columns:
        if 'MS Area' in column:
            columns_to_drop2.append(column)
    columns_to_drop2.extend(
        ['No._y', 'Identification_y', 'Mod_y', 'Site', 'M/Z', 'Charge St.', 'Missed Cleavages'])
    merged_df.drop(columns=columns_to_drop2, inplace=True)
    # 将列'RT'移动到'Mono Mass Theo.'前面
    cols = merged_df.columns.tolist()
    rt_index = cols.index('RT')
    mono_mass_theo_index = cols.index('Mono Mass Theo.')
    cols.insert(mono_mass_theo_index, cols.pop(rt_index))
    merged_df = merged_df[cols]
    # 小数点位数格式设置
    merged_df['RT'] = merged_df['RT'].map(lambda x: f'{x:.2f}' if pd.notnull(x) else x)
    merged_df['Mono Mass Theo.'] = merged_df['Mono Mass Theo.'].map(lambda x: f'{x:.4f}')
    exp_columns = [col for col in merged_df.columns if 'Mono Mass Exp.' in col]
    exp_columns2 = [col for col in merged_df.columns if '∆ ppm:' in col]
    merged_df[exp_columns] = merged_df[exp_columns].map(lambda x: f'{x:.4f}' if pd.notnull(x) else x)
    merged_df[exp_columns2] = merged_df[exp_columns2].map(lambda x: f'{x:.2f}' if pd.notnull(x) else x)
    # 字体字号设置
    styled_df = merged_df.style.set_properties(**{
        'font-family': 'Times New Roman',
        'font-size': '9pt'
    })
    # 将合并后的结果保存到一个新的Excel文件中
    filename = "3.Final_" + df_filename
    styled_df.to_excel(filename, engine='openpyxl', index=False)
    return []


def run_enzyme_process(template_path, folder, enzyme_flag):
    data_path = get_csv(folder)
    df, filename = csv_dataclean(data_path)
    merge_data(template_path, df, filename, enzyme_flag)


def main():
    print("选择肽图模板")
    template_path = get_csv(".\\1.Template")
    options = {
        1: ("Trypsin", ".\\2.Raw data with Trypsin", 1),
        2: ("Chymotrypsin", ".\\3.Raw data with Chymotrypsin", 0),
        3: "both"
    }

    try:
        switch = int(input(
            "选择运行方式：\n"
            "1. 只运行 Trypsin\n"
            "2. 只运行 Chymotrypsin\n"
            "3. Trypsin + Chymotrypsin 一起运行\n"
            "请输入数字："
        ))
    except ValueError:
        print("输入无效，请输入数字 1/2/3")
        return {"status": "error", "message": "invalid input"}

    # 单酶运行
    if switch in (1, 2):
        enzyme_name, folder, flag = options[switch]
        run_enzyme_process(template_path, folder, flag)
        return {"status": "success", "run": enzyme_name}

    # 双酶运行
    elif switch == 3:
        print("开始处理：Trypsin + Chymotrypsin")
        run_enzyme_process(template_path, ".\\2.Raw data with Trypsin", 1)
        run_enzyme_process(template_path, ".\\3.Raw data with Chymotrypsin", 0)
        return {"status": "success", "run": "both"}

    else:
        print("输入有误，请输入 1/2/3")
        return {"status": "error", "message": "invalid option"}


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print("出错啦！需重新运行，出错信息如下：")
        print(str(e))
        print("按任意键退出......")
        input()

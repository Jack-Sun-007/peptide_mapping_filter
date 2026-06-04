import os
import re
import pandas as pd


# 1. 配置区（只改这里）
# 筛选ppm范围
PPM_RANGE = (-10, 10)
# 筛选去除的杂质修饰
IMPURITIES = [
    "Al3\\+", "Ca\\+", "Ca2\\+", "Fe2\\+", "Fe3\\+",
    "Na\\+", "2x", "GasPhase", "GasPhaseH2OLoss",
    "Trisulfide", "Thioether"
]
# 常规表去除的列
DROP_KEYWORDS = ["Ratio", "%CV", "Cond.", "Avg"]
DROP_COLUMNS = [
    "Level", "∆ ppm", "Conf. Score",
    "Integration Type", "Protein", "Best ASR"
]
# Final表去除的列
FINAL_DROP = [
    "No.", "Positions", "M/Z", "Charge St.",
    "Missed Cleavages", "Mono Mass Exp.",
    "Max MS Area", "Comment"
]
# Identification列开头的数字变更成肽段
CHAIN_MAP = {"1": "LC", "2": "HC"}


# 2. 工具函数
# 二硫键位点查重，不同肽段中如果出现其中一个的重复位点，只保留MS Area最强的那条
def site_key(site):
    if pd.isna(site):
        return ()
    keys = []
    for s in site.split("/"):
        prefix = int(s.split(":")[0])
        m = re.search(r"C(\d+)", s)
        cnum = int(m.group(1)) if m else float("inf")
        keys.append((prefix, cnum))
    return tuple(sorted(keys))


# 选择原始csv文件
def select_file(folder):
    files = [f for f in os.listdir(folder) if os.path.isfile(os.path.join(folder, f))]
    for i, f in enumerate(files, 1):
        print(f"{i}. {f}")
    return os.path.abspath(os.path.join(folder, files[int(input("请选择文件序号: ")) - 1]))


# 导出带字体大小格式的excel
def style_output(df, filename):
    styled_df = df.style.set_properties(**{
        'font-family': 'Times New Roman',
        'font-size': '9pt'
    })
    styled_df.to_excel(filename, engine='openpyxl', index=False)
    return


# 3. 核心清洗流程
def clean_dataframe(path):
    df = pd.read_csv(path).dropna(subset=["Identification"])
    # 单样品兼容
    if len(df.columns) < 30:
        for new, old in {
            "Max MS Area": "MS Area",
            "MS Area:1": "MS Area",
            "∆ ppm:1": "∆ ppm",
            "ID Type:1": "ID Type",
            "Mono Mass Exp.:1": "Mono Mass Exp."
        }.items():
            df[new] = df[old]
    # ppm 过滤
    df = df[df["∆ ppm"].between(*PPM_RANGE)]
    # 保留时间过滤
    df = df[df["RT"] <= 63]
    # 删除无用列
    drop_cols = [
        c for c in df.columns
        if any(k in c for k in DROP_KEYWORDS)
        or (("RT" in c or "M/Z" in c) and ".raw" in c)
    ] + DROP_COLUMNS
    df = df.drop(columns=drop_cols)
    # 列顺序
    other, ms, ppm, mono, idt = [], [], [], [], []
    for c in df.columns:
        if "MS Area:" in c: ms.append(c)
        elif "∆ ppm:" in c: ppm.append(c)
        elif "Mono Mass Exp.:" in c: mono.append(c)
        elif "ID Type:" in c: idt.append(c)
        else: other.append(c)
    df = df[other + ms + [x for p in zip(mono, ppm) for x in p] + idt]
    # 杂质过滤
    df = df[~df["Mod"].str.contains("|".join(IMPURITIES), na=False)]
    return df


# 4. 二硫键筛选核心逻辑
def build_check_table(df, base):
    df = df[
        df["Mod"].str.contains("1ss|2ss", na=False) &
        df["ID Type"].isin(["MS2", "Both"])
    ]
    df = df.loc[df.groupby("Peptide Sequence")["Max MS Area"].idxmax()]
    # df = df.loc[df.groupby("Site")["Max MS Area"].idxmax()]
    df = df.sort_values("Max MS Area", ascending=False)
    style_output(df, f"2.check_{base}.xlsx")
    seen, keep = set(), []
    for i, r in df.iterrows():
        if pd.isna(r["Site"]):
            keep.append(i)
            continue
        s = set(r["Site"].split("/"))
        # 判断CPPC序列如果相同就排除，特定找一个漏切一个不漏切的肽段
        p = r["Peptide Sequence"].split("/")
        if len(p) == 2:
            if "CPPC" in p[0]:
                if p[0] == p[1]:
                    continue
        if not (s & seen):
            keep.append(i)
            seen |= s
    df = df.loc[keep]
    return (
        df.assign(_k=df["Site"].apply(site_key))
          .sort_values("_k")
          .drop(columns="_k")
          .reset_index(drop=True)
    )


# 5. 转换成能复制到报告中的Final表
def build_final_table(df):
    df = df.drop(columns=[c for c in df.columns if "MS Area" in c or "ID Type" in c] + FINAL_DROP)
    df = df[["Peptide Sequence", "Identification", "Site", "Mod"] +
            [c for c in df.columns if c not in ["Peptide Sequence", "Identification", "Site", "Mod"]]]
    pat = re.compile(r"(^|/)(\d+):")
    def convert(text):
        if pd.isna(text):
            return text
        text = text.split("=", 1)[0]
        return pat.sub(lambda m: f"{m.group(1)}{CHAIN_MAP.get(m.group(2), m.group(2))}:", text)

    for c in ["Identification", "Site"]:
        df[c] = df[c].apply(convert)
    # 小数点位数格式设置
    df['RT'] = df['RT'].map(lambda x: f'{x:.2f}' if pd.notnull(x) else x)
    df['Mono Mass Theo.'] = df['Mono Mass Theo.'].map(lambda x: f'{x:.4f}')
    exp_columns = [col for col in df.columns if 'Mono Mass Exp.' in col]
    exp_columns2 = [col for col in df.columns if '∆ ppm:' in col]
    df[exp_columns] = df[exp_columns].map(lambda x: f'{x:.4f}' if pd.notnull(x) else x)
    df[exp_columns2] = df[exp_columns2].map(lambda x: f'{x:.2f}' if pd.notnull(x) else x)
    return df


# 6.主入口
def main():
    path = select_file(".\\1.Raw data\\")
    base = os.path.splitext(os.path.basename(path))[0]
    df_raw = clean_dataframe(path)
    df_check = build_check_table(df_raw, base)
    df_final = build_final_table(df_check)

    style_output(df_raw, f"1.raw_{base}.xlsx")
    #style_output(df_check, f"2.check_{base}.xlsx")
    style_output(df_final, f"3.Final_{base}.xlsx")


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print("出错啦！需重新运行，出错信息如下：")
        print(str(e))
        print("按任意键退出......")
        input()

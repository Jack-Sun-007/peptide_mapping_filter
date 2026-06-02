import os
import pandas as pd
from collections import defaultdict
import PySimpleGUI as sg


# === 酶切规则 ===
enzyme_rules = {
    "Trypsin": {'type': 'C', 'cut': set('KR'), 'nocut_next': None},
    "Chymotrypsin": {'type': 'C', 'cut': set('FYWL'), 'nocut_next': None},
    "LysC": {'type': 'C', 'cut': set('K'), 'nocut_next': None},
    "AspN": {'type': 'N', 'cut': set('D'), 'nocut_next': None}
}
# 单酶切割
def enzyme_digest(seq, enzyme):
    rule = enzyme_rules[enzyme]
    cut_sites = rule['cut']
    no_cut_next = rule['nocut_next']
    cut_type = rule['type']
    peptides = []
    start = 0
    if cut_type == 'C':
        for i in range(len(seq) - 1):
            aa = seq[i]
            next_aa = seq[i + 1]
            if aa in cut_sites and (no_cut_next is None or next_aa != no_cut_next):
                peptides.append(seq[start:i + 1])
                start = i + 1
        peptides.append(seq[start:])
    elif cut_type == 'N':
        cut_positions = [0]
        for i in range(1, len(seq)):
            if seq[i] in cut_sites:
                cut_positions.append(i)
        cut_positions.append(len(seq))
        for i in range(len(cut_positions) - 1):
            peptides.append(seq[cut_positions[i]:cut_positions[i + 1]])
    return peptides


# 混合酶切：依次应用多个酶
def mixed_digest(seq, enzymes):
    peptides = [seq]
    for enzyme in enzymes:
        new_peptides = []
        for p in peptides:
            new_peptides.extend(enzyme_digest(p, enzyme))
        peptides = new_peptides
    return peptides


# 对所有链进行指定模式酶切
def digest_all(sequences, mode):
    results = []
    all_peptides = defaultdict(list)
    for chain, seq in sequences.items():
        if mode == "Trypsin_AspN":
            peptides = mixed_digest(seq, ["AspN", "Trypsin"])
        elif mode == "LysC_AspN":
            peptides = mixed_digest(seq, ["AspN", "LysC"])
        elif mode == "LysC_Chymotrypsin":
            peptides = mixed_digest(seq, ["LysC", "Chymotrypsin"])
        else:
            peptides = enzyme_digest(seq, mode)
        pos = 1
        for pep in peptides:
            start_aa = seq[pos - 1]
            end_aa = seq[pos + len(pep) - 2]
            # 单氨基酸特殊编号
            if len(pep) == 1:
                label = f"{chain}:{start_aa}{pos}"
            else:
                label = f"{chain}:{start_aa}{pos}-{end_aa}{pos + len(pep) - 1}"
            results.append({
                "Chain": chain,
                "Peptide_ID": label,
                "Sequence": pep,
                "Length": len(pep)
            })
            all_peptides[pep].append(label)
            pos += len(pep)
    # 标注重复
    for r in results:
        seq = r["Sequence"]
        if len(all_peptides[seq]) > 1:
            r["Is_Duplicate"] = "Yes"
            r["Duplicate_Peptide_IDs"] = "; ".join(all_peptides[seq])
        else:
            r["Is_Duplicate"] = ""
            r["Duplicate_Peptide_IDs"] = ""
    return results


def get_input(inputdata):
    result = {}
    keys = list(inputdata.keys())
    i = 0
    while i < len(keys):
        if inputdata[keys[i]] is True:
            # 取后面两项作为 key-value
            k = inputdata.get(keys[i + 1], None)
            v = inputdata.get(keys[i + 2], None)
            if isinstance(v, str):
                # 去除空格、制表符、换行符，并转为大写
                v = v.replace(" ", "").replace("\t", "").replace("\n", "").upper()
            if k and v:  # 确保非空字符串
                result[k] = v
            i += 3  # 跳过这三个
        elif inputdata[keys[i]] is False:
            i += 3  # False 也跳过3个
        else:
            i += 1  # 其他情况正常往后
    return result


# 主程序
def main(sequences):
    output_file = "Digestion_All_Enzymes_Analyze.xlsx"
    enzyme_modes = ["Trypsin", "Chymotrypsin", "AspN", "Trypsin_AspN", "LysC_AspN", "LysC_Chymotrypsin"]
    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        for mode in enzyme_modes:
            df = pd.DataFrame(digest_all(sequences, mode))
            df.to_excel(writer, sheet_name=mode, index=False)
    os.startfile("Digestion_All_Enzymes_Analyze.xlsx")


def GUI():
    sg.theme('GrayGrayGray')
    # 界面布局，将会按照列表顺序从上往下依次排列，二级列表中，从左往右依此排列
    layout = [[sg.Text('         链名称'), sg.Text('                 氨基酸序列')],
              [sg.Checkbox('', default=True), sg.InputText("LC", size=5), sg.Multiline(size=(40, 1))],
              [sg.Checkbox('', default=True), sg.InputText("HC", size=5), sg.Multiline(size=(40, 1))],
              [sg.Checkbox('', default=False), sg.InputText(size=5), sg.Multiline(size=(40, 1))],
              [sg.Checkbox('', default=False), sg.InputText(size=5), sg.Multiline(size=(40, 1))],
              [sg.Checkbox('', default=False), sg.InputText(size=5), sg.Multiline(size=(40, 1))],
              [sg.Text('By SHJ -- Version 1.0'),sg.Button('取消', size=10), sg.Button('运行', size=10)],
              ]
    # 创造窗口
    window = sg.Window('氨基酸理论酶切分析', layout, size=(400, 260))
    # 事件循环并获取输入值
    while True:
        event, values = window.read()
        if event in (None, '取消'):  # 如果用户关闭窗口或点击`Cancel`
            break
        if event in '运行':
            print(values)
            a = get_input(values)
            main(a)
            break
    window.close()


if __name__ == '__main__':
    try:
        GUI()
    except Exception as e:
        sg.popup_error(f'运行出错!\n报错信息：', e)

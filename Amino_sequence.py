# 氨基酸序列格式化，每行50个，每10个空格
def sequence_format(sequence):
    b = '1        '
    for i in range(len(sequence)):
        if i % 50 == 0 and i != 0:
            b += '\n' + str(i + 1)
            if 10 <= i <= 99:
                b += '       '
            elif 100 <= i <= 999:
                b += '      '
            elif 1000 <= i <= 9999:
                b += '     '
        if i % 10 == 0:
            b += sequence[i:i + 10] + ' '
    return b


# 从fasta文件提取氨基酸序列至字典中
def main():
    with open('67.fasta') as f:
        seq = {}
        for line in f:
            if line.startswith('>'):
                name = line.replace('>', '').split()[0]
                seq[name] = ''
            else:
                seq[name] += line.replace('\n', '').strip()
    c = ''
    for key, value in seq.items():
        title = key + '\n'
        seq_format = sequence_format(value) + '\n'
        c += title + seq_format
    with open('Format_sequence.txt', 'w') as d:
        d.write(c)


if __name__ == '__main__':
    main()

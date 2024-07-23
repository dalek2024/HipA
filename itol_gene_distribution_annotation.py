input_file = 'HipA_number_pergenome.tsv'
output_file = 'HipA_annotation.txt'

# 定义颜色映射
colors = {
    '1': '#FF0000',  # Red
    '2': '#00FF00',  # Green
    '3': '#0000FF',  # Blue
    '4': '#FFFF00',  # Yellow
    '5': '#FF00FF',  # Magenta
    '6': '#00FFFF',  # Cyan
}

# 打开输出文件
with open(output_file, 'w') as out:
    # 写入iTOL注释文件头部
    out.write("DATASET_COLORSTRIP\n")
    out.write("SEPARATOR COMMA\n")
    out.write("DATASET_LABEL,HipA Number Per Genome\n")
    out.write("COLOR,#000000\n")
    out.write("DATA\n")

    # 读取输入文件并处理每一行
    with open(input_file, 'r') as in_file:
        for line in in_file:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                node = parts[0]
                number = parts[2]
                color = colors.get(number)
                if color:
                    out.write(f"{node},{color}\n")

# 验证输出文件内容
with open(output_file, 'r') as file:
    generated_content = file.read()

generated_content[:500]  # 显示前500个字符以供检查

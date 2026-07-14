#!/bin/bash

# =====================================================================
# 配置区：如果你现有的脚本路径不在当前目录下，请在这里修改它们的绝对路径
# =====================================================================
GENERATE_CMD_PY="/share/nas1/renxz/pip/function/cellsplitplot/generate_cmd.py"
# =====================================================================

# 初始化变量
INPUT_FILE=""
CELLS_NPY=""
OUTPUT_DIR=""

# 打印帮助信息
usage() {
    echo "使用方法:"
    echo "  $0 --input <原始CSV> --npy <cells.npy> --output <输出目录>"
    echo "  $0 -i <原始CSV> -n <cells.npy> -o <输出目录>"
    exit 1
}

# 解析参数
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input) INPUT_FILE="$2"; shift 2 ;;
        -n|--npy) CELLS_NPY="$2"; shift 2 ;;
        -o|--output) OUTPUT_DIR="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "未知参数: $1"; usage ;;
    esac
done

# 检查参数
if [[ -z "$INPUT_FILE" ]] || [[ -z "$CELLS_NPY" ]] || [[ -z "$OUTPUT_DIR" ]]; then
    echo "错误: 缺少必填参数！"
    usage
fi

if [ ! -f "$INPUT_FILE" ] || [ ! -f "$CELLS_NPY" ]; then
    echo "错误: 输入的文件不存在，请检查路径！"
    exit 1
fi

# 创建输出目录
mkdir -p "$OUTPUT_DIR"

# -------------------------------------------------------------
# 1. 直接运行你的第一段：洗数据
# -------------------------------------------------------------
CLEANED_ALL="$OUTPUT_DIR/cleaned_all_data.csv"
echo "【步骤 1】正在调用 Awk 清洗数据..."

# 这里完全保留你原本的洗数据逻辑，直接输出到指定的输出目录
awk -F, '
BEGIN {OFS=","}
NR==1 {print $0; next}
{
    if (match($1, /cell_[0-9]+/)) {
        $1 = substr($1, RSTART, RLENGTH)
    }
    print $0
}' "$INPUT_FILE" > "$CLEANED_ALL"

echo " -> 清洗完成: $CLEANED_ALL"

# -------------------------------------------------------------
# 2. 直接调用你的第二段：generate_cmd.py
# -------------------------------------------------------------
echo "【步骤 2】正在调用 generate_cmd.py 生成绘图命令..."

# 直接调用你的 python 脚本，把清洗后的临时 CSV 路径传给它
# 并且利用 tail -n 5 捕获最后几行实际的命令（去除你的等号分割线装饰）
RAW_CMD=$(/share/nas1/dengdj/software/Miniforge/miniforge3/envs/py38/bin/python "$GENERATE_CMD_PY" \
            --cluster_csv "$CLEANED_ALL" \
            --cells_npy "$CELLS_NPY" \
            --save_dir "$OUTPUT_DIR" | grep -v "===" | sed '/^[[:space:]]*$/d')

# -------------------------------------------------------------
# 3. 自动执行绘图
# -------------------------------------------------------------
if [ -z "$RAW_CMD" ]; then
    echo "错误: 未能从 generate_cmd.py 捕获到有效的命令！"
    exit 1
fi

echo "【步骤 3】正在自动执行生成的绘图命令..."
echo -e "-------------------------------------\n$RAW_CMD\n-------------------------------------"

# 执行绘图命令
eval "$RAW_CMD"

if [ $? -eq 0 ]; then
    echo -e "\n【完成】出图成功！"
else
    echo -e "\n【错误】最终绘图脚本执行失败。"
    exit 1
fi

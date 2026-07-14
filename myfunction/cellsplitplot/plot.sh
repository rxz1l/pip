#!/bin/bash

# 默认的 Python 路径（如果外部不指定，则使用这个默认值）
DEFAULT_PYTHON="/share/nas1/dengdj/software/Miniforge/miniforge3/envs/py38/bin/python"

# --- 参数解析与帮助信息 ---
usage() {
    echo "使用方法: $0 -n <cells_npy> -c <csv_dir> -o <out_dir> [-p <python_path>]"
    echo "  -n  输入 cells.npy 文件的绝对路径"
    echo "  -c  存放 4 个 CSV 文件 (Anchor, MIA_result, RCTD, Spotlight) 的输入目录"
    echo "  -o  输出的绘图结果保存根目录"
    echo "  -p  (可选) 指定 Python 解释器路径，默认使用: $DEFAULT_PYTHON"
    exit 1
}

# 使用 getopts 解析命令行参数
while getopts "n:c:o:p:h" opt; do
    case ${opt} in
        n ) CELLS_NPY=$OPTARG ;;
        c ) CSV_DIR=$OPTARG ;;
        o ) BASE_SAVE_DIR=$OPTARG ;;
        p ) CUSTOM_PYTHON=$OPTARG ;;
        h | ? ) usage ;;
    esac
done

# 检查必填参数是否缺失
if [ -z "$CELLS_NPY" ] || [ -z "$CSV_DIR" ] || [ -z "$BASE_SAVE_DIR" ]; then
    echo "错误: 缺少必要参数！"
    usage
fi

# 确定最终使用的 Python 路径
python=${CUSTOM_PYTHON:-$DEFAULT_PYTHON}

# --- 业务逻辑开始 ---

# 定义要处理的 4 个 CSV 文件名
CSV_FILES=("Anchor.csv" "MIA_result.csv" "RCTD.csv" "Spotlight.csv")

# 循环处理每一个 CSV 文件
for csv_name in "${CSV_FILES[@]}"
do
    # 拼接完整的 CSV 路径
    FULL_CSV_PATH="${CSV_DIR}/${csv_name}"
    
    # 检查输入的 CSV 文件是否存在，不存在则跳过
    if [ ! -f "$FULL_CSV_PATH" ]; then
        echo "警告: 未找到文件 ${FULL_CSV_PATH}，跳过此文件的绘图。"
        continue
    fi

    # 提取不带后缀的文件名作为目录名（例如 Anchor, MIA_result 等）
    DIR_NAME="${csv_name%.*}"
    
    # 拼接该 CSV 专属的输出子目录
    CURRENT_SAVE_DIR="${BASE_SAVE_DIR}/${DIR_NAME}"
    
    echo "=========================================================="
    echo "正在处理文件: ${csv_name}"
    echo "输出目标目录: ${CURRENT_SAVE_DIR}"
    echo "=========================================================="
    
    # 如果输出目录不存在，则创建
    if [ ! -d "${CURRENT_SAVE_DIR}" ]; then
        mkdir -p "${CURRENT_SAVE_DIR}"
        echo "目录不存在，已成功创建: ${CURRENT_SAVE_DIR}"
    fi
    
    echo "正在分析 Cluster 并调用 cluster_plot.py 绘图..."
    
    # 提取 generate_cmd.py 打印出的核心命令并直接运行
    # 注意：运行此脚本时，请确保 generate_cmd.py 在你执行命令的当前工作目录下
    CMD=$(${python} generate_cmd.py --cluster_csv "${FULL_CSV_PATH}" --cells_npy "${CELLS_NPY}" --save_dir "${CURRENT_SAVE_DIR}" | sed -n '/\/share\/nas1\/dengdj/,/--save_dir/p')
    
    if [ -n "$CMD" ]; then
        # 执行生成的命令
        eval "$CMD"
        echo "文件 ${csv_name} 绘图完成！"
    else
        echo "错误: 未能成功生成 ${csv_name} 的绘图命令，请检查 generate_cmd.py 是否在当前目录。"
    fi
    
    echo -e "==========================================================\n"
done

echo "所有存在的 CSV 文件批量处理完毕！"

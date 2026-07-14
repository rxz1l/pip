#!/bin/bash

# ==============================================================================
# 使用说明函数
# ==============================================================================
usage() {
    echo "使用方法: $0 --input_csv <路径> --cells_npy <路径> --save_dir <路径> [--group_name <名称>]"
    echo "参数说明:"
    echo "  --input_csv  : 初始输入的原始 CSV 数据文件路径 (用于传给 split_b2c.sh)"
    echo "  --cells_npy  : 对应的 cells.npy 文件路径"
    echo "  --save_dir   : 结果输出保存的基础根目录"
    echo "  --group_name : [可选] 组名（如 AR 或 SR），默认从 input_csv 的父目录自动获取"
    exit 1
}

# ==============================================================================
# 1. 解析命令行参数
# ==============================================================================
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --input_csv) INPUT_CSV="$2"; shift ;;
        --cells_npy) CELLS_NPY="$2"; shift ;;
        --save_dir)  BASE_SAVE_DIR="$2"; shift ;;
        --group_name) GROUP_NAME="$2"; shift ;;
        -h|--help)   usage ;;
        *) echo "未知参数: $1"; usage ;;
    esac
    shift
done

# 校验核心参数是否缺失
if [ -z "$INPUT_CSV" ] || [ -z "$CELLS_NPY" ] || [ -z "$BASE_SAVE_DIR" ]; then
    echo "错误: 缺少必要参数！"
    usage
fi

# ==============================================================================
# 2. 环境配置与路径初始化
# ==============================================================================
PYTHON_ENV="/share/nas1/dengdj/software/Miniforge/miniforge3/envs/py38/bin/python"
SPLIT_SCRIPT="/share/nas1/renxz/project/DG108_s3000_Dendrobium_nobile/DZH-75551/script/cellsplit/split_b2c.sh"

# 如果未指定组名，则默认抓取输入文件所在目录的名称（对应原 AR / SR）
if [ -z "$GROUP_NAME" ]; then
    GROUP_NAME=$(basename "$(dirname "$INPUT_CSV")")
fi

# 细胞拆分脚本的 --output 参数目标目录 (即存放所有拆分出的 csv 的目录，例如 .../cellsplitplot/AR)
SPLIT_OUT_DIR="${BASE_SAVE_DIR}/${GROUP_NAME}"

echo "=========================================================="
echo " 运行环境及参数确认:"
echo " 组名 (Group):          ${GROUP_NAME}"
echo " 初始输入 CSV 文件:     ${INPUT_CSV}"
echo " 输入 NPY 文件:         ${CELLS_NPY}"
echo " 细胞拆分输出目录:      ${SPLIT_OUT_DIR}"
echo "=========================================================="

# 基础文件存在性校验
if [ ! -f "$INPUT_CSV" ]; then
    echo "错误: 初始输入的 CSV 文件不存在: ${INPUT_CSV}"
    exit 1
fi

if [ ! -f "$CELLS_NPY" ]; then
    echo "错误: 输入的 NPY 文件不存在: ${CELLS_NPY}"
    exit 1
fi

# ==============================================================================
# 3. 步骤一：执行 cell split 细胞拆分
# ==============================================================================
echo ">>> 步骤 1/2: 正在运行细胞拆分脚本 (split_b2c.sh)..."
echo "输出目录为: ${SPLIT_OUT_DIR}"

# 调用原本的第一个脚本
bash "$SPLIT_SCRIPT" \
    --input "$INPUT_CSV" \
    --output "$SPLIT_OUT_DIR"

if [ $? -eq 0 ]; then
    echo "【成功】细胞拆分步骤已完成，文件已成功生成至 ${SPLIT_OUT_DIR} 。"
else
    echo "【错误】细胞拆分脚本执行失败，请检查！"
    exit 1
fi

echo "----------------------------------------------------------"

# ==============================================================================
# 4. 步骤二：遍历拆分出来的目录，对所有 CSV 批量绘图
# ==============================================================================
echo ">>> 步骤 2/2: 开始批量绘制 ${SPLIT_OUT_DIR} 下的所有 CSV 图表..."

if [ ! -d "${SPLIT_OUT_DIR}" ]; then
    echo "错误: 拆分输出目录不存在: ${SPLIT_OUT_DIR}"
    exit 1
fi

# 优化匹配：把 cleaned_all_data.csv 直接显式写在列表最前面，后面再跟通配符
shopt -s nullglob
cluster_files=("${SPLIT_OUT_DIR}/cleaned_all_data.csv" "${SPLIT_OUT_DIR}"/cluster_*.csv)
shopt -u nullglob

if [ ${#cluster_files[@]} -eq 0 ]; then
    echo "提示: 在目录 ${SPLIT_OUT_DIR} 下未找到任何有效的 CSV 文件，跳过绘图。"
    exit 0
fi

# 循环处理目录下的每一个 csv
for full_csv_path in "${cluster_files[@]}"
do
    # 再次安全检查
    [ -f "$full_csv_path" ] || continue
    
    csv_name=$(basename "${full_csv_path}")
    cluster_name="${csv_name%.*}" # 除去后缀的文件名，如 cleaned_all_data 或 cluster_0
    
    # 最终绘图结果保存的子目录 (例如 .../cellsplitplot/AR/cleaned_all_data 或 .../cellsplitplot/AR/cluster_0)
    CURRENT_SAVE_DIR="${SPLIT_OUT_DIR}/${cluster_name}"
    
    echo "----------------------------------------------------------"
    echo "正在处理子文件: ${csv_name}"
    echo "图片保存目标目录: ${CURRENT_SAVE_DIR}"
    echo "----------------------------------------------------------"
    
    if [ ! -d "${CURRENT_SAVE_DIR}" ]; then
        mkdir -p "${CURRENT_SAVE_DIR}"
    fi
    
    # 生成并截取 python 绘图命令
    CMD=$(${PYTHON_ENV} /share/nas1/renxz/pip/function/cellsplitplot/generate_cmd.py \
        --cluster_csv "${full_csv_path}" \
        --cells_npy "${CELLS_NPY}" \
        --save_dir "${CURRENT_SAVE_DIR}" | sed -n '/\/share\/nas1\/dengdj/,/--save_dir/p')
    
    # 执行生成的命令
    if [ -n "$CMD" ]; then
        eval "$CMD"
        echo "【成功】文件 ${csv_name} 绘图完成！"
    else
        echo "【错误】未能成功生成 ${csv_name} 的绘图命令，请检查。"
    fi
done

echo "=========================================================="
echo " 整个 Pipeline 处理完成！"
echo " 拆分目录 ${SPLIT_OUT_DIR} 下的所有 CSV 已全部批量绘图完毕！"
echo "=========================================================="

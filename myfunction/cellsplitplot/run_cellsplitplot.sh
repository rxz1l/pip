#!/bin/bash

# 默认的 Python 路径（可根据需要修改）
DEFAULT_PYTHON="/share/nas1/dengdj/software/Miniforge/miniforge3/envs/py38/bin/python"

# --- 参数解析与帮助信息 ---
usage() {
    echo "使用方法: $0 -r <rds_path> -c <csv_dir> -n <cells_npy> -o <out_base_dir> [-p <python_path>]"
    echo "  -r  输入的 Seurat RDS 文件绝对路径"
    echo "  -c  包含 MIA, Anchor, Spotlight, RCTD 原始xls结果的根目录 (例如: .../result/sc_st/)"
    echo "  -n  输入的 cells.npy 文件绝对路径"
    echo "  -o  输出结果的总根目录 (脚本会自动在内部创建子文件夹)"
    echo "  -p  (可选) 指定 Python 解释器路径，默认使用: $DEFAULT_PYTHON"
    exit 1
}

# 使用 getopts 解析命令行参数
while getopts "r:c:n:o:p:h" opt; do
    case ${opt} in
        r ) RDS_PATH=$OPTARG ;;
        c ) CSV_DIR=$OPTARG ;;
        n ) CELLS_NPY=$OPTARG ;;
        o ) OUT_BASE_DIR=$OPTARG ;;
        p ) CUSTOM_PYTHON=$OPTARG ;;
        h | ? ) usage ;;
    esac
done

# 检查必填参数是否缺失
if [ -z "$RDS_PATH" ] || [ -z "$CSV_DIR" ] || [ -z "$CELLS_NPY" ] || [ -z "$OUT_BASE_DIR" ]; then
    echo "错误: 缺少必要参数！"
    usage
fi

# 确定最终使用的 Python 路径
python=${CUSTOM_PYTHON:-$DEFAULT_PYTHON}

# --- 自动定义内部中间路径 ---
# R 脚本生成的中间 CSV 存放在输出总目录下的 glass_csvs 文件夹中
MID_CSV_DIR="${OUT_BASE_DIR}/intermediate_csvs"
# 最终图片存放在输出总目录下的 plots 文件夹中
PLOT_SAVE_DIR="${OUT_BASE_DIR}/plots"


echo "=========================================================="
echo "开始运行空间转录组集成分析与绘图管线"
echo "=========================================================="

# 1. 运行 R 脚本：提取和整合 4 种算法的 Barcode 细胞类型对应表 ----
echo -e "\n[STEP 1] 正在调用 R 脚本处理分群与注释数据..."

# 确保当前目录下有之前整理好的 process_st_data.R
/share/nas1/xied/software/miniconda3/bin/Rscript /share/nas1/renxz/pip/function/cellsplitplot/sc_st_clean_data.r \
      -i "${RDS_PATH}" \
      -c "${CSV_DIR}" \
      -o "${MID_CSV_DIR}"

  


# 2. 循环读取生成的 CSV 并调用 Python 绘图 ----
echo -e "\n[STEP 2] 正在调用 Python 脚本进行批量循环绘图..."

CSV_FILES=("Anchor.csv" "MIA_result.csv" "RCTD.csv" "Spotlight.csv")

for csv_name in "${CSV_FILES[@]}"
do
    FULL_CSV_PATH="${MID_CSV_DIR}/${csv_name}"
    
    # 检查 R 脚本是否成功输出了该 CSV，不存在则跳过
    if [ ! -f "$FULL_CSV_PATH" ]; then
        echo "警告: 未找到中间文件 ${FULL_CSV_PATH}，跳过该算法的绘图。"
        continue
    fi

    # 提取算法名称作为子目录名
    DIR_NAME="${csv_name%.*}"
    CURRENT_SAVE_DIR="${PLOT_SAVE_DIR}/${DIR_NAME}"
    
    echo "----------------------------------------------------------"
    echo "正在绘图: ${csv_name}"
    echo "图片保存至: ${CURRENT_SAVE_DIR}"
    
    # 如果子输出目录不存在，则创建
    if [ ! -d "${CURRENT_SAVE_DIR}" ]; then
        mkdir -p "${CURRENT_SAVE_DIR}"
    fi
    
    # 提取 generate_cmd.py 打印出的核心命令并直接运行
    # 注意：确保 generate_cmd.py 也在当前工作目录下
    CMD=$(${python} /share/nas1/renxz/pip/function/cellsplitplot/generate_cmd.py --cluster_csv "${FULL_CSV_PATH}" --cells_npy "${CELLS_NPY}" --save_dir "${CURRENT_SAVE_DIR}" | sed -n '/\/share\/nas1\/dengdj/,/--save_dir/p')
    
    if [ -n "$CMD" ]; then
        eval "$CMD"
        echo "成功完成 ${csv_name} 的绘图任务。"
    else
        echo "错误: 未能成功生成 ${csv_name} 的绘图命令，请检查 generate_cmd.py。"
    fi
done

echo -e "\n=========================================================="
echo "🎉 所有任务执行完毕！"
echo "中间 CSV 目录: ${MID_CSV_DIR}"
echo "最终图片目录: ${PLOT_SAVE_DIR}"
echo "=========================================================="

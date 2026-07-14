import pandas as pd
import argparse
import os

def generate_cluster_command():
    # 1. 极大扩展的颜色池（共 64 种互不相同、高对比度的学术/数据可视化常用颜色）
    COLOR_POOL = [
        # 经典高对比度色板 (前20个)
        "#E64B35B2", "#4DBBD5B2", "#00A087B2", "#F39B7FB2", "#fbb612", 
        "#ba2f7b", "#8491B4B2", "#3C5488B2", "#91D1C2B2", "#DC0000B2", 
        "#7E6148B2", "#81C784", "#FFB74D", "#64B5F6", "#BA68C8", 
        "#4DB6AC", "#D4E157", "#A1887F", "#FF8A65", "#A1CE4C", 
        # 扩展高饱和度与明亮色板 (21-45)
        "#268785", "#ed1299", "#09f9f5", "#cc8e12", "#d561dd",
        "#c93f00", "#ddd53e", "#4aef7b", "#e86502", "#9ed84e",
        "#39ba30", "#6ad157", "#8249aa", "#99db27", "#e07233",
        "#ff523f", "#ce2523", "#f7aa5d", "#cebb10", "#03827f",
        "#931635", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
        # 莫兰迪/柔和及深色补色板 (46-64)
        "#1a918f", "#ff66fc", "#2927c4", "#7149af", "#57e559",
        "#8e3af4", "#f9a270", "#22547f", "#db5e92", "#edd05e",
        "#6f25e8", "#0dbc21", "#280f7a", "#6373ed", "#5b910f",
        "#7b34c1", "#0cf29a", "#d80fc1", "#51f59b"
    ]

    # 2. 参数解析
    parser = argparse.ArgumentParser(description="自动提取Cluster并生成绘图命令")
    parser.add_argument('--cluster_csv', type=str, required=True, help="输入的 cluster CSV 文件路径")
    parser.add_argument('--cells_npy', type=str, required=True, help="cells.npy 文件路径")
    parser.add_argument('--save_dir', type=str, required=True, help="结果输出目录")
    args = parser.parse_args()

    # 3. 读取 CSV 并提取去重后的 Cluster 列表
    if not os.path.exists(args.cluster_csv):
        print(f"错误: 找不到文件 {args.cluster_csv}")
        return

    # 读取并自动去掉引号
    df = pd.read_csv(args.cluster_csv, quotechar='"')
    
    # 获取去重后的群组列表
    raw_clusters = df['Cluster'].dropna().unique()
    
    # 尝试按数字大小排序，否则按字符串字母排序
    try:
        clusters = sorted(raw_clusters, key=lambda x: int(x))
    except ValueError:
        clusters = sorted(raw_clusters)

    # 4. 动态分配颜色
    color_mapping = []
    for idx, cluster in enumerate(clusters):
        # 即使 Cluster 数量超过 64 个，% 机制也能确保安全运行（循环复用）
        color = COLOR_POOL[idx % len(COLOR_POOL)]  
        color_mapping.append(f"{cluster}:{color}")
    
    cluster_and_color_str = ",".join(color_mapping)

    # 5. 拼接最终的运行命令
    command = f"""/share/nas1/dengdj/software/Miniforge/miniforge3/envs/py38/bin/python /share/nas1/fanlp/03.tools/UMICountDistPlot/cluster_plot.py \\
                --cells_path {args.cells_npy} \\
                --cluster_path {args.cluster_csv} \\
                --cluster_and_color "{cluster_and_color_str}" \\
                --background_color '#FFFFFF' \\
                --save_dir {args.save_dir}"""

    print("\n" + "="*20 + " 自动生成的命令如下 " + "="*20)
    print(command)
    print("="*60 + "\n")

if __name__ == "__main__":
    generate_cluster_command()

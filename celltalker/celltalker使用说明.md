# 1.参数说明
## 1.1参数列表
```r
/public/work/Personal/renxiaozhen/software/miniconda3/envs/R/bin/Rscript /public/work/Pipline/celltalker/celltalker_data.r -h


Options:
        -d DATA, --data=DATA
                数据路径

        -l LIGAND_RECEPTOR_PAIRS, --ligand_receptor_pairs=LIGAND_RECEPTOR_PAIRS
                配受体矩阵

        -s SPLIT, --split=SPLIT
                数据按照某一指标切分（如tissue）,若不需要切分可不输入该l参数

        -g METADATA_GROUPING, --metadata_grouping=METADATA_GROUPING
                Serurat数据中定义细胞类型的列名

        -n NUMBER_CELLS_REQUIRED, --number_cells_required=NUMBER_CELLS_REQUIRED
                进行配体/受体相互作用分析,每个分组所需的细胞数。默认为100

        -i MIN_EXPRESSION, --min_expression=MIN_EXPRESSION
                在相互作用分析中考虑配体或受体的最小表达计数,默认为1000

        -a MAX_EXPRESSION, --max_expression=MAX_EXPRESSION
                在相互作用分析中考虑配体或受体的最大表达计数,默认为2000

        -t SCRAMBLE_TIMES, --scramble_times=SCRAMBLE_TIMES
                配体/受体相互作用的随机排列次数,默认为10

        -P PLANTPHONEDB, --PlantPhoneDB=PLANTPHONEDB
                是否使用plantphoneDB进行分析，默认为不使用(TRUE/FALSE)

        -m METHOD, --method=METHOD
                plantPhoneDB的用法，默认为Average(LRscore/WeightProduct/Average/Product)

        -c MIN.PCT, --min.pct=MIN.PCT
                只测试在两个群体中任何一个群体中最小比例的细胞中检测到的基因。默认值是0.1

        -e ITERATIONS, --iterations=ITERATIONS
                如果方法设置为Average或Product，则使用排列测试通过随机打乱聚类标签来计算配体-受体相互作用得分。默认值为100

        -o OUTPUT_DIR, --output_dir=OUTPUT_DIR
                输出路径,默认为当前目录

        -h, --help
                Show this help message and exit
```
## 1.2 `-d`输入为`RDS`文件，为`seurat`对象
**RDS文件要求：**
> 1. 要包含细胞注释列，列名作为`-g`的输入。
> 2. 细胞注释列中的细胞名称不能有`_`，若有建议将`_`替换为`.`。
## 1.3 `-l`输入为`Rdata`文件，为配受体矩阵
本脚提供了<big><font color=red>rice</font></big>的配受体数据见目录下`rice_db.Rdata`。
**自定义配受体矩阵，矩阵结构如下：**
|ligand|receptor|pair|
|:----:|:----:|:----:|
|配体|受体|配体_受体|
配体受体<font color=red>命名规则</font>要与`seurat`对象一致。
# 2.输出目录
输出的目录结构为：
```bash
.
└── celltalk_output
│   ├── all
│   │   ├── all.pdf
│   │   └── top_stats_all.Rdata
│   └── Tissue
│       ├── Bud.pdf
│       ├── data_split.Rdata
│       ├── Flag.pdf
│       ├── leaf.pdf
│       ├── Root.pdf
│       ├── SAM.pdf
│       ├── Seed.pdf
│       ├── SP.pdf
│       └── ST.pdf
└── plantphone_output
    ├── all
    └── Tissue
```
该脚本会在输入的路径下创建一个`celltalk_output`文件夹,该文件夹下`all`为整体数据`celltalker`输出，其中`all.pdf`为整体弦图，`top_stats_all.Rdata`为基于整体数据分析得到数据，可以用于后续的作图；`Tissue`是以`-s`传入的分割数据列的名称创建的文件夹，其中`data_split.Rdata`为基于分隔数据分析得到的数据，可以用于后续作图，该文件夹下的`pdf`文件为根据`-s`传入列名进行分组，各个分组的弦图。`-P`为是否使用plantphoneDB进行分析，当为`TRUe`时会在输入的路径下创建一个`celltalk_output`文件夹。可以通过`-m`设置score计算方法，默认为`Average(LRscore/WeightProduct/Average/Product)`，其中`Average`计算方法与`celltalker`类似。

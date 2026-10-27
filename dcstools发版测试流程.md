# dcstools发版测试流程

## 测试分层

### 1.快速冒烟测试

MR / 提交触发，轻量级别（快速确认核心流程运行正常）

小数据集，核心参数组合，快速跑完，只检查能否正常跑完、文件是否生成，不精细比对

### 2.回归测试

发版前必跑，中等

列出的全部 case，对比 golden 基准，校验产物完整性 + 结果一致性

### 3.大容量压力测试

版本大变更才跑，比如算法重构、minibwa/mem2 升级

4G 以上基因组、超大 WGS 数据，测内存、IO、运行时长

## dcstools工具集各功能模块测试

### 1.wgs/wes功能测试

目标是能够覆盖不同基因组、倍性、数据量、比对算法参数

#### 测试数据集

| 类别  | 样本  | 说明  |
| --- | --- | --- |
| 人样本 | WES（bed输入） |     |
|     | 0.4x小数据 |     |
|     | 30x大数据 |     |
| 非人样本 | 动植物数据，小基因组 |     |
|     | 动植物数据，大基因组（4G以上吧） |     |
| 多倍体 | 土豆  |     |

**检查项**

1.  程序正常退出，exit code=0，无报错堆栈
    
2.  产物完整性：生成 BAM、索引 bai、vcf、report 文件
    
3.  文件大小校验：和 golden 基准偏差在允许阈值内（不能完全相等，时间戳、随机数会有微小差异）
    
4.  BAM 一致性校验：
    
    -   比对统计：总 reads、mapped 率、properly paired 率、插入片段分布
        
    -   抽样比对记录（可选：md5 不建议直接比对整个 BAM，会失败）
        
5.  资源监控：峰值内存、walltime（防止新版本内存暴涨）
    

### 2.mutect2变异检测测试

两种模式：Tumor-Normal 配对；单样本 no pair

#### 测试数据集

| 类别  | 样本  | 说明  |
| --- | --- | --- |
| 人   | tumor + normal pair |     |
|     | no pair |     |

**检查项**

1.  正常跑完，生成 vcf/gvcf
    
2.  变异总数、qual 分布、过滤后变异计数对比 golden
    
3.  关键变异位点抽样比对（不建议全 vcf md5）
    

### 3.jointcall联合基因分型测试

有测试脚本运行一遍

### 4.压缩测试

有测试脚本运行一遍

## 结果验证

-   需要保留之前正确版本的运行结果作为基准对比
    
-   测试运行后输出产物（关键结果文件bam，gvcf，vcf）是否存在、完整性如何、文件大小、**结果的一致性**对比
    

## 注意事项

各工具的关键参数需要覆盖，比如minibwa、mem2、ert都要运行一次

#### 测试流程目录结构

```text
dcstools-test-suite/
├── cases/                  # 测试用例配置（每个case：参数、参考基因组、样本）
│   ├── wgs_human_small
│   ├── wgs_human_large
│   ├── wes_human
│   ├── wheat_large_genome
│   ├── potato_polyploid
│   ├── mutect2_tn
│   ├── mutect2_single
│   ├── jointcall
│   └── compress
├── golden/                 # golden基准：统计结果、关键位点，可以不存原始大BAM
├── scripts/                # 校验脚本：samtools_stats_diff.py、vcf_compare.py
├── data/                   # 测试数据集配置（路径指向集群/对象存储大数据）
└── reports/                # 每次运行输出测试报告
```
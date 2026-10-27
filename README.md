# dcstools 发版测试套件

对一份 `dcs-tools` 发布包做回归：跑真实数据上的对齐、calling、基因分型，再和金标比 BAM / gVCF / VCF。一个 case 失败不会停掉后面的 case，整轮结束给一个总判定。

## 现在完善到什么程度

人、鼠、稻、薯的 WGS，人的 WES，体细胞 Mutect2，联合分型，SeqArc，以及冒烟用的小基因组，都有 case、金标和报告，可以当发版门禁。

还没收口的有四件事：

- `cycas_large_genome`、`wheat_large_genome` 脚本还在，但暂时不进 `full` / `heavy`。要跑时点名：`./run.sh cycas_large_genome wheat_large_genome`。
- `potato_polyploid` 仍在 `full` / `heavy` 里，目前没有汇总报告（数据或金标不齐时会自己跳过或报错）。
- `aligner_shared_mem` 有脚本，也跑过一次，但不在任何默认 tier 里，要按名字点名才跑。
- 30x 和 30x ERT 的 `run.log` 仍停在当时那次 32 线程跑完的 FAIL。报告里的 PASS 来自之后用单线程 genotype 重算并改过的 `case.json`，不是把整条流水线又跑了一遍。

## 跑之前

配置在 `env/cases.conf`，命令行环境变量可以覆盖。默认值就是这台机器上的路径。

| 变量 | 作用 |
|---|---|
| `DCS_HOME` | 被测的发布包，里面要有 `bin/dcs` |
| `PROXY_HOST` | license 代理，默认 `127.0.0.1:8909`。没开的话每个工具都会直接失败 |
| `DATA_ROOT` / `WGS30X_DATA` / `MOUSE_DATA` / `RICE_DATA` / `POTATO_DATA` | 各物种的 FASTQ 和参考 |
| `GOLDEN_ROOT` | 金标，一个 case 一个子目录 |
| `WORK_ROOT` | 本次产物（BAM、VCF） |
| `REPORT_ROOT` | 报告，默认仓库里的 `reports/` |
| `THREADS` | 对齐、calling 等步骤的线程，默认 32 |
| `GENOTYPER_THREADS` | genotyper 线程，默认 **1**。多线程时 QD 的随机抖动跟调度有关，同一份 gVCF 两次结果对不上 |

license 代理先起来，例如：

```bash
"$DCS_HOME/bin/dcs" licproxy --license <license 文件>
```

## 怎么运行

在套件根目录：

```bash
./run.sh --list                 # 有哪些 tier、哪些 case
make smoke                      # 打包检查 + CLI 契约 + wgs_human_small
make static                     # 只做打包和 CLI，不碰数据
make heavy                      # 30x / 鼠 / 稻 / 薯（不含苏铁、小麦）
make full                       # static + 体细胞 + WES + jointcall + SeqArc + heavy
make check                      # 只检查套件脚本自身的语法

./run.sh wgs_human_30x          # 只跑点名的 case，可写多个
./run.sh --tier full --fail-fast
```

退出码：`0` 全部通过，`1` 有检查失败，`2` 跑不起来（缺数据、没 license、工具崩溃）。

单个 case 也可以直接跑，参数只在 case 脚本上生效，`./run.sh` 不转发 `--compare-only`：

```bash
bash cases/wgs_human_small.sh
bash cases/wgs_human_30x.sh --compare-only    # 不重跑流水线，只拿现有产物对金标
bash cases/wgs_mouse.sh --update-golden       # 用这次产物换掉该 case 的金标
```

`--update-golden` 是显式动作。不要在失败时自动刷新金标，否则回归就不再抓得到回退。`make golden` 等于 `./run.sh --tier full --update-golden`。

## 测了哪些

`make full` 的内容，按报告里的分组：

**冒烟与工具契约。** 几分钟内结束。

| case | 测什么 |
|---|---|
| `bundle_audit`、`cli_contract_check` | 发布包文件和 CLI 参数是否还在（`static` tier，不进 case 卡片） |
| `wgs_human_small` | 人小基因组，FASTQ 到 VCF 的整条胚系流程 |
| `seqarc_roundtrip` | SeqArc 压缩后再解开，逐字节一致 |
| `aligner_shared_mem` | aligner 共享内存索引的创建和挂载。不在默认 tier 里 |
| `jointcall` | 多样本 gVCF 联合分型（jointcaller） |

**胚系 WGS，二倍体。** 同一行是同一套数据，三列只换对齐索引：默认 BWA、ERT（`--use-high-memory-index`）、minibwa。

| | 默认索引 | ERT | minibwa |
|---|---|---|---|
| 人 30x | `wgs_human_30x` | `wgs_human_30x_ert` | `wgs_human_30x_minibwa` |
| 小鼠 GRCm39 | `wgs_mouse` | `wgs_mouse_ert` | `wgs_mouse_minibwa` |
| 水稻 IRGSP-1.0 | `wgs_rice` | `wgs_rice_ert` | `wgs_rice_minibwa` |

非人物种不做 BQSR。minibwa 与默认索引、ERT 一样对上金标：人、鼠、稻比较 BAM、gVCF 和 genotyper VCF，马铃薯比较 BAM 和 gVCF。

**外显子组。** `wes_human`，流程与 WGS 相同，calling 限制在目标 BED 内。

**多倍体。** 马铃薯四倍体，一行三列：`wgs_potato`、`wgs_potato_ert`、`wgs_potato_minibwa`。基因分型走 jointcaller。`potato_polyploid` 是另一条 jointcaller 路径，在 heavy tier 里，目前没有汇总报告。

**体细胞。** `mutect2_tn`：肿瘤 / 正常配对，Mutect2 再加 FilterMutectCalls。不比较 BAM。

**大基因组。** `cycas_large_genome`、`wheat_large_genome` 走 ERT 索引，暂时不在 `full` / `heavy` 里。

每个胚系 case 的步骤一般是：建索引（已有则跳过）→ aligner →（人）BQSR → variantCaller 出 gVCF → genotyper 出 VCF → bam-stat / vcf-stat → report。然后检查产物大小，并把 BAM、gVCF、VCF 对上金标。

比较尺度：

- WGS / WES 的 gVCF 和 VCF 要求 INFO 与 FORMAT 逐字节一致。
- Mutect2 的两条 VCF 松一档：位点、等位基因、FILTER、基因型必须相同；QUAL / INFO / AD / DP 等注释允许抖动，但这类记录不能超过 **0.1%**。位点增删或基因型变化仍然失败。

## 怎么查看

跑完后用浏览器打开：

```text
reports/index.html
```

这是所有 case **最新一次**结果的汇总。页面上方是通过 / 失败计数，下面按上面的分组排列。胚系 WGS 和多倍体是「物种 × 索引」矩阵，点格子打开该项详情；失败项默认展开。

单次 suite 运行的页面仍在下面这两个文件里。汇总页底部的运行历史表暂时不显示。

```text
reports/suite/<run_id>/report.html
reports/suite/<run_id>/suite.json
```

单个 case 的判定和检查明细：

```text
reports/<case>/<run_id>/case.json
reports/<case>/<run_id>/run.log      # 总日志，含每步命令和工具标准输出
reports/<case>/<run_id>/<step>.out   # 该步骤工具自己的输出
```

`run.log` 写到 `=> <case>: PASS|FAIL` 才算这份运行的日志完整。事后只重算 genotype 或只改比较规则时，改的是 `case.json` 和 HTML，不会重写当时的 `run.log`。

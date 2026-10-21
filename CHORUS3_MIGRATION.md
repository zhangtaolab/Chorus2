# Chorus3 迁移方案（Python → Rust）

> 状态：规划文档，未开始开发。
> 范围约定（2026-10-17）：**首版不含 primer3 / 寡核苷酸热动力计算**（Tm、hairpin、dTm 过滤整体移除，等价于现有 `--skipdtm` 行为）；第一步目标是实现**基本的探针设计功能**（对应 Chorus2.py + ChorusGUI.py 主流程）。
> 目标平台：Windows、macOS、Linux 三平台原生运行。

---

## 1. 现状盘点：Chorus2 有什么

### 1.1 入口脚本（8 个，共约 4200 行）

| 脚本 | 功能 | 依赖的外部工具 | 迁移优先级 |
|---|---|---|---|
| `Chorus2.py` (753 行) | **核心设计管线**：k-mer 过滤 + 比对过滤 + 简单过滤 + 间距选择，输出 BED | jellyfish、bwa | **P1（第一步）** |
| `ChorusGUI.py` (894 行) | PyQt5 三页向导 GUI（程序页/参数页/运行页），QThread 跑管线 | 同上 | **P2** |
| `ChorusPBGUI.py` (613 行) | 探针浏览器：按染色体查看探针密度图（matplotlib），框选增删探针，导出 | 无（纯 pandas/matplotlib） | P3 |
| `ChorusNGSfilter.py` (352 行) | 用 NGS reads 的 k-mer 计数给全基因组打分，输出 BigWig | jellyfish | P4 |
| `ChorusNGSselect.py` (170 行) | 按分位数 + 间距 + 链交替筛选探针 | 无 | P4 |
| `ChorusNoRef.py` (746 行) | 无参考基因组流程：多物种 NGS reads → consensus 探针 | bwa、samtools、bcftools、jellyfish、pysam、pybedtools | P5 |
| `ChorusHomo.py` (505 行) | 探针在近缘物种基因组上的同源性注释 | bwa、samtools faidx | P5 |
| `ChorusDraftPrebuild.py` (91 行) | 短序列合并成 >1Mb 长序列加速搜索（N70 分隔） | 无 | P1 顺带 |

### 1.2 Choruslib 模块（约 2800 行）

| 模块 | 作用 | Rust 替代 |
|---|---|---|
| `jellyfish.py` (910 行) | 封装 jellyfish 子进程：count（canonical, -L 2）、query 交互式管道、滑窗 k-mer 打分 | **rustkmer**（库直接调用，去掉子进程） |
| `bwa.py` (444 行) | 封装 bwa 子进程：index、`mem -O 0 -B 0 -k 5`、SAM 解析（AS/XS 过滤） | **RustAlign**（v1.0 BWA-MEM 移植） |
| `prefilter.py` | N/同聚物过滤（连续 5 个相同碱基） | 内置字节扫描 |
| `primer3_filter.py` | Tm/hairpin/dTm 过滤 | **首版移除**（保留接口位，后续可接 `primer3` crate） |
| `spgenome.py` | >1G 基因组按染色体切成 <1G 子文件 | 内置（或改流式，无需落盘） |
| `loadfa.py` / `revcom.py` | FASTA 读取 / 反向互补 | needletail / 内置 |
| `pbsearch.py` / `blat.py` / `exonerate.py` | BLAT gfServer / exonerate-server 流程（**GUI 中已禁用，遗留代码**） | 不迁移 |
| `bamdepth.py` / `Baminfo.py` / `openBam.py` / `bcftools.py` | pysam 读 BAM 深度、bcftools consensus | noodles + rust-htslib 路线（P5） |
| `probecompare.py` | 多物种探针逐位比较求 consensus | 内置（纯逻辑） |
| `chorusio.py` / `subprocesspath.py` | IO 小工具 / Windows 路径空格处理 | 不需要（无子进程后消失） |

### 1.3 核心设计管线（P1 要复刻的全部逻辑）

```
输入: -g genome.fa -i input.fa -l 45 --homology 75 --step 5 -t N

1. genomesize = int(文件大小/1e6) Mb
   k = max(17, floor(log4(genomesize)) + 1)        # 经验公式
2. jellyfish count --canonical -m k -L 2 genome → kmer.jf
3. bwa index genome（拷贝到输出目录）
4. 滑窗 k-mer 过滤（按 ≤10Mb 分块，multiprocessing 并行）:
   对每条序列预计算逐位置 k-mer 计数（裁剪: 0→0, 2→1, ≥3→2），
   滑窗探针: kmerscore = sum(window[st : st+len-1])   # 注意故意排除最后一个 k-mer 起点
   score < maxkmerscore → 候选探针，右移 step；否则右移 1
   CLI 版 maxkmerscore = int(((len × homology/100) - k) × ploidy/2 + 0.5)
   GUI 版 maxkmerscore = int(len × homology/100) - k          # 两处公式不一致，见 §7
5. 候选写临时 FASTA → bwa mem -O 0 -B 0 -k 5 → SAM 流解析:
   保留 AS ≥ probe_length 且 XS < len × homology/100 的探针（取回 seq/chr/start）
6. atcg_filter: 含 N 或任一碱基连续 ≥5 → 剔除
7. [primer3 dTm 过滤 —— 首版跳过]
8. 按染色体贪心间距选择: 探针间距 > (primer长度 or 5) + probe_length 才保留
9. 输出: <input>.bed（间距选择后）、<input>_all.bed（全部）、<input>.len（染色体长度）
```

### 1.4 GUI 现状

- **ChorusGUI**：QToolBox 三页向导 —— ①程序页（jellyfish/bwa 路径检测、版本显示）②参数页（文件选择、threads/length/homology/dTm/R primer）③运行页（QTextBrowser 日志 + 进度条 + 停止按钮）。管线跑在 QThread + multiprocessing.Pool，信号回传进度。
- **ChorusPBGUI**：加载 `.bed` + `.len`，matplotlib 画每条染色体的探针密度/散点，SpanSelector 框选区域，添加/删除探针，总览 dock，另存项目。pandas 做数据层。

---

## 2. Rust 资产盘点：rustkmer 与 RustAlign

### 2.1 rustkmer v0.5.0（MIT）—— jellyfish 的完整替代

- `KmerCounter`：DashMap 并发哈希，2-bit 打包（k≤32 用 u64，33–64 用 u128），canonical 模式，rayon 并行，`-L/-U` 计数过滤（jellyfish 兼容参数名），u32 饱和计数。
- `.rkdb` 二进制库格式：排序落盘 → 二分查询（官方数据 22,700+ q/s，比未排序快 384–1526×），支持 exact/prefix/fuzzy 批量查询。
- IO 层现成：FASTA/FASTQ 流式 + gz/bz2/xz 透明解压 + mmap。
- 有 pyo3 绑定（`pyrustkmer`），可选作为过渡期给旧 Chorus2 加速的路径，但本方案直接作库依赖。
- **与 jellyfish 的差异（重要）**：
  - 计数是**精确、全内存**的（DashMap），没有 jellyfish 的 Bloom Filter 前置过滤。rustkmer 文档口径约 96 B/kmer（内存路线）；对 17-mer 的大基因组（如小麦 17 Gbp，理论 unique k-mer 数百亿）**内存会爆**。
  - `.rkdb` 与 `.jf` 格式不兼容（无所谓，内部格式）。
  - k > 64 不支持（Chorus2 的 k 由基因组大小推出，水稻 ~17、玉米 ~19、小麦 ~23，远小于 64，够用）。
- **结论**：P1 直接 `path` 依赖接入；大基因组内存问题列为风险 R-1（§8），备选方案是在 rustkmer 前加 Bloom Filter 或分染色体计数后 merge（其 streaming_merge 已具备）。

### 2.2 RustAlign（MIT）—— bwa 的替代

- **v1.0-final tag（`d147799`）是功能完整的 BWA-MEM 移植**：`index` + `mem`（SE/PE）、SAM 输出、641+ 测试通过、小基因组上索引与 C BWA 字节一致。
- 当前 HEAD（v2.0 分支）是重构骨架，`index/mem` 是 stub，但 `ksw_oracle`（Smith-Waterman 三件套）、mapq、SAM/coords 模块是真实可用代码。
- **迁移策略：从 `v1.0-final` tag 用 `git show v1.0-final:<path>` 提取代码**，按 RustAlign 仓库自己的 donors 规则（只搬代码、不做链接），清理其中的裸指针。
- Chorus2 只用了 bwa 的一个极小子集：`bwa index`、`bwa mem -O 0 -B 0 -k 5`（零错配/零 gap 罚分 → AS≈比对上的碱基数）、解析 AS/XS、`@SQ` 取染色体长度。RustAlign v1.0 全覆盖，且 SAM 解析可以直接用 noodles 或直接消费其内部 Alignment 结构（跳过 SAM 文本序列化，更快）。
- **注意**：v2 的 ksw 是标量实现（SIMD 排期在 Phase 22），探针级短序列（45 bp）比对量小，性能不是瓶颈。
- BLAT/exonerate 路线（`pbsearch.py`）在现有 GUI 中已删除，不迁移；如未来需要，可用 FM-index + ksw 组出 seed-extend 扫描器。

### 2.3 Windows 实测验证（2026-10-17，本机 Rust 1.92 / MSVC）

对两个项目做了真实构建 + 测试 + 功能验证，结论：**均可用，但有明确的适配清单**。

**rustkmer：核心路径可用 ✓**

- 库 + CLI 在 Windows 构建 ✓；`cargo test --lib` 236/237 通过，15 个集成测试文件中 14 个全绿（含 golden SHA256、round-trip、dense 并行计数）。
- **受控 k-mer 实验全部符合 jellyfish 语义**：同链 2 次→计数 2；3 次→3；正链 1 次 + 反链 1 次→canonical 合并为 2；仅 1 次→被 `-L 2` 过滤（查询无返回）；**对反向互补形式的查询同样返回正确计数（链方向透明）**——Chorus2 的正向链查询可直接复用。
- Windows 缺陷（均集中在 merge 功能，Chorus2 核心管线不使用）：
  - `tests/merge_cleanup_tests.rs` 使用 `std::os::unix::fs::symlink` 无门控，**该测试文件在 Windows 编译失败**；
  - `merge_bounded_memory_tests` 2 败（Windows 文件锁语义致 mtime 更新 PermissionDenied；`/proc/self/status` RSS 不可用）、lib 测试 1 败（readonly 目录在 Windows 不阻止建文件，测试假设失效）；
  - 两个 `#[cfg(unix)]`（默认内存上限、RSS 读取）有回退但回退为 0/估算值。
  - CI 无 Windows 矩阵。
- 接入要求：只走 count/query 路径；merge 相关 Windows 修复（约 4 处）在 Chorus3 内自修或提上游；把 Windows 加入 CI。
  **已提上游 issue**：[rustkmer#13](https://github.com/zhangtaolab/rustkmer/issues/13)（symlink 编译失败）、[#14](https://github.com/zhangtaolab/rustkmer/issues/14)（3 个 merge 测试 Windows 失败）、[#15](https://github.com/zhangtaolab/rustkmer/issues/15)（CI 加 Windows 矩阵）。

**RustAlign v1.0-final：端到端可用 ✓（3 个适配点）**

- 在 Windows 编译 ✓（98 个警告，均为 unused 类，源自 C 移植）。`v1.0-final`（`d147799`）经独立克隆验证。
- **真实管线烟雾测试通过**：40 kb 双染色体基因组 → `index`（SA-IS 路径）→ `mem -o 0 -b 0 -k 5 -t 4`（Chorus2 精确参数）→ SAM 输出逐项正确：
  - 精确探针：位置/CIGAR/NM/MD/**AS:i:45** 全对；
  - 5 错配探针：AS:i:40 = 40 匹配×1 + 5 错配×0，**与 C bwa 的 `-B 0` 打分语义一致**；
  - 双拷贝探针：**AS:i:45 + XS:i:45 + XA 备选标签**——Chorus2 的 homology 过滤（XS < len×homology%）正依赖此，行为正确；
  - 不存在探针：flag 4 未比对。
- 适配点：
  1. **`mem` 只接受 FASTQ，拒绝 FASTA**（C bwa 的 kseq 两者皆收；Chorus2 正是喂 FASTA 探针）→ 需修 `mem/read.rs` 的 bseq_read 或喂料前转换 FASTQ；
  2. index 前缀**剥掉扩展名**（`genome3.fa` → `genome3.bwt`；C bwa 是 `genome3.fa.bwt`）→ 内部使用无妨，库接入时需按其前缀规则取索引；
  3. 测试套件在 Windows：约 20 个测试因硬编码 `/tmp` 路径失败（`pac_io`×15、`try_parse_ann_file`×2、`bseq_read`×3 等，纯测试基建问题）；1 个测试（`mem::collect::tests::test_mem_collect_sorted`）向裸指针代码喂**伪造 FM-index 数据**（全零 BWT）稳定触发 ACCESS_VIOLATION——说明移植代码对畸形索引数据无防护（真实数据路径未发现问题），库化封装时应做输入校验。
- Linux 侧证据：仓库自带与 C BWA 的索引字节级一致性和 SAM parity 对比测试（CI 曾通过）。
- 接入要求：从 `v1.0-final` 提取后①修 FASTA 输入；②统一前缀规则；③给裸指针热点（rope/FM-index 遍历）加边界校验；④Windows 加入 CI。
  **已提上游 issue**：[RustAlign#1](https://github.com/zhangtaolab/RustAlign/issues/1)（mem 拒绝 FASTA 输入）、[#2](https://github.com/zhangtaolab/RustAlign/issues/2)（索引前缀剥扩展名）、[#3](https://github.com/zhangtaolab/RustAlign/issues/3)（约 20 个测试硬编码 /tmp）、[#4](https://github.com/zhangtaolab/RustAlign/issues/4)（伪造索引数据触发 UB 段错误）。

---

## 3. Chorus3 架构设计

### 3.1 Cargo Workspace 布局

```
chorus3/
├── Cargo.toml                 # [workspace]
├── crates/
│   ├── chorus-core/           # 纯算法库：管线编排、参数、探针选择、BED/len 输出
│   ├── chorus-io/             # FASTA 索引读写(.fai)、BED/len/bigwig IO
│   ├── chorus-kmer/           # 封装 rustkmer：count + 滑窗打分（含计数裁剪逻辑）
│   ├── chorus-align/          # 封装 RustAlign(v1.0)：index + mem + AS/XS 过滤
│   ├── chorus-cli/            # clap CLI：chorus3 design 子命令（P1）
│   └── chorus-gui/            # egui 应用（P2/P3）
├── rustkmer/                  # path 依赖（git submodule 或 vendored）
├── rustalign/                 # v1.0 提取源
└── .github/workflows/ci.yml   # 三平台矩阵
```

关键原则：

1. **`chorus-core` 不依赖任何 GUI/CLI crate**，管线全部实现为可取消的异步任务（进UCT via channel），CLI 和 GUI 只是两个壳。这直接复刻并改良现有 `ChorusGUI.RunChorus(QThread)` 的架构。
2. **外部子进程归零**：jellyfish/bwa 调用全部换成 rustkmer/RustAlign 库调用。好处：Windows 下不再有路径空格/`shutil.which`/`ps -A` 杀进程等一堆平台 hack（`subprocesspath.py` 整个模块消失），进度报告从"解析子进程 stdout"变成真实回调。
3. **并行模型**：rayon 全局线程池替代 multiprocessing.Pool；GUI 侧用 `crossbeam-channel` 回传 `(阶段, 进度, 消息)`，`AtomicBool`/`CancellationToken` 做停止按钮。

### 3.2 P1 管线映射（设计功能）

| Chorus2 步骤 | Chorus3 实现 | 所在 crate |
|---|---|---|
| `pyfasta.Fasta` 随机访问 | needletail 解析 + 自建 `.fai` 式索引（或 noodles-fasta） | chorus-io |
| k-mer 大小公式 | 原样复刻 `max(17, log4(size_mb)+1)` | chorus-core |
| `jellyfish count --canonical -L 2` | `rustkmer::KmerCounter::new(k).canonical()...` + 落盘 `.rkdb` | chorus-kmer |
| 10Mb 分块滑窗 + Pool 并行 | rayon `par_iter` 按块并行；计数裁剪 0/1/2 原样复刻 | chorus-kmer |
| `bwa index` | RustAlign v1.0 index 构建函数 | chorus-align |
| `bwa mem -O 0 -B 0 -k 5` + AS/XS 过滤 | RustAlign v1.0 mem（同参数），直接消费 Alignment 结构取 AS/次优分 | chorus-align |
| `atcg_filter` | 字节扫描（N 检测 + 游程 ≥5） | chorus-core |
| primer3 dTm | **跳过**（保留 `ProbeFilter` trait 的一个 no-op 实现，未来挂 primer3 crate） | chorus-core |
| 染色体贪心间距选择 | 原样复刻 | chorus-core |
| 输出 .bed/_all.bed/.len | chorus-io | chorus-io |

### 3.3 GUI 框架建议

**推荐：egui（eframe）+ egui_plot。** 理由：

| 维度 | egui/eframe | 说明 |
|---|---|---|
| 跨平台 | Windows/macOS(含 ARM)/Linux + Web | 纯 Rust，无系统 GUI SDK 依赖；0.34 起默认 wgpu 后端 |
| 与本项目的契合 | 工具型科学应用 + 流式进度更新 | immediate mode 天然适合"管线日志实时刷 + 进度条 + 取消" |
| 图表 | egui_plot（官方作者维护，折线/散点/框选） | 覆盖 ChorusPBGUI 的染色体密度图需求；BoxZoom/Drag 天然支持 |
| 后台任务 | 管线在线程池，每帧 `try_recv` channel | 比 QThread 信号更简单 |
| 成熟度 | 2026 年各调查中均为生产级首选之一；无障碍/IME 可用 | 社区活跃，文档齐全 |
| 代价 | 非原生控件外观；**CJK 字体需手动内嵌**（中文 UI 必须带字体文件，用 `egui::FontDefinitions` 注册） | 可接受 |

**备选方案对比**（若你想要其他取向）：

| 框架 | 取向 | 优点 | 不选的理由 |
|---|---|---|---|
| **iced** | Elm 架构声明式 | 纯 Rust、设计良好 | 大应用状态抽象有摩擦、图表生态弱于 egui_plot |
| **Slint** | 声明式 DSL + 原生感观 | 专利桌面应用免费、工具链好 | DSL 学习成本；科学图表仍要自绘 |
| **Tauri** | 系统 WebView + 前端 | JS 图表生态（ECharts 等）最强 | 引入 Node/前端工具链，与"整个项目 Rust 化"目标相悖；体积/依赖重 |
| **Dioxus** | React 风格 | 生态热 | 桌面端仍走 WebView，开发重心转向 Blitz，动荡 |
| **GPUI + GPUI Component** | Zed 系 | 性能好 | 裸 API 粗糙、版本混乱，尚不适合工具型应用 |

结论：**egui 是当前最稳的选择**，与你最初的意向一致；若未来要做很复杂的交互式基因组浏览器再评估 Tauri/JS 路线。注意 egui_plot 官方声明"寻找维护者"，但功能已覆盖本项目需求，风险低（必要时可 fork 或换 plotters + 自绘）。

### 3.4 GUI 信息架构（P2/P3）

- **P2 = ChorusGUI 等价物**：三页向导（参数表单 + 文件选择 `rfd` + 运行页日志/进度/取消）。dTm 控件省略。
- **P3 = ChorusPBGUI 等价物**：左栏染色体列表 + 主区 egui_plot 密度散点图（拖拽缩放、框选）+ 探针表 + 项目保存。数据层用 `polars` 或纯 `Vec<Probe>`（数据量 ≤ 百万行，纯 Vec 足够）。

---

## 4. Rust 依赖包清单

### 4.1 P1（CLI 核心设计功能）

| crate | 版本基线 | 用途 | 替代的 Python/外部依赖 | 许可证 |
|---|---|---|---|---|
| `clap`（derive） | 4.x | CLI 参数 | argparse | MIT/Apache-2.0 |
| `anyhow` + `thiserror` | 1.x / 2.x | 错误处理 | — | MIT/Apache-2.0 |
| `needletail` | 0.5+ | FASTA/FASTQ 流式解析 | pyfasta | MIT |
| `rustkmer`（path 依赖） | 0.5.0 | k-mer 计数与查询 | jellyfish 子进程 | MIT |
| RustAlign v1.0 提取代码 | — | 比对/索引 | bwa 子进程 | MIT |
| `rayon` | 1.8+ | 数据并行 | multiprocessing.Pool | MIT/Apache-2.0 |
| `memmap2` | 0.9 | 大文件 mmap | — | MIT/Apache-2.0 |
| `serde` + `serde_json` | 1.x | 参数/运行元数据持久化 | — | MIT/Apache-2.0 |
| `csv` | 1.x | BED/表格输出 | print 手写 | MIT/Unlicense |
| `indicatif` | 0.17 | CLI 进度条 | print 计数 | MIT |
| `tracing` + `tracing-subscriber` | 0.1+ | 结构化日志 | print | MIT |
| `flate2` | 1.0 | gz IO（rustkmer 已带） | — | MIT/Apache-2.0 |
| `tempfile` | 3.x | 临时文件 | os.path 手拼 | MIT/Apache-2.0 |

### 4.2 P2/P3（GUI）

| crate | 用途 | 备注 |
|---|---|---|
| `eframe`/`egui`（0.34+，wgpu 默认后端） | GUI 框架 | 中文 UI 需内嵌 CJK 字体（如 思源黑体 子集） |
| `egui_plot` | 染色体探针密度图、框选缩放 | 0.37，官方作者维护 |
| `rfd` | 原生文件/目录选择对话框 | 三平台原生（Windows IFileDialog / macOS NSPanel / GTK） |
| `crossbeam-channel` | 管线线程 → GUI 进度/日志 | 替代 Qt 信号 |
| `parking_lot` | 轻量锁 | rustkmer 已依赖 |

### 4.3 P4/P5（后期，先列出备齐）

| crate | 用途 | 替代的 Python 依赖 |
|---|---|---|
| `bigtools` | BigWig 读写（有论文、性能好） | pyBigWig |
| `bedrs` | 区间 intersect/sort | pybedtools |
| `noodles`（bam/bcf/vcf/faidx） | BAM 深度、VCF/BCF 读取、fai | pysam、samtools/bcftools 子进程 |
| （未来）`primer3` crate | Tm/hairpin/dTm 恢复 | primer3-py；GPL-2.0，引入前需决策许可问题 |
| （备选）`polars` | 仅当表格操作复杂时 | pandas（quantile 等小需求手写即可，避免重依赖） |

> 许可证结论：当前清单全部 MIT/Apache-2.0，与 Chorus2 的 MIT 完全兼容。唯一的 GPL 风险点 primer3 已按约定移出首版；未来恢复时要么接受整体转 GPL，要么改用自实现 SantaLucia 最近邻热力学（工作量可控，另立决策）。

---

## 5. 分阶段迁移计划

```
P0  Workspace 骨架 + CI（三平台矩阵：build + test + clippy + fmt）
     └─ 交付：空管线可跑通 cargo test

P1  ★ 基本设计功能（第一步，纯 CLI）
     1. chorus-io: FASTA 索引读写、BED/len 输出
     2. chorus-kmer: rustkmer 封装 + 滑窗打分（复刻计数裁剪与 off-by-one 语义）
     3. chorus-align: RustAlign v1.0 提取 + 参数映射(-O 0 -B 0 -k 5) + AS/XS 过滤
     4. chorus-core: 管线编排（k 公式、分块、间距选择、输出）+ 取消/进度
     5. chorus-cli: `chorus3 design -g ... -i ...` 子命令
     6. ★ 验收门：Golden test —— 用 TAIR10/水稻等小基因组，
        Chorus2 以 --skipdtm 运行的输出 vs Chorus3 输出，探针集合一致
     └─ 交付：三平台 `chorus3` 二进制

P2  GUI 向导（egui）：三页表单 + 运行监控 + 取消
     └─ 验收：从 GUI 完成一次完整设计并得到与 CLI 相同的输出

P3  探针浏览器 GUI：染色体密度图 + 框选编辑 + 项目保存

P4  ChorusNGSfilter/select：bigtools 写 BigWig + 分位数筛选

P5  ChorusNoRef / ChorusHomo：noodles 读 BAM/VCF、consensus 调用、多物种比较

P6  打包分发（见 §6）
```

每阶段独立可发布；P1 完成即"基本设计功能"达成。

---

## 6. 跨平台策略（Windows / macOS / Linux）

### 6.1 依赖面

- **全纯 Rust 依赖**（rustkmer、RustAlign、egui 均无 C 库），三平台只需 Rust 工具链，`cargo build --release` 即可；这是相比旧版（要求用户自装 jellyfish/bwa，Windows 还要 MSYS 环境）最大的体验改善。
- egui/eframe(wgpu)：Windows 用 D3D12/运行时自带；macOS 用 Metal；Linux 用 Vulkan/GL（主流发行版开箱即用）。Linux 的 Wayland/X11 均支持。
- 无 C 依赖也意味着交叉编译友好，CI 矩阵简单。

### 6.2 各平台要点

| 平台 | 工具链/目标 | 打包 | 注意 |
|---|---|---|---|
| Windows | `x86_64-pc-windows-msvc` | MSI + 便携 zip（cargo-packager） | 路径含空格/中文是旧版重灾区，新版全内部处理；长路径注册表开关不依赖 |
| macOS | `aarch64-apple-darwin` + `x86_64-apple-darwin` → universal2 | DMG（cargo-packager） | 分发需 Developer ID 签名 + 公证（可选，教育用户可免） |
| Linux | `x86_64-unknown-linux-gnu`（+ 可选 aarch64） | AppImage / deb / tarball | 在最老支持基线（如 ubuntu-22.04 容器）编译控制 glibc；AppImage 免装 |

### 6.3 CI/发布

- GitHub Actions 矩阵：`windows-latest / macos-latest(arm) / ubuntu-latest`，跑 `cargo fmt --check && cargo clippy -- -D warnings && cargo test`。
- 发布：`cargo-dist`（CLI 的 tarball+脚本安装器）+ `cargo-packager`（GUI 的 MSI/DMG/AppImage），tag 触发自动出三平台产物。

---

## 7. 迁移时必须统一的实现细节（Chorus2 中的不一致/坑）

这些是逐行阅读发现的差异，Chorus3 需明确选择一种语义并写进测试：

1. **step 默认值**：CLI 默认 5，GUI 硬编码 4。
2. **maxkmerscore 公式**：CLI = `int(((len×h/100) - k) × ploidy/2 + 0.5)`；GUI = `int(len×h/100) - k`。
3. **滑窗求和的 off-by-one**：`sum(count[st : st+len-1])` 故意不含最后一个 k-mer 起点——复刻时建议保留行为但注释说明（或经 golden test 验证后修正）。
4. **k-mer 计数裁剪**：query 原始计数 → `2→1, ≥3→2, 其他→0`。
5. **`getch()`（tty/termios）只工作在 Unix**——旧版 Windows CLI 的 Y/N 确认根本不可用；Chorus3 用 `std::io::stdin` 常规读取或 `--force` 参数。
6. **bwa 参数**：`-O 0 -B 0 -k 5`，`-E` 沿用默认；AS 阈值 = length，XS 阈值 = `length×homology/100`（`int()` 截断）。
7. **间距选择窗口** = `max(len(rprimer),5) + probe_length`；首条探针从 0 起算。
8. `bwa.bwareflength` 靠发一条假读段取 `@SQ`——RustAlign 直接给参考长度，此 hack 消失。

---

## 8. 风险与开放问题

| # | 风险 | 影响 | 缓解 |
|---|---|---|---|
| R-1 | rustkmer 全内存精确计数，>10 Gbp 基因组（小麦 17 Gbp）内存可能超限 | 大基因组用户 | ① 分染色体 count + streaming_merge（已有）；② 给 rustkmer 加 Bloom Filter 前置（jellyfish `--bf-size` 等价）；③ 文档标注内存需求。P1 验收用拟南芥/水稻/玉米（≤2.5 Gbp）。Windows 上 count/query 路径已实测可用（§2.3） |
| R-2 | RustAlign 代码从 v1.0 tag 提取，含裸指针需清理；v2 重构未完成 | P1 工作量 | 只提取 `mem/index/bntseq/ksw` 所需模块；641 个测试随迁作为安全网。v1.0 已在 Windows 实测端到端可用（§2.3），提取时附带修复：FASTA 输入、索引前缀规则、裸指针输入校验 |
| R-3 | golden test 对齐：bwa 与 RustAlign 在长重复区段的 MAPQ/次优命中可能有小差异 → XS 过滤边缘差异 | 输出集合差异 | 验收标准定为"探针集合 ≥99.9% 一致 + 差异逐条解释" |
| R-4 | egui_plot 维护者空缺 | P3 | 需求已覆盖；必要时 fork |
| R-5 | 中文 UI 字体内嵌（体积 +2~5 MB） | 包体积 | 字体子集化 |
| R-6 | macOS 签名/公证 | 分发体验 | 首版可免签名发布（用户右键打开）；后续补 |

## 9. 待你确认的决策点

1. **GUI 框架**：确认 egui（我推荐），还是想看 Slint/Tauri 的详细评估？
2. **rustkmer/RustAlign 接入方式**：git submodule / vendor 进 chorus3 仓库 / crates.io 发布后作普通依赖？（建议 vendor 或 submodule，保证可改）
3. **§7 中两处公式不一致**（step、maxkmerscore）：以 CLI 版为准？
4. 仓库位置：Chorus3 在当前 Chorus2 仓库的 `Chorus3` 分支上新建 `chorus3/` 目录，还是独立新仓库？

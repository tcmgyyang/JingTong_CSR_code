# 单细胞 + 双虚拟（虚拟敲除 + 虚拟过表达）— 服务器运行指引

把 18 个多源 hub 基因**定位到颈椎/椎间盘终板的具体细胞类型**，并用 **CellOracle** 做
**虚拟敲除(KO) + 虚拟过表达(OE)**，在单细胞分辨率上验证 hub 作为治疗靶点的作用方向。
对标 OMICS 的单细胞+扰动范式。

## 数据：GSE160756（人椎间盘单细胞图谱）
**GSE160756** — 91,295 个 IVD 细胞（3 髓核 NP + 2 软骨终板 CEP + 2 纤维环 AF，10X），
**含软骨终板(CEP)，与你 bulk 的 GSE153761（颈椎退变终板）组织相符**。
- ⚠️ **节段局限（务必在论文写明）**：经穷尽检索，**无颈椎椎间盘/终板 scRNA 公开数据**；
  GSE160756 节段未标注、按惯例多为腰椎/尸体。软骨终板的细胞类型生物学在颈/腰高度共通，
  此处仅作**"hub 基因的细胞来源"组织参考**，节段差异列为局限——这比"换疾病数据集"温和得多。
- 文件为**每样本一个 `.loom`**（`GSM..._umi_h{NP,CEP,AF}_n.loom`）。
- 我已本地下载 `GSE160756_RAW.tar`(433M) 并解压到本目录 `GSE160756/`（.loom 已 gunzip）。
  你直接把 `GSE160756/` 连同脚本上传服务器即可。

## 环境（建议 conda，服务器/内存≥32G）
```bash
conda create -n oracle python=3.10 -y && conda activate oracle
pip install scanpy leidenalg igraph
pip install celloracle            # 依赖较多；如失败见 CellOracle 官方 conda 配方
```

## 运行（两步）
```bash
# 1) 注释：读 GSE160756/ 下的 7 个 .loom（自动加 NP/CEP/AF 标签）→QC→聚类→marker注释→UMAP/hub图
py 01_annotate.py GSE160756/ .       # 目录内含 *.loom；脚本也支持 10x mtx 目录 / .h5ad

# 2) CellOracle GRN + 双虚拟（KO + OE）
py 02_celloracle_dualvirtual.py adata_annotated.h5ad .
```
> .loom 需为解压状态（非 .loom.gz）。01_annotate.py 会按文件名 hNP/hCEP/hAF 自动标注 tissue。

## 关键生物学说明（脚本已正确处理）
CellOracle 扰动的是**调控因子(TF)**。18 hub 里 **JUN / MYC / ESR1 / HIF1A** 是 TF，可做
KO/OE；**MMP9 / IL1B / IL6 / TNF / CCL2 / CXCL8** 等是**效应基因**（非 GRN 调控层），
用**细胞定位**（step 1 的 dotplot/UMAP）呈现，不强行虚拟扰动。脚本自动取
`18hub ∩ oracle调控基因` 做双虚拟，名单写入 `perturbable_hubs.csv`。

## 产出
- `adata_annotated.h5ad`、`umap_clusters.png`、`umap_hub_expr.png`、`dotplot_hub_dotplot.png`
- `perturb_figs/<gene>_KO.png`、`<gene>_OE.png`（UMAP 上的扰动向量场）
- `dualvirtual_summary.csv`（每个可扰动 hub 的 KO/OE 总位移；位移越大=对细胞命运影响越大）
- `GRN_degree*`（GRN 中心性）

## 写进论文
- hub 基因的**细胞来源**：预期炎症基因(IL1B/IL6/TNF/CXCL8)→免疫/巨噬+退变软骨细胞；
  MMP9/MMP→肥大-成骨样软骨细胞。
- **双虚拟**：虚拟敲除 TF hub（如 JUN/HIF1A）若把细胞推离"退变态"，支持其为干预靶点；
  虚拟过表达作反向佐证。务必措辞为**计算预测**（in silico），非湿实验。

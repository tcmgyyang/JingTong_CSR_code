# 分子对接重做 + 扩充方案（含结构纠错）

## ⚠️ 0. 为什么要重做：结构纠错 QC
逐一核对原对接所用 PDB 的真实身份（序列比对 + RCSB 核实）：

| 靶点 | 原用 PDB | **真实身份** | 判定 |
|---|---|---|---|
| AKT1 | 8R5K | **FKBP5/FKBP51**（FK1 域，肽基脯氨酰异构酶） | ❌ 错，重做 |
| IL6 | 5GW9 | **G-CSF**（粒细胞集落刺激因子，序列已比对吻合） | ❌ 错，重做 |
| IL1B | 8RYS | Interleukin-1β（1.16 Å） | ✅ 正确，保留 |

→ **AKT1、IL6 的对接与基于同结构的 MD 全部作废**，须用下列正确结构重做；IL1B 保留。

## 1. 验证过的正确/新增结构（均已逐一开 RCSB 核实）
| 靶点 | 用途 | PDB | 说明 | 对接口袋 | 重对接对照配体 |
|---|---|---|---|---|---|
| AKT1 | 重做 | **4EKL** | 人 AKT1 激酶域, 2.0 Å, ATP 竞争抑制剂 GDC-0068 | ATP 口袋 | 0RF (GDC-0068) |
| IL6 | 重做 | **1ALU** | 人 IL-6, 1.9 Å, 四螺旋束 | 盲对接(Site II/III) | 无 |
| IL1B | 保留 | 8RYS | 人 IL-1β, 1.16 Å | (已完成) | — |
| MMP9 | 新增 | **1GKC** | 人 MMP-9 催化域, 2.3 Å, 催化锌+羟肟酸 NFH | 催化锌/S1' | NFH |
| TP53 | 新增(可选) | **2OCJ** | 人 p53 核心域(96–289), 2.05 Å | 盲对接 | 无 |

> **MMP9** 是新增首选：它既是 11 个 hub 基因之一，又是 ML ROC=0.74 的血液诊断基因，且为经典黄酮可成药酶——把"对接 ↔ ML"闭环。
> **TP53** 是 ML 最佳(ROC=0.78)但**对接困难**（野生型核心域无深口袋，属转录因子/PPI 界面）。建议：要么用 Y220C 可成药突变体口袋(如 2VUK)，要么按盲对接报告并注明为探索性；优先把 MMP9 做扎实。
> **IL6** 也是细胞因子 PPI 界面、无经典小分子口袋，盲对接结果偏软，需在文中如实说明。

## 2. 对接盒（已用共结晶配体算好，见 config_*.txt）
| 靶点 | center (x,y,z) | size | 模式 |
|---|---|---|---|
| AKT1(4EKL) | 28.03, 5.22, 10.89 | 22³ | ATP 口袋(配体 0RF) |
| MMP9(1GKC) | 65.84, 30.71, 117.75 | 22³ | 催化位点(配体 NFH) |
| IL6(1ALU) | 2.67, −19.98, 9.04 | 30³ | 盲对接 |
| TP53(2OCJ) | 2.65, −0.06, 31.41 | 30³ | 盲对接 |

## 3. 操作流程（用你 docking/ 里的 AutoDock Vina）
**依赖**：MGLTools(AutoDockTools) 或 Open Babel（受体/配体转 PDBQT）。本机未装 obabel，请在服务器装：`conda install -c conda-forge openbabel mgltools`。

### 3.1 受体准备（每个 PDB）
```bash
# 删水、删无关配体；MMP9 必须保留催化锌 ZN！加氢、加 Gasteiger 电荷 -> pdbqt
# MGLTools 路线：
pythonsh prepare_receptor4.py -r 4EKL.pdb -o AKT1_receptor.pdbqt -A hydrogens -U waters
pythonsh prepare_receptor4.py -r 1GKC.pdb -o MMP9_receptor.pdbqt -A hydrogens -U waters   # 之后手动把 ZN 行加回 pdbqt
pythonsh prepare_receptor4.py -r 1ALU.pdb -o IL6_receptor.pdbqt  -A hydrogens -U waters
pythonsh prepare_receptor4.py -r 2OCJ.pdb -o TP53_receptor.pdbqt -A hydrogens -U waters
```
### 3.2 配体准备（你已有 quercetin.sdf / kaempferol.sdf 在 分子对接/）
```bash
obabel quercetin.sdf  -O quercetin.pdbqt  -p 7.4 --partialcharge gasteiger
obabel kaempferol.sdf -O kaempferol.pdbqt -p 7.4 --partialcharge gasteiger
```
### 3.3 运行对接 + 阳性对照
```bash
bash run_docking.sh        # 见同目录脚本：4 靶点 × 2 配体；并对 AKT1/MMP9 做重对接 RMSD<2Å 对照
```

## 4. 对接完成后
- 取每个复合物结合能最低构象 → 更新 heat.csv / 结合能热图（沿用你 `分子对接结果/对接热能代码.R`）。
- 关键结合残基/氢键：把对接复合物存成 pdb（蛋白+最佳配体，配体改名 UNK），跑
  `py 分子对接结果/binding_residues.py`（已写好）自动生成残基表。
- 正确的对接构象同时作为**分子动力学起始结构**（替换之前 AKT1/IL6 的错误复合物）。

## 5. 结合残基表现状
`分子对接结果/binding_residues_table.csv` 中**仅 IL1B 行有效**（AKT1/IL6 两行来自错误蛋白，已作废）。
AKT1/IL6/MMP9/TP53 的残基表在重做对接后用上面脚本一键生成。

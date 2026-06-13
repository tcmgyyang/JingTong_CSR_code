# 药物虚拟筛选（聚焦版）— 成分×靶点 对接矩阵 服务器运行指引

把现有的"挑几个对接"升级为**系统对接矩阵**：**43 个颈痛颗粒活性成分 × 5 个 hub 受体
(AKT1/IL1B/IL6/MMP9/TP53)**，复用 `../../docking_redo` 里已纠错的结构与盒子，得到
亲和力矩阵 + 热图，系统排出最强的"成分-靶点"对（聚焦本方剂机制，不做外部老药重定位）。

## 前置（在 docking_redo 里先备好受体）
受体 pdbqt 需先准备好：`docking_redo/<HUB>_receptor.pdbqt`
（按 `docking_redo/run_docking.sh` 里的 prepare_receptor 流程；结构与盒子已就绪：
AKT1=4EKL、IL1B=8RYS、IL6=1ALU、MMP9=1GKC、TP53=2OCJ）。

## 环境（服务器）
```bash
# AutoDock Vina 1.2+ (vina 命令)、Open Babel、rdkit、pandas、seaborn
conda install -c conda-forge vina openbabel rdkit pandas seaborn -y
```
Windows 上仓库内已带 `docking/vina.exe`，脚本会自动找到。

## 运行
```bash
py dock_matrix.py
```
流程：读 `../swisstarget_input.csv` 的 44 个成分 SMILES → Open Babel 生成 3D pdbqt →
对 5 个受体逐一 Vina 对接（exhaustiveness=16）→ 取最优构象亲和力 → 矩阵 + 热图。

## 产出
- `docking_affinity_matrix.csv`（44×5 亲和力，kcal/mol）
- `../figures/Fig_docking_matrix.png/.pdf`（热图，越负=结合越强=颜色越深；最强成分排上方）
- 终端打印 ≤ −8 kcal/mol 的强结合对（候选重点对接，可挑去做 MD）

## 写进论文
- 用矩阵热图替代/补充原来的单个对接图，体现**系统性筛选**。
- 报告每个 hub 的最强成分（预期黄酮类槲皮素/山奈酚对 AKT1/MMP9 结合强），
  与 PPI/ML/MR 的 hub 优先级交叉印证；选 top 对接对进入 MD（与 AKT1/IL6/IL1B 的 MD 衔接）。
- 注意 TP53 是 ML/MR 标志但非≥2源药物 hub，矩阵里保留作参考行。

## 备注
- 配体 pdbqt 默认用 Open Babel `--gen3d`；如需更规范可改用 **meeko**（Vina 官方配体准备）。
- 矩阵=220 次对接×exhaustiveness16，建议服务器多核；可调低 exhaustiveness 先粗筛。

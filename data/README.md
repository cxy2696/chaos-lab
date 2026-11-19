# 数据目录

本目录的数据全部由 [`scripts/generate_data.sh`](../scripts/generate_data.sh) 从本仓库的 `chaos_cli` 生成，不是手工编辑的。参数、行数、SHA256 与对应命令见 [`docs/data-provenance.md`](../docs/data-provenance.md)。

## 文件

| 文件 | 行数 | 内容 |
|---|---:|---|
| `lorenz_trajectory.csv` | 6000 | Lorenz $t,x,y,z$，RK4 $h=0.01$，丢弃前 5 个时间单位 |
| `lorenz_series.csv` | 8000 | Lorenz 同参数，丢弃前 50 个时间单位，用于嵌入/GP/Rosenstein |
| `rossler_trajectory.csv` | 8000 | Rössler $a=b=0.2,c=5.7$ |
| `chua_trajectory.csv` | 10000 | Chua 双涡卷 |
| `three_body_figure_eight.csv` | 7000 | 8 字三体， $h=0.001$ |
| `double_pendulum.csv` | 5000 | 双摆角度/角速度 |
| `henon_orbit.csv` | 5000 | Hénon 映射轨道 |
| `logistic_orbit.csv` | 2000 | Logistic $r=4$ 轨道 |
| `logistic_bifurcation.csv` | 120000 | $r\in[2.5,4]$，每个参数保留 80 个后期迭代 |
| `logistic_lyapunov.csv` | 1500 | $\lambda(r)$ 扫描 |
| `lorenz_poincare.csv` | 240 | $z=27$ 下穿截面 |
| `lorenz_gp.csv` | 60 | 关联积分与局部斜率 |
| `lorenz_rosenstein.csv` | 301 | 平均发散曲线 |
| `lorenz_spectrum.csv` | 4097 | 单边功率谱 |
| `lorenz_rqa_matrix.csv` | 794×794 | 递归矩阵（无表头） |
| `lorenz_zeroone.csv` | 64 | $K_c$ 样本 |
| `lorenz_sindy.csv` | 10 | SINDy 系数表 |
| `lorenz_sync_error.csv` | 8001 | 同步误差曲线 |
| `feigenbaum.csv` | 8 | 超稳定周期参数与 $\delta$ 估计 |

## 校验与再生

```bash
cd ..
./scripts/generate_data.sh          # 重新生成 CSV
python3 scripts/make_metadata.py   # 重新计算 SHA256 与 metrics.json
shasum -a 256 -c data/checksums.sha256   # macOS
# 或：sha256sum -c data/checksums.sha256
```

> `make_metadata.py` 与 `plot_figures.py` 只依赖 Python 标准库 / Pillow；仓库不把 matplotlib 作为运行依赖。

## 注意

- 这些是**数值实验数据**，不是观测数据；不要把它们当作物理实验的原始测量。
- 混沌系统的长期轨迹对步长/舍入误差敏感；数据用于图示与算法验证，不用于高精度定量结论。
- 若修改了生成参数，请同步更新 `docs/data-provenance.md` 和 `data/metrics.json`，并重新生成校验和。

# 数据溯源与数值验证

本页回答三个问题：

1. 仓库里的数据是怎么来的？
2. 文档中的数值断言对应哪条命令、哪个文件？
3. 哪些是文献值，哪些是本仓库实测值？

理论文档是 [`docs/chaos-systems.md`](chaos-systems.md)；算法与 API 是 [`docs/algorithms.md`](algorithms.md)。

---

## 1. 数据生成链

```
docs/chaos-systems.md（理论、公式、引用 [n]）
        ├─ include/chaos/*.hpp + src/*.cpp（算法）
        ├─ apps/chaos_cli.cpp（可复现命令）
        └─ tests/*（数值断言）
             └─ data/*.csv（由 CLI 生成）
                  └─ docs/figures/*.png（由 scripts/plot_figures.py 读取数据）
```

生成脚本：[`scripts/generate_data.sh`](../scripts/generate_data.sh)。它使用的二进制是当前 `./build/chaos_cli`；编译参数与 `CMakeLists.txt` 一致（Release，`-O3` 由 CMake 默认）。

环境（本次生成）：

| 项 | 值 |
|---|---|
| 机器 | macOS arm64 (Apple Silicon) |
| 编译器 | Apple clang 21.0.0 |
| CMake | 4.5.0 |
| 构建类型 | Release |
| 库版本 | 0.1.0 |

其他平台可能得到末位不同的浮点结果；测试容差按物理/数值意义设定，不要求逐位一致。

---

## 2. 数据文件与命令

| 文件 | 行数 | 生成命令（省略 `./build/chaos_cli`） |
|---|---:|---|
| `lorenz_trajectory.csv` | 6000 | `simulate --system lorenz --n 6000 --h 0.01 --transient 5` |
| `lorenz_series.csv` | 8000 | `simulate --system lorenz --n 8000 --h 0.01 --transient 50` |
| `rossler_trajectory.csv` | 8000 | `simulate --system rossler --n 8000 --h 0.01 --transient 20` |
| `chua_trajectory.csv` | 10000 | `simulate --system chua --n 10000 --h 0.01 --transient 20` |
| `three_body_figure_eight.csv` | 7000 | `simulate --system three_body --n 7000 --h 0.001` |
| `double_pendulum.csv` | 5000 | `simulate --system double_pendulum --n 5000 --h 0.005` |
| `henon_orbit.csv` | 5000 | `simulate --map henon --n 5000 --transient 1000` |
| `logistic_orbit.csv` | 2000 | `simulate --map logistic --n 2000 --transient 1000` |
| `logistic_bifurcation.csv` | 120000 | `bifurcation --map logistic --pmin 2.5 --pmax 4 --nparams 1500 --transient 1000 --keep 80` |
| `logistic_lyapunov.csv` | 1500 | `lyapunov-scan --map logistic --pmin 2.5 --pmax 4 --nparams 1500 --transient 1000 --n 20000` |
| `lorenz_poincare.csv` | 240 | `poincare --system lorenz --coord 2 --value 27 --direction -1 --h 0.005 --transient 20 --t 200` |
| `lorenz_gp.csv` | 60 | `dimension --input data/lorenz_series.csv --column 1 --m 3 --tau 17 --theiler 50 --nr 60` |
| `lorenz_rosenstein.csv` | 301 | `rosenstein --input data/lorenz_series.csv --column 1 --m 3 --tau 17 --theiler 100 --max-k 300 --dt 0.01` |
| `lorenz_spectrum.csv` | 4097 | `spectrum --input data/lorenz_series.csv --column 1 --dt 0.01` |
| `lorenz_rqa_matrix.csv` | 794×794 | `rqa --input data/lorenz_series.csv --column 1 --stride 10 --m 3 --tau 3 --theiler 5 --rate 0.03 --matrix-out ...` |
| `lorenz_zeroone.csv` | 64 | `zeroone --input data/lorenz_series.csv --column 1 --stride 10` |
| `lorenz_sindy.csv` | 10 | `sindy --system lorenz --n 5000 --transient 5000 --h 0.002 --degree 2 --threshold 0.1` |
| `lorenz_sync_error.csv` | 8001 | `sync --t 40 --h 0.005` |
| `feigenbaum.csv` | 8 | `scripts/feigenbaum.py` |

校验和：[`data/checksums.sha256`](../data/checksums.sha256)。`data/metrics.json` 是上述实测值的机器可读汇总。

---

## 3. “断言 → 来源 → 复现”表

### 3.1 文献值与本仓库实测值的区别

| 断言 | 文献值 | 本仓库实测 | 类型 | 复现/测试 |
|---|---|---|---|---|
| Logistic $r=4$ 的 $\lambda$ | $\ln2=0.693147$ | $0.693149$ | 数值 | `test_core.cpp::logistic_lyapunov_known_values`； $10^6$ 次迭代 |
| Feigenbaum $\delta$ | $4.669201609$ | $4.669157$ | 数值 | `scripts/feigenbaum.py`；超稳定周期 $2^k$ |
| Logistic 倍周期累积点 | $r_\infty=3.569945672$ | 由 $\lambda$ 变号与周期窗口定位 | 数值 | `logistic_lyapunov.csv`、`logistic_bifurcation.csv` |
| Lorenz $\lambda$ 谱 | $(0.906,0,-14.57)$ [9,17] | $(0.9101,-0.0043,-14.5724)$ | 数值 | `lyapunov --system lorenz --h 0.01 --transient 50 --duration 500` |
| Lorenz 维数 | $\approx2.06$ [9,24] | $D_2=2.047$； $D_{KY}=2.062$ | 数值 | `dimension ... --m 3 --tau 17`； $D_{KY}$ 由谱计算 |
| Hénon $\lambda$ | $(0.419,-1.623)$ [11,17] | $(0.4206,-1.6246)$ | 数值 | `lyapunov --map henon --n 100000` |
| Hénon $D_2$ | $\approx1.21$ [24] | 自动 $1.156$；手动尺度 $1.211$ | 数值 | `test_analysis.cpp`（两种设置） |
| Rössler $\lambda_1$ | $\approx0.071$（ $a=b=0.2,c=5.7$） | $0.0790$ | 数值 | `lyapunov --system rossler` |
| Chua 双涡卷 | 正 $\lambda_1$ [13] | $0.4247$ | 数值 | `lyapunov --system chua` |
| 三体 8 字周期 | $T=6.32591398$ [14] | 一个周期后回到初值 $2.0\times10^{-4}$ | 数值 | `test_core.cpp::three_body_figure_eight_periodic` |
| 0-1 检验 Lorenz | $K\approx1$ [28] | 降采样后 $0.9905$；原采样 $-0.274$ | 数值/采样陷阱 | `zeroone ... --stride 10`；见 §5 |
| 同步 | 误差指数收敛 [33] | $t\approx7.6$ 后 $<10^{-6}$，最终机器精度 | 数值 | `sync --t 40 --h 0.005` |

### 3.2 文档中的历史/理论断言

| 断言 | 出处 |
|---|---|
| 三体问题中复杂轨道（1890） | [1] |
| Lorenz 1963 对流模型与蝴蝶效应 | [2] |
| 集合预报/可预报性 | [3,40,41] |
| 准周期路径（Ruelle–Takens–Newhouse） | [4,47] |
| 周期三意味着混沌 | [5] |
| Logistic 映射的倍周期图景 | [6] |
| Feigenbaum 普适常数 | [7,8] |
| Lorenz 系统的谱/维数/参数图景 | [9] |
| Hénon 吸引子 | [11] |
| Rössler 方程 | [12] |
| Chua 双涡卷 | [13] |
| 8 字三体 | [14] |
| Takens 嵌入 | [20,21] |
| AMI / FNN / GP | [22,23,24] |
| Wolf / Rosenstein / Benettin | [16,17,19] |
| RQA | [25,26] |
| 替代数据 | [27] |
| 0-1 检验 | [28] |
| Kaplan–Yorke 维数 | [29] |
| OGY / Pyragas / Pecora–Carroll | [31,32,33] |
| SINDy / RC / DMD | [36,37,38] |
| 混沌图像加密的风险 | [45] |
| 渗流/重整化群 | [46,50] |

---

## 4. 数值协议

### 4.1 积分

- Lorenz/Rössler/Chua：RK4， $h=0.01$；瞬态 50（Lorenz/Rössler/Chua 谱）。
- 三体：RK4， $h=0.002$；能量漂移 $<10^{-12}$。
- 双摆：RK4， $h=0.001$；20 秒能量漂移 $<10^{-9}$。
- 辛积分器测试：Verlet/Yoshida4， $h=0.01$，10000 步；Yoshida 能量漂移 $<10^{-9}$。

### 4.2 Lyapunov

- Benettin：解析 Jacobian；每步重正交化；瞬态 30–50，统计窗口 120–500。
- Hénon：1000 次瞬态，100000 次迭代；每步 QR。
- Rosenstein：Lorenz $m=3,\tau=17$、Theiler=100、`max_k=300`、 $dt=0.01$；自动拟合窗口 $k\in[148,297]$。
- Hénon Rosenstein： $m=2,\tau=1$、Theiler=10、手动拟合 $k=0..10$。

### 4.3 嵌入与 GP

- Lorenz： $m=3$， $\tau=17$（AMI 第一极小），Theiler=50，60 个半径，拟合分位 0.2–0.8。
- Hénon： $m=2,\tau=1$，Theiler=10；自动范围给出 1.156，手动尺度给出 1.211。
- GP 自动半径：取距离分布的 0.005%–0.3% 分位（至少 100 对），这只是**起点**；必须检查局部斜率平台。

### 4.4 0-1 检验

- Logistic：4000 点， $c$ 在 $(\pi/5,4\pi/5)$ 随机，64 个 $c$，取中位数。
- Lorenz：原始 $dt=0.01$ 给 $K\approx-0.274$；`--stride 10`（ $dt=0.1$）给 $K\approx0.9905$。这不是库的“不稳定”，而是 0-1 检验对过采样的已知敏感性；见理论文档 §5.5。
- 正弦： $K\approx-0.002$。

### 4.5 RQA

- Lorenz：先 `--stride 10` 降采样， $m=3,\tau=3$，Theiler=5，目标 RR=3%， $\epsilon=3.20195$。
- 结果：RR=2.96%，DET=0.691，LAM=0.293，ENTR=2.039 nat，TT=2.044，Lmax=73。
- 注意：DET/LAM 依赖 $\epsilon$ 和 Theiler；不要跨不同设置直接比较。

---

## 5. 已知差异与解释

| 现象 | 解释 |
|---|---|
| Lorenz Benettin $0.910$ vs 文献 $0.906$ | 参数、步长、瞬态、统计窗口；差约 0.4% |
| Rosenstein $0.855$ vs Benettin $0.910$ | 有限数据、拟合区间、嵌入；差约 6% |
| Hénon GP 自动 $1.156$ vs 文献 $1.21$ | 自动半径偏小；手动尺度区可到 1.211 |
| 0-1 原采样 $K\approx0$ | 过采样；降采样后 $K\approx0.9905$ |
| 三体周期返回 $2\times10^{-4}$ | 周期常数截断 + 混沌敏感性；不是积分失败 |
| Chua 谱和 $-3.80$ | 分段线性区域的平均散度；不同参数/区域会变 |

---

## 6. 校验与再生

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
./scripts/generate_data.sh
/Users/xinyi/.dsh/.../python3 scripts/make_metadata.py   # 或任意带 Pillow 的 Python
shasum -a 256 -c data/checksums.sha256
```

`data/metrics.json` 中的每个值都能在本文第 3 节找到命令；如果修改了算法或参数，请同时更新测试容差、metrics 和本文。

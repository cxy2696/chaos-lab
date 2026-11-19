# 实验手册

每个实验都给出：目标 → 命令 → 期望输出 → 需要观察/思考的点。所有命令在仓库根目录执行；先构建：

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

---

## 实验 0：先确认数值与测试基线

```bash
./build/chaos_cli version
./build/chaos_cli models
./build/chaos_tests
```

期望：`chaos_tests` 输出 `33 passed, 0 failed`（不同平台可能有末位差异，但断言应通过）。测试覆盖：

- RK4/RK45 精度、辛积分器能量；
- 解析 Jacobian 与有限差分一致；
- Lorenz/Hénon/猫映射 Lyapunov 谱；
- Logistic 已知值；
- 三体 8 字周期、双摆能量、Chua 有界；
- AMI/FNN/GP/Rosenstein/0-1/FFT/替代数据/RQA/Poincaré/分岔；
- OGY/Pyragas/同步/SINDy；
- C API。

如果失败，先看具体断言与 `docs/data-provenance.md` 的容差说明，不要急着改算法。

---

## 实验 1：Logistic 分岔图与 Lyapunov 指数

### 目标

观察倍周期级联、周期窗口和 $\lambda$ 变号；验证 $r=4$ 时 $\lambda=\ln2$。

### 命令

```bash
./build/chaos_cli bifurcation --map logistic --pmin 2.5 --pmax 4.0 \
  --nparams 1500 --transient 1000 --keep 80 --out data/logistic_bifurcation.csv
./build/chaos_cli lyapunov-scan --map logistic --pmin 2.5 --pmax 4.0 \
  --nparams 1500 --transient 1000 --n 20000 --out data/logistic_lyapunov.csv
```

### 期望

- 分岔图从 $r=3$ 开始周期 1→2→4→8；
- $r_\infty\approx3.5699457$ 附近进入混沌；
- $r\approx3.8284$ 有周期 3 窗口；
- $\lambda(r)$ 在 $r=4$ 接近 $0.693$。

`data/metrics.json` 记录 $\lambda(4)=0.693149$；Feigenbaum $\delta$ 的估计见实验 1b。

### 1b：Feigenbaum 常数

```bash
python3 scripts/feigenbaum.py
cat data/feigenbaum.csv
```

该脚本找超稳定周期 $f^{2^k}(1/2)=1/2$ 的参数 $s_k$，用

$$
\delta_k=\frac{s_{k-1}-s_{k-2}}{s_k-s_{k-1}}
$$

估计。期望最后一个 $\delta\approx4.669157$（参考 $4.669201609$）。

### 思考

- 为什么 $\lambda$ 在周期窗口内会掉到负值？
- 分岔图里“空白带”是什么？与周期窗口/混沌带的关系？
- 改变 `--keep`、`--transient` 对图有什么影响？哪些是真实结构，哪些是采样伪影？

---

## 实验 2：Lorenz 系统与 Lyapunov 谱

### 目标

用 Benettin 算法得到完整谱，检查与散度、Kaplan–Yorke 维数的关系。

### 命令

```bash
./build/chaos_cli lyapunov --system lorenz --h 0.01 --transient 50 --duration 500
./build/chaos_cli simulate --system lorenz --n 6000 --h 0.01 --transient 5 \
  --out data/lorenz_trajectory.csv
```

### 期望

```
system lorenz exponents: 0.910064 -0.004271 -14.572357  sum=-13.666564
```

- $\lambda_1>0$，$\lambda_2\approx0$，$\lambda_3<0$；
- 和 $\approx-(\sigma+1+\beta)=-13.6667$；
- $D_{KY}=2+\lambda_1/|\lambda_3|\approx2.062$。

### 收敛性检查

```bash
for h in 0.02 0.01 0.005; do
  ./build/chaos_cli lyapunov --system lorenz --h $h --transient 50 --duration 500
done
```

观察 $\lambda_1$ 随 $h$ 的变化。若变化很大，说明步长不够；若统计窗口太短，指数会抖动。可把 `duration` 增到 2000 看是否稳定。

### 思考

- 为什么 $\lambda_2$ 接近 0？连续自治流的零指数对应什么方向？
- 数值积分误差会不会破坏“和 = 散度”？为什么？
- 双轨迹法（测试中的 `lyapunov_max_two_trajectories`）与 Benettin 的 $\lambda_1$ 是否一致？

---

## 实验 3：从一维观测重构 Lorenz 吸引子

### 目标

只用 $x(t)$，用 AMI 选 $\tau$、FNN 选 $m$、GP 估维数、Rosenstein 估 $\lambda_{\max}$。

### 命令

```bash
# data/lorenz_series.csv 已由实验 2 的 simulate 生成（或直接使用仓库数据）
./build/chaos_cli ami --input data/lorenz_series.csv --column 1 --max-tau 60
./build/chaos_cli fnn --input data/lorenz_series.csv --column 1 --tau 17 --max-m 8
./build/chaos_cli dimension --input data/lorenz_series.csv --column 1 \
  --m 3 --tau 17 --theiler 50 --nr 60 --out data/lorenz_gp.csv
./build/chaos_cli rosenstein --input data/lorenz_series.csv --column 1 \
  --m 3 --tau 17 --theiler 100 --max-k 300 --dt 0.01 --out data/lorenz_rosenstein.csv
```

### 期望

| 量 | 结果 |
|---|---|
| AMI 第一极小 | $\tau=17$ |
| 自相关 $1/e$ | 约 28 |
| FNN 建议 | $m=3$ |
| GP | $D_2=2.0466$，$R^2=0.9999$ |
| Rosenstein | $\lambda=0.855$，拟合 $k\in[148,297]$ |

### 关键观察

- 把 `--m` 改成 2 或 4，看 $D_2$ 和 FNN 比例如何变；
- 把 `--tau` 改成 5 或 30，看 GP 标度区是否还平；
- 打开 `data/lorenz_gp.csv` 的 `local_slope` 列，确认平台；不要只看一个 $D_2$ 数字；
- Rosenstein 的自动窗口对 Hénon 可能偏保守，可手动指定 `--fit-start/--fit-end`。

### 思考

- Takens 定理要求 $m\ge2d+1$；这里 $d\approx2.06$，为什么 $m=3$ 就够用（提示：泛型条件、有限数据、$D_2$ 不是状态维数）？
- 如果只用 2000 个点，GP 的标度区会怎样变化？
- 如果加 1% 噪声，AMI/FNN/GP 哪个最先失效？

---

## 实验 4：控制、同步与 SINDy

### 4a. OGY 控制

运行测试中的 OGY 用例：

```bash
./build/chaos_tests ogy
```

期望：logistic $r=3.8$、$x^*=0.736842$，200000 步内锁定，尾段误差 $<10^{-8}$。

手工改参数：把 `tolerance` 从 0.2 降到 0.05，观察激活次数与锁定时间；把 `dp_max` 降到 0.005，观察是否还能锁定。这能说明 OGY 的“局部线性化 + 执行器饱和”限制。

### 4b. Pyragas

```bash
./build/chaos_tests pyragas
```

期望：延迟 $d=1$、$K=-0.5$ 时，logistic 不动点稳定。尝试 $K=+0.5$，观察是否失稳；用线性化特征方程解释。

### 4c. Pecora–Carroll 同步

```bash
./build/chaos_cli sync --t 40 --h 0.005 --out data/lorenz_sync_error.csv
```

期望：`synchronized=true`，`settling_time≈7.58`，最终误差 0（机器精度）。把响应初值改成 `--response0 "5,-5"`，看同步时间是否变长。

### 4d. SINDy

```bash
./build/chaos_cli sindy --system lorenz --n 5000 --transient 5000 --h 0.002 \
  --degree 2 --threshold 0.1 --out data/lorenz_sindy.csv
```

期望输出：

```
d x0/dt = -10*x0 + 10*x1
d x1/dt = 28*x0 - 1*x1 - 1*x0*x2
d x2/dt = -2.66667*x2 + 1*x0*x1
```

把 `--threshold` 改成 0.01 或 1.0，看是否出现虚假项或漏项；把 `--n` 降到 500，看系数误差如何增大。

---

## 实验 5：0-1 检验与采样陷阱

### 命令

```bash
./build/chaos_cli zeroone --input data/lorenz_series.csv --column 1 --stride 1
./build/chaos_cli zeroone --input data/lorenz_series.csv --column 1 --stride 10
```

### 期望

- `--stride 1`（$dt=0.01$）：$K\approx-0.274$，**误判为规则**；
- `--stride 10`（$dt=0.1$）：$K\approx0.9905$，正确判为混沌。

### 解释

0-1 检验的 $M_c(n)$ 需要看到扩散项随 $n$ 线性增长。过采样时，强周期分量主导累积和，$n_{\max}=N/10$ 覆盖的时间尺度太短，$M_c$ 看起来有界。**先按相关时间降采样，再做检验**。

### 思考

- 对 Logistic 映射，`--stride` 有意义吗？（映射本身是离散的，$dt=1$。）
- 对 Rössler，合适的 stride 是多少？用自相关/AMI 估计相关时间。
- 加噪声后 $K$ 如何变化？噪声太强时 0-1 检验会失效吗？

---

## 实验 6：RQA 与递归图

### 命令

```bash
./build/chaos_cli rqa --input data/lorenz_series.csv --column 1 --stride 10 \
  --m 3 --tau 3 --theiler 5 --rate 0.03 --matrix-out data/lorenz_rqa_matrix.csv
```

### 期望

```
epsilon=3.20195 RR=0.0296232 DET=0.691293 LAM=0.293427 ENTR=2.03944 TT=2.04445 Lmax=73
```

### 观察

- 改变 `--rate`（0.01、0.05、0.1），看 DET/LAM 如何变化；
- 改变 `--theiler`，看对角线结构是否被时间近邻污染；
- 对周期正弦和噪声各跑一次，比较 DET（规则信号高、噪声低）。

### 思考

- 为什么 RR 要固定？跨数据集比较 DET 时，$\epsilon$ 不同会带来什么偏差？
- LAM/TT 对层流/阵发性信号有什么意义？
- 递归图的对角线长度与 Lyapunov 时间有什么关系？

---

## 实验 7（扩展）：自己加一个系统

1. 在 `include/chaos/systems.hpp` 声明 struct，实现 `f` 与 `jacobian`；
2. 在 `src/systems.cpp` 实现，并加入 `ode_catalog()` 与 `default_initial_state()`；
3. 在 `tests/test_core.cpp` 加 Jacobian 对照 + 一个已知不变量/周期；
4. 用 `chaos_cli simulate/lyapunov` 跑通；
5. 在 `docs/chaos-systems.md` 的模型表加一行，注明参数与出处。

建议从 **Duffing 振子**、**标准映射**、**Kuramoto–Sivashinsky（PDE，需离散化）** 中选一个。加系统时优先写解析 Jacobian，否则 Lyapunov 谱会慢很多。

---

## 附：实验记录模板

```text
日期/提交：
命令：
参数：
观测（数值/图）：
与文献/上次的差异：
可能原因：
下一步：
```

把实验记录放在自己的分支/文档里，不要直接改 `data/` 的生成结果；需要新数据时改 `scripts/generate_data.sh` 并重新生成校验和。

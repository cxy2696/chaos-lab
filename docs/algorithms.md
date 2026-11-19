# 算法库 API 与数值约定

本文是代码的使用说明，不重复理论推导；理论见 [`chaos-systems.md`](chaos-systems.md)。所有头文件位于 [`include/chaos/`](../include/chaos/)。

## 1. 构建与安装

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
./build/chaos_cli models
```

只构建库：

```bash
cmake -S . -B build -DCHAOS_BUILD_TESTS=OFF -DCHAOS_BUILD_CLI=OFF
cmake --build build -j
```

安装：

```bash
cmake --install build --prefix /usr/local
```

要求：C++17、CMake ≥ 3.20；无第三方运行时依赖。C 接口只需要 C99 头文件。

## 2. 核心类型

```cpp
namespace chaos {
using State = std::vector<double>;

struct OdeSystem {
    std::string name, description;
    int dim;
    std::vector<std::string> param_names;
    std::vector<double> params;
    std::function<State(const State&, double)> rhs;
    std::function<Matrix(const State&, double)> jacobian; // 可选
};

struct MapSystem {
    std::string name, description;
    int dim;
    std::vector<std::string> param_names, params;
    std::function<State(const State&)> step;
    std::function<Matrix(const State&)> jacobian; // 可选
};
}
```

`Matrix` 是行主序的 `rows × cols` 稠密矩阵，仅用于小维度（三维 Jacobian、SINDy 正规方程）。不要拿它做大规模线性代数。

## 3. 模型与工厂

每个模型是一个轻量 struct，既能直接调用方程，也能转成 `OdeSystem`/`MapSystem`：

```cpp
Lorenz lorenz{10, 28, 8.0/3};
State dx = lorenz.f({1,1,1});
Matrix J = lorenz.jacobian({1,1,1});
OdeSystem sys = lorenz.ode();
```

工厂用于 CLI/运行时选择：

```cpp
OdeSystem sys = make_ode_system("lorenz");                       // 默认参数
OdeSystem sys2 = make_ode_system("lorenz", {10, 28, 2.7});       // 覆盖
MapSystem henon = make_map_system("henon", {1.4, 0.3});
State x0 = default_initial_state("lorenz");
```

可用名称：

| 类型 | 名称 | 参数（顺序） |
|---|---|---|
| ODE | `lorenz` | `sigma,rho,beta` |
| ODE | `rossler` | `a,b,c` |
| ODE | `chua` | `alpha,beta,m0,m1` |
| ODE | `double_pendulum` | `m1,m2,l1,l2,g` |
| ODE | `three_body` | `m1,m2,m3,G` |
| 映射 | `logistic` | `r` |
| 映射 | `henon` | `a,b` |
| 映射 | `cat` | — |

`make_ode_system` 只接受**完整**参数列表；空列表表示默认值。参数错误会抛 `std::invalid_argument`。

## 4. 积分器

### 4.1 RK4

```cpp
State rk4_step(const OdeSystem& sys, const State& x, double t, double h);
template<class RHS> State rk4_step_t(const RHS& rhs, const State& x, double t, double h);
```

模板版本避免 `std::function` 开销，适合在 Poincaré/分岔扫描的紧循环中使用。

### 4.2 Dormand–Prince 5(4)

```cpp
struct Rk45Options {
    double rtol = 1e-8, atol = 1e-10;
    double h_min = 1e-14, h_max = 0.0;  // 0 => 不额外限制
    std::size_t max_steps = 20000000;
};
struct Rk45Result {
    bool ok; std::string message;
    State x; double t, h_last;
    std::size_t accepted, rejected;
};
Rk45Result rk45_integrate(const OdeSystem&, State x, double t0, double t1,
                          double h0, const Rk45Options& = {});
```

误差按分量相对/绝对混合范数控制；成功时 `ok=true`，`x` 在 `t1`。只支持 `t1 > t0`。

### 4.3 辛积分器

```cpp
using Accel = std::function<State(const State&)>;
void verlet_step(const Accel& a, State& q, State& v, double h);
void yoshida4_step(const Accel& a, State& q, State& v, double h);
```

`a(q)` 返回 $q''$。单位质量、 separable Hamiltonian $H=|v|^2/2+V(q)$。Yoshida4 是 velocity Verlet 的三段对称组合，四阶、辛。

### 4.4 能量

```cpp
double hamiltonian_energy(const State& q, const State& v,
                          const std::function<double(const State&)>& potential);
```

## 5. Lyapunov 指数

### 5.1 连续系统：Benettin

```cpp
struct LyapunovOptions {
    double h = 0.01;
    double transient = 100.0;  // 丢弃
    double duration = 1000.0;  // 统计窗口
    std::size_t reorthonormalize_every = 1;
    int n_exponents = 0;       // 0 => 全部
};
struct LyapunovSpectrum {
    std::vector<double> exponents;  // 降序（实现末尾排序；QR 原始顺序可能不同）
    double integration_time;
    std::size_t steps, reorthonormalizations;
};
LyapunovSpectrum lyapunov_spectrum(const OdeSystem&, State x0,
                                   const LyapunovOptions& = {});
```

- 单位：**每单位时间**。
- 无解析 Jacobian 时用中心差分；此时每步 4 次 Jacobian（RK4 阶段），较慢。
- `transient` 与 `duration` 以时间计，内部按 `h` 取整步。
- 对耗散系统，指数和应接近平均散度；可作为健全性检查。

### 5.2 离散映射

```cpp
LyapunovSpectrum lyapunov_spectrum_map(const MapSystem&, State x0,
                                       std::size_t transient, std::size_t n,
                                       std::size_t reorth = 1, int n_exponents = 0);
```

单位：**每次迭代**。与 ODE 的区别是时间步固定为 1；不要把它当成“每秒”。

### 5.3 双轨迹最大指数

```cpp
double lyapunov_max_two_trajectories(const OdeSystem&, State x0, double d0,
                                     double h, double transient, double duration);
```

返回每单位时间。`d0` 建议 $10^{-8}\sim10^{-10}$（太大进入非线性区，太小被舍入噪声吞掉）。

### 5.4 Rosenstein（标量序列）

```cpp
struct RosensteinResult {
    double lambda;          // 每单位时间
    double intercept, r2;
    double dt;
    std::size_t fit_begin, fit_end;  // 索引进 k/divergence
    std::vector<double> k, divergence;
};
RosensteinResult rosenstein_lyapunov(const std::vector<double>& series,
                                     int m, int tau, int theiler, int max_k,
                                     double dt = 1.0,
                                     int fit_begin = -1, int fit_end = -1);
```

- `series` 等间隔采样，`dt` 是采样间隔；`tau`、`theiler` 以样本计。
- **返回的是 $\lambda$，不是斜率**。实现内部用 `slope/dt`；不要在外面再除以 `tau*dt`。
- 自动拟合窗口会找一个“长、正斜率、近似线性”的区间；对 Hénon 这类快速饱和的数据，自动窗口可能偏保守，建议手动看 `divergence` 曲线并同时指定 `fit_begin/fit_end`。两者要么都省略（自动），要么同时给出；只给一个会抛异常。
- 近邻搜索是 $O(N^2)$；N 很大时需 KD-tree（本库未实现）。

## 6. 相空间重构与非线性指标

### 6.1 延迟嵌入

```cpp
Embedding delay_embed(const std::vector<double>&, int m, int tau);
// emb.n = N - (m-1)tau; emb.at(i,j)
```

### 6.2 AMI 与自相关

```cpp
std::vector<double> average_mutual_information(const std::vector<double>&,
                                               int max_tau, int bins = 32);
int first_ami_minimum(const std::vector<double>& ami);
std::vector<double> autocorrelation(const std::vector<double>&, int max_lag);
int autocorrelation_1e(const std::vector<double>& ac);
```

AMI 使用等宽分箱估计；分箱数影响曲线，建议 `bins=16..64` 做敏感性检查。`first_ami_minimum` 返回第一个内部局部极小，找不到返回 -1。

### 6.3 FNN

```cpp
struct FnnResult {
    int suggested_m;
    std::vector<double> false_ratio;  // 索引 m 表示“用 m 维时仍需 m+1 的比例”
    std::vector<int> valid_points;
};
FnnResult false_nearest_neighbors(const std::vector<double>&, int max_m, int tau,
                                  int theiler = 0, double rtol = 15.0,
                                  double atol = 0.0);
```

判据：$|x_{i+m\tau}-x_{j+m\tau}|/\|\mathbf X_i-\mathbf X_j\|_m>rtol$ 记为假近邻。复杂度 $O(N^2\cdot max\_m)$。

### 6.4 关联维（GP）

```cpp
struct CorrelationResult {
    double d2, r2, r_min, r_max;
    std::vector<double> radii, log_r, log_c, local_slope;
    std::size_t fit_begin, fit_end;
};
CorrelationResult correlation_dimension(const Embedding&, double r_min=0,
                                        double r_max=0, int n_r=40,
                                        int theiler=0, double fit_lo=0.25,
                                        double fit_hi=0.75);
```

- `r_min/r_max` 为 0 时自动选分位数范围；自动范围只是起点，**必须看 `local_slope` 是否平台**。
- 复杂度 $O(N^2 n_r)$；N=8000 时几秒量级。
- Theiler 窗用于排除时间近邻；对连续流建议取 10–100 个样本，依采样率而定。

## 7. 时间序列工具

### 7.1 功率谱

```cpp
struct Spectrum { std::vector<double> frequency, power; double df; std::size_t nfft; };
Spectrum power_spectrum(const std::vector<double>&, double dt,
                        bool detrend = true, bool hann = true, bool zero_pad = true);
```

默认去均值 + Hann 窗 + 零填充到 2 的幂；`zero_pad=false` 时用精确长度的 DFT（非 2 的幂走 $O(n^2)$ 路径，适合需要严格比较相位的场景）。

### 7.2 替代数据

```cpp
std::vector<double> phase_randomized_surrogate(const std::vector<double>&, Rng&);
```

保持均值与功率谱（在精确 DFT 意义下），随机化相位；**不保持振幅分布**。更严格的 AAFT/IAAFT 未实现。

### 7.3 0-1 检验

```cpp
struct ZeroOneResult { double K; std::vector<double> K_values; };
ZeroOneResult zero_one_test(const std::vector<double>&, int n_c = 64,
                            std::uint64_t seed = ...);
```

$K\approx0$ 规则，$K\approx1$ 混沌。**采样间隔至关重要**：过采样可能让混沌序列给出 $K\approx0$（见理论文档 §5.5）。先降采样到与相关时间相当的间隔。

### 7.4 RQA

```cpp
std::vector<std::uint8_t> recurrence_matrix(const Embedding&, double epsilon,
                                            int theiler = 0);
double epsilon_for_recurrence_rate(const Embedding&, double rate, int theiler = 0);
struct RqaResult { double epsilon, recurrence_rate, determinism, laminarity,
                   entropy, trapping_time, longest_diagonal, mean_diagonal;
                   std::size_t points; };
RqaResult rqa(const Embedding&, double epsilon, int theiler = 0, int min_line = 2);
```

- `RR` 不计主对角线；`DET`/`LAM` 也只统计非对角递归点。
- 跨数据集比较时固定 `m, tau, epsilon, theiler, min_line`，或固定目标 RR。
- `ENTR` 以 nat 为单位（自然对数）。

## 8. 分岔与 Poincaré 截面

```cpp
struct BifurcationPoint { double parameter, value; };
std::vector<BifurcationPoint> bifurcation_scan_1d(
    const std::function<double(double,double)>& map, double x0,
    double p_min, double p_max, int n_params, std::size_t transient,
    std::size_t keep, std::size_t stride = 1);

struct PoincareOptions {
    int coordinate = 2; double value = 0; int direction = 0;  // +1/-1/0
    double h = 0.005, t_transient = 0, t_max = 1000;
    std::size_t max_points = 20000;
};
std::vector<State> poincare_section(const OdeSystem&, State x0,
                                    const PoincareOptions&);
```

Poincaré 截面用 RK4 单步 + 线性插值定位穿越点；`direction=+1` 上穿，`-1` 下穿，`0` 任意。插值误差 $O(h^2)$，需要高精度时减小 `h` 或用事件检测。

## 9. 控制、同步、SINDy

### 9.1 OGY

```cpp
struct OgyOptions {
    double x_star, p0, lambda_u, g;
    double tolerance = 0.1, dp_max = 0.05;
    std::size_t steps = 100000;
};
OgyResult ogy_control_1d(const std::function<double(double,double)>& map,
                         double x0, const OgyOptions&);
```

`lambda_u = ∂f/∂x`，`g = ∂f/∂p` 在目标点取值。控制律 $\delta p=-\lambda_u(x-x^*)/g$，仅在邻域内、且受 `dp_max` 限制。

### 9.2 Pyragas（映射）

```cpp
std::vector<double> pyragas_map_control(const std::function<double(double)>& map,
                                        double x0, int delay, double gain,
                                        std::size_t steps);
```

$x_{n+1}=f(x_n)+K(x_{n-d}-x_n)$。`delay` 以迭代步计。连续 ODE 版本未实现。

### 9.3 Pecora–Carroll

```cpp
struct SyncResult { std::vector<double> t, error; double final_error,
                    max_error_after_transient, settling_time; bool synchronized; };
SyncResult pecora_carroll_lorenz(const Lorenz&, const State& drive0,
                                 const State& response0, double h = 0.005,
                                 double t_transient = 5, double t_total = 40);
```

`drive0` 三维，`response0` 两维（$y,z$）。同步判据是最终误差 $<10^{-6}$；`settling_time` 是最后一次进入 $10^{-6}$ 后的连续时间。

### 9.4 SINDy

```cpp
struct SindyOptions {
    int degree = 2;
    bool include_constant = true, include_sine = false, include_cosine = false;
    double threshold = 0.1, ridge = 1e-10;
    int max_iter = 20;
};
struct SindyModel {
    std::vector<SindyTerm> terms;  // name + evaluate
    Matrix coefficients;           // terms x dim
    State predict_derivative(const State&) const;
    std::string pretty() const;
};
SindyModel sindy_fit(const std::vector<State>& X,
                     const std::vector<State>& dX,
                     const SindyOptions& = {});
std::vector<State> finite_difference_derivatives(const std::vector<State>&, double dt);
```

STLSQ：列归一化 → 岭正则最小二乘 → 阈值化 → 迭代。`threshold` 是绝对系数阈值；量纲变化大的系统应先归一化状态。候选库目前是多项式（次数可调）+ 可选 sin/cos。

## 10. C API

头文件 [`include/chaos/chaos_c.h`](../include/chaos/chaos_c.h)，链接 `libchaos`（C++ 运行时）。

```c
const char* chaos_version(void);
const char* chaos_last_error(void);
int chaos_lorenz_trajectory(double sigma, double rho, double beta,
                            const double x0[3], double h, size_t n,
                            double* out_t, double* out_xyz);
int chaos_lorenz_lyapunov(double sigma, double rho, double beta,
                          const double x0[3], double h, double transient,
                          double duration, double out_exponents[3]);
double chaos_logistic_lyapunov(double r, size_t transient, size_t n);
int chaos_correlation_dimension(const double* series, size_t n, int m, int tau,
                                int theiler, double* out_d2);
int chaos_zero_one_test(const double* series, size_t n, int n_c, double* out_k);
int chaos_rosenstein_lambda(const double* series, size_t n, int m, int tau,
                            int theiler, int max_k, double dt, double* out_lambda);
```

失败返回非 0，错误信息由 `chaos_last_error()` 提供；`chaos_logistic_lyapunov` 失败返回 NaN。

## 11. CLI 速查

```bash
chaos_cli models
chaos_cli simulate --system lorenz --n 5000 --h 0.01 --transient 50 --out lorenz.csv
chaos_cli simulate --map henon --n 5000 --transient 1000 --out henon.csv
chaos_cli lyapunov --system lorenz --h 0.01 --transient 50 --duration 500
chaos_cli lyapunov --map henon --n 100000
chaos_cli bifurcation --map logistic --pmin 2.5 --pmax 4 --nparams 1500 \
  --transient 1000 --keep 80 --out bif.csv
chaos_cli lyapunov-scan --map logistic --nparams 1500 --n 20000 --out lyap.csv
chaos_cli poincare --system lorenz --coord 2 --value 27 --direction -1 \
  --h 0.005 --transient 20 --t 200 --out section.csv
chaos_cli ami --input series.csv --column 1 --max-tau 60
chaos_cli fnn --input series.csv --column 1 --tau 17 --max-m 8
chaos_cli dimension --input series.csv --column 1 --m 3 --tau 17 --theiler 50
chaos_cli rosenstein --input series.csv --column 1 --m 3 --tau 17 --theiler 100 \
  --max-k 300 --dt 0.01
chaos_cli zeroone --input series.csv --column 1 --stride 10
chaos_cli rqa --input series.csv --column 1 --stride 10 --m 3 --tau 3 --rate 0.03
chaos_cli spectrum --input series.csv --column 1 --dt 0.01 --out psd.csv
chaos_cli sindy --system lorenz --n 5000 --transient 5000 --h 0.002
chaos_cli sync --t 40 --h 0.005
```

CSV 约定：`--column` 是 **0 基**列号；数据文件通常 `t,x0,...`，所以标量序列用 `--column 1`。首行若不能解析为数字会被当作表头跳过。

## 12. 复杂度与扩展

| 操作 | 复杂度 | 备注 |
|---|---|---|
| RK4 一步 | $O(d)$ 次 RHS | 模板版更快 |
| Benettin | $O(N d^2 k)$ + Jacobian | 解析 Jacobian 快很多 |
| FNN | $O(N^2 m)$ | 大 N 需 KD-tree |
| GP | $O(N^2 n_r)$ | 大 N 需树/采样 |
| Rosenstein | $O(N^2+N K)$ | 大 N 需 KD-tree |
| RQA | $O(N^2)$ 内存 | N 不宜过大 |
| SINDy | $O(N p^2)$ | $p$ 为候选库大小 |

扩展一个新 ODE：实现 `f`/`jacobian`、写 `ode()`、加入 `ode_catalog()` 与 `default_initial_state()`；测试里加 Jacobian 对照和已知不变量。扩展一个新的指标：优先在 `tests/` 里加一个能区分“已知规则/混沌”的断言，再写文档。

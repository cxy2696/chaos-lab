# 参考文献与溯源

本页是 `docs/chaos-systems.md`、`docs/algorithms.md` 与 README 中 `[n]` 编号的出处。DOI 均于 2026-11-19 通过 Crossref 核对；BibTeX 位于 [`references/references.bib`](../references/references.bib)。若 DOI 未来失效，可用标题 + 期刊 + 卷页检索。

## 基础、历史与教材

1. **Poincaré, H. (1890).** *Sur le problème des trois corps et les équations de la dynamique.* Acta Mathematica 13, 1–270. [期刊卷](https://link.springer.com/journal/11511/volumes-and-issues/13-1) · BibTeX `poincare1890`  
   三体问题中复杂轨道的起点；混沌思想的历史源头之一。
2. **Lorenz, E. N. (1963).** *Deterministic Nonperiodic Flow.* J. Atmos. Sci. 20(2), 130–141. [doi:10.1175/1520-0469(1963)020<0130:DNF>2.0.CO;2](https://doi.org/10.1175/1520-0469(1963)020%3C0130:DNF%3E2.0.CO;2) · `lorenz1963`  
   三维对流模型；蝴蝶效应的经典来源。
3. **Lorenz, E. N. (1965).** *A study of the predictability of a 28-variable atmospheric model.* Tellus 17(3), 321–333. [doi:10.1111/j.2153-3490.1965.tb01424.x](https://doi.org/10.1111/j.2153-3490.1965.tb01424.x) · `lorenz1965`  
   把初值敏感性变成可操作的“可预报性”研究，集合预报的思想来源。
4. **Ruelle, D., & Takens, F. (1971).** *On the nature of turbulence.* Comm. Math. Phys. 20(3), 167–192. [doi:10.1007/BF01646553](https://doi.org/10.1007/BF01646553) · `ruelle1971`  
   准周期路径：不动点 → 极限环 → 二维环面 → 奇怪吸引子。
5. **Li, T.-Y., & Yorke, J. A. (1975).** *Period Three Implies Chaos.* Amer. Math. Monthly 82(10), 985–992. [doi:10.1080/00029890.1975.11994008](https://doi.org/10.1080/00029890.1975.11994008) · `li1975`  
   “混沌”一词进入数学文献；周期 3 蕴含任意周期与不可数混沌集。
6. **May, R. M. (1976).** *Simple mathematical models with very complicated dynamics.* Nature 261, 459–467. [doi:10.1038/261459a0](https://doi.org/10.1038/261459a0) · `may1976`  
   Logistic 映射的普及性综述；倍周期与混沌的参数图景。
7. **Feigenbaum, M. J. (1978).** *Quantitative universality for a class of nonlinear transformations.* J. Stat. Phys. 19(1), 25–52. [doi:10.1007/BF01020332](https://doi.org/10.1007/BF01020332) · `feigenbaum1978`  
   δ ≈ 4.6692 与 α ≈ 2.5029 的普适常数。
8. **Feigenbaum, M. J. (1979).** *The universal metric properties of nonlinear transformations.* J. Stat. Phys. 21(6), 669–706. [doi:10.1007/BF01107909](https://doi.org/10.1007/BF01107909) · `feigenbaum1979`  
   普适性的重整化群解释与函数方程。
9. **Sparrow, C. (1982).** *The Lorenz Equations: Bifurcations, Chaos, and Strange Attractors.* Springer, Applied Mathematical Sciences 41. [doi:10.1007/978-1-4612-5767-7](https://doi.org/10.1007/978-1-4612-5767-7) · `sparrow1982`  
   Lorenz 系统的系统分析；谱、维数、参数图景的经典参考。
10. **Strogatz, S. H. (2015).** *Nonlinear Dynamics and Chaos.* 2nd ed., Westview Press. ISBN 978-0813349107. · `strogatz2015`  
    首选入门教材；本文许多定义与例题的表述参照该书。

## 经典模型

11. **Hénon, M. (1976).** *A two-dimensional mapping with a strange attractor.* Comm. Math. Phys. 50(1), 69–77. [doi:10.1007/BF01608556](https://doi.org/10.1007/BF01608556) · `henon1976`
12. **Rössler, O. E. (1976).** *An equation for continuous chaos.* Phys. Lett. A 57(5), 397–398. [doi:10.1016/0375-9601(76)90101-8](https://doi.org/10.1016/0375-9601(76)90101-8) · `rossler1976`
13. **Chua, L. O., Komuro, M., & Matsumoto, T. (1986).** *The double scroll family.* IEEE Trans. Circuits Syst. 33(11), 1072–1118. [doi:10.1109/TCS.1986.1085869](https://doi.org/10.1109/TCS.1986.1085869) · `chua1986`
14. **Chenciner, A., & Montgomery, R. (2000).** *A Remarkable Periodic Solution of the Three-Body Problem in the Case of Equal Masses.* Ann. of Math. 152(3), 881–901. [doi:10.2307/2661357](https://doi.org/10.2307/2661357) · `chenciner2000`  
    8 字形三体编舞；本库 `three_body` 的初值与周期来源。
15. **Arnold, V. I., & Avez, A. (1968).** *Ergodic Problems of Classical Mechanics.* W. A. Benjamin. · `arnold1968`  
    环面线性映射（猫映射）与遍历理论的经典背景。

## 数值方法与非线性时间序列

16. **Benettin, G., Galgani, L., Giorgilli, A., & Strelcyn, J.-M. (1980).** *Lyapunov Characteristic Exponents for smooth dynamical systems and for Hamiltonian systems; a method for computing all of them. Part 1: Theory.* Meccanica 15(1), 9–20. [doi:10.1007/BF02128236](https://doi.org/10.1007/BF02128236) · `benettin1980`  
    沿变分方程积分 + 周期 QR/Gram–Schmidt 重正交化，即本库 `lyapunov_spectrum` 的算法基础。
17. **Wolf, A., Swift, J. B., Swinney, H. L., & Vastano, J. A. (1985).** *Determining Lyapunov exponents from a time series.* Physica D 16(3), 285–317. [doi:10.1016/0167-2789(85)90011-9](https://doi.org/10.1016/0167-2789(85)90011-9) · `wolf1985`  
    从标量时间序列估计 Lyapunov 指数的经典方法；需要嵌入与演化。
18. **Sano, M., & Sawada, Y. (1985).** *Measurement of the Lyapunov Spectrum from a Chaotic Time Series.* Phys. Rev. Lett. 55(10), 1082–1085. [doi:10.1103/PhysRevLett.55.1082](https://doi.org/10.1103/PhysRevLett.55.1082) · `sano1985`  
    从时间序列重构局部 Jacobian 以估计谱。
19. **Rosenstein, M. T., Collins, J. J., & De Luca, C. J. (1993).** *A practical method for calculating largest Lyapunov exponents from small data sets.* Physica D 65(1–2), 117–134. [doi:10.1016/0167-2789(93)90009-P](https://doi.org/10.1016/0167-2789(93)90009-P) · `rosenstein1993`  
    最近邻轨道平均发散；本库 `rosenstein_lyapunov`。
20. **Takens, F. (1981).** *Detecting strange attractors in turbulence.* In *Dynamical Systems and Turbulence, Warwick 1980*, LNM 898, 366–381. Springer. [doi:10.1007/BFb0091924](https://doi.org/10.1007/BFb0091924) · `takens1981`  
    延迟嵌入定理：m ≥ 2d+1（泛型）时重构吸引子与原吸引子微分同胚。
21. **Sauer, T., Yorke, J. A., & Casdagli, M. (1991).** *Embedology.* J. Stat. Phys. 65(3–4), 579–616. [doi:10.1007/BF01053745](https://doi.org/10.1007/BF01053745) · `sauer1991`  
    嵌入的严格处理与噪声下的推广。
22. **Fraser, A. M., & Swinney, H. L. (1986).** *Independent coordinates for strange attractors from mutual information.* Phys. Rev. A 33(2), 1134–1140. [doi:10.1103/PhysRevA.33.1134](https://doi.org/10.1103/PhysRevA.33.1134) · `fraser1986`  
    用平均互信息（AMI）第一极小值选延迟 τ。
23. **Kennel, M. B., Brown, R., & Abarbanel, H. D. I. (1992).** *Determining embedding dimension for phase-space reconstruction using a geometrical construction.* Phys. Rev. A 45(6), 3403–3411. [doi:10.1103/PhysRevA.45.3403](https://doi.org/10.1103/PhysRevA.45.3403) · `kennel1992`  
    假近邻法（FNN）选嵌入维数 m。
24. **Grassberger, P., & Procaccia, I. (1983).** *Measuring the strangeness of strange attractors.* Physica D 9(1–2), 189–208. [doi:10.1016/0167-2789(83)90298-1](https://doi.org/10.1016/0167-2789(83)90298-1) · `gp1983`  
    关联积分 C(r) ∝ r^{D₂} 与关联维数；本库 `correlation_dimension`。
25. **Eckmann, J.-P., Kamphorst, S. O., & Ruelle, D. (1987).** *Recurrence Plots of Dynamical Systems.* Europhys. Lett. 4(9), 973–977. [doi:10.1209/0295-5075/4/9/004](https://doi.org/10.1209/0295-5075/4/9/004) · `eckmann1987`
26. **Marwan, N., Romano, M. C., Thiel, M., & Kurths, J. (2007).** *Recurrence plots for the analysis of complex systems.* Phys. Rep. 438(5–6), 237–329. [doi:10.1016/j.physrep.2006.11.001](https://doi.org/10.1016/j.physrep.2006.11.001) · `marwan2007`  
    递归图与 RQA 的 RR/DET/LAM/ENTR/TT 定义。
27. **Theiler, J., Eubank, S., Longtin, A., Galdrikian, B., & Farmer, J. D. (1992).** *Testing for nonlinearity in time series: the method of surrogate data.* Physica D 58(1–4), 77–94. [doi:10.1016/0167-2789(92)90102-S](https://doi.org/10.1016/0167-2789(92)90102-S) · `theiler1992`  
    替代数据检验：先排除线性随机过程，再谈“非线性/混沌”。
28. **Gottwald, G. A., & Melbourne, I. (2004).** *A new test for chaos in deterministic systems.* Proc. R. Soc. A 460, 603–611. [doi:10.1098/rspa.2003.1183](https://doi.org/10.1098/rspa.2003.1183) · `gottwald2004`  
    0-1 检验；K≈0 规则、K≈1 混沌。
29. **Kaplan, J. L., & Yorke, J. A. (1979).** *Chaotic behavior of multidimensional difference equations.* LNM 730, 204–227. Springer. [doi:10.1007/BFb0064319](https://doi.org/10.1007/BFb0064319) · `kaplan1979`  
    Kaplan–Yorke 维数 D_KY = j + Σλ_i/|λ_{j+1}|，用于把 Lyapunov 谱与分形维数联系起来。
30. **Kantz, H., & Schreiber, T. (2004).** *Nonlinear Time Series Analysis.* 2nd ed., Cambridge University Press. ISBN 978-0521529020. · `kantz2004`  
    时间序列方法的系统参考；也明确讨论噪声、非平稳与替代数据。

## 控制、同步与数据驱动建模

31. **Ott, E., Grebogi, C., & Yorke, J. A. (1990).** *Controlling chaos.* Phys. Rev. Lett. 64(11), 1196–1199. [doi:10.1103/PhysRevLett.64.1196](https://doi.org/10.1103/PhysRevLett.64.1196) · `ogy1990`  
    OGY：利用不稳定周期轨道，微小参数扰动稳定混沌。
32. **Pyragas, K. (1992).** *Continuous control of chaos by self-controlling feedback.* Phys. Lett. A 170(6), 421–428. [doi:10.1016/0375-9601(92)90745-8](https://doi.org/10.1016/0375-9601(92)90745-8) · `pyragas1992`  
    延迟反馈 u(t)=K[x(t−T)−x(t)]；不需要目标轨道的解析形式。
33. **Pecora, L. M., & Carroll, T. L. (1990).** *Synchronization in chaotic systems.* Phys. Rev. Lett. 64(8), 821–824. [doi:10.1103/PhysRevLett.64.821](https://doi.org/10.1103/PhysRevLett.64.821) · `pecora1990`  
    驱动–响应同步；本库 `pecora_carroll_lorenz`。
34. **Cuomo, K. M., & Oppenheim, A. V. (1993).** *Circuit implementation of synchronized chaos with applications to communications.* Phys. Rev. Lett. 71(1), 65–68. [doi:10.1103/PhysRevLett.71.65](https://doi.org/10.1103/PhysRevLett.71.65) · `cuomo1993`
35. **Kocarev, L., & Parlitz, U. (1995).** *General Approach for Chaotic Synchronization with Applications to Communication.* Phys. Rev. Lett. 74(25), 5028–5031. [doi:10.1103/PhysRevLett.74.5028](https://doi.org/10.1103/PhysRevLett.74.5028) · `kocarev1995`
36. **Brunton, S. L., Proctor, J. L., & Kutz, J. N. (2016).** *Discovering governing equations from data by sparse identification of nonlinear dynamical systems.* PNAS 113(15), 3932–3937. [doi:10.1073/pnas.1517384113](https://doi.org/10.1073/pnas.1517384113) · `brunton2016`  
    SINDy：稀疏回归 + 候选函数库；本库 `sindy_fit`。
37. **Pathak, J., Hunt, B. R., Girvan, M., Lu, Z., & Ott, E. (2018).** *Model-Free Prediction of Large Spatiotemporally Chaotic Systems from Data: A Reservoir Computing Approach.* Phys. Rev. Lett. 120, 024102. [doi:10.1103/PhysRevLett.120.024102](https://doi.org/10.1103/PhysRevLett.120.024102) · `pathak2018`
38. **Schmid, P. J. (2010).** *Dynamic mode decomposition of numerical and experimental data.* J. Fluid Mech. 656, 5–28. [doi:10.1017/S0022112010001217](https://doi.org/10.1017/S0022112010001217) · `schmid2010`
39. **Bradley, E., & Kantz, H. (2015).** *Nonlinear time-series analysis revisited.* Chaos 25(9), 097610. [doi:10.1063/1.4917289](https://doi.org/10.1063/1.4917289) · `bradley2015`  
    对“在噪声、非平稳数据上直接套非线性指标”的清醒警告；金融/生理数据尤其应读。

## 应用、普适性与警告

40. **Toth, Z., & Kalnay, E. (1997).** *Ensemble Forecasting at NCEP and the Breeding Method.* Mon. Wea. Rev. 125(12), 3297–3319. [doi:10.1175/1520-0493(1997)125<3297:EFANAT>2.0.CO;2](https://doi.org/10.1175/1520-0493(1997)125%3C3297:EFANAT%3E2.0.CO;2) · `toth1997`
41. **Kalnay, E. (2002).** *Atmospheric Modeling, Data Assimilation and Predictability.* Cambridge University Press. [doi:10.1017/CBO9780511802270](https://doi.org/10.1017/CBO9780511802270) · `kalnay2002`
42. **Belbruno, E., & Miller, J. K. (1993).** *Sun-perturbed Earth-to-moon transfers with ballistic capture.* J. Guidance Control Dyn. 16(4), 770–775. [doi:10.2514/3.21079](https://doi.org/10.2514/3.21079) · `belbruno1993`  
    低能量转移与弹道捕获；Hiten 任务的理论背景之一。
43. **Aref, H. (1984).** *Stirring by chaotic advection.* J. Fluid Mech. 143, 1–21. [doi:10.1017/S0022112084001233](https://doi.org/10.1017/S0022112084001233) · `aref1984`
44. **Goldberger, A. L., Amaral, L. A. N., Hausdorff, J. M., Ivanov, P. Ch., Peng, C.-K., & Stanley, H. E. (2002).** *Fractal dynamics in physiology: Alterations with disease and aging.* PNAS 99(Suppl 1), 2466–2472. [doi:10.1073/pnas.012582499](https://doi.org/10.1073/pnas.012582499) · `goldberger2002`
45. **Solak, E., Çokal, C., Yildiz, O. T., & Biyikoğlu, T. (2010).** *Cryptanalysis of Fridrich's Chaotic Image Encryption.* Int. J. Bifurcation Chaos 20(5), 1405–1413. [doi:10.1142/S0218127410026563](https://doi.org/10.1142/S0218127410026563) · `solak2010`  
    混沌图像加密被攻破的代表性案例；不要用于真实安全。
46. **Stauffer, D., & Aharony, A. (1994).** *Introduction to Percolation Theory.* 2nd ed., Taylor & Francis. ISBN 978-0748400270. · `stauffer1994`
47. **Newhouse, S., Ruelle, D., & Takens, F. (1978).** *Occurrence of strange Axiom A attractors near quasi periodic flows on T^m, m ≥ 3.* Comm. Math. Phys. 64(1), 35–40. [doi:10.1007/BF01940759](https://doi.org/10.1007/BF01940759) · `newhouse1978`
48. **Ott, E. (2002).** *Chaos in Dynamical Systems.* 2nd ed., Cambridge University Press. [doi:10.1017/CBO9780511803260](https://doi.org/10.1017/CBO9780511803260) · `ott2002`  
    进阶教材；控制、同步、阴影引理与统计性质的系统讨论。
49. **Pomeau, Y., & Manneville, P. (1980).** *Intermittent transition to turbulence in dissipative dynamical systems.* Comm. Math. Phys. 74(2), 189–197. [doi:10.1007/BF01197757](https://doi.org/10.1007/BF01197757) · `pomeau1980`
50. **Wilson, K. G. (1971).** *Renormalization Group and Critical Phenomena. I. Renormalization Group and the Kadanoff Scaling Picture.* Phys. Rev. B 4(9), 3174–3183. [doi:10.1103/PhysRevB.4.3174](https://doi.org/10.1103/PhysRevB.4.3174) · `wilson1971`  
    重整化群与临界指数普适性的基础；与 Feigenbaum 普适性的对照。

---

## 编号使用约定

- 正文中的 `[n]` 只指向本页编号；不把网页博客当作首要来源。
- 数值断言若来自本仓库计算，会在 `docs/data-provenance.md` 的“断言 → 来源 → 复现命令”表中单独标注，不冒充文献值。
- 书籍页码/版次以 BibTeX 为准；不同版次页码可能不同。

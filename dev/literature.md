# Summary of Literature on Spatial and Spatiotemporal ARCH Models

This document provides a technical summary of spatial and spatiotemporal autoregressive conditional heteroscedasticity (ARCH) models, incorporating explicit mathematical formulations.

---

## 1. Core Model Definition
For a univariate stochastic process $\{Y(s) : s \in D\}$, the general approach to modeling spatial ARCH-type processes is given by:
$$Y = \text{diag}(h)^{1/2}\varepsilon$$
where $h = (h(s_1), \dots, h(s_n))'$ is the vector of conditional variances and $\varepsilon = (\varepsilon(s_1), \dots, \varepsilon(s_n))'$ is a noise component (typically i.i.d. with mean 0 and variance 1).

## 2. Model Variants

### A. Spatial ARCH (spARCH)
The volatility $h_O$ is defined as:
$$h_O = \alpha \mathbf{1}_n + \rho W \text{diag}(Y)Y$$
Alternatively, in terms of squared observations $Y^{(2)} = (Y(s_1)^2, \dots, Y(s_n)^2)'$:
$$Y^{(2)} = \alpha (I - A)^{-1} \varepsilon^{(2)}, \quad \text{where } A = \rho \text{diag}(\varepsilon^{(2)})W$$

### B. Spatial GARCH (spGARCH)
Incorporating past conditional variances, the model is:
$$h = \alpha \mathbf{1}_n + W_1 Y^{(2)} + W_2 h$$
The conditional variance $h$ can be specified as:
$$h = (I - W_1 E - W_2)^{-1} \alpha$$
where $E = \text{diag}(\varepsilon(s_1)^2, \dots, \varepsilon(s_n)^2)$.

### C. Exponential Spatial ARCH (E-spARCH)
To prevent negative variances without strict parameter constraints, the logarithm of the volatility is modeled:
$$\ln h_E = \alpha \mathbf{1}_n + \rho W g_b(\varepsilon)$$
where $g_b(\varepsilon) = (\ln |\varepsilon(s_1)|^b, \dots, \ln |\varepsilon(s_n)|^b)'$. This can be rewritten as:
$$\ln h_E = S(\alpha \mathbf{1}_n + \rho b W \ln |Y|) \quad \text{with } S = (I + \frac{1}{2}\rho b W)^{-1}$$

### D. SARspARCH Model
Used as an error model for spatial regression:
$$Y = \lambda B Y + X\beta + u, \quad u = \text{diag}(h)^{1/2}\varepsilon$$
Here, $\lambda B$ captures spatial dependence in the mean, while the spatial ARCH model in $h$ captures spatial heteroscedasticity.

## 3. Estimation Strategy
Parameters $\vartheta = (\rho, \lambda, \alpha)'$ are typically estimated using a non-linear least-squares approach (NLSE) or Quasi-Maximum Likelihood (QML). For the logarithmic transformation $H = \ln(Y^{(2)}) - E(\ln(\varepsilon^{(2)}))$, the estimators minimize:
$$Q_n(\vartheta) = \frac{1}{n} \sum_{i=1}^n \left( H_i - \tau^{-1}(\alpha c_i(\lambda) + \rho d_i(\lambda)' \tilde{\gamma}(H + c\mathbf{1})) \right)^2$$
where $\tau(x) = f(\exp(x))$ and $\tilde{\gamma}$ is a function relating the log-volatility to past squared observations.

## 4. Spatiotemporal Extensions
For spatiotemporal processes $Y_t(s)$, the model can be represented as:
$$\ln h_t(s_i) = \sum_{l=1}^p \sum_{j=1}^n \rho_{l0} m_{l,ij} \ln Y_t^2(s_j) + \gamma_0 \ln Y_{t-1}^2(s_i) + \sum_{l=1}^p \sum_{j=1}^n \delta_{l0} m_{l,ij} \ln Y_{t-1}^2(s_i) + x_t(s_i)' \beta_0 + \mu_0(s_i) + \alpha_{t0}$$
This incorporates spatial ($m_{l,ij}$), temporal ($\gamma_0$), and spatiotemporal ($\delta_{l0}$) dependencies explicitly.

---
*Sources: Otto & Schmid (2023), Otto (2019), and Otto et al. (2025).*

# **Robustifying Marginal Linear Models for Correlated Responses Using a Constructive Multivariate Huber Distribution**

## **Overview**
We present a novel approach to analyzing correlated data through a robust marginal model incorporating a multivariate Huber distribution. This model provides robustness against outliers and features a **tuning parameter** to control the degree of robustness, where larger values approximate normality and smaller values enhance robustness. 

Unlike prior models with constant tuning parameters (often fixed near 2), our approach allows **subject-specific tuning parameters**, enabling analysts to adjust the influence of individual subjects based on their observations. Parameters are estimated using the exact likelihood function via the **Hamiltonian Monte Carlo (HMC)** algorithm implemented in Stan. 

Additionally, the **modified Cholesky decomposition** is employed to construct the multivariate density function, transforming constrained dispersion matrix parameters into two vectors of unconstrained parameters for easier interpretation and estimation. The model supports **flexible covariance structures**, including compound symmetry, autoregressive, and random coefficients.

---

## **Features**
- **Robust Modeling**: Introduces robustness against outliers for correlated responses through a constructive multivariate Huber distribution.
- **Flexible Parameterization**: Employs the modified Cholesky decomposition to estimate covariance matrices efficiently.
- **Subject-Varying Robustness**: Allows tuning parameters to vary across subjects, enabling individualized control over robustness.
- **Hamiltonian Monte Carlo Sampling**: Uses the RStan library for computational efficiency and exact likelihood estimation.
- **Diagnostics and Model Selection**: Facilitates model evaluation and selection with robust diagnostic measures.

---

## 📄 Reference to the Paper
This repository implements the methodology introduced in the following paper:

**"Robustifying Marginal Linear Models for Correlated Responses Using a Constructive Multivariate Huber Distribution"**  
*Raziyeh Mohammadi¹, Iraj Kazemi²*  

¹ Duke-NUS Medical School, National University of Singapore, Singapore  
² Department of Statistics, Faculty of Mathematics & Statistics, University of Isfahan, Iran  

📖 **Published in:** *Statistical Analysis and Modeling*  
📅 **First Published:** 31 January 2025  
🔗 **DOI:** [https://doi.org/10.1002/sam.70011](https://doi.org/10.1002/sam.70011)  
📚 **Volume:** 18, **Issue:** 1, **Article:** e70011 (February 2025)
---


Contact
For any questions or feedback, feel free to reach out at:

📧 Email: raziyeh.mohammadi@duke-nus.edu.sg

---

## **Installation**
To install the necessary packages, run the following command in R:

```R
install.packages(c("rstan", "tidyr", "dplyr"))
---





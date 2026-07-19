$\Sigma^{(t)}$: the covariance we are trying to estimate at iteration $t$ of EM algorithm.
$\beta^{(t)}$: choice covariate vector.
$Z_i$: latent utiilty of observation $i$.
$\mu_i^{(t)} = E[Z_i | X_i, Y_i, \beta^{(t)}, \Sigma^{(t)}]$: conditional mean of latent utility given data and parameters at iteration $t$ of EM algorithm.
$S_i^{(t)}=Var(Z_i | X_i, Y_i \beta^{(t)}, \Sigma^{(t)})$ : conditional covariance of latent utility 
$X_i$: covariate $i$.
$\hat{S}^{(t+1)}= E[S|(X_i, Y_i)_{i=1}^n, \beta^{(t)}, \Sigma^{(t)}] = \frac{1}{n}\sum_{i=1}^n S_i^{(t)} + (\mu_i^{(t)} - X_i\beta) (\mu_i^{(t)} - X_i\beta)^T$: conditional sample covariance

To update $\Sigma^{(t)}$ we solve
$\max_{\Sigma} -\log\det(\Sigma) - tr(\Sigma^{-1} \hat{S})$, such that $tr(\Sigma) = m$.

Do a eigendecomposition of both $\Sigma$ and $\hat{S}$, and have them share eigen vectors so that:
$\hat{S} = U diag(s_1, \ldots, s_m) U^T$
$\Sigma = U diag(\sigma_1, \ldots, \sigma_m) U^T$

Then the optimization problem is
$\max_{\sigma_i} \sum_{i=1}^m ( \log \sigma_i + \frac{s_i}{\sigma_i})$, such that $\sum_i^m \sigma_i = p$.

For Lagrange multiplier $\lambda$ define
$$
\sigma_i(\lambda) = \begin{cases}
	\lambda \ne 0 & \frac{2s_i}{1 + \sqrt{1 + 4\lambda s_i}} 
	\\
	\lambda = 0 & \sigma_i = s_i
\end{cases}
$$
Solve the equation 
$$
\sum_i^m \sigma_i (\lambda) = m
$$
for $\lambda$, get $\sigma_i(\lambda)$ for $i = 1, \ldots, m$, and plug those into the decomposition to get the maximized $\Sigma$. 


Cunningham et al 
$N$: number of observations.
$m$: number of choices not including the base case.

$y_i \in \{0, \ldots, m\}$: choice of observation $i$.
$z_i \in \mathbb{R}^m$: relative choice utilities of observation $i$.
$A \in \mathbb{R}^{m\times m}$: constraint matrix of obs $i$ utilities, given choice $y_i$. 
$A_j$: $j$th constraint of observation $i$, $p \times 1$ vector.
$e_j$: $j$th unit vector.

$\Sigma$: approximate covariance for obs $i$.
$\mu$: approximate mean for obs $i$.
$\tau_j$: site $j$ precision for obs $i$.
$\eta_j$: site $j$ mean for obs $i$.
$\epsilon_j$: site $j$ damping parameter

For a matrix $M$, $M_k$ denotes the $k$th column of $M$.

constraint matrix $A(y) = -I_{m-1} + 1_{m-1} e_y^T + e_ye_y^T$.

For the non base case, the EP approximation involves computing $a_j(y)^T \Sigma a_j(y)$, where $a_j(y) = e_y - e_j$, $j \ne y$, is a contrast constraint, the $j$th row of $A(y)$. Then define 
$s_{yj}^2 = a_j(y)^T \Sigma a_j(y) = \Sigma_{yy} - 2\Sigma_{y j} + \Sigma_{jj}$, 
$g_j(y) = \Sigma A(y)_j = \Sigma_y - \Sigma_j$. 
$d_{yj} = \mu_y - \mu_j$.

Meanwhile in the base case where $y=0$, $A(y)=I_m$, so
$s_{yj}^2 = a_j(y)^T \Sigma a_j(y) = \Sigma_{jj}$ and 
$g_j(y) = \Sigma a_j(y) = \Sigma_j$.  
$d_{yj} = \mu_j$.

```pseudo
\begin{algorithm}
\caption{EPMNP: Expectation propagation for multinomial probit moment approximation}
\begin{algorithmic}
	\Input TMVN mean $\mu_0$, covariance $\Sigma_0$, non base choice $y \in \{1, \ldots, m\}$, 
	\Output approximate mean $\mu$ and covariance $\Sigma$
	
	\State $\mu \gets \mu_0$
	\State $\Sigma \gets \Sigma_0$
	
	\While{not converged}
		\For{$j=1:p$}
				\State 
					$
						\tau^{/j} = (
							s_{yj}^2
						)^{-1} 
						- \tau_j
					$ 
				\State 
					$
						\eta^{/j} = 
							(\mu_y - \mu_j)
							(s_{yj}^2)^{-1}
						- \eta_j
					$
				
				\State 
					$\hat{\mu} = \mathbb{E}[Q^{/j}]$
				
				\State
					${\hat{\sigma}^2} = \mathbb{Var}[Q^{/j}]$
				\Comment{cavity density diagonal covariance}
				
				\State
					$
						\Delta \tau_j = 
							\epsilon_j
							(
								{\hat{\sigma}^2_j} -
								\tau^{/j}
							)^{-1} - \tau_j
					$ 
				
				\State
					$
						\Delta \eta_j = 
							\epsilon_j 
							(
								\hat{\mu}_j {\hat{\sigma}^2_j}^{-1} - 
								\eta^{/j} - \eta_j
							)
					$
				\State
					$ \tau_j += \Delta \tau_j$
				\State
					$ \eta_j += \Delta \tau_j$
				\State 
					$
						g_j = \Sigma_y - \Sigma_j
					$
				\State
					$
						\Sigma -= \Delta \tau_j(1 + \Delta \tau_j s_{yj}^2)^{-1} g_j g_j^T
					$
				\State
					$
						\mu += 
						(
						\Delta \eta_j - \Delta \tau_j (\mu_y - \mu_j)) 
						(1 + \Delta \tau_j s_{yj}^2)^{-1} g_j
					$
		\EndFor
	\EndWhile
\end{algorithmic}
\end{algorithm}
```



Pakman and Paninski (2014)

$N$: number of observations.
$m$: number of choices not including the base case.

$y_i \in \{0, \ldots, m\}$: choice of observation $i$.
$z_i \in \mathbb{R}^m$: relative choice utilities of observation $i$.
$e_j$: $j$th unit vector.

$\mu_0$: untruncated mean for obs $i$.
$\Sigma_0$: untruncated covariance for obs $i$.
$P_0 = \Sigma_0^{-1}$: untruncated precision for obs $i$.
$R$: upper triangular Cholesky factor of the precision, $P_0 = R^T R$.
$x = R(z-\mu_0)$: utility in whitened coordinates.
$v$: Hamiltonian velocity in whitened coordinates.
$T$: Hamiltonian integration time.
$B$: number of burn-in transitions.
$L$: number of retained Monte Carlo samples.

For reference identification, choice $0$ is the base choice and $z_j$ is the utility of non-base choice $j$ relative to the base choice. The observed choice defines the truncation region

$$
\mathcal{C}(y) =
\begin{cases}
\{z : z_j \leq 0,\ j=1,\ldots,m\}, & y=0,\\
\{z : z_y \geq 0,\ z_y-z_j \geq 0,\ j\neq y\}, & y\in\{1,\ldots,m\}.
\end{cases}
$$

Write every constraint as $a_j(y)^Tz \geq 0$. The reference-identified constraint normals are

$$
a_j(0) = -e_j
$$

for the base choice, and

$$
a_j(y) =
\begin{cases}
e_y, & j=y,\\
e_y-e_j, & j\neq y
\end{cases}
$$

for a non-base choice. Thus the constraints do not need to be stored as an arbitrary dense matrix.

Since $z=\mu_0+R^{-1}x$, constraint $j$ becomes

$$
f_j^Tx+g_j \geq 0,
\qquad
f_j=R^{-T}a_j(y),
\qquad
g_j=a_j(y)^T\mu_0.
$$

In whitened coordinates the unconstrained target is standard normal. With a freshly sampled velocity $v\sim\mathcal{N}(0,I_m)$, its exact Hamiltonian trajectory is

$$
x(t)=v\sin(t)+x\cos(t),
\qquad
v(t)=v\cos(t)-x\sin(t).
$$

For constraint $j$, define

$$
\alpha_j=f_j^Tv,
\qquad
\beta_j=f_j^Tx,
\qquad
\rho_j=\sqrt{\alpha_j^2+\beta_j^2},
\qquad
\phi_j=\operatorname{atan2}(-\alpha_j,\beta_j).
$$

The wall-hit equation is

$$
f_j^Tx(t)+g_j
=
\rho_j\cos(t+\phi_j)+g_j
=0.
$$

When $\rho_j>|g_j|$, its candidate hit times over one period are

$$
t_{j,\pm}
=
\left[
\pm\arccos\left(-\frac{g_j}{\rho_j}\right)-\phi_j
\right]_{2\pi},
$$

where $[\cdot]_{2\pi}$ maps a time into $[0,2\pi)$. Select the smallest positive candidate at which the particle is leaving the feasible half-space,

$$
\frac{d}{dt}\left(f_j^Tx(t)+g_j\right)
=
\alpha_j\cos(t)-\beta_j\sin(t)
<0.
$$

At a collision with wall $j$, reflect the velocity in the wall normal:

$$
v \gets v-2\frac{f_j^Tv}{f_j^Tf_j}f_j.
$$

The implementation uses $T=\pi/2$ and $B=30$. A feasible initial point in the original utility coordinates is

$$
z^{(0)}=
\begin{cases}
-1_m, & y=0,\\
2e_y-1_m, & y\in\{1,\ldots,m\},
\end{cases}
$$

which is transformed to $x^{(0)}=R(z^{(0)}-\mu_0)$.

```pseudo
\begin{algorithm}
\caption{HMCMNP: Exact HMC moment estimation for reference-identified multinomial probit}
\begin{algorithmic}
	\Input TMVN mean $\mu_0$, precision $P_0$, observed choice $y\in\{0,\ldots,m\}$, retained sample count $L$
	\Output Monte Carlo mean $\hat{\mu}$ and covariance $\hat{\Sigma}$

	\State $R \gets \operatorname{chol}(P_0)$ such that $P_0=R^TR$
	\For{$j=1:m$}
		\State Construct $a_j(y)$ from the reference-identified choice constraints
		\State $f_j \gets R^{-T}a_j(y)$
		\State $g_j \gets a_j(y)^T\mu_0$
	\EndFor

	\If{$y=0$}
		\State $z^{(0)} \gets -1_m$
	\Else
		\State $z^{(0)} \gets 2e_y-1_m$
	\EndIf
	\State $x \gets R(z^{(0)}-\mu_0)$
	\State $\hat{\mu} \gets 0_m$
	\State $C \gets 0_{m\times m}$

	\For{$s=1:(B+L)$}
		\Repeat
			\State $x_{\mathrm{start}} \gets x$
			\State $v \sim \mathcal{N}(0,I_m)$
			\State $t_{\mathrm{rem}} \gets T$

			\While{$t_{\mathrm{rem}}>0$}
				\For{$j=1:m$}
					\State Compute the smallest positive outgoing wall-hit time $t_j$
				\EndFor
				\State $k \gets \arg\min_j t_j$
				\If{$t_k<t_{\mathrm{rem}}$}
					\State $t_* \gets t_k$
				\Else
					\State $t_* \gets t_{\mathrm{rem}}$
				\EndIf
				\State $(x,v) \gets \left(v\sin(t_*)+x\cos(t_*),\ v\cos(t_*)-x\sin(t_*)\right)$
				\State $t_{\mathrm{rem}} \gets t_{\mathrm{rem}}-t_*$

				\If{$t_k=t_*$ and $t_{\mathrm{rem}}>0$}
					\State $v \gets v-2(f_k^Tv)(f_k^Tf_k)^{-1}f_k$
				\Else
					\State \textbf{break}
				\EndIf
			\EndWhile
			\If{any $f_j^Tx+g_j<0$}
				\State $x \gets x_{\mathrm{start}}$
			\EndIf
		\Until{$f_j^Tx+g_j\geq0$ for every $j$}
		\Comment{On numerical failure, restore the previous state and resample $v$}

		\If{$s>B$}
			\State $\ell \gets s-B$
			\State $z^{(\ell)} \gets \mu_0+R^{-1}x$
			\State $\delta \gets z^{(\ell)}-\hat{\mu}$
			\State $\hat{\mu} \gets \hat{\mu}+\delta/\ell$
			\State $C \gets C+\delta(z^{(\ell)}-\hat{\mu})^T$
		\EndIf
	\EndFor

	\State $\hat{\Sigma} \gets C/(L-1)$
	\State \Return $(\hat{\mu},\hat{\Sigma})$
\end{algorithmic}
\end{algorithm}
```

The online covariance update is the multivariate Welford recurrence. It gives the ordinary sample covariance of the retained Markov chain states without storing the full $L\times m$ sample matrix. Because the samples form a Markov chain, $L$ is the number of retained draws rather than the effective sample size.

Reference: Pakman, A. and Paninski, L. (2014). Exact Hamiltonian Monte Carlo for Truncated Multivariate Gaussians. *Journal of Computational and Graphical Statistics*, 23(2), 518--542. https://doi.org/10.1080/10618600.2013.788448

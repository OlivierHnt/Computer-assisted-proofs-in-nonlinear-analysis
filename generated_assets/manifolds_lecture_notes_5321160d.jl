### A Pluto.jl notebook ###
# v0.19.46

#> [frontmatter]
#> chapter = 4
#> order = 1
#> title = "Lecture notes: invariant manifolds of equilibria"
#> tags = ["module4"]
#> layout = "layout.jlhtml"

using Markdown
using InteractiveUtils

# ╔═╡ d0e623ee-b096-4a27-977d-dc32567d6020
using PlutoTeachingTools, PlutoUI

# ╔═╡ 018ecc45-8638-4a59-b561-efb086bdc751
using RadiiPolynomial, LinearAlgebra

# ╔═╡ 7fc40507-eda3-474d-a454-04e9173a7adb
html"""<style>
main {
    max-width: 1000px;
    margin-left: auto;
    margin-right: auto;
    text-align: justify;
}
"""

# ╔═╡ 584fad48-cadf-4c57-a4c8-4de54f53ebbe
md"""
# Local stable manifolds of equilibria
"""

# ╔═╡ d5a510a3-c518-47ed-96bc-7bb22e3b08b5
md"""
This part of the tutorial is concerned with the computation of local stable manifolds (the computation of local unstable manifolds is done similarly). Our focus is the development of efficient high order approximation methods --with validated error bounds--
for these manifolds.

To be more specific, consider a differential equation $\dot{x} = f(x)$ with associated flow $\varphi\colon \mathbb{R} \times \mathbb{R}^n \to \mathbb{R}^n$ and assume that $\tilde{x}$ is an equilibrium.
The *stable set* for $\tilde{x}$ is defined to be

$\begin{align}
W^s(\tilde{x}) = \{x\in\mathbb{R}^n : \lim_{t\to \infty} \varphi(t,x) = \tilde{x} \}.
\end{align}$

At this level of generality there is very little that can be said about $W^s(\tilde{x})$, however a natural step is to simplify the problem by restricting the dynamics to a neighborhood of $\tilde{x}$. Hence, let us define a more localized notion of a stable set. Given a neighborhood $U$ of $\tilde{x}$ the associated *local stable manifold* is given by

$\begin{align}
W^s_{\text{loc}}(\tilde{x}) = W^s(\tilde{x},U) := \{x \in W^s(\tilde{x}) \mid \varphi([0,\infty),x) \subset U \}.
\end{align}$"""

# ╔═╡ f73b88cb-19e2-4d50-a45e-fbe90fb691cd
md"""
## The Parameterization Method
"""

# ╔═╡ cf717e9b-26d8-449d-9f3b-d1d5d0bd2bb5
md"""
### Assumptions
"""

# ╔═╡ c33dc650-3f94-11ef-398a-8bbc4a2b69b8
md"""

**For sake of simplicity of the presentation**, we make the following assumptions:

-  $\tilde{x}$ is an equilibrium solution of $\dot{x}=f(x)$
-  there are exactly $m$ stable eigenvalues $\{\lambda_1, \ldots, \lambda_m \}$, they are real and there exist $m$ linearly independent associated eigenvector $\xi_i \in \mathbb{R}^n$

Denote
$\begin{align}
\Lambda =
\begin{pmatrix}
 \lambda_1 &  \ldots & 0  \\
 \vdots & \ddots  & \vdots  \\
 0  & \ldots  & \lambda_m
\end{pmatrix} \in M_m(\mathbb{R})
\quad \text{and} \quad A_0 = [\xi_1 | \ldots | \xi_m]
\end{align}$
an $n \times m$ matrix whose columns are the associated (stable) eigenvectors.

The linearized equation restricted to the stable subspace is $\dot{\theta} = \Lambda \theta$ ($\theta \in \mathbb{R}^m$) with associated flow $e^{\Lambda t}$.
"""

# ╔═╡ e0c480c1-ab86-4af4-9027-75ebb4fba960
md"""
### The conjugacy relation
"""

# ╔═╡ 1ea9c840-369e-48e6-a87b-7e4bdaa1abb0
md"""
Denote by $B^m := \{ \theta = (\theta_1,\dots,\theta_m) : |\theta_j| \le 1 \}$ the unit closed $\ell^\infty$ ball in $\mathbb{R}^m$.

To solve for the parameterization $P:B^m \to \mathbb{R}^n$ of the local stable manifold, we impose

(1) $P(0) = \tilde{x}$ (equilibria are mapped to one another)

(2) $DP(0) = A_0$ ($P$ is tangent to the stable eigenspace $\mathbb{E}^s$ at $0$})

(3) $\varphi (t, P(\theta)) = P( e^{\Lambda t} \theta)$, for all $\theta \in B^m$ (**conjugacy relation**)
"""

# ╔═╡ 1bfcc76a-266a-4784-b76b-3871445e603e
Foldable("""Why such parameterization provides a local stable manifold?""", md"""
We get that $P(B^m) \subset W^s_{\text{loc}}(\tilde{x})$, since for $\theta \in B^m$,

$\begin{align}
\lim_{t \to \infty} \varphi (t, P(\theta) ) = \lim_{t \to \infty}
	P( e^{\Lambda t} \theta)
			= P\left(\lim_{t \to \infty} e^{\Lambda t} \theta \right)
		    = P(0)
		    = \tilde{x}.
\end{align}$
""")

# ╔═╡ 3e3a3a60-f76f-401d-b265-a99e713d98c5
md"""
### The parameterization lemma
"""

# ╔═╡ 8f3e96dc-0cdc-411e-8f94-e9e47488abf3
md"""

*Remark.* The utility of (3) above is limited by the appearance of the flow $\varphi$, which is only known implicitly. Let us introduce a more practical infinitesimal version of the conjugacy relation, which requires the following important result."""

# ╔═╡ 56bbfe4c-1f28-47e5-9a16-e521f1d55350
Markdown.MD(Markdown.Admonition("tip", "Parameterization Lemma",[md"""
Let $P \colon B^m  \to \mathbb{R}^n$ be a function satisfying (1) and (2). Then $P$ satisfies (3) $\Longleftrightarrow$ $P$ solves the *partial* differential equation

$\begin{align}
\text{(4)} \quad DP(\theta) \Lambda \theta = f(P(\theta)), \quad \theta \in B^m.
\end{align}$"""]))

# ╔═╡ be501ea1-e0d3-4a4b-8f9e-19135d21aa9d
Foldable("""Proof""", md"""
($\Longleftarrow$) Assume that the assumption (4) above holds and let $\theta \in B^m$. Let $\gamma(t) := P(e^{\Lambda t} \theta)$. Then, $\gamma(0) = P(\theta)$ and

$\begin{align}
\gamma'(t) = \frac{d}{dt}  P\left( e^{\Lambda t} \theta \right)
   	= D P\left( e^{\Lambda t} \theta \right) \Lambda e^{\Lambda t} \theta
 	= f\left( P\left( e^{\Lambda t} \theta \right)\right)
 	= f ( \gamma(t)),
\end{align}$

that is $\gamma$ solves the IVP $x'=f(x), x(0) = P(\theta)$. By uniqueness, $\varphi( t, P(\theta)) = \gamma(t) = P(e^{\Lambda t} \theta)$.

($\Longrightarrow$) Suppose that $P$ solves (3). Differentiate on both sides with respect to $t$

$\begin{align}
f(\varphi(t,P(\theta))) = DP( e^{\Lambda t} \theta) \Lambda e^{\Lambda t} \theta, \quad \text{and take } t=0. \quad \square
\end{align}$
""")

# ╔═╡ b5d97ade-bcde-4c7a-ae78-697e73c10a0c
md"""
From the Parameterization Lemma, to solve for the parameterization $P:B^m \to \mathbb{R}^n$, we impose

(1) $P(0) = \tilde{x}$ (equilibria are mapped to one another)

(2) $DP(0) = A_0$ ($P$ is tangent to the stable eigenspace $\mathbb{E}^s$ at $0$)

(4) $DP(\theta) \Lambda \theta = f(P(\theta))$, $\forall~\theta \in B^m$ (invariance equation)

The invariance equation can be written more explicitly, for $\theta = (\theta_1, \ldots, \theta_m) \in B^m$, as

$\begin{align}
(4) \quad \lambda_1 \theta_1 \frac{\partial}{\partial \theta_1} P(\theta_1, \ldots, \theta_m)
+ \ldots + \lambda_m \theta_m \frac{\partial}{\partial \theta_m} P(\theta_1, \ldots, \theta_m)
= f(P(\theta_1, \ldots, \theta_m))
\end{align}$
"""

# ╔═╡ 3ebed1fa-923d-4a0e-9d80-f2919972b1ff
md"""
### Taylor series representation
"""

# ╔═╡ 62868ea5-3018-4e1a-a71e-89ede6e443c2
md"""
We use a Taylor series representation of the form

$\begin{align}
P(\theta) = \sum_{|\alpha|=0}^\infty \begin{pmatrix} (a_1)_\alpha \\ (a_2)_\alpha \\ \vdots \\ (a_n)_\alpha \end{pmatrix}
\theta^\alpha
\end{align}$

where $\alpha = (\alpha_1,\dots,\alpha_m) \in \mathbb{N}^m$, $|\alpha| = \alpha_1 + \dots + \alpha_m$, $\theta=(\theta_1,\dots,\theta_m)$ and $\theta^\alpha = \theta_1^{\alpha_1} \cdots \theta_m^{\alpha_m}$.

Denote $a = (a_1,\dots,a_n)$, where each $a_i = ((a_i)_\alpha)_{|\alpha| \ge 0}$. Plugging the above Taylor series in equations (1), (2) and (4):

(1*) $\begin{pmatrix} (a_1)_{(0,\dots,0)} \\ \vdots \\ (a_n)_{(0,\dots,0)} \end{pmatrix} = \tilde{x}$, (for $|\alpha|=0$: the equilibrium solution)

(2*) $\begin{pmatrix} (a_1)_{e_j} \\ \vdots \\ (a_n)_{e_j} \end{pmatrix} = \xi_j$, $j=1,\dots,m$
 (for $|\alpha|=1$: the $m$ stable eigenvectors)

(4*) $(\alpha \cdot \lambda) (a_i)_{\alpha} = (\phi_i(a))_{\alpha}$, $i=1,\dots,n$ (for $|\alpha| \ge 2$)

where

$\begin{align}
\alpha \cdot \lambda = \alpha_1 \lambda_1 +\cdots + \alpha_m \lambda_m \quad \text{and} \quad f_i(P(\theta))= \sum_{|\alpha| \ge 0} (\phi_i(a))_\alpha \theta^\alpha
\end{align}$
"""

# ╔═╡ f5a85447-f9a0-4f56-b2bd-49da4d74326c
md"""
### The Zero Finding Problem, Banach Space and Norm
"""

# ╔═╡ 31104635-eb81-42b0-b412-b26a3635d25f
md"""
For each $i=1,\dots,n$, define $F_i(a)=\{ (F_i(a))_\alpha \}_{|\alpha| \ge 0}$ component-wise by

$\begin{align}
(F_i(a))_\alpha := \begin{cases}
(a_i)_{(0,\dots,0)}  - \tilde{x}_i,&  |\alpha| = 0, \\
(a_i)_{e_j}  - (\xi_{j})_i, & |\alpha| = 1, {\rm~with~} \alpha=e_j {\rm~for~}j=1,\dots,m,  \\
(a_i)_{\alpha} - \frac{1}{\alpha \cdot \lambda}(\phi_i(a))_{\alpha}, & |\alpha| \ge 2.
\end{cases}
\end{align}$

For each $i=1,\dots,n$, let $c_i = \{ (c_i)_\alpha\}_{|\alpha| \ge 0}$ component-wise as

$\begin{align}
(c_i)_\alpha =
\begin{cases}
\tilde{x}_i, &  {\rm~if~} |\alpha| = 0, \\
\xi_{j}, &  {\rm~for~} |\alpha| = 1 {\rm~with~} \alpha=e_j {\rm~for~} j=1,\dots,m,  \\
0, & {\rm~if~} |\alpha| \ge 2
\end{cases}
\end{align}$

and let $c := (c_1,\dots,c_n)$. Define the diagonal operator $K$ acting on a sequence $b = \{ b_\alpha\}_{|\alpha| \ge 0}$ as

$\begin{align}
(K b)_\alpha =
\begin{cases}
0, &   {\rm~if~} |\alpha| = 0, \\
0, &   {\rm~if~} |\alpha| = 1,  \\
\frac{1}{\alpha \cdot \lambda} b_\alpha, &  {\rm~if~} |\alpha| \ge 2.
\end{cases}
\end{align}$

Define the operator $\mathcal{K}$ acting on $a=(a_1,\dots,a_n)$ as
$\begin{align}
\mathcal{K} a = (K a_1,K a_2,\dots,K a_n).
\end{align}$

With the above notation, for each $i=1,\dots,n$,

$\begin{align}
F_i(a) = a_i - K \phi_i(a) - c_i.
\end{align}$

Letting $F(a) = (F_1(a),\dots,F_n(a))$ and $\phi(a) = (\phi_1(a),\dots,\phi_n(a))$, one can densely write

$\begin{align}
F(a) = a - \mathcal{K} \phi(a) - c
\end{align}$

Define the Banach space $\mathcal{X} = (\ell^1)^n$, where

$\begin{align}
\ell^1 := \left\{ c = (c_\alpha)_{|\alpha| \ge 0} : c_\alpha \in \mathbb{R}, \quad \|c\|_1 := \sum_{|\alpha|=0}^\infty |c_\alpha| < \infty  \right\}.
\end{align}$

The norm on $\mathcal{X}$ is given by

$\begin{align}
\|a\|_{\mathcal{X}} = \sum_{i=1}^n \| a_i \|_1
\end{align}$
"""

# ╔═╡ 96bb7a72-c10b-45d7-aa64-4b61e488249a
md"""
###### Code that implements the operator $K$
"""

# ╔═╡ 3752c9c6-dd41-4194-aaf2-1cd1c2160665
function K(λ, N)
    op = zeros(typeof(λ), Taylor(N), Taylor(N))
    for n = 2:N
        op[n,n] = inv(λ * n)
    end
    return op
end

# ╔═╡ dbaf57aa-8a89-476e-b21f-56fe3d8696db
K(-2.0,5)

# ╔═╡ 97729076-0443-4740-9b0a-6e64cc75b567
Markdown.MD(Markdown.Admonition("note", "Exercise (Compact operator)",[md"""
Show that the operator $\mathcal{K}:\mathcal{X} \to \mathcal{X}$ is compact.
"""]))

# ╔═╡ 1c12da5f-7c44-49c2-8077-d93d9de9f32e
md"""
We now formalize how a solution to the above zero finding problem leads to local stable manifolds of equilibria of our differential equation.
"""

# ╔═╡ fc8635b9-2867-4580-a9ca-fd2383c0efd1
Markdown.MD(Markdown.Admonition("lemma", "Lemma (A solution to the zero finding problem yields a parameterization of the local stable manifold)",[md""" Assume that
-  $\tilde{x}$ is an equilibrium solution of $\dot{x}=f(x)$;
-  there are exactly $m$ stable eigenvalues $\{\lambda_1, \ldots, \lambda_m \}$, they are real and there exist $m$ linearly independent associated eigenvector $\xi_i \in \mathbb{R}^n$.

If $a \in \mathcal{X}$ is a solution to $F(a)=0$, then the corresponding Taylor series expansion

$\begin{align}
P(\theta) = \sum_{|\alpha|=0}^\infty \begin{pmatrix} (a_1)_\alpha \\ (a_2)_\alpha \\ \vdots \\ (a_n)_\alpha \end{pmatrix}
\theta^\alpha
\end{align}$

converges on $B^m := \{ \theta = (\theta_1,\dots,\theta_m) : |\theta_j| \le 1 \}$ and the function $P:B^m \to \mathbb{R}^n$ provides a parameterization of the local stable manifold $W^s_{\text{loc}}(\tilde{x})$, that is

$\begin{align}
P(B^m) = W^s_{\text{loc}}(\tilde{x}).
\end{align}$
"""]))

# ╔═╡ a1f8e3df-15ff-4c95-b88b-717901b67c1d
md"""
Rather than presenting the bounds for the Newton-Kantorovich Theorem in general, let us rather focus on an example. Specifically, we compute a rigorous parameterization of a one-dimensional local stable manifold of an equilibrium solution in the Lorenz system.
"""

# ╔═╡ d9604e75-0e9e-44b6-bf2e-885f01a96bdc
md"""
## A one-dimensional local stable manifold in the Lorenz system
"""

# ╔═╡ 62719aa9-bfda-4913-8b74-f3694bfc33c0
md"""
The Lorenz system is given by

$\begin{align}
\dot x = f(x) := \begin{pmatrix}
\sigma (x_2 - x_1) \\
\rho x_1 - x_2 - x_1 x_3 \\
- \beta x_3 + x_1 x_2
\end{pmatrix}
\end{align}$

where $\sigma,\rho, \beta>0$ are parameters.

It is straightforward (**Exercise**) to verify that if $\rho>1$, there are two nontrivial equilibria given by
$\left\{
\begin{pmatrix} -\sqrt{\beta(\rho-1)} \\ -\sqrt{\beta(\rho-1)} \\ \rho-1 \end{pmatrix},
\quad
\begin{pmatrix} \sqrt{\beta(\rho-1)} \\ \sqrt{\beta(\rho-1)} \\ \rho-1 \end{pmatrix} \right\}$.
"""

# ╔═╡ fc381b29-6671-4364-a9df-b35db07af45b
md"""
We make the following **assumptions**:
-  Assume that $\rho>1$, and let us choose the equilibrium solution to be
$\begin{align}
\tilde x := \begin{pmatrix} \sqrt{\beta(\rho-1)} \\ \sqrt{\beta(\rho-1)} \\ \rho-1 \end{pmatrix}.
\end{align}$
-  Assume that there is a unique stable real eigenvalue $\lambda<0$ with associated eigenvector $\xi \in \mathbb{R}^3$. Hence, the dimension of $W^s_{\text{loc}}(\tilde{x})$ is one, that is $m=1$.
"""

# ╔═╡ 194ee83b-7442-42f6-862a-2f5bd1c8b710
Markdown.MD(Markdown.Admonition("warning", "Remark",[md""" For $\rho>24.75$, $\sigma=10$ and $\beta=8/3$, $Df(\tilde{x})$ has one real negative eigenvalues and two complex conjugate eigenvalues with positive real parts (**Exercise**). In this case, the local stable manifold $W^s_{\text{loc}}(\tilde{x})$ is one-dimensional.
"""]))

# ╔═╡ 9829987c-2bf5-4b33-8ae1-6b80d43ed978
md"""
Under the above assumptions, $m=1$, and the map $F$ is given by

$\begin{align}
F(a) = a - \mathcal{K} \phi(a) - c, \quad \quad \quad {\rm where } \quad \phi(a):=\begin{pmatrix}
\sigma (a_2 - a_1) \\
\rho a_1 - a_2 - a_1*a_3 \\
- \beta a_3 + a_1*a_2
\end{pmatrix}
\end{align}$

where $a_i*a_j$ is the **Cauchy product** defined by

$\begin{align}
(a_i*a_j)_n = \sum_{k=0}^n (a_i)_k (a_j)_{n-k}
\end{align}$

The following result will be useful when performing the nonlinear analysis.
"""

# ╔═╡ 1f8fba92-d023-4444-b4e3-18c751f678c5
Markdown.MD(Markdown.Admonition("note", "Exercise (Banach algebra)",[md"""
Prove that for all $a_1,a_2 \in \ell^1$, $\|a_1*a_2\|_1 \le \|a_1\|_1 \|a_2\|_1$.
"""]))

# ╔═╡ 11aa1d55-1a2c-48cf-a7f5-30ff16e2e92a
md"""
Note that in the case of a one-dimensional manifold, the sequences are indexed over a single index running over $\mathbb{N}$, as opposed to the more general multi-index case of an $m$-dimensional stable manifold. In the 1D case,

$\begin{align}
\ell^1 = \left\{ c = (c_n)_{n \ge 0} : c_n \in \mathbb{R}, \quad \|c\|_1 := \sum_{n=0}^\infty |c_n| < \infty  \right\}
\end{align}$
"""

# ╔═╡ 668058ef-e66b-4fd4-ad67-2a8389589e86
Markdown.MD(Markdown.Admonition("note", "Exercise",[md""" Show that for all $a = (a_1,a_2,a_3) \in \mathcal{X} := (\ell^1)^3$, $F(a) \in \mathcal{X}$. That shows that $F:\mathcal{X} \to \mathcal{X}$.
"""]))

# ╔═╡ 0562b884-0653-4d5b-8a66-e35467e6693e
md"""
#### Evaluating the operator norm
"""

# ╔═╡ a47a4e2a-371e-4244-8574-493f00d456df
Markdown.MD(Markdown.Admonition("note", "Exercise (characterization of operator norms on product of sequence spaces)",[md"""

Consider the Banach space $\mathcal{X}=(\ell^1)^n$ with norm

$\begin{align}
\|a\|_{\mathcal{X}} = \sum_{i=1}^n \| a_i \|_1.
\end{align}$

Let $\Gamma \in B(\mathcal{X})$ be a bounded linear operator acting on $\mathcal{X} = (\ell^1)^n$ and represented by

$\begin{align}
\Gamma = \begin{pmatrix}
\Gamma_{1,1} & \cdots & \Gamma_{1,n} \\
\vdots & \ddots & \vdots \\
\Gamma_{n,1} & \cdots & \Gamma_{n,n}
\end{pmatrix}
\end{align}$

where each $\Gamma_{i,j}$ is a bounded linear operator on $\ell^1$. More explicitly, $\Gamma$ acts on an element $a = (a_1,\dots,a_n) \in \mathcal{X}$ as

$\begin{align}
\Gamma a = \left( (\Gamma a)_1,\dots,(\Gamma a)_n \right) = \left( \sum_{j=1}^n \Gamma_{1,j} a_j , \dots, \sum_{j=1}^n \Gamma_{n,j} a_j \right),
\end{align}$

where$(\Gamma_{i,j} a_j)_{n_1} = \sum_{n_2 \ge 0 } (\Gamma_{i,j})_{n_1,n_2} (a_j)_{n_2}$. Then, show that

$\begin{align}
\boxed{
  \|\Gamma\|_{B(\mathcal{X})} = \max_{j=1,\dots,n} \sum_{i=1}^n \| \Gamma_{i,j} \|_{B(\ell^1)}
}
\end{align}$

where

$\begin{align}
  \|\Gamma_{i,j} \|_{B(\ell^1)} = \sup_{n_2 \ge 0} \sum_{n_1 \ge 0} |(\Gamma_{i,j})_{n_1,n_2}| .
\end{align}$
"""]))

# ╔═╡ beca7238-ca9b-409a-8d5a-240ca720f58a
md"""
Here is an example of how to compute the operator norm of $\pi^{\le 5}\mathcal{K}$
"""

# ╔═╡ 5e6482db-9b76-4bf0-bc58-0a6f779c4556
begin
	mat = K(-2.0, 5)
	M = zeros(Taylor(5)^3, Taylor(5)^3) # The operator K
	component(M, 1, 1) .= mat
	component(M, 2, 2) .= mat
	component(M, 3, 3) .= mat
	ell1 = Ell1() #no entry means that the weight is 1
	X = NormedCartesianSpace((ell1, ell1, ell1), ell1) # each space has l1 and the outer norm is l1
	opnorm(M, X) # this is equivalent to opnorm(M, 1)
end

# ╔═╡ 8c0ac2a2-b3db-455a-8a6d-3de6728b61d3
Markdown.MD(Markdown.Admonition("note", "Exercise",[md"""
Show theoretically that $\|\mathcal{K}\|_{B(\mathcal{X})} \le \frac{1}{2|\lambda|}$
"""]))

# ╔═╡ 7cfbe32b-d677-4611-97d9-ed775c60e48e
md"""
#### Finite dimensional projection
"""

# ╔═╡ 13f3b4e4-0409-49a2-a544-0c34ff9289db
md"""
Computing an approximate solution to $F=0$ requires first considering a finite dimensional projection. We need the projection operator on a finite number of Taylor coefficients ($N+1$ coefficients) and its complement:

$\begin{align}
(\pi^{\leq N} b)_n := \begin{cases}
b_n, & n=0,\dots,N, \\
0, & n > N,
\end{cases}
& \qquad\qquad
(\pi^{>N} b)_n := \begin{cases}
0, & n=0,\dots,N, \\
b_n, & n > N.
\end{cases}
\end{align}$

Now, given three sequences $a_1,a_2,a_3$ and for $a=(a_1,a_2,a_3)$, consider the finite dimensional projection on a finite number of Taylor coefficients ($3(N+1)$ coefficients) and its complement:

$\begin{align}
\pi^{\leq N} a = \left( \pi^{\leq N} a_1,\pi^{\leq N} a_2,\pi^{\leq N} a_3 \right) \qquad {\rm and} \qquad
\pi^{> N} a = \left( \pi^{> N} a_1,\pi^{> N} a_2,\pi^{> N} a_3 \right)
\end{align}$

so that $a=\pi^{\leq N}a+\pi^{>N}a$ and $\pi^{\leq N}\pi^{>N}=0=\pi^{>N}\pi^{\leq N}$.

The range $\pi^{\leq N} \mathcal{X}$ is finite dimensional, and the restriction of $\pi^{\leq N} F$ to $\pi^{\leq N} \mathcal{X}$, which we denote by $F^{\leq N}$, is what we will work with in the computer. An implementation is given below, as well as an implementation of its derivative.
"""

# ╔═╡ 40156798-59d5-4afe-a054-29c753a51e46
md"""
## Numerics
"""

# ╔═╡ 42ad70a3-114e-4cd3-8dd1-2429ac80bfc2
md"""
##### Implementation code for the finite dimensional projection $F^{\leq N}$
"""

# ╔═╡ cdbbb2e1-56bc-4beb-93d3-2db1a2aa71d5
function F(a, x, λ, ξ, par, N)
    a₁, a₂, a₃ = eachcomponent(a)

    σ, ρ, β = par

    F₁ = project(a₁ - K(λ, N) * (σ*(a₂ - a₁)), Taylor(N))
    F₁[0] -= x[1]
    F₁[1] -= ξ[1]

    F₂ = project(a₂ - K(λ, N) * (a₁*(ρ - a₃) - a₂), Taylor(N))
    F₂[0] -= x[2]
    F₂[1] -= ξ[2]

    F₃ = project(a₃ - K(λ, N) * (a₁*a₂ - β*a₃), Taylor(N))
    F₃[0] -= x[3]
    F₃[1] -= ξ[3]

    return Sequence(Taylor(N)^3, [coefficients(F₁) ; coefficients(F₂) ; coefficients(F₃)])
end

# ╔═╡ aaf278f6-b8ee-4094-9173-ffde2cb0b461
md"""
##### Implementation code for the Jacobian matrix $DF^{\leq N}$
"""

# ╔═╡ d77a192e-4be0-4987-8023-8d91b22d8d9f
function DF(a, λ, par, N)
    a₁, a₂, a₃ = eachcomponent(a)

    σ, ρ, β = par

    D₁F₁ = I - K(λ, N) * (-σ)
    D₁F₁[0,:] .= 0; D₁F₁[1,:] .= 0; D₁F₁[0,0] = D₁F₁[1,1] = 1
    D₂F₁ =   - K(λ, N) * σ
    D₃F₁ = zeros(Taylor(N), Taylor(N))

    D₁F₂ =   - K(λ, N) * project(Multiplication(ρ - a₃), Taylor(N), Taylor(N))
    D₂F₂ = I - K(λ, N) * (-1)
    D₂F₂[0,:] .= 0; D₂F₂[1,:] .= 0; D₂F₂[0,0] = D₂F₂[1,1] = 1
    D₃F₂ =   - K(λ, N) * project(Multiplication(-a₁),    Taylor(N), Taylor(N))

    D₁F₃ =   - K(λ, N) * project(Multiplication(a₂), Taylor(N), Taylor(N))
    D₂F₃ =   - K(λ, N) * project(Multiplication(a₁), Taylor(N), Taylor(N))
    D₃F₃ = I - K(λ, N) * (-β)
    D₃F₃[0,:] .= 0; D₃F₃[1,:] .= 0; D₃F₃[0,0] = D₃F₃[1,1] = 1

    return LinearOperator(Taylor(N)^3, Taylor(N)^3,
        [coefficients(D₁F₁) coefficients(D₂F₁) coefficients(D₃F₁);
         coefficients(D₁F₂) coefficients(D₂F₂) coefficients(D₃F₂);
         coefficients(D₁F₃) coefficients(D₂F₃) coefficients(D₃F₃)])
end

# ╔═╡ 7410e377-715c-4464-a83d-d7cd48c4d1ae
md"""
From now on, we fix the parameters in the Lorenz system to be the *classical* ones, namely
$\begin{align}
\rho=28, \quad \sigma=10 \quad \text{and} \quad \beta=8/3.
\end{align}$
"""

# ╔═╡ fbeac0c7-8fcb-4e4b-9662-9f1d48b23fb3
begin
σ, ρ, β = 10.0, 28.0, 8/3
par = [σ, ρ, β]
tσ, tρ, tβ = interval(10.0), interval(28.0), interval(8)/interval(3)
tpar = [tσ, tρ, tβ]
end

# ╔═╡ ae88c826-f486-4786-9447-affcc466d1b5
md"""
For these parameters, there is a unique real negative eigenvalue $\lambda$ of $Df(\tilde{x})$ with associated eigenvector $\xi \in \mathbb{R}^3$ (see the above Remark).
"""

# ╔═╡ 76aa085c-d7ab-4e2f-9648-89b3eebe11cd
md"""
### Rigorous computation of the stable eigenvalue $\lambda$ and eigenvector $\xi \in \mathbb{R}^3$
"""

# ╔═╡ 585b99f6-fad0-44b1-b1ca-f7aa6f8f52f8
md"""
While we could work hard and compute the eigenpair $(\lambda,\xi)$ by hand, we resort to a computer-assisted proof to obtain a rigorous enclosure. Here is a code that performs this task.
"""

# ╔═╡ 400fe5e0-0dfd-4d82-baf9-b66d931133a8
md"""
##### Numerical approximation and rigorous enclosure of the equilibrium solution $\tilde{x}$
"""

# ╔═╡ e0f75d0e-555d-4df7-8077-9871a6d7b197
begin
bx = [sqrt(β*(ρ-1)), sqrt(β*(ρ-1)), ρ-1] #Numerical approximation
tx = [sqrt(tβ*(tρ-interval(1))), sqrt(tβ*(tρ-interval(1))), tρ-interval(1)] #Rigorous enclosure
end

# ╔═╡ deb48c67-19f1-481e-b2de-1765e449e6ae
md"""
##### Jacobian matrix for the Lorenz system
"""

# ╔═╡ c106822b-2c0c-4759-84fb-9ddba91fe295
function Df(u, par)
    σ, ρ, β = par
    u₁, u₂, u₃ = u
    return [-σ     σ         zero(u₃)
             ρ-u₃ -one(u₂)  -u₁
             u₂    u₁       -β]
end

# ╔═╡ f75c26b4-8e6f-4520-b1a8-f30c20b9a194
md"""
##### Computation of an approximate eigenpair
"""

# ╔═╡ 05a4dd64-b5c9-4014-a18c-f26938e02ebd
begin
approx_eigenvalues, approx_eigenvectors = eigen(Df(bx, par))
bλ, bξ = real(approx_eigenvalues[1]), real(approx_eigenvectors[:,1])
end

# ╔═╡ 3cb601bf-43ea-430c-950d-402b0929e515
md"""
##### Zero finding problem for an eigenpair (see Part I of the tutorial)
"""

# ╔═╡ 2a76f46f-c203-4028-80ba-e7bc1f481076
function F_eig(X, M, ū)
	λ = X[1]
	u = X[2:end]
	return Sequence(
		[sum(u .* conj(ū)) - 1;
		 M*u - λ*u])
end

# ╔═╡ 8a802d5f-e2f5-4e03-b0d9-b7cd0c081d2c
md"""
##### The Jacobian of the zero finding map for an eigenpair
"""

# ╔═╡ 5ef06ffe-d2f3-4cd2-946a-df8fd78f863f
function DF_eig(X, M, ū)
	λ = X[1]
	u = X[2:end]
	n = length(X) - 1
	temp = zeros(eltype(X), n+1, n+1)
	temp[2:end,1] = -u
	temp[1,2:end] = ū'
	temp[2:end,2:end] = M - λ*I
	return LinearOperator(temp)
end

# ╔═╡ e4a015b0-e7e6-4841-b97f-4a845f462196
md"""
##### A Newton-Kantorovich Theorem to obtain a computer-assisted proof for an eigenpair
"""

# ╔═╡ c3d0d2ae-5711-451b-946d-6cdca682fe06
function validate_eigenpair(λ, u, M)
	X = [λ ; u]
	ū = u
	iA = interval.(inv(DF_eig(X, mid.(M), ū)))
	iX = interval.(X)
	Y = norm(iA * F_eig(iX, M, ū), 1)
	Z₁ = opnorm(I - iA * DF_eig(iX, M, ū), 1)
	Z₂ = opnorm(iA, 1)
	return interval_of_existence(Y, Z₁, Z₂, Inf)
end

# ╔═╡ 4ff1b050-3036-4e0b-b9db-b07b4d17327c
md"""
###### The rigorous enclosure of the eigenpair
"""

# ╔═╡ ec3cd93f-877d-419e-a608-cef748306927
begin
err_eig = inf(validate_eigenpair(bλ, bξ, Df(tx, tpar)))
tλ = interval(bλ, err_eig; format = :midpoint)
tξ = interval.(bξ, err_eig; format = :midpoint)
end

# ╔═╡ 6e715b86-4604-49c3-8df1-929921e38f9a
md"""
##### We introduce a rescaling for the eigenvector (will impact the decay of the Taylor coefficients)
"""

# ╔═╡ ffe5218b-88ca-42e0-8a0d-074ab2e6d898
begin
eigenvector_scaling = 60;
bξ_scaled = eigenvector_scaling*bξ
tξ_scaled = eigenvector_scaling*tξ
end

# ╔═╡ 98c9ae61-0a84-48c0-8773-0fe9bbd260ac
md"""
#### Compute an numerical approximation of the Taylor coefficients with Newton's method
"""

# ╔═╡ 5838bc82-6896-49a7-bc2c-b700fb3a3eb8
begin

N = 150 #Choose a Taylor order

initialguess = zeros(Taylor(N)^3) # The initial guess consists of the linear approximation of the manifold
for i = 1:3
    component(initialguess, i)[0] = bx[1]
    component(initialguess, i)[1] = bξ[1]
end

a, _ = newton(a -> (F(a, bx, bλ, bξ_scaled, par, N), DF(a, bλ, par, N)), initialguess)

end

# ╔═╡ 3c23bf6e-fd79-4770-83c3-385832eddac6
begin
	using OrdinaryDiffEq
	using GLMakie
	function lorenz(du, u, p, t)
	    du[1] = 10.0(u[2] - u[1])
	    du[2] = u[1] * (28.0 - u[3]) - u[2]
	    du[3] = u[1] * u[2] - (8 / 3) * u[3]
	end
	u0 = [0.0; 10.0; 10.0]
	tspan = (0.0, 50.0)
	prob = ODEProblem(lorenz, u0, tspan)
	sol = solve(prob, Tsit5())
	fig = Figure()
	ax = Axis3(fig[1,1]; aspect = :data, azimuth = 3π/5)
	lines!(ax, [Point3f(sol(t)) for t = LinRange(tspan..., 10_001)]; color = :royalblue1)
	lines!(ax, [Point3f(a(t)) for t = LinRange(-1, 1, 101)]; color = :salmon, linewidth = 3)
	meshscatter!(ax, [Point3f(bx), Point3f(-bx[1], -bx[2], bx[3]), Point3f(0,0,0)]; markersize = 1, color = :palegreen3)
	fig
end

# ╔═╡ 0bff34a8-2ad3-43d8-9f25-60a926f5b5c6
md"""
##### Plot the image of the parameterization
"""

# ╔═╡ c4363fdf-6500-4830-bc42-9777d616ac34
md"""
### Choosing the approximate inverse $A$ for the Fréchet derivative $DF(\bar{a})$
"""

# ╔═╡ e72d43a1-b804-455a-bff0-8b1b715ddcee
md"""
The derivative of $F(a)=a-\mathcal{K}\phi(a)-c$ can now be written compactly as

$\begin{align}
DF(a) = I - \mathcal{K} D\phi(a)
\end{align}$

acting on an element $b=(b_1,b_2,b_3) \in \mathcal{X}$ as

$\begin{align}
DF(a)b = b  - \mathcal{K} D\phi(a) b.
\end{align}$
"""

# ╔═╡ 24e1a5af-569e-4ae8-a3d7-19b7962db48b
md"""
We use the projection to split $A$ into a finite part and a tail part:

$\begin{align}A := A^{\le N}\pi^{\le N} + I^{>N}\pi^{>N},\end{align}$

where $A^{\le N}$ is a linear map (to be chosen below) on $\pi^{\le N} \mathcal{X}$ and $I^{>N}$ is the identity on $\pi^{> N} \mathcal{X}$. Hence by construction one may also write $A$ as:

$\begin{align}A = \pi^{\le N}A^{\le N}\pi^{\le N} + \pi^{>N}I^{>N}\pi^{>N}.\end{align}$
"""

# ╔═╡ 922f5fa4-f9db-4eae-bd81-dc3dd525d212
md"""
We compute a numerical inverse of the Jacobian $DF^{\le N}(\bar{a})$ and we denote this inverse matrix by $A^{\le N}$.
Moreover, we compute the operator norm of the operator $A$.
"""

# ╔═╡ bb5ea822-f29b-499c-bebb-53208cf42eb1
begin
A_finite = interval.( inv(DF(a, bλ, par, N)) )
opnorm_A_finite = opnorm(A_finite, 1)
opnormA = max(opnorm_A_finite, 1)
end

# ╔═╡ 0d75d077-53c3-4f03-99d3-af796d973d66
md"""
## The Newton-Kantorovich proof
"""

# ╔═╡ 9433b6f1-3460-48c4-b587-9009823327c0
md"""
We are now in the Newton-Kantorovich context and thus need to derive the appropriate bounds.
"""

# ╔═╡ 4741e38b-0c29-4f02-a540-c543f7b3f234
md"""
The following code is used to compute the operator $\pi^{> N}$ acting on a single sequence or a list of sequences
"""

# ╔═╡ 656ae53a-1dbc-4abc-8837-234762db80c0
function tail(x, N)
    # get the tail of a Taylor series
    y = zeros(eltype(x), space(x))
    for n = N+1:order(x)
        y[n] = x[n]
    end
    return y
end

# ╔═╡ 9fe5120a-824a-4395-afe8-325a36ab9a07
md"""
### The bound $Y$
"""

# ╔═╡ eb2ae1fd-5062-436b-ac50-eb3b169bca42
md"""
Recall that that bound $Y$ satisfies $\| A F(\bar{a}) \|_{\mathcal{X}} \le Y$. Note that $A F(\bar{a}) = \pi^{\le 2N} A F(\bar{a}) + \stackrel{=0}{\overbrace{\pi^{> 2N} A F(\bar{a})}}$, since $(\bar{a}_1*\bar{a}_2)_n=(\bar{a}_1*\bar{a}_3)_n=0$ for $n>2N$. Hence,

$\begin{align}
A F(\bar{a}) &= \pi^{\le 2N} A F(\bar{a})  \\
& = \pi^{\le N} A F(\bar{a})+ (\pi^{\le 2N} - \pi^{\le N}) A F(\bar{a}) \\
& = \pi^{\le N} A^{\le N} F^{\le N}(\bar{a}) + (\pi^{\le 2 N}-\pi^{\le N}) (\bar{a} - \mathcal{K} \phi(\bar{a})-c) \\
& = \pi^{\le N} A^{\le N} F^{\le N}(\bar{a}) - (\pi^{\le 2 N}-\pi^{\le N}) \mathcal{K} \phi(\bar{a})
\end{align}$

The code just below computes with interval arithmetic the bound $Y$ satisfying

$\begin{align}
\| A F(\bar{a}) \|_{\mathcal{X}}
\le \| \pi^{\le N} A^{\le N} F^{\le N}(\bar{a}) \|_{\mathcal{X}} + \frac{1}{|\lambda| (N+1)} \| (\pi^{\le 2 N}-\pi^{\le N}) \phi(\bar{a}) \|_{\mathcal{X}} \le Y.
\end{align}$
"""

# ╔═╡ 3021e7c0-508a-4a95-9c16-098fa35b65a0
begin
full_F = F(interval.(a), tx, tλ, tξ_scaled, par, 2N)
Y = norm(A_finite * full_F, 1) + sum(i -> norm(tail(component(full_F, i), N), 1), 1:3) / (abs(tλ)*(N+1))
end

# ╔═╡ fc18e3ff-4a39-4aba-a2d2-94521035013e
md"""
### The bound $Z_1$
"""

# ╔═╡ 40a4ef18-8cf6-46bf-bf84-a7c89d898f51
md"""
We begin with three Exercises that will be useful in deriving the $Z_1$ bound.
"""

# ╔═╡ 454efc82-4918-4c60-8b9d-fb630302beec
Markdown.MD(Markdown.Admonition("note", "Exercise 1 for Z_1",[md""" Show that $\| \pi^{>N} \mathcal{K} \|_{B(\mathcal{X})} \le \frac{1}{|\lambda|(N+1)}$.
"""]))

# ╔═╡ c52ec039-8ef2-41fc-a9cd-34fcef45f735
Markdown.MD(Markdown.Admonition("note", "Exercise 2 for Z_1",[md""" Show that
$\begin{align} \| D\phi(\bar{a}) \|_{B(\mathcal{X})}
= \max_{i=1,2,3} \sum_{j=1}^3 \| \partial_i \phi_j(\bar{a})\|_1 \end{align}
= \max \left\{
\sigma + \|\rho-a_3\|_1 + \|a_2\|_1,
\sigma + 1 + \|a_1\|_1,
\|a_1\|_1 + \beta
\right\}.$
"""]))

# ╔═╡ 7229c43e-dcda-436c-86fc-1e626fe43daf
Foldable("""Hint""", md"""
Note that

$\begin{align}
D\phi(a) = \begin{pmatrix}
-\sigma I & \sigma I & 0 \\
\rho I - \mathcal{M}_{a_3} & -I & - \mathcal{M}_{a_1} \\
\mathcal{M}_{a_2} & \mathcal{M}_{a_1} & - \beta I
\end{pmatrix}
\end{align},$

where
$\mathcal{M}_{a_i} \in B(\ell^1)$ is the multiplication operator defined by $\mathcal{M}_{a_i} b = a_i*b$ with norm $\| \mathcal{M}_{a_i} \|_{B(\ell^1)} = \| a_i\|_1$.
""")

# ╔═╡ bc448153-fe09-4a05-b857-5b9fd51dcdae
Markdown.MD(Markdown.Admonition("note", "Exercise 3 for Z_1",[md""" Show that
$\pi^{\le N} D\phi(\bar{a})\pi^{> N}=0$.
"""]))

# ╔═╡ 1d3d35a7-c7a8-419e-b751-a144b86bbe1b
Foldable("""Hint""", md"""
Note that each none blocks of $D\phi(\bar{a})$ is lower triangular since of the Cauchy products.
""")

# ╔═╡ 002c7c0a-2b38-4a31-9e08-a424fa5679b8
md"""
Recall that the bound $Z_1$ satisfies
$\begin{align}
\| I - A DF(\bar{a}) \|_{B(\mathcal{X})} \le Z_1.
\end{align}$

We use a similar splitting as in the case of periodic orbits using Fourier series (see Part 3 of the tutorial), and observe that

$\begin{align}
\|\Gamma\|_{B(\mathcal{X})} = \max \Big( \|\Gamma \pi^{\le N}\|_{B(\mathcal{X})}, \|\Gamma \pi^{> N}\|_{B(\mathcal{X})} \Big).
\end{align}$

We compute separately two bounds $M_1,M_2$ satisfying $\|(I - A DF(\bar{a})) \pi^{\le N}\|_{B(\mathcal{X})} \le M_1$ and $\| (I - A DF(\bar{a})) \pi^{> N}\|_{B(\mathcal{X})} \le M_2$.

Recall that $A = A^{\le N}\pi^{\le N} + I^{>N}\pi^{>N}$.
"""

# ╔═╡ 9f022768-29a7-4ce9-bf9a-ceb04ea0adf2
md"""
###### Computing $M_1$
"""

# ╔═╡ 090c74fc-9bf7-47c0-af5f-887d6bf43ccd
md"""
Note first that

$\begin{align}
\|(I - A DF(\bar{a})) \pi^{\le N}\|_{B(\mathcal{X})} &\le
\|(I - A^{\le N}\pi^{\le N} DF(\bar{a})) \pi^{\le N}\|_{B(\mathcal{X})}
 + \| \pi^{>N} DF(\bar{a}) \pi^{\le N}\|_{B(\mathcal{X})},
\end{align}$

where $\pi^{>N} DF(\bar{a}) \pi^{\le N} = \pi^{>N} (I-\mathcal{K} D\phi(\bar{a})) \pi^{\le N}
= \pi^{>N} \pi^{\le N} - \pi^{>N} \mathcal{K} D\phi(\bar{a}) \pi^{\le N}
= - \pi^{>N} \mathcal{K} D\phi(\bar{a}) \pi^{\le N}$, and hence by Exercises 1 and 2:

$\begin{align}
\| \pi^{>N} DF(\bar{a}) \pi^{\le N}\|_{B(\mathcal{X})} & \le
\| \pi^{>N} \mathcal{K} \|_{B(\mathcal{X})} \| \pi^{>N} D\phi(\bar{a}) \|_{B(\mathcal{X})} \\
& \le \frac{1}{|\lambda|(N+1)} \| D\phi(\bar{a}) \|_{B(\mathcal{X})} \\
& \le \frac{1}{|\lambda|(N+1)} \max \left\{
\sigma + \|\rho-a_3\|_1 + \|a_2\|_1,
\sigma + 1 + \|a_1\|_1,
\|a_1\|_1 + \beta
\right\}
\end{align}$

Using interval arithmetic, one can compute an upper bound $M_1$ such that

$\begin{align}
\|(I - A DF(\bar{a})) \pi^{\le N}\|_{B(\mathcal{X})} &\le
\|(I - A^{\le N}\pi^{\le N} DF(\bar{a})) \pi^{\le N}\|_{B(\mathcal{X})} +
\frac{1}{|\lambda|(N+1)} \max \left\{
\sigma + \|\rho-a_3\|_1 + \|a_2\|_1,
\sigma + 1 + \|a_1\|_1,
\|a_1\|_1 + \beta
\right\}\le M_1
\end{align}$
"""

# ╔═╡ 9e8d66d9-ad6e-42b5-a24f-6d0a9060f688
md"""
###### Computing $M_2$
"""

# ╔═╡ 6297e819-a105-4451-bbb7-c2241ff77f51
md"""
Now, for $M_2$:

$\begin{align}
\|(I - A DF(\bar{a})) \pi^{> N}\|_{B(\mathcal{X})} & =
\|(I - (A^{\le N}\pi^{\le N} + I^{>N}\pi^{>N}) DF(\bar{a})) \pi^{> N}\|_{B(\mathcal{X})} \\
&\le
\| A^{\le N}\stackrel{=0}{\overbrace{\pi^{\le N} DF(\bar{a})\pi^{> N}}}\|_{B(\mathcal{X})}
 + \| \pi^{> N}-\pi^{>N} DF(\bar{a}) \pi^{> N}\|_{B(\mathcal{X})},
\end{align}$

since $\pi^{\le N} DF(\bar{a})\pi^{> N} = \pi^{\le N} (I- \mathcal{K} D\phi(\bar{a}))\pi^{> N}
=\pi^{\le N}\pi^{> N} - \pi^{\le N} D\phi(\bar{a}))\pi^{> N} = - \pi^{\le N} D\phi(\bar{a}))\pi^{> N}=0$ by Exercise 3.

Now, using a similar idea as for $M_1$, one can show that

$\begin{align}
\| \pi^{> N}-\pi^{>N} DF(\bar{a}) \pi^{> N}\|_{B(\mathcal{X})} \le \frac{1}{|\lambda|(N+1)} \max \left\{
\sigma + \|\rho-a_3\|_1 + \|a_2\|_1,
\sigma + 1 + \|a_1\|_1,
\|a_1\|_1 + \beta
\right\}\end{align}$.
"""

# ╔═╡ fb8527a5-4fa9-40d3-b0af-2cce8cacac3d
md"""
###### Computing $Z_1$
"""

# ╔═╡ 7e60324b-f65e-4d08-b642-13542b3637f3
md"""
Since clearly $M_2 \le M_1$, then we let $Z_1=M_1$. Then, using interval arithmetic, one can compute an upper bound $Z_1$ such that

$\begin{align}
\| I - A DF(\bar{a}) \|_{B(\mathcal{X})}
 &\le
\|(I - A^{\le N}\pi^{\le N} DF(\bar{a})) \pi^{\le N}\|_{B(\mathcal{X})} +
\frac{1}{|\lambda|(N+1)} \max \left\{
\sigma + \|\rho-a_3\|_1 + \|a_2\|_1,
\sigma + 1 + \|a_1\|_1,
\|a_1\|_1 + \beta
\right\}\le Z_1
\end{align}$
which we do below.
"""

# ╔═╡ c5d63b4c-7823-4122-a81f-9e3a18a634ff
begin
a₁, a₂, a₃ = eachcomponent(interval.(a))

opnorm_Dϕ = max(
    tσ + norm(ρ - a₃, 1) + norm(a₂, 1),
    tσ + 1               + norm(a₁, 1),
         norm(a₁, 1)     + β)

Z₁ = opnorm(I - A_finite * DF(interval.(a), tλ, tpar, N), 1) + opnorm_Dϕ / (abs(tλ)*(N+1))
end

# ╔═╡ c0ce2c85-eaae-4b08-8ef7-2a9297b10e05
md"""
### The bound $Z_2$
"""

# ╔═╡ 6f346c67-879a-4427-afd4-6acdeb17af21
md"""
This concerns the bound on the operator $A[DF(a)-DF(\bar{a})]$. Since this is a Lipschitz bound on the first derivative, one approach is to bound the second derivative. However, here we just do a direct computation to see that
$A[DF(a)-DF(\bar{a})]b = - A \mathcal{K} [D\phi(a)-D\phi(\bar{a})]b$.

Note that

$\begin{align}
(D\phi(a) - D\phi(\bar{a}))b = \begin{pmatrix}
0 & 0 & 0 \\
0 - (\mathcal{M}_{a_3}-\mathcal{M}_{\bar{a}_3}) & 0 & - (\mathcal{M}_{a_1}-\mathcal{M}_{\bar{a}_1}) \\
(\mathcal{M}_{a_2}-\mathcal{M}_{\bar{a}_2}) & (\mathcal{M}_{a_1}-\mathcal{M}_{\bar{a}_1}) & 0
\end{pmatrix} \begin{pmatrix} b_1 \\ b_2 \\ b_3 \end{pmatrix}
= \begin{pmatrix} 0 \\ -(a_3-\bar{a}_3)*b_1 -(a_1-\bar{a}_1)*b_3 \\ (a_2-\bar{a}_2)*b_1 + (a_1-\bar{a}_1)*b_2\end{pmatrix}
\end{align},$

hence by the Banach algebra property and the definition of the operator norm we arrive at

$\begin{align}
\|A[DF(a)-DF(\bar{a})]b\|_{\mathcal{X}} &\le \|A\|_{B(\mathcal{X})} \| \mathcal{K} \|_{B(\mathcal{X})}
(
\| a_3-\bar{a}_3 \|_1 \|b_1\|_1 + \| a_1-\bar{a}_1\|_1 \|b_3\|_1 + \|a_2-\bar{a}_2\|_1 \|b_1\|_1 + \|a_1-\bar{a}_1\|_1 \|b_2\|_1
)
\\
& \le \|A\|_{B(\mathcal{X})} \| \mathcal{K} \|_{B(\mathcal{X})} \|a-\bar{a}\|_{\mathcal{X}} \|b\|_{\mathcal{X}}\\
& = \frac{1}{2|\lambda|} \|A\|_{B(\mathcal{X})} \|a-\bar{a}\|_{\mathcal{X}} \|b\|_{\mathcal{X}}.
\end{align}$

Computing $Z_2$ boils down to compute (or bound) $\frac{1}{2|\lambda|} \|A\|_{B(\mathcal{X})}$.
"""

# ╔═╡ 6688fe5d-0019-485d-9059-bb5bd2be7d6e
Z₂ = (2/abs(tλ)) * opnormA

# ╔═╡ 5d0b6b50-0170-4f8b-ad5a-0dda212a00e1
md"""
## Finishing the CAP
"""

# ╔═╡ a9f597fe-4058-47a8-8a39-8d6b36594a96
md"""
We now evaluate the radii polynomial to finish the Newton-Kantorovich proof.
"""

# ╔═╡ 1082f423-448c-4d74-8baa-76b2e6508084
begin
r_star = Inf
r = interval_of_existence(Y, Z₁, Z₂, r_star)
end

# ╔═╡ f82c6a7a-fffb-4a7b-97df-bfafd2072f5c
md"""
## Your main task now: the following Exercise.
"""

# ╔═╡ 9a0d1fd7-e8b6-4737-b764-95f917442335
Markdown.MD(Markdown.Admonition("note", "Exercise",[md"""Using Taylor series, compute rigorously (that is with the Newton-Kantorovich approach) a parameterization of order $N=20$ of the stable and unstable manifold of the steady state $(0,0)$ for the equation

$\begin{align}
\dot x = \begin{pmatrix} \dot x_1  \\ \dot x_2 \end{pmatrix} =
\begin{pmatrix}
x_1+6 x_2 + x_1 x_2
\\
4 x_1 + 3 x_2 - x_1^2
\end{pmatrix}.
\end{align}$
"""]))

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
GLMakie = "e9467ef8-e4e7-5192-8a1a-b1aee30e663a"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
OrdinaryDiffEq = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
RadiiPolynomial = "f2081a94-c849-46b6-8dc9-07bb90ed72a9"
ToeplitzMatrices = "c751599d-da0a-543b-9d20-d0a503d91d24"

[compat]
GLMakie = "~0.10.9"
OrdinaryDiffEq = "~6.87.0"
Plots = "~1.40.5"
PlutoTeachingTools = "~0.2.15"
PlutoUI = "~0.7.59"
RadiiPolynomial = "~0.8.12"
ToeplitzMatrices = "~0.8.4"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.5"
manifest_format = "2.0"
project_hash = "ce363e2746a58f11b3b70d064fa12cdf20c8a8d7"

[[deps.ADTypes]]
git-tree-sha1 = "99a6f5d0ce1c7c6afdb759daa30226f71c54f6b0"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.7.1"
weakdeps = ["ChainRulesCore", "EnzymeCore"]

    [deps.ADTypes.extensions]
    ADTypesChainRulesCoreExt = "ChainRulesCore"
    ADTypesEnzymeCoreExt = "EnzymeCore"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "LinearAlgebra", "MacroTools", "Markdown", "Test"]
git-tree-sha1 = "f61b15be1d76846c0ce31d3fcfac5380ae53db6a"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.37"

    [deps.Accessors.extensions]
    AccessorsAxisKeysExt = "AxisKeys"
    AccessorsIntervalSetsExt = "IntervalSets"
    AccessorsStaticArraysExt = "StaticArrays"
    AccessorsStructArraysExt = "StructArrays"
    AccessorsUnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    Requires = "ae029012-a4dd-5104-9daa-d747884805df"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "6a55b747d1812e699320963ffde36f1ebdda4099"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.0.4"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AdaptivePredicates]]
git-tree-sha1 = "7d5da5dd472490d048b081ca1bda4a7821b06456"
uuid = "35492f91-a3bd-45ad-95db-fcad7dcfedb7"
version = "1.1.1"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e81c509d2c8e49592413bfb0bb3b08150056c79d"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.1"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "d57bd3762d308bded22c3b82d033bff85f6195c6"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.4.0"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "3640d077b6dafd64ceb8fd5c1ec76f7ca53bcf76"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.16.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = "CUDSS"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceSparseArraysExt = "SparseArrays"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra"]
git-tree-sha1 = "0dd7edaff278e346eb0ca07a7e75c9438408a3ce"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "1.10.3"
weakdeps = ["SparseArrays"]

    [deps.ArrayLayouts.extensions]
    ArrayLayoutsSparseArraysExt = "SparseArrays"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Automa]]
deps = ["PrecompileTools", "TranscodingStreams"]
git-tree-sha1 = "014bc22d6c400a7703c0f5dc1fdc302440cf88be"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "1.0.4"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "16351be62963a67ac4083f748fdb3cca58bfd52f"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.7"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "f21cfd4950cb9f0587d5067e69405ad2acd27b87"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.6"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9e2a6b69137e6969bab0152632dcb3bc108c8bdd"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+1"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "PrecompileTools", "Static"]
git-tree-sha1 = "5a97e67919535d6841172016c9530fd69494e5ec"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.2.6"

[[deps.CRC32c]]
uuid = "8bf52ea8-c179-5cab-976a-9e18b702a9bc"

[[deps.CRlibm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e329286945d0cfc04456972ea732551869af1cfc"
uuid = "4e9b3aee-d8a1-5a3d-ad8b-7d824db253f0"
version = "1.0.1+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a2f1c8c668c8e3cb4cca4e57a8efdb09067bb3fd"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.0+2"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "71acdbf594aab5bbb2cec89b208c41b4c411e49f"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.24.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CloseOpenIntervals]]
deps = ["Static", "StaticArrayInterface"]
git-tree-sha1 = "05ba0d07cd4fd8b7a39541e31a7b0254704ea581"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.13"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "c0216e792f518b39b22212127d4a84dc31e4e386"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.3.5"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "b8fe8546d52ca154ac556809e10c75e6e7430ac8"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.5"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON", "Test"]
git-tree-sha1 = "61c5334f33d91e570e1d0c3eb5465835242582c4"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "4b270d6465eb21ae89b732182c20dc165f8bf9f2"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.25.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "362a287c3aa50601b0bc359053d5c2468f0e7ce0"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.11"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "b1c55339b7c6c350ee89f2c1604299660525b248"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.15.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConcreteStructs]]
git-tree-sha1 = "f749037478283d372048690eb3b5f92a79432b34"
uuid = "2569d6c7-a4a2-43d3-a901-331e8e4be471"
version = "0.2.3"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "ea32b83ca4fefa1768dc84e504cc0a94fb1ab8d1"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.4.2"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d8a9c0b6ac2d9081bf76324b39c78ca3ce4f0c98"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.6"
weakdeps = ["IntervalSets", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "fcbb72b032692610bfbdb15018ac16a36cf2e406"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.1"

[[deps.DSP]]
deps = ["Compat", "FFTW", "IterTools", "LinearAlgebra", "Polynomials", "Random", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "f7f4319567fe769debfcf7f8c03d8da1dd4e2fb0"
uuid = "717857b8-e6f2-59f4-9121-6e50c889abd2"
version = "0.7.9"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fc173b380865f70627d7dd1190dc2fce6cc105af"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.14.10+0"

[[deps.DelaunayTriangulation]]
deps = ["AdaptivePredicates", "EnumX", "ExactPredicates", "Random"]
git-tree-sha1 = "46f12daa85e5acc0ea5d5f9f8c3f1fc679e0f7e5"
uuid = "927a84f5-c5f4-47a5-9785-b46e178433df"
version = "1.2.0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DiffEqBase]]
deps = ["ArrayInterface", "ConcreteStructs", "DataStructures", "DocStringExtensions", "EnumX", "EnzymeCore", "FastBroadcast", "FastClosures", "ForwardDiff", "FunctionWrappers", "FunctionWrappersWrappers", "LinearAlgebra", "Logging", "Markdown", "MuladdMacro", "Parameters", "PreallocationTools", "PrecompileTools", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "Setfield", "SparseArrays", "Static", "StaticArraysCore", "Statistics", "Tricks", "TruncatedStacktraces"]
git-tree-sha1 = "d1e8a4642e28b0945bde6e2e1ac569b9e0abd728"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.151.5"

    [deps.DiffEqBase.extensions]
    DiffEqBaseCUDAExt = "CUDA"
    DiffEqBaseChainRulesCoreExt = "ChainRulesCore"
    DiffEqBaseDistributionsExt = "Distributions"
    DiffEqBaseEnzymeExt = ["ChainRulesCore", "Enzyme"]
    DiffEqBaseGeneralizedGeneratedExt = "GeneralizedGenerated"
    DiffEqBaseMPIExt = "MPI"
    DiffEqBaseMeasurementsExt = "Measurements"
    DiffEqBaseMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    DiffEqBaseReverseDiffExt = "ReverseDiff"
    DiffEqBaseTrackerExt = "Tracker"
    DiffEqBaseUnitfulExt = "Unitful"

    [deps.DiffEqBase.weakdeps]
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    GeneralizedGenerated = "6b9d7cbe-bcb9-11e9-073f-15a7a543e2eb"
    MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.DifferentiationInterface]]
deps = ["ADTypes", "Compat", "DocStringExtensions", "FillArrays", "LinearAlgebra", "PackageExtensionCompat", "SparseArrays", "SparseMatrixColorings"]
git-tree-sha1 = "9b23f9a816790b8ab9914c3c86321a546e92cbe7"
uuid = "a0c0ee7d-e4b9-4e03-894e-1c5f64a51d63"
version = "0.5.17"

    [deps.DifferentiationInterface.extensions]
    DifferentiationInterfaceChainRulesCoreExt = "ChainRulesCore"
    DifferentiationInterfaceDiffractorExt = "Diffractor"
    DifferentiationInterfaceEnzymeExt = "Enzyme"
    DifferentiationInterfaceFastDifferentiationExt = "FastDifferentiation"
    DifferentiationInterfaceFiniteDiffExt = "FiniteDiff"
    DifferentiationInterfaceFiniteDifferencesExt = "FiniteDifferences"
    DifferentiationInterfaceForwardDiffExt = "ForwardDiff"
    DifferentiationInterfacePolyesterForwardDiffExt = "PolyesterForwardDiff"
    DifferentiationInterfaceReverseDiffExt = "ReverseDiff"
    DifferentiationInterfaceSymbolicsExt = "Symbolics"
    DifferentiationInterfaceTapirExt = "Tapir"
    DifferentiationInterfaceTrackerExt = "Tracker"
    DifferentiationInterfaceZygoteExt = ["Zygote", "ForwardDiff"]

    [deps.DifferentiationInterface.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Diffractor = "9f5e2b26-1114-432f-b630-d3fe2085c51c"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    FastDifferentiation = "eb9bf01b-bf85-4b60-bf87-ee5de06c00be"
    FiniteDiff = "6a86dc24-6348-571c-b903-95158fe2bd41"
    FiniteDifferences = "26cc04aa-876d-5657-8c51-4c34ba976000"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    PolyesterForwardDiff = "98d1487c-24ca-40b6-b7ab-df2af84e126b"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
    Tapir = "07d77754-e150-4737-8c94-cd238a1fb45b"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "e6c693a0e4394f8fda0e51a5bdf5aef26f8235e9"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.111"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[deps.EnzymeCore]]
git-tree-sha1 = "8f205a601760f4798a10f138c3940f0451d95188"
uuid = "f151be2c-9106-41f4-ab19-57ee4f262869"
version = "0.7.8"
weakdeps = ["Adapt"]

    [deps.EnzymeCore.extensions]
    AdaptExt = "Adapt"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8e9441ee83492030ace98f9789a654a6d0b1f643"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+0"

[[deps.ExactPredicates]]
deps = ["IntervalArithmetic", "Random", "StaticArrays"]
git-tree-sha1 = "b3f2ff58735b5f024c392fde763f29b057e4b025"
uuid = "429591f6-91af-11e9-00e2-59fbe8cec110"
version = "2.2.8"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "dcb08a0d93ec0b1cdc4af184b26b591e9695423a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.10"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c6317308b9dc757616f0b5cb379db10494443a7"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.2+0"

[[deps.ExponentialUtilities]]
deps = ["Adapt", "ArrayInterface", "GPUArraysCore", "GenericSchur", "LinearAlgebra", "PrecompileTools", "Printf", "SparseArrays", "libblastrampoline_jll"]
git-tree-sha1 = "8e18940a5ba7f4ddb41fe2b79b6acaac50880a86"
uuid = "d4d017d3-3776-5f7e-afef-a10c40355c18"
version = "1.26.1"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.Expronicon]]
deps = ["MLStyle", "Pkg", "TOML"]
git-tree-sha1 = "fc3951d4d398b5515f91d7fe5d45fc31dccb3c9b"
uuid = "6b7a57c9-7cc1-4fdf-b7f5-e857abae3636"
version = "0.8.5"

[[deps.Extents]]
git-tree-sha1 = "81023caa0021a41712685887db1fc03db26f41f5"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.4"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "4820348781ae578893311153d69049a93d05f39d"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.8.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FastBroadcast]]
deps = ["ArrayInterface", "LinearAlgebra", "Polyester", "Static", "StaticArrayInterface", "StrideArraysCore"]
git-tree-sha1 = "ab1b34570bcdf272899062e1a56285a53ecaae08"
uuid = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
version = "0.3.5"

[[deps.FastClosures]]
git-tree-sha1 = "acebe244d53ee1b461970f8910c235b259e772ef"
uuid = "9aa1b823-49e4-5ca5-8b0f-3971ec8bab6a"
version = "0.3.2"

[[deps.FastLapackInterface]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "cbf5edddb61a43669710cbc2241bc08b36d9e660"
uuid = "29a986be-02c6-4525-aec4-84b980013641"
version = "2.0.4"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "82d8afa92ecf4b52d78d869f038ebfb881267322"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.3"

[[deps.FilePaths]]
deps = ["FilePathsBase", "MacroTools", "Reexport", "Requires"]
git-tree-sha1 = "919d9412dbf53a2e6fe74af62a73ceed0bce0629"
uuid = "8fc22ac5-c921-52a6-82fd-178b2807b824"
version = "0.8.3"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates"]
git-tree-sha1 = "7878ff7172a8e6beedd1dea14bd27c3c6340d361"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.22"
weakdeps = ["Mmap", "Test"]

    [deps.FilePathsBase.extensions]
    FilePathsBaseMmapExt = "Mmap"
    FilePathsBaseTestExt = "Test"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "0653c0a2396a6da5bc4766c43041ef5fd3efbe57"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.11.0"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Setfield", "SparseArrays"]
git-tree-sha1 = "f9219347ebf700e77ca1d48ef84e4a82a6701882"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.24.0"

    [deps.FiniteDiff.extensions]
    FiniteDiffBandedMatricesExt = "BandedMatrices"
    FiniteDiffBlockBandedMatricesExt = "BlockBandedMatrices"
    FiniteDiffStaticArraysExt = "StaticArrays"

    [deps.FiniteDiff.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "db16beca600632c95fc8aca29890d83788dd8b23"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.96+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "cf0fe81336da9fb90944683b8c41984b08793dad"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.36"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "907369da0f8e80728ab49c1c7e09327bf0d6d999"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.1.1"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "5c1d8ae0efc6c2e7b1fc502cbe25def8f661b7bc"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.2+0"

[[deps.FreeTypeAbstraction]]
deps = ["ColorVectorSpace", "Colors", "FreeType", "GeometryBasics"]
git-tree-sha1 = "2493cdfd0740015955a8e46de4ef28f49460d8bc"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.10.3"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1ed150b39aebcc805c26b93a8d0122c940f64ce2"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.14+0"

[[deps.FunctionWrappers]]
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[deps.FunctionWrappersWrappers]]
deps = ["FunctionWrappers"]
git-tree-sha1 = "b104d487b34566608f8b4e1c39fb0b10aa279ff8"
uuid = "77dc65aa-8811-40c2-897b-53d922fa7daf"
version = "0.1.3"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW]]
deps = ["GLFW_jll"]
git-tree-sha1 = "7ed24cfc4cb29fb10c0e8cca871ddff54c32a4c3"
uuid = "f7f18e0c-5ee9-5ccd-a5bf-e8befd85ed98"
version = "3.4.3"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "3f74912a156096bd8fdbef211eff66ab446e7297"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.0+0"

[[deps.GLMakie]]
deps = ["ColorTypes", "Colors", "FileIO", "FixedPointNumbers", "FreeTypeAbstraction", "GLFW", "GeometryBasics", "LinearAlgebra", "Makie", "Markdown", "MeshIO", "ModernGL", "Observables", "PrecompileTools", "Printf", "ShaderAbstractions", "StaticArrays"]
git-tree-sha1 = "f644c26d8b4d0d9d6ac6ff1f7fefe7b8f70c3e92"
uuid = "e9467ef8-e4e7-5192-8a1a-b1aee30e663a"
version = "0.10.9"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "ec632f177c0d990e64d955ccc1b8c04c485a0950"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.6"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "3e527447a45901ea392fe12120783ad6ec222803"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.6"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "182c478a179b267dd7a741b6f8f4c3e0803795d6"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.6+0"

[[deps.GenericSchur]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "af49a0851f8113fcfae2ef5027c6d49d0acec39b"
uuid = "c145ed77-6b09-5dd9-b285-bf645a82121e"
version = "0.5.4"

[[deps.GeoFormatTypes]]
git-tree-sha1 = "59107c179a586f0fe667024c5eb7033e81333271"
uuid = "68eda718-8dee-11e9-39e7-89f7f65f511f"
version = "0.4.2"

[[deps.GeoInterface]]
deps = ["Extents", "GeoFormatTypes"]
git-tree-sha1 = "5921fc0704e40c024571eca551800c699f86ceb4"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.3.6"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "Extents", "GeoInterface", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "b62f2b2d76cee0d61a2ef2b3118cd2a3215d3134"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.11"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "7c82e6a6cd34e9d935e9aa4051b66c6ff3af59ba"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.80.2+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "ebd18c326fa6cee1efb7da9a3b45cf69da2ed4d9"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.11.2"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "fc713f007cff99ff9e50accba6373624ddd33588"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.11.0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "d1d712be3164d61d1fb98e7ce9bcbc6cc06b45ed"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.8"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "8e070b599339d622e9a081d17230d74a5c473293"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.17"

[[deps.HypergeometricFunctions]]
deps = ["LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "7c4195be1649ae622304031ed46a2f4df989f1eb"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.24"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "2e4520d67b0cef90865b3ef727594d2a58e0e1f8"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.11"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "eb49b82c172811fd2c86759fa0553a2221feb909"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.7"

[[deps.ImageCore]]
deps = ["ColorVectorSpace", "Colors", "FixedPointNumbers", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "PrecompileTools", "Reexport"]
git-tree-sha1 = "b2a7eaa169c13f5bcae8131a83bc30eff8f71be0"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.10.2"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs"]
git-tree-sha1 = "437abb322a41d527c197fa800455f79d414f0a3c"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.8"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "355e2b974f2e3212a75dfb60519de21361ad3cb7"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.9"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0936ba688c6d201805a83da835b55c61a180db52"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.11+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "14eb2b542e748570b56446f4c50fbfb2306ebc45"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2024.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "88a101217d7cb38a7b481ccd50d21876e1d1b0e0"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.15.1"
weakdeps = ["Unitful"]

    [deps.Interpolations.extensions]
    InterpolationsUnitfulExt = "Unitful"

[[deps.IntervalArithmetic]]
deps = ["CRlibm_jll", "MacroTools", "RoundingEmulator"]
git-tree-sha1 = "433b0bb201cd76cb087b017e49244f10394ebe9c"
uuid = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
version = "0.22.14"
weakdeps = ["DiffRules", "ForwardDiff", "RecipesBase"]

    [deps.IntervalArithmetic.extensions]
    IntervalArithmeticDiffRulesExt = "DiffRules"
    IntervalArithmeticForwardDiffExt = "ForwardDiff"
    IntervalArithmeticRecipesBaseExt = "RecipesBase"

[[deps.IntervalSets]]
git-tree-sha1 = "dba9ddf07f77f60450fe5d2e2beb9854d9a49bd0"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.10"
weakdeps = ["Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.InverseFunctions]]
git-tree-sha1 = "2787db24f4e03daf859c6509ff87764e4182f7d1"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.16"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.Isoband]]
deps = ["isoband_jll"]
git-tree-sha1 = "f9b6d97355599074dc867318950adaa6f9946137"
uuid = "f1662d9f-8043-43de-a69a-05efc1cc6ff4"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "a53ebe394b71470c7f97c2e7e170d51df21b17af"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.7"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "fa6d0bcff8583bac20f1ffa708c3913ca605c611"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.5"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c84a835e1a09b289ffcd2271bf2a337bbdda6637"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.3+0"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "a6adc2dcfe4187c40dc7c2c9d2128e326360e90a"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.32"

[[deps.KLU]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse_jll"]
git-tree-sha1 = "07649c499349dad9f08dde4243a4c597064663e9"
uuid = "ef3ab10e-7fda-4108-b977-705223b18434"
version = "0.6.0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "7d703202e65efa1369de1279c162b915e245eed1"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.9"

[[deps.Krylov]]
deps = ["LinearAlgebra", "Printf", "SparseArrays"]
git-tree-sha1 = "267dad6b4b7b5d529c76d40ff48d33f7e94cb834"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.9.6"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "170b660facf5df5de098d866564877e119141cbd"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.2+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d986ce2d884d49126836ea94ed5bfb0f12679713"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "70c5da094887fd2cae843b8db33920bac4b6f07d"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.2+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "5b0d630f3020b82c0775a51d05895852f8506f50"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.4"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "a9eaadb366f5493a5654e843864c13d8b107548c"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.17"

[[deps.LazyArrays]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "MacroTools", "SparseArrays"]
git-tree-sha1 = "507b423197fdd9e77b74aa2532c0a05eb7eb4004"
uuid = "5078a376-72f3-5289-bfd5-ec5146d43c02"
version = "2.2.0"

    [deps.LazyArrays.extensions]
    LazyArraysBandedMatricesExt = "BandedMatrices"
    LazyArraysBlockArraysExt = "BlockArrays"
    LazyArraysBlockBandedMatricesExt = "BlockBandedMatrices"
    LazyArraysStaticArraysExt = "StaticArrays"

    [deps.LazyArrays.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockArrays = "8e7c35d0-a365-5155-bbbb-fb81a777f24e"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll"]
git-tree-sha1 = "9fd170c4bbfd8b935fdc5f8b7aa33532c991a673"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.11+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fbb1f2bef882392312feb1ede3615ddc1e9b99ed"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.49.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0c4f9c4f1a50d8f35048fa0532dabbadf702f81e"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.40.1+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "2da088d113af58221c52828a80378e16be7d037a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.5.1+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "5ee6203157c120d79034c748a2acba45b82b8807"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.40.1+0"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "e4c3be53733db1051cc15ecf573b1042b3a712a1"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.3.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LinearSolve]]
deps = ["ArrayInterface", "ChainRulesCore", "ConcreteStructs", "DocStringExtensions", "EnumX", "FastLapackInterface", "GPUArraysCore", "InteractiveUtils", "KLU", "Krylov", "LazyArrays", "Libdl", "LinearAlgebra", "MKL_jll", "Markdown", "PrecompileTools", "Preferences", "RecursiveFactorization", "Reexport", "SciMLBase", "SciMLOperators", "Setfield", "SparseArrays", "Sparspak", "StaticArraysCore", "UnPack"]
git-tree-sha1 = "6c5e4555ac2bc449a28604e184f759d18fc08420"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "2.34.0"

    [deps.LinearSolve.extensions]
    LinearSolveBandedMatricesExt = "BandedMatrices"
    LinearSolveBlockDiagonalsExt = "BlockDiagonals"
    LinearSolveCUDAExt = "CUDA"
    LinearSolveCUDSSExt = "CUDSS"
    LinearSolveEnzymeExt = ["Enzyme", "EnzymeCore"]
    LinearSolveFastAlmostBandedMatricesExt = ["FastAlmostBandedMatrices"]
    LinearSolveHYPREExt = "HYPRE"
    LinearSolveIterativeSolversExt = "IterativeSolvers"
    LinearSolveKernelAbstractionsExt = "KernelAbstractions"
    LinearSolveKrylovKitExt = "KrylovKit"
    LinearSolveMetalExt = "Metal"
    LinearSolvePardisoExt = "Pardiso"
    LinearSolveRecursiveArrayToolsExt = "RecursiveArrayTools"

    [deps.LinearSolve.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockDiagonals = "0a1fb500-61f7-11e9-3c65-f5ef3456f9f0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    FastAlmostBandedMatrices = "9d29842c-ecb8-4973-b1e9-a27b1157504e"
    HYPRE = "b5ffcf37-a2bd-41ab-a3da-4bd9bc8ad771"
    IterativeSolvers = "42fd0dbc-a981-5370-80f2-aaf504508153"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    KrylovKit = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    Pardiso = "46dd5b70-b6fb-5a00-ae2d-e8fea33afaf2"
    RecursiveArrayTools = "731186ca-8d62-57ce-b412-fbd966d074cd"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "a2d09619db4e765091ee5c6ffe8872849de0feea"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.28"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "c1dd6d7978c12545b4179fb6153b9250c96b0075"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.3"

[[deps.LoopVectorization]]
deps = ["ArrayInterface", "CPUSummary", "CloseOpenIntervals", "DocStringExtensions", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "PrecompileTools", "SIMDTypes", "SLEEFPirates", "Static", "StaticArrayInterface", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "8084c25a250e00ae427a379a5b607e7aed96a2dd"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.171"
weakdeps = ["ChainRulesCore", "ForwardDiff", "SpecialFunctions"]

    [deps.LoopVectorization.extensions]
    ForwardDiffExt = ["ChainRulesCore", "ForwardDiff"]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "eeaedcf337f33c039f9f3a209a8db992deefd7e9"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "2.4.8"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "f046ccd0c6db2832a9f639e2c669c6fe867e5f4f"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2024.2.0+0"

[[deps.MLStyle]]
git-tree-sha1 = "bc38dff0548128765760c79eb7388a4b37fae2c8"
uuid = "d8e11817-5142-5d16-987a-aa16d5891078"
version = "0.4.17"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.Makie]]
deps = ["Animations", "Base64", "CRC32c", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "Contour", "Dates", "DelaunayTriangulation", "Distributions", "DocStringExtensions", "Downloads", "FFMPEG_jll", "FileIO", "FilePaths", "FixedPointNumbers", "Format", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageIO", "InteractiveUtils", "IntervalSets", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MacroTools", "MakieCore", "Markdown", "MathTeXEngine", "Observables", "OffsetArrays", "Packing", "PlotUtils", "PolygonOps", "PrecompileTools", "Printf", "REPL", "Random", "RelocatableFolders", "Scratch", "ShaderAbstractions", "Showoff", "SignedDistanceFields", "SparseArrays", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "TriplotBase", "UnicodeFun", "Unitful"]
git-tree-sha1 = "204f06860af9008fa08b3a4842f48116e1209a2c"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.21.9"

[[deps.MakieCore]]
deps = ["ColorTypes", "GeometryBasics", "IntervalSets", "Observables"]
git-tree-sha1 = "b0e2e3473af351011e598f9219afb521121edd2b"
uuid = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
version = "0.8.6"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.MappedArrays]]
git-tree-sha1 = "2dab0221fe2b0f2cb6754eaa743cc266339f527e"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.2"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "UnicodeFun"]
git-tree-sha1 = "e1641f32ae592e415e3dbae7f4a188b5316d4b62"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.6.1"

[[deps.MaybeInplace]]
deps = ["ArrayInterface", "LinearAlgebra", "MacroTools"]
git-tree-sha1 = "54e2fdc38130c05b42be423e90da3bade29b74bd"
uuid = "bb5d69b7-63fc-4a16-80bd-7e42200c7bdb"
version = "0.1.4"
weakdeps = ["SparseArrays"]

    [deps.MaybeInplace.extensions]
    MaybeInplaceSparseArraysExt = "SparseArrays"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.MeshIO]]
deps = ["ColorTypes", "FileIO", "GeometryBasics", "Printf"]
git-tree-sha1 = "dc182956229ff16d5a4d90a562035e633bd2561d"
uuid = "7269a6da-0436-5bbc-96c2-40638cbb6118"
version = "0.4.12"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.ModernGL]]
deps = ["Libdl"]
git-tree-sha1 = "b76ea40b5c0f45790ae09492712dd326208c28b2"
uuid = "66fc600b-dfda-50eb-8b99-91cfa97b1301"
version = "1.1.7"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.MuladdMacro]]
git-tree-sha1 = "cac9cc5499c25554cba55cd3c30543cff5ca4fab"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.4"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "a0b464d183da839699f4c79e7606d9d186ec172c"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.3"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.NonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "ConcreteStructs", "DiffEqBase", "FastBroadcast", "FastClosures", "FiniteDiff", "ForwardDiff", "LazyArrays", "LineSearches", "LinearAlgebra", "LinearSolve", "MaybeInplace", "PrecompileTools", "Preferences", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SimpleNonlinearSolve", "SparseArrays", "SparseDiffTools", "StaticArraysCore", "SymbolicIndexingInterface", "TimerOutputs"]
git-tree-sha1 = "bcd8812e751326ff1d4b2dd50764b93df51f143b"
uuid = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
version = "3.14.0"

    [deps.NonlinearSolve.extensions]
    NonlinearSolveBandedMatricesExt = "BandedMatrices"
    NonlinearSolveFastLevenbergMarquardtExt = "FastLevenbergMarquardt"
    NonlinearSolveFixedPointAccelerationExt = "FixedPointAcceleration"
    NonlinearSolveLeastSquaresOptimExt = "LeastSquaresOptim"
    NonlinearSolveMINPACKExt = "MINPACK"
    NonlinearSolveNLSolversExt = "NLSolvers"
    NonlinearSolveNLsolveExt = "NLsolve"
    NonlinearSolveSIAMFANLEquationsExt = "SIAMFANLEquations"
    NonlinearSolveSpeedMappingExt = "SpeedMapping"
    NonlinearSolveSymbolicsExt = "Symbolics"
    NonlinearSolveZygoteExt = "Zygote"

    [deps.NonlinearSolve.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    FastLevenbergMarquardt = "7a0df574-e128-4d35-8cbd-3d84502bf7ce"
    FixedPointAcceleration = "817d07cb-a79a-5c30-9a31-890123675176"
    LeastSquaresOptim = "0fc2ff8b-aaa3-5acd-a817-1944a5e08891"
    MINPACK = "4854310b-de5a-5eb6-a2a5-c1dee2bd17f9"
    NLSolvers = "337daf1e-9722-11e9-073e-8b9effe078ba"
    NLsolve = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
    SIAMFANLEquations = "084e46ad-d928-497d-ad5e-07fa361a48c4"
    SpeedMapping = "f1835b91-879b-4a3f-a438-e4baacf14412"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
git-tree-sha1 = "1a27764e945a152f7ca7efa04de513d473e9542e"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.14.1"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "327f53360fdb54df7ecd01e96ef1983536d1e633"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.2"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "8292dd5c8a38257111ada2174000a33745b06d4e"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.2.4+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "38cb508d080d21dc1128f7fb04f20387ed4c0af4"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.3"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a028ee3cb5641cccc4c24e90c36b0a4f7707bdf5"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.14+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.OrdinaryDiffEq]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "EnumX", "ExponentialUtilities", "FastBroadcast", "FastClosures", "FillArrays", "FiniteDiff", "ForwardDiff", "FunctionWrappersWrappers", "IfElse", "InteractiveUtils", "LineSearches", "LinearAlgebra", "LinearSolve", "Logging", "MacroTools", "MuladdMacro", "NonlinearSolve", "Polyester", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SciMLStructures", "SimpleNonlinearSolve", "SimpleUnPack", "SparseArrays", "SparseDiffTools", "Static", "StaticArrayInterface", "StaticArrays", "TruncatedStacktraces"]
git-tree-sha1 = "a8b2d333cd90562b58b977b4033739360b37fb1f"
uuid = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
version = "6.87.0"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "949347156c25054de2db3b166c52ac4728cbad65"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.31"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "67186a2bc9a90f9f85ff3cc8277868961fb57cbd"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.3"

[[deps.PackageExtensionCompat]]
git-tree-sha1 = "fb28e33b8a95c4cee25ce296c817d89cc2e53518"
uuid = "65ce6f38-6b18-4e1d-a461-8949797d7930"
version = "1.0.2"
weakdeps = ["Requires", "TOML"]

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "ec3edfe723df33528e085e632414499f26650501"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.5.0"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "cb5a2ab6763464ae0f19c86c56c63d4a2b0f5bda"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.52.2+0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "35621f10a7531bc8fa58f74610b1bfb70a3cfc6b"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.43.4+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "6e55c6841ce3411ccb3457ee52fc48cb698d6fb0"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.2.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "7b1a9df27f072ac4c9c7cbe5efb198489258d1f5"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.1"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "082f0c4b70c202c37784ce4bfbc33c9f437685bf"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.5"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoHooks]]
deps = ["InteractiveUtils", "Markdown", "UUIDs"]
git-tree-sha1 = "072cdf20c9b0507fdd977d7d246d90030609674b"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0774"
version = "0.0.5"

[[deps.PlutoLinks]]
deps = ["FileWatching", "InteractiveUtils", "Markdown", "PlutoHooks", "Revise", "UUIDs"]
git-tree-sha1 = "8f5fa7056e6dcfb23ac5211de38e6c03f6367794"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0420"
version = "0.1.6"

[[deps.PlutoTeachingTools]]
deps = ["Downloads", "HypertextLiteral", "LaTeXStrings", "Latexify", "Markdown", "PlutoLinks", "PlutoUI", "Random"]
git-tree-sha1 = "5d9ab1a4faf25a62bb9d07ef0003396ac258ef1c"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.2.15"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "ab55ee1510ad2af0ff674dbcced5e94921f867a9"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.59"

[[deps.Polyester]]
deps = ["ArrayInterface", "BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "ManualMemory", "PolyesterWeave", "Static", "StaticArrayInterface", "StrideArraysCore", "ThreadingUtilities"]
git-tree-sha1 = "6d38fea02d983051776a856b7df75b30cf9a3c1f"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.7.16"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "645bed98cd47f72f67316fd42fc47dee771aefcd"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.2.2"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "RecipesBase", "Requires", "Setfield", "SparseArrays"]
git-tree-sha1 = "1a9cfb2dc2c2f1bd63f1906d72af39a79b49b736"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "4.0.11"

    [deps.Polynomials.extensions]
    PolynomialsChainRulesCoreExt = "ChainRulesCore"
    PolynomialsFFTWExt = "FFTW"
    PolynomialsMakieCoreExt = "MakieCore"
    PolynomialsMutableArithmeticsExt = "MutableArithmetics"

    [deps.Polynomials.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
    MakieCore = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
    MutableArithmetics = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "ForwardDiff"]
git-tree-sha1 = "6c62ce45f268f3f958821a1e5192cf91c75ae89c"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.24"

    [deps.PreallocationTools.extensions]
    PreallocationToolsReverseDiffExt = "ReverseDiff"

    [deps.PreallocationTools.weakdeps]
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "8f6bc219586aef8baf0ff9a5fe16ee9c70cb65e4"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.10.2"

[[deps.PtrArrays]]
git-tree-sha1 = "77a42d78b6a92df47ab37e177b2deac405e1c88f"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.2.1"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "18e8f4d1426e965c7b532ddd260599e1510d26ce"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.0"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "492601870742dcd38f233b23c3ec629628c1d724"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.7.1+1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "1d587203cf851a51bf1ea31ad7ff89eff8d625ea"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.0"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.RadiiPolynomial]]
deps = ["IntervalArithmetic", "LinearAlgebra", "Printf", "Reexport", "SparseArrays"]
git-tree-sha1 = "89f57ab86310e5ca7009cb236441505ba6b3242a"
uuid = "f2081a94-c849-46b6-8dc9-07bb90ed72a9"
version = "0.8.12"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "b034171b93aebc81b3e1890a036d13a9c4a9e3e0"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "3.27.0"

    [deps.RecursiveArrayTools.extensions]
    RecursiveArrayToolsFastBroadcastExt = "FastBroadcast"
    RecursiveArrayToolsForwardDiffExt = "ForwardDiff"
    RecursiveArrayToolsMeasurementsExt = "Measurements"
    RecursiveArrayToolsMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    RecursiveArrayToolsReverseDiffExt = ["ReverseDiff", "Zygote"]
    RecursiveArrayToolsSparseArraysExt = ["SparseArrays"]
    RecursiveArrayToolsTrackerExt = "Tracker"
    RecursiveArrayToolsZygoteExt = "Zygote"

    [deps.RecursiveArrayTools.weakdeps]
    FastBroadcast = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.RecursiveFactorization]]
deps = ["LinearAlgebra", "LoopVectorization", "Polyester", "PrecompileTools", "StrideArraysCore", "TriangularSolve"]
git-tree-sha1 = "6db1a75507051bc18bfa131fbc7c3f169cc4b2f6"
uuid = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
version = "0.2.23"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Revise]]
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "Pkg", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "85ddd93ea15dcd8493400600e09104a9e94bb18d"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.5.15"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e60724fd3beea548353984dc61c943ecddb0e29a"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.3+0"

[[deps.RoundingEmulator]]
git-tree-sha1 = "40b9edad2e5287e05bd413a38f61a8ff55b9557b"
uuid = "5eaf0fd0-dfba-4ccb-bf02-d820a40db705"
version = "0.2.1"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "04c968137612c4a5629fa531334bb81ad5680f00"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.13"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "2803cab51702db743f3fda07dd1745aadfbf43bd"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.5.0"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "456f610ca2fbd1c14f5fcf31c6bfadc55e7d66e0"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.43"

[[deps.SciMLBase]]
deps = ["ADTypes", "Accessors", "ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "Expronicon", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "PrecompileTools", "Preferences", "Printf", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "SciMLStructures", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "8001043f80051c86f264fd6e936d97e6b9eff401"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "2.52.0"

    [deps.SciMLBase.extensions]
    SciMLBaseChainRulesCoreExt = "ChainRulesCore"
    SciMLBaseMakieExt = "Makie"
    SciMLBasePartialFunctionsExt = "PartialFunctions"
    SciMLBasePyCallExt = "PyCall"
    SciMLBasePythonCallExt = "PythonCall"
    SciMLBaseRCallExt = "RCall"
    SciMLBaseZygoteExt = "Zygote"

    [deps.SciMLBase.weakdeps]
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    PartialFunctions = "570af359-4316-4cb7-8c74-252c00c2016b"
    PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
    PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
    RCall = "6f49c342-dc21-5d91-9882-a32aef131414"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SciMLOperators]]
deps = ["Accessors", "ArrayInterface", "DocStringExtensions", "LinearAlgebra", "MacroTools"]
git-tree-sha1 = "e39c5f217f9aca640c8e27ab21acf557a3967db5"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "0.3.10"
weakdeps = ["SparseArrays", "StaticArraysCore"]

    [deps.SciMLOperators.extensions]
    SciMLOperatorsSparseArraysExt = "SparseArrays"
    SciMLOperatorsStaticArraysCoreExt = "StaticArraysCore"

[[deps.SciMLStructures]]
deps = ["ArrayInterface"]
git-tree-sha1 = "25514a6f200219cd1073e4ff23a6324e4a7efe64"
uuid = "53ae85a6-f571-4167-b2af-e1d143709226"
version = "1.5.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.ShaderAbstractions]]
deps = ["ColorTypes", "FixedPointNumbers", "GeometryBasics", "LinearAlgebra", "Observables", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "79123bc60c5507f035e6d1d9e563bb2971954ec8"
uuid = "65257c39-d410-5151-9873-9b3e5be5013e"
version = "0.4.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SignedDistanceFields]]
deps = ["Random", "Statistics", "Test"]
git-tree-sha1 = "d263a08ec505853a5ff1c1ebde2070419e3f28e9"
uuid = "73760f76-fbc4-59ce-8f25-708e95d2df96"
version = "0.4.0"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.SimpleNonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "ConcreteStructs", "DiffEqBase", "DiffResults", "DifferentiationInterface", "FastClosures", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "MaybeInplace", "PrecompileTools", "Reexport", "SciMLBase", "Setfield", "StaticArraysCore"]
git-tree-sha1 = "4d7a7c177bcb4c6dc465f09db91bfdb28c578919"
uuid = "727e6d20-b764-4bd8-a329-72de5adea6c7"
version = "1.12.0"

    [deps.SimpleNonlinearSolve.extensions]
    SimpleNonlinearSolveChainRulesCoreExt = "ChainRulesCore"
    SimpleNonlinearSolveReverseDiffExt = "ReverseDiff"
    SimpleNonlinearSolveTrackerExt = "Tracker"
    SimpleNonlinearSolveZygoteExt = "Zygote"

    [deps.SimpleNonlinearSolve.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SimpleUnPack]]
git-tree-sha1 = "58e6353e72cde29b90a69527e56df1b5c3d8c437"
uuid = "ce78b400-467f-4804-87d8-8f486da07d0a"
version = "1.1.0"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "2da10356e31327c7096832eb9cd86307a50b1eb6"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.SparseDiffTools]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "Compat", "DataStructures", "FiniteDiff", "ForwardDiff", "Graphs", "LinearAlgebra", "PackageExtensionCompat", "Random", "Reexport", "SciMLOperators", "Setfield", "SparseArrays", "StaticArrayInterface", "StaticArrays", "Tricks", "UnPack", "VertexSafeGraphs"]
git-tree-sha1 = "c9e5d7ee75cf6a1ca3a22c9a6a4ef451792cf62b"
uuid = "47a9eef4-7e08-11e9-0b38-333d64bd3804"
version = "2.20.0"

    [deps.SparseDiffTools.extensions]
    SparseDiffToolsEnzymeExt = "Enzyme"
    SparseDiffToolsPolyesterExt = "Polyester"
    SparseDiffToolsPolyesterForwardDiffExt = "PolyesterForwardDiff"
    SparseDiffToolsSymbolicsExt = "Symbolics"
    SparseDiffToolsZygoteExt = "Zygote"

    [deps.SparseDiffTools.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    Polyester = "f517fe37-dbe3-4b94-8317-1923a5111588"
    PolyesterForwardDiff = "98d1487c-24ca-40b6-b7ab-df2af84e126b"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SparseMatrixColorings]]
deps = ["ADTypes", "Compat", "DataStructures", "DocStringExtensions", "LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "996dff77d814c45c3f2342fa0113e4ad31e712e8"
uuid = "0a514795-09f3-496d-8182-132a7b665d35"
version = "0.4.0"

[[deps.Sparspak]]
deps = ["Libdl", "LinearAlgebra", "Logging", "OffsetArrays", "Printf", "SparseArrays", "Test"]
git-tree-sha1 = "342cf4b449c299d8d1ceaf00b7a49f4fbc7940e7"
uuid = "e56a9233-b9d6-4f03-8d0f-1825330902ac"
version = "0.3.9"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "2f5d4697f21388cbe1ff299430dd169ef97d7e14"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.4.0"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "d2fdac9ff3906e27f7a618d47b676941baa6c80c"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.8.10"

[[deps.StaticArrayInterface]]
deps = ["ArrayInterface", "Compat", "IfElse", "LinearAlgebra", "PrecompileTools", "Static"]
git-tree-sha1 = "96381d50f1ce85f2663584c8e886a6ca97e60554"
uuid = "0d7ed370-da01-4f52-bd93-41d350b8b718"
version = "1.8.0"
weakdeps = ["OffsetArrays", "StaticArrays"]

    [deps.StaticArrayInterface.extensions]
    StaticArrayInterfaceOffsetArraysExt = "OffsetArrays"
    StaticArrayInterfaceStaticArraysExt = "StaticArrays"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "eeafab08ae20c62c44c8399ccb9354a04b80db50"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.7"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "5cf7606d6cef84b543b483848d4ae08ad9832b21"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.3"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "cef0472124fab0695b58ca35a77c6fb942fdab8a"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.1"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.StrideArraysCore]]
deps = ["ArrayInterface", "CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface", "ThreadingUtilities"]
git-tree-sha1 = "f35f6ab602df8413a50c4a25ca14de821e8605fb"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.5.7"

[[deps.StructArrays]]
deps = ["ConstructionBase", "DataAPI", "Tables"]
git-tree-sha1 = "f4dc295e983502292c4c3f951dbb4e985e35b3be"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.18"
weakdeps = ["Adapt", "GPUArraysCore", "SparseArrays", "StaticArrays"]

    [deps.StructArrays.extensions]
    StructArraysAdaptExt = "Adapt"
    StructArraysGPUArraysCoreExt = "GPUArraysCore"
    StructArraysSparseArraysExt = "SparseArrays"
    StructArraysStaticArraysExt = "StaticArrays"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.SymbolicIndexingInterface]]
deps = ["Accessors", "ArrayInterface", "RuntimeGeneratedFunctions", "StaticArraysCore"]
git-tree-sha1 = "c9fce29fb41a10677e24f74421ebe31220b81ad0"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.28"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "598cd7c1f68d1e205689b1c2fe65a9f85846f297"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "eda08f7e9818eb53661b3deb74e3159460dfbc27"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.2"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "ProgressMeter", "SIMD", "UUIDs"]
git-tree-sha1 = "bc7fd5c91041f44636b2c134041f7e5263ce58ae"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.10.0"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "5a13ae8a41237cff5ecf34f73eb1b8f42fff6531"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.24"

[[deps.ToeplitzMatrices]]
deps = ["AbstractFFTs", "DSP", "FillArrays", "LinearAlgebra"]
git-tree-sha1 = "05a042dcb3dabaedb4f1c20de0932c34a0fcee76"
uuid = "c751599d-da0a-543b-9d20-d0a503d91d24"
version = "0.8.4"
weakdeps = ["StatsBase"]

    [deps.ToeplitzMatrices.extensions]
    ToeplitzMatricesStatsBaseExt = "StatsBase"

[[deps.TranscodingStreams]]
git-tree-sha1 = "60df3f8126263c0d6b357b9a1017bb94f53e3582"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.0"
weakdeps = ["Random", "Test"]

    [deps.TranscodingStreams.extensions]
    TestExt = ["Test", "Random"]

[[deps.TriangularSolve]]
deps = ["CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "LoopVectorization", "Polyester", "Static", "VectorizationBase"]
git-tree-sha1 = "be986ad9dac14888ba338c2554dcfec6939e1393"
uuid = "d5829a12-d9aa-46ab-831f-fb7c9ab06edf"
version = "0.2.1"

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

[[deps.TriplotBase]]
git-tree-sha1 = "4d4ed7f294cda19382ff7de4c137d24d16adc89b"
uuid = "981d1d27-644d-49a2-9326-4793e63143c3"
version = "0.1.0"

[[deps.TruncatedStacktraces]]
deps = ["InteractiveUtils", "MacroTools", "Preferences"]
git-tree-sha1 = "ea3e54c2bdde39062abf5a9758a23735558705e1"
uuid = "781d530d-4396-4725-bb49-402e4bee1e77"
version = "1.4.0"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "dd260903fdabea27d9b6021689b3cd5401a57748"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.20.0"
weakdeps = ["ConstructionBase", "InverseFunctions"]

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "975c354fcd5f7e1ddcc1f1a23e6e091d99e99bc8"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.4"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "e7f5b81c65eb858bed630fe006837b935518aca5"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.70"

[[deps.VertexSafeGraphs]]
deps = ["Graphs"]
git-tree-sha1 = "8351f8d73d7e880bfc042a8b6922684ebeafb35c"
uuid = "19fa3120-7c27-5ec5-8db8-b0b0aa330d6f"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "7558e29847e99bc3f04d6569e82d0f5c54460703"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+1"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "93f43ab61b16ddfb2fd3bb13b3ce241cafb0e6c9"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.31.0+0"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c1a7aa6219628fcd757dede0ca95e245c5cd9511"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.0.0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "d9717ce3518dc68a99e6b96300813760d887a01d"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.1+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "a54ee957f4c86b526460a720dbc882fa5edcbefc"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.41+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ac88fb95ae6447c8dda6a5503f3bafd496ae8632"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.4.6+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "326b4fea307b0b39892b3e85fa451692eda8d46c"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.1+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "3796722887072218eabafb494a13c963209754ce"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.4+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "d2d1a5c49fae4ba39983f63de6afcbea47194e85"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.6+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "47e45cd78224c53109495b3e324df0c37bb61fbe"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.11+0"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "bcd466676fef0878338c61e655629fa7bbc69d8e"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.0+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "730eeca102434283c50ccf7d1ecdadf521a765a4"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "04341cb870f29dcd5e39055f895c39d016e18ccd"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.4+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "330f955bc41bb8f5270a369c473fc4a5a4e4d3cb"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e678132f07ddb5bfa46857f0d7620fb9be675d3b"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.6+0"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "gperf_jll"]
git-tree-sha1 = "431b678a28ebb559d224c0b6b6d01afce87c51ba"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.9+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a68c9655fbe6dfcab3d972808f1aafec151ce3f8"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.43.0+0"

[[deps.gperf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3516a5630f741c9eecb3720b1ec9d8edc3ecc033"
uuid = "1a1c6b14-54f6-533d-8383-74cd7377aa70"
version = "3.1.1+0"

[[deps.isoband_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51b5eeb3f98367157a7a12a1fb0aa5328946c03c"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.3+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1827acba325fdcdf1d2647fc8d5301dd9ba43a9d"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.9.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "141fe65dc3efabb0b1d5ba74e91f6ad26f84cc22"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.11.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "ad50e5b90f222cfe78aa3d5183a20a12de1322ce"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.18.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d7015d2e18a5fd9a4f47de711837e980519781a4"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.43+1"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "libpng_jll"]
git-tree-sha1 = "d4f63314c8aa1e48cd22aa0c17ed76cd1ae48c3c"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.3+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "814e154bdb7be91d78b6802843f76b6ece642f11"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.6+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7d0ea0f4895ef2f5cb83645fa689e52cb55cf493"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2021.12.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9c304562909ab2bab0262639bd4f444d7bc2be37"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+1"
"""

# ╔═╡ Cell order:
# ╟─7fc40507-eda3-474d-a454-04e9173a7adb
# ╠═d0e623ee-b096-4a27-977d-dc32567d6020
# ╠═018ecc45-8638-4a59-b561-efb086bdc751
# ╟─584fad48-cadf-4c57-a4c8-4de54f53ebbe
# ╟─d5a510a3-c518-47ed-96bc-7bb22e3b08b5
# ╟─f73b88cb-19e2-4d50-a45e-fbe90fb691cd
# ╟─cf717e9b-26d8-449d-9f3b-d1d5d0bd2bb5
# ╟─c33dc650-3f94-11ef-398a-8bbc4a2b69b8
# ╟─e0c480c1-ab86-4af4-9027-75ebb4fba960
# ╟─1ea9c840-369e-48e6-a87b-7e4bdaa1abb0
# ╟─1bfcc76a-266a-4784-b76b-3871445e603e
# ╟─3e3a3a60-f76f-401d-b265-a99e713d98c5
# ╟─8f3e96dc-0cdc-411e-8f94-e9e47488abf3
# ╟─56bbfe4c-1f28-47e5-9a16-e521f1d55350
# ╟─be501ea1-e0d3-4a4b-8f9e-19135d21aa9d
# ╟─b5d97ade-bcde-4c7a-ae78-697e73c10a0c
# ╟─3ebed1fa-923d-4a0e-9d80-f2919972b1ff
# ╟─62868ea5-3018-4e1a-a71e-89ede6e443c2
# ╟─f5a85447-f9a0-4f56-b2bd-49da4d74326c
# ╟─31104635-eb81-42b0-b412-b26a3635d25f
# ╟─96bb7a72-c10b-45d7-aa64-4b61e488249a
# ╠═3752c9c6-dd41-4194-aaf2-1cd1c2160665
# ╠═dbaf57aa-8a89-476e-b21f-56fe3d8696db
# ╟─97729076-0443-4740-9b0a-6e64cc75b567
# ╟─1c12da5f-7c44-49c2-8077-d93d9de9f32e
# ╟─fc8635b9-2867-4580-a9ca-fd2383c0efd1
# ╟─a1f8e3df-15ff-4c95-b88b-717901b67c1d
# ╟─d9604e75-0e9e-44b6-bf2e-885f01a96bdc
# ╟─62719aa9-bfda-4913-8b74-f3694bfc33c0
# ╟─fc381b29-6671-4364-a9df-b35db07af45b
# ╟─194ee83b-7442-42f6-862a-2f5bd1c8b710
# ╟─9829987c-2bf5-4b33-8ae1-6b80d43ed978
# ╟─1f8fba92-d023-4444-b4e3-18c751f678c5
# ╟─11aa1d55-1a2c-48cf-a7f5-30ff16e2e92a
# ╟─668058ef-e66b-4fd4-ad67-2a8389589e86
# ╟─0562b884-0653-4d5b-8a66-e35467e6693e
# ╟─a47a4e2a-371e-4244-8574-493f00d456df
# ╟─beca7238-ca9b-409a-8d5a-240ca720f58a
# ╠═5e6482db-9b76-4bf0-bc58-0a6f779c4556
# ╟─8c0ac2a2-b3db-455a-8a6d-3de6728b61d3
# ╟─7cfbe32b-d677-4611-97d9-ed775c60e48e
# ╟─13f3b4e4-0409-49a2-a544-0c34ff9289db
# ╟─40156798-59d5-4afe-a054-29c753a51e46
# ╟─42ad70a3-114e-4cd3-8dd1-2429ac80bfc2
# ╠═cdbbb2e1-56bc-4beb-93d3-2db1a2aa71d5
# ╟─aaf278f6-b8ee-4094-9173-ffde2cb0b461
# ╠═d77a192e-4be0-4987-8023-8d91b22d8d9f
# ╟─7410e377-715c-4464-a83d-d7cd48c4d1ae
# ╠═fbeac0c7-8fcb-4e4b-9662-9f1d48b23fb3
# ╟─ae88c826-f486-4786-9447-affcc466d1b5
# ╟─76aa085c-d7ab-4e2f-9648-89b3eebe11cd
# ╟─585b99f6-fad0-44b1-b1ca-f7aa6f8f52f8
# ╟─400fe5e0-0dfd-4d82-baf9-b66d931133a8
# ╠═e0f75d0e-555d-4df7-8077-9871a6d7b197
# ╟─deb48c67-19f1-481e-b2de-1765e449e6ae
# ╠═c106822b-2c0c-4759-84fb-9ddba91fe295
# ╟─f75c26b4-8e6f-4520-b1a8-f30c20b9a194
# ╠═05a4dd64-b5c9-4014-a18c-f26938e02ebd
# ╟─3cb601bf-43ea-430c-950d-402b0929e515
# ╠═2a76f46f-c203-4028-80ba-e7bc1f481076
# ╟─8a802d5f-e2f5-4e03-b0d9-b7cd0c081d2c
# ╠═5ef06ffe-d2f3-4cd2-946a-df8fd78f863f
# ╟─e4a015b0-e7e6-4841-b97f-4a845f462196
# ╠═c3d0d2ae-5711-451b-946d-6cdca682fe06
# ╟─4ff1b050-3036-4e0b-b9db-b07b4d17327c
# ╠═ec3cd93f-877d-419e-a608-cef748306927
# ╟─6e715b86-4604-49c3-8df1-929921e38f9a
# ╠═ffe5218b-88ca-42e0-8a0d-074ab2e6d898
# ╟─98c9ae61-0a84-48c0-8773-0fe9bbd260ac
# ╠═5838bc82-6896-49a7-bc2c-b700fb3a3eb8
# ╟─0bff34a8-2ad3-43d8-9f25-60a926f5b5c6
# ╟─3c23bf6e-fd79-4770-83c3-385832eddac6
# ╟─c4363fdf-6500-4830-bc42-9777d616ac34
# ╟─e72d43a1-b804-455a-bff0-8b1b715ddcee
# ╟─24e1a5af-569e-4ae8-a3d7-19b7962db48b
# ╟─922f5fa4-f9db-4eae-bd81-dc3dd525d212
# ╠═bb5ea822-f29b-499c-bebb-53208cf42eb1
# ╟─0d75d077-53c3-4f03-99d3-af796d973d66
# ╟─9433b6f1-3460-48c4-b587-9009823327c0
# ╟─4741e38b-0c29-4f02-a540-c543f7b3f234
# ╠═656ae53a-1dbc-4abc-8837-234762db80c0
# ╟─9fe5120a-824a-4395-afe8-325a36ab9a07
# ╟─eb2ae1fd-5062-436b-ac50-eb3b169bca42
# ╠═3021e7c0-508a-4a95-9c16-098fa35b65a0
# ╟─fc18e3ff-4a39-4aba-a2d2-94521035013e
# ╟─40a4ef18-8cf6-46bf-bf84-a7c89d898f51
# ╟─454efc82-4918-4c60-8b9d-fb630302beec
# ╟─c52ec039-8ef2-41fc-a9cd-34fcef45f735
# ╟─7229c43e-dcda-436c-86fc-1e626fe43daf
# ╟─bc448153-fe09-4a05-b857-5b9fd51dcdae
# ╟─1d3d35a7-c7a8-419e-b751-a144b86bbe1b
# ╟─002c7c0a-2b38-4a31-9e08-a424fa5679b8
# ╟─9f022768-29a7-4ce9-bf9a-ceb04ea0adf2
# ╟─090c74fc-9bf7-47c0-af5f-887d6bf43ccd
# ╟─9e8d66d9-ad6e-42b5-a24f-6d0a9060f688
# ╟─6297e819-a105-4451-bbb7-c2241ff77f51
# ╟─fb8527a5-4fa9-40d3-b0af-2cce8cacac3d
# ╟─7e60324b-f65e-4d08-b642-13542b3637f3
# ╠═c5d63b4c-7823-4122-a81f-9e3a18a634ff
# ╟─c0ce2c85-eaae-4b08-8ef7-2a9297b10e05
# ╟─6f346c67-879a-4427-afd4-6acdeb17af21
# ╠═6688fe5d-0019-485d-9059-bb5bd2be7d6e
# ╟─5d0b6b50-0170-4f8b-ad5a-0dda212a00e1
# ╟─a9f597fe-4058-47a8-8a39-8d6b36594a96
# ╠═1082f423-448c-4d74-8baa-76b2e6508084
# ╟─f82c6a7a-fffb-4a7b-97df-bfafd2072f5c
# ╟─9a0d1fd7-e8b6-4737-b764-95f917442335
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

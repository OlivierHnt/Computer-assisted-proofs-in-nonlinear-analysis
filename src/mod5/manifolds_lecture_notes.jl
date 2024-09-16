### A Pluto.jl notebook ###
# v0.19.46

#> [frontmatter]
#> chapter = 5
#> order = 1
#> title = "Lecture notes: invariant manifolds of equilibria"
#> tags = ["module5"]
#> layout = "layout.jlhtml"

using Markdown
using InteractiveUtils

# ╔═╡ d0e623ee-b096-4a27-977d-dc32567d6020
using PlutoTeachingTools, PlutoUI # packages for the notebook

# ╔═╡ 018ecc45-8638-4a59-b561-efb086bdc751
using RadiiPolynomial, LinearAlgebra

# ╔═╡ 7fc40507-eda3-474d-a454-04e9173a7adb
html"""<style>
main {
    margin-left: auto;
    margin-right: auto;
    text-align: justify;
}
"""

# ╔═╡ 4412bd33-8b3f-41de-9e1a-a8f1187c793f
TableOfContents(title="Table of Contents"; indent=true, depth=4, aside=true)

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

# ╔═╡ 0bff34a8-2ad3-43d8-9f25-60a926f5b5c6
# md"""
# ##### Plot the image of the parameterization
# """

# ╔═╡ 3c23bf6e-fd79-4770-83c3-385832eddac6
# begin
# 	using OrdinaryDiffEq
# 	using GLMakie
# 	function lorenz(du, u, p, t)
# 	    du[1] = 10.0(u[2] - u[1])
# 	    du[2] = u[1] * (28.0 - u[3]) - u[2]
# 	    du[3] = u[1] * u[2] - (8 / 3) * u[3]
# 	end
# 	u0 = [0.0; 10.0; 10.0]
# 	tspan = (0.0, 50.0)
# 	prob = ODEProblem(lorenz, u0, tspan)
# 	sol = solve(prob, Tsit5())
# 	fig = Figure()
# 	ax = Axis3(fig[1,1]; aspect = :data, azimuth = 3π/5)
# 	lines!(ax, [Point3f(sol(t)) for t = LinRange(tspan..., 10_001)]; color = :royalblue1)
# 	lines!(ax, [Point3f(a(t)) for t = LinRange(-1, 1, 101)]; color = :salmon, linewidth = 3)
# 	meshscatter!(ax, [Point3f(bx), Point3f(-bx[1], -bx[2], bx[3]), Point3f(0,0,0)]; markersize = 1, color = :palegreen3)
# 	fig
# end

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
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
RadiiPolynomial = "f2081a94-c849-46b6-8dc9-07bb90ed72a9"

[compat]
PlutoTeachingTools = "~0.2.15"
PlutoUI = "~0.7.60"
RadiiPolynomial = "~0.8.13"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.5"
manifest_format = "2.0"
project_hash = "4f12b0bdeca87c29c864818987e167b160e4de90"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.CRlibm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e329286945d0cfc04456972ea732551869af1cfc"
uuid = "4e9b3aee-d8a1-5a3d-ad8b-7d824db253f0"
version = "1.0.1+0"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "7eee164f122511d3e4e1ebadb7956939ea7e1c77"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.3.6"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

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

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IntervalArithmetic]]
deps = ["CRlibm_jll", "MacroTools", "RoundingEmulator"]
git-tree-sha1 = "fe30dec78e68f27fc416901629c6e24e9d5f057b"
uuid = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
version = "0.22.16"

    [deps.IntervalArithmetic.extensions]
    IntervalArithmeticDiffRulesExt = "DiffRules"
    IntervalArithmeticForwardDiffExt = "ForwardDiff"
    IntervalArithmeticIntervalSetsExt = "IntervalSets"
    IntervalArithmeticLinearAlgebraExt = "LinearAlgebra"
    IntervalArithmeticRecipesBaseExt = "RecipesBase"

    [deps.IntervalArithmetic.weakdeps]
    DiffRules = "b552c78f-8df3-52c6-915a-8e097449b14b"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "f389674c99bfcde17dc57454011aa44d5a260a40"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.6.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "2984284a8abcfcc4784d95a9e2ea4e352dd8ede7"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.36"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "ce5f5621cac23a86011836badfedf664a612cee4"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.5"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

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

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "c2b5e92eaf5101404a58ce9c6083d595472361d6"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "3.0.2"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

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
git-tree-sha1 = "eba4810d5e6a01f612b948c9fa94f905b49087b0"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.60"

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

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.RadiiPolynomial]]
deps = ["IntervalArithmetic", "LinearAlgebra", "Printf", "Reexport", "SparseArrays"]
git-tree-sha1 = "8442e84088a316034b2b9d8128d6af0ac4ab4fad"
uuid = "f2081a94-c849-46b6-8dc9-07bb90ed72a9"
version = "0.8.13"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Revise]]
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "7b7850bb94f75762d567834d7e9802fc22d62f9c"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.5.18"

[[deps.RoundingEmulator]]
git-tree-sha1 = "40b9edad2e5287e05bd413a38f61a8ff55b9557b"
uuid = "5eaf0fd0-dfba-4ccb-bf02-d820a40db705"
version = "0.2.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Tricks]]
git-tree-sha1 = "7822b97e99a1672bfb1b49b668a6d46d58d8cbcb"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.9"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╟─7fc40507-eda3-474d-a454-04e9173a7adb
# ╠═d0e623ee-b096-4a27-977d-dc32567d6020
# ╠═018ecc45-8638-4a59-b561-efb086bdc751
# ╟─4412bd33-8b3f-41de-9e1a-a8f1187c793f
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

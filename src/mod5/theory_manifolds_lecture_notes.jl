### A Pluto.jl notebook ###
# v0.19.46

#> [frontmatter]
#> chapter = 5
#> order = 1
#> title = "Parameterization method"
#> tags = ["lecture", "module5"]
#> layout = "layout.jlhtml"

using Markdown
using InteractiveUtils

# ╔═╡ 12bc1703-41cc-44d2-8b56-cd42012cc8e1
using PlutoTeachingTools, PlutoUI # packages for the notebook

# ╔═╡ fc651998-8263-475c-aece-4e257cd6be70
using RadiiPolynomial, LinearAlgebra

# ╔═╡ be25c236-6e65-47c7-a36c-c0a6eb594f57
html"""
<style>
main {
    margin-right: auto;
    text-align: justify;
}
</style>
"""

# ╔═╡ 5482d7a2-df89-4ca9-b3f7-5fac7967f6e1
TableOfContents(title = "Table of Contents"; indent = true, depth = 4, aside = true)

# ╔═╡ e507fbe8-e81c-4386-b33b-b0482815f8c1
md"""
# Local stable manifolds of equilibria
"""

# ╔═╡ 5c7f76e1-ba35-4cb7-9fdc-cfbf05a9ab13
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

# ╔═╡ 51622c96-c674-4079-ad9d-3554b6562ed4
md"""
## The Parameterization Method
"""

# ╔═╡ a707c66a-8e9f-454f-a0ad-bcc4d628d82d
md"""
### Assumptions
"""

# ╔═╡ 1e145c27-90cd-43f1-b8b2-4ae19b8b9ec3
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

# ╔═╡ aa2f4387-e3d9-48b1-bee1-af8c80dade1f
md"""
### The conjugacy relation
"""

# ╔═╡ c60b9061-3143-46b4-997f-c5c072b81da5
md"""
Denote by $B^m := \{ \theta = (\theta_1,\dots,\theta_m) : |\theta_j| \le 1 \}$ the unit closed $\ell^\infty$ ball in $\mathbb{R}^m$.

To solve for the parameterization $P:B^m \to \mathbb{R}^n$ of the local stable manifold, we impose

(1) $P(0) = \tilde{x}$ (equilibria are mapped to one another)

(2) $DP(0) = A_0$ ($P$ is tangent to the stable eigenspace $\mathbb{E}^s$ at $0$})

(3) $\varphi (t, P(\theta)) = P( e^{\Lambda t} \theta)$, for all $\theta \in B^m$ (**conjugacy relation**)
"""

# ╔═╡ e4163286-ceda-4277-8da9-5a9394b15e54
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

# ╔═╡ 6d00a202-4d6b-497e-b423-76e3323c4f24
md"""
### The parameterization lemma
"""

# ╔═╡ 07dcb89c-9926-432e-94d4-44cb98c2c8a1
md"""

*Remark.* The utility of (3) above is limited by the appearance of the flow $\varphi$, which is only known implicitly. Let us introduce a more practical infinitesimal version of the conjugacy relation, which requires the following important result."""

# ╔═╡ 68c49689-e740-4ddb-bc67-6f811e8b891d
Markdown.MD(Markdown.Admonition("tip", "Parameterization Lemma",[md"""
Let $P \colon B^m  \to \mathbb{R}^n$ be a function satisfying (1) and (2). Then $P$ satisfies (3) $\Longleftrightarrow$ $P$ solves the *partial* differential equation

$\begin{align}
\text{(4)} \quad DP(\theta) \Lambda \theta = f(P(\theta)), \quad \theta \in B^m.
\end{align}$"""]))

# ╔═╡ 4d46da29-6ec9-4d0b-b28e-843c1774b2c4
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

# ╔═╡ 4e9479fe-2ddd-4ea5-a8dd-93ef889aef33
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

# ╔═╡ a85db154-e650-40bd-8191-175a28a5de02
md"""
### Taylor series representation
"""

# ╔═╡ ac318edc-cca6-456a-b4ac-7ece084089c5
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

# ╔═╡ da95dea7-9ec5-46f4-beea-9a22cf5ca6ee
md"""
### The Zero Finding Problem, Banach Space and Norm
"""

# ╔═╡ 510a3921-abe3-414f-ad17-c165415d8e40
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

# ╔═╡ 32179842-b1fc-4cbb-bb98-a39b6ebce9ee
md"""
###### Code that implements the operator $K$
"""

# ╔═╡ f1bb1b2a-c2b5-48a3-897a-bcee907333e3
function K(λ, N)
    op = zeros(typeof(λ), Taylor(N), Taylor(N))
    for n = 2:N
        op[n,n] = inv(λ * n)
    end
    return op
end

# ╔═╡ ccb1e680-934b-46b0-ad8f-4bd9516010c6
K(-2.0,5)

# ╔═╡ 3cfbcab7-c3cc-45ce-99cc-b9e800bec98b
Markdown.MD(Markdown.Admonition("note", "Exercise (Compact operator)",[md"""
Show that the operator $\mathcal{K}:\mathcal{X} \to \mathcal{X}$ is compact.
"""]))

# ╔═╡ b9ad918e-6b48-4c14-a2cc-ead87c667070
md"""
We now formalize how a solution to the above zero finding problem leads to local stable manifolds of equilibria of our differential equation.
"""

# ╔═╡ 0548b47b-cc62-4730-8293-26b2afe3f793
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

# ╔═╡ ce8ea7a5-dea3-4cda-9dea-eab93c54b24d
md"""
Rather than presenting the bounds for the Newton-Kantorovich Theorem in general, let us rather focus on an example. Specifically, we compute a rigorous parameterization of a one-dimensional local stable manifold of an equilibrium solution in the Lorenz system.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
RadiiPolynomial = "f2081a94-c849-46b6-8dc9-07bb90ed72a9"

[compat]
PlutoTeachingTools = "~0.3.0"
PlutoUI = "~0.7.60"
RadiiPolynomial = "~0.8.13"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.5"
manifest_format = "2.0"
project_hash = "e8d35381806b79379583eb0dc6132f4e83486b98"

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
deps = ["Downloads", "HypertextLiteral", "Latexify", "Markdown", "PlutoLinks", "PlutoUI"]
git-tree-sha1 = "e2593782a6b53dc5176058d27e20387a0576a59e"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.3.0"

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
# ╠═12bc1703-41cc-44d2-8b56-cd42012cc8e1
# ╠═fc651998-8263-475c-aece-4e257cd6be70
# ╟─be25c236-6e65-47c7-a36c-c0a6eb594f57
# ╟─5482d7a2-df89-4ca9-b3f7-5fac7967f6e1
# ╟─e507fbe8-e81c-4386-b33b-b0482815f8c1
# ╟─5c7f76e1-ba35-4cb7-9fdc-cfbf05a9ab13
# ╟─51622c96-c674-4079-ad9d-3554b6562ed4
# ╟─a707c66a-8e9f-454f-a0ad-bcc4d628d82d
# ╟─1e145c27-90cd-43f1-b8b2-4ae19b8b9ec3
# ╟─aa2f4387-e3d9-48b1-bee1-af8c80dade1f
# ╟─c60b9061-3143-46b4-997f-c5c072b81da5
# ╟─e4163286-ceda-4277-8da9-5a9394b15e54
# ╟─6d00a202-4d6b-497e-b423-76e3323c4f24
# ╟─07dcb89c-9926-432e-94d4-44cb98c2c8a1
# ╟─68c49689-e740-4ddb-bc67-6f811e8b891d
# ╟─4d46da29-6ec9-4d0b-b28e-843c1774b2c4
# ╟─4e9479fe-2ddd-4ea5-a8dd-93ef889aef33
# ╟─a85db154-e650-40bd-8191-175a28a5de02
# ╟─ac318edc-cca6-456a-b4ac-7ece084089c5
# ╟─da95dea7-9ec5-46f4-beea-9a22cf5ca6ee
# ╟─510a3921-abe3-414f-ad17-c165415d8e40
# ╟─32179842-b1fc-4cbb-bb98-a39b6ebce9ee
# ╠═f1bb1b2a-c2b5-48a3-897a-bcee907333e3
# ╠═ccb1e680-934b-46b0-ad8f-4bd9516010c6
# ╟─3cfbcab7-c3cc-45ce-99cc-b9e800bec98b
# ╟─b9ad918e-6b48-4c14-a2cc-ead87c667070
# ╟─0548b47b-cc62-4730-8293-26b2afe3f793
# ╟─ce8ea7a5-dea3-4cda-9dea-eab93c54b24d
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

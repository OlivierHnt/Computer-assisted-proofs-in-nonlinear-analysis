### A Pluto.jl notebook ###
# v0.20.1

#> [frontmatter]
#> homework_number = 3
#> order = 3
#> title = "Roots"
#> tags = ["module2", "homeworks"]
#> layout = "layout.jlhtml"

using Markdown
using InteractiveUtils

# ╔═╡ 7fc40507-eda3-474d-a454-04e9173a7adb
html"""
<style>
main {
    max-width: 1000px;
    margin-left: auto;
    margin-right: auto;
    text-align: justify;
}
</style>
"""

# ╔═╡ e2ba1e46-84af-4deb-8042-43d9b7c78194
md"""
In this exercise, we will formulate a zero-finding problem to compute rigorously the square root of a Fourier series.

Similarly to the case of Taylor series, we discretize the function space of periodic functions via hte bi-sided sequence space

```math
\ell^1_{\mathbb{Z}} \overset{\text{def}}{=} \left\{ a \in \mathbb{C}^\mathbb{Z} \, : \, \| a \|_1 \overset{\text{def}}{=} \sum_{k \in \mathbb{Z}} |a_k| < \infty \right\}.
```
"""

# ╔═╡ 3819f407-5f4f-4c2a-bec9-fe0409384ffb
md"""
**1.** Prove that if $A \in B(\ell^1)$, then $\| A \|_{B(\ell^1)} \le \max_{l \in \mathbb{Z}} \sum_{k \in \mathbb{Z}} |A_{k,l}|$ (in fact the inequality is an equality).
"""

# ╔═╡ 4c9f56dd-b79d-4e7e-8b1d-38e91f6632df
md"""
**2.** If $u(t) = \sum_{k \in \mathbb{Z}} a_k e^{ikt}$ and $v(t) = \sum_{k \in \mathbb{Z}} b_k e^{ikt}$, then for their product we have $u(t) v(t) = \sum_{k \in \mathbb{Z}} (a*b)_k e^{ikt}$ where

```math
(a*b)_k \overset{\text{def}}{=} \sum_{l \in \mathbb{Z}} a_{k-l} b_l, \qquad k \in \mathbb{Z}.
```

Prove that $\ell^1_{\mathbb{Z}}$ together with $*$ forms a unital Banach algebra, i.e. $\|a*b\|_1 \le \|a\|_1 \|b\|_1$.
"""

# ╔═╡ 2c0bf040-7984-49ca-8860-bf090631ff22
md"""
**3.** Given $b \in \ell^1_{\mathbb{Z}}$, consider $F : \ell^1_{\mathbb{Z}} \to \ell^1_{\mathbb{Z}}$ given by

```math
F(a) \overset{\text{def}}{=} a*a - b.
```

Noting that the multiplication operator associated with $b \in \ell^1_\mathbb{Z}$ has the structure of a [Toepliz matrix](https://en.wikipedia.org/wiki/Toeplitz_matrix), specifically

```math
M_b =
\begin{pmatrix}
\ddots & \ddots & \ddots \\
\ddots & b_0    & b_{-1} & b_{-2} \\
\ddots & b_1    & b_0    & b_{-1} & \ddots \\
       & b_2    & b_1    & b_0    & \ddots \\
       &        & \ddots & \ddots & \ddots
\end{pmatrix},
```

use the Newton-Kantorovich theorem to obtain a rigorous enclosure of the square root of $v(t) = 2 + \cos(t)$.
"""

# ╔═╡ 15be0e8a-408b-4db2-af7e-15261f54238e
md"""
**4.** Define a suitable $F = 0$ problem for the cubic root of a function using $\ell^1_\mathbb{Z}$, and derive the bounds needed to apply the Newton-Kantorovich theorem in that context.
"""

# ╔═╡ 4ea373ca-e44b-49cf-9e6c-74a5043dee79
md"""
**5.** Find the cubic root of $v(t) = 2 + \cos (t)$.
"""

# ╔═╡ Cell order:
# ╟─7fc40507-eda3-474d-a454-04e9173a7adb
# ╟─e2ba1e46-84af-4deb-8042-43d9b7c78194
# ╟─3819f407-5f4f-4c2a-bec9-fe0409384ffb
# ╟─4c9f56dd-b79d-4e7e-8b1d-38e91f6632df
# ╟─2c0bf040-7984-49ca-8860-bf090631ff22
# ╟─15be0e8a-408b-4db2-af7e-15261f54238e
# ╟─4ea373ca-e44b-49cf-9e6c-74a5043dee79

### A Pluto.jl notebook ###
# v0.19.47

#> [frontmatter]
#> homework_number = 4
#> order = 4
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

Similarly to the case of Taylor series, we discretize the function space of (absolutely continuous) periodic functions via the two-sided sequence space

```math
\ell^1_{\mathbb{Z}} \overset{\text{def}}{=} \left\{ x \in \mathbb{C}^\mathbb{Z} \, : \, \| x \|_1 \overset{\text{def}}{=} \sum_{k \in \mathbb{Z}} |x_k| < \infty \right\}.
```
"""

# ╔═╡ 4c9f56dd-b79d-4e7e-8b1d-38e91f6632df
md"""
**1.** Let $\omega > 0$.
If $u(t) = \sum_{k \in \mathbb{Z}} x_k e^{i k \omega t}$ and $v(t) = \sum_{k \in \mathbb{Z}} y_k e^{i k \omega t}$, then $u(t) v(t) = \sum_{k \in \mathbb{Z}} (x*y)_k e^{i k \omega t}$ where

```math
(x*y)_k \overset{\text{def}}{=} \sum_{l \in \mathbb{Z}} x_{k-l} y_l, \qquad k \in \mathbb{Z}.
```

Prove that $\ell^1_{\mathbb{Z}}$ together with $*$ satisfies the unital Banach algebra property $\|x * y\|_1 \le \|x\|_1 \|y\|_1$.
"""

# ╔═╡ 3819f407-5f4f-4c2a-bec9-fe0409384ffb
md"""
**2.** Show that $\| M \|_{\mathscr{B}(\ell^1_\mathbb{Z})} = \max_{l \in \mathbb{Z}} \sum_{k \in \mathbb{Z}} |M_{k,l}|$ for all $M \in \mathscr{B}(\ell^1_\mathbb{Z})$.

Moreover, note that the multiplication operator $\mathcal{M}_x : \ell^1_\mathbb{Z} \to \ell^1_\mathbb{Z}$ associated with $x \in \ell^1_\mathbb{Z}$ has the structure of a [Toepliz matrix](https://en.wikipedia.org/wiki/Toeplitz_matrix), specifically

```math
\mathcal{M}_x =
\begin{pmatrix}
\ddots & \ddots & \ddots \\
\ddots & x_0    & x_{-1} & x_{-2} \\
\ddots & x_1    & x_0    & x_{-1} & \ddots \\
       & x_2    & x_1    & x_0    & \ddots \\
       &        & \ddots & \ddots & \ddots
\end{pmatrix},
```

For all $x \in \ell^1_\mathbb{Z}$, show that $\| \mathcal{M}_x \|_{\mathscr{B}(\ell^1_\mathbb{Z})} = \| x \|_1$.
"""

# ╔═╡ 2c0bf040-7984-49ca-8860-bf090631ff22
md"""
**3.** Given $y \in \ell^1_\mathbb{Z}$, consider $F : \ell^1_\mathbb{Z} \to \ell^1_\mathbb{Z}$ given by

```math
F(x) \overset{\text{def}}{=} x*x - y.
```

Derive the bounds needed to apply the Newton-Kantorovich theorem.
Find the square root of $v(t) = 2 + \cos(t)$.
Why is the obtained Fourier series real-valued?
"""

# ╔═╡ 15be0e8a-408b-4db2-af7e-15261f54238e
md"""
**4.** Define a suitable $F = 0$ problem for the cubic root of an absolutely converging Fourier series, and derive the bounds needed to apply the Newton-Kantorovich theorem in that context.
Find the cubic root of $v(t) = 2 + \cos(t)$.
Why is the obtained Fourier series real-valued?
"""

# ╔═╡ Cell order:
# ╟─7fc40507-eda3-474d-a454-04e9173a7adb
# ╟─e2ba1e46-84af-4deb-8042-43d9b7c78194
# ╟─4c9f56dd-b79d-4e7e-8b1d-38e91f6632df
# ╟─3819f407-5f4f-4c2a-bec9-fe0409384ffb
# ╟─2c0bf040-7984-49ca-8860-bf090631ff22
# ╟─15be0e8a-408b-4db2-af7e-15261f54238e

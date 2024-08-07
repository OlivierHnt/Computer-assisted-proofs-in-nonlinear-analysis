### A Pluto.jl notebook ###
# v0.19.45

#> [frontmatter]
#> homework_number = 1
#> order = 1.5
#> title = "Operator norms"
#> tags = ["module1", "homeworks"]
#> layout = "layout.jlhtml"

using Markdown
using InteractiveUtils

# ╔═╡ 7fc40507-eda3-474d-a454-04e9173a7adb
html"""<style>
main {
    max-width: 1000px;
    margin-left: auto;
    margin-right: auto;
    text-align: justify;
}
"""

# ╔═╡ 15be0e8a-408b-4db2-af7e-15261f54238e
md"""
**1.** Consider the 1-norm on $\mathbb{R}^d$: $\Vert x\Vert_1 = \sum_{i=1}^d \vert x_i\vert$, and let $A$ be a $d\times d$ matrix of real numbers. Show that the operator norm induced by the 1-norm,

$\begin{align}
\Vert A\Vert_1 := \sup_{\substack{x\in\mathbb{R}^d \\ x\neq 0}} \frac{\Vert Ax\Vert_1}{\Vert x\Vert_1},
\end{align}$

can be expressed as

$\begin{align}
\left\Vert A\right\Vert_1 = \max_{1\leq j\leq d} \sum_{i=1}^d \vert A_{i,j}\vert.
\end{align}$
"""

# ╔═╡ 24407233-72ba-4341-84ec-c1a1148c3ff0
md"""
**2.** Let us generalize this result to infinite dimensional spaces (this will be useful for later parts of the tutorial). Consider a sequence of positive weights $\omega_k$, $k\in\mathbb{N}$, and the sequence space

$\begin{align}
\ell^1_\omega := \left\{ x \in \mathbb{R}^{\mathbb{N}},\ \Vert x\Vert_\omega := \sum_{i\in\mathbb{N}} \vert x_i\vert \omega_i < \infty \right\}.
\end{align}$

Given an *infinite matrix* $A=\left(A_{i,j}\right)_{i,j\in\mathbb{N}}$ representing a linear operator on $\ell^1_\omega$, derive a formula for the operator norm of $A$.
"""

# ╔═╡ Cell order:
# ╟─7fc40507-eda3-474d-a454-04e9173a7adb
# ╟─15be0e8a-408b-4db2-af7e-15261f54238e
# ╟─24407233-72ba-4341-84ec-c1a1148c3ff0

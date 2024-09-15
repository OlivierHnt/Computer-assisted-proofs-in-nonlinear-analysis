### A Pluto.jl notebook ###
# v0.19.40

#> [frontmatter]
#> homework_number = 1
#> order = 1
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
Consider the $1$-norm on $\mathbb{R}^d$: $\Vert x \Vert_1 := \sum_{i=1}^d \vert x_i \vert$, and let $A$ be a $d \times d$ matrix of real numbers. Show that the operator norm induced by the $1$-norm,

$\begin{align}
\Vert A\Vert_1 := \sup_{\substack{x\in\mathbb{R}^d \\ x\neq 0}} \frac{\Vert Ax\Vert_1}{\Vert x\Vert_1},
\end{align}$

can be expressed as

$\begin{align}
\left\Vert A\right\Vert_1 = \max_{1\leq j\leq d} \sum_{i=1}^d \vert A_{i,j}\vert.
\end{align}$
"""

# ╔═╡ Cell order:
# ╟─7fc40507-eda3-474d-a454-04e9173a7adb
# ╟─15be0e8a-408b-4db2-af7e-15261f54238e

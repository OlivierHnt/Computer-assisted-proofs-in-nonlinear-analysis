### A Pluto.jl notebook ###
# v0.20.1

#> [frontmatter]
#> homework_number = 1
#> order = 1
#> title = "Operator norms"
#> tags = ["module1", "homeworks"]
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

# ╔═╡ 15be0e8a-408b-4db2-af7e-15261f54238e
md"""
Consider the $1$-norm on $\mathbb{R}^d$

```math
\| x \|_1 \overset{\text{def}}{=} \sum_{i=1}^d | x_i |,
```

and let $A$ be a $d \times d$ matrix of real numbers.

Show that the operator norm induced by the $1$-norm, that is

```math
\| A \|_1 \overset{\text{def}}{=} \sup_{\substack{x \in \mathbb{R}^d \\ x\ne 0}} \frac{\| Ax \|_1}{\| x \|_1},
```

satisfies

```math
\| A \|_1 = \max_{1 \le j \le d} \sum_{i=1}^d | A_{i,j} |.
```
"""

# ╔═╡ Cell order:
# ╟─7fc40507-eda3-474d-a454-04e9173a7adb
# ╟─15be0e8a-408b-4db2-af7e-15261f54238e

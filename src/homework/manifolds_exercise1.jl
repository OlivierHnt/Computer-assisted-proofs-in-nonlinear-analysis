### A Pluto.jl notebook ###
# v0.19.47

#> [frontmatter]
#> homework_number = 1
#> order = 1
#> title = "Stable and unstable manifolds"
#> tags = ["module5", "homeworks"]
#> layout = "layout.jlhtml"

using Markdown
using InteractiveUtils

# ╔═╡ 5c656cec-3f36-434d-a85e-2af114df6483
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

# ╔═╡ 372f61b4-7842-11ef-3b13-f5944a4e171d
md"""
Using Taylor series, compute rigorously (that is with the Newton-Kantorovich approach) a parameterization of order $N = 20$ of the stable and unstable manifolds of the equilibrium $(0,0)$ for the equation

```math
\dot x = \begin{pmatrix} \dot x_1  \\ \dot x_2 \end{pmatrix} =
\begin{pmatrix}
x_1+6 x_2 + x_1 x_2
\\
4 x_1 + 3 x_2 - x_1^2
\end{pmatrix}.
```
"""

# ╔═╡ Cell order:
# ╟─5c656cec-3f36-434d-a85e-2af114df6483
# ╟─372f61b4-7842-11ef-3b13-f5944a4e171d

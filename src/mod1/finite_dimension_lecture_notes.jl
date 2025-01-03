### A Pluto.jl notebook ###
# v0.19.47

#> [frontmatter]
#> chapter = 1
#> order = 0.5
#> title = "Finite-dimensional problems"
#> tags = ["lecture", "module1"]
#> layout = "layout.jlhtml"

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ a5169b4b-f0d4-4f66-8198-0aec5f8e135a
html"""
<style>
main {
    margin-right: auto;
    text-align: justify;
}
</style>
"""

# ╔═╡ d0e623ee-b096-4a27-977d-dc32567d6020
using PlutoTeachingTools, PlutoUI # packages for the notebook

# ╔═╡ 2661bfc9-e398-41ed-87d9-c78f05da64cb
using Plots

# ╔═╡ b3845641-1537-4a27-8550-1eff30900a6b
TableOfContents(title = "Table of Contents"; indent = true, depth = 4, aside = true)

# ╔═╡ ed6f44f4-3696-494a-b8e3-b30184f5bb06
md"# Motivating example: a chaotic dynamical system"

# ╔═╡ c33dc650-3f94-11ef-398a-8bbc4a2b69b8
md"""
Consider the sequence defined by

$\begin{align}
x_{n+1} = \mu x_n (1-x_n),
\end{align}$

where $\mu$ is a parameter in $[0,4]$ and the initial condition $x_1$ is in $[0,1]$.
This is known as the [logistic map](https://en.wikipedia.org/wiki/Logistic_map).

The slider below can be used to vary the value of $\mu$.
"""

# ╔═╡ 9e0bfed7-7ccc-42b2-9b1a-15420d85c356
@bind μ Slider(0:0.1:4; default = 3.9)

# ╔═╡ e3ca081f-0e6a-4dcf-9436-ca4f34467379
μ

# ╔═╡ 79138e2e-6123-4923-8329-e5c17499d785
begin
	x1 = 0.4
	maxiter = 50
	xs = zeros(maxiter+1)
	xs[1] = x1
	for n = 1:maxiter
		xs[n+1] = μ*xs[n]*(1-xs[n])
	end
	plot(1:maxiter+1, xs; marker = (:circle, 5), legend = false)
	xlabel!("n")
	ylims!(-0.1, 1.1)
end

# ╔═╡ 1d9286fa-647a-4743-806f-9cea5aab6be1
md"""
As $\mu$ increases, the dynamics become progressively more complex.
This behaviour arises from a cascade of [period doubling bifurcations](https://en.wikipedia.org/wiki/Period-doubling_bifurcation), ultimately leading to [chaos](https://en.wikipedia.org/wiki/Chaos_theory).
A well-known result for discrete dynamical systems states that the existence of a period-3 orbit implies chaos.
"""

# ╔═╡ 1e56532e-9886-4a41-9ba0-1a619a85e014
md"""
!!! theorem "Theorem"
	Let $f : [0,1] \to [0,1]$ be continuous, and consider the dynamical system defined by $x_{n+1} = f(x_n)$.
	If there exists an orbit of period-3, i.e. $x_1 \in [0,1]$ such that $x_1 \ne x_2 \ne x_3$ and $x_4 = x_1$, then the system is chaotic.
	In particular, there exist orbits of any period.[^Sha64][^LY75]
"""

# ╔═╡ 8d6c631f-832a-48a9-891f-650f53a0ac6c
md"""
Let us fix a value of $\mu$ for which the dynamics seem complicated, say $\mu = 3.9$.
Our goal is to demonstrate the existence of a period-3 orbit.
We will begin by numerically approximating such an orbit and then rigorously establish the existence of a true period-3 orbit nearby.

While there are various methods to prove the existence of a period-3 orbit (both with and without computer assistance), the approach we present will extend naturally to more challenging infinite-dimensional problems.
"""

# ╔═╡ 429e5261-d7d0-47a3-94b5-8e7905a5bd53
md"# Contraction mapping and the Newton-Kantorovich theorem"

# ╔═╡ a27b45ca-77d1-408e-98e9-6b5248aca85e
md"""
Given a problem and an approximate solution $\bar{x}$, our goal is to define a fixed-point operator $T$ such that
- the solutions to our problem are in one-to-one correspondence with fixed-points of $T$,
-  $T$ is a contraction on a small neighborhood of $\bar{x}$.

The specific form of $T$ may vary depending on the nature of the problem, but we can outline a general strategy for constructing such an operator.
First, we provide sufficient conditions that can be verified in practice to demonstrate that an operator $T$ is contracting on a small neighborhood of $\bar{x}$.
Although this first module focuses on finite-dimentional problems, we state these conditions in the more general setting of a [Banach space](https://en.wikipedia.org/wiki/Banach_space), as they will play an important role in the next modules.
"""

# ╔═╡ a801e2c7-c907-40d2-9976-bba19c590b2b
md"""
!!! theorem "Theorem (Contraction mapping)"
	Let $\mathcal{X}$ be a Banach space, $\bar{x} \in \mathcal{X}$, and $T : \mathcal{X} \to \mathcal{X}$ a continuously differentiable map.
	Assume that there exist a constant $Y$ and a non-decreasing map $Z : [0, \infty) \to [0, \infty)$ such that

	```math
	\begin{align}
	\| T(\bar{x}) - \bar{x} \| &\le Y, \\
	\| DT(x) \| &\le Z(\| x-\bar{x} \|), \qquad \forall x \in \mathcal{X}.
	\end{align}
	```

	If there exists $r > 0$ such that

	```math
	\begin{align}
	Y + \int_0^r Z(s) \, ds &\le r, \\
	Z(r) &< 1,
	\end{align}
	```

	then $T$ has a unique fixed point $\tilde{x}$ such that $\| \tilde{x} - \bar{x} \| \le r$.
"""

# ╔═╡ 6caaf267-3015-4a66-975f-ec1c3a6f86f8
md"""
!!! note "Remark"
	- The bound $Z$ is actually only needed locally, i.e. we only need $Z(s)$ for $s \le r$.
	  Therefore, one can fix a priori some $R > 0$, and weaken the assumption by only asking for $Z$ to be defined on $[0, R]$ and to satisfy $\| DT(x) \| \le Z(\| x-\bar{x} \|)$ for all $x \in \mathcal{X}$ such that $\| x-\bar{x} \| \le R$.
	  Of course, we are then only allowed to consider $r \in (0, R]$.
	- In practice, studying carefully how $\| DT(x) \|$ depends on $\| x-\bar{x} \|$ helps us to get better bounds.
	  However, we can sometimes get away with coarser estimates.
	  Specifically, we derive a constant $Z_*$ such that $\| DT(x) \| \le Z_*$ for all $x \in \mathcal{X}$ such that $\| x-\bar{x} \| \le R$, and then use $Z(s) = Z_*$ for all $s \in [0, R]$.
	  The conditions on $r$ then simplify, and, if $Z_* < 1$, we can take any $r \ge \frac{Y}{1-Z_*}$.
"""

# ╔═╡ 77380230-67f2-4794-bb51-39ceeac7bef1
md"""
For a dynamical system $x_{n+1} = f(x_n)$, a natural fixed-point operator for period-3 orbits is given by $T(x) \overset{\text{def}}{=} f^3(x)$, which represents $f$ composed with itself three times.
However, in general such a $T$ has no reason to be contracting.

To develop a more general procedure for constructing a fixed-point problem, one can start by finding a map $F$ whose isolated zeros are solutions to our problem.
For this example, we consider the mapping $F : \mathbb{R}^3 \to \mathbb{R}^3$ given by

```math
F(x_1, x_2, x_3) \overset{\text{def}}{=}
\begin{pmatrix}
f(x_1) - x_2 \\
f(x_2) - x_3 \\
f(x_3) - x_1
\end{pmatrix}.
```

Given an approximate solution $\bar{x}$, one can consider the Newton-like operator $T(x) \overset{\text{def}}{=} x - DF(\bar{x})^{-1} F(x)$.
This operator should be contracting in a neighborhood of $\bar{x}$, as $DT(\bar{x}) = 0$.

For finite-dimensional problems of moderate size, computing the inverse of $DF(\bar{x})$ is possible.
However, in infinite-dimensional problems, working with $DF(\bar{x})^{-1}$ can be challenging.
To address this, we often consider an approximate inverse $A \approx DF(\bar{x})^{-1}$ and define the fixed-point operator $T(x) \overset{\text{def}}{=} x - AF(x)$.
In the following results, we will operate under this framework, where $A$ is not specified and could be $DF(\bar{x})^{-1}$ or another suitable approximation.
"""

# ╔═╡ 08dab1b4-c78d-4e85-89b6-cf9b485030e0
md"""
!!! corollary "Corollary (Newton-Kantorovich)"
	Let $\mathcal{X}$ and $\mathcal{Y}$ be two Banach spaces, $\bar{x} \in \mathcal{X}$, $F : \mathcal{X} \to \mathcal{Y}$ a continuously differentiable map, $A : \mathcal{Y} \to \mathcal{X}$ am injective linear map, and $R > 0$.
	Assume that there exist constants $Y$, $Z_1$ and $Z_2$ such that

	```math
	\begin{align}
	\| AF(\bar{x}) \| &\le Y, \\
	\| I-ADF(\bar{x}) \| &\le Z_1, \\
	\| A(DF(x)-DF(\bar{x})) \| &\le Z_2 \| x-\bar{x} \|, \qquad \forall x \in \mathcal{X} \text{ such that } \| x-\bar{x} \| \le R.
	\end{align}
	```

	If there exists $r \in (0,R]$ such that

	```math
	\begin{align}
	Y + Z_1 r + \frac{1}{2} Z_2 r^2 &\le r, \\
	Z_1 + Z_2 r &< 1,
	\end{align}
	```

	then $F$ has a unique zero $\tilde{x}$ such that $\| \tilde{x} - \bar{x} \| \le r$.
"""

# ╔═╡ 3dc209e6-d7c0-464b-a477-0c46ede9483e
md"""
This corollary is a simplified version of the [Newton-Kantorovich Theorem](https://en.wikipedia.org/wiki/Kantorovich_theorem) [^Ort68].
Variations of this result can be found throughout the literature and are commonly used in computer-assisted proofs (CAPs) (see e.g. [^BL15] [^NPW19] and the references therein).
"""

# ╔═╡ f0d6ab5b-13cf-4408-aaae-93f959ce29c4
md"# Interval arithmetic"

# ╔═╡ 24ddc2c5-7976-4c67-8525-f6bb5d2e04a5
md"""
When we want to apply the contraction mapping theorem with $\bar{x}$ as an approximate solution obtained using the computer, we also need to compute values for the bounds $Y$ and $Z$.
However, these calculations typically rely on floating-point arithmetic, which introduces rounding errors.
For instance, if we execute `Y = norm(A * F(x̄), 1)` on the computer, we cannot be certain that $Y$ is an upper bound for $\|A F(\bar{x})\|_1$.

[Interval arithmetic](https://en.wikipedia.org/wiki/Interval_arithmetic) provides a robust method to control rounding errors and ensuring guaranteed results from computer calculations.
In this section, we offer a brief description of interval arithmetic; for more detailed discussions, please refer to [^Moo79] and [^Tuc11].
"""

# ╔═╡ 3b123092-1241-4a2a-bb61-36c7f68a0b35
md"## Overview"

# ╔═╡ 66ca51ad-0f2e-4f3c-8dd0-e66649b224d9
md"""
The central concept of interval arithmetic is to represent any real number by an interval, whose endpoints are floating-point numbers.
More precisely, for $a \in \mathbb{R}$, we define the interval $[\underline{a}, \overline{a}]$ where $\underline{a}, \overline{a}$ are floating-point numbers satisfying $\underline{a} \le a \le \overline{a}$.
While we give up the hope of representing numbers exactly, we recover guaranteed information: the real number $a$ is contained within the interval stored in the computer.

The rules of interval arithmetic ensure that this containment property is preserved for arithmetic operations.
For instance, consider two intervals $[a]$ and $[b]$ of the form $[\underline{a}, \overline{a}]$ and $[\underline{b}, \overline{b}]$, which enclose the reals numbers $a$ and $b$ respectively.
When we perform the operation $[c] = [a] + [b]$, the resulting interval $[c] = [\underline{c}, \overline{c}]$ is computed as follows:
- for $\underline{c}$, we take $\underline{a} + \underline{b}$ **rounded downward**, and
- for $\overline{c}$, we take $\overline{a}+\overline{b}$ **rounded upward**.
In particular, the real number $c = a + b$ is contained within the interval $[c]$.

Similarly, if $[c]$ is an interval $[\underline{c}, \overline{c}]$ and we compute $[d] = e^{[c]}$, the resulting interval $[d]$ contains $e^c$.

In this course, we use the Julia library [IntervalArithmetic](https://github.com/JuliaIntervals/IntervalArithmetic.jl), though many other interval arithmetic libraries are available in different programming languages.
Note that IntervalArithmetic is automatically included when using [RadiiPolynomial](https://github.com/OlivierHnt/RadiiPolynomial.jl).
This library provides standard arithmetic operations as well as implementations of elementary functions (such as `exp`, `log`, `cos`, etc.) that comply with the rules of interval arithmetic.
"""

# ╔═╡ 249bc3fc-3e7a-4dde-b526-06ddd77093aa
using IntervalArithmetic # not necessary if you did `using RadiiPolynomial`

# ╔═╡ f8eab26f-c893-4d2e-b4e8-6b59f33cbc9c
a = interval(2, 4)

# ╔═╡ 556163e0-80d2-4f4e-b198-b8e0922438db
cos(exp(sqrt(a)))

# ╔═╡ 97fadca1-fcf9-46b8-aa07-c615ca0deb7f
md"""
!!! note "Remark"
	Even when doing computer-assisted proofs, it’s not necessary to use interval arithmetic for every computation.
	For finding an approximate solution $\bar{x}$, floating-point calculations are perfectly fine.
	It’s only when computing the bounds for the contraction theorem, and verifying the radii polynomial inequalities, that interval arithmetic must be used.
"""

# ╔═╡ 0bf1cce8-6caf-40d5-93a6-b495b764d943
md"## Beware of typed-in floating-point numbers"

# ╔═╡ c09b05d4-0325-486a-a433-6b32fbb7dfc7
md"""
When you execute a command like `x = 0.1` in Julia, the variable `x` is assigned the closest floating-point number approximating $1/10$.
Since $1/10$ cannot be represented exactly in base 2, the value stored in `x` is not precisely equal to the typed-in number `0.1`.
In fact, it turns out to be slightly larger than $1/10$, but it could just as easily have been smaller (as seen with `x = 0.3`).
For more information on this topic, you can refer to this [webpage](https://0.30000000000000004.com).
"""

# ╔═╡ 2d68d26d-6e10-407f-8227-515ddadd9599
x = 0.1

# ╔═╡ 78f17262-7ffc-4d74-835c-f17863817b9a
x > 1//10

# ╔═╡ c44fd902-7655-43b5-82cc-53c2ecc4f77b
ix = interval(1)/interval(10)

# ╔═╡ 4c15bb79-2591-4f2a-9243-ff811de70df7
in_interval(1//10, ix)

# ╔═╡ 68243410-d2df-4606-82a5-53d515403e40
md"# Back to period-3 implies chaos"

# ╔═╡ 2698cfa7-5b3c-4f5a-9be4-b2c811e7a4ce
using RadiiPolynomial

# ╔═╡ 9a0e21bd-da27-48e6-9a08-b20b3cb805bd
md"""
Let’s recap what we have covered so far:
- Proving the existence of a period-3 orbit is sufficient to establish chaotic dynamics.
- We reformulated the task of finding a period-3 orbit as solving an equation $F = 0$.
- We derived sufficient conditions, in terms of bound estimates $Y, Z_1, Z_2$, to determine the radius of a ball, centered at a numerical approximation $\bar{x}$, that contains a (locally unique) true zero $\tilde{x}$ of $F$.
- We introduced interval arithmetic to rigorously estimate the bounds $Y$, $Z_1$, and $Z_2$ with computer assistance.

Now, we walk through the implementation of a computer-assisted proof to demonstrate a period-3 orbit for $\mu = 3.9$.

We start by defining the logistic map in Julia:
"""

# ╔═╡ b3121e23-03c4-46f2-89b7-c6b39f8efdad
f(x, μ) = μ * x * (1 - x) # logistic map

# ╔═╡ aeda376d-b1d1-4613-9964-05c24daf4b9a
md"""
and the zero-finding problem:
"""

# ╔═╡ c31f9913-0f16-4e74-b395-97038706be3c
function F(x, μ)
    v = zeros(eltype(x), 3)
    v[1] = f(x[1], μ) - x[2]
	v[2] = f(x[2], μ) - x[3]
	v[3] = f(x[3], μ) - x[1]
	return v
end

# ╔═╡ 129c8b50-1ba3-427e-a306-5be57c9755a3
md"""
To apply the Newton-Kantorovich theorem we also need the first derivate:
"""

# ╔═╡ d2f64f7e-2acb-4bf6-a961-3e2d8a87c409
Df(x, μ) = μ * (1 - 2x) # derivative of the logistic map w.r.t. x

# ╔═╡ 8591da5a-0b9a-4bff-8adc-fe1d3d9fca09
function DF(x, μ) # derivative of the zero-finding problem w.r.t. x
	M = zeros(eltype(x), 3, 3)
	M[1,1] = Df(x[1], μ)
	M[1,2] = -1
	M[2,2] = Df(x[2], μ)
	M[2,3] = -1
	M[3,1] = -1
	M[3,3] = Df(x[3], μ)
	return M
end

# ╔═╡ fa0be782-c09e-45b7-9cbb-d0b3161c3317
md"## Numerical zero"

# ╔═╡ cf53d067-ec00-4c06-b4ae-65ed7aa19d34
md"""
We apply [Newton's method](https://en.wikipedia.org/wiki/Newton%27s_method) to $F$ using the `newton` function provided in the RadiiPolynomial library.
"""

# ╔═╡ ef6afc18-f779-4370-8ccc-0a9494341166
initial_guess = [0.5, 1.0, 0.1]

# ╔═╡ a71e7ddf-e4d7-4200-abfb-b248883200e1
F_DF(x) = (F(x, μ), DF(x, μ)) # returns `F` and `DF` at `μ = 3.9`

# ╔═╡ b1ff2f19-31b7-4064-9ac0-e600ccd937d9
x_bar, success = newton(F_DF, initial_guess)

# ╔═╡ 8cc5d9d9-6630-43eb-85b7-2e658da6ec03
md"""
The `newton` function returns two outputs: the result of the iterations and a value indicating whether the method converged -- `true` if it did (as it does here) and `false` otherwise.
"""

# ╔═╡ 81d97321-ff79-42d6-adc6-3d66f91c438b
md"""
###### Initializing Newton's method

Finding a good starting point for Newton's method is not always easy.
Here we were a bit lucky that a simple guess worked well.
For other parameter values, however, we may need to take a closer look at the iterations of the map $f$.
"""

# ╔═╡ a2c9c052-1628-4ccd-8ac6-e2b16e0fc687
md"## Constructing $A$"

# ╔═╡ 32b5d274-3d0f-41aa-b1cf-85dfe40d0f76
md"""
As mentioned earlier, $A$ doesn’t need to be exactly $DF(\bar{x})^{-1}$; it only needs to be a reasonably close approximation.
For finite-dimensional problems, this inverse can be efficiently computed using built-in numerical algorithms.
"""

# ╔═╡ 16ec745d-0353-471f-b177-c25c5440e860
A = inv(DF(x_bar, μ))

# ╔═╡ 8ba3efdb-bff7-4211-bcd0-1537ae0420ed
md"## Verifying the contraction"

# ╔═╡ d9117644-f54f-42ce-94df-5933e57d765c
md"""
We now verify that $T(x) \overset{\text{def}}{=} x - A F(x)$ is a contraction around $\bar{x} \in \mathbb{R}^3$.
To rigorously compute the bounds $Y, Z_1, Z_2$, we enclose all relevant data into intervals.
"""

# ╔═╡ 787a3624-7a81-4f9f-9bca-f28e914a8f07
# since 39/10 is not representable as a floating-point number
iμ = interval(39)/interval(10)

# ╔═╡ a182a2c2-2e20-44f5-92ba-e85805050170
ix_bar = interval(x_bar)

# ╔═╡ f6759116-c48c-4a54-b073-8c2a0adc97cf
iA = interval(A)

# ╔═╡ 761a7427-6561-4754-8815-b1bf491a17da
md"##### Computing $Y$"

# ╔═╡ 7a037025-57ef-4b62-a880-8a2538096a9e
Y = norm(iA * F(ix_bar, iμ), 1)

# ╔═╡ 7a451708-9427-425a-93e0-0e66cad2c80b
md"##### Computing $Z_1$"

# ╔═╡ b7d36207-5f2c-41de-b142-7e8ee8700f79
Z₁ = opnorm(I - iA * DF(ix_bar, iμ), 1)

# ╔═╡ 7c1876f9-ccb7-4084-841f-316b1c1f6308
md"##### Computing $Z_2$"

# ╔═╡ 1a196987-95f2-48e0-9648-66373dce0ed5
md"""
```math
\|A(DF(x) - DF(\bar{x}))\|_1 \le \|A\|_1 \|DF(x) - DF(\bar{x})\|_1 = 2 \mu \|A\|_1 \| x - \bar{x} \|_1, \qquad \forall x \in \mathbb{R}^3.
```

So $Z_2 \ge 2 \mu \|A\|_1$ is a constant.
In other words, we have a Lipschitz control for $DF$ on the whole space $\mathbb{R}^3$, and we can freely choose $R = \infty$.
"""

# ╔═╡ 0b1f112c-d8b0-4903-966c-8ac6f4be69dc
R = Inf

# ╔═╡ 6f38f2e6-170c-4fea-b260-e55ca3fc3849
Z₂ = interval(2) * iμ * opnorm(iA, 1)

# ╔═╡ f7db0f90-473b-4265-8918-1b7acb08ebd2
md"##### Finishing the CAP"

# ╔═╡ afadd46e-458b-4975-925c-2add34ec75af
md"""
To determine the values of $r$ that satisfy the radii polynomial inequalities, we use the `interval_of_existence` function from RadiiPolynomial.
This function returns an interval of existence, within which all contained values satisfy the contraction conditions.
"""

# ╔═╡ 980ecce4-f02d-4976-b953-587b31009adb
ie = interval_of_existence(Y, Z₁, Z₂, R)

# ╔═╡ 958265d4-f604-47a4-8195-64ba21275f6d
md"""
The largest value in the interval represents the maximum radius of the ball centered at $\bar{x}$ within which the true zero $\tilde{x}$ of $F$ is proven to be unique.
On the other hand, the smallest value provides the sharpest error bound.
"""

# ╔═╡ 9eb669e1-f191-46cf-a3cb-b441b7d7f9e5
r = inf(ie)

# ╔═╡ c1a3ff98-5846-4601-87f7-b94b4f1fdf17
md"""
In fact, we cannot yet conclude that $\tilde{x}$ is a zero of $F$ since we did not check that $A$ was injective.
This property actually comes as a by-product of the contraction argument.
Indeed, $\|I - A DF(\bar{x})\| \le Z_1 < 1$ ensures that $A DF(\bar{x})$ is invertible, which implies that $A$ is surjective.
Since $A$ is a square matrix, it follows that $A$ must also be injective.

Finally, we verify that the true zero $\tilde{x}$ of $F$ is a period-3 orbit, that is $\tilde{x}_1 \ne \tilde{x}_2 \ne \tilde{x}_3$.
Since we have found $r$ such that $\| \tilde{x} - \bar{x} \|_1 \le r$, it follows that $\tilde{x}_j \in [\bar{x}_j - r, \bar{x}_j + r]$ for $j = 1, 2, 3$.
We then confirm that the three intervals $[\bar{x}_j - r, \bar{x}_j + r]$, for $j = 1, 2, 3$, are mutually disjoint:
"""

# ╔═╡ 611dc902-dd8d-4ec8-a2c0-1e9ac949a128
enclosure_x_tilde_1 = interval(ix_bar[1], r; format = :midpoint)

# ╔═╡ 61e71daf-8a6d-4c27-bccd-91f5e4db992c
enclosure_x_tilde_2 = interval(ix_bar[2], r; format = :midpoint)

# ╔═╡ 5f7e039e-b1a7-40a2-a076-ee56215da40c
enclosure_x_tilde_3 = interval(ix_bar[3], r; format = :midpoint)

# ╔═╡ 0f3d6fcf-2794-4d09-ae35-9dc5db3d0fa8
isdisjoint_interval(enclosure_x_tilde_1, enclosure_x_tilde_2, enclosure_x_tilde_3)

# ╔═╡ 1a9d9d5b-22bd-4fea-b480-817d80c63f8e
md"""
###### A final note on interval arithmetic

In the code provided in this notebook we used the Julia interval arithmetic library [IntervalArithmetic.jl](https://juliaintervals.github.io/IntervalArithmetic.jl/stable/). We took a bit of liberty by taking advantage of automatic conversion of floats and integers to intervals, when variables of these types are combined in elementary operations. The interval arithmetic library is not entirely happy about this:
"""

# ╔═╡ 26a665d6-77c1-4bb4-b2df-04e05a80ee80
isguaranteed(Y), isguaranteed(Z₁), isguaranteed(Z₂)

# ╔═╡ 0a2208b6-f90c-4587-978f-33643fd4bdba
md"# References"

# ╔═╡ b8d56ba6-01da-4604-8dad-3e63ec203fd4
md"""
[^BL15]: J. B. van den Berg and J.-P. Lessard. Rigorous numerics in dynamics. *Notices Of The American Mathematical Society*, 62(9), 2015.
[^LY75]: T.-Y. Li and J. A. Yorke. Period three implies chaos. *The American Mathematical Monthly*, 82(10):985--992, 1975.
[^Moo79]: R. E. Moore. *Methods and applications of interval analysis*. SIAM, 1979
[^Ort68]: J. M. Ortega. The Newton-Kantorovich theorem. *The American Mathematical Monthly*, 75(6):658--660, 1968.
[^NPW19]: M. T. Nakao, M. Plum, and Y. Watanabe. *Numerical Verification Methods and Computer-Assisted Proofs for Partial Differential Equations*. Springer Singapore, 2019.
[^Sha64]: A. N. Sharkowskii. Co-existence of the cycles of a continuous mapping of the line into itself. *Ukrainian Mathematical Journal*, 16(1), 1964.
[^Tuc11]: W. Tucker. *Validated numerics: a short introduction to rigorous computations*. Princeton University Press, 2011.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
IntervalArithmetic = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
RadiiPolynomial = "f2081a94-c849-46b6-8dc9-07bb90ed72a9"

[compat]
IntervalArithmetic = "~0.22.19"
Plots = "~1.40.9"
PlutoTeachingTools = "~0.3.1"
PlutoUI = "~0.7.60"
RadiiPolynomial = "~0.8.15"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.1"
manifest_format = "2.0"
project_hash = "06e17a9b7d31c4ac898dcef66dc31f92b0175889"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "8873e196c2eb87962a2048b3b8e08946535864a1"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+2"

[[deps.CRlibm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e329286945d0cfc04456972ea732551869af1cfc"
uuid = "4e9b3aee-d8a1-5a3d-ad8b-7d824db253f0"
version = "1.0.1+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "009060c9a6168704143100f36ab08f06c2af4642"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.2+1"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "7eee164f122511d3e4e1ebadb7956939ea7e1c77"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.3.6"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "bce6804e5e6044c6daab27bb533d1295e4a2e759"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.6"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "c785dfb1b3bfddd1da557e861b919819b82bbe5b"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.27.1"

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

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "64e15186f0aa277e174aa81798f7eb8598e0157e"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "ea32b83ca4fefa1768dc84e504cc0a94fb1ab8d1"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.4.2"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fc173b380865f70627d7dd1190dc2fce6cc105af"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.14.10+0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8e9441ee83492030ace98f9789a654a6d0b1f643"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+0"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "d36f682e590a83d63d1c7dbd287573764682d12a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.11"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "cc5231d52eb1771251fbd37171dbc408bcc8a1b6"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.4+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "53ebe7511fa11d33bec688a9178fac4e49eeee00"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.2"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

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

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "fa8e19f44de37e225aa0f1695bc223b05ed51fb4"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.3+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1ed150b39aebcc805c26b93a8d0122c940f64ce2"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.14+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "532f9126ad901533af1d4f5c198867227a7bb077"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.0+1"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "ee28ddcd5517d54e417182fec3886e7412d3926f"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.8"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "f31929b9e67066bee48eec8b03c0df47d31a74b3"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.8+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "b36c7e110080ae48fdef61b0c31e6b17ada23b33"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.82.2+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "01979f9b37367603e2848ea225918a3b3861b606"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+1"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "PrecompileTools", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "6c4d6a1babbbee6f283b3da64ac895f0a3bfbc96"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.11"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "401e4f3f30f43af2c8478fc008da50096ea5240f"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.3.1+0"

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
version = "1.11.0"

[[deps.IntervalArithmetic]]
deps = ["CRlibm_jll", "LinearAlgebra", "MacroTools", "RoundingEmulator"]
git-tree-sha1 = "24c095b1ec7ee58b936985d31d5df92f9b9cfebb"
uuid = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
version = "0.22.19"

    [deps.IntervalArithmetic.extensions]
    IntervalArithmeticDiffRulesExt = "DiffRules"
    IntervalArithmeticForwardDiffExt = "ForwardDiff"
    IntervalArithmeticIntervalSetsExt = "IntervalSets"
    IntervalArithmeticRecipesBaseExt = "RecipesBase"

    [deps.IntervalArithmetic.weakdeps]
    DiffRules = "b552c78f-8df3-52c6-915a-8e097449b14b"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "71b48d857e86bf7a1838c4736545699974ce79a2"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.9"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "be3dc50a92e5a386872a493a10050136d4703f9b"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.6.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "25ee0be4d43d0269027024d75a24c24d6c6e590c"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.4+0"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "10da5154188682e5c0726823c2b5125957ec3778"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.38"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "170b660facf5df5de098d866564877e119141cbd"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.2+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "36bdbc52f13a7d1dcb0f3cd694e01677a515655b"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.0.0+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "78211fb6cbc872f77cad3fc0b6cf647d923f4929"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "854a9c268c43b77b0a27f22d7fab8d33cdb3a731"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.2+1"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

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
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll"]
git-tree-sha1 = "8be878062e0ffa2c3f67bb58a595375eda5de80b"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.11.0+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c6ce1e19f3aec9b59186bdf06cdf3c4fc5f5f3e6"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.50.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "61dfdba58e585066d8bce214c5a51eaa0539f269"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+1"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "84eef7acd508ee5b3e956a2ae51b05024181dee0"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.40.2+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "b404131d06f7886402758c9ce2214b636eb4d54a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "edbf5309f9ddf1cab25afc344b1e8150b7c832f9"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.40.2+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

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
version = "1.11.0"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f02b56007b064fbfddb4c9cd60161b6dd0f40df3"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.1.0"

[[deps.LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "688d6d9e098109051ae33d126fcfc88c4ce4a021"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "3.1.0"

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
version = "1.11.0"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

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
git-tree-sha1 = "7493f61f55a6cce7325f197443aa80d32554ba10"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.15+1"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6703a85cb3781bd5909d48730a67205f3f31a575"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.3+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e127b609fb9ecba6f201ba7ab753d5a605d53801"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.54.1+0"

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
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "41031ef3a1be6f5bbbf3e8073f210556daeae5ca"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.3.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "3ca9a356cd2e113c420f2c13bea19f8d3fb1cb18"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.3"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "dae01f8c2e069a683d3a6e17bbae5070ab94786f"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.9"

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
deps = ["Downloads", "HypertextLiteral", "Latexify", "Markdown", "PlutoLinks", "PlutoUI"]
git-tree-sha1 = "8252b5de1f81dc103eb0293523ddf917695adea1"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.3.1"

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
version = "1.11.0"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "492601870742dcd38f233b23c3ec629628c1d724"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.7.1+1"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll"]
git-tree-sha1 = "e5dd466bf2569fe08c91a2cc29c1003f4797ac3b"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.7.1+2"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "1a180aeced866700d4bebc3120ea1451201f16bc"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.7.1+1"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "729927532d48cf79f49070341e1d918a65aba6b0"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.7.1+1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.RadiiPolynomial]]
deps = ["IntervalArithmetic", "LinearAlgebra", "Printf", "Reexport", "SparseArrays"]
git-tree-sha1 = "891c7a0c4355fdef85cfe09abcd8c7b0d896b0e6"
uuid = "f2081a94-c849-46b6-8dc9-07bb90ed72a9"
version = "0.8.15"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

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
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "470f48c9c4ea2170fd4d0f8eb5118327aada22f5"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.6.4"

[[deps.RoundingEmulator]]
git-tree-sha1 = "40b9edad2e5287e05bd413a38f61a8ff55b9557b"
uuid = "5eaf0fd0-dfba-4ccb-bf02-d820a40db705"
version = "0.2.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "83e6cce8324d49dfaf9ef059227f91ed4441a8e5"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.2"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

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

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

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
version = "1.11.0"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

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
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "d95fe458f26209c66a187b1114df96fd70839efd"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.21.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "975c354fcd5f7e1ddcc1f1a23e6e091d99e99bc8"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.4"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
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

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "a2fccc6559132927d4c5dc183e3e01048c6dcbd6"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.5+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "a54ee957f4c86b526460a720dbc882fa5edcbefc"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.41+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "15e637a697345f6743674f1322beefbc5dcd5cfc"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.6.3+0"

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
git-tree-sha1 = "2b0e27d52ec9d8d483e2ca0b72b3cb1a8df5c27a"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+1"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "02054ee01980c90297412e4c809c8694d7323af3"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+1"

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
git-tree-sha1 = "fee57a273563e273f0f53275101cd41a8153517a"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+1"

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
git-tree-sha1 = "b9ead2d2bdb27330545eb14234a2e300da61232e"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+1"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "555d1076590a6cc2fdee2ef1469451f872d8b41b"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.6+1"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "gperf_jll"]
git-tree-sha1 = "431b678a28ebb559d224c0b6b6d01afce87c51ba"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.9+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6e50f145003024df4f5cb96c7fce79466741d601"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.56.3+0"

[[deps.gperf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3516a5630f741c9eecb3720b1ec9d8edc3ecc033"
uuid = "1a1c6b14-54f6-533d-8383-74cd7377aa70"
version = "3.1.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1827acba325fdcdf1d2647fc8d5301dd9ba43a9d"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.9.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e17c115d55c5fbb7e52ebedb427a0dca79d4484e"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.2+0"

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
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a22cf860a7d27e4f3498a0fe0811a7957badb38"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.3+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "ad50e5b90f222cfe78aa3d5183a20a12de1322ce"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.18.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "b70c870239dc3d7bc094eb2d6be9b73d27bef280"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.44+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "490376214c4721cdaca654041f635213c6165cb3"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+2"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "814e154bdb7be91d78b6802843f76b6ece642f11"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.6+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

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
# ╟─a5169b4b-f0d4-4f66-8198-0aec5f8e135a
# ╠═d0e623ee-b096-4a27-977d-dc32567d6020
# ╠═2661bfc9-e398-41ed-87d9-c78f05da64cb
# ╟─b3845641-1537-4a27-8550-1eff30900a6b
# ╟─ed6f44f4-3696-494a-b8e3-b30184f5bb06
# ╟─c33dc650-3f94-11ef-398a-8bbc4a2b69b8
# ╟─9e0bfed7-7ccc-42b2-9b1a-15420d85c356
# ╠═e3ca081f-0e6a-4dcf-9436-ca4f34467379
# ╟─79138e2e-6123-4923-8329-e5c17499d785
# ╟─1d9286fa-647a-4743-806f-9cea5aab6be1
# ╟─1e56532e-9886-4a41-9ba0-1a619a85e014
# ╟─8d6c631f-832a-48a9-891f-650f53a0ac6c
# ╟─429e5261-d7d0-47a3-94b5-8e7905a5bd53
# ╟─a27b45ca-77d1-408e-98e9-6b5248aca85e
# ╟─a801e2c7-c907-40d2-9976-bba19c590b2b
# ╟─6caaf267-3015-4a66-975f-ec1c3a6f86f8
# ╟─77380230-67f2-4794-bb51-39ceeac7bef1
# ╟─08dab1b4-c78d-4e85-89b6-cf9b485030e0
# ╟─3dc209e6-d7c0-464b-a477-0c46ede9483e
# ╟─f0d6ab5b-13cf-4408-aaae-93f959ce29c4
# ╟─24ddc2c5-7976-4c67-8525-f6bb5d2e04a5
# ╟─3b123092-1241-4a2a-bb61-36c7f68a0b35
# ╟─66ca51ad-0f2e-4f3c-8dd0-e66649b224d9
# ╠═249bc3fc-3e7a-4dde-b526-06ddd77093aa
# ╠═f8eab26f-c893-4d2e-b4e8-6b59f33cbc9c
# ╠═556163e0-80d2-4f4e-b198-b8e0922438db
# ╟─97fadca1-fcf9-46b8-aa07-c615ca0deb7f
# ╟─0bf1cce8-6caf-40d5-93a6-b495b764d943
# ╟─c09b05d4-0325-486a-a433-6b32fbb7dfc7
# ╠═2d68d26d-6e10-407f-8227-515ddadd9599
# ╠═78f17262-7ffc-4d74-835c-f17863817b9a
# ╠═c44fd902-7655-43b5-82cc-53c2ecc4f77b
# ╠═4c15bb79-2591-4f2a-9243-ff811de70df7
# ╟─68243410-d2df-4606-82a5-53d515403e40
# ╠═2698cfa7-5b3c-4f5a-9be4-b2c811e7a4ce
# ╟─9a0e21bd-da27-48e6-9a08-b20b3cb805bd
# ╠═b3121e23-03c4-46f2-89b7-c6b39f8efdad
# ╟─aeda376d-b1d1-4613-9964-05c24daf4b9a
# ╠═c31f9913-0f16-4e74-b395-97038706be3c
# ╟─129c8b50-1ba3-427e-a306-5be57c9755a3
# ╠═d2f64f7e-2acb-4bf6-a961-3e2d8a87c409
# ╠═8591da5a-0b9a-4bff-8adc-fe1d3d9fca09
# ╟─fa0be782-c09e-45b7-9cbb-d0b3161c3317
# ╟─cf53d067-ec00-4c06-b4ae-65ed7aa19d34
# ╠═ef6afc18-f779-4370-8ccc-0a9494341166
# ╠═a71e7ddf-e4d7-4200-abfb-b248883200e1
# ╠═b1ff2f19-31b7-4064-9ac0-e600ccd937d9
# ╟─8cc5d9d9-6630-43eb-85b7-2e658da6ec03
# ╟─81d97321-ff79-42d6-adc6-3d66f91c438b
# ╟─a2c9c052-1628-4ccd-8ac6-e2b16e0fc687
# ╟─32b5d274-3d0f-41aa-b1cf-85dfe40d0f76
# ╠═16ec745d-0353-471f-b177-c25c5440e860
# ╟─8ba3efdb-bff7-4211-bcd0-1537ae0420ed
# ╟─d9117644-f54f-42ce-94df-5933e57d765c
# ╠═787a3624-7a81-4f9f-9bca-f28e914a8f07
# ╠═a182a2c2-2e20-44f5-92ba-e85805050170
# ╠═f6759116-c48c-4a54-b073-8c2a0adc97cf
# ╟─761a7427-6561-4754-8815-b1bf491a17da
# ╠═7a037025-57ef-4b62-a880-8a2538096a9e
# ╟─7a451708-9427-425a-93e0-0e66cad2c80b
# ╠═b7d36207-5f2c-41de-b142-7e8ee8700f79
# ╟─7c1876f9-ccb7-4084-841f-316b1c1f6308
# ╟─1a196987-95f2-48e0-9648-66373dce0ed5
# ╠═0b1f112c-d8b0-4903-966c-8ac6f4be69dc
# ╠═6f38f2e6-170c-4fea-b260-e55ca3fc3849
# ╟─f7db0f90-473b-4265-8918-1b7acb08ebd2
# ╟─afadd46e-458b-4975-925c-2add34ec75af
# ╠═980ecce4-f02d-4976-b953-587b31009adb
# ╟─958265d4-f604-47a4-8195-64ba21275f6d
# ╠═9eb669e1-f191-46cf-a3cb-b441b7d7f9e5
# ╟─c1a3ff98-5846-4601-87f7-b94b4f1fdf17
# ╠═611dc902-dd8d-4ec8-a2c0-1e9ac949a128
# ╠═61e71daf-8a6d-4c27-bccd-91f5e4db992c
# ╠═5f7e039e-b1a7-40a2-a076-ee56215da40c
# ╠═0f3d6fcf-2794-4d09-ae35-9dc5db3d0fa8
# ╟─1a9d9d5b-22bd-4fea-b480-817d80c63f8e
# ╠═26a665d6-77c1-4bb4-b2df-04e05a80ee80
# ╟─0a2208b6-f90c-4587-978f-33643fd4bdba
# ╟─b8d56ba6-01da-4604-8dad-3e63ec203fd4
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

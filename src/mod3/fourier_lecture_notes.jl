### A Pluto.jl notebook ###
# v0.19.45

#> [frontmatter]
#> chapter = 3
#> order = 1
#> title = "Lecture notes: periodic orbit problems"
#> tags = ["module3"]
#> layout = "layout.jlhtml"

using Markdown
using InteractiveUtils

# ╔═╡ d0e623ee-b096-4a27-977d-dc32567d6020
using PlutoTeachingTools, PlutoUI

# ╔═╡ 2661bfc9-e398-41ed-87d9-c78f05da64cb
using RadiiPolynomial, Plots

# ╔═╡ 018ecc45-8638-4a59-b561-efb086bdc751
using LinearAlgebra,ToeplitzMatrices

# ╔═╡ 7fc40507-eda3-474d-a454-04e9173a7adb
html"""<style>
main {
    max-width: 1000px;
    margin-left: auto;
    margin-right: auto;
    text-align: justify;
}
"""

# ╔═╡ d5a510a3-c518-47ed-96bc-7bb22e3b08b5
md"""
We continue to do CAPs to find solutions of ODEs, but rather than looking at an initial value problem, we now look for *periodic* solutions. The setting is again infinite dimensional, with some small differences:
- we will use **Fourier series** rather than Taylor series
- the Fourier coefficients are two-sided sequences of complex numbers
- the convolution product is given by an infinite series
- the problem of finding a closed orbit is global in the sense that there is no equivalent of "short time" integration
"""

# ╔═╡ f73b88cb-19e2-4d50-a45e-fbe90fb691cd
md"""
# Simple example: a forced, damped, nonlinear pendulum
"""

# ╔═╡ cf717e9b-26d8-449d-9f3b-d1d5d0bd2bb5
md"""
## Introduction and setup
"""

# ╔═╡ c33dc650-3f94-11ef-398a-8bbc4a2b69b8
md"""
We consider the second order ODE

$\begin{align}
u''(t) + \beta_1 u'(t) + \beta_2 u(t) - u(t)^2 = \beta_3 \cos(t)
\end{align}$

where $\beta_1, \beta_2, \beta_3$ are parameters, and we aim to find a periodic solution $u(t)$.
"""

# ╔═╡ 58892dec-c89c-4f1f-9f97-700bc27ef338
md"""
Since the (minimal) period of the forcing is $2\pi$, we will look for a periodic solution of period $2\pi$.
"""

# ╔═╡ 1bfcc76a-266a-4784-b76b-3871445e603e
Foldable("""Why is this a "simple" example?""", md"""
A more natural problem to look at is to find a periodic orbit in a system of ODEs. In such an *autonomous* system of ODEs the period of a periodic solution is unknown a priori, and determining the period $\tau$ is thus part of the problem, in addition to finding the unknown $\tau$-periodic function $u(t)$. Furthermore, when $u(t)$ is a periodic solution, so is $u(t-t_0)$ for any time shift $t_0$. In particular, periodic solutions are not isolated unless an additional phase condition is imposed. An additional complication is that for a system of equations the bookkeeping is more involved than for a single equation.

Although these difficulties can be overcome, the notation can get in the way of understanding. Since we want to focus here on how to deal with periodic functions, we choose to look at a scalar, nonautonomous problem where the period is given.
""")

# ╔═╡ 3e3a3a60-f76f-401d-b265-a99e713d98c5
md"""
### Fourier series
"""

# ╔═╡ 8f3e96dc-0cdc-411e-8f94-e9e47488abf3
md"""
We substitute the Fourier series 

$\begin{align}
u(t) = \sum_{n\in\mathbb{Z}} a_n e^{int}
\end{align}$

in the ODE to find the infinite set of equations ($n \in \mathbb{Z}$)

$\begin{align}
\lambda_n a_n + (a*a)_n + c_n = 0.
\end{align}$

where
$\lambda_n:= n^2-in \beta_1-\beta_2$, while
$c$ represents the Fourier transform of $\beta_3 \cos(t)$,
 given by

$\begin{align}
c_n := \begin{cases} \frac{\beta_3}{2} & n= \pm 1\\ 0 & \text{otherwise} \end{cases}
\end{align}$
"""

# ╔═╡ 7c9ca144-02f1-4941-bb0b-becb96b4fc8a
md"""
### Convolution
"""

# ╔═╡ b5d97ade-bcde-4c7a-ae78-697e73c10a0c
md"""
The convolution product is given by 

$\begin{align}
(a*b)_n := \sum_{k\in\mathbb{Z}} a_k b_{n-k} 
\end{align}$
"""

# ╔═╡ 47e39d6a-231e-40d1-906e-eeafe1b1f807
Markdown.MD(Markdown.Admonition("tip", "Lemma",[md"""When $u(t) = \sum_{n\in\mathbb{Z}} a_n e^{int}$
and $\tilde{u}(t) = \sum_{n\in\mathbb{Z}} \tilde{a}_n e^{int}$ then for their product we have $u(t)\tilde{u}(t) = \sum_{n\in\mathbb{Z}} (a*\tilde{a})_n e^{int}$.
 """]))

# ╔═╡ ac3e3cf5-e5c8-4863-8c7b-ccd5f58a3073
md"""
When $a_n=0$ for $|n|>N_a$ and $b_n=0$ for $|n|>N_b$ then the convolution can be computed. 
An implementation is below.
"""

# ╔═╡ 4a8b8e9e-9562-4cb5-a029-856127e84b98
function convolution(a, b, nab=-1)
    # convolution for 1D Fourier; nab is size of output (default is all nonzero coefficients)
    na = (length(a)-1)÷2
    nb = (length(b)-1)÷2
	if nab<0
		nab = na + nb
	end
	
    # Initialize the result array with zeros
    result = zeros(typeof(a[1]),2*nab+1)
    
    # Perform the convolution
    for ia in -na:na
		minib=max(-nab-ia,-nb)
		maxib=min(nab-ia,nb)
        for ib in minib:maxib 
            result[nab+ia+ib+1] += a[na+ia+1] * b[nb+ib+1]
        end
    end
    
    return result
end

# ╔═╡ 312d9592-1be0-4626-8773-51cc8604dd6a
md"""
### Norm
"""

# ╔═╡ c4d9c4f0-fd89-4a2a-96a7-17c4ba6dd747
md"""
We will work in the sequence space 

$\begin{align}
X=X_\nu=\{ a \in \mathbb{C}^\mathbb{Z} : \|a\|_X:= \sum_{n\in\mathbb{Z}} \nu^{|n|} |a_n| < \infty \}.
\end{align}$

We will choose the value of $\nu$ later.
"""

# ╔═╡ 766164f3-eef3-4798-bf7c-afdca52d3ba0
Markdown.MD(Markdown.Admonition("tip", "Lemma",[md"""When $a\in X_\nu$ for some $\nu \geq 1$ then the Fourier series $\sum_{n\in\mathbb{Z}} a_n e^{int}$ converges uniformly. When $\nu>1$ then term-by-term differentiation up to any order of the Fourier series is justified."""]))

# ╔═╡ bbf2e943-8c5c-469f-a2f3-5f8b758f6638
Foldable("""Is X with this norm a Banach space?""", md"""**Exercise**: Yes""")

# ╔═╡ f91abfa7-8a8f-4204-95d9-637c622dd720
md"""
The product and the norm play together nicely: they give the space $X$ the structure of a Banach algebra, as expressed by the following lemma.
"""

# ╔═╡ 8f6431ed-346c-4c51-83c4-ea7ec12a8d60
Markdown.MD(Markdown.Admonition("tip", "Lemma (Banach algebra property)",[md"""For any $\nu \geq 1$ and $a,b\in X=X_\nu$ we have $\|a*b\|_X\leq \|a\|_X \|b\|_X$."""]))

# ╔═╡ c5cfe80e-1728-4393-b837-187e3e7d49a9
Foldable("""Proof""", md"""**Exercise** (rearranging series and using triangle inequality)""")

# ╔═╡ 39d13194-bca9-47c1-a686-2ea5f9d62289
md"""
The operator norm (recall: $\|\Gamma\|_{B(X)}$ is the smallest number such that $\|\Gamma a\|_X\leq \|\Gamma\|_{B(X)} \|a\|_X$ for all $a\in X$) has a particularly nice representation for weighted $l^1$ spaces such as $X_\nu$.
"""

# ╔═╡ 7e78028c-f70b-4474-8bee-783ad7d99d56
Markdown.MD(Markdown.Admonition("tip", "Lemma (characterisation of operator norm)",[md"""
Let $\Gamma$ be a bounded operator on $X$, represented by $(\Gamma a)_n=\sum_{k\in\mathbb{Z}} \Gamma_{nk} a_k$. 
Then 

$\begin{align}
  \|\Gamma\|_{B(X)} = \sup_{k\in\mathbb{Z}} \frac{1}{\nu^{|k|}} \sum_{n\in\mathbb{Z}} \nu^{|n|} |\Gamma_{nk}| .
\end{align}$

Let $e_k$ denote the basis vector given by $(e_k)_n=\delta_{kn}$ for any $k,n \in \mathbb{Z}$. Then 

$\begin{align}
  \|\Gamma\|_{B(X)} = \sup_{k\in\mathbb{Z}} \frac{\|\Gamma e_k\|_X}{\|e_k\|_X}.
\end{align}$
"""]))

# ╔═╡ c7a10c87-a95d-4b0b-b63d-f9546e97bb36
md"""
Code for the (operator) norm is given below for the case of finitely represented elements/operators in/on $X$.
"""

# ╔═╡ efb766e3-e41f-4113-ae6e-c4dd7766c784
function nunorm(B, nu)
    # weighted l1-norm; works for both vectors and matrices B
    n1 = (size(B,1)-1)÷2 
    n2 = (size(B,2)-1)÷2

	weights=[nu^(abs(i)) for i in -n1:n1]'
	invweights=[nu^(-abs(i)) for i in -n2:n2]'

	result=maximum((weights*abs.(B)).*invweights)
    
	return result
end

# ╔═╡ 0f04a415-b4c1-4b05-a35a-6d3029decb12
maximum([1.9414153246843213 1.9290258804863343 1.4084075531504734])

# ╔═╡ 31a59e3d-32af-4328-b3bc-cc434728429d
md"""
## Zero finding problem
"""

# ╔═╡ 8d7fb15d-32d3-4fbb-9d06-c8ec1ed7d00e
md"""
We define the zero finding problem $F(a)=0$ on $X$ by ($n \in \mathbb{Z}$)

$\begin{align}
F_n(a) := a_n + \lambda_n^{-1} [(a*a)+c]_n .
\end{align}$
"""

# ╔═╡ b2defc79-da78-47c8-bb78-44cb9000ff58
Foldable("""Why divide by lambda?""", md"""
Dividing by $\lambda_n$ is a choice that simplifies the linear term, which is beneficial for some of the algebra. In terms of estimates it merely shifts the difficulty to the nonlinear term. Working with the alternative $\tilde{F}_n(a) := \lambda_n a_n + [(a*a)+b]_n$ is also a perfectly valid choice.
""")

# ╔═╡ 9ab6ce17-0637-4008-be29-e0bf4fda4287
Markdown.MD(Markdown.Admonition("note", "Conjugate symmetry",[md"""
We denote the *complex conjugate* of $z\in \mathbb{C}$ by $\overline{z}$. For real-valued functions $u(t)$ the Fourier coefficients satisfy $a_{-n} = \overline{a_n}$ for all $n\in \mathbb{Z}$. And the other way around: if $a_{-n} =\overline{a_n}$ for all $n\in \mathbb{Z}$, then the Fourier series is real-valued.

For $a\in X$, we denote by $a^\dagger \in X$ the *conjugate* given by $(a^\dagger)_n = \overline{a_{-n}}$. With this notation, $a=a^\dagger$ implies the Fourier series is real-valued.

We note that both $\lambda_{-n} =\overline{\lambda_n}$ and $c_{-n} =\overline{c_n}$. Additionally, $a^\dagger * b^\dagger = (a * b)^\dagger$. Hence the zero finding problem is conjugate symmetric:

$\begin{align} F(a^\dagger)=F(a)^\dagger . \end{align}$
"""]))

# ╔═╡ 9002b709-ec9c-4063-a3a8-7b0c00834095
Foldable("""Exercise""", md"""
Prove that $a^\dagger * b^\dagger = (a * b)^\dagger$ and $F(a^\dagger)=F(a)^\dagger$.
""")

# ╔═╡ c2f4618c-dd22-477c-8fe7-324c82d3ff46
md"""
We now formalize that the zero finding problem leads to periodic solutions of our differential equation.
"""

# ╔═╡ 7d79ed91-7cbb-48b4-966a-d3165f2a3d89
Markdown.MD(Markdown.Admonition("tip", "Lemma",[md"""Let $\nu>1$. When $a\in X$ with $a^\dagger=a$ solves $F(a)=0$, then the corresponding Fourier series $u(t)=\sum_{n\in\mathbb{Z}} a_n e^{int}$ is a real-valued periodic solution of the differential equation $u''(t) + \beta_1 u'(t) + \beta_2 u(t) - u(t)^2 = \beta_3 \cos(t)$."""]))

# ╔═╡ 66deaecf-e201-40ec-88d3-99574195253a
md"""
## Finite dimensional projection
"""

# ╔═╡ c32d9839-7cb7-4018-bab1-b2211f05b71e
md"""
We need the projection operator on a finite number of modes ($2N+1$ coefficents) and its complement:

$\begin{align}
(\pi^{\leq N} a)_n := \begin{cases}
a_n & |n| \leq N \\
0 & |n| > N
\end{cases}
& \qquad\qquad
(\pi^{>N} a)_n := \begin{cases}
0 & |n| \leq N \\
a_n & |n| > N
\end{cases}
\end{align}$

so that $a=\pi^{\leq N}a+\pi^{>N}a$ and $\pi^{\leq N}\pi^{>N}=0=\pi^{>N}\pi^{\leq N}$.
"""

# ╔═╡ bd76555d-aa30-4313-b4c4-c3f4b9000194
Foldable("""Why is this so complicated?""", md"""
The strange looking superscripts $\leq N$ and $>N$ are in fact meant to indicate clearly which index values are involved.

The projections provide a sleek way to do the algebra that essentially splits the problem into a finite dimensional part handled by the computer and an infinite dimensional "tail" part handled entirely by human analysis. It avoids writing down *infinite* matrices or formulas with a lot of indices.
""")

# ╔═╡ 740167dc-5bcd-48c0-b1f0-a7bf5e1b6cb7
md"""
The range $\pi^{\leq N}X$ is finite dimensional, and 
the restriction of $\pi^{\leq N} F$ to $\pi^{\leq N} X$, which we denote by $F^{\leq N}$, is what we will work with in the computer. An implementation is given below"""

# ╔═╡ 026aec24-4ed0-4189-af4f-28f56e6964ef
function F(a, beta, Nf=-1)
	# default length of output is same as input
    Na = (length(a) - 1) ÷ 2
	if Nf<0
		Nf = Na
	end

	# deal with sizes
	Naf=min(Na,Nf);
	a1=zeros(typeof(a[1]),2*Nf+1);
	a1[Nf+1-Naf:Nf+1+Naf]=a[Na+1-Naf:Na+1+Naf]

	# the nonlinear term
	nonlinear=convolution(a,a,Nf)
	# add the forcing term
	if Nf>0
	    nonlinear[Nf] += beta[3]/2
		nonlinear[Nf+2] += beta[3]/2
	end

	# linear terms
    lambdainv = 1 ./ [n^2 - 1im * beta[1] * n - beta[2] for n in -Nf:Nf]
    result = a1 + lambdainv.*nonlinear

    return result
end

# ╔═╡ 06f3d467-af63-440d-b486-8da7f67c314f
md"""It is convenient to introduce the diagonal linear operator 

$\begin{align}
(\Lambda a)_n := \lambda_n a_n
\end{align}$

When $\beta_1 \neq 0$ then the inverse $\Lambda^{-1}$ is well-defined as a bounded operator on $X$. 
"""

# ╔═╡ f5cc03ca-09eb-4419-9d51-78c956ee73d9
md"""
The derivative of $F(a)=a+\Lambda^{-1}[a*a+c]$ can now be written compactly as 

$\begin{align}
DF(a)b = b + 2 \Lambda^{-1} [a*b].
\end{align}$
"""

# ╔═╡ c5706263-a08a-4742-b376-b980232d6191
md"""
An implementation of the Jacobian $DF^{\leq N}$ of the truncated problem is given below.
"""

# ╔═╡ e5b2628d-2e64-47ca-b594-a0e4cb8ea64a
function DF(a, beta, Ndf=-1)
	# Jacobian; default size of output is same as input
    Na = (length(a) - 1) ÷ 2
	if Ndf<0
		Ndf = Na
	end

	# deal with nonlinear term in terms of a Toeplitz matrix
	Naf=min(Na,2*Ndf);
	a1c=zeros(typeof(a[1]),2*Ndf+1);
	a1r=zeros(typeof(a[1]),2*Ndf+1);
	a1c[1:Naf+1]=a[Na+1:Na+1+Naf]
	a1r[1:Naf+1]=a[Na+1:-1:Na+1-Naf]
    Ta=Toeplitz(a1c,a1r);

	# linear terms
    lambda = [n^2 - 1im * beta[1] * n - beta[2] for n in -Ndf:Ndf]
 	Lambdainv=Diagonal( 1 ./ lambda );
    result = I(2*Ndf+1)+2*Lambdainv*Ta;    
	
    return result
end

# ╔═╡ c4363fdf-6500-4830-bc42-9777d616ac34
md"""
## Choosing the approximate inverse of the Jacobian
"""

# ╔═╡ 24e1a5af-569e-4ae8-a3d7-19b7962db48b
md"""
We also use the project to split $A$ into a finite part and a tail part:

$\begin{align}A:=A^{\leq N}\pi^{\leq N} + I^{>N}\pi^{>N}, \end{align}$

where $A^{\leq N}$ is a linear map (to be chosen below) on $\pi^{\leq N}X$ and $I^{>N}$ is the identity on
$\pi^{> N}X$. Hence by construction one may also write this is block diagonal form:

$\begin{align}A=\pi^{\leq N}A^{\leq N}\pi^{\leq N} + \pi^{> N}I^{>N}\pi^{>N}, \end{align}$
"""

# ╔═╡ 77b91c27-ddb0-4115-b0d0-a78f2729dadc
Foldable("""Is A injective?""", md"""
The operator $A$ is injective on $X$ if and only $A^{\leq N}$ is injective on $\pi^{\leq N} X$. The operator $A^{\leq N}$ has a (finite) matrix representation, hence injectivity can be proved by computer.

One may also observe that $A=I+\pi^{\leq N}(A^{\leq N}-I^{\leq N})\pi^{\leq N}$, hence $A$ is of the form identity plus compacy (and thus a Fredholm operator of index $0$).
""")

# ╔═╡ 51e35c7e-77f2-4911-914c-7f1a242c3b91
md"""
## Numerics
"""

# ╔═╡ 4a9daa68-86ed-4327-98fc-03d436885d2e
md"""
We consider the finite truncation $F^{\leq N}$ and use a Newton algorithm to find a point $\bar{a} \in \pi^{\leq N} X$ such that $\pi^{\leq N} F (\bar{a}) \approx 0$. If the truncation dimensions $N$ is large enough we may hope that also 
$\pi^{> N} F (\bar{a}) \approx 0$, since we expect the coefficients of $\bar{a}$ to decrease for large $n$.
"""

# ╔═╡ 3842ad29-ae9c-4743-b487-e76d94fb01cf
md"""
We choose parameter values
"""

# ╔═╡ 1688b3d8-8228-45ee-8729-22265a86484a
# in interval form
beta= [ I"0.1"
	    I"4.0"
	    I"1.0" ]

# ╔═╡ aa150fa3-1729-43b1-bb55-73c741221760
md"""
and truncation dimension
"""

# ╔═╡ e1f6ef9b-5818-415a-8945-97b65288a433
N=10

# ╔═╡ d0858eee-5086-4754-966f-de2940183573
md"""
To simplify some of the expressions appearing in the bounds, we will for convenience assume that $N \geq \sqrt{\beta_3}$.
"""

# ╔═╡ a89606ee-b46f-4656-9cb9-d89850646593
println(N >= sup(sqrt(beta[3])))

# ╔═╡ 192507bd-8285-4581-8a6d-f9ff1b2c4793
md"""
We then run Newton iterations.
"""

# ╔═╡ 301bf3ed-b180-4c44-a92f-0f9f4c2509ed
begin
	initialdata = zeros(typeof(0.0im),2*N+1)   # Fourier coefficients are complex-valued
	initialdata[N+1] += 0.0 	
	initialdata=(initialdata+conj(reverse(initialdata)))/2

	Fbeta=mid.(beta)   # floats
	
	a0, success = newton(a -> (F(a, Fbeta), DF(a, Fbeta)), initialdata)
	a0=(a0+conj(reverse(a0)))/2    # symmetry
end

# ╔═╡ aa320281-42ce-4035-8c55-1f44661737e2
md"""
We plot the approximate periodic solution
"""

# ╔═╡ d5800e45-8809-46d6-9f41-7565e5b1cbfd
begin
	Na= (length(a0)-1) ÷ 2
	Nplot=201;
	t=LinRange(0, 2*pi, Nplot)
	u=zeros(typeof(0.0im),Nplot)
	for n=-Na:Na
		u += a0[Na+1+n]*exp.(im*n*t)
	end
	u=real(u)
	plot(t, u, legend=false)
	xlabel!("t")
	ylabel!("u")
	xlims!(0, 2*pi)
end

# ╔═╡ 5f3d91ee-14c2-4b5e-8e2a-a607d086eb51
Foldable("""Initializing Newton's method""", md"""Finding a good starting point for Newton's method is not always easy. Here we were a bit lucky that a simple guess suffices. For other parameter values one may need to study the solution of the ODE in more detail using numerical integration.""")

# ╔═╡ 936423fd-0d78-41da-8b49-5fca3c125547
md"""
Next we compute a numerical inverse of the Jacobian $DF^{\leq N}(\bar{a})$ and we denote this inverse matrix by $A^{\leq N}$.
"""

# ╔═╡ e37811fa-cd67-424e-8707-25ad8534d293
AN=inv(DF(a0,Fbeta))

# ╔═╡ 0b4e2571-cbec-4a4d-aac6-dd8b32e3f474
Markdown.MD(Markdown.Admonition("note", "Remark on truncation dimension",
[md"""
It is not necessary to choose the finite truncation dimensions for $\bar{a}\in \pi^{\leq N} X$ and $A^{\leq N}$ equal. 
For purposes of exposition we do not introduce two different truncation dimension parameters here, but essentially the truncation dimension for $\bar{a}\in \pi^{\leq N} X$ controls the size of the residue, while the truncation dimension for the approximate inverse $A^{\leq N}$ controls the contractivity of the fixed point operator. These can in principle be  controled rather independently.
"""]))

# ╔═╡ 0d75d077-53c3-4f03-99d3-af796d973d66
md"""
## The Newton-Kantorovich proof
"""

# ╔═╡ b6cfde9d-357b-43e9-82ec-f88e8cb7ba8e
md"""
We choose a value for the weight $\nu$ in the Banach space norm.
"""

# ╔═╡ 671babb1-4b04-4d6a-8dba-79f807114d63
nu=1.1

# ╔═╡ 83f772ea-aeec-49a3-baa2-d2fd8771c0a7
md"""
Since we want to prove something we need to resort to interval arithmetic. Hence we will need to convert $\bar{a}$, $\beta$ and $A$ to interval-valued variables. 
"""

# ╔═╡ 1ace9001-93ae-46f2-98e6-7b371221beb3
begin
	Ia0=interval(a0)
	Ibeta=beta
	IAN=interval(AN)
	Inu=interval(nu)
	# for testing without intervals:	
	 # Ia0=a0
	 # Ibeta=Fbeta
	 # IA=AN
	 # Inu=nu
end

# ╔═╡ 9433b6f1-3460-48c4-b587-9009823327c0
md"""
We are now in the Newton-Kantorovich context and thus need to derive the appropriate bounds. 
"""

# ╔═╡ 9fe5120a-824a-4395-afe8-325a36ab9a07
md"""
### Computable expression for $Y$
"""

# ╔═╡ eb2ae1fd-5062-436b-ac50-eb3b169bca42
md"""
This concerns the bound on the residue $A F(\bar{a})$. Determining the residue requires only a finite computation. Indeed, since $(\bar{a}*\bar{a})_n=0$ for $|n|>2N$ (why?) we infer that $\pi^{>2N} F(\bar{a})=0$ and it can be computed explicitly. In turn this implies that $A F(\bar{a})= \pi^{\leq 2N} A \pi^{\leq 2N} F(\bar{a})$. We conclude that $\|A F(\bar{a})\|_X= \|\pi^{\leq 2N} A \pi^{\leq 2N} F(\bar{a})\|_X $ requires a finite computation only, which we can perform with interval arithmetic to get a rigorous upper bound $Y \geq \|A F(\bar{a})\|_X$.
"""

# ╔═╡ 3021e7c0-508a-4a95-9c16-098fa35b65a0
function boundY(a, beta, A, nu)
	# deal with sizes
    Na = (length(a) - 1) ÷ 2
    NA = (size(A,1)- 1) ÷ 2
	NAfa=max(NA,2*Na)

	# compute residue
	Fa=F(a,beta,2*Na);
	AFa=zeros(typeof(a[1]),2*NAfa+1)
	AFa[NAfa+1-2*Na:NAfa+1+2*Na]=Fa;
    AFa[NAfa+1-NA:NAfa+1+NA]=A*AFa[NAfa+1-NA:NAfa+1+NA];

	result=nunorm(AFa,nu)

    return result
end

# ╔═╡ 9bf0e55e-de5f-4ec9-92e9-5169d3d14302
IY=boundY(Ia0,Ibeta,IAN,Inu)

# ╔═╡ fc18e3ff-4a39-4aba-a2d2-94521035013e
md"""
### Computable expression for $Z_1$
"""

# ╔═╡ 002c7c0a-2b38-4a31-9e08-a424fa5679b8
md"""
This concerns the bound on the operator $I-A DF(\bar{a})$. In the analysis we will split the bound in a part which is computable by "brute force" and an estimate of the tail. The triangle inequality is often helpful for sch arguments, but in our case the characterisation of the operator norm gives an additional tool. In particular we will use that for any $M \in \mathbb{N}$ we have that 

$\begin{align}
\|\Gamma\|_{B(X)} = \max \{ \|\Gamma \pi^{\leq M}\|_{B(X)} , \|\Gamma \pi^{> M}\|_{B(X)} \}.
\end{align}$

We will make a suitable choice for $M$ later.
"""

# ╔═╡ 907f8fb5-859b-4751-b483-2006ba358d5c
Markdown.MD(Markdown.Admonition("note", "Finite bandwidth",[md"""
It follows from the expression for $DF(a)$ that $DF(\bar{a})$ is an finite bandwidth operator of width $N+1$.
In particular, for any $M \in \mathbb{N}$ we have $\pi^{> M+N} DF(\bar{a}) \pi^{\leq M} =0$ and for any $M\geq N$ we have that $\pi^{\leq M-N} DF(\bar{a}) \pi^{> M} =0$.

It follows from the block-diagonal structure of $A$ that this implies that 

$\begin{align}
A DF(\bar{a}) \pi^{> 2N} =  (A^{\leq N}\pi^{\leq N} + I^{>N} \pi^{>N}) DF(\bar{a}) \pi^{> 2N} =   \pi^{>N} DF(\bar{a}) \pi^{> 2N},
\end{align}$

and this expression does not involve $A^{\leq N}$.

Furthermore, it also follows from the block-diagonal structure of $A$ and finite bandwidth of $DF(\bar{a})$ that 

$\begin{align}
A DF(\bar{a}) \pi^{\leq 2N} = A \pi^{\leq 3N} DF(\bar{a}) \pi^{\leq 2N} = \pi^{\leq 3N} A DF(\bar{a}) \pi^{\leq 2N},
\end{align}$

which can be represented as a matrix and requires only a finite computation to evaluate.
"""]))

# ╔═╡ 3d28b669-3294-47fb-8221-688ec3ca5491
Foldable("""Exercise""", md"""Provide the details proving the statements in the box above.""")

# ╔═╡ d35592ba-c0cf-45fe-94c1-c58c3c0bc939
md"""
By choosing $M=2N$ we can thus split the computation of the norm in two parts:

$\begin{align}
\|I-ADF(\bar{a})\|_{B(X)} = \max \{ \|\pi^{\leq 3N} [I-A DF(\bar{a})] \pi^{\leq 2N}\|_{B(X)} , \| [I-ADF(\bar{a})] \pi^{>2N} \|_{B(X)} \}.
\end{align}$

The latter term is computable thanks to the following lemma.
"""

# ╔═╡ d71cba6e-4f62-48d0-aa17-061e1f6bd744
Markdown.MD(Markdown.Admonition("tip", "Lemma",[md"""
Final lemma:
Assume that $N\geq \sqrt{\beta_3}$. 
Then

$\begin{align}
\| [I-ADF(\bar{a})]\pi^{>2N} \|_{B(X)} = 
\max \left\{ \frac{\|[I-ADF(\bar{a})]e_{2N+1}\|_X}{\|e_{2N+1}\|_X} ,  
\frac{\|[I-ADF(\bar{a})]e_{-(2N+1)}\|_X}{\|e_{-(2N+1)}\|_X} \right\}.
\end{align}$
"""]))

# ╔═╡ ac18ecec-5d4e-4d4d-99f1-d5cefcb1355d
Foldable("""Proof""", md"""**Exercise**. Hint: for $|k| \geq 2N+1$
we have $[I-ADF(\bar{a})]e_k= -2\Lambda^{-1} \bar{a} * e_k$ and it follows from monotonicity of $\lambda_n$ in $|n|$ for $|n|\geq \sqrt{\beta_3}$ that $\frac{\|\Lambda^{-1} \bar{a} * e_k \|_X}{\|e_{k}\|_X}$ decreases monotonically in $|k|$ for $|k|\geq 2N+1$. """)

# ╔═╡ b6ae8c03-57aa-4a45-9d6a-9886421ee934
md"""
Collecting and summarizing the above results, we conclude that 

$\begin{align}
\|I-ADF(\bar{a})\|_{B(X)} = \|\pi^{\leq 3N+1}[I-ADF(\bar{a})]\pi^{\leq 2N+1}\|_{B(X)}.
\end{align}$

The latter operator can be represented by a finite matrix and we can compute the operator norm $\|I-ADF(\bar{a})\|_{B(X)}$ using interval arithmetic to obtain the bound $Z_1$.
"""

# ╔═╡ a18bb200-3166-4c23-bcb9-dee71ff5035a
Foldable("""Do we need to determine the norm of I-ADF(a̅) exactly?""", md"""No, we just need an estimate. But we need a pretty good estimate, since we can only succeed with the contraction argument if $Z_1<1$. A naive estimate like $\|I-ADF(\bar{a})\|_{B(X)} \leq \|I\|_{B(X)}+\|A\|_{B(X)}\|DF(\bar{a})\|_{B(X)}$ will certainly not work out. Indeed, we need to take advantage of the cancellations occuring in $I-ADF(\bar{a})$ due to $A$ being an approximate inverse of $DF(\bar{a})$.

We are often willing to spend quite some computer power on getting a rather sharp bound.
For the current problem we could also get away with a somewhat rougher bound,
such as (provided $N\geq \sqrt{\beta_3}$) 

$\begin{align}
Z_1=\max\left\{ \|\pi^{\leq 3N}[I- ADF(\bar{a})]\pi^{\leq 2N}\|_{B(X)} , \frac{2  \|\bar{a}\|_X}{|\lambda_{N+1}|} \right\} .
\end{align}$

We note that sharpness of the $Z_2$-bound is less of an issue, as it plays a considerably less important role in the radii polynomial inequality. 
""")

# ╔═╡ 4a57ce24-2872-4ef1-ada1-73a64ab4f17a
function boundZ1(a, beta, A, nu)
    Na = (length(a) - 1) ÷ 2
    NA = (size(A,1)- 1) ÷ 2

	# deal with sizes
	Ncol = max(NA,ceil(Integer,sup(sqrt(beta[2])))) + Na + 1
	Nrow = Ncol + Na

	# compute A*DF(a)
	Adfa = DF(a,beta,Nrow)
	Adfa[Nrow+1-NA:Nrow+1+NA,:] = A*Adfa[Nrow+1-NA:Nrow+1+NA,:]

	B = I(2*Nrow+1)-Adfa;
	B = B[:,Nrow+1-Ncol:Nrow+1+Ncol]
	result=nunorm(B,nu)
		
    return result
end

# ╔═╡ 47ef9015-eae8-4b29-807f-d223754701d3
IZ1=boundZ1(Ia0,Ibeta,IAN,Inu)

# ╔═╡ abdca76c-832e-42aa-98fc-8a6c2c03965e
Foldable("""Injectivity of A revisited.""", md"""**Exercise**: Prove that $\|I-A DF(\bar{a})\|_{B(X)}<1$ implies that $A$ is surjective on $X$.

**Corollary**: The block-diagonal structure of $A$ implies that it follows from surjectivity of $A$ onto $X$ that $A^{\leq N}$ is surjective as an operator on $X^{\leq N}$. Since $X^{\leq N}$ is finite dimensional, this implies that $A^{\leq N}$  is injective on $X^{\leq N}$. We conclude that $A$ is also injective on $X$.
""") 

# ╔═╡ c0ce2c85-eaae-4b08-8ef7-2a9297b10e05
md"""
### Computable expression for $Z_2$
"""

# ╔═╡ 6f346c67-879a-4427-afd4-6acdeb17af21
md"""
This concerns the bound on the operator $A[DF(a)-DF(\bar{a})]$. We see that 
$A[DF(a)-DF(\bar{a})]b =2A\Lambda^{-1}(a-\bar{a})*b$, hence by the Banach algebra property and the definition of the operator norm we arrive at

$\begin{align}
\|A[DF(a)-DF(\bar{a})]b\|_X \leq 2 \|A\Lambda^{-1}\|_{B(X)} \|a-\bar{a}\|_X \|b\|_X.
\end{align}$

We thus need to compute (or bound) $2\|A\Lambda^{-1}\|_{B(X)}$.
"""

# ╔═╡ 1c2cf614-faba-4fd4-8ceb-2a5d1373d93a
md"""
Since the operator $\Lambda^{-1}$ is diagonal, it has a natural restriction to $\pi^{\leq N} X$. The next lemma shows that the operator norm of $\|A\Lambda^{-1}\|_{B(X)}$ requires only a finite computation.
"""

# ╔═╡ 62cec2b5-e408-46a4-947a-ebc58745c1ec
Markdown.MD(Markdown.Admonition("tip", "Lemma",[md"""Assume $N\geq\sqrt{\beta_3}$. Then 
$\|A\Lambda^{-1}\|_{B(X)} = \max\{ \|A^{\leq N} \Lambda^{-1}  \|_{B(\pi^{\leq N} X)} , |\lambda_{N+1}|^{-1}  \}$
"""]))

# ╔═╡ 84478c6d-e157-40dc-987b-80c7f4d11a68
Foldable("""Proof""", md"""**Exercise**. Hint: consider $\|A\Lambda^{-1}e_n\|_X$ for $|n|\leq N$ and for $|n|>N$ separately.""")

# ╔═╡ 6688fe5d-0019-485d-9059-bb5bd2be7d6e
function boundZ2(a, beta, A, nu)
	NA = (size(A,1)- 1) ÷ 2

	# bound on finite projection
    lambda = [n^2 - 1im * beta[1] * n - beta[2] for n in -NA:NA]
 	Lambdainv=Diagonal( 1 ./ lambda );
	normALinv=nunorm(A*Lambdainv,nu)

	#tail bound
	Nmax=max(NA+1,ceil(Integer,sup(sqrt(beta[2])))+1)
	abslambda=[abs(n^2 - 1im * beta[1] * n - beta[2]) for n in NA+1:Nmax]
    tailbound=1/minimum(abslambda)

	result=2*max(normALinv,tailbound);
	
    return result
end

# ╔═╡ 0d7146b3-3015-416b-aae8-697ac0850603
IZ2=boundZ2(Ia0,Ibeta,IAN,Inu)

# ╔═╡ 5d0b6b50-0170-4f8b-ad5a-0dda212a00e1
md"""
## Finishing the CAP
"""

# ╔═╡ a9f597fe-4058-47a8-8a39-8d6b36594a96
md"""
We now evaluate the radii polynomial to finish the proof.
"""

# ╔═╡ fceb89d1-127a-41b0-b0ec-ec58763a9caa
begin
	# set proper bounds and absorb case of non-interval-arithmetic
	Y=interval(sup(IY))
	Z1=interval(sup(IZ1))
	Z2=interval(sup(IZ2))
end

# ╔═╡ 60b43603-9937-4869-b5fc-f74c1c18e710
r=interval_of_existence(Y, Z1, Z2, Inf)

# ╔═╡ ee6195ff-70b0-46f9-9911-c17bac67aa82
md"""
## Interval arithmetic in Julia 
"""

# ╔═╡ a16e9d49-7fd2-4594-8b81-f8c758b7efd1
md"""
This does not look great:
"""

# ╔═╡ 64a95b2a-f747-4de1-aca1-c54e43119f64
begin
	[isguaranteed(IY)
	 isguaranteed(IZ1)
	 isguaranteed(IZ2)]
end

# ╔═╡ b0612642-f710-496b-be28-c12d5f867e20
md"""
## Symmetry
"""

# ╔═╡ eb1f2c71-60fc-4ee5-9716-f18dd5eb0d2d
md"""
We made sure that the numerical approximation was conjugate symmetric: $\bar{a}^\dagger=\bar{a}$.
"""

# ╔═╡ 6df91aa3-ad36-4b9e-ba4a-4dfb46eb1d27
	println(conj(reverse(a0))==a0)

# ╔═╡ 49eae38b-786f-45fc-b331-3b01c32363a1
md"""
Hence the following lemma finishes our proof.
"""

# ╔═╡ 35df9eed-5f22-4e26-b119-ea7354d2c762
Markdown.MD(Markdown.Admonition("tip", "Lemma (symmetry of the solution)",
[md"""
Assume that $\bar{a}^\dagger=\bar{a}$ and that the asummptions of the Newton-Kantorovich theorem were satisfied for some $r=r_0>0$. Then the zero $a^*$ of $F$ such that $\|a^*-\bar{a}\|_X \leq r_0$ satisfies $(a^*)^\dagger=a^*$.
"""]))

# ╔═╡ bae570fe-8219-4cf2-b763-fcd9f0f02735
md"""
*Proof*: First, based on the symmetry $F(a^\dagger)=F(a)^\dagger$ of the zero finding problem, we see that $F(a^*)=0$ implies that $F((a^*)^\dagger)=(F(a^*))^\dagger=0$, hence $(a^*)^\dagger$ is a zero of $F$.  Second, the norm on $X$ is such that $\|a^\dagger\|_X=\|a\|_X$  for all $a\in X$. We find that 

$\begin{align} 
\|(a^*)^\dagger-\bar{a}\|_X=\|(a^*)^\dagger-\bar{a}^\dagger\|_X 
=\|(a^*-\bar{a})^\dagger\|_X = \|a^*-\bar{a}\|_X \leq r_0.
\end{align}$

Hence $(a^*)^\dagger$ is a zero of $F$ inside the ball of radius $r_0$ around $\bar{a}$, and since this zero is unique by the Newton-Kantorovich theorem, we conclude that $(a^*)^\dagger=a^*$.
"""

# ╔═╡ 09922158-6ac9-410e-90ac-5240a9056086
md"""
## Conclusion
"""

# ╔═╡ 3cc0fc2c-9512-41d6-9798-4e8cbec94c9e
md"""This means that we have proven that there is a periodic solution of the differential equation within distance"""

# ╔═╡ 670dd539-75b4-45b5-b9f3-de4354945648
r₀=inf(r)

# ╔═╡ f4f432f8-ea80-4169-94ad-ee3289e1a93f
md"""
of the approximation depicted in the figure below (with distance measured in the max-norm on $[0,2\pi]$). 
"""

# ╔═╡ 26d06cfe-24cc-41ec-87d2-8e489865b5b8
begin
	plot(t, u, legend=false, width=2)
	xlabel!("t")
	ylabel!("u")
	xlims!(0, 2*pi)
	middle=(maximum(u)+minimum(u))/2
	annotate!(pi,middle, Plots.text("Capified",50, :magenta, rotation = 30 ))
end

# ╔═╡ e40d36d5-93a6-4b76-bccd-c629f621b22b
Foldable("""Error bound""", md"""
**Exercise**: Prove that the error between the solution $u^*(t)=\sum_{n\in \mathbb{Z}} a^*_n e^{int}$ and the approximation $\bar{u}(t)=\sum_{|n| \leq N} \bar{a}_n e^{int}$ is indeed controled by

$\begin{align}
\max_{t \in [0,2\pi]} | u^*(t) - \bar{u}(t)| \leq r_0
\end{align}$
""")

# ╔═╡ 194e7df2-a234-45c7-8760-8c49a0d6e651


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
RadiiPolynomial = "f2081a94-c849-46b6-8dc9-07bb90ed72a9"
ToeplitzMatrices = "c751599d-da0a-543b-9d20-d0a503d91d24"

[compat]
Plots = "~1.40.5"
PlutoTeachingTools = "~0.2.15"
PlutoUI = "~0.7.59"
RadiiPolynomial = "~0.8.12"
ToeplitzMatrices = "~0.8.4"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.4"
manifest_format = "2.0"
project_hash = "db3d94c0ad35dd4b82a364910d1cea9158e423fd"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

    [deps.AbstractFFTs.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

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

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9e2a6b69137e6969bab0152632dcb3bc108c8bdd"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+1"

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

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

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

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

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
git-tree-sha1 = "dcb08a0d93ec0b1cdc4af184b26b591e9695423a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.10"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c6317308b9dc757616f0b5cb379db10494443a7"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.2+0"

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

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "0653c0a2396a6da5bc4766c43041ef5fd3efbe57"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.11.0"

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

    [deps.FillArrays.weakdeps]
    PDMats = "90014a1f-27ba-587c-ab20-58faa44d9150"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

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
git-tree-sha1 = "5c1d8ae0efc6c2e7b1fc502cbe25def8f661b7bc"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.2+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1ed150b39aebcc805c26b93a8d0122c940f64ce2"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.14+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "3f74912a156096bd8fdbef211eff66ab446e7297"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.0+0"

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

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "14eb2b542e748570b56446f4c50fbfb2306ebc45"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2024.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IntervalArithmetic]]
deps = ["CRlibm_jll", "MacroTools", "RoundingEmulator"]
git-tree-sha1 = "433b0bb201cd76cb087b017e49244f10394ebe9c"
uuid = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
version = "0.22.14"

    [deps.IntervalArithmetic.extensions]
    IntervalArithmeticDiffRulesExt = "DiffRules"
    IntervalArithmeticForwardDiffExt = "ForwardDiff"
    IntervalArithmeticRecipesBaseExt = "RecipesBase"

    [deps.IntervalArithmetic.weakdeps]
    DiffRules = "b552c78f-8df3-52c6-915a-8e097449b14b"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

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

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

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

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

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

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

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

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

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
version = "0.3.23+4"

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

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "cb5a2ab6763464ae0f19c86c56c63d4a2b0f5bda"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.52.2+0"

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

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "492601870742dcd38f233b23c3ec629628c1d724"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.7.1+1"

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
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "Pkg", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "85ddd93ea15dcd8493400600e09104a9e94bb18d"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.5.15"

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

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

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

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "2f5d4697f21388cbe1ff299430dd169ef97d7e14"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.4.0"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

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

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

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

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

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
version = "5.8.0+1"

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
# ╠═2661bfc9-e398-41ed-87d9-c78f05da64cb
# ╠═018ecc45-8638-4a59-b561-efb086bdc751
# ╟─d5a510a3-c518-47ed-96bc-7bb22e3b08b5
# ╟─f73b88cb-19e2-4d50-a45e-fbe90fb691cd
# ╟─cf717e9b-26d8-449d-9f3b-d1d5d0bd2bb5
# ╟─c33dc650-3f94-11ef-398a-8bbc4a2b69b8
# ╟─58892dec-c89c-4f1f-9f97-700bc27ef338
# ╟─1bfcc76a-266a-4784-b76b-3871445e603e
# ╟─3e3a3a60-f76f-401d-b265-a99e713d98c5
# ╟─8f3e96dc-0cdc-411e-8f94-e9e47488abf3
# ╟─7c9ca144-02f1-4941-bb0b-becb96b4fc8a
# ╟─b5d97ade-bcde-4c7a-ae78-697e73c10a0c
# ╟─47e39d6a-231e-40d1-906e-eeafe1b1f807
# ╟─ac3e3cf5-e5c8-4863-8c7b-ccd5f58a3073
# ╠═4a8b8e9e-9562-4cb5-a029-856127e84b98
# ╟─312d9592-1be0-4626-8773-51cc8604dd6a
# ╟─c4d9c4f0-fd89-4a2a-96a7-17c4ba6dd747
# ╟─766164f3-eef3-4798-bf7c-afdca52d3ba0
# ╟─bbf2e943-8c5c-469f-a2f3-5f8b758f6638
# ╟─f91abfa7-8a8f-4204-95d9-637c622dd720
# ╟─8f6431ed-346c-4c51-83c4-ea7ec12a8d60
# ╟─c5cfe80e-1728-4393-b837-187e3e7d49a9
# ╟─39d13194-bca9-47c1-a686-2ea5f9d62289
# ╟─7e78028c-f70b-4474-8bee-783ad7d99d56
# ╟─c7a10c87-a95d-4b0b-b63d-f9546e97bb36
# ╠═efb766e3-e41f-4113-ae6e-c4dd7766c784
# ╠═0f04a415-b4c1-4b05-a35a-6d3029decb12
# ╟─31a59e3d-32af-4328-b3bc-cc434728429d
# ╟─8d7fb15d-32d3-4fbb-9d06-c8ec1ed7d00e
# ╟─b2defc79-da78-47c8-bb78-44cb9000ff58
# ╟─9ab6ce17-0637-4008-be29-e0bf4fda4287
# ╟─9002b709-ec9c-4063-a3a8-7b0c00834095
# ╟─c2f4618c-dd22-477c-8fe7-324c82d3ff46
# ╟─7d79ed91-7cbb-48b4-966a-d3165f2a3d89
# ╟─66deaecf-e201-40ec-88d3-99574195253a
# ╟─c32d9839-7cb7-4018-bab1-b2211f05b71e
# ╟─bd76555d-aa30-4313-b4c4-c3f4b9000194
# ╟─740167dc-5bcd-48c0-b1f0-a7bf5e1b6cb7
# ╠═026aec24-4ed0-4189-af4f-28f56e6964ef
# ╟─06f3d467-af63-440d-b486-8da7f67c314f
# ╟─f5cc03ca-09eb-4419-9d51-78c956ee73d9
# ╟─c5706263-a08a-4742-b376-b980232d6191
# ╠═e5b2628d-2e64-47ca-b594-a0e4cb8ea64a
# ╟─c4363fdf-6500-4830-bc42-9777d616ac34
# ╟─24e1a5af-569e-4ae8-a3d7-19b7962db48b
# ╟─77b91c27-ddb0-4115-b0d0-a78f2729dadc
# ╟─51e35c7e-77f2-4911-914c-7f1a242c3b91
# ╟─4a9daa68-86ed-4327-98fc-03d436885d2e
# ╟─3842ad29-ae9c-4743-b487-e76d94fb01cf
# ╠═1688b3d8-8228-45ee-8729-22265a86484a
# ╟─aa150fa3-1729-43b1-bb55-73c741221760
# ╠═e1f6ef9b-5818-415a-8945-97b65288a433
# ╟─d0858eee-5086-4754-966f-de2940183573
# ╠═a89606ee-b46f-4656-9cb9-d89850646593
# ╟─192507bd-8285-4581-8a6d-f9ff1b2c4793
# ╠═301bf3ed-b180-4c44-a92f-0f9f4c2509ed
# ╟─aa320281-42ce-4035-8c55-1f44661737e2
# ╠═d5800e45-8809-46d6-9f41-7565e5b1cbfd
# ╟─5f3d91ee-14c2-4b5e-8e2a-a607d086eb51
# ╟─936423fd-0d78-41da-8b49-5fca3c125547
# ╠═e37811fa-cd67-424e-8707-25ad8534d293
# ╟─0b4e2571-cbec-4a4d-aac6-dd8b32e3f474
# ╟─0d75d077-53c3-4f03-99d3-af796d973d66
# ╟─b6cfde9d-357b-43e9-82ec-f88e8cb7ba8e
# ╠═671babb1-4b04-4d6a-8dba-79f807114d63
# ╟─83f772ea-aeec-49a3-baa2-d2fd8771c0a7
# ╠═1ace9001-93ae-46f2-98e6-7b371221beb3
# ╟─9433b6f1-3460-48c4-b587-9009823327c0
# ╟─9fe5120a-824a-4395-afe8-325a36ab9a07
# ╟─eb2ae1fd-5062-436b-ac50-eb3b169bca42
# ╠═3021e7c0-508a-4a95-9c16-098fa35b65a0
# ╠═9bf0e55e-de5f-4ec9-92e9-5169d3d14302
# ╟─fc18e3ff-4a39-4aba-a2d2-94521035013e
# ╟─002c7c0a-2b38-4a31-9e08-a424fa5679b8
# ╟─907f8fb5-859b-4751-b483-2006ba358d5c
# ╟─3d28b669-3294-47fb-8221-688ec3ca5491
# ╟─d35592ba-c0cf-45fe-94c1-c58c3c0bc939
# ╟─d71cba6e-4f62-48d0-aa17-061e1f6bd744
# ╟─ac18ecec-5d4e-4d4d-99f1-d5cefcb1355d
# ╟─b6ae8c03-57aa-4a45-9d6a-9886421ee934
# ╟─a18bb200-3166-4c23-bcb9-dee71ff5035a
# ╠═4a57ce24-2872-4ef1-ada1-73a64ab4f17a
# ╠═47ef9015-eae8-4b29-807f-d223754701d3
# ╟─abdca76c-832e-42aa-98fc-8a6c2c03965e
# ╟─c0ce2c85-eaae-4b08-8ef7-2a9297b10e05
# ╟─6f346c67-879a-4427-afd4-6acdeb17af21
# ╟─1c2cf614-faba-4fd4-8ceb-2a5d1373d93a
# ╟─62cec2b5-e408-46a4-947a-ebc58745c1ec
# ╟─84478c6d-e157-40dc-987b-80c7f4d11a68
# ╠═6688fe5d-0019-485d-9059-bb5bd2be7d6e
# ╠═0d7146b3-3015-416b-aae8-697ac0850603
# ╟─5d0b6b50-0170-4f8b-ad5a-0dda212a00e1
# ╟─a9f597fe-4058-47a8-8a39-8d6b36594a96
# ╠═fceb89d1-127a-41b0-b0ec-ec58763a9caa
# ╠═60b43603-9937-4869-b5fc-f74c1c18e710
# ╟─ee6195ff-70b0-46f9-9911-c17bac67aa82
# ╟─a16e9d49-7fd2-4594-8b81-f8c758b7efd1
# ╠═64a95b2a-f747-4de1-aca1-c54e43119f64
# ╟─b0612642-f710-496b-be28-c12d5f867e20
# ╟─eb1f2c71-60fc-4ee5-9716-f18dd5eb0d2d
# ╠═6df91aa3-ad36-4b9e-ba4a-4dfb46eb1d27
# ╟─49eae38b-786f-45fc-b331-3b01c32363a1
# ╟─35df9eed-5f22-4e26-b119-ea7354d2c762
# ╟─bae570fe-8219-4cf2-b763-fcd9f0f02735
# ╟─09922158-6ac9-410e-90ac-5240a9056086
# ╟─3cc0fc2c-9512-41d6-9798-4e8cbec94c9e
# ╠═670dd539-75b4-45b5-b9f3-de4354945648
# ╟─f4f432f8-ea80-4169-94ad-ee3289e1a93f
# ╟─26d06cfe-24cc-41ec-87d2-8e489865b5b8
# ╟─e40d36d5-93a6-4b76-bccd-c629f621b22b
# ╠═194e7df2-a234-45c7-8760-8c49a0d6e651
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

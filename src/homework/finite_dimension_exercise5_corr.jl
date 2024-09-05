### A Pluto.jl notebook ###
# v0.19.43

#> [frontmatter]
#> homework_number = 5
#> order = 5.5
#> title = "Rigorous control of the entire spectrum - correction"
#> tags = ["module1", "homeworks"]
#> layout = "layout.jlhtml"

using Markdown
using InteractiveUtils

# ╔═╡ 9099f41e-6239-4f7f-a3ec-82e7d4787f8f
using PlutoTeachingTools

# ╔═╡ 2661bfc9-e398-41ed-87d9-c78f05da64cb
using RadiiPolynomial, LinearAlgebra

# ╔═╡ 7fc40507-eda3-474d-a454-04e9173a7adb
html"""<style>
main {
    max-width: 1000px;
    margin-left: auto;
    margin-right: auto;
    text-align: justify;
}
"""

# ╔═╡ 45e10c91-f161-43a5-9c9e-5b4dca6a8e53
md"""
Exercise 3 can be used to enclose eigenvalues of a given matrix one by one. We now present an alternate strategy to enclose the entire spectrum at once (but not the corresponding eigenvectors), which can sometimes also be adapted in infinite dimension. This strategy relies on the Gershgorin circle theorem, which we recall below.
"""

# ╔═╡ a6e349ff-ba01-48c4-91a8-5400faa05346
Markdown.MD(Markdown.Admonition("tip", "Theorem", [md"For any matrix $A=\left(A_{i,j}\right)_{1\leq i,j\leq d}$, its spectrum $\sigma(A)$ is controlled as follows.

$\begin{align}
\sigma(A) \subset \bigcup_{i=1}^d D\left(A_{i,i},\, \sum_{j\neq i} \vert A_{i,j}\vert\right),
\end{align}$

where $D(z,r)$ denotes the closed disk of center $z$ and radius $r$ in the complex plane. Moreover, if $I \subset \{1,\ldots,d\}$ is such that $\bigcup_{i \in I} D\left(A_{i,i},\, \sum_{j \ne i} \vert A_{i,j} \vert\right)$ is disjoint from $\bigcup_{i \notin I} D\left(A_{i,i},\, \sum_{j \ne i} \vert A_{i,j} \vert\right)$, then $\bigcup_{i \in I} D\left(A_{i,i},\, \sum_{j \ne i} \vert A_{i,j} \vert\right)$ contains exactly $\vert I \vert$ eigenvalues."]))

# ╔═╡ d01c5817-c4b6-4502-81f6-51ae5715117b
md"""
**1.** Using the Gershgorin circle theorem, get as tight as possible rigorous enclosures of all eigenvalues of $W_{3}$ (defined in the third exercise)).
"""

# ╔═╡ 02939955-c9aa-4152-8e08-f31e1e0c0e9c
Foldable("Hint",
md"You may first compute numerically a matrix $P$ of approximate eigenvectors of $W_3$, then rigorously compute $\tilde W_{3} = P^{-1}W_3 P$, and finally apply the Gershgoring circle theorem to $\tilde W_{3}$."
)

# ╔═╡ 257d92b2-ccc4-4351-b3ac-ff53f6f2d54d
function my_rigorous_inv(A)
	B = inv(A)
	δ = opnorm(I - LinearOperator(interval.(A))*LinearOperator(interval.(B)), 1)
	if sup(δ) < 1
		r = δ/(1-δ) * opnorm(LinearOperator(interval.(B)), 1)
	else
		println("Unable to rigorously invert")
		r = NaN
	end
	return interval(B, sup(r)*ones(size(B)); format = :midpoint)
end

# ╔═╡ 773d5d67-2819-4225-add4-f7913fd6d296
function gershgorin(M)
	eigenvalues, P = eigen(M)
	iP = interval.(P)
	iPinv = my_rigorous_inv(P)
	M̃ = iPinv * M * iP
	centers = mid.(diag(M̃))
	temp = M̃
	for n = 1:size(M,1)
		temp[n,n] = 0
	end
	radii = sum(abs.(temp); dims=2) + radius.(diag(M̃))
	return interval(centers, sup.(radii); format = :midpoint)
end

# ╔═╡ 0db57db5-c43a-451a-9ebb-eccd799514ef
function W(N)
	M = zeros(2N+1, 2N+1)
	for i = 1:2N+1
		M[i,i] = abs(N - i + 1)
		if i+1 ≤ 2N+1
			M[i,i+1] = 1
		end
		if i-1 ≥ 1
			M[i,i-1] = 1
		end
	end
	return M
end

# ╔═╡ 6339a6ca-985b-4cc0-891f-72b943b5b812
N = 3

# ╔═╡ ad19ec24-3779-453d-9315-ce70db1b9b4e
spectrum = gershgorin(W(N))

# ╔═╡ ab8f264a-c1a4-40d5-a38c-dc1a5ef398a6
function count_nb_neg_eig(spectrum)
	nb_neg = sum(isstrictless.(spectrum, interval(0)))
	nb_pos = sum(isstrictless.(interval(0), spectrum))
	if nb_neg + nb_pos == length(spectrum)
		println("We have proven that there is exactly ",nb_neg," negative eigenvalue(s)")
	else
		println("We have proven that there are at least ",nb_neg," negative eigenvalue(s) and at least ",nb_pos," positive eigenvalue(s)")
	end
end

# ╔═╡ c308fcbf-d094-4f5b-9cb5-23c17a6c6b73
count_nb_neg_eig(spectrum)

# ╔═╡ bc052916-db0e-43a2-bbb0-76b917d6f638
md"""
**2.** Try to prove that $W_{1000}$ has exactly one eigenvalue with negative real part.
"""

# ╔═╡ 3686e5b2-daac-470d-b9e3-4bd5706e5894
Foldable("Hint",md"You do not need to numerically diagonalize all of $W_{1000}$: for most rows, the corresponding Gershgorin disk already lies in the left half of the complex plane. ")

# ╔═╡ a75b3d18-e03e-4fc6-a5e0-0c03579e53a9
function gershgorin_wilkinson_cheap(N, d)
	Md = W(d)
	eigenvalues, P = eigen(Md)
	iPd = interval.(Matrix(I(2*d+3)))
	iPd[2:2*d+2,2:2*d+2] = interval.(P)
	iPdinv = interval.(Matrix(I(2*d+3)))
	iPdinv[2:2*d+2,2:2*d+2] = my_rigorous_inv(P)
	M̃ = interval.(W(N))
	M̃[N-(d+1).+(1:2*d+3),N-(d+1).+(1:2*d+3)] = iPdinv * M̃[N-(d+1).+(1:2*d+3),N-(d+1).+(1:2*d+3)] * iPd
	centers = mid.(diag(M̃))
	temp = M̃
	for n = 1:2*N+1
		temp[n,n] = 0
	end
	radii = sum(abs.(temp); dims=2) + radius.(diag(M̃))
	return interval(centers, sup.(radii); format = :midpoint)
end

# ╔═╡ eb4abb7a-9605-4099-ba30-716af7e2d173
spectrum_approx = gershgorin_wilkinson_cheap(1000, 5)

# ╔═╡ cdcdea86-7463-47d2-9416-930f816caad5
count_nb_neg_eig(spectrum_approx)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
RadiiPolynomial = "f2081a94-c849-46b6-8dc9-07bb90ed72a9"

[compat]
PlutoTeachingTools = "~0.2.15"
RadiiPolynomial = "~0.8.12"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.2"
manifest_format = "2.0"
project_hash = "49d29ae341b6f6390131461e857efcdfb8528a3c"

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
git-tree-sha1 = "c0216e792f518b39b22212127d4a84dc31e4e386"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.3.5"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.0+0"

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
git-tree-sha1 = "8b72179abc660bfab5e28472e019392b97d0985c"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.4"

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

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "a6adc2dcfe4187c40dc7c2c9d2128e326360e90a"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.32"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "e0b5cd21dc1b44ec6e64f351976f961e6f31d6c4"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.3"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
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
git-tree-sha1 = "eeaedcf337f33c039f9f3a209a8db992deefd7e9"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "2.4.8"

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
git-tree-sha1 = "ab55ee1510ad2af0ff674dbcced5e94921f867a9"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.59"

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
git-tree-sha1 = "89f57ab86310e5ca7009cb236441505ba6b3242a"
uuid = "f2081a94-c849-46b6-8dc9-07bb90ed72a9"
version = "0.8.12"

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

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

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
# ╟─9099f41e-6239-4f7f-a3ec-82e7d4787f8f
# ╠═2661bfc9-e398-41ed-87d9-c78f05da64cb
# ╟─45e10c91-f161-43a5-9c9e-5b4dca6a8e53
# ╟─a6e349ff-ba01-48c4-91a8-5400faa05346
# ╟─d01c5817-c4b6-4502-81f6-51ae5715117b
# ╟─02939955-c9aa-4152-8e08-f31e1e0c0e9c
# ╠═773d5d67-2819-4225-add4-f7913fd6d296
# ╠═257d92b2-ccc4-4351-b3ac-ff53f6f2d54d
# ╠═0db57db5-c43a-451a-9ebb-eccd799514ef
# ╠═6339a6ca-985b-4cc0-891f-72b943b5b812
# ╠═ad19ec24-3779-453d-9315-ce70db1b9b4e
# ╠═ab8f264a-c1a4-40d5-a38c-dc1a5ef398a6
# ╠═c308fcbf-d094-4f5b-9cb5-23c17a6c6b73
# ╟─bc052916-db0e-43a2-bbb0-76b917d6f638
# ╟─3686e5b2-daac-470d-b9e3-4bd5706e5894
# ╠═a75b3d18-e03e-4fc6-a5e0-0c03579e53a9
# ╠═eb4abb7a-9605-4099-ba30-716af7e2d173
# ╠═cdcdea86-7463-47d2-9416-930f816caad5
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

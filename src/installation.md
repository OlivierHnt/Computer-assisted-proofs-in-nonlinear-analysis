---
title: "Software installation"
tags: ["welcome"]
order: 1
layout: "md.jlmd"
---

$(
    begin
        # these special elements will automatically update to read the latest Julia version (cf. the JavaScript snippet at the bottom of this page)

        version = html"<auto-julia-version>1.10.0</auto-julia-version>"
        pkg_version = html"<auto-julia-version short>1.10</auto-julia-version>"

        nothing
    end
)

## Step 1: Install Julia

Download the current stable release of Julia on [https://julialang.org/downloads](https://julialang.org/downloads), using the correct version for your operating system (Linux x86, Mac, Windows, etc).

## Step 2: Install RadiiPolynomial

Next we will install [RadiiPolynomial](https://github.com/OlivierHnt/RadiiPolynomial.jl), the library for computer-assisted proofs that we will be using during the course. RadiiPolynomial is a Julia package designed for computer-assisted proofs relying on spectral methods and the contraction mapping theorem.

Open the **Julia REPL**. This is the command-line interface to Julia: there you type a command, and when you press ENTER, it runs, and you see the result.

To install RadiiPolynomial, we want to run a _package manager command_. To switch from _Julia mode_ to _Pkg mode_, type `]` (closing square bracket) at the <code><span style="color: #399746;">julia></span></code> prompt:

<pre><code>
<span style="color: #399746;">julia></span> ]

<span style="color: #1f83ff;">(&#64;v$(pkg_version)) pkg></span>
</code></pre>

The line turns blue and the prompt changes to <code><span style="color: #1f83ff;">pkg></span></code>, telling you that you are now in _package manager mode_. This mode allows you to do operations on **packages** (also called libraries).

Run the following (case sensitive) command to *add* (install) the package to your system by downloading it from the internet (you should only need to do this *once* for each installation of Julia):

<pre><code>
<span style="color: #1f83ff;">(&#64;v$(pkg_version)) pkg></span> add RadiiPolynomial
</code>
</pre>





<script defer>
const run = f => f();
run(async () => {
const versions = await (await fetch(`https://julialang-s3.julialang.org/bin/versions.json`)).json()
const sortby = v => v.split("-")[0].split(".").map(parseFloat).reduce((a,b) => a*10000 + b)
const version_names = Object.keys(versions).sort((a,b) => sortby(a) - sortby(b)).reverse()
const stable = version_names.find(v => versions[v].stable)
console.log({stable})
const pkg_stable = /\\d+\\.\\d+/.exec(stable)[0]
document.querySelectorAll("auto-julia-version").forEach(el => {
    console.log(el)
    el.innerText = el.getAttribute("short") == null ? stable : pkg_stable
})
});
</script>

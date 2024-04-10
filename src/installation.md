---
title: "Software installation"
tags: ["welcome"]
order: 2
layout: "md.jlmd"
---

$(
    begin
        # these special elements will automatically update to read the latest Julia version. See the JavaScript snippet at the bottom of this page to see how it works!

        version = html"<auto-julia-version>1.10.0</auto-julia-version>"
        pkg_version = html"<auto-julia-version short>1.10</auto-julia-version>"

        nothing
    end
)

# Install Julia & RadiiPolynomial

## Step 1: Install Julia

Go to [https://julialang.org/downloads](https://julialang.org/downloads) and download the current stable release, using the correct version for your operating system (Linux x86, Mac, Windows, etc).

## Step 2: Run Julia

After installing, **make sure that you can run Julia**. On some systems, this means searching for the Julia program installed on your computer; in others, it means running the command `julia` in a terminal.

*Make sure that you are able to launch Julia before proceeding!*

## Step 3: Install RadiiPolynomial

Next we will install [RadiiPolynomial](https://github.com/OlivierHnt/RadiiPolynomial.jl), the library for computer-assisted proofs that we will be using during the course. RadiiPolynomial is a Julia package designed for rigorous spectral methods.

Open the **Julia REPL**. This is the command-line interface to Julia, similar to the previous screenshot.
There you type a command, and when you press ENTER, it runs, and you see the result.

To install RadiiPolynomial, we want to run a _package manager command_. To switch from _Julia_ mode to _Pkg_ mode, type `]` (closing square bracket) at the `julia>` prompt:

<pre><code>
julia> ]

(&#64;v$(pkg_version)) pkg>
</code></pre>

The line turns blue and the prompt changes to `pkg>`, telling you that you are now in _package manager mode_. This mode allows you to do operations on **packages** (also called libraries).

Then, run the following (case sensitive) command to *add* (install) the package to your system by downloading it from the internet (you should only need to do this *once* for each installation of Julia):

<pre><code>
(&#64;v$(pkg_version)) pkg> add RadiiPolynomial
</code></pre>

You can now close the terminal.

## Step 4: Running RadiiPolynomial

Repeat the following steps whenever you want to work on a project.

Start the Julia REPL, like you did during the setup. In the REPL, type:
```julia
julia> using RadiiPolynomial
```





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

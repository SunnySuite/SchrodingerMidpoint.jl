Companion code to the paper:

* D. Dahlbom et al., _Geometric integration of classical spin dynamics via a mean-field Schrödinger equation_, [arXiv:2204.07563](https://arxiv.org/abs/2204.07563).

To install this package from within the Julia terminal, enter
```
using Pkg
Pkg.develop(url="https://github.com/SunnySuite/SchrodingerMidpoint.jl")
```
This will `git clone` the package into the `.julia/dev/` directory. Source file modifications will be automatically picked up by Julia.

To run the package tests,
```
Pkg.test("SchrodingerMidpoint")
```

To generate paper figures,
```
using SchrodingerMidpoint

# Each of these commands should bring up a plot window
fig1()
fig2()
fig3()
fig4()
```

---
Authors: Kipton Barros, David Dahlbom
# AUrisim-julia
# An Oritatami systems simulator written in Julia

There are two notebooks now, a Jupyter notebook and a Pluto.jl notebook.

Both require Julia: 
https://julialang.org/

The Jupyter notebook requires Jupyter: 
https://jupyter.org/

The Pluto notebook uses the Pluto.jl package in Julia, so does not need external installation, but you need to install it along with the other package requirements.
I will mostly focus on the Pluto notebook from now on, probably, because of its reactive nature and because the source file in the end is a valid Julia file, as opposed to the jupyter notebooks.

To use the Oritatami-pluto notebook, you need to install the packages below, after installing Julia:
 - Pluto
 - PlutoUI
 - Plots
 - Colors
 - ColorSchemes
 - LightGraphs
 - GraphPlot

To install the packages, you need to run
> using Pkg

Then, you can install the package called PackageName with the command
> Pkg.add("PackageName")

After installing the packages, the following two commands run Pluto, where you can open the notebook **Oritatami-pluto.jl** and start folding :)
> import Pluto
> Pluto.run()


Julia packages used in the Jupyter notebook (can be installed from the notebook by uncommenting the relevant lines in the first cell):
- Plots
- Makie (through a backend, e.g. WGLMakie)
- Colors
- ColorSchemes


For now there is no documentation, but it will come when I have more time.


![Image of BinaryCounter](https://github.com/szfazekas/AUrisim-julia/blob/main/counter1k.gif)

# AUrisim-julia
# An Oritatami systems simulator written in Julia

The goal of this project is to create an easy to use interactive simulator for Oritatami systems.
These systems were introduced by Cody Geary, Olivier Meunier, Nicolas Schabanel and Shinnosuke Seki as a mathematical model of engineering shapes by exploiting the cotransciptional of RNA. The papers that started the field are [Geary14](https://science.sciencemag.org/content/345/6198/799.abstract) and [GMSS16](https://drops.dagstuhl.de/opus/volltexte/2016/6456/pdf/LIPIcs-MFCS-2016-43.pdf). For a recent paper, see, e.g., the one on the  [infinite binary counter](https://link.springer.com/chapter/10.1007/978-3-030-38919-2_46), which can be loaded from the file bincount.auri.txt.

There are two notebooks now, a Jupyter notebook and a Pluto.jl notebook.

Both require Julia: 
https://julialang.org/

The Jupyter notebook requires Jupyter: 
https://jupyter.org/

The Pluto notebook uses the Pluto.jl package in Julia, so does not need external installation, but you need to install it along with the other package requirements.
I will mostly focus on the Pluto notebook from now on, probably, because of its reactive nature and because the source file in the end is a valid Julia file, as opposed to the Jupyter notebooks. In Pluto I found it much easier to set up the interactive controls, which is a key in this project. After opening the notebook and all cells finished loading you should see something like below (with fewer beads in the conformation).
![Image of Pluto notebook](https://github.com/szfazekas/AUrisim-julia/blob/main/Screenshot1_pluto.png)


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

-------------------------------------------------------

Julia packages used in the Jupyter notebook (can be installed from the notebook by uncommenting the relevant lines in the first cell):
- Plots
- Makie (through a backend, e.g. WGLMakie)
- Colors
- ColorSchemes


For now there is no documentation, but it will come when I have more time.


![Image of BinaryCounter](https://github.com/szfazekas/AUrisim-julia/blob/main/counter1k.gif)

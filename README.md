# AMLStateSpace
This repository contains code to generate figures for the paper:
> Henrik Lilljebjörn, Pablo Peña-Martínez, Hanna Thorsson, Rasmus Henningsson, Marianne Rissler, Niklas Landberg, Noelia Puente-Moncada, Sofia von Palffy, Vendela Rissler, Petr Stanek, Jonathan Desponds, Xiangfu Zhong, Gunnar Juliusson, Vladimir Lazarevic, Sören Lehmann, Magnus Fontes, Helena Ågerstam, Carl Sandén, Christina Orsmark-Pietras, Thoas Fioretos. "The cellular state space of AML unveils novel NPM1 subtypes with distinct clinical outcomes and immune evasion properties".

## General plots
Most plots were generated in R.
See [R/README.md](R/README.md) for details on how to run the code to create them.

## Projection plots
The Projection plots were generated in [Julia](www.julialang.org) using the [SingleCellProjections.jl](https://github.com/BioJulia/SingleCellProjections.jl) package.
See [ProjectionPlots/README.md](ProjectionPlots/README.md) for details on how to run the code to create these plots.

## Data
The data files needed to generate the plots are available for download at https://doi.org/10.17044/scilifelab.23715648.
(For the projection plots, the required files can be downloaded automatically, see [ProjectionPlots/README.md](ProjectionPlots/README.md) for more info.)

# ProjectionPlots

Following these instructions will ensure that you use the right version of Julia and of the all Julia packages.
This is important to guarantee reproducibility.

The original figures were produced on a Linux computer.
(If you have problems reproducing the results, consider using Linux.)

Note that these plots were generated using an older version of the [SingleCellProjections.jl](https://github.com/BioJulia/SingleCellProjections.jl) package. For new projects we recommend using the latest version, which is available in the Julia registry.


## Installation
Install Julia 1.10 using [juliaup](https://julialang.org/install/), by running
```
juliaup add 1.10
```
in a terminal.

## Running
In this folder, run
```
julia +1.10 --project
```
in a terminal to start `Julia 1.10` and load the right project.

In the Julia REPL, run
```julia-repl
julia> using ProjectionPlots
```
to load the package.

To download the sample files needed into the `data/samples` folder, just run:
```julia-repl
julia> download_samples()
```

Alternatively, you can download all the `.h5` files manually from https://doi.org/10.17044/scilifelab.23715648 into the `data/samples` folder or, if you want to use a custom location for the data, run:
```julia-repl
julia> set_samples_path("path/to/folder/with/AMLXYZ.h5")
```
The setting will be persisted, so you only need to do this once.

Now, you can create figures by calling:
```julia-repl
julia> figure2a()

julia> figure3b()

julia> supplementary_figure6()

julia> supplementary_figure7()

julia> supplementary_figure8()

julia> supplementary_figure11b()
```
(The source code for these functions are located in `src/main_AML.jl`.)

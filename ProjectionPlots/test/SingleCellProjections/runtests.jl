using ProjectionPlots.SingleCellProjections
using Test
using StableRNGs
using StaticArrays
using Statistics

import .SingleCellProjections: BarnesHutTree, build!


include("MatrixExpressions/runtests.jl")
include("SingleCell10x/runtests.jl")
include("SCTransform/runtests.jl")


@testset "SingleCellProjections.jl" begin
    include("test_barnes_hut.jl")
end

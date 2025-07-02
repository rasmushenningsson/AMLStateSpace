module ProjectionPlots

include("LazyCache/LazyCache.jl")
include("SingleCellProjections/SingleCellProjections.jl")


using .SingleCellProjections
using .SingleCellProjections.SingleCell10x
using .LazyCache
using LinearAlgebra
using Statistics
using SHA
using DataFrames
using CSV
using StableRNGs
using SparseArrays
using StaticArrays
using PlotlyJS
using Colors
using DelimitedFiles
using JLD2
using Preferences

using .SingleCellProjections.SCTransform
import .SCTransform: scparams

export
	set_samples_path,
	set_cache_path,
	get_samples_path,
	get_annotations_path,
	get_cache_path,
	figure2a,
	figure3b,
	supplementary_figure6,
	supplementary_figure7,
	supplementary_figure8,
	supplementary_figure11b


include("utils.jl")
include("colormaps.jl")
include("hex.jl")
include("grids.jl")
include("processdata.jl")
include("plotting.jl")
include("main.jl")
include("serialization.jl")

include("main_AML.jl")

end

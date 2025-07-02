module SingleCellProjections

export
	DataMatrix,
	ProjectionModel,
	FilterModel,
	LogTransformModel,
	TFIDFTransformModel,
	SCTransformModel,
	NormalizationModel,
	SVDModel,
	NearestNeighborModel,
	ObsAnnotationModel,
	VarCountsFractionModel,
	project,
	var_coordinates,
	obs_coordinates,
	load10x,
	loadh5ad,
	load_counts,
	merge_counts,
	filter_matrix,
	filter_var,
	filter_obs,
	covariate,
	designmatrix,
	normalization_model,
	normalize_matrix,
	svd,
	pma,
	update_matrix,
	logtransform,
	tf_idf_transform,
	knn_adjacency_matrix,
	force_layout,
	var_to_obs!,
	var_to_obs,
	var_to_obs_table,
	var_counts_fraction!,
	differentialexpression

using LinearAlgebra
import LinearAlgebra: svd

using SparseArrays
using ThreadedSparseArrays
using Statistics

using HDF5, H5Zblosc
using DataFrames # TODO: get rid of this dependency?

using Missings

using Random

# TODO: Make dependency optional
using PrincipalMomentAnalysis
import PrincipalMomentAnalysis: PMA, pma, simplices2kernelmatrixroot

using NearestNeighbors
using StaticArrays

using Distributions



include("MatrixExpressions/MatrixExpressions.jl")
using .MatrixExpressions


include("utils.jl")
include("table_utils.jl")

# SCTransform
include("SCTransform/SCTransform.jl")
import .SCTransform: scparams, sctransform

include("bilinear.jl")
include("sctransformsparse.jl")

# 10x FileIO
include("SingleCell10x/SingleCell10x.jl")
using .SingleCell10x


include("implicitsvd.jl")
include("implicitpma.jl")


include("barnes_hut.jl")
include("force_layout.jl")
include("embed.jl")

include("h5ad.jl")

include("lowrank.jl")
include("projectionmodels.jl")
include("datamatrix.jl")
include("subset_expression.jl")
include("filter.jl")
include("load.jl")
include("transform.jl")
include("normalize.jl")
include("reduce.jl")
include("annotate.jl")
include("counts_fraction.jl")

include("differentialexpression.jl")

include("precompile.jl")

end

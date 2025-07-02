# Utility function to avoid excessive serialization
discard_old_models!(data::DataMatrix) = deleteat!(data.models, 1:lastindex(data.models)-1)


lazyloaddata(samplePaths) = load_counts(load10x, string.(samplePaths); lazy=true, lazy_merge=true)
LazyCache.version(::typeof(lazyloaddata)) = v"0.1.0"

loaddata(lazy_merged) = (@info "Loading matrices"; @time load_counts(lazy_merged))
LazyCache.version(::typeof(loaddata)) = v"0.1.0"


loadannot(annotationPaths) = [DataFrame(CSV.File(p; delim=getdelim(p), stringtype=String, silencewarnings=true), copycols=true) for p in string.(annotationPaths)]
LazyCache.version(::typeof(loadannot)) = v"0.1.1"



function extendcellannotations(data::DataMatrix,annots;
                               on = nothing)
	cells = copy(data.obs)
	for df in annots
		df = df isa DataFrame ? df : DataFrame(df) # to support NamedTuples in annots
		leftjoin!(cells, df; on=something(on, hasproperty(df, :barcode) ? [:barcode,:sampleName] : [:sampleName]))
	end
	DataMatrix(data.matrix, data.var, cells)
end
function extendcellannotations(data::DataMatrix, source::DataMatrix; on=[:id,:barcode,:sampleName])
	cells = copy(data.obs)
	leftjoin!(cells, source.obs; on)
	DataMatrix(data.matrix, data.var, cells)
end
LazyCache.version(::typeof(extendcellannotations)) = v"0.2.0"


update_obs(data::DataMatrix, obs::DataFrame) = update_matrix(data, data.matrix; var=:keep, obs)
update_obs(data::DataMatrix, obs) = update_obs(data, DataFrame(obs)) # intended for NamedTuple (or other Tables) args
LazyCache.version(::typeof(update_obs)) = v"0.1.0"



function compute_scparams(countData, feature_type)
	@info "Computing SCTransform parameters";
	@time model = SCTransformModel(countData; var_filter = :feature_type => isequal(feature_type))
	model
end
LazyCache.version(::typeof(compute_scparams)) = v"0.4.0"


function compute_sctransform(countData, sctModel; kwargs...)
	@info "Computing SCTransform"
	@time transformed = project(countData, sctModel; kwargs...)
	empty!(transformed.models)
	transformed
end
LazyCache.version(::typeof(compute_sctransform)) = v"0.4.0"



function create_designmatrix(data, args...; kwargs...)
	@info "Creating Design Matrix"
	c = (x->covariate(x...)).(args)
	@time designmatrix(data, c...; kwargs...)
end
function create_designmatrix(data, model::NormalizationModel)
	@info "Creating Design Matrix"
	@time designmatrix(data, model.covariates)
end
LazyCache.version(::typeof(create_designmatrix)) = v"0.1.0"

function compute_normalization_model(args...; kwargs...)
	@info "Creating Normalization Model"
	@time NormalizationModel(args...; kwargs...)
end
LazyCache.version(::typeof(compute_normalization_model)) = v"0.1.0"

function apply_normalization(args...; kwargs...)
	@info "Applying Normalization"
	@time normalized = project(args...; kwargs...)
	empty!(normalized.models) # avoid excessive serialization
	normalized
end
LazyCache.version(::typeof(apply_normalization)) = v"0.1.0"


function sign_stabilize!(F::SVD)
	# make it so that the largest value (in absolute terms) is always positive
	ext = vec(extrema(F.U,dims=1))
	flip = last.(ext) .< .-first.(ext)
	F.U[:,flip] .*= -1
	F.Vt[flip,:] .*= -1
	F
end

function compute_svd(data; nsv, subspacedims=8nsv, niter=2, seed)
	@info "Computing SVD"
	rng = StableRNG(seed)
	@time reduced = svd(data; nsv, subspacedims, niter, rng)
	sign_stabilize!(reduced.matrix) # together with seed this should make it reproducible
	discard_old_models!(reduced) # avoid excessive serialization
	reduced
end
LazyCache.version(::typeof(compute_svd)) = v"0.5.1"


function compute_adjacency_matrix(data; k, kwargs...)
	@info "Computing adjacency matrix (k=$k)"
	@time knn_adjacency_matrix(obs_coordinates(data); k, kwargs...)
end
LazyCache.version(::typeof(compute_adjacency_matrix)) = v"0.1.0"

function compute_force_layout(reduced, adj; seed, kwargs...)
	@info "Computing force layout"
	rng = StableRNG(seed)
	@time fl = force_layout(reduced; adj, rng, obs=:keep, kwargs...)
	empty!(fl.models) # avoid excessive serialization
	fl
end
LazyCache.version(::typeof(compute_force_layout)) = v"0.2.0"


function compute_svd_project(data, base::DataMatrix{<:SVD}, args...; kwargs...)
	@info "Projecting sample points (SVD)"
	@time reduced = project(data,base,args...;kwargs...)
	empty!(reduced.models) # avoid excessive serialization
	reduced
end
LazyCache.version(::typeof(compute_svd_project)) = v"0.6.0"

function compute_nn_project(data, base, args...; kwargs...)
	@info "Computing nearest neighbor projection"
	adj = Ref(spzeros(Bool,0,0))
	@time projected = project(data, base, args...; adj_out = adj, kwargs...)
	empty!(projected.models) # avoid excessive serialization
	(;adj, data=projected)
end
LazyCache.version(::typeof(compute_nn_project)) = v"0.2.2"




function create_grid(coords, subdivisions)
	bounds = extrema(coords; dims=2)
	HexagonGrid(bounds, subdivisions)
end
create_grid(data::DataMatrix, args...; kwargs...) = create_grid(obs_coordinates(data), args...; kwargs...)
LazyCache.version(::typeof(create_grid)) = v"1.0.2"

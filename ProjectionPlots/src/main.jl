rotmatrix2d(α) = [cos(α) -sin(α); sin(α) cos(α)]
function maxbounds(points,α)
	ext = extrema(rotmatrix2d(α)*points; dims=2)
	maximum(last.(ext).-first.(ext))
end

# find the rotation of a data set that leads to best usage of (square) plotting area.
function optimalrotation(coords; niter=100)
	if size(coords,1)==2 # only support 2d for now
		val,α = minimum(α->(maxbounds(coords,α),α), range(0,stop=pi/2,length=niter))
		rotmatrix2d(α)
	else
		R = I(size(coords,1))
	end
end
optimalrotation(data::DataMatrix; kwargs...) = optimalrotation(obs_coordinates(data); kwargs...)

apply_rotation(R, data::DataMatrix) = update_matrix(data, R*data.matrix; var=:keep, obs=:keep)
apply_rotation(R, nt::NamedTuple) = apply_rotation(R, nt.data)
LazyCache.version(::typeof(apply_rotation)) = v"0.1.0"


maketitle(p,d) = (get(p,:titlefun,(p,d)->p.plotTitle)(p,d))


function get_ncells_string(d)
	rotEmbedding = haskey(d,:rotSampleEmbedding) ? d.rotSampleEmbedding : d.rotForceLayout
	coords = obs_coordinates(rotEmbedding[])

	cellCount = string(size(coords,2))
	cellCount = reverse(join(map(i->cellCount[i], Iterators.partition(length(cellCount):-1:1,3)), " ")) # add thousands separators to string!
end


function load_sample_data(;samplePath, additionalCellAnnot, sampleIDs, cellfilterfun, update_obs_fun=nothing)
	sampleIDs = sort(sampleIDs) # for increased caching stability
	basePaths = checksummedpath.(joinpath.(samplePath,string.(sampleIDs,".h5")))

	lazyCountData = lazy(lazyloaddata, basePaths)
	lazyCountDataAnn = lazy(extendcellannotations, lazyCountData, additionalCellAnnot)
	if update_obs_fun !== nothing
		obs = NamedTuple(pairs(eachcol(update_obs_fun(lazyCountDataAnn[].obs))))
		lazyCountDataAnn = lazy(update_obs, lazyCountDataAnn, obs)
	end

	countData = lazy(loaddata, lazyCountData)
	countDataAnn = lazy(extendcellannotations, countData, lazyCountDataAnn)

	# NB: By computing cellInd directly (not lazily), we get caching by value in the filtering step, which means
	# that we avoid expensive recomputations if annotation files have changed, but the mask is actually the same!
	cellInd = findall(cellfilterfun, eachrow(lazyCountDataAnn[].obs))

	lazyCountData, lazyCountDataAnn, countData, countDataAnn, cellInd
end



function transform_counts(;lazyCountDataAnn, countData, cellInd, featureType, elim_args)
	sctModel = cached(compute_scparams, countData, featureType)

	transformed = lazy(compute_sctransform, countData, sctModel; post_obs_filter=cellInd)
	transformedAnn = lazy(extendcellannotations, transformed, lazyCountDataAnn) # could be moved out from this function


	# Extract only those annotations needed for the design matrix - and cache them by value.
	# This prevents unnecessary invalidations when annotations change, but the design matrix does not.
	annForElim = lazyCountDataAnn[].obs[:, vcat("barcode","sampleName",first.(elim_args)...)] # keep only columns actually used for eliminated factors
	annForElim = NamedTuple(pairs(eachcol(annForElim))) # convert to NamedTuple that can be checksummed
	transformedAnnForElim = lazy(extendcellannotations, transformed, (annForElim,))
	dm = lazy(create_designmatrix, transformedAnnForElim, elim_args...)




	# ...but not to normalize
	normalizationModel = cached(compute_normalization_model, transformed, dm)
	normalized = lazy(apply_normalization, transformed, normalizationModel, dm)

	sctModel, transformed, transformedAnn, normalizationModel, normalized
end

function transform_counts2(;lazySampleCountDataAnn, sampleCountData, sctModel, sampleCellInd, normalizationModel, elim_args)
	transformed2 = lazy(compute_sctransform, sampleCountData, sctModel; post_obs_filter=sampleCellInd)

	transformedAnn2 = lazy(extendcellannotations, transformed2, lazySampleCountDataAnn) # could be moved out from this function

	# Extract only those annotations needed for the design matrix - and cache them by value.
	# This prevents unnecessary invalidations when annotations change, but the design matrix does not.
	annForElim2 = lazySampleCountDataAnn[].obs[:, vcat("barcode","sampleName",first.(elim_args)...)] # keep only columns actually used for eliminated factors
	annForElim2 = NamedTuple(pairs(eachcol(annForElim2))) # convert to NamedTuple that can be checksummed
	transformedAnnForElim2 = lazy(extendcellannotations, transformed2, (annForElim2,))
	dm2 = lazy(create_designmatrix, transformedAnnForElim2, normalizationModel)



	normalized2 = lazy(apply_normalization, transformed2, normalizationModel, dm2) # doesn't need annotations

	transformed2, transformedAnn2, normalized2
end


function force(;reduced, force_args)
	forceAdj = cached(compute_adjacency_matrix, reduced; k=force_args.k)
	force_args = Base.structdiff(force_args, NamedTuple{(:k,)}) # This is the "recommended" way until structdiff is exported from Julia base
	forceLayout = cached(compute_force_layout, reduced, forceAdj; force_args...)
	forceAdj, forceLayout
end
function nn_projection_model(pre, post; k)
	lazy(NearestNeighborModel, "force_layout", pre, post;
	                           k, var="Force Layout Dim ", obs=:keep)
end




function add_cellcount_overlay!(traces, layout, p, d; color="black")
	rotEmbedding = haskey(d,:rotSampleEmbedding) ? d.rotSampleEmbedding : d.rotForceLayout
	coords = obs_coordinates(rotEmbedding[])

	cellCount = get_ncells_string(d)

	ext = extrema(obs_coordinates(d.rotForceLayout[]); dims=2)
	pos = first.(ext).*[0.88, 0.95] .+ last.(ext).*[0.12, 0.05]

	add_text_overlay!(traces, layout, pos, string("Cell count=",cellCount); fontSize=p.overlayFontSize, color)
end

_grid_x_text_offset(g::HexagonGrid) = sqrt(3√3/2)*g.side

function add_heatmap_size_overlay!(traces, layout, p, d; color="black", cutoffs, radii, max_cells, pixel_color=RGB(0.5,0.5,0.5), pos=[0.8, 0.99])
	@assert length(cutoffs)==length(radii)

	# Get rid of unused cutoffs (to avoid displaying when upper bound being below lower bound)
	max_ind = findlast(<=(max_cells), cutoffs)
	if max_ind !== nothing
		cutoffs = cutoffs[1:max_ind]
		radii = radii[1:max_ind]
	end


	baseCoords = obs_coordinates(d.rotForceLayout[])
	@assert size(baseCoords,1)==2
	ext = extrema(baseCoords; dims=2)
	abspos = first.(ext).*(1.0.-pos) .+ last.(ext).*pos
	y_offset = -(last(ext[2])-first(ext[2]))*0.04 * get(p,:nCellsYSep,1.0)

	grid = d.grid[]
	shapes = [plotly_shape(grid, abspos[1], abspos[2]+y_offset*(i-1), r; fillcolor=pixel_color, opacity=1.0, line_width=0) for (i,r) in enumerate(radii)]


	_relayout!(layout; shapes=vcat(layout.shapes, shapes))

	for (i,c) in enumerate(cutoffs)
		if i < length(cutoffs)
			c2 = cutoffs[i+1]-1
		elseif max_cells !== nothing
			c2 = max_cells
		else
			c2 = "?"
		end
		t = string(c, '–', c2, " cells") # NB: en-dash
		add_text_overlay!(traces, layout, abspos+[_grid_x_text_offset(grid), y_offset*(i-1)], t; fontSize=p.overlayFontSize, color, xanchor="left")
	end

	add_text_overlay!(traces, layout, abspos + [-_grid_x_text_offset(grid)/2,-y_offset], "Number of cells"; fontSize=p.overlayFontSize, color, xanchor="left")

	traces, layout
end




function add_overlays!(traces, layout, p, d; celltype=p.overlay_celltype, color="black", embedding=:force)
	if embedding==:force
		rotEmbedding = d.rotForceLayoutAnn[]
		coords = obs_coordinates(rotEmbedding)
		cellAnnot = rotEmbedding.obs

		if celltype !== nothing
			traces,layout = add_categorical_overlay!(traces, layout, coords; text=cellAnnot[!,celltype], color, p.category_overlay_args...)
		end
	end
	traces, layout
end

function plot_categorical(p,d; annotationName::String, colorDict, outName="categorical_$annotationName", embedding=:force, legendTitle=annotationName, order=nothing, bgAnnotationName=annotationName)
	if embedding == :force
		rotEmbedding = haskey(d,:rotSampleEmbedding) ? d.rotSampleEmbeddingAnn[] : d.rotForceLayoutAnn[]
		baseEmbedding = d.rotForceLayout[]
	elseif embedding == :pca
		@assert !haskey(d,:rotSampleEmbedding) "PCA projection not yet supported." # TEMP workaround
		baseEmbedding = d.reduced[]
		rotEmbedding = d.reducedAnn[]
	end
	cellAnnot = rotEmbedding.obs

	colorBy = coalesce.(cellAnnot[!,annotationName], "NA")

	baseCoords = obs_coordinates(baseEmbedding)[1:p.ndim,:]
	sampleCoords = obs_coordinates(rotEmbedding)[1:p.ndim,:]


	traces,layout = plotscattercategorical(sampleCoords; colorby=colorBy, colorDict, title=maketitle(p,d), p.axisPrefix, p.markerSize, legendTitle, p.fontSize, order, p.axis_kwargs)

	baseEmbedding.matrix != rotEmbedding.matrix && prepend!(traces, bg_scatter(baseCoords; p.bgMarkerSize, p.bgColor, p.bgOpacity))
	(add_overlays!(traces, layout, p, d; embedding, celltype=bgAnnotationName)..., outName)
end
plot_categorical(annotationName;kwargs...) = (p,d)->plot_categorical(p,d;annotationName,kwargs...)




function plot_numerical_heatmap(p,d; annotationName, f=mean, logtransform=false, logoffset=1, show_size=true, legendTitle=annotationName, cmin=nothing, cmax=nothing)
	rotEmbedding = haskey(d,:rotSampleEmbedding) ? d.rotSampleEmbeddingAnn : d.rotForceLayoutAnn

	baseCoords = obs_coordinates(d.rotForceLayout[])
	sampleCoords = obs_coordinates(rotEmbedding[])
	colorby = rotEmbedding[].obs[:,annotationName]

	extra = (;)
	if show_size
		cutoffs, radii, size_by, size_fun = default_heatmap_sizing(d.grid[], size(sampleCoords,2); cellcountmin=p.heatmap_mincells)
		extra = (;size_by, size_fun)
	end

	heatmapData, pointInds = setupheatmapdata(d.grid[], sampleCoords, colorby, f; cellcountmin=p.heatmap_mincells, logtransform, logoffset, extra...)
	traces,layout = plotheatmap(d.grid[], baseCoords, sampleCoords, heatmapData, pointInds; p.bgMarkerSize, p.bgColor, p.bgOpacity, p.markerSize, title=maketitle(p,d), p.axisPrefix, colorScale=p.color_scales.numerical, cmin, cmax, colorbarTitle=legendTitle, p.fontSize, p.axis_kwargs)
	if show_size
		traces,layout = add_heatmap_size_overlay!(traces,layout,p,d; cutoffs, radii, max_cells=maximum(heatmapData.NbrPoints), p.heatmap_size_overlay_args...)
	end

	(add_overlays!(traces, layout, p, d)..., "heatmap_numerical_$annotationName")
end
plot_numerical_heatmap(annotationName; kwargs...) = (p,d)->plot_numerical_heatmap(p, d; annotationName, kwargs...)



function plot_cellcount(p,d; normalize=true, logtransform=true, cmin=nothing, cmax=nothing, cellcountmin=0, outName="cellcount")
	baseCoords = obs_coordinates(d.rotForceLayout[])
	sampleCoords = obs_coordinates(d.rotSampleEmbedding[])
	N = size(sampleCoords,2)

	if normalize
		cellWeights = fill(1e4/N, N)
		colorbarTitle = "Normalized<br>cell count"
	else
		cellWeights = fill(1.0, N)
		colorbarTitle = "Cell count"
	end
	logtransform && (colorbarTitle = string("log₁₀<br>",colorbarTitle))

	heatmapData, pointInds = setupheatmapdata(d.grid[], sampleCoords, cellWeights, sum; logtransform, cellcountmin)

	traces,layout = plotheatmap(d.grid[], baseCoords, sampleCoords, heatmapData, pointInds; p.bgMarkerSize, p.bgColor, p.bgOpacity, p.markerSize, title=maketitle(p,d), p.axisPrefix, colorScale=p.color_scales.cell_counts, colorbarTitle, cmin, cmax, p.fontSize, p.axis_kwargs)
	(add_overlays!(traces, layout, p, d)..., outName)
end
plot_cellcount(; kwargs...) = (p,d)->plot_cellcount(p,d; kwargs...)



function het2mutationfreq(mutData; min_reads=0)
	mut = first.(mutData)
	wt = last.(mutData)

	mut = coalesce.(mut, 0)
	wt = coalesce.(wt, 0)


	mutTot = sum(mut)
	wtTot = sum(wt)

	mutTot+wtTot < min_reads && return missing

	# Model 1 - based on total number of reads
	f = mutTot / max(1, mutTot+wtTot) # put to zero if no data - use missing values instead???
	mutFreq = min(2f,1)

	# Model 2 - based on number of cells with at least one mutated transcript detected
	mutFreq2 = count(!iszero, mut) / length(mut) # frequency of cells with at least one mutated transcript - more conservative

	# take the max of the two models and convert to %
	100*max(mutFreq,mutFreq2)
end
het2mutationfreq(;kwargs...) = mutData->het2mutationfreq(mutData; kwargs...)


function plot_mutation_heatmap(p,d; show_size=true)
	baseCoords = obs_coordinates(d.rotForceLayout[])
	sampleCoords = obs_coordinates(d.rotSampleEmbedding[])
	cellAnnot = d.rotSampleEmbeddingAnn[].obs

	colorby = tuple.(cellAnnot[!,"heterozygous.mutSum"], cellAnnot[!,"heterozygous.wtSum"])

	extra = (;)
	if show_size
		cutoffs, radii, size_by, size_fun = default_heatmap_sizing(d.grid[], size(sampleCoords,2); cellcountmin=p.heatmap_mincells)
		extra = (;size_by, size_fun)
	end

	# TODO: Make parameter
	missing_color = colorant"#A08968"


	heatmapData, pointInds = setupheatmapdata(d.grid[], sampleCoords, colorby, het2mutationfreq(;min_reads=3); cellcountmin=p.heatmap_mincells, extra...)
	traces,layout = plotheatmap(d.grid[], baseCoords, sampleCoords, heatmapData, pointInds; cmin=0.0, cmax=100.0, p.bgMarkerSize, p.bgColor, p.bgOpacity, p.markerSize, title=maketitle(p,d), p.axisPrefix, colorScale=p.color_scales.mutations, p.fontSize, missing_color, p.axis_kwargs)

	if show_size
		traces,layout = add_heatmap_size_overlay!(traces,layout,p,d; cutoffs, radii, max_cells=maximum(heatmapData.NbrPoints), p.heatmap_size_overlay_args...)
	end
	(add_overlays!(traces, layout, p, d)..., "mutation_frequency")
end





function setupmain(p)
	additionalCellAnnot = lazy(loadannot, checksummedpath.(p.annotationPaths))
	lazyCountData, lazyCountDataAnn, countData, countDataAnn, cellInd = load_sample_data(;p.samplePath, additionalCellAnnot, sampleIDs=p.baseIDs, p.cellfilterfun, update_obs_fun=get(p,:update_obs_fun,nothing))
	sctModel, transformed, transformedAnn, normalizationModel, normalized = transform_counts(;lazyCountDataAnn, countData, cellInd, p.featureType, p.elim_args)
	reduced = cached(compute_svd, normalized; p.pca_args...)
	forceAdj, forceLayout = force(;reduced, p.force_args)

	R = lazy(optimalrotation, forceLayout)
	p.rotation !== nothing && (R = lazy(*, p.rotation, R)) # Custom rotation (if specified)
	rotForceLayout = lazy(apply_rotation, R, forceLayout)

	# Annotated versions of the data matrices
	normalizedAnn = lazy(extendcellannotations, normalized, lazyCountDataAnn)
	reducedAnn = lazy(extendcellannotations, reduced, lazyCountDataAnn)
	rotForceLayoutAnn = lazy(extendcellannotations, rotForceLayout, lazyCountDataAnn)

	grid = lazy(create_grid, rotForceLayout, p.subdivisions)

	(;additionalCellAnnot, lazyCountData, lazyCountDataAnn, countData, countDataAnn, cellInd, sctModel, transformed, transformedAnn, normalizationModel, normalized, normalizedAnn, reduced, reducedAnn, forceAdj, forceLayout, R, rotForceLayout, rotForceLayoutAnn, grid)
end


function setupembedding(p, d)
	lazySampleCountData, lazySampleCountDataAnn, sampleCountData, sampleCountDataAnn, sampleCellInd = load_sample_data(;p.samplePath, d.additionalCellAnnot, p.sampleIDs, cellfilterfun=p.cellfilterfun_sample, update_obs_fun=get(p,:update_obs_fun,nothing))

	transformed2, transformedAnn2, normalized2 = transform_counts2(;lazySampleCountDataAnn, sampleCountData, d.sctModel, sampleCellInd, d.normalizationModel, p.elim_args)
	
	reduced2 = cached(compute_svd_project, normalized2, d.reduced)
	nnProjectionModel = nn_projection_model(d.reduced, d.forceLayout; p.embed_args...)
	sampleEmbedding = cached(compute_nn_project, reduced2, nnProjectionModel)


	rotSampleEmbedding = lazy(apply_rotation, d.R, sampleEmbedding)

	# Annotated versions of the data matrices
	normalizedAnn2 = lazy(extendcellannotations, normalized2, lazySampleCountDataAnn)
	reducedAnn2 = lazy(extendcellannotations, reduced2, lazySampleCountDataAnn)
	rotSampleEmbeddingAnn = lazy(extendcellannotations, rotSampleEmbedding, lazySampleCountDataAnn)

	(;lazySampleCountData, lazySampleCountDataAnn, sampleCountData, sampleCountDataAnn, sampleCellInd, transformed2, transformedAnn2, normalized2, normalizedAnn2, reduced2, reducedAnn2, sampleEmbedding, rotSampleEmbedding, rotSampleEmbeddingAnn)
end





let cache_initialized = false
	global function initcache()
		cache_initialized && return
		LazyCache.init(expanduser(get_cache_path()))
		cache_initialized = true
	end
end

function main(p; setup=setupmain, setupproj=setupembedding)
	initcache()
	outPath = p.outPath
	fnPrefix = ""
	d = setup(p)

	if p.plotmode == :groups
		d2 = setupproj(p,d)
		d = (;d..., d2...)
		pt = replace(p.plotTitle,":"=>"") # ':' can appear in fusion gene names, remove ':' for file system compatibility
		fnPrefix = string(pt,'_')
	end

	if p.savePlots
		!isdir(outPath) && mkpath(outPath)
	end

	for plotfun in p.plotfuns
		x = plotfun(p,d)
		if x !== nothing
			traces,layout,plotName = x
			pl = Plot(traces,layout)
			p.displayPlots && display(plot(pl))
			if p.savePlots
				for fmt in p.plotFormats
					fn = joinpath(outPath,string(fnPrefix, plotName, '.', fmt))
					@info "Saving figure \"$fn\""
					@time savefig(pl, fn; width=p.plotWidth, height=p.plotHeight, scale=get(p,:plotScale,1))
				end
			end
		end
	end
end

function main_samples(p; setup=setupmain, setupproj=setupembedding)
	@assert length(p.sampleIDs)==length(p.plotTitle)
	initcache()
	d = setup(p)

	for (id,sampleName) in zip(p.sampleIDs,p.plotTitle)
		pSample = (; p..., plotTitle=sampleName, sampleIDs=[id])
		d2 = setupproj(pSample,d)

		sampleFilename = replace(sampleName,":"=>"") # ':' can appear in fusion gene names, remove ':' for file system compatibility

		outPath = p.outPath

		dSample = (;d..., d2...)

		if p.savePlots
			!isdir(outPath) && mkpath(outPath)
		end

		for plotfun in p.plotfuns
			x = plotfun(pSample,dSample)
			if x !== nothing
				traces,layout,plotName = x
				pl = Plot(traces,layout)
				p.displayPlots && display(plot(pl))
				if p.savePlots
					for fmt in p.plotFormats
						fn = joinpath(outPath,string(sampleFilename,'_',plotName, '.', fmt))
						@info "Saving figure \"$fn\""
						@time savefig(pl, fn; width=p.plotWidth, height=p.plotHeight, scale=get(p,:plotScale,1))
					end
				end
			end
		end
	end
end

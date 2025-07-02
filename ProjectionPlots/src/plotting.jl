# Workaround for PlotlyBase.jl where relayout! destroys the margin attribute.
_getmargin(sp::PlotlyJS.SyncPlot) = _getmargin(sp.plot)
_getmargin(p::Plot) = _getmargin(p.layout)
_getmargin(l::Layout) = l.margin
_relayout!(x, update; kwargs...) = relayout!(x, update; margin=get(update, :margin, _getmargin(x)), kwargs...)
_relayout!(x; kwargs...) = relayout!(x; margin=_getmargin(x), kwargs...)


function createlayout(ndim; title="", axisPrefix="dim", legendTitle=nothing, colorScale=nothing, colorbarTitle=nothing, cmin=nothing, cmax=nothing, fontSize=nothing, axis_kwargs=(;), kwargs...)
	margin = attr(l=0,t=5,r=0,b=65,pad=0)

	args1 = Dict{Symbol,Any}(pairs((;
	                                 title = attr(text=title,x=0.07,y=0.03), # Updated - Title at the bottom with space for two lines (lower smaller font size)
	                                 plot_bgcolor = colorant"white",
	                                 paper_bgcolor = colorant"white",
	                                 hovermode = "closest",
	                                 margin,
	                                 kwargs...)))

	if isempty(axisPrefix)
		args1[:xaxis] = attr(;mirror=true, showgrid=false, showticklabels=false, zeroline=false, constrain="domain", axis_kwargs...)
		args1[:yaxis] = attr(;mirror=true, showgrid=false, showticklabels=false, zeroline=false, constrain="domain", scaleanchor='x', axis_kwargs...) # scaleanchor used to set aspect ratio to 1
	else
		args1[:xaxis] = attr(;title="$(axisPrefix) 1", showgrid=false, showticklabels=false, zeroline=false, constrain="domain", axis_kwargs...)
		args1[:yaxis] = attr(;title="$(axisPrefix) 2", showgrid=false, showticklabels=false, zeroline=false, constrain="domain", scaleanchor='x', axis_kwargs...) # scaleanchor used to set aspect ratio to 1
	end

	isnothing(legendTitle) || (args1[:legend] = attr(title_text=legendTitle, itemsizing="constant"))

	args2 = Dict()
	cs = isnothing(colorScale) ? (;) : (;colorscale=colorScale)
	ct = isnothing(colorbarTitle) ? (;) : (;colorbar=attr(title=attr(text=colorbarTitle)))
	cmin = isnothing(cmin) ? (;) : (; cmin)
	cmax = isnothing(cmax) ? (;) : (; cmax)

	caxis = attr(; cs..., ct..., cmin..., cmax...)
	isempty(caxis) || (args2[:coloraxis] = caxis)

	isnothing(fontSize) || (args2[:font] = attr(size=fontSize))

	Layout(;args1..., args2...)
end


_scatter(V; kwargs...) = size(V,2)==2 ? scatter(;x=V[:,1], y=V[:,2], kwargs...) : scatter3d(;x=V[:,1], y=V[:,2], z=V[:,3], kwargs...)



project_shapes(args...) = false,nothing
function project_shapes(::Val{:sphere}, c1, c2, radius)
	d = c2-c1
	dist = norm(d)
	dist >= radius*0.99 && return false, d # early out
	d .*= (radius-dist) ./ dist
	true, d
end
function project_shapes(::Val{:box}, c1, c2, extents)
	d = c2-c1
	s = sign.(d)
	d .= abs.(d)
	any(d .>= 0.99.*extents) && return false, d

	d .= (extents .- d)/2
	dim = argmin(d)
	# d[1:dim-1]   .= 0
	# d[dim+1:end] .= 0
	d[1:dim-1]   .= 0.5*d[dim] # move slightly in the
	d[dim+1:end] .= 0.5*d[dim] # other directions to!
	d .*= s
	true, d
end

# Treat each text center as a sphere/box and make sure they do not overlap
function layout_categorical(centers, shape, extents; maxiter=1000)
	centers = copy(centers)
	N = size(centers,2)
	noverlap = zeros(Int,N)
	v = zeros(size(centers))

	for iter=1:maxiter
		noverlap .= 0
		v .= 0
		for i=1:N-1
			for j=i+1:N
				overlap,d = project_shapes(shape, centers[:,i], centers[:,j], extents)
				if overlap
					# @show d
					noverlap[i] += 1
					noverlap[j] += 1
					v[:,i] .-= d
					v[:,j] .+= d
				end
			end
		end
		all(iszero, noverlap) && break
		# @show noverlap
		noverlap .= max.(1,noverlap) # avoid div by zero
		v ./= noverlap'
		centers .+= v
	end
	centers
end


function categorical_setup(points; text, minCells, shape, extents, offsets)
	df = DataFrame(;text=string.(coalesce.(text,"NA")))
	centers = Float64[]
	texts = String[]
	for sub in groupby(df,:text)
		tb = sub[1,:text]
		ismissing(tb) && continue
		if size(sub,1)<minCells
			@info "Removed \"$tb\" overlay text due to too few cells ($(size(sub,1))<$minCells)"
			continue
		end
		ind = parentindices(sub)[1]
		append!(centers, vec(mean(points[:,ind],dims=2)))
		push!(texts, tb)
	end
	centers = reshape(centers, size(points,1), :)

	if offsets !== nothing
		centers2 = copy(centers)
		for (j,t) in enumerate(texts)
			centers2[:,j] += get(offsets, t, zeros(size(points,1)))
		end

	elseif shape !== nothing && extents !== nothing
		centers2 = layout_categorical(centers, Val(shape), extents)
	else
		centers2 = copy(centers)
	end


	texts,centers,centers2
end


function add_categorical_overlay!(traces, layout, points; text, fontSize, minCells=0, color="black", shape=:none, extents=nothing, offsets=nothing)
	texts,centers,centers2 = categorical_setup(points; text, minCells, shape, extents, offsets)

	annotations = [attr(x=centers[1,i], y=centers[2,i], ax=centers2[1,i], ay=centers2[2,i], axref='x', ayref='y', text=t, font=attr(;size=fontSize, color), showarrow=true) for (i,t) in enumerate(texts)]
	_relayout!(layout; annotations=vcat(get(layout.fields,:annotations,[]), annotations))

	traces, layout
end

function add_text_overlay!(traces, layout, pos, text; fontSize, color="black", kwargs...)
	annotations = [attr(; x=pos[1], y=pos[2], text, font=attr(;size=fontSize, color), showarrow=false, kwargs...)]
	_relayout!(layout; annotations=vcat(get(layout.fields,:annotations,[]), annotations))
	traces, layout
end



function _plotscattercategorical(points; colorby, colorDict, markerSize, opacity=1.0, showlegend=true, order=nothing)
	missingKeys = setdiff(unique(colorby), keys(colorDict))
	@assert isempty(missingKeys) "No color given for keys: $(join(missingKeys,','))"

	df = DataFrame(;colorby)
	traces = GenericTrace[]

	order === nothing && (order = unique(df[!,:colorby]))
	subs = groupby(df,:colorby)
	subIDs = [sub[1,:colorby] for sub in subs]
	@assert isempty(setdiff(subIDs, order)) "$(join(setdiff(subIDs, order),',')) missing from order"
	subOrder = indexin(order, subIDs)

	for subInd in subOrder
		subInd === nothing && continue
		sub = subs[subInd]
		cb = sub[1,:colorby]
		ind = parentindices(sub)[1]
		push!(traces, _scatter(points[:,ind]'; mode="markers", marker=attr(;color=colorDict[cb], size=markerSize, line_width=0), opacity, name=string(cb), showlegend))
	end
	traces
end

function plotscattercategorical(points; colorby, colorDict, markerSize, legendTitle, order=nothing, kwargs...)
	traces = _plotscattercategorical(points; colorby, colorDict, markerSize, showlegend=legendTitle!==nothing, order)
	traces, createlayout(size(points,1); legendTitle, kwargs...)
end




function bg_scatter(background_points; bgMarkerSize, bgColor, bgOpacity)
	colorby = fill("",size(background_points,2))
	colorDict = Dict(""=>bgColor)
	_plotscattercategorical(background_points; colorby, colorDict, markerSize=bgMarkerSize, opacity=bgOpacity, showlegend=false)
end


function default_heatmap_sizing(grid::HexagonGrid, N; cellcountmin=0)
	size_by = fill(nothing,N)
	cutoffs = ceil.(Int, 5.0.^(-1:1) * 3N * grid.volume_frac)
	cutoffs = max.(cellcountmin, cutoffs)

	n = length(cutoffs)
	radii = (1:n)/n

	k = findlast(==(first(cutoffs)), cutoffs)
	cutoffs = cutoffs[k:end]
	radii = radii[k:end]

	size_fun = v->(i=searchsortedlast(cutoffs, length(v)); i==0 ? 0.0 : radii[i])
	cutoffs, radii, size_by, size_fun
end



function heatmap_aggregate(::Val{ndim}, pointInds, colorby, f, logtransform=false, logoffset=1) where ndim
	# compute block aggregates
	df = DataFrame(pointInd=reinterpret(NTuple{ndim,Int},vec(pointInds)), colorby=colorby)
	df = combine(groupby(df,:pointInd), nrow=>:NbrPoints, :colorby=>f=>:aggregate)
	logtransform && (df.aggregate = log10.(df.aggregate.+logoffset))
	df
end
heatmap_aggregate(pointInds::AbstractArray, args...) = heatmap_aggregate(Val(size(pointInds,1)), pointInds, args...)


heatmap2centers(df, grid) = (reshape(reinterpret(Int,df.pointInd),length(grid.offset),:).-0.5) .* grid.side .+ grid.offset # convert cube indices to cube centers

function heatmap2pointvalues(::Val{ndim}, df, pointInds) where ndim
	dfPoints = DataFrame(pointInd=reinterpret(NTuple{ndim,Int},vec(pointInds)))
	leftjoin!(dfPoints, df; on=:pointInd).aggregate
end
heatmap2pointvalues(df, pointInds) = heatmap2pointvalues(Val(size(pointInds,1)), df, pointInds)


function polygon_svgpath(points)
	@assert size(points,1)==2
	io = IOBuffer()
	for i in 1:size(points,2)
		print(io, i==1 ? 'M' : 'L', ' ', points[1,i], ',', points[2,i], ' ')
	end
	print(io, "z")
	String(take!(io))
end


function hexagon_svgpath(x, y, side)
	dx1 = side*√3/2
	dy1 = side/2
	dy2 = side
	# NB: for some reason, relative coordinates do not work, use absolute coordinates
	"M $(x+dx1),$(y+dy1) L $(x+dx1),$(y-dy1) L $(x),$(y-dy2) L $(x-dx1),$(y-dy1) L $(x-dx1),$(y+dy1) L $(x),$(y+dy2) z"
end
plotly_rect(x, y, r; kwargs...) = rect(x-r, x+r, y-r, y+r; kwargs...)

function plotly_shape(g::HexagonGrid, x, y, r; kwargs...)
	path(hexagon_svgpath(x, y, r*g.side); kwargs...)
end



function heatmap2shapes(df, grid::HexagonGrid, cmin, cmax, colormap; missing_color, show_block_id=false, fontSize=8)
	@assert grid_ndims(grid)==2 # only supports 2d

	pointInd = reshape(reinterpret(Int,df.pointInd), length(grid.offset), :)
	p = grid_centers(grid, pointInd)

	# map color values to [0,1] interval
	c = df.aggregate

	if cmin === nothing
		cmin = all(ismissing,c) ? 0.0 : minimum(skipmissing(c))
	end
	if cmax === nothing
		cmax = all(ismissing,c) ? 1.0 : maximum(skipmissing(c))
	end

	c = (c.-cmin)./(cmax.-cmin)
	clamp!(c,0,1)

	c = coalesce.(lookup.(Ref(colormap),c), missing_color)

	r = hasproperty(df, :aggregate_size) ? df.aggregate_size : fill(1.0,size(df,1))

	shapes = [plotly_shape(grid, p[1,i], p[2,i], r[i]; fillcolor=color(c[i]), opacity=alpha(c[i]), line_width=0) for i in 1:size(p,2)]

	annotations = []
	if show_block_id
		annotations = [attr(; x=p[1,i], y=p[2,i], text="$(pointInd[1,i]),$(pointInd[2,i])", font=attr(;size=fontSize), showarrow=false) for i in 1:size(p,2)]
	end

	shapes, annotations
end


function setupheatmapdata(grid::HexagonGrid, points, colorby, f; size_by=nothing, size_fun=nothing, logtransform=false, logoffset=1, cellcountmin=0)
	pointInds = grid_indices(grid, points)

	heatmapData = heatmap_aggregate(pointInds, colorby, f, logtransform, logoffset)
	if size_by !== nothing
		heatmapData2 = heatmap_aggregate(pointInds, size_by, size_fun)
		heatmapData.aggregate_size = heatmapData2.aggregate

		filter!(row->row.aggregate_size>0, heatmapData) # skip blocks with size 0
	end

	nCellsTotal = sum(heatmapData.NbrPoints)

	filter!(row->row.NbrPoints>=cellcountmin, heatmapData) # skip blocks with too few cells

	nCellsFiltered = sum(heatmapData.NbrPoints)
	@info "$nCellsFiltered/$nCellsTotal cells remaining in heatmap after filtering"

	heatmapData, pointInds
end



function plotheatmap(grid, background_points, points, heatmapData, pointInds; bgMarkerSize, bgColor, bgOpacity, markerSize, title, colorbarTitle=nothing, cmin=nothing, cmax=nothing, colorScale, missing_color=RGB(1.0,0.0,1.0), kwargs...)

	traces = bg_scatter(background_points; bgMarkerSize, bgColor, bgOpacity)
	layout = createlayout(size(points,1); title, colorbarTitle, colorScale, cmin, cmax, kwargs...)

	# Enable this to filter out "pixels" with missing values
	# filter!(row->row.aggregate!==missing, heatmapData)

	if !all(ismissing,heatmapData.aggregate)
		if cmax !== nothing
			dataMax = maximum(skipmissing(heatmapData.aggregate))
			cmax<dataMax && @warn "Maximum color value: $dataMax clamped to $cmax for plot \"$title\" - \"$colorbarTitle\""
		end
		if cmin !== nothing
			dataMin = minimum(skipmissing(heatmapData.aggregate))
			cmin>dataMin && @warn "Minimum color value: $dataMin clamped to $cmin for plot \"$title\" - \"$colorbarTitle\""
		end
	end

	shapes,annotations = heatmap2shapes(heatmapData, grid, cmin, cmax, colorScale; missing_color)

	# To get colorbar to show
	if colorbarTitle !== nothing
		!all(ismissing,heatmapData.aggregate) && @show extrema(skipmissing(heatmapData.aggregate))
		c1 = cmin===nothing ? minimum(heatmapData.aggregate) : cmin
		c2 = cmax===nothing ? maximum(heatmapData.aggregate) : cmax
		h = _scatter(zeros(2,2); mode="markers", opacity=0.0, marker=attr(;color=[c1,c2], coloraxis="coloraxis", name=""), showlegend=false)
		push!(traces,h)
	end

	_relayout!(layout; shapes, annotations=vcat(get(layout.fields,:annotations,[]), annotations))

	traces, layout
end

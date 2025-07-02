struct HexagonGrid
	side::Float64
	offset::Vector{Float64}
	volume_frac::Float64
end
function HexagonGrid(bounds, subdivisions)
	@assert size(bounds,1) == 2
	bounding_side = maximum(last.(bounds).-first.(bounds))
	side = bounding_side / (subdivisions*sqrt(3√3/2)) # NB: normalization factor to make hexagons have the same area as squares above
	offset = vec((last.(bounds).+first.(bounds))./2)
	bounding_volume = bounding_side^2
	volume = 3√3/2 * side^2
	HexagonGrid(side, offset, volume/bounding_volume)
end

grid_ndims(g::HexagonGrid) = length(g.offset)
grid_indices(g::HexagonGrid, points) = roundhex(cartesian2hex((points.-g.offset)./g.side))
grid_centers(g::HexagonGrid, indices) = hex2cartesian(indices).*g.side .+ g.offset

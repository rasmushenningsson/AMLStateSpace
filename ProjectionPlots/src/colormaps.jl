function scientificcolourmap(name::String; rev::Bool=false)
	p = joinpath("ScientificColourMaps6", "$name.txt")
	@assert isfile(p) "Couldn't find colormap \"$name\""
	A = readdlm(p)
	rev && (A=reverse(A;dims=1))
	x = range(0.0;stop=1.0,length=size(A,1))
	[(t,RGB{Float64}(r...)) for (t,r) in zip(x,eachrow(A))]
end

function lookup(colormap::Vector{Tuple{Float64,T}}, t::Real) where {T<:Colorant}
	i2 = searchsortedfirst(colormap, t; by=first)
	i2>length(colormap) && return last(colormap[end])
	i2==1 && return last(colormap[1])
	i1 = i2-1
	α = (t-first(colormap[i1])) / (first(colormap[i2])-first(colormap[i1]))
	weighted_color_mean(α, last(colormap[i2]), last(colormap[i1])) # NB: order because α=1 means first color in weighted_color_mean
end
lookup(::Any, ::Missing) = missing

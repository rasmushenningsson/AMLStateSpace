# See https://www.redblobgames.com/grids/hexagons/ for handy formulas

cartesian2hex(points) = [√3/3 -1/3; 0 2/3] * points
hex2cartesian(points) = [√3 √3/2; 0 3/2] * points


function _roundhex(point::SVector{2})
	q,r = point
	s = -(q+r)

	iq = round(Int,q)
	ir = round(Int,r)
	is = round(Int,s)

	dq = abs(q-iq)
	dr = abs(r-ir)
	ds = abs(s-is)

	if dq > dr && dq > ds
		iq = -(ir+is)
	elseif dr > ds
		ir = -(iq+is)
	# else
	# 	is = -(ir+iq)
	end
	@SVector [iq,ir]
end

function roundhex(points::AbstractMatrix{T}) where T
	@assert size(points,1)==2
	# Work with SVectors but return a Matrix
	p = reinterpret(SVector{2,T},vec(points))
	ip = Matrix{Int}(undef, size(points))
	ip2 = reinterpret(SVector{2,Int},vec(ip))
	ip2 .= _roundhex.(p)
	ip
end

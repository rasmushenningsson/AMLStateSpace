function checksum(fn::String, fv::VersionNumber, args...; kwargs...)
	fn[1]=='#' && return "" # lazy() with anonymous functions cannot get stable checksums

	io = IOBuffer()
	print(io, fn)
	print(io, "___SEP___")
	print(io, fv)
	for a in args
		a isa Result && isempty(checksum(a)) && return "" # propagation of anonymous
		checksum(io, a)
		print(io, "___SEP___")
	end
	print(io, "___ARGSEP___")
	for (k,v) in kwargs
		print(io,k)
		print(io, "___SEP___")
		checksum(io, v)
		print(io, "___SEP___")
	end
	bytes2hex(sha256(take!(io)))
end


function checksum(x)
	io = IOBuffer()
	checksum(io, x)
	String(take!(io))
end


checksum(x::Result) = x.ch
checksum(io::IO, x::Result) = print(io, checksum(x))

checksum(io::IO, x::Missing) = print(io, "__missing__")
checksum(io::IO, x::String) = bytes2hex(io,sha256(x))
checksum(io::IO, x::AbstractString) = checksum(io, string(x))
checksum(io::IO, x::Symbol) = (print(io, "Symbol___TYPESEP___"), bytes2hex(io,sha256(string(x))))
checksum(io::IO, x::AbstractVector{<:UInt8}) = bytes2hex(io,sha256(x))
checksum(io::IO, x::T) where T<:Number = (print(io,T,"___TYPESEP___"); bytes2hex(io,reinterpret(UInt8,[x])))

function checksum(io::IO, x::Tuple)
	print(io, length(x))
	for v in x
		print(io, "___TUPSEP___")
		checksum(io, v)
	end
end

function checksum(io::IO, x::AbstractArray)
	print(io, size(x))
	for v in x
		print(io, "___ARRSEP___")
		checksum(io, v)
	end
end

function checksum(io::IO, nt::NamedTuple)
	print(io, "__NamedTuple__")
	for (name,col) in pairs(nt)
		checksum(io, name)
		checksum(io, col)
	end
end

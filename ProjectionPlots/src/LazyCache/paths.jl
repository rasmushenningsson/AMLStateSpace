struct TimestampedFilePath{T}
	path::T
	timestamp::Float64
end
function TimestampedFilePath(path)
	@assert isfile(path) "File \"$path\" not found."
	TimestampedFilePath(path, mtime(path))
end
checksum(io::IO, x::TimestampedFilePath) = print(io, x.path, "0x", reinterpret(UInt64,x.timestamp))

struct ChecksummedFilePath{T}
	path::T
	timestamp::Float64
	checksum::String
end
function ChecksummedFilePath(x::TimestampedFilePath)
	t = mtime(x.path)
	@assert t==x.timestamp "File \"$(x.path)\" modified during analysis."
	ChecksummedFilePath(x.path, t, bytes2hex(open(sha256,x.path)))
end
ChecksummedFilePath(path) = ChecksummedFilePath(TimestampedFilePath(path))
version(::ChecksummedFilePath) = v"1.0.0"
checksum(io::IO, x::ChecksummedFilePath) = print(io, x.checksum)
function Base.string(x::ChecksummedFilePath)
	@assert mtime(x.path)==x.timestamp "File \"$(x.path)\" modified during analysis."
	string(x.path)
end
Base.show(io::IO, x::ChecksummedFilePath) = print(io, "ChecksummedFilePath(", x.checksum[1:min(6,end)], " \"", x.path, "\")")

(checksummedpath(path::TimestampedFilePath{T}; cache=true)::ChecksummedFilePath{T}) where T = cache ? cached(ChecksummedFilePath, path)[] : ChecksummedFilePath(path)
checksummedpath(path; kwargs...) = checksummedpath(TimestampedFilePath(path); kwargs...)

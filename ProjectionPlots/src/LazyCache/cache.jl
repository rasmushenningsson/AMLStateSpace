struct Cache
	dir::String
	contentsLock::ReentrantLock
	checksumLocks::Dict{String,ReentrantLock} # checksum -> checksum-specific lock
end
Cache(dir::String) = Cache(dir,ReentrantLock(),Dict{String,ReentrantLock}())

function lockchecksum(f, cache::Cache, ch::String)
	local checksumLock::ReentrantLock
	checksumLock = lock(cache.contentsLock) do
		get!(()->ReentrantLock(), cache.checksumLocks, ch)
	end
	lock(checksumLock) do
		f()
	end
end
filename(cache::Cache, ch::String) = joinpath(cache.dir, string(ch,".jld2"))


let cache::Union{Cache,Nothing} = nothing
	global function init(dir::String)
		isdir(dir) || mkdir(dir)
		cache = Cache(dir)
		nothing
	end
	global function getcache()::Cache
		@assert cache!=nothing "Please run LazyCache.init(mydir) before using cache."
		cache
	end
end



struct CachedValue end

struct LazyValue{F,S,T}
	f::F
	args::S
	kwargs::T
	useCache::Bool
end

mutable struct Result
	fname::String # name of function that produced the result, only for printout
	ch::String # checksum of input - can be "" (for lazy() used on anonymous functions)
	value::Any
end
Result(fname, ch) = Result(fname, ch,CachedValue())
Base.show(io::IO, r::Result) = print(io, "Result(", r.ch[1:min(6,end)], " ", r.fname, "->", nameof(typeof(r.value)), ')')



function Base.getindex(r::Result)
	if r.value isa LazyValue
		args2, kwargs2 = unpackargs(r.value.args, r.value.kwargs)
		fn = filename(getcache(), r.ch)
		lockchecksum(getcache(), r.ch) do
			useCache = r.value.useCache
			if useCache && isfile(fn)
				@warn "Unexpected recomputation of cached value, ignoring new value. Function: $(r.fname)."
				r.value = CachedValue() # fall through to cache loading below.
			else
				r.value = r.value.f(args2...; kwargs2...)
				useCache && save(fn, "value", r.value; compress=true)
			end
		end
	end

	r.value !== CachedValue() && return r.value
	cache = getcache()
	fn = filename(cache, r.ch)
	lockchecksum(cache, r.ch) do
		@assert isfile(fn) "Value unexpectedly missing in cache. (Checksum: $r.ch)"
		r.value = load(fn, "value")
		touch(fn) # update file timestamp (in case we want to remove old cached files later)
	end
	return r.value
end


getvalue(x) = x
getvalue(r::Result) = r[]
getvalue(a::AbstractArray{<:Result}) = getvalue.(a)

unpackargs(args, kwargs) = (getvalue.(args), map(getvalue,values(kwargs)))


function cached(f, args...; function_name_=string(nameof(f)), function_version_=version(f), kwargs...)
	cache = getcache()
	ch = checksum(function_name_, function_version_, args...; kwargs...)
	isempty(ch) && throw(ArgumentError("Anonymous functions cannot be used for caching"))
	fn = filename(cache, ch)
	fname = string(nameof(f))
	if lockchecksum(()->isfile(fn), cache, ch)
		return Result(fname, ch, CachedValue())
	end
	return Result(fname, ch, LazyValue(f,args,kwargs,true))
end
function lazy(f, args...; function_name_=string(nameof(f)), function_version_=version(f), kwargs...)
	ch = checksum(function_name_, function_version_, args...; kwargs...)
	Result(string(nameof(f)), ch, LazyValue(f,args,kwargs,false))
end


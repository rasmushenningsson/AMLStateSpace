module LazyCache

using SHA
using FileIO
using JLD2

export cached, lazy, checksummedpath

# these methods are not exported, but designed to be extended by the user
version(::Any) = v"0.0.0"
function checksum end


include("cache.jl")
include("checksum.jl")
include("paths.jl")



end

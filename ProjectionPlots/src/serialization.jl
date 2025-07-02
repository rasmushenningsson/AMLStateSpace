# Custom serialization of DataFrames since
# 1. That is more stable over time
# 2. We don't use anything fancy, e.g. Metadata
struct SerializedDataFrame
    nt::NamedTuple
end
JLD2.writeas(::Type{DataFrame}) = SerializedDataFrame
JLD2.wconvert(::Type{SerializedDataFrame}, df::DataFrame) = SerializedDataFrame(NamedTuple(pairs(eachcol(df))))
JLD2.rconvert(::Type{DataFrame}, s::SerializedDataFrame) = DataFrame(s.nt)

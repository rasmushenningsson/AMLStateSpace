const Index = Union{AbstractVector, Colon}

"""
	LowRank

A matrix decomposition `UVᵀ` where each row of `U` represents a variable and each column of `Vᵀ` represents a sample.
Intended for situations where the product is low rank, i.e. `size(U,2)==size(Vt,1)` is small.
"""
struct LowRank{T1,T2}
	U::T1
	Vt::T2
end

Base.size(F::LowRank) = (size(F.U,1),size(F.Vt,2))
Base.size(F::LowRank, dim::Integer) = dim<2 ? size(F.U,dim) : size(F.Vt,dim)

const LowRankFactorization = Union{SVD,PMA,LowRank}


innersize(F::LowRankFactorization) = length(F.S)
innersize(F::LowRank) = size(F.U,2)


var_coordinates(F::LowRankFactorization) = F.U
obs_coordinates(F::Union{SVD,PMA}) = Diagonal(F.S)*F.Vt
obs_coordinates(F::LowRank) = F.Vt


function _subsetmatrix(F::Union{SVD,PMA}, I::Index, J::Index)
	U = F.U[I,:]
	Vt = F.Vt[:,J]
	lmul!(Diagonal(F.S), Vt)
	LowRank(U, Vt)
end

_subsetmatrix(F::LowRank, I::Index, J::Index) = LowRank(F.U[I,:], F.Vt[:,J])


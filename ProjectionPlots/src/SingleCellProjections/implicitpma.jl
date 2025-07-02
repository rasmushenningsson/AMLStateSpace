function _implicitpma(A, samplekernelroot::AbstractMatrix; kwargs...)
	prod(size(A))==0 && return PMA(zeros(0,0), zeros(0), zeros(0,0), zeros(0,0))
	B = matrixproduct(A, samplekernelroot)
	F = implicitsvd(B; kwargs...)

	U = F.U
	S = F.S
	V = A'U
	V ./= S'
	PMA(U,S,Matrix(V'),Matrix(F.V))
end

"""
	implicitpma(A, G::SimplexGraph; nsv=3, subspacedims=8nsv, niter=2)

Computes the Principal Moment Analysis of the implicitly given matrix `A` (variables × samples) using the sample simplex graph `G`.
"""
implicitpma(A, G::SimplexGraph; kwargs...) = _implicitpma(A, simplices2kernelmatrixroot(G; simplify=false); kwargs...)
implicitpma(A, G::AbstractMatrix{Bool}; kwargs...) = implicitpma(A, SimplexGraph(G); kwargs...)

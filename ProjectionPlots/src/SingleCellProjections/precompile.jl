using SnoopPrecompile

@precompile_setup begin
	test_path = joinpath(@__DIR__, "../../test")
	h5_path = joinpath(test_path, "SingleCellProjections/SingleCell10x/data/500_PBMC_3p_LT_Chromium_X_50genes/filtered_feature_bc_matrix.h5")
	@precompile_all_calls begin
		counts = load10x(h5_path)
		counts2 = load_counts([h5_path,h5_path]; sample_names=["a","b"])
		transformed = sctransform(counts; use_cache=false, verbose=false)
		centered = normalize_matrix(transformed)
		reduced = svd(centered; nsv=10)
		fl = force_layout(reduced; ndim=3, k=100)
	end
end
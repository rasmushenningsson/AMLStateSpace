# normal samples
const samplesNormal = ["NBM4-MNC", "NBM7-MNC", "NBM8-MNC", "NBM10-MNC", "NBM11-MNC", "NBM10-CD34", "NBM11-CD34", "NBM8-CD34"]

# AML samples
const samplesAML = ["AML10","AML104","AML105","AML110","AML111","AML117","AML123","AML124","AML126","AML136","AML138","AML151","AML155","AML157","AML161","AML172","AML21","AML23","AML24","AML25","AML27","AML28","AML32","AML33","AML34","AML37","AML48","AML55","AML61","AML62","AML7","AML79","AML80","AML83","AML85D","AML85R","AML9","AML97"]
const samplesTP53 = ["AML37","AML48","AML62","AML80","AML83","AML85D","AML85R"] # with AML85 Relapse
const samplesNPM1 = ["AML79","AML27","AML28","AML25","AML7","AML9","AML24","AML33","AML97","AML104","AML105","AML136"]
const samplesMR = ["AML32","AML21","AML110","AML138","AML10","AML55","AML111","AML126","AML151","AML157","AML172"]
const samplesCBFBMYH11 = ["AML23","AML123","AML124"] # inv16
const samplesRUNX1RUNX1T1 = ["AML61","AML117","AML161"] # t(8;21)(q22;q22)
const samplesOther = ["AML155", "AML34"]

# NPM1 groups
const samplesNPM1_ClassA = ["AML24","AML25","AML27","AML28","AML79","AML104","AML136"]
const samplesNPM1_ClassB = ["AML7","AML9","AML33","AML97","AML105"]

const samplesCluster_C3a = ["AML24", "AML25", "AML27", "AML28", "AML79", "AML104"]
const samplesCluster_C3b = ["AML7", "AML9", "AML33", "AML97", "AML105"]

const samplesComplete = union(samplesAML,samplesNormal)

function update_vaf_tags(x)
	x = replace(x, "<u>"=>"<span style='text-decoration:underline'>")
	x = replace(x, "</u>"=>"</span>")
	x = replace(x, "<c>"=>"<span style='color:grey'>")
	x = replace(x, "</c>"=>"</span>")
end
function vafstrings() # See above for how tags are used
	d = Dict{String,String}(
		"AML10"  => "<i>JAK2</i> (47%), <i>ASXL1</i> (43%), <b><u><i>RUNX1</i> (41%)</u></b>, <b><u><i>SRSF2</i> (41%)</u></b>, <b><u><i>IDH2</i> (39%)</u></b>",
		"AML104" => "<b><u><i>NPM1</i> (59%)</u></b>, <i>FLT3</i>-ITD (54%), <b><u><i>DNMT3A</i> (46%)</u></b>, <b><u><i>IDH2</i> (43%)</u></b>, <b><c><i>PCMTD1</i> (26%)</c></b>",
		"AML105" => "<b><u><i>DNMT3A</i> (43%)</u></b>, <b><u><i>NPM1</i> (37%)</u></b>, <i>FLT3</i>-ITD (31%)",
		"AML110" => "<b><c><u><i>MKI67</i> (46%)</u></c></b>, <i>STAG2</i> (45%), <b><u><i>SRSF2</i> (40%)</u></b>, <b><u><i>IDH2</i> (37%)</u></b>",
		"AML111" => "<i>BCOR</i> (59%), <b><u><i>RUNX1</i> (40%)</u></b>, <b><c><u><i>PAM16</i> (32%)</u></c></b>, <b><c><u><i>NARS</i> (31%)</u></c></b>",
		"AML117" => "<i>RUNX1</i>::<i>RUNX1T1</i>, <b><c><u><i>LRRFIP1</i> (50%)</u></c></b>, <b><c><u><i>TSPOAP1</i> (50%)</u></c></b>, <b><i>NRAS</i> (15%)</b>",
		"AML123" => "<b><i>CBFB</i>::<i>MYH11</i></b>, <i>SETBP1</i> (38%)",
		"AML124" => "<b><i>CBFB</i>::<i>MYH11</i></b>, <b><u><i>KIT</i> (43%)</u></b>, <b><c><i>RTCA</i> (36%)</c></b>, <b><c><i>COPS4</i> (8%)</c></b>",
		"AML126" => "<b><u><i>CEBPA</i> (38%)</u></b>, <b><i>RAD21</i> (32%)</b>, <b><u><i>IDH2</i> (29%)</u></b>, <b><u><i>CBFB</i> (26%)</u></b>, <i>ASXL1</i> (17%)",
		"AML136" => "<b><u><i>NPM1</i> (40%)</u></b>, <b><u><i>IDH2</i> (36%)</u></b>, <b><i>EFHD2</i> (17%)</b>",
		"AML138" => "<i>ZRSR2</i> (89%), <i>EZH2</i> (54%), <i>TET2</i> (41%), <b><c><u><i>GRHPR</i> (39%)</u></c></b>, <b><c><u><i>PRPF38B</i> (35%)</u></c></b>, <b><u><i>RUNX1</i> (25%)</u></b>",
		"AML144" => "2x<i>EZH2</i> (52%;51%), <b><c><u><i>FOXP1</i> (46%)</u></c></b>, <i>ASXL1</i> (41%), <b><u><i>RUNX1</i> (39%)</u></b>, <b><c><u><i>BRCC3</i> (39%)</u></c></b>,  2x<i>TET2</i> (39%;35%), <i>BCOR</i> (37%)",
		"AML151" => "<b><i>RUNX1</i> (100%)</b>, <b><c><i>SOX4</i> (81%) </c></b>, <i>TET2</i> (98%), <i>PHF6</i> (95%), <b><u><i>SRSF2</i> (53%)</u></b>, <i>IDH1</i> (38%)",
		"AML154" => "<i>SF3B1</i> (47%), <i>PTPN11</i> (40%), <i>RUNX1</i> (37%), <i>BCOR</i> (25%)", # Not present in data
		"AML155" => "<i>KMT2A</i>::<i>MLLT3</i>, <b><u><i>FLT3</i> (27%)</u></b>, <b><u><i>TP53</i> (27%)</u></b>",
		"AML157" => "<i>BCOR</i> (88%), <b><c><u><i>MS4A6A</i> (54%)</u></c></b>, <b><u><i>ETV6</i> (53%)</u></b>, <b><u><i>SRSF2</i> (53%)</u></b>, 2x<i>RUNX1</i> (43%;42%), <i>NRAS</i> (41%)",
		"AML161" => "<i>RUNX1</i>::<i>RUNX1T1</i>, <i>FLT3</i>-ITD (45%), <b><c><u><i>TTC1</i> (40%)</u></c></b>, <i>RAD21</i> (37%), <b><c><u><i>FGFR1OP2</i> (31%)</u></c></b>, <i>KIT</i> (9%), <b><i>KRAS</i> (8%)</b>",
		"AML172" => "<i>FLT3</i>-ITD (67%), <b><u><i>RUNX1</i> (44%)</u></b>, <b><c><u><i>FYB1</i> (43%)</u></c></b>, <b><c><u><i>SEC31A</i> (43%)</u></c></b>, <i>ASXL1</i> (38%), <i>BCORL1</i> (7%)",
		"AML21"  => "<i>STAG2</i> (88%), <i>BCOR</i> (84%), <b><u><i>IDH2</i> (46%)</u></b>, <b><u><i>SRSF2</i> (33%)</u></b>",
		"AML23"  => "<b><i>CBFB</i>::<i>MYH11</i></b>, <b><u><i>NRAS</i> (43%)</u></b>, <i>SF3B1</i> (12%)",
		"AML24"  => "<b><u><i>NPM1</i> (50%)</u></b>, <b><u><i>IDH2</i> (46%)</u></b>, <i>FLT3</i>-ITD (18%)",
		"AML25"  => "<i>TET2</i> (91%), <b><c><u><i>MCM2</i> (50%)</u></c></b>, <b><c><u><i>DOCK8</i> (48%)</u></c></b>, <b><u><i>NPM1</i> (46%)</u></b>, <b><c><u><i>CD37</i> (40%)</u></c></b>",
		"AML27"  => "<i>FLT3</i>-ITD (59%), <b><u><i>IDH1</i> (52%)</u></b>, <b><c><u><i>UBR2</i> (51%)</u></c></b>, <b><u><i>NPM1</i> (40%)</u></b>, <b><c><i>IK</i> (13%)</c></b>",
		"AML28"  => "2x<i>TET2</i> (48%;47%), <b><c><u><i>METTL5</i> (42%)</u></c></b>, <b><c><u><i>GSTK1</i> (41%)</u></c></b>, <i>FLT3</i>-ITD (40%), <b><c><u><i>MAPRE2</i> (38%)</u></c></b>, <b><u><i>NPM1</i> (36%)</u></b>",
		"AML32"  => "<i>ASXL1</i> (54%), <b><u><i>SRSF2</i> (48%)</u></b>, 2x<b><u><i>RUNX1</i></u></b> (46%;<b><u>26%</u></b>), <b><c><u><i>LAMTOR5</i> (37%)</u></c></b>, <i>STAG2</i> (30%)",
		"AML33"  => "<b><u><i>NPM1</i> (58%)</u></b>, <i>DNMT3A</i> (53%), <i>FLT3</i>-ITD (48%), <b><c><u><i>MRPL52</i> (48%)</u></c></b>, <b><u><i>IDH1</i> (45%)</u></b>",
		"AML34"  => "<i>DNMT3A</i> (46%), <i>NRAS</i> (43%), <i>CEBPA</i> (42%), <b><u><i>IDH1</i> (41%)</u></b>, <b><c><i>PRIM2</i> (16%)</c></b>, <b><c><i>SRRM1</i> (8%)</c></b>",
		"AML35"  => "<b><u><i>RUNX1</i> (60%)</u></b>, <i>DNMT3A</i> (43%), <b><u><i>IDH2</i> (35%)</u></b>, <i>ASXL1</i> (29%)",
		"AML37"  => "<b><u><i>TP53</i> (62%)</u></b>, <b><i>FLT3</i> (25%)</b>, <b><c><i>RPL28</i> (17%)</c></b>, <b><c><i>EGR1</i> (13%)</c></b>",
		"AML48"  => "<b><c><i>SMG1</i> (64%)</c></b>, <b><c><u><i>DOCK8</i> (56%)</u></c></b>, 2x<b><u><i>TP53</i></u></b> (<b><u>44%</u></b>;41%)",
		"AML55"  => "<i>ASXL1</i> (39%), <i>EZH2</i> (36%), <b><c><u><i>PFN1</i> (35%)</u></c></b>, <i>PTPN11</i> (26%), <b><i>SRSF2</i> (26%)</b>",
		"AML61"  => "<i>RUNX1</i>::<i>RUNX1T1</i>, <b><c><u><i>SCP2</i> (46%)</u></c></b>, <b><c><u><i>BCL7C</i> (44%)</u></c></b>",
		"AML62"  => "<b><i>TP53</i> (69%)</b>, <i>DNMT3A</i> (44%), <i>NF1</i> (44%), <b><c><u><i>MRTFA</i> (30%)</u></c></b>, <b><c><u><i>TNRC18</i> (20%)</u></c></b>",
		"AML7"   => "<i>DNMT3A</i> (38%), <b><u><i>NPM1</i> (38%)</u></b>, 2x<i>WT1</i> (19%;19%), <b><i>IDH2</i> (11%)</b>, <b><i>KRAS</i> (6%)</b>",
		"AML70"  => "<b><i>EZH2</i> (99%)</b>, <b><u><i>RUNX1</i> (54%)</u></b>, <i>ASXL1</i> (53%) , 2x<i>TET2</i> (49%;39%)",
		"AML79"  => "<i>FLT3</i>-ITD (63%), <b><u><i>NPM1</i> (48%)</u></b>, <b><u><i>IDH2</i> (44%)</u></b>, <b><c><u><i>ENTPD4</i> (44%)</u></c></b>, <b><c><u><i>NHEJ1</i> (36%)</u></c></b>",
		"AML80"  => "<b><i>TP53</i> (94%)</b>, <b><c><u><i>BANF1</i> (45%)</u></c></b>, <b><c><u><i>R3HDM4</i> (41%)</u></c></b>, <b><c><u><i>PRKCSH</i> (33%)</u></c></b>",
		"AML83"  => "<b><i>TP53</i> (88%)</b>, <i>NF1</i> (87%), <b><c><u><i>MCM4</i> (58%)</u></c></b>, <i>TET2</i> (41%), <b><c><u><i>ERBIN</i> (28%)</u></c></b>, <b><c><i>CTNNB1</i> (17%)</c></b>",
		"AML85D" => "<b><u><i>TP53</i> (53%)</u></b>, <i>DNMT3A</i> (38%), <b><c><i>SEPHS2</i> (25%)</c></b>",
		"AML85R" => "<b><i>TP53</i> (97%)</b>, <i>DNMT3A</i> (71%), <b><c><i>SEPHS2</i> (25%)</c></b>",
		"AML9"   => "<i>DNMT3A</i> (40%), <b><u><i>NPM1</i> (39%)</u></b>, 2x<i>TET2</i> (34%;29%), <b><c><u><i>ECHDC1</i> (31%)</u></c></b>",
		"AML97"  => "<b><u><i>IDH1</i> (46%)</u></b>, <i>DNMT3A</i> (42%), <b><u><i>NPM1</i> (33%)</u></b>, <b><c><u><i>BCAP31</i> (31%)</u></c></b>, <i>NRAS</i> (25%)",
		"NPM1 Class I"   => string("n=", length(samplesNPM1_ClassA)),
		"NPM1 Class II"  => string("n=", length(samplesNPM1_ClassB)),
		"Cluster C3a"    => string("n=", length(samplesCluster_C3a)),
		"Cluster C3b"    => string("n=", length(samplesCluster_C3b)),
		"NPM1"           => string("n=", length(samplesNPM1)),
		"AML-MR"         => string("n=", length(samplesMR)),
		"TP53"           => string("n=", length(samplesTP53)),
		"CBFB::MYH11"    => string("n=", length(samplesCBFBMYH11)),
		"RUNX1::RUNX1T1" => string("n=", length(samplesRUNX1RUNX1T1)),
	)
	Dict((k=>update_vaf_tags(v) for (k,v) in d))
end



blue_gradient() = [(0.0,RGB(19/255,43/255,67/255)), (1.0,RGB(86/255,177/255,247/255))]




l4order() = ["HSC","LMPP","GMP","Monocytes","Megakaryocytic cells","Erythroid cells","Dendritic cells","NK-cells","T-cells","B-cells"]

colordict_celltype_l4() = Dict(
	"HSC"                  => colorant"#66b266",
	"Monocytes"            => colorant"#fa9fb5",
	"GMP"                  => colorant"#008000",
	"Megakaryocytic cells" => colorant"#e5687e",
	"LMPP"                 => colorant"#756bb1",
	"Erythroid cells"      => colorant"#a7241d",
	"B-cells"              => colorant"#a64ca6",
	"T-cells"              => colorant"#7fb4e5",
	"NK-cells"             => colorant"#196fbe",
	"Dendritic cells"      => colorant"#ff4d00",
)



function setup_args(; plotmode::Symbol=:groups, aml_immature_only=false)
	@assert plotmode in (:base, :groups, :samples)

	ndim = 2
	pca_args = (;nsv=40, seed=1234, niter=3)
	embed_args = (;k=10)

	# - Eliminated Factors -
	elim_args = (("percent.mt",:numerical),)

	# - Cell Filtering -
	# Remove cells with percent.mt Unknown or >=15
	cellfilterfun = row -> row["percent.mt"]!==missing && row["percent.mt"]<15


	if aml_immature_only
		cellfilterfun_sample = row -> row["percent.mt"]!==missing && row["percent.mt"]<15 && row["celltype.aml"]=="AML Immature"
	else
		cellfilterfun_sample = row -> row["percent.mt"]!==missing && row["percent.mt"]<15 # Use me
	end


	# - Colors -
	color_scales = (; score=scientificcolourmap("lajolla"),
	                  heatmap_score=scientificcolourmap("lajolla"),
	                  cell_counts=scientificcolourmap("lajolla"),
	                  numerical=scientificcolourmap("lajolla"),
	                  mutations=blue_gradient(),
	               )

	force_args = (;k=160, ndim, niter=400, link_distance=40, link_strength=0.05, charge=40, seed=9560) # Latest

	# Custom rotation
	rotation = [1 0; 0 -1]*rotmatrix2d(π/4) # 9560

	bgColor = RGB(0.85,0.85,0.85) # gray - for mutation plots
	bgOpacity = 1.0

	bgMarkerSize = 3
	markerSize = 6
	fontSize = 13
	overlayFontSize = 22
	heatmap_mincells = plotmode==:samples ? 3 : 5
	subdivisions = 40

	cmin = (; score=0,   score_mean=95,  score_q90=95)
	cmax = (; score=450, score_mean=310, score_q90=310)

	offsetDict = Dict{String,Vector}(
		"HSC"                  => [1250,200],
		"Monocytes"            => [600,-1250],
		"GMP"                  => [300,1150],
		"Megakaryocytic cells" => [2200,500],
		"LMPP"                 => [-600,800],
		"Erythroid cells"      => [1400,-900],
		"B-cells"              => [-1000,-850],
		"T-cells"              => [1050,-1250],
		"NK-cells"             => [-1200,50],
		"Dendritic cells"      => [-1000,2600],
	)
	category_overlay_args = (;fontSize=30, minCells=50, offsets=offsetDict) # Use me

	heatmap_size_overlay_args = (;pos=[0.05, 0.10])

	axis_kwargs = (;showline=true, linewidth=4, linecolor="black")
	# axis_kwargs = (;showline=true, linewidth=4, linecolor="black", showgrid=true, nticks=40, showticklabels=true, gridcolor="black")

	plotFormats=["png"]
	# plotFormats=["pdf"]


	titlefun = (p,d)->begin
		t = "$(p.plotTitle), $(get_ncells_string(d)) cells."
		t = string("<span style='font-size:20pt'>", t, "</span>")

		t2 = get(vafstrings(),p.plotTitle,"")

		t = isempty(t2) ? t : string(t,"<br>",t2)
		t = string("<span style='color:black'>", t, "</span>")
	end

	samplePath = get_samples_path()
	annotationPaths = joinpath.(get_annotations_path(), ["scRNA_AML.csv", "mutation_status.csv"])

	args = (; baseIDs=samplesNormal,
	          annotationPaths, samplePath, featureType="Gene Expression", ndim, pca_args, force_args, embed_args, elim_args,
	          cellfilterfun, cellfilterfun_sample,
	          plotmode,
	          subdivisions, heatmap_mincells,
	          color_scales, cmin, cmax, category_overlay_args,
	          heatmap_size_overlay_args,
	          axisPrefix="",
	          axis_kwargs,
	          overlay_celltype="celltype.l4",
	          titlefun,
	          rotation,
	          bgMarkerSize, bgColor, bgOpacity,
	          markerSize, lineWidth=0, fontSize, overlayFontSize,
	          displayPlots=false, savePlots=true, plotFormats, plotWidth=1024, plotHeight=1024,
	          plotScale=2, # Not used for paper
	       )
end


# Figure parts
function nbm_force_plot(outPath)
	plotfuns = [
		plot_categorical("celltype.l4"; colorDict=colordict_celltype_l4(), outName="celltypebroad", legendTitle=nothing, order=l4order())
	]
	main((; setup_args(plotmode=:base)..., outPath, plotfuns, plotTitle="NBM"))
end
function sample_mutation_plots(outPath, sampleIDs)
	plotfuns = [plot_mutation_heatmap]
	main_samples((;setup_args(plotmode=:samples)..., outPath, plotfuns, sampleIDs, plotTitle=sampleIDs))
end

function group_mutation_plot(outPath, plotTitle, sampleIDs)
	plotfuns = [plot_mutation_heatmap]
	main((;setup_args()..., outPath, plotfuns, plotTitle, sampleIDs))
end

function aml_immature_counts_plot(outPath)
	plotfuns = [
		plot_cellcount(;normalize=false, cmin=0.0, logtransform=false, cellcountmin=100, outName="cellcount_linear_min100cells")
	]
	main((;setup_args(aml_immature_only=true)..., outPath, plotfuns, plotTitle="AML Immature", sampleIDs=samplesAML))
end



# Figures
function figure2a()
	outPath = "out/figure2a"
	nbm_force_plot(outPath)
	sample_mutation_plots(outPath, ["AML61", "AML124", "AML83", "AML28", "AML9", "AML151", "AML172"])
end

function figure3b()
	outPath = "out/figure3b"
	aml_immature_counts_plot(outPath)
end


function supplementary_figure6()
	outPath = "out/supplementary_figure6"
	nbm_force_plot(outPath)
	group_mutation_plot(outPath, "TP53",           samplesTP53)
	group_mutation_plot(outPath, "NPM1",           samplesNPM1)
	group_mutation_plot(outPath, "AML-MR",         samplesMR)
	group_mutation_plot(outPath, "CBFB::MYH11",    samplesCBFBMYH11)
	group_mutation_plot(outPath, "RUNX1::RUNX1T1", samplesRUNX1RUNX1T1)
end

function supplementary_figure7()
	outPath = "out/supplementary_figure7"
	sample_mutation_plots(outPath, vcat(samplesNPM1, samplesMR))
end

function supplementary_figure8()
	outPath = "out/supplementary_figure8"
	sample_mutation_plots(outPath, vcat(samplesTP53, samplesCBFBMYH11, samplesRUNX1RUNX1T1, samplesOther))
end

function supplementary_figure11b()
	outPath = "out/supplementary_figure11b"
	group_mutation_plot(outPath, "Cluster C3a", samplesCluster_C3a)
	group_mutation_plot(outPath, "Cluster C3b", samplesCluster_C3b)
end

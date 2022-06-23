### A Pluto.jl notebook ###
# v0.19.8

using Markdown
using InteractiveUtils

# ╔═╡ c82bb3e6-e052-11ec-3392-eba4a4b464ba
# ╠═╡ show_logs = false
begin
	let
		using Pkg
		Pkg.activate(mktempdir())
		Pkg.Registry.update()
		Pkg.add("PlutoUI")
		Pkg.add("Statistics")
		Pkg.add("StatsBase")
		Pkg.add("ImageMorphology")
		Pkg.add("ImageFiltering")
		Pkg.add("CSV")
		Pkg.add("DataFrames")
		Pkg.add("GLM")
		Pkg.add(url="https://github.com/JuliaHealth/DICOM.jl")
		Pkg.add(url="https://github.com/Dale-Black/DICOMUtils.jl")
		Pkg.add(url="https://github.com/Dale-Black/PhantomSegmentation.jl")
		Pkg.add(url="https://github.com/Dale-Black/CalciumScoring.jl")
		Pkg.add("CairoMakie")
		Pkg.add("HypothesisTests")
		Pkg.add("Colors")
		Pkg.add("MLJBase")
		Pkg.add("AlgebraOfGraphics")
	end
	
	using PlutoUI
	using Statistics
	using StatsBase: quantile!, rmsd
	using ImageMorphology
	using ImageFiltering
	using CSV
	using DataFrames
	using GLM
	using DICOM
	using DICOMUtils
	using PhantomSegmentation
	using CalciumScoring
	using CairoMakie
	using HypothesisTests
	using Colors
	using MLJBase
	using AlgebraOfGraphics
end

# ╔═╡ 2df207ed-cf9c-4d7d-8354-4da14f93276c
TableOfContents()

# ╔═╡ 4cf61487-3913-4e11-970e-4ca31a5ffc8d
md"""
## Integrated
"""

# ╔═╡ 57d46998-368a-4916-8380-ee49d5473a49
path_integrated = "/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/output_new/integrated_scoring";

# ╔═╡ 2c960bd8-ae64-453f-b29f-275bf5263774
df_i = CSV.read(string(path_integrated, "/full.csv"), DataFrame);

# ╔═╡ 7f4fae09-7916-4a98-a102-7f861900c457
df_i_low, df_i_normal = groupby(df_i, :DENSITY);

# ╔═╡ c987e55d-9ffd-40a5-9d64-1653d4000837
df_i_low_small, df_i_low_medium, df_i_low_large = groupby(df_i_low, :SIZE);

# ╔═╡ 0c367cb1-d959-4c10-b31c-18d3cb396255
df_i_normal_small, df_i_normal_medium, df_i_normal_large = groupby(df_i_normal, :SIZE);

# ╔═╡ a58beea4-b186-44c6-bf75-98d368f7ebe4
md"""
## Agatston
"""

# ╔═╡ 4fd5b649-cea1-4d1f-95ac-35aa12f8fbff
path_agat = "/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/output_new/agatston";

# ╔═╡ 3fa3748a-22d6-49d8-9888-a749dca99ce2
df_a = CSV.read(string(path_agat, "/full.csv"), DataFrame);

# ╔═╡ c38b9994-7698-4c81-b65b-534b994d8ff0
df_a_low, df_a_normal = groupby(df_a, :DENSITY);

# ╔═╡ 79952ee5-3ccd-4d85-95b4-f7610c90d63d
df_a_low_small, df_a_low_medium, df_a_low_large = groupby(df_a_low, :SIZE);

# ╔═╡ 7223ae32-7241-463f-be13-2d2ce6bc7de7
df_a_normal_small, df_a_normal_medium, df_a_normal_large = groupby(df_a_normal, :SIZE);

# ╔═╡ f39a9ed5-c897-4b90-8494-bf231312970e
md"""
## Spatially Weighted
"""

# ╔═╡ 32623778-9a19-4318-81c6-1fc7aa0157fe
path_swcs = "/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/output_new/swcs";

# ╔═╡ 317f8c26-7b39-45c3-b9d3-02925d4b0514
df_s = CSV.read(string(path_swcs, "/full.csv"), DataFrame);

# ╔═╡ 06deefa7-9170-4e2f-ac49-d6dc65c0af76
df_s_low, df_s_normal = groupby(df_s, :DENSITY);

# ╔═╡ 5a4cc121-fd02-49bf-bdf5-07b71b48ce19
df_s_low_small, df_s_low_medium, df_s_low_large = groupby(df_s_low, :SIZE);

# ╔═╡ d5d8e9da-eba8-4788-b787-95b3a2b41329
df_s_normal_small, df_s_normal_medium, df_s_normal_large = groupby(df_s_normal, :SIZE);

# ╔═╡ 310c1999-77a6-432e-850c-4111411cceb0
md"""
# Figures
"""

# ╔═╡ 3a8c6f3c-0fa6-477c-85e2-7739cf2537ae
md"""
## Zero CAC Scores (FIG)
"""

# ╔═╡ 932e1a11-ec52-48bc-af27-bb3e09e6b27f
md"""
#### SWCS
"""

# ╔═╡ dc8805a5-2b4c-47fe-88c8-0f4434ddd5cf
df_s_80, df_s_100, df_s_120, df_s_135, = groupby(df_s, :scan);

# ╔═╡ 2c78dcc4-03c7-49b5-ac45-30c198671761
begin
	mean_swcs_80, std_swcs_80= mean(df_s_80[!, :swcs_bkg]), std(df_s_80[!, :swcs_bkg])
	mean_swcs_100, std_swcs_100= mean(df_s_100[!, :swcs_bkg]), std(df_s_100[!, :swcs_bkg])
	mean_swcs_120, std_swcs_120= mean(df_s_120[!, :swcs_bkg]), std(df_s_120[!, :swcs_bkg])
	mean_swcs_135, std_swcs_135 = mean(df_s_135[!, :swcs_bkg]), std(df_s_135[!, :swcs_bkg])
end;

# ╔═╡ d3e13dcd-201a-49ab-86f7-a7d3320088c2
begin
	array_s_80 = Array(df_s_80[!, 12:14])
	array_s_100 = Array(df_s_100[!, 12:14])
	array_s_120 = Array(df_s_120[!, 12:14])
	array_s_135 = Array(df_s_135[!, 12:14])
end;

# ╔═╡ 4b97d892-b0db-4f75-9239-9da215f905f0
begin
	num_zeroCAC_80 = length(findall(x -> x < mean_swcs_80 - std_swcs_80, array_s_80))
	num_zeroCAC_100 = length(findall(x -> x < mean_swcs_100 - std_swcs_100, array_s_100))
	num_zeroCAC_120 = length(findall(x -> x < mean_swcs_120 - std_swcs_120, array_s_120))
	num_zeroCAC_135 = length(findall(x -> x < mean_swcs_135 - std_swcs_135, array_s_135))

	total_zero_s = num_zeroCAC_80 + num_zeroCAC_100 + num_zeroCAC_120 + num_zeroCAC_135
	total_cac = length(array_s_80) * 4
end;

# ╔═╡ 01d903d8-fca2-4317-a2d7-93429f06d47c
md"""
#### Agatston
"""

# ╔═╡ 97a612d5-44e7-43f1-b61e-8518c7e6aa28
array_a = hcat(df_a[!, 7], df_a[!, 9], df_a[!, 11]);

# ╔═╡ a31dfad2-e420-4c41-b805-227778a39fc9
num_zero_a = length(findall(x -> x == 0, array_a));

# ╔═╡ 86d5f7fd-6a1b-4fdc-9e06-6cd66f8a169c
md"""
#### Integrated
"""

# ╔═╡ aabe8e7b-faed-499d-aea2-fd44836a3cfb
df_i_80, df_i_100, df_i_120, df_i_135, = groupby(df_i, :scan);

# ╔═╡ dd644657-51e9-4b27-a186-6efba0bb83f8
begin
	mean_i_80, std_i_80= mean(df_i_80[!, :mass_bkg]), std(df_i_80[!, :mass_bkg])
	mean_i_100, std_i_100= mean(df_i_100[!, :mass_bkg]), std(df_i_100[!, :mass_bkg])
	mean_i_120, std_i_120= mean(df_i_120[!, :mass_bkg]), std(df_i_120[!, :mass_bkg])
	mean_i_135, std_i_135 = mean(df_i_135[!, :mass_bkg]), std(df_i_135[!, :mass_bkg])
end

# ╔═╡ d79d6478-0837-40f0-a7a5-270bf44d6166
begin
	array_i_80 = hcat(df_i_80[!, 7], df_i_80[!, 9], df_i_80[!, 11])
	array_i_100 = hcat(df_i_100[!, 7], df_i_100[!, 9], df_i_100[!, 11])
	array_i_120 = hcat(df_i_120[!, 7], df_i_120[!, 9], df_i_120[!, 11])
	array_i_135 = hcat(df_i_135[!, 7], df_i_135[!, 9], df_i_135[!, 11])
end;

# ╔═╡ 2b5a8dca-faae-4163-ab14-ae4f98b2a7cc
begin
	num_zeroCAC_80_i = length(findall(x -> x < mean_i_80 - std_i_80, array_i_80))
	num_zeroCAC_100_i = length(findall(x -> x < mean_i_100 - std_i_100, array_i_100))
	num_zeroCAC_120_i = length(findall(x -> x < mean_i_120 - std_i_120, array_i_120))
	num_zeroCAC_135_i = length(findall(x -> x < mean_i_135 - std_i_135, array_i_135))

	total_zero_i = num_zeroCAC_80_i + num_zeroCAC_100_i + num_zeroCAC_120_i + num_zeroCAC_135_i
end;

# ╔═╡ b702760a-9d9b-4d2a-992c-37ed4f8bcfdf
let
	f = Figure()
	f[1, 1]

	colors = Makie.wong_colors()
	labels = ["Number Zero CAC (SWCS)", "Number Zero CAC (Agatston)"]
	elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]

	##-- TOP --##
	axtop = Axis(f[1, 1], xticks = (1:3, ["Integrated", "Spatially Weighted", "Agatston"]))
	ylims!(axtop, low=0, high=100)
	axtop.yticks = [0, 25, 50, 75, 100]
	table = [1, 2, 3]
	heights1 = [(total_zero_i / total_cac) * 100, (total_zero_s / total_cac) * 100, (num_zero_a / total_cac) * 100]
	barplot!(axtop, table, heights1, color=colors[1:3], bar_labels=:y)
	hidedecorations!(axtop, ticklabels=false, ticks=false)
	
	Label(f[1:end, 1:end, Top()], "Zero CAC Scores", valign = :center, padding = (0, 0, 0, 0), textsize=25)
	Label(f[1:end, 1:end, Left()], "% false-negative zero CAC scores", valign = :center, padding = (0, 50, 0, 0), rotation=π/2, textsize=17)

	save("/Users/daleblack/Google Drive/Research/2022-AAPM/zero_cac.png", f)
	f
end

# ╔═╡ f4b51729-d53a-40e2-919b-5d4562cffbc0
total_zero_i, total_zero_s, num_zero_a

# ╔═╡ da2b55dc-2fe9-468f-a53e-baf232ca79cd
total_cac

# ╔═╡ 7dfc24a4-e006-45f4-b5b9-977a7c3c0b7c
md"""
## Linear Regression Mass Score (FIG)
"""

# ╔═╡ d04ce8a0-8224-49c7-a057-f1c7d0d65ef8
# let
# 	f = Figure()

# 	##-- A --##
# 	axtop = Axis(f[1, 1])
	
# 	df = df_i_normal
# 	scatter!(axtop, df[!, :ground_truth_mass_large], df[!, :calculated_mass_large])
# 	errorbars!(axtop, df[!, :ground_truth_mass_large], df[!, :calculated_mass_large], rms(df[!, :ground_truth_mass_large], df[!, :calculated_mass_large]))
# 	scatter!(axtop, df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium])
# 	errorbars!(axtop, df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium], rms(df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium]))
# 	scatter!(axtop, df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], color=:red)
# 	errorbars!(axtop, df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], rms(df[!, :ground_truth_mass_small], df[!, :calculated_mass_small]))
# 	lines!(axtop, [-1000, 1000], [-1000, 1000],)
# 	lines!(axtop, collect(1:1000), pred_i_norm, label = L"y=%$(trunc(co1[2]; digits=3))x + %$(trunc(co1[1]; digits=3))", linestyle=:dashdot)
# 	# axislegend(; position=(0,1))
# 	Textbox(
# 		f[1, 1], 
# 		placeholder = "y = $(trunc(co1[2]; digits=3))x + $(trunc(co1[1]; digits=3)) \nr=$(trunc(r2_1; digits=3)) \nRMSE: $(trunc(rms_values1[1]; digits=3)) \nRMSD: $(trunc(rms_values1[2]; digits=3))", 
# 		tellheight = false,
#         tellwidth = false,
# 		boxcolor=:white,
# 		halign=:left,
# 		valign=:top,
# 		textsize=12
# 	)

# 	xlims!(axtop, low=0, high=200)
# 	ylims!(axtop, low=0, high=200)
# 	axtop.xticks = [0, 50, 100, 150, 200]
# 	axtop.yticks = [0, 50, 100, 150, 200]
# 	axtop.xlabel = "Known Mass (mg)"
# 	axtop.ylabel = "Calculated Mass (mg)"
# 	axtop.title = "Integrated (Normal-Density)"
# 	hidedecorations!(axtop, ticklabels=false, ticks=false)

# 	##-- B --##
# 	axbottom = Axis(f[2, 1])
	
# 	df2 = df_i_low
# 	scatter!(axbottom, df2[!, :ground_truth_mass_large], df2[!, :calculated_mass_large])
# 	errorbars!(axbottom, df2[!, :ground_truth_mass_large], df2[!, :calculated_mass_large], rms(df2[!, :ground_truth_mass_large], df2[!, :calculated_mass_large]))
# 	scatter!(axbottom, df2[!, :ground_truth_mass_medium], df2[!, :calculated_mass_medium])
# 	errorbars!(axbottom, df2[!, :ground_truth_mass_medium], df2[!, :calculated_mass_medium], rms(df2[!, :ground_truth_mass_medium], df2[!, :calculated_mass_medium]))
# 	scatter!(axbottom, df2[!, :ground_truth_mass_small], df2[!, :calculated_mass_small], color=:red)
# 	errorbars!(axbottom, df2[!, :ground_truth_mass_small], df2[!, :calculated_mass_small], rms(df2[!, :ground_truth_mass_small], df2[!, :calculated_mass_small]))
# 	lines!(axbottom, [-1000, 1000], [-1000, 1000])
# 	lines!(axbottom, collect(1:1000), pred_i_low, label = L"y=%$(trunc(co2[2]; digits=3))x + %$(trunc(co2[1]; digits=3))", linestyle=:dashdot)
# 	# axislegend(; position=(0,1))
# 	Textbox(
# 		f[2, 1], 
# 		placeholder = "y = $(trunc(co2[2]; digits=3))x + $(trunc(co2[1]; digits=3)) \nr=$(trunc(r2_2; digits=3)) \nRMSE: $(trunc(rms_values2[1]; digits=3)) \nRMSD: $(trunc(rms_values2[2]; digits=3))", 
# 		tellheight = false,
#         tellwidth = false,
# 		boxcolor=:white,
# 		halign=:left,
# 		valign=:top,
# 		textsize=12
# 	)

# 	xlims!(axbottom, low=0, high=25)
# 	ylims!(axbottom, low=0, high=25)
# 	axbottom.xticks = [0, 5, 10, 15, 20, 25]
# 	axbottom.yticks = [0, 5, 10, 15, 20, 25]
# 	axbottom.xlabel = "Known Mass (mg)"
# 	axbottom.ylabel = "Calculated Mass (mg)"
# 	axbottom.title = "Integrated (Low-Density)"
# 	hidedecorations!(axbottom, ticklabels=false, ticks=false)

# 	##-- C --##
# 	axtopright = Axis(f[1, 2])
	
# 	df3 = df_a_normal
# 	scatter!(axtopright, df3[!, :ground_truth_mass_large], df3[!, :calculated_mass_large])
# 	errorbars!(axtopright, df3[!, :ground_truth_mass_large], df3[!, :calculated_mass_large], rms(df3[!, :ground_truth_mass_large], df3[!, :calculated_mass_large]))
# 	scatter!(axtopright, df3[!, :ground_truth_mass_medium], df3[!, :calculated_mass_medium])
# 	errorbars!(axtopright, df3[!, :ground_truth_mass_medium], df3[!, :calculated_mass_medium], rms(df3[!, :ground_truth_mass_medium], df3[!, :calculated_mass_medium]))
# 	scatter!(axtopright, df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small], color=:red)
# 	errorbars!(axtopright, df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small], rms(df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small]))
# 	lines!(axtopright, [-1000, 1000], [-1000, 1000])
# 	lines!(axtopright, collect(1:1000), pred_a_norm, label = L"y=%$(trunc(co3[2]; digits=3))x + %$(trunc(co3[1]; digits=3))", linestyle=:dashdot)
# 	# axislegend(; position=(0,1))
# 	Textbox(
# 		f[1, 2], 
# 		placeholder = "y = $(trunc(co3[2]; digits=3))x + $(trunc(co3[1]; digits=3)) \nr=$(trunc(r2_3; digits=3)) \nRMSE: $(trunc(rms_values3[1]; digits=3)) \nRMSD: $(trunc(rms_values3[2]; digits=3))", 
# 		tellheight = false,
#         tellwidth = false,
# 		boxcolor=:white,
# 		halign=:left,
# 		valign=:top,
# 		textsize=12
# 	)
	
# 	xlims!(axtopright, low=0, high=200)
# 	ylims!(axtopright, low=0, high=200)
# 	axtopright.xticks = [0, 50, 100, 150, 200]
# 	axtopright.yticks = [0, 50, 100, 150, 200]
# 	axtopright.xlabel = "Known Mass (mg)"
# 	axtopright.ylabel = "Calculated Mass (mg)"
# 	axtopright.title = "Agatston (Normal-Density)"
# 	hidedecorations!(axtopright, ticklabels=false, ticks=false)

# 	##-- D --##
# 	axbottomright = Axis(f[2, 2])
	
# 	df4 = df_a_low
# 	sc1 = scatter!(axbottomright, df4[!, :ground_truth_mass_large], df4[!, :calculated_mass_large])
# 	errorbars!(axbottomright, df4[!, :ground_truth_mass_large], df4[!, :calculated_mass_large], rms(df4[!, :ground_truth_mass_large], df4[!, :calculated_mass_large]))
# 	sc2 = scatter!(axbottomright, df4[!, :ground_truth_mass_medium], df4[!, :calculated_mass_medium])
# 	errorbars!(axbottomright, df4[!, :ground_truth_mass_medium], df4[!, :calculated_mass_medium], rms(df4[!, :ground_truth_mass_medium], df4[!, :calculated_mass_medium]))
# 	sc3 = scatter!(axbottomright, df4[!, :ground_truth_mass_small], df4[!, :calculated_mass_small], color=:red)
# 	errorbars!(axbottomright, df4[!, :ground_truth_mass_small], df4[!, :calculated_mass_small], rms(df4[!, :ground_truth_mass_small], df4[!, :calculated_mass_small]))
# 	ln1 = lines!(axbottomright, [-1000, 1000], [-1000, 1000])
# 	ln2 = lines!(axbottomright, collect(1:1000), pred_a_low, label = L"y=%$(trunc(co4[2]; digits=3))x %$(trunc(co4[1]; digits=3))", linestyle=:dashdot)
# 	# axislegend(; position=(0,1))
# 	Textbox(
# 		f[2, 2], 
# 		placeholder = "y = $(trunc(co4[2]; digits=3))x + $(trunc(co4[1]; digits=3)) \nr=$(trunc(r2_4; digits=3)) \nRMSE: $(trunc(rms_values4[1]; digits=3)) \nRMSD: $(trunc(rms_values4[2]; digits=3))", 
# 		tellheight = false,
#         tellwidth = false,
# 		boxcolor=:white,
# 		halign=:left,
# 		valign=:top,
# 		textsize=12
# 	)

# 	xlims!(axbottomright, low=0, high=25)
# 	ylims!(axbottomright, low=0, high=25)
# 	axbottomright.xticks = [0, 5, 10, 15, 20, 25]
# 	axbottomright.yticks = [0, 5, 10, 15, 20, 25]
# 	axbottomright.xlabel = "Known Mass (mg)"
# 	axbottomright.ylabel = "Calculated Mass (mg)"
# 	axbottomright.title = "Agatston (Low-Density)"
# 	hidedecorations!(axbottomright, ticklabels=false, ticks=false)

# 	##-- LABELS --##

# 	f[1:2, 3] = Legend(f, [sc1, sc2, sc3, ln1, ln2], ["Large Inserts", "Medium Inserts", "Small Inserts", "Unity", "Fitted Line"], framevisible = false)

	
# 	for (label, layout) in zip(["A", "B", "C", "D"], [f[1,1], f[2,1], f[1,2], f[2,2]])
# 	    Label(layout[1, 1, TopLeft()], label,
# 	        textsize = 25,
# 	        padding = (0, 60, 25, 0),
# 	        halign = :right)
# 	end

# 	f
# end

# ╔═╡ 33e21324-b4f8-4b5a-b759-dc38057a612d
md"""
#### Normal Density
"""

# ╔═╡ 1ef8e352-7e89-461f-9ac8-8a39cf2c14d7
let
	df = df_i_normal
	gt_array = vec(hcat(df[!, :ground_truth_mass_large], df[!, :ground_truth_mass_medium], df[!, :ground_truth_mass_small]))
	calc_array = vec(hcat(df[!, :calculated_mass_large], df[!, :calculated_mass_medium], df[!, :calculated_mass_small]))
	data = DataFrame(
		X = gt_array,
		Y= calc_array
	)
	global model_i_normal
	model_i_normal = lm(@formula(Y ~ X), data)
	global r2_1
	r2_1 = GLM.r2(model_i_normal)
	global rms_values1
	rms_values1 = [
		rms(data[!, :X], data[!, :Y]),
		rmsd(data[!, :Y], GLM.predict(model_i_normal))
	]
end

# ╔═╡ ece4a68f-1f3b-4140-bc34-064cd3953763
begin
	newX1 = DataFrame(X=collect(1:1000));
	pred_i_norm = GLM.predict(model_i_normal, newX1)
end

# ╔═╡ 731edfcf-6c6b-4bbc-b23f-74be29d80513
co1 = coef(model_i_normal)

# ╔═╡ c86c3d16-ace0-4de2-a7fb-267ea2925302
let
	df = df_a_normal
	gt_array = vec(hcat(df[!, :ground_truth_mass_large], df[!, :ground_truth_mass_medium], df[!, :ground_truth_mass_small]))
	calc_array = vec(hcat(df[!, :calculated_mass_large], df[!, :calculated_mass_medium], df[!, :calculated_mass_small]))
	data = DataFrame(
		X = gt_array,
		Y= calc_array
	)
	global model_a_normal
	model_a_normal = lm(@formula(Y ~ X), data)
	global r2_3
	r2_3 = GLM.r2(model_a_normal)
	global rms_values3
	rms_values3 = [
		rms(data[!, :X], data[!, :Y]),
		rmsd(data[!, :Y], GLM.predict(model_a_normal))
	]
end

# ╔═╡ 36718970-38e9-4497-b9dd-d9947b1717e4
begin
	newX3 = DataFrame(X=collect(1:1000));
	pred_a_norm = GLM.predict(model_a_normal, newX3)
end

# ╔═╡ e7bb990c-3414-480a-b9f4-1b532bc006c9
co3 = coef(model_a_normal)

# ╔═╡ 42019ced-5c19-4a6f-8e98-c5bdb9e0ce1f
let
	f = Figure()

	##-- A --##
	axtop = Axis(f[1, 1])
	
	df = df_i_normal
	scatter!(axtop, df[!, :ground_truth_mass_large], df[!, :calculated_mass_large])
	errorbars!(axtop, df[!, :ground_truth_mass_large], df[!, :calculated_mass_large], rms(df[!, :ground_truth_mass_large], df[!, :calculated_mass_large]))
	scatter!(axtop, df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium])
	errorbars!(axtop, df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium], rms(df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium]))
	scatter!(axtop, df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], color=:red)
	errorbars!(axtop, df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], rms(df[!, :ground_truth_mass_small], df[!, :calculated_mass_small]))
	lines!(axtop, [-1000, 1000], [-1000, 1000],)
	lines!(axtop, collect(1:1000), pred_i_norm, linestyle=:dashdot)
	Textbox(
		f[1, 1], 
		placeholder = "y = $(trunc(co1[2]; digits=3))x + $(trunc(co1[1]; digits=3)) \nr = $(trunc(r2_1; digits=3)) \nRMSE: $(trunc(rms_values1[1]; digits=3)) \nRMSD: $(trunc(rms_values1[2]; digits=3))", 
		tellheight = false,
        tellwidth = false,
		boxcolor=:white,
		halign=:left,
		valign=:top,
		textsize=12
	)

	xlims!(axtop, low=0, high=200)
	ylims!(axtop, low=0, high=200)
	axtop.xticks = [0, 50, 100, 150, 200]
	axtop.yticks = [0, 50, 100, 150, 200]
	axtop.xlabel = "Known Mass (mg)"
	axtop.ylabel = "Calculated Mass (mg)"
	axtop.title = "Integrated (Normal-Density)"
	hidedecorations!(axtop, ticklabels=false, ticks=false, label=false)
	
	##-- B --##
	axtopright = Axis(f[2, 1])
	
	df3 = df_a_normal
	sc1=scatter!(axtopright, df3[!, :ground_truth_mass_large], df3[!, :calculated_mass_large])
	errorbars!(axtopright, df3[!, :ground_truth_mass_large], df3[!, :calculated_mass_large], rms(df3[!, :ground_truth_mass_large], df3[!, :calculated_mass_large]))
	sc2=scatter!(axtopright, df3[!, :ground_truth_mass_medium], df3[!, :calculated_mass_medium])
	errorbars!(axtopright, df3[!, :ground_truth_mass_medium], df3[!, :calculated_mass_medium], rms(df3[!, :ground_truth_mass_medium], df3[!, :calculated_mass_medium]))
	sc3=scatter!(axtopright, df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small], color=:red)
	errorbars!(axtopright, df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small], rms(df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small]))
	ln1=lines!(axtopright, [-1000, 1000], [-1000, 1000])
	ln2=lines!(axtopright, collect(1:1000), pred_a_norm, linestyle=:dashdot)
	Textbox(
		f[2, 1], 
		placeholder = "y = $(trunc(co3[2]; digits=3))x + $(trunc(co3[1]; digits=3)) \nr = $(trunc(r2_3; digits=3)) \nRMSE: $(trunc(rms_values3[1]; digits=3)) \nRMSD: $(trunc(rms_values3[2]; digits=3))", 
		tellheight = false,
        tellwidth = false,
		boxcolor=:white,
		halign=:left,
		valign=:top,
		textsize=12
	)
	
	xlims!(axtopright, low=0, high=200)
	ylims!(axtopright, low=0, high=200)
	axtopright.xticks = [0, 50, 100, 150, 200]
	axtopright.yticks = [0, 50, 100, 150, 200]
	axtopright.xlabel = "Known Mass (mg)"
	axtopright.ylabel = "Calculated Mass (mg)"
	axtopright.title = "Agatston (Normal-Density)"
	hidedecorations!(axtopright, ticklabels=false, ticks=false, label=false)
	
	##-- LABELS --##

	f[1:2, 2] = Legend(f, [sc1, sc2, sc3, ln1, ln2], ["Large Inserts", "Medium Inserts", "Small Inserts", "Unity", "Fitted Line"], framevisible = false)

	
	for (label, layout) in zip(["A", "B"], [f[1,1], f[2,1]])
	    Label(layout[1, 1, TopLeft()], label,
	        textsize = 25,
	        padding = (0, 60, 25, 0),
	        halign = :right)
	end

	save("/Users/daleblack/Google Drive/Research/2022-AAPM/linear_reg_norm.png", f)
	f
end

# ╔═╡ cf59cd3f-026c-4321-a1ac-de217177b52e
md"""
#### Low Density
"""

# ╔═╡ 0ce89635-b285-4ad7-a936-134eaac981ca
let
	df = df_i_low
	gt_array = vec(hcat(df[!, :ground_truth_mass_large], df[!, :ground_truth_mass_medium], df[!, :ground_truth_mass_small]))
	calc_array = vec(hcat(df[!, :calculated_mass_large], df[!, :calculated_mass_medium], df[!, :calculated_mass_small]))
	data = DataFrame(
		X = gt_array,
		Y= calc_array
	)
	global model_i_low
	model_i_low = lm(@formula(Y ~ X), data)
	global r2_2
	r2_2 = GLM.r2(model_i_low)
	global rms_values2
	rms_values2 = [
		rms(data[!, :X], data[!, :Y]),
		rmsd(data[!, :Y], GLM.predict(model_i_low))
	]
end

# ╔═╡ 6a163966-b24c-4e6a-b2a1-795bbf54bc36
begin
	newX2 = DataFrame(X=collect(1:1000));
	pred_i_low = GLM.predict(model_i_low, newX2)
end

# ╔═╡ 27e191ad-b04c-4e48-9763-128637fdb4b6
co2 = coef(model_i_low)

# ╔═╡ d4c89efd-fb4b-4553-8920-847c393cb4bc
let
	df = df_a_low
	gt_array = vec(hcat(df[!, :ground_truth_mass_large], df[!, :ground_truth_mass_medium], df[!, :ground_truth_mass_small]))
	calc_array = vec(hcat(df[!, :calculated_mass_large], df[!, :calculated_mass_medium], df[!, :calculated_mass_small]))
	data = DataFrame(
		X = gt_array,
		Y= calc_array
	)
	global model_a_low
	model_a_low = lm(@formula(Y ~ X), data)
	global r2_4
	r2_4 = GLM.r2(model_a_low)
	global rms_values4
	rms_values4 = [
		rms(data[!, :X], data[!, :Y]),
		rmsd(data[!, :Y], GLM.predict(model_a_low))
	]
end

# ╔═╡ 5ccbe974-afb9-454f-bd1a-88879a507565
begin
	newX4 = DataFrame(X=collect(1:1000));
	pred_a_low = GLM.predict(model_a_low, newX4)
end

# ╔═╡ d43c412e-22b5-4406-a883-3410641b5b16
co4 = coef(model_a_low)

# ╔═╡ 2c325dab-f159-48e0-8ee2-f0ce2dd5ea6a
let
	f = Figure()
	##-- B --##
	axbottom = Axis(f[1, 1])
	
	df2 = df_i_low
	sc1 = scatter!(axbottom, df2[!, :ground_truth_mass_large], df2[!, :calculated_mass_large])
	errorbars!(axbottom, df2[!, :ground_truth_mass_large], df2[!, :calculated_mass_large], rms(df2[!, :ground_truth_mass_large], df2[!, :calculated_mass_large]))
	sc2 = scatter!(axbottom, df2[!, :ground_truth_mass_medium], df2[!, :calculated_mass_medium])
	errorbars!(axbottom, df2[!, :ground_truth_mass_medium], df2[!, :calculated_mass_medium], rms(df2[!, :ground_truth_mass_medium], df2[!, :calculated_mass_medium]))
	sc3 = scatter!(axbottom, df2[!, :ground_truth_mass_small], df2[!, :calculated_mass_small], color=:red)
	errorbars!(axbottom, df2[!, :ground_truth_mass_small], df2[!, :calculated_mass_small], rms(df2[!, :ground_truth_mass_small], df2[!, :calculated_mass_small]))
	ln1 = lines!(axbottom, [-1000, 1000], [-1000, 1000])
	ln2 = lines!(axbottom, collect(1:1000), pred_i_low, linestyle=:dashdot)
	Textbox(
		f[1, 1], 
		placeholder = "y = $(trunc(co2[2]; digits=3))x + $(trunc(co2[1]; digits=3)) \nr = $(trunc(r2_2; digits=3)) \nRMSE: $(trunc(rms_values2[1]; digits=3)) \nRMSD: $(trunc(rms_values2[2]; digits=3))", 
		tellheight = false,
        tellwidth = false,
		boxcolor=:white,
		halign=:left,
		valign=:top,
		textsize=12
	)

	xlims!(axbottom, low=0, high=25)
	ylims!(axbottom, low=0, high=25)
	axbottom.xticks = [0, 5, 10, 15, 20, 25]
	axbottom.yticks = [0, 5, 10, 15, 20, 25]
	axbottom.xlabel = "Known Mass (mg)"
	axbottom.ylabel = "Calculated Mass (mg)"
	axbottom.title = "Integrated (Low-Density)"
	hidedecorations!(axbottom, ticklabels=false, ticks=false, label=false)
	
	
	##-- D --##
	axbottomright = Axis(f[2, 1])
	
	df4 = df_a_low
	sc1 = scatter!(axbottomright, df4[!, :ground_truth_mass_large], df4[!, :calculated_mass_large])
	errorbars!(axbottomright, df4[!, :ground_truth_mass_large], df4[!, :calculated_mass_large], rms(df4[!, :ground_truth_mass_large], df4[!, :calculated_mass_large]))
	sc2 = scatter!(axbottomright, df4[!, :ground_truth_mass_medium], df4[!, :calculated_mass_medium])
	errorbars!(axbottomright, df4[!, :ground_truth_mass_medium], df4[!, :calculated_mass_medium], rms(df4[!, :ground_truth_mass_medium], df4[!, :calculated_mass_medium]))
	sc3 = scatter!(axbottomright, df4[!, :ground_truth_mass_small], df4[!, :calculated_mass_small], color=:red)
	errorbars!(axbottomright, df4[!, :ground_truth_mass_small], df4[!, :calculated_mass_small], rms(df4[!, :ground_truth_mass_small], df4[!, :calculated_mass_small]))
	ln1 = lines!(axbottomright, [-1000, 1000], [-1000, 1000])
	ln2 = lines!(axbottomright, collect(1:1000), pred_a_low, linestyle=:dashdot)
	Textbox(
		f[2, 1], 
		placeholder = "y = $(trunc(co4[2]; digits=3))x + $(trunc(co4[1]; digits=3)) \nr = $(trunc(r2_4; digits=3)) \nRMSE: $(trunc(rms_values4[1]; digits=3)) \nRMSD: $(trunc(rms_values4[2]; digits=3))", 
		tellheight = false,
        tellwidth = false,
		boxcolor=:white,
		halign=:left,
		valign=:top,
		textsize=12
	)

	xlims!(axbottomright, low=0, high=25)
	ylims!(axbottomright, low=0, high=25)
	axbottomright.xticks = [0, 5, 10, 15, 20, 25]
	axbottomright.yticks = [0, 5, 10, 15, 20, 25]
	axbottomright.xlabel = "Known Mass (mg)"
	axbottomright.ylabel = "Calculated Mass (mg)"
	axbottomright.title = "Agatston (Low-Density)"
	hidedecorations!(axbottomright, ticklabels=false, ticks=false, label=false)

	

	##-- LABELS --##

	f[1:2, 2] = Legend(f, [sc1, sc2, sc3, ln1, ln2], ["Large Inserts", "Medium Inserts", "Small Inserts", "Unity", "Fitted Line"], framevisible = false)

	
	for (label, layout) in zip(["A", "B"], [f[1,1], f[2,1]])
	    Label(layout[1, 1, TopLeft()], label,
	        textsize = 25,
	        padding = (0, 60, 25, 0),
	        halign = :right)
	end

	save("/Users/daleblack/Google Drive/Research/2022-AAPM/linear_reg_low.png", f)
	f
end

# ╔═╡ 69a5a71a-b05f-4c71-bfb2-4d085e48c645
md"""
#### SWCS
"""

# ╔═╡ d98ac8e2-0d07-4194-8bd6-2644a27f9b48
array_s = Array(hcat(df_s[!, end-3:end-1]));

# ╔═╡ 2a790b87-d76e-4a66-8add-4c776b0f1a93
vec_s = vec(array_s);

# ╔═╡ ea7e4439-7c29-4ae9-8038-0dd760cf9545
md"""
#### Agatston
"""

# ╔═╡ 6e0c83a5-45b2-4966-a7e4-721b2b403ec7
vec_a = vec(array_a);

# ╔═╡ 1abbc7c8-46cd-42ad-bf2a-2f2f4f883f79
md"""
#### Integrated
"""

# ╔═╡ 6372b044-a4c1-4a67-888e-1d473f3d6a87
array_i = hcat(df_i[!, 7], df_i[!, 9], df_i[!, 11]);

# ╔═╡ 1bca72fc-adb9-4d54-8e30-8b23ba4da329
vec_i = vec(array_i);

# ╔═╡ 4fde5c38-784d-4461-ab16-525ffb678a42
md"""
#### Known mass
"""

# ╔═╡ 3e9a9103-5d84-44de-aeab-9adce3adb876
array_k = hcat(df_a[!, 6], df_a[!, 8], df_a[!, 10]);

# ╔═╡ 49a9cc2f-bedd-4992-99d0-08f8001c5ba0
vec_k = vec(array_k);

# ╔═╡ ee93ecc4-a964-4b0d-8673-da2b2f1430b7
md"""
#### Correlations
"""

# ╔═╡ b3d14f1b-d538-415a-bb8d-df9abb36be70
pearson_a = cor(hcat(vec_a, vec_k));

# ╔═╡ e7c525f4-e9da-424e-9483-e2d97832e19d
pearson_i = cor(hcat(vec_i, vec_k));

# ╔═╡ 7ddde041-6e01-4f4c-b930-c2b5143fd51f
pearson_s = cor(hcat(vec_s, vec_k));

# ╔═╡ c05610e1-ecc6-4391-b6cf-964bb4696362
p_a = EqualVarianceTTest(vec_a, vec_k);

# ╔═╡ 8fd6011f-24b4-43a8-809e-371ef33ceedc
p_i = EqualVarianceTTest(vec_i, vec_k);

# ╔═╡ e89f1021-7a7a-4d32-82b3-dd93d2823388
p_s = EqualVarianceTTest(vec_s, vec_k);

# ╔═╡ 870ba06e-5f39-4200-a29d-9a771396c9f8
p_corr_coef = [
	(pearson_a[2], "p = ", pvalue(p_a))
	(pearson_i[2], "p = ", pvalue(p_i))
	(pearson_s[2], "p = ", pvalue(p_s))
];

# ╔═╡ 70dcfdd9-b6c2-4b1c-8b6f-78724050acaf
let
	df = df_a
	gt_array = vec(hcat(df[!, :ground_truth_mass_large], df[!, :ground_truth_mass_medium], df[!, :ground_truth_mass_small]))
	calc_array = vec(hcat(df[!, :calculated_mass_large], df[!, :calculated_mass_medium], df[!, :calculated_mass_small]))
	data = DataFrame(
		X = gt_array,
		Y= calc_array
	)
	model = lm(@formula(Y ~ X), data)
	global r_squared_a
	r_squared_a = r2(model)
end;

# ╔═╡ ad15a980-a7a3-4118-9152-6903bf24b605
let
	df = df_i
	gt_array = vec(hcat(df[!, :ground_truth_mass_large], df[!, :ground_truth_mass_medium], df[!, :ground_truth_mass_small]))
	calc_array = vec(hcat(df[!, :calculated_mass_large], df[!, :calculated_mass_medium], df[!, :calculated_mass_small]))
	data = DataFrame(
		X = gt_array,
		Y= calc_array
	)
	model = lm(@formula(Y ~ X), data)
	global r_squared_i
	r_squared_i = r2(model)
end;

# ╔═╡ b11aadd1-6124-4155-9ce1-c2f162053bf9
r_squared = [
	r_squared_a,
	r_squared_i,
	"NA"
]

# ╔═╡ 636cfdaa-4b53-426b-9465-aa348f4bee7f
RMS_values = [
	rms(vec_a, vec_k),
	rms(vec_i, vec_k),
	"NA"
]

# ╔═╡ 4dcbfbdb-d620-438b-b4c5-c288563d63a8
t_test_rmse1 = OneSampleTTest(Float64.(RMS_values[1:2]));

# ╔═╡ 0f59b081-4fcd-43ba-9d79-83af3ae2b3ac
pvalue(t_test_rmse1)

# ╔═╡ 20dc54a6-6351-440f-bc32-9b95e962f970
DataFrame(
	Pearson_Corr_Coeff = p_corr_coef,
	R_Squared = r_squared,
	RMS_values = RMS_values
)

# ╔═╡ 7ec78f71-6b20-4c8c-a9da-a216404bee72
md"""
## Reproducibility (FIG)
"""

# ╔═╡ 387fafb4-4a14-482f-8a8e-4214217f82a5
md"""
#### Integrated (Repeated)
"""

# ╔═╡ 53c1b176-e2f7-4cb9-baa9-26d61ab8c18f
path_integrated_r = "/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/output_repeated/integrated_scoring";

# ╔═╡ 89b28ac5-dd69-4812-84fd-64b54606d146
df_i_r = CSV.read(string(path_integrated_r, "/full.csv"), DataFrame);

# ╔═╡ 4f5328e7-66dd-433c-a7f6-5dc461c83a4c
df_i_low_r, df_i_normal_r = groupby(df_i_r, :DENSITY);

# ╔═╡ 6e515ee8-2836-4285-a0a9-131e1e7b4b7f
df_i_low_small_r, df_i_low_medium_r, df_i_low_large_r = groupby(df_i_low_r, :SIZE);

# ╔═╡ 9ae3b4cb-fac5-4168-868d-724de9b0c9d2
df_i_normal_small_r, df_i_normal_medium_r, df_i_normal_large_r = groupby(df_i_normal_r, :SIZE);

# ╔═╡ 4b82a62c-9f2c-4e61-b4de-f61fa8292dad
let
	df_r = df_i_r
	df = df_i
	r_array = vec(hcat(df_r[!, :calculated_mass_large], df_r[!, :calculated_mass_medium], df_r[!, :calculated_mass_small]))
	calc_array = vec(hcat(df[!, :calculated_mass_large], df[!, :calculated_mass_medium], df[!, :calculated_mass_small]))
	data = DataFrame(
		X = r_array,
		Y= calc_array
	)
	global model_ir
	model_ir = lm(@formula(Y ~ X), data)
	global r2_1r
	r2_1r = GLM.r2(model_ir)
	global rms_values1r
	rms_values1r = [
		rms(data[!, :X], data[!, :Y]),
		rmsd(data[!, :Y], GLM.predict(model_ir))
	]
end

# ╔═╡ 5fb6004b-e9b4-4da5-961b-51170ab8e57e
begin
	newX1r = DataFrame(X=collect(1:1000));
	pred_ir= GLM.predict(model_ir, newX1r)
end

# ╔═╡ 437ae6df-aeb3-457e-ae3a-3c3701e53f34
co1r = coef(model_ir)

# ╔═╡ 5248fac3-4dc2-4579-a52b-8d04ebf45b07
md"""
#### Agatston (Repeated)
"""

# ╔═╡ fcd83805-06e2-4cda-aa56-0c13c69424d8
path_agat_r = "/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/output_repeated/agatston";

# ╔═╡ 5f2eb2d2-84c4-4160-a92a-f181b4126450
df_a_r = CSV.read(string(path_agat_r, "/full.csv"), DataFrame);

# ╔═╡ 869de9cc-6ab0-4e0b-ad1f-59db1fd2f69f
df_a_low_r, df_a_normal_r = groupby(df_a_r, :DENSITY);

# ╔═╡ 43b3ba9b-247c-41cd-968e-3a27736426b0
df_a_low_small_r, df_a_low_medium_r, df_a_low_large_r = groupby(df_a_low_r, :SIZE);

# ╔═╡ b16090b5-a232-4e71-9148-b71296efa999
df_a_normal_small_r, df_a_normal_medium_r, df_a_normal_large_r= groupby(df_a_normal_r, :SIZE);

# ╔═╡ a3996a98-d95f-4a51-8776-d6286f5c8ed5
let
	df_r = df_a_r
	df = df_a
	r_array = vec(hcat(df_r[!, :calculated_mass_large], df_r[!, :calculated_mass_medium], df_r[!, :calculated_mass_small]))
	calc_array = vec(hcat(df[!, :calculated_mass_large], df[!, :calculated_mass_medium], df[!, :calculated_mass_small]))
	data = DataFrame(
		X = r_array,
		Y= calc_array
	)
	global model_ar
	model_ar = lm(@formula(Y ~ X), data)
	global r2_2r
	r2_2r = GLM.r2(model_ar)
	global rms_values2r
	rms_values2r = [
		rms(data[!, :X], data[!, :Y]),
		rmsd(data[!, :Y], GLM.predict(model_ar))
	]
end

# ╔═╡ 664e7901-31bc-4c49-8d6a-e770e0553c31
begin
	newX2r = DataFrame(X=collect(1:1000));
	pred_ar= GLM.predict(model_ar, newX2r)
end

# ╔═╡ f278a992-1527-466c-909b-8ff83fed0a61
co2r = coef(model_ar)

# ╔═╡ c1d820ac-d0ca-4c18-b61d-849402ae647e
md"""
#### SWCS (Repeated)
"""

# ╔═╡ 594b8d0e-f65c-4fe3-b571-e2ee68a9ab5f
path_swcs_r = "/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/output_repeated/swcs";

# ╔═╡ 1d208fd8-b847-4afe-84b5-3a553f50a858
df_s_r = CSV.read(string(path_swcs_r, "/full.csv"), DataFrame);

# ╔═╡ a9f68d7e-7818-47d5-ba1f-965634225b30
df_s_low_r, df_s_normal_r = groupby(df_s_r, :DENSITY);

# ╔═╡ f2177e1f-2a55-4f81-9a69-ab985ae3b7c2
df_s_low_small_r, df_s_low_medium_r, df_s_low_large_r = groupby(df_s_low_r, :SIZE);

# ╔═╡ 35bc90e2-9d70-4daf-a825-b462831c5bf6
df_s_normal_small_r, df_s_normal_medium_r, df_s_normal_large_r = groupby(df_s_normal_r, :SIZE);

# ╔═╡ a385671b-af3a-480b-bb0c-18d1c13cf943
let
	df_r = df_s_r
	df = df_s
	r_array = vec(hcat(df_r[!, :calculated_swcs_large], df_r[!, :calculated_swcs_medium], df_r[!, :calculated_swcs_small]))
	calc_array = vec(hcat(df[!, :calculated_swcs_large], df[!, :calculated_swcs_medium], df[!, :calculated_swcs_small]))
	data = DataFrame(
		X = r_array,
		Y= calc_array
	)
	global model_sr
	model_sr = lm(@formula(Y ~ X), data)
	global r2_3r
	r2_3r = GLM.r2(model_sr)
	global rms_values3r
	rms_values3r = [
		rms(data[!, :X], data[!, :Y]),
		rmsd(data[!, :Y], GLM.predict(model_sr))
	]
end

# ╔═╡ fd610b17-74b8-44e1-bcac-5ac910db05e5
begin
	newX3r = DataFrame(X=collect(1:1000));
	pred_sr= GLM.predict(model_sr, newX3r)
end

# ╔═╡ 8f6d3806-cc69-4eb2-97bb-cf0cd0e0b565
co3r = coef(model_sr)

# ╔═╡ a85f777a-76f2-4c64-9973-ea9dec245600
let
	f = Figure()

	##-- A --##
	axtop = Axis(f[1, 1])
	
	df = df_i
	df_r = df_i_r
	scatter!(axtop, df_r[!, :calculated_mass_large], df[!, :calculated_mass_large])
	scatter!(axtop, df_r[!, :calculated_mass_medium], df[!, :calculated_mass_medium])
	scatter!(axtop, df_r[!, :calculated_mass_small], df[!, :calculated_mass_small], color=:red)
	lines!(axtop, [-1000, 1000], [-1000, 1000], label="Unity")
	lines!(axtop, collect(1:1000), pred_ir, linestyle=:dashdot)
	Textbox(
		f[1, 1], 
		placeholder = "y = $(trunc(co1r[2]; digits=3))x + $(trunc(co1r[1]; digits=3)) \nr = $(trunc(r2_1r; digits=3)) \nRMSE: $(trunc(rms_values1r[1]; digits=3)) \nRMSD: $(trunc(rms_values1r[2]; digits=3))", 
		tellheight = false,
        tellwidth = false,
		boxcolor=:white,
		halign=:left,
		valign=:top,
		textsize=12
	)

	xlims!(axtop, low=0, high=200)
	ylims!(axtop, low=0, high=200)
	axtop.xticks = [0, 50, 100, 150, 200]
	axtop.yticks = [0, 50, 100, 150, 200]
	axtop.xlabel = "Mass (mg)"
	axtop.ylabel = "Mass (mg)"
	axtop.title = "Integrated"
	hidedecorations!(axtop, ticklabels=false, ticks=false, label=false)

	##-- B --##
	axtopright = Axis(f[2, 1])
	
	df3 = df_a
	df3_r = df_a_r
	scatter!(axtopright, df3_r[!, :calculated_mass_large], df3[!, :calculated_mass_large], label="Large Inserts")
	scatter!(axtopright, df3_r[!, :calculated_mass_medium], df3[!, :calculated_mass_medium], label="Medium Inserts")
	scatter!(axtopright, df3_r[!, :calculated_mass_small], df3[!, :calculated_mass_small], label="Small Inserts", color=:red)
	lines!(axtopright, [-1000, 1000], [-1000, 1000], label="Unity")
	lines!(axtopright, collect(1:1000), pred_ar, linestyle=:dashdot)
	Textbox(
		f[2, 1], 
		placeholder = "y = $(trunc(co2r[2]; digits=3))x + $(trunc(co2r[1]; digits=3)) \nr = $(trunc(r2_2r; digits=3)) \nRMSE: $(trunc(rms_values2r[1]; digits=3)) \nRMSD: $(trunc(rms_values2r[2]; digits=3))", 
		tellheight = false,
        tellwidth = false,
		boxcolor=:white,
		halign=:left,
		valign=:top,
		textsize=12
	)

	xlims!(axtopright, low=0, high=200)
	ylims!(axtopright, low=0, high=200)
	axtopright.xticks = [0, 50, 100, 150, 200]
	axtopright.yticks = [0, 50, 100, 150, 200]
	axtopright.xlabel = "Mass (mg)"
	axtopright.ylabel = "Mass (mg)"
	axtopright.title = "Agatston"
	hidedecorations!(axtopright, ticklabels=false, ticks=false, label=false)

	# ##-- C --##
	axbottomright = Axis(f[3, 1])
	
	df4 = df_s
	df4_r = df_s_r
	scatter!(axbottomright, df4_r[!, :calculated_swcs_large], df4[!, :calculated_swcs_large], label="Large Inserts")
	scatter!(axbottomright, df4_r[!, :calculated_swcs_medium], df4[!, :calculated_swcs_medium], label="Medium Inserts")
	scatter!(axbottomright, df4_r[!, :calculated_swcs_small], df4[!, :calculated_swcs_small], label="Small Inserts", color=:red)
	lines!(axbottomright, [-1000, 1000], [-1000, 1000], label="Unity")
	lines!(axbottomright, collect(1:1000), pred_ar, linestyle=:dashdot)
	Textbox(
		f[3, 1], 
		placeholder = "y = $(trunc(co3r[2]; digits=3))x + $(trunc(co3r[1]; digits=3)) \nr = $(trunc(r2_3r; digits=3)) \nRMSE: $(trunc(rms_values3r[1]; digits=3)) \nRMSD: $(trunc(rms_values3r[2]; digits=3))", 
		tellheight = false,
        tellwidth = false,
		boxcolor=:white,
		halign=:left,
		valign=:top,
		textsize=12
	)

	xlims!(axbottomright, low=0, high=500)
	ylims!(axbottomright, low=0, high=500)
	axbottomright.xlabel = "SWCS"
	axbottomright.ylabel = "SWCS"
	axbottomright.title = "Spatially Weighted"
	hidedecorations!(axbottomright, ticklabels=false, ticks=false, label=false)

	##-- LABELS --##
	f[2, 2] = Legend(f, axbottomright, framevisible = false)

	
	for (label, layout) in zip(["A", "B", "C"], [f[1,1], f[2,1], f[3,1]])
	    Label(layout[1, 1, TopLeft()], label,
	        textsize = 25,
	        padding = (0, -10, 5, 0),
	        halign = :right)
	end
	
	save("/Users/daleblack/Google Drive/Research/2022-AAPM/reprod.png", f)
	f
end

# ╔═╡ 3e361d54-19f3-48d0-8e85-7351ec2f3335
md"""
## Motion Blurring (FIG)
"""

# ╔═╡ f8c2e340-8aa4-4af6-b4c8-43cb33bef721
md"""
#### Integrated
"""

# ╔═╡ 620ad304-cb77-4a7c-afc7-76a49d220cfd
df_ii = CSV.read(string(path_integrated, "/full.csv"), DataFrame);

# ╔═╡ 663a9cab-3235-486b-97c2-23a4a3c86c10
df_ii0, df_ii05, df_ii1, df_ii15, df_ii2 = groupby(df_ii, :blur);

# ╔═╡ 909f4f32-cfe6-4c82-bfa0-cea38a8de5ca
let
	df = df_ii
	gt_array = vec(hcat(df[!, :ground_truth_mass_large], df[!, :ground_truth_mass_medium], df[!, :ground_truth_mass_small]))
	calc_array = vec(hcat(df[!, :calculated_mass_large], df[!, :calculated_mass_medium], df[!, :calculated_mass_small]))
	data = DataFrame(
		X = gt_array,
		Y= calc_array
	)
	global model_ii
	model_ii = lm(@formula(Y ~ X), data)
	global r2_ii
	r2_ii = GLM.r2(model_ii)
	global rms_valuesii
	rms_valuesii = [
		rms(data[!, :X], data[!, :Y]),
		rmsd(data[!, :Y], GLM.predict(model_ii))
	]
end

# ╔═╡ 9ffe97f5-624d-459e-a064-be4e45123b85
begin
	newXii = DataFrame(X=collect(1:1000));
	pred_ii = GLM.predict(model_ii, newXii)
end

# ╔═╡ 67f78b69-7202-4128-8918-ee65fac779bb
coii = coef(model_ii)

# ╔═╡ 56efa0dd-c7dc-4364-a840-26ceb052826c
md"""
#### Agatston
"""

# ╔═╡ bd7e30d2-916e-49b5-905e-7caa2d9080b4
df_aa = CSV.read(string(path_agat, "/full2.csv"), DataFrame);

# ╔═╡ 8e381b33-f2f5-476f-9a9f-d0b7ec1bb774
let
	df = df_aa
	gt_array = vec(hcat(df[!, :ground_truth_mass_large], df[!, :ground_truth_mass_medium], df[!, :ground_truth_mass_small]))
	calc_array = vec(hcat(df[!, :calculated_mass_large], df[!, :calculated_mass_medium], df[!, :calculated_mass_small]))
	data = DataFrame(
		X = gt_array,
		Y= calc_array
	)
	global model_aa
	model_aa = lm(@formula(Y ~ X), data)
	global r2_aa
	r2_aa = GLM.r2(model_aa)
	global rms_valuesaa
	rms_valuesaa = [
		rms(data[!, :X], data[!, :Y]),
		rmsd(data[!, :Y], GLM.predict(model_aa))
	]
end

# ╔═╡ 41946bb8-fe10-4086-831d-3a68878fe944
begin
	newXaa = DataFrame(X=collect(1:1000));
	pred_aa = GLM.predict(model_aa, newXaa)
end

# ╔═╡ ca8f3808-2173-4a6e-b248-bafc0473bbdb
coaa = coef(model_aa)

# ╔═╡ f00cfcc7-9ecc-4a84-a23d-2a30845be8f8
let
	f = Figure()

	##-- A --##
	axtop = Axis(f[1, 1])
	
	df = df_ii
	scatter!(axtop, df[!, :ground_truth_mass_large], df[!, :calculated_mass_large])
	errorbars!(axtop, df[!, :ground_truth_mass_large], df[!, :calculated_mass_large], rms(df[!, :ground_truth_mass_large], df[!, :calculated_mass_large]))
	scatter!(axtop, df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium])
	errorbars!(axtop, df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium], rms(df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium]))
	scatter!(axtop, df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], color=:red)
	errorbars!(axtop, df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], rms(df[!, :ground_truth_mass_small], df[!, :calculated_mass_small]))
	lines!(axtop, [-1000, 1000], [-1000, 1000],)
	lines!(axtop, collect(1:1000), pred_ii, linestyle=:dashdot)
	Textbox(
		f[1, 1], 
		placeholder = "y = $(trunc(coii[2]; digits=3))x + $(trunc(coii[1]; digits=3)) \nr = $(trunc(r2_ii; digits=3)) \nRMSE: $(trunc(rms_valuesii[1]; digits=3)) \nRMSD: $(trunc(rms_valuesii[2]; digits=3))", 
		tellheight = false,
        tellwidth = false,
		boxcolor=:white,
		halign=:left,
		valign=:top,
		textsize=12
	)

	xlims!(axtop, low=0, high=200)
	ylims!(axtop, low=0, high=200)
	axtop.xticks = [0, 50, 100, 150, 200]
	axtop.yticks = [0, 50, 100, 150, 200]
	axtop.xlabel = "Known Mass (mg)"
	axtop.ylabel = "Calculated Mass (mg)"
	axtop.title = "Integrated"
	hidedecorations!(axtop, ticklabels=false, ticks=false, label=false)
	
	##-- B --##
	axtopright = Axis(f[2, 1])
	
	df3 = df_aa
	sc1=scatter!(axtopright, df3[!, :ground_truth_mass_large], df3[!, :calculated_mass_large])
	errorbars!(axtopright, df3[!, :ground_truth_mass_large], df3[!, :calculated_mass_large], rms(df3[!, :ground_truth_mass_large], df3[!, :calculated_mass_large]))
	sc2=scatter!(axtopright, df3[!, :ground_truth_mass_medium], df3[!, :calculated_mass_medium])
	errorbars!(axtopright, df3[!, :ground_truth_mass_medium], df3[!, :calculated_mass_medium], rms(df3[!, :ground_truth_mass_medium], df3[!, :calculated_mass_medium]))
	sc3=scatter!(axtopright, df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small], color=:red)
	errorbars!(axtopright, df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small], rms(df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small]))
	ln1=lines!(axtopright, [-1000, 1000], [-1000, 1000])
	ln2=lines!(axtopright, collect(1:1000), pred_aa, linestyle=:dashdot)
	Textbox(
		f[2, 1], 
		placeholder = "y = $(trunc(coaa[2]; digits=3))x + $(trunc(coaa[1]; digits=3)) \nr = $(trunc(r2_aa; digits=3)) \nRMSE: $(trunc(rms_valuesaa[1]; digits=3)) \nRMSD: $(trunc(rms_valuesaa[2]; digits=3))", 
		tellheight = false,
        tellwidth = false,
		boxcolor=:white,
		halign=:left,
		valign=:top,
		textsize=12
	)
	
	xlims!(axtopright, low=0, high=200)
	ylims!(axtopright, low=0, high=200)
	axtopright.xticks = [0, 50, 100, 150, 200]
	axtopright.yticks = [0, 50, 100, 150, 200]
	axtopright.xlabel = "Known Mass (mg)"
	axtopright.ylabel = "Calculated Mass (mg)"
	axtopright.title = "Agatston"
	hidedecorations!(axtopright, ticklabels=false, ticks=false, label=false)
	
	##-- LABELS --##

	f[1:2, 2] = Legend(f, [sc1, sc2, sc3, ln1, ln2], ["Large Inserts", "Medium Inserts", "Small Inserts", "Unity", "Fitted Line"], framevisible = false)

	
	for (label, layout) in zip(["A", "B"], [f[1,1], f[2,1]])
	    Label(layout[1, 1, TopLeft()], label,
	        textsize = 25,
	        padding = (0, 60, 25, 0),
	        halign = :right)
	end

	# save("/Users/daleblack/Google Drive/Research/2022-AAPM/linear_reg_norm.png", f)
	f
end

# ╔═╡ bd115bb7-065c-41af-b332-4b701bd6120d
df_aa0, df_aa05, df_aa1, df_aa15, df_aa2 = groupby(df_aa, :blur);

# ╔═╡ Cell order:
# ╠═c82bb3e6-e052-11ec-3392-eba4a4b464ba
# ╠═2df207ed-cf9c-4d7d-8354-4da14f93276c
# ╟─4cf61487-3913-4e11-970e-4ca31a5ffc8d
# ╠═57d46998-368a-4916-8380-ee49d5473a49
# ╠═2c960bd8-ae64-453f-b29f-275bf5263774
# ╠═7f4fae09-7916-4a98-a102-7f861900c457
# ╠═c987e55d-9ffd-40a5-9d64-1653d4000837
# ╠═0c367cb1-d959-4c10-b31c-18d3cb396255
# ╟─a58beea4-b186-44c6-bf75-98d368f7ebe4
# ╠═4fd5b649-cea1-4d1f-95ac-35aa12f8fbff
# ╠═3fa3748a-22d6-49d8-9888-a749dca99ce2
# ╠═c38b9994-7698-4c81-b65b-534b994d8ff0
# ╠═79952ee5-3ccd-4d85-95b4-f7610c90d63d
# ╠═7223ae32-7241-463f-be13-2d2ce6bc7de7
# ╟─f39a9ed5-c897-4b90-8494-bf231312970e
# ╠═32623778-9a19-4318-81c6-1fc7aa0157fe
# ╠═317f8c26-7b39-45c3-b9d3-02925d4b0514
# ╠═06deefa7-9170-4e2f-ac49-d6dc65c0af76
# ╠═5a4cc121-fd02-49bf-bdf5-07b71b48ce19
# ╠═d5d8e9da-eba8-4788-b787-95b3a2b41329
# ╟─310c1999-77a6-432e-850c-4111411cceb0
# ╟─3a8c6f3c-0fa6-477c-85e2-7739cf2537ae
# ╟─b702760a-9d9b-4d2a-992c-37ed4f8bcfdf
# ╟─932e1a11-ec52-48bc-af27-bb3e09e6b27f
# ╠═dc8805a5-2b4c-47fe-88c8-0f4434ddd5cf
# ╠═2c78dcc4-03c7-49b5-ac45-30c198671761
# ╠═d3e13dcd-201a-49ab-86f7-a7d3320088c2
# ╠═4b97d892-b0db-4f75-9239-9da215f905f0
# ╟─01d903d8-fca2-4317-a2d7-93429f06d47c
# ╠═97a612d5-44e7-43f1-b61e-8518c7e6aa28
# ╠═a31dfad2-e420-4c41-b805-227778a39fc9
# ╟─86d5f7fd-6a1b-4fdc-9e06-6cd66f8a169c
# ╠═aabe8e7b-faed-499d-aea2-fd44836a3cfb
# ╠═dd644657-51e9-4b27-a186-6efba0bb83f8
# ╠═d79d6478-0837-40f0-a7a5-270bf44d6166
# ╠═2b5a8dca-faae-4163-ab14-ae4f98b2a7cc
# ╠═f4b51729-d53a-40e2-919b-5d4562cffbc0
# ╠═da2b55dc-2fe9-468f-a53e-baf232ca79cd
# ╟─7dfc24a4-e006-45f4-b5b9-977a7c3c0b7c
# ╟─42019ced-5c19-4a6f-8e98-c5bdb9e0ce1f
# ╟─2c325dab-f159-48e0-8ee2-f0ce2dd5ea6a
# ╟─d04ce8a0-8224-49c7-a057-f1c7d0d65ef8
# ╟─33e21324-b4f8-4b5a-b759-dc38057a612d
# ╠═1ef8e352-7e89-461f-9ac8-8a39cf2c14d7
# ╠═ece4a68f-1f3b-4140-bc34-064cd3953763
# ╠═731edfcf-6c6b-4bbc-b23f-74be29d80513
# ╠═c86c3d16-ace0-4de2-a7fb-267ea2925302
# ╠═36718970-38e9-4497-b9dd-d9947b1717e4
# ╠═e7bb990c-3414-480a-b9f4-1b532bc006c9
# ╟─cf59cd3f-026c-4321-a1ac-de217177b52e
# ╠═0ce89635-b285-4ad7-a936-134eaac981ca
# ╠═6a163966-b24c-4e6a-b2a1-795bbf54bc36
# ╠═27e191ad-b04c-4e48-9763-128637fdb4b6
# ╠═d4c89efd-fb4b-4553-8920-847c393cb4bc
# ╠═5ccbe974-afb9-454f-bd1a-88879a507565
# ╠═d43c412e-22b5-4406-a883-3410641b5b16
# ╟─69a5a71a-b05f-4c71-bfb2-4d085e48c645
# ╠═d98ac8e2-0d07-4194-8bd6-2644a27f9b48
# ╠═2a790b87-d76e-4a66-8add-4c776b0f1a93
# ╟─ea7e4439-7c29-4ae9-8038-0dd760cf9545
# ╠═6e0c83a5-45b2-4966-a7e4-721b2b403ec7
# ╟─1abbc7c8-46cd-42ad-bf2a-2f2f4f883f79
# ╠═6372b044-a4c1-4a67-888e-1d473f3d6a87
# ╠═1bca72fc-adb9-4d54-8e30-8b23ba4da329
# ╟─4fde5c38-784d-4461-ab16-525ffb678a42
# ╠═3e9a9103-5d84-44de-aeab-9adce3adb876
# ╠═49a9cc2f-bedd-4992-99d0-08f8001c5ba0
# ╟─ee93ecc4-a964-4b0d-8673-da2b2f1430b7
# ╠═b3d14f1b-d538-415a-bb8d-df9abb36be70
# ╠═e7c525f4-e9da-424e-9483-e2d97832e19d
# ╠═7ddde041-6e01-4f4c-b930-c2b5143fd51f
# ╠═c05610e1-ecc6-4391-b6cf-964bb4696362
# ╠═8fd6011f-24b4-43a8-809e-371ef33ceedc
# ╠═e89f1021-7a7a-4d32-82b3-dd93d2823388
# ╠═870ba06e-5f39-4200-a29d-9a771396c9f8
# ╠═70dcfdd9-b6c2-4b1c-8b6f-78724050acaf
# ╠═ad15a980-a7a3-4118-9152-6903bf24b605
# ╠═b11aadd1-6124-4155-9ce1-c2f162053bf9
# ╠═636cfdaa-4b53-426b-9465-aa348f4bee7f
# ╠═4dcbfbdb-d620-438b-b4c5-c288563d63a8
# ╠═0f59b081-4fcd-43ba-9d79-83af3ae2b3ac
# ╠═20dc54a6-6351-440f-bc32-9b95e962f970
# ╟─7ec78f71-6b20-4c8c-a9da-a216404bee72
# ╟─a85f777a-76f2-4c64-9973-ea9dec245600
# ╟─387fafb4-4a14-482f-8a8e-4214217f82a5
# ╠═53c1b176-e2f7-4cb9-baa9-26d61ab8c18f
# ╠═89b28ac5-dd69-4812-84fd-64b54606d146
# ╠═4f5328e7-66dd-433c-a7f6-5dc461c83a4c
# ╠═6e515ee8-2836-4285-a0a9-131e1e7b4b7f
# ╠═9ae3b4cb-fac5-4168-868d-724de9b0c9d2
# ╠═4b82a62c-9f2c-4e61-b4de-f61fa8292dad
# ╠═5fb6004b-e9b4-4da5-961b-51170ab8e57e
# ╠═437ae6df-aeb3-457e-ae3a-3c3701e53f34
# ╟─5248fac3-4dc2-4579-a52b-8d04ebf45b07
# ╠═fcd83805-06e2-4cda-aa56-0c13c69424d8
# ╠═5f2eb2d2-84c4-4160-a92a-f181b4126450
# ╠═869de9cc-6ab0-4e0b-ad1f-59db1fd2f69f
# ╠═43b3ba9b-247c-41cd-968e-3a27736426b0
# ╠═b16090b5-a232-4e71-9148-b71296efa999
# ╠═a3996a98-d95f-4a51-8776-d6286f5c8ed5
# ╠═664e7901-31bc-4c49-8d6a-e770e0553c31
# ╠═f278a992-1527-466c-909b-8ff83fed0a61
# ╟─c1d820ac-d0ca-4c18-b61d-849402ae647e
# ╠═594b8d0e-f65c-4fe3-b571-e2ee68a9ab5f
# ╠═1d208fd8-b847-4afe-84b5-3a553f50a858
# ╠═a9f68d7e-7818-47d5-ba1f-965634225b30
# ╠═f2177e1f-2a55-4f81-9a69-ab985ae3b7c2
# ╠═35bc90e2-9d70-4daf-a825-b462831c5bf6
# ╠═a385671b-af3a-480b-bb0c-18d1c13cf943
# ╠═fd610b17-74b8-44e1-bcac-5ac910db05e5
# ╠═8f6d3806-cc69-4eb2-97bb-cf0cd0e0b565
# ╟─3e361d54-19f3-48d0-8e85-7351ec2f3335
# ╟─f00cfcc7-9ecc-4a84-a23d-2a30845be8f8
# ╟─f8c2e340-8aa4-4af6-b4c8-43cb33bef721
# ╠═620ad304-cb77-4a7c-afc7-76a49d220cfd
# ╠═663a9cab-3235-486b-97c2-23a4a3c86c10
# ╠═909f4f32-cfe6-4c82-bfa0-cea38a8de5ca
# ╠═9ffe97f5-624d-459e-a064-be4e45123b85
# ╠═67f78b69-7202-4128-8918-ee65fac779bb
# ╟─56efa0dd-c7dc-4364-a840-26ceb052826c
# ╠═bd7e30d2-916e-49b5-905e-7caa2d9080b4
# ╠═8e381b33-f2f5-476f-9a9f-d0b7ec1bb774
# ╠═41946bb8-fe10-4086-831d-3a68878fe944
# ╠═ca8f3808-2173-4a6e-b248-bafc0473bbdb
# ╠═bd115bb7-065c-41af-b332-4b701bd6120d

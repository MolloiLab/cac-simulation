### A Pluto.jl notebook ###
# v0.19.4

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

# ╔═╡ 5144e2a2-91e2-49ab-9dcf-df57ed767a36
df_i

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
## Figure 1. Zero CAC Scores
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
end

# ╔═╡ d3e13dcd-201a-49ab-86f7-a7d3320088c2
begin
	array_s_80 = Array(df_s_80[!, 12:14])
	array_s_100 = Array(df_s_100[!, 12:14])
	array_s_120 = Array(df_s_120[!, 12:14])
	array_s_135 = Array(df_s_135[!, 12:14])
end;

# ╔═╡ 4b97d892-b0db-4f75-9239-9da215f905f0
begin
	num_zeroCAC_80 = length(findall(x -> x < mean_swcs_80, array_s_80))
	num_zeroCAC_100 = length(findall(x -> x < mean_swcs_100, array_s_100))
	num_zeroCAC_120 = length(findall(x -> x < mean_swcs_120, array_s_120))
	num_zeroCAC_135 = length(findall(x -> x < mean_swcs_135, array_s_135))

	total_zero_s = num_zeroCAC_80 + num_zeroCAC_100 + num_zeroCAC_120 + num_zeroCAC_135
	total_cac = length(array_s_80) * 4
end;

# ╔═╡ f3cfe6db-d3f5-43c9-988e-25f332a5998a
total_cac

# ╔═╡ 8e4fc2cf-16fd-43f2-84f1-8339b6181baf
total_zero_s

# ╔═╡ 01d903d8-fca2-4317-a2d7-93429f06d47c
md"""
#### Agatston
"""

# ╔═╡ 97a612d5-44e7-43f1-b61e-8518c7e6aa28
array_a = hcat(df_a[!, 7], df_a[!, 9], df_a[!, 11]);

# ╔═╡ a31dfad2-e420-4c41-b805-227778a39fc9
num_zero_a = length(findall(x -> x == 0, array_a))

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
	num_zeroCAC_80_i = length(findall(x -> x < mean_i_80, array_i_80))
	num_zeroCAC_100_i = length(findall(x -> x < mean_i_100, array_i_100))
	num_zeroCAC_120_i = length(findall(x -> x < mean_i_120, array_i_120))
	num_zeroCAC_135_i = length(findall(x -> x < mean_i_135, array_i_135))

	total_zero_i = num_zeroCAC_80_i + num_zeroCAC_100_i + num_zeroCAC_120_i + num_zeroCAC_135_i
end;

# ╔═╡ f4b51729-d53a-40e2-919b-5d4562cffbc0
total_zero_i, total_zero_s, num_zero_a

# ╔═╡ b702760a-9d9b-4d2a-992c-37ed4f8bcfdf
let
	f = Figure()
	ga = f[1, 1] = GridLayout()
	# gb = f[2, 1] = GridLayout()

	colors = Makie.wong_colors()
	labels = ["Number Zero CAC (SWCS)", "Number Zero CAC (Agatston)"]
	elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]

	##-- TOP --##
	axtop = Axis(ga[1, 1], xticks = (1:3, ["Integrated", "Spatially Weighted", "Agatston"]))
	# hidexdecorations!(axtop, grid = false)
	# axtop.ylabel = ""
	ylims!(axtop, low=0, high=100)
	# axtop.xticks = [0, 25, 50, 100]
	axtop.yticks = [0, 25, 50, 75, 100]
	table = [1, 2, 3]
	heights1 = [(total_zero_i / total_cac) * 100, (total_zero_s / total_cac) * 100, (num_zero_a / total_cac) * 100]
	barplot!(axtop, table, heights1, color=colors[1:3], bar_labels=:y)
	
	Label(ga[1:end, 1:end, Top()], "Zero CAC Scores", valign = :center, padding = (0, 0, 0, 0), textsize=25)
	Label(ga[1:end, 1:end, Left()], "% Zero CAC Scores", valign = :center, padding = (0, 50, 0, 0), rotation=π/2, textsize=17)

	# for (label, layout) in zip(["A"], [ga])
	#     Label(layout[1, 1, TopLeft()], label,
	#         textsize = 25,
	#         padding = (0, 60, 25, 0),
	#         halign = :right)
	# end

	# Legend(f[1, 2], elements, labels, framevisible = false)
	
	f
end

# ╔═╡ b1a054d6-863f-4675-bb59-01196565bee6
md"""
## Figure 2. Pearson Correlation Coefficient
"""

# ╔═╡ 69a5a71a-b05f-4c71-bfb2-4d085e48c645
md"""
#### SWCS
"""

# ╔═╡ d98ac8e2-0d07-4194-8bd6-2644a27f9b48
array_s = Array(hcat(df_s[!, end-3:end-1]));

# ╔═╡ 2a790b87-d76e-4a66-8add-4c776b0f1a93
vec_s = vec(array_s)

# ╔═╡ ea7e4439-7c29-4ae9-8038-0dd760cf9545
md"""
#### Agatston
"""

# ╔═╡ 6e0c83a5-45b2-4966-a7e4-721b2b403ec7
vec_a = vec(array_a)

# ╔═╡ 1abbc7c8-46cd-42ad-bf2a-2f2f4f883f79
md"""
#### Integrated
"""

# ╔═╡ 6372b044-a4c1-4a67-888e-1d473f3d6a87
array_i = hcat(df_i[!, 7], df_i[!, 9], df_i[!, 11]);

# ╔═╡ 1bca72fc-adb9-4d54-8e30-8b23ba4da329
vec_i = vec(array_i)

# ╔═╡ 4fde5c38-784d-4461-ab16-525ffb678a42
md"""
#### Known mass
"""

# ╔═╡ 3e9a9103-5d84-44de-aeab-9adce3adb876
array_k = hcat(df_a[!, 6], df_a[!, 8], df_a[!, 10]);

# ╔═╡ 49a9cc2f-bedd-4992-99d0-08f8001c5ba0
vec_k = vec(array_k)

# ╔═╡ ee93ecc4-a964-4b0d-8673-da2b2f1430b7
md"""
#### Correlations
"""

# ╔═╡ 7ddde041-6e01-4f4c-b930-c2b5143fd51f
pearson_s = cor(hcat(vec_s, vec_k))

# ╔═╡ b3d14f1b-d538-415a-bb8d-df9abb36be70
pearson_a = cor(hcat(vec_a, vec_k))

# ╔═╡ e7c525f4-e9da-424e-9483-e2d97832e19d
pearson_i = cor(hcat(vec_i, vec_k))

# ╔═╡ e89f1021-7a7a-4d32-82b3-dd93d2823388
p_s = EqualVarianceTTest(vec_s, vec_k)

# ╔═╡ 8fd6011f-24b4-43a8-809e-371ef33ceedc
p_i = EqualVarianceTTest(vec_i, vec_k)

# ╔═╡ c05610e1-ecc6-4391-b6cf-964bb4696362
p_a = EqualVarianceTTest(vec_a, vec_k)

# ╔═╡ 3b3269ef-1e0a-4335-b0a1-e14223670682
let
	f = Figure()
	ga = f[1, 1] = GridLayout()

	colors = Makie.wong_colors()
	labels = ["Number Zero CAC (SWCS)", "Number Zero CAC (Agatston)"]
	elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]

	##-- TOP --##
	axtop = Axis(ga[1, 1], xticks = (1:3, ["Integrated", "Spatially Weighted", "Agatston"]))
	ylims!(axtop, low=0, high=1.1)
	axtop.yticks = [0, 0.25, 0.50, 0.75, 1.0]
	table = [1, 2, 3]
	heights1 = [pearson_i[1, 2], pearson_s[1, 2], pearson_a[1, 2]]
	barplot!(axtop, table, heights1, color=colors[1:3], bar_labels=:y)
	
	Label(ga[1:end, 1:end, Top()], "Correlation to Known Calcium Mass", valign = :center, padding = (0, 0, 0, 0), textsize=25)
	Label(ga[1:end, 1:end, Left()], "Pearson Correlation", valign = :center, padding = (0, 50, 0, 0), rotation=π/2, textsize=17)
	
	f
end

# ╔═╡ 7dfc24a4-e006-45f4-b5b9-977a7c3c0b7c
md"""
## Figure 3. Linear Regression Mass Score
"""

# ╔═╡ d04ce8a0-8224-49c7-a057-f1c7d0d65ef8
let
	f = Figure()
	
	ga = f[1, 1] = GridLayout()
	gb = f[2, 1] = GridLayout()
	gcd = f[1:2, 2] = GridLayout()
	gc = gcd[1, 1] = GridLayout()
	gd = gcd[2, 1] = GridLayout()
	ge = f[1:2, 3] = GridLayout()

	##-- A --##
	axtop = Axis(ga[1, 1])
	
	df = df_i_normal
	scatter!(axtop, df[!, :ground_truth_mass_large], df[!, :calculated_mass_large], label="Large Inserts")
	scatter!(axtop, df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium], label="Medium Inserts")
	scatter!(axtop, df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], label="Small Inserts", color=:red)
	lines!(axtop, [-1000, 1000], [-1000, 1000], label="Unity")

	xlims!(axtop, low=0, high=200)
	ylims!(axtop, low=0, high=200)
	axtop.xticks = [0, 50, 100, 150, 200]
	axtop.yticks = [0, 50, 100, 150, 200]
	axtop.xlabel = "Known Mass (mg)"
	axtop.ylabel = "Calculated Mass (mg)"
	# Label(ga[1:end, 1:end, Top()], "Integrated", valign = :bottom, padding = (0, 0, 120, 0), textsize = 25)
	# Label(ga[1:end, 1:end, Left()], "High Density", valign = :center, padding = (0, 120, 0, 0), rotation = pi/2, textsize=25)

	##-- B --##
	axbottom = Axis(gb[1, 1])
	
	df2 = df_i_low
	scatter!(axbottom, df2[!, :ground_truth_mass_large], df2[!, :calculated_mass_large], label="Large Inserts")
	scatter!(axbottom, df2[!, :ground_truth_mass_medium], df2[!, :calculated_mass_medium], label="Medium Inserts")
	scatter!(axbottom, df2[!, :ground_truth_mass_small], df2[!, :calculated_mass_small], label="Small Inserts", color=:red)
	lines!(axbottom, [-1000, 1000], [-1000, 1000], label="Unity")

	xlims!(axbottom, low=0, high=25)
	ylims!(axbottom, low=0, high=25)
	axbottom.xticks = [0, 5, 10, 15, 20, 25]
	axbottom.yticks = [0, 5, 10, 15, 20, 25]
	axbottom.xlabel = "Known Mass (mg)"
	axbottom.ylabel = "Calculated Mass (mg)"
	# Label(gb[1:end, 1:end, Left()], "Low Density", valign = :center, padding = (0, 120, 0, 0), rotation = pi/2, textsize=25)

	##-- C --##
	axtopright = Axis(gc[1, 1])
	
	df3 = df_a_normal
	scatter!(axtopright, df3[!, :ground_truth_mass_large], df3[!, :calculated_mass_large], label="Large Inserts")
	scatter!(axtopright, df3[!, :ground_truth_mass_medium], df3[!, :calculated_mass_medium], label="Medium Inserts")
	scatter!(axtopright, df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small], label="Small Inserts", color=:red)
	lines!(axtopright, [-1000, 1000], [-1000, 1000], label="Unity")

	xlims!(axtopright, low=0, high=200)
	ylims!(axtopright, low=0, high=200)
	axtopright.xticks = [0, 50, 100, 150, 200]
	axtopright.yticks = [0, 50, 100, 150, 200]
	axtopright.xlabel = "Known Mass (mg)"
	axtopright.ylabel = "Calculated Mass (mg)"
	# Label(gc[1:end, 1:end, Top()], "Agatston", valign = :bottom, padding = (0, 0, 120, 0), textsize = 25)

	##-- D --##
	axbottomright = Axis(gd[1, 1])
	
	df4 = df_a_low
	scatter!(axbottomright, df4[!, :ground_truth_mass_large], df4[!, :calculated_mass_large], label="Large Inserts")
	scatter!(axbottomright, df4[!, :ground_truth_mass_medium], df4[!, :calculated_mass_medium], label="Medium Inserts")
	scatter!(axbottomright, df4[!, :ground_truth_mass_small], df4[!, :calculated_mass_small], label="Small Inserts", color=:red)
	lines!(axbottomright, [-1000, 1000], [-1000, 1000], label="Unity")

	xlims!(axbottomright, low=0, high=25)
	ylims!(axbottomright, low=0, high=25)
	axbottomright.xticks = [0, 5, 10, 15, 20, 25]
	axbottomright.yticks = [0, 5, 10, 15, 20, 25]
	axbottomright.xlabel = "Known Mass (mg)"
	axbottomright.ylabel = "Calculated Mass (mg)"

	##-- LABELS --##

	f[1:2, 3] = Legend(f, axbottomright, framevisible = false)

	
	for (label, layout) in zip(["A", "B", "C", "D"], [ga, gb, gc, gd])
	    Label(layout[1, 1, TopLeft()], label,
	        textsize = 25,
	        padding = (0, 60, 25, 0),
	        halign = :right)
	end

	# colsize!(f.layout, 1, Auto(0.5))
	# rowsize!(gcd, 1, Auto(1.5))

	f
end

# ╔═╡ 1ef8e352-7e89-461f-9ac8-8a39cf2c14d7
let
	df = df_i_normal
	df2 = df_i_low
	# gt_array = Array(df[!, :ground_truth_mass_large])
	gt_array = vec(hcat(df[!, :ground_truth_mass_large], df[!, :ground_truth_mass_medium], df[!, :ground_truth_mass_small], df2[!, :ground_truth_mass_large], df2[!, :ground_truth_mass_medium], df2[!, :ground_truth_mass_small]))
	# calc_array = Array(df[!, :calculated_mass_large])
	calc_array = vec(hcat(df[!, :calculated_mass_large], df[!, :calculated_mass_medium], df[!, :calculated_mass_small], df2[!, :calculated_mass_large], df2[!, :calculated_mass_medium], df2[!, :calculated_mass_small]))
	data = DataFrame(
		X = gt_array,
		Y= calc_array
	)
	model = lm(@formula(Y ~ X), data)
	r2 = r2(model)
end

# ╔═╡ c86c3d16-ace0-4de2-a7fb-267ea2925302
let
	df = df_a_normal
	df2 = df_a_low
	# gt_array = Array(df[!, :ground_truth_mass_large])
	gt_array = vec(hcat(df[!, :ground_truth_mass_large], df[!, :ground_truth_mass_medium], df[!, :ground_truth_mass_small], df2[!, :ground_truth_mass_large], df2[!, :ground_truth_mass_medium], df2[!, :ground_truth_mass_small]))
	# calc_array = Array(df[!, :calculated_mass_large])
	calc_array = vec(hcat(df[!, :calculated_mass_large], df[!, :calculated_mass_medium], df[!, :calculated_mass_small], df2[!, :calculated_mass_large], df2[!, :calculated_mass_medium], df2[!, :calculated_mass_small]))
	data = DataFrame(
		X = gt_array,
		Y= calc_array
	)
	model = lm(@formula(Y ~ X), data)
	r2 = r2(model)
end

# ╔═╡ d4c89efd-fb4b-4553-8920-847c393cb4bc
let
	dfs = df_s_normal
	dfs2 = df_s_low
	
	dfk = df_a_normal
	dfk2 = df_a_low
	
	gt_array = vec(hcat(dfk[!, :ground_truth_mass_large], dfk[!, :ground_truth_mass_medium], dfk[!, :ground_truth_mass_small], dfk2[!, :ground_truth_mass_large], dfk2[!, :ground_truth_mass_medium], dfk2[!, :ground_truth_mass_small]))
	calc_array = vec(hcat(dfs[!, :calculated_swcs_large], dfs[!, :calculated_swcs_medium], dfs[!, :calculated_swcs_small], dfs2[!, :calculated_swcs_large], dfs2[!, :calculated_swcs_medium], dfs2[!, :calculated_swcs_small]))
	
	data = DataFrame(
		X = gt_array,
		Y= calc_array
	)
	model = lm(@formula(Y ~ X), data)
	r2 = r2(model)
end

# ╔═╡ dd07e78f-1f65-4510-891d-ce4f56fdac05
md"""
## Figure 4. RMSE Mass Score
"""

# ╔═╡ e0a3618d-9686-41a9-932e-0381df126a2d
begin
	agat_arr = [agat_normal_small_rms, agat_low_small_rms, agat_normal_medium_rms, agat_low_medium_rms, agat_normal_large_rms, agat_low_large_rms]
	int_arr = [int_normal_small_rmsd, int_low_small_rms, int_normal_medium_rmsd, int_low_medium_rms, int_normal_large_rmsd, int_low_large_rms]
end;

# ╔═╡ 787ba4f2-8747-48c6-a03e-705ee8e07aa3
total_rms_agat = mean(agat_arr)

# ╔═╡ de4d6989-969a-4a4d-a089-7252031df010
total_rms_int = mean(int_arr)

# ╔═╡ 73dd7766-befd-41e9-bdfb-0632e9500326
val_rms = EqualVarianceTTest(agat_arr, int_arr)

# ╔═╡ 8f4de118-39af-44b9-a408-3f5026c33fe5
let
	f = Figure()
	ga = f[1, 1] = GridLayout()
	gb = f[2, 1] = GridLayout()

	colors = Makie.wong_colors()
	labels = ["Agatston", "Integrated"]
	elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]

	##-- TOP --##
	axtop = Axis(ga[1, 1])
	hidexdecorations!(axtop, grid = false)
	axtop.ylabel = "RMSD"
	labels = ["Agatston", "Integrated"]
	table = [1, 1, 2, 2, 3, 3]
	grp = [1, 2, 1, 2, 1, 2]
	heights1 = [agat_normal_small_rms, int_normal_small_rmsd, agat_normal_medium_rms, int_normal_medium_rmsd, agat_normal_large_rms, int_normal_medium_rmsd]
	barplot!(axtop, table, heights1, dodge=grp, color=colors[grp])
	# Label(ga[1:end, 1:end, Left()], "High Density", valign = :center, padding = (0, 120, 0, 0), rotation = pi/2, textsize=25)


	##-- BOTTOM --##
	axbottom = Axis(gb[1, 1], xticks = (1:3, ["Small Patient", "Medium Patient", "Large Patient"]))
	axbottom.ylabel = "RMSD"	
	heights2 = [agat_low_small_rms, int_low_small_rms, agat_low_medium_rms, int_low_medium_rms, agat_low_large_rms, int_low_medium_rms]
	barplot!(axbottom, table, heights2, dodge=grp, color=colors[grp])
	linkaxes!(axtop, axbottom)
	# Label(gb[1:end, 1:end, Left()], "Low Density", valign = :center, padding = (0, 120, 0, 0), rotation = pi/2, textsize=25)

	for (label, layout) in zip(["A", "B"], [ga, gb])
	    Label(layout[1, 1, TopLeft()], label,
	        textsize = 25,
	        padding = (0, 60, 25, 0),
	        halign = :right)
	end

	Legend(f[1:2, 2], elements, labels, framevisible = false)
	
	f
end

# ╔═╡ 7ec78f71-6b20-4c8c-a9da-a216404bee72
md"""
## Figure 5. Reproducibility
"""

# ╔═╡ 387fafb4-4a14-482f-8a8e-4214217f82a5
md"""
#### Integrated (Repeated)
"""

# ╔═╡ 53c1b176-e2f7-4cb9-baa9-26d61ab8c18f
path_integrated_r = "/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/output_repeated/integrated_scoring";

# ╔═╡ 89b28ac5-dd69-4812-84fd-64b54606d146
df_i_r = CSV.read(string(path_integrated_r, "/full.csv"), DataFrame);

# ╔═╡ 8ade536a-6ed2-481d-8b4b-cb6e1d57aca4
df_i_r

# ╔═╡ 4f5328e7-66dd-433c-a7f6-5dc461c83a4c
df_i_low_r, df_i_normal_r = groupby(df_i_r, :DENSITY);

# ╔═╡ 6e515ee8-2836-4285-a0a9-131e1e7b4b7f
df_i_low_small_r, df_i_low_medium_r, df_i_low_large_r = groupby(df_i_low_r, :SIZE);

# ╔═╡ 9ae3b4cb-fac5-4168-868d-724de9b0c9d2
df_i_normal_small_r, df_i_normal_medium_r, df_i_normal_large_r = groupby(df_i_normal_r, :SIZE);

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

# ╔═╡ a85f777a-76f2-4c64-9973-ea9dec245600
let
	f = Figure()
	
	ga = f[1, 1] = GridLayout()
	gb = f[2, 1] = GridLayout()
	gc = f[1, 2] = GridLayout()
	gd = f[2, 2] = GridLayout()
	ge = f[1:2, 3] = GridLayout()

	##-- A --##
	axtop = Axis(ga[1, 1])
	
	df = df_i
	df_r = df_i_r
	scatter!(axtop, df_r[!, :calculated_mass_large], df[!, :calculated_mass_large], label="Large Inserts")
	scatter!(axtop, df_r[!, :calculated_mass_medium], df[!, :calculated_mass_medium], label="Medium Inserts")
	scatter!(axtop, df_r[!, :calculated_mass_small], df[!, :calculated_mass_small], label="Small Inserts", color=:red)
	lines!(axtop, [-1000, 1000], [-1000, 1000], label="Unity")

	xlims!(axtop, low=0, high=200)
	ylims!(axtop, low=0, high=200)
	axtop.xticks = [0, 50, 100, 150, 200]
	axtop.yticks = [0, 50, 100, 150, 200]
	axtop.xlabel = "Calculated Mass 1 (mg)"
	axtop.ylabel = "Calculated Mass 2 (mg)"

	# ##-- B --##
	axbottomright = Axis(gb[1, 1])
	
	df4 = df_s
	df4_r = df_s_r
	scatter!(axbottomright, df4_r[!, :calculated_swcs_large], df4[!, :calculated_swcs_large], label="Large Inserts")
	scatter!(axbottomright, df4_r[!, :calculated_swcs_medium], df4[!, :calculated_swcs_medium], label="Medium Inserts")
	scatter!(axbottomright, df4_r[!, :calculated_swcs_small], df4[!, :calculated_swcs_small], label="Small Inserts", color=:red)
	lines!(axbottomright, [-1000, 1000], [-1000, 1000], label="Unity")

	xlims!(axbottomright, low=0, high=500)
	ylims!(axbottomright, low=0, high=500)
	# axbottomright.xticks = [0, 5, 10, 15, 20, 25]
	# axbottomright.yticks = [0, 5, 10, 15, 20, 25]
	axbottomright.xlabel = "SWCS 1"
	axbottomright.ylabel = "SWCS 2"

	##-- C --##
	axtopright = Axis(gc[1, 1])
	
	df3 = df_a
	df3_r = df_a_r
	scatter!(axtopright, df3_r[!, :calculated_mass_large], df3[!, :calculated_mass_large], label="Large Inserts")
	scatter!(axtopright, df3_r[!, :calculated_mass_medium], df3[!, :calculated_mass_medium], label="Medium Inserts")
	scatter!(axtopright, df3_r[!, :calculated_mass_small], df3[!, :calculated_mass_small], label="Small Inserts", color=:red)
	lines!(axtopright, [-1000, 1000], [-1000, 1000], label="Unity")

	xlims!(axtopright, low=0, high=200)
	ylims!(axtopright, low=0, high=200)
	axtopright.xticks = [0, 50, 100, 150, 200]
	axtopright.yticks = [0, 50, 100, 150, 200]
	axtopright.xlabel = "Calculated Mass 1 (mg)"
	axtopright.ylabel = "Calculated Mass 2 (mg)"

	##-- D --##

	##-- LABELS --##
	f[1:2, 3] = Legend(f, axbottomright, framevisible = false)

	
	for (label, layout) in zip(["A", "B", "C", "D"], [ga, gb, gc, gd])
	    Label(layout[1, 1, TopLeft()], label,
	        textsize = 25,
	        padding = (0, 60, 25, 0),
	        halign = :right)
	end

	# colsize!(f.layout, 1, Auto(0.5))
	# rowsize!(gcd, 1, Auto(1.5))

	f
end

# ╔═╡ 4f19aa51-ba28-4433-bc75-dc21e55da969
let
	df = df_i
	df2 = df_i_r
	# gt_array = Array(df[!, :ground_truth_mass_large])
	arr1 = vec(hcat(df[!, :calculated_mass_large], df[!, :calculated_mass_medium], df[!, :calculated_mass_small]))
	arr2 = vec(hcat(df2[!, :calculated_mass_large], df2[!, :calculated_mass_medium], df2[!, :calculated_mass_small]))
	data = DataFrame(
		X = arr1,
		Y= arr2
	)
	model = lm(@formula(Y ~ X), data)
	r2 = GLM.r2(model)
end

# ╔═╡ 6b61bd3d-4d4b-4304-9591-06f8338b75cc
let
	df = df_s
	df2 = df_s_r
	arr1 = vec(hcat(df[!, :calculated_swcs_large], df[!, :calculated_swcs_medium], df[!, :calculated_swcs_small]))
	arr2 = vec(hcat(df2[!, :calculated_swcs_large], df2[!, :calculated_swcs_medium], df2[!, :calculated_swcs_small]))
	data = DataFrame(
		X = arr1,
		Y= arr2
	)
	model = lm(@formula(Y ~ X), data)
	r2 = GLM.r2(model)
end

# ╔═╡ 84d5724b-de74-40e8-a5c7-e22b982e99d7
let
	df = df_a
	df2 = df_a_r
	arr1 = vec(hcat(df[!, :calculated_mass_large], df[!, :calculated_mass_medium], df[!, :calculated_mass_small]))
	arr2 = vec(hcat(df2[!, :calculated_mass_large], df2[!, :calculated_mass_medium], df2[!, :calculated_mass_small]))
	data = DataFrame(
		X = arr1,
		Y= arr2
	)
	model = lm(@formula(Y ~ X), data)
	r2 = GLM.r2(model)
end

# ╔═╡ 40f443f8-6e6b-4de3-9c2e-b70599640c5d
md"""
# TEMPLATE FIGURE
"""

# ╔═╡ f33bd1c6-009a-4071-9cb5-93f4f71c5415
let
	f = Figure()
	
	ga = f[1, 1] = GridLayout()
	gb = f[2, 1] = GridLayout()
	gcd = f[1:2, 2] = GridLayout()
	gc = gcd[1, 1] = GridLayout()
	gd = gcd[2, 1] = GridLayout()

	##-- A --##
	axtop = Axis(ga[1, 1])
	axmain = Axis(ga[2, 1], xlabel = "before", ylabel = "after")
	axright = Axis(ga[2, 2])
	
	linkyaxes!(axmain, axright)
	linkxaxes!(axmain, axtop)
	
	labels = ["treatment", "placebo", "control"]
	data = randn(3, 100, 2) .+ [1, 3, 5]
	
	for (label, col) in zip(labels, eachslice(data, dims = 1))
	    scatter!(axmain, col, label = label)
	    density!(axtop, col[:, 1])
	    density!(axright, col[:, 2], direction = :y)
	end
	
	ylims!(axtop, low = 0)
	xlims!(axright, low = 0)
	axmain.xticks = 0:3:9
	leg = Legend(ga[1, 2], axmain)
	hidedecorations!(axtop, grid = false)
	hidedecorations!(axright, grid = false)
	leg.tellheight = true
	colgap!(ga, 10)
	rowgap!(ga, 10)
	
	Label(ga[1, 1:2, Top()], "Stimulus ratings", valign = :bottom,
    padding = (0, 0, 5, 0))

	##-- B --##
	xs = LinRange(0.5, 6, 50)
	ys = LinRange(0.5, 6, 50)
	data1 = [sin(x^1.5) * cos(y^0.5) for x in xs, y in ys] .+ 0.1 .* randn.()
	data2 = [sin(x^0.8) * cos(y^1.5) for x in xs, y in ys] .+ 0.1 .* randn.()
	
	ax1, hm = contourf(gb[1, 1], xs, ys, data1,
	    levels = 6)
	ax1.title = "Histological analysis"
	contour!(ax1, xs, ys, data1, levels = 5, color = :black)
	hidexdecorations!(ax1)
	
	ax2, hm2 = contourf(gb[2, 1], xs, ys, data2,
	    levels = 6)
	contour!(ax2, xs, ys, data2, levels = 5, color = :black)

	cb = Colorbar(gb[1:2, 2], hm, label = "cell group")
	low, high = extrema(data1)
	edges = range(low, high, length = 7)
	centers = (edges[1:6] .+ edges[2:7]) .* 0.5
	cb.ticks = (centers, string.(1:6))
	cb.alignmode = Mixed(right = 0)
	colgap!(gb, 10)
	rowgap!(gb, 10)

	##-- C --##
	# ...

	##-- D --##
	axs = [Axis(gd[row, col]) for row in 1:3, col in 1:2]
	hidedecorations!.(axs, grid = false, label = false)
	
	for row in 1:3, col in 1:2
	    xrange = col == 1 ? (0:0.1:6pi) : (0:0.1:10pi)
	
	    eeg = [sum(sin(pi * rand() + k * x) / k for k in 1:10)
	        for x in xrange] .+ 0.1 .* randn.()
	
	    lines!(axs[row, col], eeg, color = (:black, 0.5))
	end
	
	axs[3, 1].xlabel = "Day 1"
	axs[3, 2].xlabel = "Day 2"
	Label(gd[1, :, Top()], "EEG traces", valign = :bottom,
	    font = "TeX Gyre Heros Bold",
	    padding = (0, 0, 5, 0))
	rowgap!(gd, 10)
	colgap!(gd, 10)
	for (i, label) in enumerate(["sleep", "awake", "test"])
	    Box(gd[i, 3], color = :gray90)
	    Label(gd[i, 3], label, rotation = pi/2, tellheight = false)
	end
	colgap!(gd, 2, 0)
	n_day_1 = length(0:0.1:6pi)
	n_day_2 = length(0:0.1:10pi)
	
	colsize!(gd, 1, Auto(n_day_1))
	colsize!(gd, 2, Auto(n_day_2))

	##-- LABELS --##
	for (label, layout) in zip(["A", "B", "C", "D"], [ga, gb, gc, gd])
	    Label(layout[1, 1, TopLeft()], label,
	        textsize = 26,
	        font = "TeX Gyre Heros Bold",
	        padding = (0, 5, 5, 0),
	        halign = :right)
	end

	# colsize!(f.layout, 1, Auto(0.5))
	# rowsize!(gcd, 1, Auto(1.5))

	f
end

# ╔═╡ Cell order:
# ╠═c82bb3e6-e052-11ec-3392-eba4a4b464ba
# ╠═2df207ed-cf9c-4d7d-8354-4da14f93276c
# ╟─4cf61487-3913-4e11-970e-4ca31a5ffc8d
# ╠═57d46998-368a-4916-8380-ee49d5473a49
# ╠═2c960bd8-ae64-453f-b29f-275bf5263774
# ╠═7f4fae09-7916-4a98-a102-7f861900c457
# ╠═c987e55d-9ffd-40a5-9d64-1653d4000837
# ╠═0c367cb1-d959-4c10-b31c-18d3cb396255
# ╠═5144e2a2-91e2-49ab-9dcf-df57ed767a36
# ╠═8ade536a-6ed2-481d-8b4b-cb6e1d57aca4
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
# ╟─932e1a11-ec52-48bc-af27-bb3e09e6b27f
# ╠═dc8805a5-2b4c-47fe-88c8-0f4434ddd5cf
# ╠═2c78dcc4-03c7-49b5-ac45-30c198671761
# ╠═d3e13dcd-201a-49ab-86f7-a7d3320088c2
# ╠═4b97d892-b0db-4f75-9239-9da215f905f0
# ╠═f3cfe6db-d3f5-43c9-988e-25f332a5998a
# ╠═8e4fc2cf-16fd-43f2-84f1-8339b6181baf
# ╟─01d903d8-fca2-4317-a2d7-93429f06d47c
# ╠═97a612d5-44e7-43f1-b61e-8518c7e6aa28
# ╠═a31dfad2-e420-4c41-b805-227778a39fc9
# ╟─86d5f7fd-6a1b-4fdc-9e06-6cd66f8a169c
# ╠═aabe8e7b-faed-499d-aea2-fd44836a3cfb
# ╠═dd644657-51e9-4b27-a186-6efba0bb83f8
# ╠═d79d6478-0837-40f0-a7a5-270bf44d6166
# ╠═2b5a8dca-faae-4163-ab14-ae4f98b2a7cc
# ╠═f4b51729-d53a-40e2-919b-5d4562cffbc0
# ╟─b702760a-9d9b-4d2a-992c-37ed4f8bcfdf
# ╟─b1a054d6-863f-4675-bb59-01196565bee6
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
# ╠═7ddde041-6e01-4f4c-b930-c2b5143fd51f
# ╠═b3d14f1b-d538-415a-bb8d-df9abb36be70
# ╠═e7c525f4-e9da-424e-9483-e2d97832e19d
# ╠═e89f1021-7a7a-4d32-82b3-dd93d2823388
# ╠═8fd6011f-24b4-43a8-809e-371ef33ceedc
# ╠═c05610e1-ecc6-4391-b6cf-964bb4696362
# ╟─3b3269ef-1e0a-4335-b0a1-e14223670682
# ╟─7dfc24a4-e006-45f4-b5b9-977a7c3c0b7c
# ╟─d04ce8a0-8224-49c7-a057-f1c7d0d65ef8
# ╠═1ef8e352-7e89-461f-9ac8-8a39cf2c14d7
# ╠═c86c3d16-ace0-4de2-a7fb-267ea2925302
# ╠═d4c89efd-fb4b-4553-8920-847c393cb4bc
# ╟─dd07e78f-1f65-4510-891d-ce4f56fdac05
# ╠═e0a3618d-9686-41a9-932e-0381df126a2d
# ╠═787ba4f2-8747-48c6-a03e-705ee8e07aa3
# ╠═de4d6989-969a-4a4d-a089-7252031df010
# ╠═73dd7766-befd-41e9-bdfb-0632e9500326
# ╟─8f4de118-39af-44b9-a408-3f5026c33fe5
# ╟─7ec78f71-6b20-4c8c-a9da-a216404bee72
# ╟─387fafb4-4a14-482f-8a8e-4214217f82a5
# ╠═53c1b176-e2f7-4cb9-baa9-26d61ab8c18f
# ╠═89b28ac5-dd69-4812-84fd-64b54606d146
# ╠═4f5328e7-66dd-433c-a7f6-5dc461c83a4c
# ╠═6e515ee8-2836-4285-a0a9-131e1e7b4b7f
# ╠═9ae3b4cb-fac5-4168-868d-724de9b0c9d2
# ╟─5248fac3-4dc2-4579-a52b-8d04ebf45b07
# ╠═fcd83805-06e2-4cda-aa56-0c13c69424d8
# ╠═5f2eb2d2-84c4-4160-a92a-f181b4126450
# ╠═869de9cc-6ab0-4e0b-ad1f-59db1fd2f69f
# ╠═43b3ba9b-247c-41cd-968e-3a27736426b0
# ╠═b16090b5-a232-4e71-9148-b71296efa999
# ╟─c1d820ac-d0ca-4c18-b61d-849402ae647e
# ╠═594b8d0e-f65c-4fe3-b571-e2ee68a9ab5f
# ╠═1d208fd8-b847-4afe-84b5-3a553f50a858
# ╠═a9f68d7e-7818-47d5-ba1f-965634225b30
# ╠═f2177e1f-2a55-4f81-9a69-ab985ae3b7c2
# ╠═35bc90e2-9d70-4daf-a825-b462831c5bf6
# ╟─a85f777a-76f2-4c64-9973-ea9dec245600
# ╠═4f19aa51-ba28-4433-bc75-dc21e55da969
# ╠═6b61bd3d-4d4b-4304-9591-06f8338b75cc
# ╠═84d5724b-de74-40e8-a5c7-e22b982e99d7
# ╟─40f443f8-6e6b-4de3-9c2e-b70599640c5d
# ╠═f33bd1c6-009a-4071-9cb5-93f4f71c5415

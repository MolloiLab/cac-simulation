### A Pluto.jl notebook ###
# v0.19.2

using Markdown
using InteractiveUtils

# ╔═╡ 0b4fbf14-cb07-11ec-2bf5-cf802f56d04b
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
end

# ╔═╡ 08788f31-c822-4d0b-8eca-6e83c4cc9f12
TableOfContents()

# ╔═╡ 9aabc4c6-36ed-47da-9e7c-630d505ea8e2
md"""
## Integrated
"""

# ╔═╡ 021474f3-9d65-4f8d-b758-c6488ba94d10
path_integrated = "/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/output/integrated_scoring";

# ╔═╡ 509cb793-4d6d-47ee-a8a2-edf018e73fe5
df_i = CSV.read(string(path_integrated, "/full.csv"), DataFrame);

# ╔═╡ 0cf92498-d46c-4604-a582-269f46b3ae93
df_i_low, df_i_normal = groupby(df_i, :DENSITY);

# ╔═╡ 642ff265-d962-41bc-a35f-af0cb054d870
df_i_low_small, df_i_low_medium, df_i_low_large = groupby(df_i_low, :SIZE);

# ╔═╡ 2e78a50c-3798-4502-bc58-becdfab06a16
df_i_normal_small, df_i_normal_medium, df_i_normal_large = groupby(df_i_normal, :SIZE);

# ╔═╡ 198bbf97-7931-4acc-8748-2d6f364b38c0
md"""
### Normal Density
"""

# ╔═╡ 631664c6-0d4d-4793-89a9-f9097efa0986
let
	f = Figure()
	ax = Axis(f[1, 1])

	df = df_i_normal_small

	scatter!(df[!, :ground_truth_mass_large], df[!, :calculated_mass_large], label="Large Inserts")
	scatter!(df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium], label="Medium Inserts")
	scatter!(df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], label="Small Inserts", color=:red)
	lines!([0, 100], [0, 100], label="Unity")
	
	ax.title = "Mass vs Known Mass (Small Patient)"
	ax.ylabel = "Calculated Mass (mg)"
	ax.xlabel = "Known Mass (mg)"
	ax.xticks = [0, 25, 50, 75, 100]
	ax.yticks = [0, 25, 50, 75, 100]

	xlims!(ax, 0, 100)
	ylims!(ax, 0, 100)
	
	f[1, 2] = Legend(f, ax, framevisible = false)
	f
end

# ╔═╡ 79325025-f1f2-425e-bb6c-f52972b4768b
let
	f = Figure()
	ax = Axis(f[1, 1])

	df = df_i_normal_medium

	scatter!(df[!, :ground_truth_mass_large], df[!, :calculated_mass_large], label="Large Inserts")
	scatter!(df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium], label="Medium Inserts")
	scatter!(df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], label="Small Inserts", color=:red)
	lines!([0, 100], [0, 100], label="Unity")
	
	ax.title = "Mass vs Known Mass (Medium Patient)"
	ax.ylabel = "Calculated Mass (mg)"
	ax.xlabel = "Known Mass (mg)"
	ax.xticks = [0, 25, 50, 75, 100]
	ax.yticks = [0, 25, 50, 75, 100]

	xlims!(ax, 0, 100)
	ylims!(ax, 0, 100)
	
	f[1, 2] = Legend(f, ax, framevisible = false)
	f
end

# ╔═╡ 1cdfd787-0957-4933-bef8-aa6a7084325c
let
	f = Figure()
	ax = Axis(f[1, 1])

	df = df_i_normal_large

	scatter!(df[!, :ground_truth_mass_large], df[!, :calculated_mass_large], label="Large Inserts")
	scatter!(df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium], label="Medium Inserts")
	scatter!(df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], label="Small Inserts", color=:red)
	lines!([0, 100], [0, 100], label="Unity")
	
	ax.title = "Mass vs Known Mass (Large Patient)"
	ax.ylabel = "Calculated Mass (mg)"
	ax.xlabel = "Known Mass (mg)"
	ax.xticks = [0, 25, 50, 75, 100]
	ax.yticks = [0, 25, 50, 75, 100]

	xlims!(ax, 0, 100)
	ylims!(ax, 0, 100)
	
	f[1, 2] = Legend(f, ax, framevisible = false)
	f
end

# ╔═╡ 51a79c37-dbbd-41af-b3ca-47cb67fa0ae9
let
	df = df_i_normal_small
	global int_normal_small_rms
	int_normal_small_rms = (rmsd(df[!, :calculated_mass_large], df[!, :ground_truth_mass_large]) + rmsd(df[!, :calculated_mass_medium], df[!, :ground_truth_mass_medium]) + rmsd(df[!, :calculated_mass_small], df[!, :ground_truth_mass_small]))/ 3
end

# ╔═╡ d114d2db-d813-4aff-b993-5df344e2317f
let
	df = df_i_normal_medium
	global int_normal_medium_rms
	int_normal_medium_rms = (rmsd(df[!, :calculated_mass_large], df[!, :ground_truth_mass_large]) + rmsd(df[!, :calculated_mass_medium], df[!, :ground_truth_mass_medium]) + rmsd(df[!, :calculated_mass_small], df[!, :ground_truth_mass_small]))/ 3
end

# ╔═╡ 42f95b3b-016b-4bd1-a139-0bdfbdbbae38
let
	df = df_i_normal_large
	global int_normal_large_rms
	int_normal_large_rms = (rmsd(df[!, :calculated_mass_large], df[!, :ground_truth_mass_large]) + rmsd(df[!, :calculated_mass_medium], df[!, :ground_truth_mass_medium]) + rmsd(df[!, :calculated_mass_small], df[!, :ground_truth_mass_small])) / 3
end

# ╔═╡ 3e321850-eb7c-4620-99b4-60025921c10b
let
	tbl = [1, 2, 3]
    height = [int_normal_small_rms, int_normal_medium_rms, int_normal_large_rms]
	barplot(tbl, height, axis = (xticks = (1:3, ["small", "medium", "large"]),))
end

# ╔═╡ Cell order:
# ╠═0b4fbf14-cb07-11ec-2bf5-cf802f56d04b
# ╠═08788f31-c822-4d0b-8eca-6e83c4cc9f12
# ╟─9aabc4c6-36ed-47da-9e7c-630d505ea8e2
# ╠═021474f3-9d65-4f8d-b758-c6488ba94d10
# ╠═509cb793-4d6d-47ee-a8a2-edf018e73fe5
# ╠═0cf92498-d46c-4604-a582-269f46b3ae93
# ╠═642ff265-d962-41bc-a35f-af0cb054d870
# ╠═2e78a50c-3798-4502-bc58-becdfab06a16
# ╟─198bbf97-7931-4acc-8748-2d6f364b38c0
# ╟─631664c6-0d4d-4793-89a9-f9097efa0986
# ╟─79325025-f1f2-425e-bb6c-f52972b4768b
# ╟─1cdfd787-0957-4933-bef8-aa6a7084325c
# ╠═51a79c37-dbbd-41af-b3ca-47cb67fa0ae9
# ╠═d114d2db-d813-4aff-b993-5df344e2317f
# ╠═42f95b3b-016b-4bd1-a139-0bdfbdbbae38
# ╠═3e321850-eb7c-4620-99b4-60025921c10b

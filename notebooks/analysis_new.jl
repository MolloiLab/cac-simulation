### A Pluto.jl notebook ###
# v0.19.4

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

# ╔═╡ a09621a9-ea7a-4722-9b96-136ca1114edf
md"""
### -- Small Inserts
"""

# ╔═╡ 2425b571-6931-4c97-9998-7e382cd5b941
let
	f = Figure()
	ax = Axis(f[1, 1])

	df = df_i_normal_large
	scatter!(df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], label="Small Inserts", color=:red)
	lines!([-100, 100], [-100, 100], label="Unity")
	
	ax.title = "Mass vs Known Mass (Large Patient)"
	ax.ylabel = "Calculated Mass (mg)"
	ax.xlabel = "Known Mass (mg)"

	xlims!(ax, -1, 10)
	ylims!(ax, -1, 10)
	
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

# ╔═╡ b0a5ec8a-3630-43fe-aeb4-2ae2d5f66d55
df_i_rms_normal = DataFrame(
	small_rms = int_normal_small_rms,
	medium_rms = int_normal_medium_rms,
	large_rms = int_normal_large_rms
)

# ╔═╡ 3e321850-eb7c-4620-99b4-60025921c10b
let
	f = Figure()
	df = df_i_rms_normal
	
	ax = Makie.Axis(
		f[1, 1],
		xticks = (1:3, ["Small Patient", "Medium Patient", "Large Patient"]),
	)
	ax.title = "RMS (Large)"
	ax.ylabel = "RMSD"
	
	barplot!(ax, [1], df[!, :small_rms])
	barplot!(ax, [2], df[!, :medium_rms])
	barplot!(ax, [3], df[!, :large_rms])

	f
end

# ╔═╡ 3776dcd4-e1eb-4f0a-ab80-5f770ae66d83
md"""
### Low Density
"""

# ╔═╡ 31933102-a867-4d9a-9af2-5cdc8597025d
let
	f = Figure()
	ax = Axis(f[1, 1])

	df = df_i_low_small

	scatter!(df[!, :ground_truth_mass_large], df[!, :calculated_mass_large], label="Large Inserts")
	scatter!(df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium], label="Medium Inserts")
	scatter!(df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], label="Small Inserts", color=:red)
	lines!([0, 100], [0, 100], label="Unity")
	
	ax.title = "Mass vs Known Mass (Small Patient)"
	ax.ylabel = "Calculated Mass (mg)"
	ax.xlabel = "Known Mass (mg)"
	ax.xticks = [0, 5, 10, 15, 20]
	ax.yticks = [0, 5, 10, 15, 20]

	xlims!(ax, 0, 20)
	ylims!(ax, 0, 20)
	
	f[1, 2] = Legend(f, ax, framevisible = false)
	f
end

# ╔═╡ 5c663e87-ad90-4f73-8a6f-ca62a522dab9
let
	f = Figure()
	ax = Axis(f[1, 1])

	df = df_i_low_medium

	scatter!(df[!, :ground_truth_mass_large], df[!, :calculated_mass_large], label="Large Inserts")
	scatter!(df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium], label="Medium Inserts")
	scatter!(df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], label="Small Inserts", color=:red)
	lines!([0, 100], [0, 100], label="Unity")
	
	ax.title = "Mass vs Known Mass (Medium Patient)"
	ax.ylabel = "Calculated Mass (mg)"
	ax.xlabel = "Known Mass (mg)"
	ax.xticks = [0, 5, 10, 15, 20]
	ax.yticks = [0, 5, 10, 15, 20]

	xlims!(ax, 0, 20)
	ylims!(ax, 0, 20)
	
	f[1, 2] = Legend(f, ax, framevisible = false)
	f
end

# ╔═╡ 21f9fc96-3560-4d9e-bbf8-947fd6d3a7c8
let
	f = Figure()
	ax = Axis(f[1, 1])

	df = df_i_low_large

	scatter!(df[!, :ground_truth_mass_large], df[!, :calculated_mass_large], label="Large Inserts")
	scatter!(df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium], label="Medium Inserts")
	scatter!(df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], label="Small Inserts", color=:red)
	lines!([0, 100], [0, 100], label="Unity")
	
	ax.title = "Mass vs Known Mass (Large Patient)"
	ax.ylabel = "Calculated Mass (mg)"
	ax.xlabel = "Known Mass (mg)"
	ax.xticks = [0, 5, 10, 15, 20]
	ax.yticks = [0, 5, 10, 15, 20]

	xlims!(ax, 0, 20)
	ylims!(ax, 0, 20)
	
	f[1, 2] = Legend(f, ax, framevisible = false)
	f
end

# ╔═╡ 4ec7c0c9-9bfb-4ef1-bf47-33b17bb4c3b6
md"""
### -- Small Inserts
"""

# ╔═╡ 4a65be29-d4e7-43e9-bba3-cf80304ab5f3
let
	f = Figure()
	ax = Axis(f[1, 1])

	df = df_i_low_large
	scatter!(df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], label="Small Inserts", color=:red)
	lines!([-100, 100], [-100, 100], label="Unity")
	
	ax.title = "Mass vs Known Mass (Large Patient)"
	ax.ylabel = "Calculated Mass (mg)"
	ax.xlabel = "Known Mass (mg)"

	xlims!(ax, 0, 1)
	ylims!(ax, -0.25, 1)
	
	f[1, 2] = Legend(f, ax, framevisible = false)
	f
end

# ╔═╡ abc35a1a-4aad-4a3b-9f27-5ef8b3f6d8a3
let
	df = df_i_low_small
	global int_low_small_rms
	int_low_small_rms = (rmsd(df[!, :calculated_mass_large], df[!, :ground_truth_mass_large]) + rmsd(df[!, :calculated_mass_medium], df[!, :ground_truth_mass_medium]) + rmsd(df[!, :calculated_mass_small], df[!, :ground_truth_mass_small]))/ 3
end

# ╔═╡ 083eb4f4-a60f-4dd6-9ec6-fdf6ff15a697
let
	df = df_i_low_medium
	global int_low_medium_rms
	int_low_medium_rms = (rmsd(df[!, :calculated_mass_large], df[!, :ground_truth_mass_large]) + rmsd(df[!, :calculated_mass_medium], df[!, :ground_truth_mass_medium]) + rmsd(df[!, :calculated_mass_small], df[!, :ground_truth_mass_small]))/ 3
end

# ╔═╡ 5aefde17-282b-4b0c-a08e-75e17f59c26f
let
	df = df_i_low_large
	global int_low_large_rms
	int_low_large_rms = (rmsd(df[!, :calculated_mass_large], df[!, :ground_truth_mass_large]) + rmsd(df[!, :calculated_mass_medium], df[!, :ground_truth_mass_medium]) + rmsd(df[!, :calculated_mass_small], df[!, :ground_truth_mass_small])) / 3
end

# ╔═╡ dae44dae-9d36-4591-8f64-34bc6d531524
df_i_rms_low = DataFrame(
	small_rms = int_low_small_rms,
	medium_rms = int_low_medium_rms,
	large_rms = int_low_large_rms
)

# ╔═╡ 3c1510f9-b62b-48f1-9670-fbbb0d2bbc78
let
	f = Figure()
	df = df_i_rms_low
	
	ax = Makie.Axis(
		f[1, 1],
		xticks = (1:3, ["Small Patient", "Medium Patient", "Large Patient"]),
	)
	ax.title = "RMS (Large)"
	ax.ylabel = "RMSD"
	
	barplot!(ax, [1], df[!, :small_rms])
	barplot!(ax, [2], df[!, :medium_rms])
	barplot!(ax, [3], df[!, :large_rms])

	f
end

# ╔═╡ 14391d86-a2e2-4a6a-98f5-eef12be4502f
md"""
## Agatston
"""

# ╔═╡ 67780e42-059d-40fe-8116-f9140909cc0a
path_agat = "/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/output/agatston";

# ╔═╡ 8144a8c1-1f9a-4de1-8bf2-55d443147540
df_a = CSV.read(string(path_agat, "/full.csv"), DataFrame);

# ╔═╡ a20354d4-37d0-4e2d-a5f4-581a302fde54
df_a_low, df_a_normal = groupby(df_a, :DENSITY);

# ╔═╡ 6b923331-0784-48cf-b4ca-9684ecb15335
df_a_low_small, df_a_low_medium, df_a_low_large = groupby(df_a_low, :SIZE);

# ╔═╡ 817510ce-847c-4c17-95c5-b5445f38e379
df_a_normal_small, df_a_normal_medium, df_a_normal_large = groupby(df_a_normal, :SIZE);

# ╔═╡ dfb96fe4-2cb7-449d-bb68-b5e20299f4e3
df_a_low

# ╔═╡ c3a13ff7-72e3-4233-8530-36dab53e4c94
md"""
### Normal Density
"""

# ╔═╡ fdc99587-f3e5-42cd-a387-c36d0adf8eed
let
	f = Figure()
	ax = Axis(f[1, 1])

	df = df_a_normal_small

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

# ╔═╡ 7ae72022-136c-4594-87b1-94d2102e59dd
let
	f = Figure()
	ax = Axis(f[1, 1])

	df = df_a_normal_medium

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

# ╔═╡ ba3b4f2e-d0e0-47c9-ba61-4741cf296333
let
	f = Figure()
	ax = Axis(f[1, 1])

	df = df_a_normal_large

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

# ╔═╡ 6bbbae82-ebcc-455f-a035-9049291ce94f
let
	df = df_a_normal_small
	global agat_normal_small_rms
	agat_normal_small_rms = (rmsd(df[!, :calculated_mass_large], df[!, :ground_truth_mass_large]) + rmsd(df[!, :calculated_mass_medium], df[!, :ground_truth_mass_medium]) + rmsd(df[!, :calculated_mass_small], df[!, :ground_truth_mass_small]))/ 3
end

# ╔═╡ 9c5d2ad9-3a12-4816-ae2b-e11621f8a783
let
	df = df_a_normal_medium
	global agat_normal_medium_rms
	agat_normal_medium_rms = (rmsd(df[!, :calculated_mass_large], df[!, :ground_truth_mass_large]) + rmsd(df[!, :calculated_mass_medium], df[!, :ground_truth_mass_medium]) + rmsd(df[!, :calculated_mass_small], df[!, :ground_truth_mass_small]))/ 3
end

# ╔═╡ 2ab0eeea-79b2-4d18-bc85-91f415a5f3c6
let
	df = df_a_normal_large
	global agat_normal_large_rms
	agat_normal_large_rms = (rmsd(df[!, :calculated_mass_large], df[!, :ground_truth_mass_large]) + rmsd(df[!, :calculated_mass_medium], df[!, :ground_truth_mass_medium]) + rmsd(df[!, :calculated_mass_small], df[!, :ground_truth_mass_small])) / 3
end

# ╔═╡ fc4c940c-ed2d-45f9-9989-eae7979992a9
df_a_rms = DataFrame(
	normal_small_rms = agat_normal_small_rms,
	normal_medium_rms = agat_normal_medium_rms,
	normal_large_rms = agat_normal_large_rms
)

# ╔═╡ 6830ae86-54d5-4238-84d4-4f81541b4714
let
	f = Figure()
	df = df_a_rms
	
	ax = Makie.Axis(
		f[1, 1],
		xticks = (1:3, ["Small Patient", "Medium Patient", "Large Patient"])
	)
	ax.title = "RMS"
	ax.ylabel = "RMSD"
	barplot!(ax, [1], df[!, :normal_small_rms])
	barplot!(ax, [2], df[!, :normal_medium_rms])
	barplot!(ax, [3], df[!, :normal_large_rms])
	f
end

# ╔═╡ ff238353-1281-40a3-9484-5f7d3177b7d6
md"""
### Low Density
"""

# ╔═╡ 72d80487-acce-4b8d-afa3-28f30ad884ab
let
	f = Figure()
	ax = Axis(f[1, 1])

	df = df_a_low_small

	scatter!(df[!, :ground_truth_mass_large], df[!, :calculated_mass_large], label="Large Inserts")
	scatter!(df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium], label="Medium Inserts")
	scatter!(df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], label="Small Inserts", color=:red)
	lines!([0, 100], [0, 100], label="Unity")
	
	ax.title = "Mass vs Known Mass (Small Patient)"
	ax.ylabel = "Calculated Mass (mg)"
	ax.xlabel = "Known Mass (mg)"
	ax.xticks = [0, 5, 10, 15, 20]
	ax.yticks = [0, 5, 10, 15, 20]

	xlims!(ax, 0, 20)
	ylims!(ax, 0, 20)
	
	f[1, 2] = Legend(f, ax, framevisible = false)
	f
end

# ╔═╡ 60001e02-4def-4d06-a5ff-03f65687f251
let
	f = Figure()
	ax = Axis(f[1, 1])

	df = df_a_low_medium

	scatter!(df[!, :ground_truth_mass_large], df[!, :calculated_mass_large], label="Large Inserts")
	scatter!(df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium], label="Medium Inserts")
	scatter!(df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], label="Small Inserts", color=:red)
	lines!([0, 100], [0, 100], label="Unity")
	
	ax.title = "Mass vs Known Mass (Medium Patient)"
	ax.ylabel = "Calculated Mass (mg)"
	ax.xlabel = "Known Mass (mg)"
	ax.xticks = [0, 5, 10, 15, 20]
	ax.yticks = [0, 5, 10, 15, 20]

	xlims!(ax, 0, 20)
	ylims!(ax, 0, 20)
	
	f[1, 2] = Legend(f, ax, framevisible = false)
	f
end

# ╔═╡ 514ca198-2b3c-40c4-be80-634bfbf70e57
let
	f = Figure()
	ax = Axis(f[1, 1])

	df = df_a_low_large

	scatter!(df[!, :ground_truth_mass_large], df[!, :calculated_mass_large], label="Large Inserts")
	scatter!(df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium], label="Medium Inserts")
	scatter!(df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], label="Small Inserts", color=:red)
	lines!([0, 100], [0, 100], label="Unity")
	
	ax.title = "Mass vs Known Mass (Large Patient)"
	ax.ylabel = "Calculated Mass (mg)"
	ax.xlabel = "Known Mass (mg)"
	ax.xticks = [0, 5, 10, 15, 20]
	ax.yticks = [0, 5, 10, 15, 20]

	xlims!(ax, 0, 20)
	ylims!(ax, 0, 20)
	
	f[1, 2] = Legend(f, ax, framevisible = false)
	f
end

# ╔═╡ 5b758329-4f92-432e-9ebd-f6881479d9f1
let
	df = df_a_low_small
	global agat_low_small_rms
	agat_low_small_rms = (rmsd(df[!, :calculated_mass_large], df[!, :ground_truth_mass_large]) + rmsd(df[!, :calculated_mass_medium], df[!, :ground_truth_mass_medium]) + rmsd(df[!, :calculated_mass_small], df[!, :ground_truth_mass_small]))/ 3
end

# ╔═╡ 06ff4fee-8f8e-4714-9275-ee922737960d
let
	df = df_a_low_medium
	global agat_low_medium_rms
	agat_low_medium_rms = (rmsd(df[!, :calculated_mass_large], df[!, :ground_truth_mass_large]) + rmsd(df[!, :calculated_mass_medium], df[!, :ground_truth_mass_medium]) + rmsd(df[!, :calculated_mass_small], df[!, :ground_truth_mass_small]))/ 3
end

# ╔═╡ 85fcfae8-f03f-4540-97b1-b62889bbd4c6
let
	df = df_a_low_large
	global agat_low_large_rms
	agat_low_large_rms = (rmsd(df[!, :calculated_mass_large], df[!, :ground_truth_mass_large]) + rmsd(df[!, :calculated_mass_medium], df[!, :ground_truth_mass_medium]) + rmsd(df[!, :calculated_mass_small], df[!, :ground_truth_mass_small])) / 3
end

# ╔═╡ 1316f1e8-1654-43db-840a-33fd91414f57
df_a_rms_low = DataFrame(
	small_rms = agat_low_small_rms,
	medium_rms = agat_low_medium_rms,
	large_rms = agat_low_large_rms
)

# ╔═╡ f3938476-0ce7-4fe1-be05-2bf6d1aa193b
let
	f = Figure()
	df = df_a_rms_low
	
	ax = Makie.Axis(
		f[1, 1],
		xticks = (1:3, ["Small Patient", "Medium Patient", "Large Patient"])
	)
	ax.title = "RMS"
	ax.ylabel = "RMSD"
	barplot!(ax, [1], df[!, :small_rms])
	barplot!(ax, [2], df[!, :medium_rms])
	barplot!(ax, [3], df[!, :large_rms])
	f
end

# ╔═╡ 02ea26a4-b582-4837-bc67-1a34d356a452
md"""
## Spatially Weighted
"""

# ╔═╡ 20445589-aeb7-4725-a23c-1f5a697cca75
path_swcs = "/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/output/swcs";

# ╔═╡ ef9c691b-6326-4044-9345-587c6ca2604a
df_s = CSV.read(string(path_swcs, "/full.csv"), DataFrame);

# ╔═╡ 4cfaa7fd-be47-4a5c-b7ce-fe8a03f284e9
df_s_low, df_s_normal = groupby(df_s, :DENSITY);

# ╔═╡ 487ab047-c95c-4342-87e2-d62894041650
df_s_low_small, df_s_low_medium, df_s_low_large = groupby(df_s_low, :SIZE);

# ╔═╡ 7a3047d1-0332-44fe-8316-252a5bdea26a
df_s_normal_small, df_s_normal_medium, df_s_normal_large = groupby(df_s_normal, :SIZE);

# ╔═╡ 1e769016-54e3-4859-9bc2-f315ee0c31d6
md"""
### Normal Density
"""

# ╔═╡ d11d45d7-8336-4b95-bfb5-68cd27d9d6d3
let
	f = Figure()
	ax = Axis(f[1, 1])

	df = df_s_normal_small

	scatter!(df[!, :ground_truth_mass_large], df[!, :calculated_swcs_large], label="Large Inserts")
	scatter!(df[!, :ground_truth_mass_medium], df[!, :calculated_swcs_medium], label="Medium Inserts")
	scatter!(df[!, :ground_truth_mass_small], df[!, :calculated_swcs_small], label="Small Inserts", color=:red)
	# lines!([0, 100], [0, 100], label="Unity")
	
	ax.title = "Mass vs Known Mass (Small Patient)"
	ax.ylabel = "Calculated Mass (mg)"
	ax.xlabel = "Known Mass (mg)"
	# ax.xticks = [0, 25, 50, 75, 100]
	# ax.yticks = [0, 25, 50, 75, 100]

	xlims!(ax, 0, 150)
	ylims!(ax, 0, 400)
	
	f[1, 2] = Legend(f, ax, framevisible = false)
	f
end

# ╔═╡ a7977d83-0502-4c0a-a8f1-57d637620aa4
let
	f = Figure()
	ax = Axis(f[1, 1])

	df = df_s_normal_medium

	scatter!(df[!, :ground_truth_mass_large], df[!, :calculated_swcs_large], label="Large Inserts")
	scatter!(df[!, :ground_truth_mass_medium], df[!, :calculated_swcs_medium], label="Medium Inserts")
	scatter!(df[!, :ground_truth_mass_small], df[!, :calculated_swcs_small], label="Small Inserts", color=:red)
	# lines!([0, 100], [0, 100], label="Unity")
	
	ax.title = "Mass vs Known Mass (Medium Patient)"
	ax.ylabel = "Calculated Mass (mg)"
	ax.xlabel = "Known Mass (mg)"
	# ax.xticks = [0, 25, 50, 75, 100]
	# ax.yticks = [0, 25, 50, 75, 100]

	xlims!(ax, 0, 150)
	ylims!(ax, 0, 400)
	
	f[1, 2] = Legend(f, ax, framevisible = false)
	f
end

# ╔═╡ 80ccbf94-47a6-42d2-9896-76882e8ed746
let
	f = Figure()
	ax = Axis(f[1, 1])

	df = df_s_normal_large

	scatter!(df[!, :ground_truth_mass_large], df[!, :calculated_swcs_large], label="Large Inserts")
	scatter!(df[!, :ground_truth_mass_medium], df[!, :calculated_swcs_medium], label="Medium Inserts")
	scatter!(df[!, :ground_truth_mass_small], df[!, :calculated_swcs_small], label="Small Inserts", color=:red)
	# lines!([0, 100], [0, 100], label="Unity")
	
	ax.title = "Mass vs Known Mass (Large Patient)"
	ax.ylabel = "Calculated Mass (mg)"
	ax.xlabel = "Known Mass (mg)"
	# ax.xticks = [0, 25, 50, 75, 100]
	# ax.yticks = [0, 25, 50, 75, 100]
	
	xlims!(ax, 0, 150)
	ylims!(ax, 0, 400)
	
	f[1, 2] = Legend(f, ax, framevisible = false)
	f
end

# ╔═╡ 848db47c-c578-4fd4-85ce-db926e369854
md"""
### -- Small Inserts
"""

# ╔═╡ 3a2ca3b6-8112-406e-8efb-6015ad5400bb
let
	f = Figure()
	ax = Axis(f[1, 1])

	df = df_s_normal_large
	scatter!(df[!, :ground_truth_mass_small], df[!, :calculated_swcs_small], label="Small Inserts", color=:red)
	
	ax.title = "Mass vs Known Mass (Large Patient)"
	ax.ylabel = "Calculated Mass (mg)"
	ax.xlabel = "Known Mass (mg)"

	xlims!(ax, 0, 10)
	ylims!(ax, 0, 70)
	
	f[1, 2] = Legend(f, ax, framevisible = false)
	f
end

# ╔═╡ c56ab225-b627-4b5f-9f71-e5cc37b19bd3
md"""
### Low Density
"""

# ╔═╡ f471c3aa-2c6d-43ab-bcfd-6841beef961b
let
	f = Figure()
	ax = Axis(f[1, 1])

	df = df_s_low_small

	scatter!(df[!, :ground_truth_mass_large], df[!, :calculated_swcs_large], label="Large Inserts")
	scatter!(df[!, :ground_truth_mass_medium], df[!, :calculated_swcs_medium], label="Medium Inserts")
	scatter!(df[!, :ground_truth_mass_small], df[!, :calculated_swcs_small], label="Small Inserts", color=:red)
	# lines!([0, 100], [0, 100], label="Unity")
	
	ax.title = "Mass vs Known Mass (Small Patient)"
	ax.ylabel = "Calculated Mass (mg)"
	ax.xlabel = "Known Mass (mg)"
	# ax.xticks = [0, 25, 50, 75, 100]
	# ax.yticks = [0, 25, 50, 75, 100]

	xlims!(ax, 0, 30)
	ylims!(ax, 0, 200)
	
	f[1, 2] = Legend(f, ax, framevisible = false)
	f
end

# ╔═╡ caf5c810-30e9-45ca-9103-8b83c3c9fdb4
let
	f = Figure()
	ax = Axis(f[1, 1])

	df = df_s_low_medium

	scatter!(df[!, :ground_truth_mass_large], df[!, :calculated_swcs_large], label="Large Inserts")
	scatter!(df[!, :ground_truth_mass_medium], df[!, :calculated_swcs_medium], label="Medium Inserts")
	scatter!(df[!, :ground_truth_mass_small], df[!, :calculated_swcs_small], label="Small Inserts", color=:red)
	# lines!([0, 100], [0, 100], label="Unity")
	
	ax.title = "Mass vs Known Mass (Medium Patient)"
	ax.ylabel = "Calculated Mass (mg)"
	ax.xlabel = "Known Mass (mg)"
	# ax.xticks = [0, 25, 50, 75, 100]
	# ax.yticks = [0, 25, 50, 75, 100]

	xlims!(ax, 0, 30)
	ylims!(ax, 0, 200)
	
	f[1, 2] = Legend(f, ax, framevisible = false)
	f
end

# ╔═╡ 7e9151dd-61c9-4e62-aec8-75c3be9e7caa
let
	f = Figure()
	ax = Axis(f[1, 1])

	df = df_s_low_large

	scatter!(df[!, :ground_truth_mass_large], df[!, :calculated_swcs_large], label="Large Inserts")
	scatter!(df[!, :ground_truth_mass_medium], df[!, :calculated_swcs_medium], label="Medium Inserts")
	scatter!(df[!, :ground_truth_mass_small], df[!, :calculated_swcs_small], label="Small Inserts", color=:red)
	# lines!([0, 100], [0, 100], label="Unity")
	
	ax.title = "Mass vs Known Mass (Large Patient)"
	ax.ylabel = "Calculated Mass (mg)"
	ax.xlabel = "Known Mass (mg)"
	# ax.xticks = [0, 25, 50, 75, 100]
	# ax.yticks = [0, 25, 50, 75, 100]
	
	xlims!(ax, 0, 30)
	ylims!(ax, 0, 200)
	
	f[1, 2] = Legend(f, ax, framevisible = false)
	f
end

# ╔═╡ 635fa6ec-67d1-42e9-b169-f4193a18ba3e
md"""
### -- Small Inserts
"""

# ╔═╡ da3b005a-1fe1-4928-9100-1d7601fb64e2
let
	f = Figure()
	ax = Axis(f[1, 1])

	df = df_s_low_large
	scatter!(df[!, :ground_truth_mass_small], df[!, :calculated_swcs_small], label="Small Inserts", color=:red)
	
	ax.title = "Mass vs Known Mass (Large Patient)"
	ax.ylabel = "Calculated Mass (mg)"
	ax.xlabel = "Known Mass (mg)"

	xlims!(ax, 0, 1)
	ylims!(ax, 0, 45)
	
	f[1, 2] = Legend(f, ax, framevisible = false)
	f
end

# ╔═╡ ef6aefe6-b83c-43ef-97fe-f3dc26e0584b
md"""
## Comparisons
"""

# ╔═╡ cab17220-a118-439e-9a88-1c1282671dfa
let
	f = Figure()
	
	ax = Makie.Axis(
		f[1, 1],
		xticks = (1:3, ["Small Patient", "Medium Patient", "Large Patient"]),
	)
	ax.title = "RMS Agatston vs Integrated (Normal Density)"
	ax.ylabel = "RMSD"
	labels = ["Agatston", "Integrated"]

	table = [1, 1, 2, 2, 3, 3]
	grp = [1, 2, 1, 2, 1, 2]
	heights = [agat_normal_small_rms, int_normal_small_rms, agat_normal_medium_rms, int_normal_medium_rms, agat_normal_large_rms, int_normal_medium_rms]
	barplot!(ax, table, heights, dodge=grp, color=grp)

	# Legend(f[1,2], heights, labels)

	f
end

# ╔═╡ 35d60a2a-9db2-4520-a192-e37e4360226f
let
	f = Figure()
	
	ax = Makie.Axis(
		f[1, 1],
		xticks = (1:3, ["Small Patient", "Medium Patient", "Large Patient"]),
	)
	ax.title = "RMS Agatston vs Integrated (Low Density)"
	ax.ylabel = "RMSD"
	labels = ["Agatston", "Integrated"]

	table = [1, 1, 2, 2, 3, 3]
	grp = [1, 2, 1, 2, 1, 2]
	heights = [agat_low_small_rms, int_low_small_rms, agat_low_medium_rms, int_low_medium_rms, agat_low_large_rms, int_low_medium_rms]
	barplot!(ax, table, heights, dodge=grp, color=grp)

	# Legend(f[1,2], heights, labels)

	f
end

# ╔═╡ 9dc00f8c-cf9e-4c2b-a36e-323c121110c1
md"""
### Normal Density
"""

# ╔═╡ 04285fb0-2921-42c3-8192-4e9edab4ed20
md"""
### --SWCS
"""

# ╔═╡ 961f3fa6-255f-481b-bf6e-1a4e84c534d2
let
	f = Figure()
	ax = Axis(f[1, 1])

	df1 = df_s_normal_large
	df2 = df_s_normal_medium
	df3 = df_s_normal_small
	scatter!(df1[!, :ground_truth_mass_small], df1[!, :calculated_swcs_small], label="Large Patient", color=:red)
	scatter!(df2[!, :ground_truth_mass_small], df2[!, :calculated_swcs_small], label="Medium Patient", color=:blue)
	scatter!(df3[!, :ground_truth_mass_small], df3[!, :calculated_swcs_small], label="Small Patient", color=:green)
	
	ax.title = "Small Inserts"
	ax.ylabel = "SWCS"
	ax.xlabel = "Known Mass (mg)"

	xlims!(ax, 0, 7)
	ylims!(ax, 0, 60)
	
	f[1, 2] = Legend(f, ax, framevisible = false)
	f
end

# ╔═╡ cc501607-a04f-4bc5-a620-7ba888d15e20
md"""
### --Integrated Score
"""

# ╔═╡ 2ad10af0-28bd-4ce0-a879-dce357ac926c
let
	f = Figure()
	ax = Axis(f[1, 1])

	df1 = df_i_normal_large
	df2 = df_i_normal_medium
	df3 = df_i_normal_small
	scatter!(df1[!, :ground_truth_mass_small], df1[!, :calculated_mass_small], label="Large Patient", color=:red)
	scatter!(df2[!, :ground_truth_mass_small], df2[!, :calculated_mass_small], label="Medium Patient", color=:blue)
	scatter!(df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small], label="Small Patient", color=:green)
	
	ax.title = "Small Inserts"
	ax.ylabel = "Integrated Mass Score"
	ax.xlabel = "Known Mass (mg)"

	xlims!(ax, 0, 7)
	ylims!(ax, 0, 7)
	
	f[1, 2] = Legend(f, ax, framevisible = false)
	f
end

# ╔═╡ 2a1bf4ad-984e-4edc-855d-97e6968bb459
md"""
### --Agatston Score
"""

# ╔═╡ 2e53194a-b88b-4723-85cc-585335e95e46
let
	f = Figure()
	ax = Axis(f[1, 1])

	df1 = df_a_normal_large
	df2 = df_a_normal_medium
	df3 = df_a_normal_small
	scatter!(df1[!, :ground_truth_mass_small], df1[!, :calculated_mass_small], label="Large Patient", color=:red)
	scatter!(df2[!, :ground_truth_mass_small], df2[!, :calculated_mass_small], label="Medium Patient", color=:blue)
	scatter!(df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small], label="Small Patient", color=:green)
	
	ax.title = "Small Inserts"
	ax.ylabel = "Agatston Mass Score"
	ax.xlabel = "Known Mass (mg)"

	xlims!(ax, 0, 7)
	ylims!(ax, 0, 15)
	
	f[1, 2] = Legend(f, ax, framevisible = false)
	f
end

# ╔═╡ 35ceff8f-7ca4-41f3-94f3-7ef17c1a2695
md"""
### Low Density
"""

# ╔═╡ 942a1d3a-74d0-45ce-a103-be1183a04d22
md"""
### --SWCS
"""

# ╔═╡ 651ee7a0-ed7d-4941-9ac0-3d88c046695f
let
	f = Figure()
	ax = Axis(f[1, 1])

	df1 = df_s_low_large
	df2 = df_s_low_medium
	df3 = df_s_low_small
	scatter!(df1[!, :ground_truth_mass_small], df1[!, :calculated_swcs_small], label="Small Inserts", color=:red)
	scatter!(df2[!, :ground_truth_mass_small], df2[!, :calculated_swcs_small], label="Small Inserts", color=:blue)
	scatter!(df3[!, :ground_truth_mass_small], df3[!, :calculated_swcs_small], label="Small Inserts", color=:green)
	
	ax.title = "Small Inserts (All Patient Size, Low Intensities)"
	ax.ylabel = "SWCS"
	ax.xlabel = "Known Mass (mg)"

	xlims!(ax, 0, 1)
	ylims!(ax, 0, 45)
	
	f[1, 2] = Legend(f, ax, framevisible = false)
	f
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
# ╟─a09621a9-ea7a-4722-9b96-136ca1114edf
# ╟─2425b571-6931-4c97-9998-7e382cd5b941
# ╠═51a79c37-dbbd-41af-b3ca-47cb67fa0ae9
# ╠═d114d2db-d813-4aff-b993-5df344e2317f
# ╠═42f95b3b-016b-4bd1-a139-0bdfbdbbae38
# ╠═b0a5ec8a-3630-43fe-aeb4-2ae2d5f66d55
# ╟─3e321850-eb7c-4620-99b4-60025921c10b
# ╟─3776dcd4-e1eb-4f0a-ab80-5f770ae66d83
# ╟─31933102-a867-4d9a-9af2-5cdc8597025d
# ╟─5c663e87-ad90-4f73-8a6f-ca62a522dab9
# ╟─21f9fc96-3560-4d9e-bbf8-947fd6d3a7c8
# ╟─4ec7c0c9-9bfb-4ef1-bf47-33b17bb4c3b6
# ╟─4a65be29-d4e7-43e9-bba3-cf80304ab5f3
# ╠═abc35a1a-4aad-4a3b-9f27-5ef8b3f6d8a3
# ╠═083eb4f4-a60f-4dd6-9ec6-fdf6ff15a697
# ╠═5aefde17-282b-4b0c-a08e-75e17f59c26f
# ╠═dae44dae-9d36-4591-8f64-34bc6d531524
# ╟─3c1510f9-b62b-48f1-9670-fbbb0d2bbc78
# ╟─14391d86-a2e2-4a6a-98f5-eef12be4502f
# ╠═67780e42-059d-40fe-8116-f9140909cc0a
# ╠═8144a8c1-1f9a-4de1-8bf2-55d443147540
# ╠═a20354d4-37d0-4e2d-a5f4-581a302fde54
# ╠═6b923331-0784-48cf-b4ca-9684ecb15335
# ╠═817510ce-847c-4c17-95c5-b5445f38e379
# ╠═dfb96fe4-2cb7-449d-bb68-b5e20299f4e3
# ╟─c3a13ff7-72e3-4233-8530-36dab53e4c94
# ╟─fdc99587-f3e5-42cd-a387-c36d0adf8eed
# ╟─7ae72022-136c-4594-87b1-94d2102e59dd
# ╟─ba3b4f2e-d0e0-47c9-ba61-4741cf296333
# ╠═6bbbae82-ebcc-455f-a035-9049291ce94f
# ╠═9c5d2ad9-3a12-4816-ae2b-e11621f8a783
# ╠═2ab0eeea-79b2-4d18-bc85-91f415a5f3c6
# ╠═fc4c940c-ed2d-45f9-9989-eae7979992a9
# ╟─6830ae86-54d5-4238-84d4-4f81541b4714
# ╟─ff238353-1281-40a3-9484-5f7d3177b7d6
# ╟─72d80487-acce-4b8d-afa3-28f30ad884ab
# ╟─60001e02-4def-4d06-a5ff-03f65687f251
# ╟─514ca198-2b3c-40c4-be80-634bfbf70e57
# ╠═5b758329-4f92-432e-9ebd-f6881479d9f1
# ╠═06ff4fee-8f8e-4714-9275-ee922737960d
# ╠═85fcfae8-f03f-4540-97b1-b62889bbd4c6
# ╠═1316f1e8-1654-43db-840a-33fd91414f57
# ╟─f3938476-0ce7-4fe1-be05-2bf6d1aa193b
# ╟─02ea26a4-b582-4837-bc67-1a34d356a452
# ╠═20445589-aeb7-4725-a23c-1f5a697cca75
# ╠═ef9c691b-6326-4044-9345-587c6ca2604a
# ╠═4cfaa7fd-be47-4a5c-b7ce-fe8a03f284e9
# ╠═487ab047-c95c-4342-87e2-d62894041650
# ╠═7a3047d1-0332-44fe-8316-252a5bdea26a
# ╟─1e769016-54e3-4859-9bc2-f315ee0c31d6
# ╟─d11d45d7-8336-4b95-bfb5-68cd27d9d6d3
# ╟─a7977d83-0502-4c0a-a8f1-57d637620aa4
# ╟─80ccbf94-47a6-42d2-9896-76882e8ed746
# ╟─848db47c-c578-4fd4-85ce-db926e369854
# ╟─3a2ca3b6-8112-406e-8efb-6015ad5400bb
# ╟─c56ab225-b627-4b5f-9f71-e5cc37b19bd3
# ╟─f471c3aa-2c6d-43ab-bcfd-6841beef961b
# ╟─caf5c810-30e9-45ca-9103-8b83c3c9fdb4
# ╟─7e9151dd-61c9-4e62-aec8-75c3be9e7caa
# ╟─635fa6ec-67d1-42e9-b169-f4193a18ba3e
# ╟─da3b005a-1fe1-4928-9100-1d7601fb64e2
# ╟─ef6aefe6-b83c-43ef-97fe-f3dc26e0584b
# ╟─cab17220-a118-439e-9a88-1c1282671dfa
# ╟─35d60a2a-9db2-4520-a192-e37e4360226f
# ╟─9dc00f8c-cf9e-4c2b-a36e-323c121110c1
# ╟─04285fb0-2921-42c3-8192-4e9edab4ed20
# ╟─961f3fa6-255f-481b-bf6e-1a4e84c534d2
# ╟─cc501607-a04f-4bc5-a620-7ba888d15e20
# ╟─2ad10af0-28bd-4ce0-a879-dce357ac926c
# ╟─2a1bf4ad-984e-4edc-855d-97e6968bb459
# ╟─2e53194a-b88b-4723-85cc-585335e95e46
# ╟─35ceff8f-7ca4-41f3-94f3-7ef17c1a2695
# ╟─942a1d3a-74d0-45ce-a103-be1183a04d22
# ╟─651ee7a0-ed7d-4941-9ac0-3d88c046695f

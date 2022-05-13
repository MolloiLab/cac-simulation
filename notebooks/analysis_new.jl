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

# ╔═╡ 6258d79a-3564-40d0-b7c5-c05cd1fa7691
df_i_low

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

# ╔═╡ d1e626e7-688a-4d8e-bb3b-53dedbfc4d47
df_a_low

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

# ╔═╡ a4df3923-548b-4818-b614-7d819012c5a7
md"""
## Integrated New
"""

# ╔═╡ d74881c2-a303-4d4a-8410-7be4a02ca1a1
df_i2 = CSV.read(string(path_integrated, "/full_new.csv"), DataFrame);

# ╔═╡ 6dbf443f-78ed-4e0e-b9a5-8f2449920760
df_i2_low, df_i2_normal = groupby(df_i2, :DENSITY);

# ╔═╡ 27aff353-492d-40ae-bf41-e00109b93b82
df_i2_low_small, df_i2_low_medium, df_i2_low_large = groupby(df_i2_low, :SIZE);

# ╔═╡ c1035609-c201-4ad3-aae1-c51f9daf7111
df_i2_normal_small, df_i2_normal_medium, df_i2_normal_large = groupby(df_i2_normal, :SIZE);

# ╔═╡ 481303d3-af6d-4a2a-9ade-f66ea00803b3
md"""
### Normal Density
"""

# ╔═╡ bfa68345-1fa6-4972-927a-f401d7ea6ece
let
	f = Figure()
	ax = Axis(f[1, 1])

	df = df_i2_normal_small

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

# ╔═╡ aa87744a-fe7e-432a-a916-1630218b390c
let
	f = Figure()
	ax = Axis(f[1, 1])

	df = df_i2_normal_medium

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

# ╔═╡ 9fccc335-1884-40eb-b3a2-0b3bac3ee727
let
	f = Figure()
	ax = Axis(f[1, 1])

	df = df_i2_normal_large

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

# ╔═╡ e8ee1296-f7c7-4d5f-9716-9eaa0c923435
md"""
### -- Small Inserts
"""

# ╔═╡ 255a5a5c-f9cb-4f80-b7a1-19cfb8e8f7ab
let
	f = Figure()
	ax = Axis(f[1, 1])

	df = df_i2_normal_large
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

# ╔═╡ 1178556d-1f11-426b-a3a1-213a9678a513
let
	df = df_i2_normal_small
	global int2_normal_small_rms
	int2_normal_small_rms = (rmsd(df[!, :calculated_mass_large], df[!, :ground_truth_mass_large]) + rmsd(df[!, :calculated_mass_medium], df[!, :ground_truth_mass_medium]) + rmsd(df[!, :calculated_mass_small], df[!, :ground_truth_mass_small]))/ 3
end

# ╔═╡ 555b9ad1-2625-494f-aa01-4cb25afa7c60
let
	df = df_i2_normal_medium
	global int2_normal_medium_rms
	int2_normal_medium_rms = (rmsd(df[!, :calculated_mass_large], df[!, :ground_truth_mass_large]) + rmsd(df[!, :calculated_mass_medium], df[!, :ground_truth_mass_medium]) + rmsd(df[!, :calculated_mass_small], df[!, :ground_truth_mass_small]))/ 3
end

# ╔═╡ 5a4b88cb-2fa3-44d4-9bcc-dd2e75a7db3a
let
	df = df_i2_normal_large
	global int2_normal_large_rms
	int2_normal_large_rms = (rmsd(df[!, :calculated_mass_large], df[!, :ground_truth_mass_large]) + rmsd(df[!, :calculated_mass_medium], df[!, :ground_truth_mass_medium]) + rmsd(df[!, :calculated_mass_small], df[!, :ground_truth_mass_small])) / 3
end

# ╔═╡ 8eef26bf-97f0-428a-8a00-5d1f9691b1ee
df_i2_rms_normal = DataFrame(
	small_rms = int2_normal_small_rms,
	medium_rms = int2_normal_medium_rms,
	large_rms = int2_normal_large_rms
)

# ╔═╡ 6663d028-42de-45b6-870a-3fe8d51a9cda
let
	f = Figure()
	df = df_i2_rms_normal
	
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

# ╔═╡ 8aa6961e-04e6-41b3-afa3-0e6cbd3ec897
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
	heights = [int_normal_small_rms, int2_normal_small_rms, int_normal_medium_rms, int2_normal_medium_rms, int_normal_medium_rms, int2_normal_medium_rms]
	barplot!(ax, table, heights, dodge=grp, color=grp)

	# Legend(f[1,2], heights, labels)

	f
end

# ╔═╡ c104d140-e410-4f0b-b33e-b3828b5c03e5
md"""
### Low Density
"""

# ╔═╡ ff0bb9ee-9d2a-478f-b7bc-20cf0565f017
let
	f = Figure()
	ax = Axis(f[1, 1])

	df = df_i2_low_small

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

# ╔═╡ 22c8739d-b262-4b15-bc43-04f790cb469e
let
	f = Figure()
	ax = Axis(f[1, 1])

	df = df_i2_low_medium

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

# ╔═╡ d815abc5-f695-4be3-b615-51be6c02b442
let
	f = Figure()
	ax = Axis(f[1, 1])

	df = df_i2_low_large

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

# ╔═╡ e5261f46-dce6-45cc-baa8-39c389e190ad
md"""
### -- Small Inserts
"""

# ╔═╡ fa349c10-3cd9-47dd-a290-129328b70d32
let
	f = Figure()
	ax = Axis(f[1, 1])

	df = df_i2_low_large
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

# ╔═╡ 5537da9f-0361-4f16-b681-2ec251251034
let
	df = df_i2_low_small
	global int2_low_small_rms
	int2_low_small_rms = (rmsd(df[!, :calculated_mass_large], df[!, :ground_truth_mass_large]) + rmsd(df[!, :calculated_mass_medium], df[!, :ground_truth_mass_medium]) + rmsd(df[!, :calculated_mass_small], df[!, :ground_truth_mass_small]))/ 3
end

# ╔═╡ 3fa49994-e141-4f70-a076-da43874d24cb
let
	df = df_i2_low_medium
	global int2_low_medium_rms
	int2_low_medium_rms = (rmsd(df[!, :calculated_mass_large], df[!, :ground_truth_mass_large]) + rmsd(df[!, :calculated_mass_medium], df[!, :ground_truth_mass_medium]) + rmsd(df[!, :calculated_mass_small], df[!, :ground_truth_mass_small]))/ 3
end

# ╔═╡ eb0cb4a2-3387-4eee-8e9e-6c691971f8a5
let
	df = df_i2_low_large
	global int2_low_large_rms
	int2_low_large_rms = (rmsd(df[!, :calculated_mass_large], df[!, :ground_truth_mass_large]) + rmsd(df[!, :calculated_mass_medium], df[!, :ground_truth_mass_medium]) + rmsd(df[!, :calculated_mass_small], df[!, :ground_truth_mass_small])) / 3
end

# ╔═╡ 5f33c24c-70ea-421e-a7ec-357328a337e5
df_i2_rms_low = DataFrame(
	small_rms = int_low_small_rms,
	medium_rms = int_low_medium_rms,
	large_rms = int_low_large_rms
)

# ╔═╡ 524361c2-d710-474f-aa84-74b19564d811
let
	f = Figure()
	df = df_i2_rms_low
	
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

# ╔═╡ cbff013e-48a6-462e-90a1-35bcdb78e574
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
	heights = [int_low_small_rms, int2_low_small_rms, int_low_medium_rms, int2_low_medium_rms, int_low_medium_rms, int2_low_medium_rms]
	barplot!(ax, table, heights, dodge=grp, color=grp)

	# Legend(f[1,2], heights, labels)

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
# ╠═6258d79a-3564-40d0-b7c5-c05cd1fa7691
# ╠═d1e626e7-688a-4d8e-bb3b-53dedbfc4d47
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
# ╟─ef6aefe6-b83c-43ef-97fe-f3dc26e0584b
# ╟─cab17220-a118-439e-9a88-1c1282671dfa
# ╟─35d60a2a-9db2-4520-a192-e37e4360226f
# ╟─a4df3923-548b-4818-b614-7d819012c5a7
# ╠═d74881c2-a303-4d4a-8410-7be4a02ca1a1
# ╠═6dbf443f-78ed-4e0e-b9a5-8f2449920760
# ╠═27aff353-492d-40ae-bf41-e00109b93b82
# ╠═c1035609-c201-4ad3-aae1-c51f9daf7111
# ╟─481303d3-af6d-4a2a-9ade-f66ea00803b3
# ╟─bfa68345-1fa6-4972-927a-f401d7ea6ece
# ╟─aa87744a-fe7e-432a-a916-1630218b390c
# ╟─9fccc335-1884-40eb-b3a2-0b3bac3ee727
# ╟─e8ee1296-f7c7-4d5f-9716-9eaa0c923435
# ╟─255a5a5c-f9cb-4f80-b7a1-19cfb8e8f7ab
# ╠═1178556d-1f11-426b-a3a1-213a9678a513
# ╠═555b9ad1-2625-494f-aa01-4cb25afa7c60
# ╠═5a4b88cb-2fa3-44d4-9bcc-dd2e75a7db3a
# ╠═8eef26bf-97f0-428a-8a00-5d1f9691b1ee
# ╟─6663d028-42de-45b6-870a-3fe8d51a9cda
# ╟─8aa6961e-04e6-41b3-afa3-0e6cbd3ec897
# ╟─c104d140-e410-4f0b-b33e-b3828b5c03e5
# ╟─ff0bb9ee-9d2a-478f-b7bc-20cf0565f017
# ╟─22c8739d-b262-4b15-bc43-04f790cb469e
# ╟─d815abc5-f695-4be3-b615-51be6c02b442
# ╟─e5261f46-dce6-45cc-baa8-39c389e190ad
# ╟─fa349c10-3cd9-47dd-a290-129328b70d32
# ╠═5537da9f-0361-4f16-b681-2ec251251034
# ╠═3fa49994-e141-4f70-a076-da43874d24cb
# ╠═eb0cb4a2-3387-4eee-8e9e-6c691971f8a5
# ╠═5f33c24c-70ea-421e-a7ec-357328a337e5
# ╟─524361c2-d710-474f-aa84-74b19564d811
# ╟─cbff013e-48a6-462e-90a1-35bcdb78e574

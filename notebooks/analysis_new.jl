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
	global int_normal_small_rmsd
	int_normal_small_rmsd = (rmsd(df[!, :calculated_mass_large], df[!, :ground_truth_mass_large]) + rmsd(df[!, :calculated_mass_medium], df[!, :ground_truth_mass_medium]) + rmsd(df[!, :calculated_mass_small], df[!, :ground_truth_mass_small]))/ 3
end

# ╔═╡ d114d2db-d813-4aff-b993-5df344e2317f
let
	df = df_i_normal_medium
	global int_normal_medium_rmsd
	int_normal_medium_rmsd = (rmsd(df[!, :calculated_mass_large], df[!, :ground_truth_mass_large]) + rmsd(df[!, :calculated_mass_medium], df[!, :ground_truth_mass_medium]) + rmsd(df[!, :calculated_mass_small], df[!, :ground_truth_mass_small]))/ 3
end

# ╔═╡ 42f95b3b-016b-4bd1-a139-0bdfbdbbae38
let
	df = df_i_normal_large
	global int_normal_large_rmsd
	int_normal_large_rmsd = (rmsd(df[!, :calculated_mass_large], df[!, :ground_truth_mass_large]) + rmsd(df[!, :calculated_mass_medium], df[!, :ground_truth_mass_medium]) + rmsd(df[!, :calculated_mass_small], df[!, :ground_truth_mass_small])) / 3
end

# ╔═╡ b0a5ec8a-3630-43fe-aeb4-2ae2d5f66d55
df_i_rmsd_normal = DataFrame(
	small_rmsd = int_normal_small_rmsd,
	medium_rmsd = int_normal_medium_rmsd,
	large_rmsd = int_normal_large_rmsd
)

# ╔═╡ 3e321850-eb7c-4620-99b4-60025921c10b
let
	f = Figure()
	df = df_i_rmsd_normal
	
	ax = Makie.Axis(
		f[1, 1],
		xticks = (1:3, ["Small Patient", "Medium Patient", "Large Patient"]),
	)
	ax.title = "RMS (Large)"
	ax.ylabel = "RMSD"
	
	barplot!(ax, [1], df[!, :small_rmsd])
	barplot!(ax, [2], df[!, :medium_rmsd])
	barplot!(ax, [3], df[!, :large_rmsd])

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
	
	ax.title = "Spatially Weighted (Small Patient)"
	ax.ylabel = "SWCS"
	ax.xlabel = "Known Mass (mg)"
	# ax.xticks = [0, 25, 50, 75, 100]
	# ax.yticks = [0, 25, 50, 75, 100]

	xlims!(ax, 0, 150)
	ylims!(ax, 0, 500)
	
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
	
	ax.title = "Spatially Weighted (Medium Patient)"
	ax.ylabel = "SWCS"
	ax.xlabel = "Known Mass (mg)"
	# ax.xticks = [0, 25, 50, 75, 100]
	# ax.yticks = [0, 25, 50, 75, 100]

	xlims!(ax, 0, 150)
	ylims!(ax, 0, 500)
	
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
	# lines!([0, 500], [0, 500], label="Unity")
	
	ax.title = "Spatially Weighted (Large Patient)"
	ax.ylabel = "SWCS"
	ax.xlabel = "Known Mass (mg)"
	# ax.xticks = [0, 25, 50, 75, 100]
	# ax.yticks = [0, 25, 50, 75, 100]
	
	xlims!(ax, 0, 150)
	ylims!(ax, 0, 500)
	
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

	# xlims!(ax, 0, 10)
	# ylims!(ax, 0, 70)
	
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

	
	colors = Makie.wong_colors()
	labels = ["Agatston", "Integrated"]
	elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]
	title = "Groups"
	
	ax.title = "RMS Agatston vs Integrated (Normal Density)"
	ax.ylabel = "RMSD"
	labels = ["Agatston", "Integrated"]

	table = [1, 1, 2, 2, 3, 3]
	grp = [1, 2, 1, 2, 1, 2]
	heights = [agat_normal_small_rms, int_normal_small_rmsd, agat_normal_medium_rms, int_normal_medium_rmsd, agat_normal_large_rms, int_normal_medium_rmsd]
	barplot!(ax, table, heights, dodge=grp, color=colors[grp])

	Legend(f[1, 2], elements, labels, title)

	f
end

# ╔═╡ 35d60a2a-9db2-4520-a192-e37e4360226f
let
	f = Figure()
	
	ax = Makie.Axis(
		f[1, 1],
		xticks = (1:3, ["Small Patient", "Medium Patient", "Large Patient"]),
	)


	colors = Makie.wong_colors()
	labels = ["Agatston", "Integrated"]
	elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]
	title = "Groups"
	
	ax.title = "RMS Agatston vs Integrated (Low Density)"
	ax.ylabel = "RMSD"
	
	table = [1, 1, 2, 2, 3, 3]
	grp = [1, 2, 1, 2, 1, 2]
	heights = [agat_low_small_rms, int_low_small_rms, agat_low_medium_rms, int_low_medium_rms, agat_low_large_rms, int_low_medium_rms]
	barplot!(ax, table, heights, dodge=grp, color=colors[grp])
	
	Legend(f[1, 2], elements, labels, title)

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

	xlims!(ax, 0, 8)
	ylims!(ax, 0, 200)
	
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
	scatter!(df1[!, :ground_truth_mass_small], df1[!, :calculated_swcs_small], label="Large Patient", color=:red)
	scatter!(df2[!, :ground_truth_mass_small], df2[!, :calculated_swcs_small], label="Medium Patient", color=:blue)
	scatter!(df3[!, :ground_truth_mass_small], df3[!, :calculated_swcs_small], label="Small Patient", color=:green)
	
	ax.title = "Small Inserts"
	ax.ylabel = "SWCS"
	ax.xlabel = "Known Mass (mg)"

	xlims!(ax, 0, 1)
	ylims!(ax, 0, 100)
	
	f[1, 2] = Legend(f, ax, framevisible = false)
	f
end

# ╔═╡ 72dca916-fcec-41e6-928f-f4cd524ba49a
md"""
### --Integrated Score
"""

# ╔═╡ f12f6cd7-4cfc-46ed-bdd3-10186547e6f0
let
	f = Figure()
	ax = Axis(f[1, 1])

	df1 = df_i_low_large
	df2 = df_i_low_medium
	df3 = df_i_low_small
	scatter!(df1[!, :ground_truth_mass_small], df1[!, :calculated_mass_small], label="Large Patient", color=:red)
	scatter!(df2[!, :ground_truth_mass_small], df2[!, :calculated_mass_small], label="Medium Patient", color=:blue)
	scatter!(df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small], label="Small Patient", color=:green)
	
	ax.title = "Small Inserts"
	ax.ylabel = "Integrated Mass"
	ax.xlabel = "Known Mass (mg)"

	xlims!(ax, 0, 1)
	ylims!(ax, -0.1, 1)
	
	f[1, 2] = Legend(f, ax, framevisible = false)
	f
end

# ╔═╡ a3c8b60e-dad5-4eb2-ac85-91e49fdfc0f8
md"""
### --Agatston Score
"""

# ╔═╡ 1565f83e-a82c-4a70-97d3-efc7c88fff45
let
	f = Figure()
	ax = Axis(f[1, 1])

	df1 = df_a_low_large
	df2 = df_a_low_medium
	df3 = df_a_low_small
	scatter!(df1[!, :ground_truth_mass_small], df1[!, :calculated_mass_small], label="Large Patient", color=:red)
	scatter!(df2[!, :ground_truth_mass_small], df2[!, :calculated_mass_small], label="Medium Patient", color=:blue)
	scatter!(df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small], label="Small Patient", color=:green)
	
	ax.title = "Small Inserts"
	ax.ylabel = "Agatston Mass"
	ax.xlabel = "Known Mass (mg)"

	xlims!(ax, 0, 1)
	ylims!(ax, -0.1, 1)
	
	f[1, 2] = Legend(f, ax, framevisible = false)
	f
end

# ╔═╡ 990071b7-d417-463e-ae39-8f1591df0d7f
md"""
# Figures
"""

# ╔═╡ 19016ebb-93c3-470b-93ae-4930e04a9a29
md"""
## Figure 1. Zero CAC Scores
"""

# ╔═╡ c58fd175-a769-44ad-9e0b-c5e99b3e5bf7
md"""
#### SWCS
"""

# ╔═╡ a03ee3a8-1c80-47e1-b9a9-f6245e2259ae
df_s_80, df_s_100, df_s_120, df_s_135, = groupby(df_s, :scan);

# ╔═╡ f26ebcf5-bb40-4b9a-aaeb-66fcd8e28d66
begin
	mean_swcs_80, std_swcs_80= mean(df_s_80[!, :swcs_bkg]), std(df_s_80[!, :swcs_bkg])
	mean_swcs_100, std_swcs_100= mean(df_s_100[!, :swcs_bkg]), std(df_s_100[!, :swcs_bkg])
	mean_swcs_120, std_swcs_120= mean(df_s_120[!, :swcs_bkg]), std(df_s_120[!, :swcs_bkg])
	mean_swcs_135, std_swcs_135 = mean(df_s_135[!, :swcs_bkg]), std(df_s_135[!, :swcs_bkg])
end

# ╔═╡ 0eeaab7f-6dca-4096-82f0-421e0f5d3815
begin
	array_s_80 = Array(df_s_80[!, 12:14])
	array_s_100 = Array(df_s_100[!, 12:14])
	array_s_120 = Array(df_s_120[!, 12:14])
	array_s_135 = Array(df_s_135[!, 12:14])
end;

# ╔═╡ 4c2cda1f-4bfb-4ac0-9b8b-387045647148
begin
	num_zeroCAC_80 = length(findall(x -> x < mean_swcs_80 - std_swcs_80, array_s_80))
	num_zeroCAC_100 = length(findall(x -> x < mean_swcs_100 - std_swcs_100, array_s_100))
	num_zeroCAC_120 = length(findall(x -> x < mean_swcs_120 - std_swcs_120, array_s_120))
	num_zeroCAC_135 = length(findall(x -> x < mean_swcs_135 - std_swcs_135, array_s_135))

	total_zero_s = num_zeroCAC_80 + num_zeroCAC_100 + num_zeroCAC_120 + num_zeroCAC_135
	total_cac = length(array_s_80) * 4
end;

# ╔═╡ e522127d-aa19-46c9-9617-2b2c8574455d
total_cac

# ╔═╡ 9c595b59-39a2-42c6-a1fc-081ca6a91df1
total_zero_s

# ╔═╡ 90573c0e-ee3a-4b75-9d92-1677e644df08
md"""
#### Agatston
"""

# ╔═╡ 5e0ecc11-1ce7-4542-ab45-4c1da3484b32
array_a = hcat(df_a[!, 7], df_a[!, 9], df_a[!, 11]);

# ╔═╡ 89246a65-ef7f-4697-875c-eb0da534d557
num_zero_a = length(findall(x -> x == 0, array_a))

# ╔═╡ 9be8af6c-c35d-4824-8041-9a49c763ebf1
md"""
#### Integrated
"""

# ╔═╡ 9541d8c2-1b5b-4fe2-bec4-754a5a431b40
df_i_80, df_i_100, df_i_120, df_i_135, = groupby(df_i, :scan);

# ╔═╡ 89f6ae83-2c91-456e-807c-9585e3f029e6
begin
	mean_i_80, std_i_80= mean(df_i_80[!, :mass_bkg]), std(df_i_80[!, :mass_bkg])
	mean_i_100, std_i_100= mean(df_i_100[!, :mass_bkg]), std(df_i_100[!, :mass_bkg])
	mean_i_120, std_i_120= mean(df_i_120[!, :mass_bkg]), std(df_i_120[!, :mass_bkg])
	mean_i_135, std_i_135 = mean(df_i_135[!, :mass_bkg]), std(df_i_135[!, :mass_bkg])
end

# ╔═╡ ad26cc0b-d5f7-45fb-bb30-a9c9a5901c34
begin
	array_i_80 = hcat(df_i_80[!, 7], df_i_80[!, 9], df_i_80[!, 11])
	array_i_100 = hcat(df_i_100[!, 7], df_i_100[!, 9], df_i_100[!, 11])
	array_i_120 = hcat(df_i_120[!, 7], df_i_120[!, 9], df_i_120[!, 11])
	array_i_135 = hcat(df_i_135[!, 7], df_i_135[!, 9], df_i_135[!, 11])
end;

# ╔═╡ 16db3a02-7614-4bcb-9eb7-43f695cfb81e
begin
	num_zeroCAC_80_i = length(findall(x -> x < mean_i_80 - std_i_80, array_i_80))
	num_zeroCAC_100_i = length(findall(x -> x < mean_i_100 - std_i_100, array_i_100))
	num_zeroCAC_120_i = length(findall(x -> x < mean_i_120 - std_i_120, array_i_120))
	num_zeroCAC_135_i = length(findall(x -> x < mean_i_135 - std_i_135, array_i_135))

	total_zero_i = num_zeroCAC_80_i + num_zeroCAC_100_i + num_zeroCAC_120_i + num_zeroCAC_135_i
end;

# ╔═╡ 05d6a990-55de-48e4-b620-eaf326e3a208
total_zero_i

# ╔═╡ bd0b16ef-ea0f-4028-9d1c-3cf7e69d7c43


# ╔═╡ 54d07958-b2c9-41d7-973a-44bef3067b93
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

# ╔═╡ 1f08025c-f326-4d7f-b031-4dad002daa14
md"""
## Figure 2. Pearson Correlation Coefficient
"""

# ╔═╡ ce58f6f6-c0e4-41f5-a77b-5f7a01351446
md"""
#### SWCS
"""

# ╔═╡ 7c46403c-21fd-41d2-80c5-e42168cf7c72
array_s = Array(hcat(df_s[!, end-3:end-1]));

# ╔═╡ c1f2bfb6-5860-4af2-9926-b6acd00189d7
vec_s = vec(array_s)

# ╔═╡ 87ec1484-9e5b-41e5-a7a1-94d4803188b4
md"""
#### Agatston
"""

# ╔═╡ 3897a45c-e7d2-4406-9a4e-a4dd427ce5cb
vec_a = vec(array_a)

# ╔═╡ 1f6ff3a6-489e-4149-9eab-d9473a9f15ff
md"""
#### Integrated
"""

# ╔═╡ 2829e8a3-7eab-4c2b-9745-b01c5362428a
vec_i = vec(array_i)

# ╔═╡ 2ac67a93-383d-4420-a5b7-0a10a0bfa3d6
md"""
#### Known mass
"""

# ╔═╡ 00389fea-d2da-49a4-b5c4-daca3f79857f
array_k = hcat(df_a[!, 6], df_a[!, 8], df_a[!, 10]);

# ╔═╡ 7627848b-4421-4002-afce-313244b56697
vec_k = vec(array_k)

# ╔═╡ 0ae2a51f-ec55-4fc9-a2b6-a68d75745045
md"""
#### Correlations
"""

# ╔═╡ 98dbcf02-7703-402d-a307-65aa7a0c653f
pearson_s = cor(hcat(vec_s, vec_k))

# ╔═╡ 06033041-bf24-4f8f-b7da-1642d27c31a3
pearson_a = cor(hcat(vec_a, vec_k))

# ╔═╡ 48882170-9362-4990-b001-8f616395f92f
pearson_i = cor(hcat(vec_i, vec_k))

# ╔═╡ 0130ad0f-9aad-4f2c-88b1-83f25ff5fcfa
p_s = EqualVarianceTTest(vec_s, vec_k)

# ╔═╡ 37f9c0f0-0be3-4628-941f-f26e0dd2ece2
p_i = EqualVarianceTTest(vec_i, vec_k)

# ╔═╡ 0194aa4a-e363-45b2-a37a-78cc7e9be9b0
p_a = EqualVarianceTTest(vec_a, vec_k)

# ╔═╡ a2d49d55-a715-4623-87f3-e6d9079e1c13
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

# ╔═╡ 80c3eb0e-20ed-46a7-a112-8562d79eb514
md"""
## Figure 3. Linear Regression Mass Score
"""

# ╔═╡ 547dad9f-2856-45e6-aa66-eddb2041f82a
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

# ╔═╡ 2a84d850-0d54-4a0a-9be3-d03a7fa16ba7
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

# ╔═╡ b90e187e-e5a7-4e1a-8944-1b8e5e6fb96f
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

# ╔═╡ f66f1f95-f35f-414f-b4d6-193822f16376
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

# ╔═╡ 4d13c232-a4ae-4193-9697-9c93c1446af2
# let
# 	df = df_i_low
# 	# gt_array = Array(df[!, :ground_truth_mass_large])
# 	gt_array = vec(hcat(df[!, :ground_truth_mass_large], df[!, :ground_truth_mass_medium], df[!, :ground_truth_mass_small]))
# 	# calc_array = Array(df[!, :calculated_mass_large])
# 	calc_array = vec(hcat(df[!, :calculated_mass_large], df[!, :calculated_mass_medium], df[!, :calculated_mass_small]))
# 	data = DataFrame(
# 		X = gt_array,
# 		Y= calc_array
# 	)
# 	model = lm(@formula(Y ~ X), data)
# 	r2 = r2(model)
# end

# ╔═╡ 37d0e6f9-2c73-43bd-88cb-a20583483abc
# let
# 	df = df_a_low
# 	# gt_array = Array(df[!, :ground_truth_mass_large])
# 	gt_array = vec(hcat(df[!, :ground_truth_mass_large], df[!, :ground_truth_mass_medium], df[!, :ground_truth_mass_small]))
# 	# calc_array = Array(df[!, :calculated_mass_large])
# 	calc_array = vec(hcat(df[!, :calculated_mass_large], df[!, :calculated_mass_medium], df[!, :calculated_mass_small]))
# 	data = DataFrame(
# 		X = gt_array,
# 		Y= calc_array
# 	)
# 	model = lm(@formula(Y ~ X), data)
# 	r2 = r2(model)
# end

# ╔═╡ 94dc1c51-561c-4404-a91a-0721f4065f66
md"""
## Figure 4. RMSE Mass Score
"""

# ╔═╡ c19df90c-41ab-4dd7-accb-ab059ddb8195
begin
	agat_arr = [agat_normal_small_rms, agat_low_small_rms, agat_normal_medium_rms, agat_low_medium_rms, agat_normal_large_rms, agat_low_large_rms]
	int_arr = [int_normal_small_rmsd, int_low_small_rms, int_normal_medium_rmsd, int_low_medium_rms, int_normal_large_rmsd, int_low_large_rms]
end;

# ╔═╡ 46b4d4aa-b9bc-4557-a713-52ce2191e8d3
total_rms_agat = mean(agat_arr)

# ╔═╡ 9fb19a68-7fe4-4117-8892-206f1fa4eaa1
total_rms_int = mean(int_arr)

# ╔═╡ c0e8f5af-c0e2-4c7f-9c94-3d580ce2520f
val_rms = EqualVarianceTTest(agat_arr, int_arr)

# ╔═╡ e1e24e06-7063-4851-afd2-b25dd414b4d7
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

# ╔═╡ 19049b16-bacc-49b9-8533-546deea85775
md"""
## TEMPLATE FIGURE
"""

# ╔═╡ ea678826-dc11-4c22-9b84-17b21ac59002
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
# ╟─72dca916-fcec-41e6-928f-f4cd524ba49a
# ╟─f12f6cd7-4cfc-46ed-bdd3-10186547e6f0
# ╟─a3c8b60e-dad5-4eb2-ac85-91e49fdfc0f8
# ╟─1565f83e-a82c-4a70-97d3-efc7c88fff45
# ╟─990071b7-d417-463e-ae39-8f1591df0d7f
# ╟─19016ebb-93c3-470b-93ae-4930e04a9a29
# ╟─c58fd175-a769-44ad-9e0b-c5e99b3e5bf7
# ╠═a03ee3a8-1c80-47e1-b9a9-f6245e2259ae
# ╠═f26ebcf5-bb40-4b9a-aaeb-66fcd8e28d66
# ╠═0eeaab7f-6dca-4096-82f0-421e0f5d3815
# ╠═4c2cda1f-4bfb-4ac0-9b8b-387045647148
# ╠═e522127d-aa19-46c9-9617-2b2c8574455d
# ╠═9c595b59-39a2-42c6-a1fc-081ca6a91df1
# ╟─90573c0e-ee3a-4b75-9d92-1677e644df08
# ╠═5e0ecc11-1ce7-4542-ab45-4c1da3484b32
# ╠═89246a65-ef7f-4697-875c-eb0da534d557
# ╟─9be8af6c-c35d-4824-8041-9a49c763ebf1
# ╠═9541d8c2-1b5b-4fe2-bec4-754a5a431b40
# ╠═89f6ae83-2c91-456e-807c-9585e3f029e6
# ╠═ad26cc0b-d5f7-45fb-bb30-a9c9a5901c34
# ╠═16db3a02-7614-4bcb-9eb7-43f695cfb81e
# ╠═05d6a990-55de-48e4-b620-eaf326e3a208
# ╠═bd0b16ef-ea0f-4028-9d1c-3cf7e69d7c43
# ╟─54d07958-b2c9-41d7-973a-44bef3067b93
# ╟─1f08025c-f326-4d7f-b031-4dad002daa14
# ╟─ce58f6f6-c0e4-41f5-a77b-5f7a01351446
# ╠═7c46403c-21fd-41d2-80c5-e42168cf7c72
# ╠═c1f2bfb6-5860-4af2-9926-b6acd00189d7
# ╟─87ec1484-9e5b-41e5-a7a1-94d4803188b4
# ╠═3897a45c-e7d2-4406-9a4e-a4dd427ce5cb
# ╟─1f6ff3a6-489e-4149-9eab-d9473a9f15ff
# ╠═2829e8a3-7eab-4c2b-9745-b01c5362428a
# ╟─2ac67a93-383d-4420-a5b7-0a10a0bfa3d6
# ╠═00389fea-d2da-49a4-b5c4-daca3f79857f
# ╠═7627848b-4421-4002-afce-313244b56697
# ╟─0ae2a51f-ec55-4fc9-a2b6-a68d75745045
# ╠═98dbcf02-7703-402d-a307-65aa7a0c653f
# ╠═06033041-bf24-4f8f-b7da-1642d27c31a3
# ╠═48882170-9362-4990-b001-8f616395f92f
# ╠═0130ad0f-9aad-4f2c-88b1-83f25ff5fcfa
# ╠═37f9c0f0-0be3-4628-941f-f26e0dd2ece2
# ╠═0194aa4a-e363-45b2-a37a-78cc7e9be9b0
# ╟─a2d49d55-a715-4623-87f3-e6d9079e1c13
# ╟─80c3eb0e-20ed-46a7-a112-8562d79eb514
# ╟─547dad9f-2856-45e6-aa66-eddb2041f82a
# ╠═2a84d850-0d54-4a0a-9be3-d03a7fa16ba7
# ╠═b90e187e-e5a7-4e1a-8944-1b8e5e6fb96f
# ╠═f66f1f95-f35f-414f-b4d6-193822f16376
# ╟─4d13c232-a4ae-4193-9697-9c93c1446af2
# ╟─37d0e6f9-2c73-43bd-88cb-a20583483abc
# ╟─94dc1c51-561c-4404-a91a-0721f4065f66
# ╠═c19df90c-41ab-4dd7-accb-ab059ddb8195
# ╠═46b4d4aa-b9bc-4557-a713-52ce2191e8d3
# ╠═9fb19a68-7fe4-4117-8892-206f1fa4eaa1
# ╠═c0e8f5af-c0e2-4c7f-9c94-3d580ce2520f
# ╟─e1e24e06-7063-4851-afd2-b25dd414b4d7
# ╟─19049b16-bacc-49b9-8533-546deea85775
# ╟─ea678826-dc11-4c22-9b84-17b21ac59002

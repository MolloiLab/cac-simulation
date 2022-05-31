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

# ╔═╡ bab5cbcc-2193-49c1-8baa-f1e3d83ebaeb
md"""
### Normal Density
"""

# ╔═╡ 05074eb9-996b-40e9-8c66-c56562d54f84
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

# ╔═╡ 1874695d-4aba-443c-b71e-c071d49025ef
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

# ╔═╡ a97dbbed-31d2-44a8-b7f2-1d1b59e3f493
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

# ╔═╡ a4c255ee-1f96-466a-87a6-6963f03da0b3
md"""
### -- Small Inserts
"""

# ╔═╡ 3a6a0850-62c0-48c9-a793-901267ea7912
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

# ╔═╡ 7f103d79-215d-4a5c-9ee3-1afc90eae60d
let
	df = df_i_normal_small
	global int_normal_small_rmsd
	int_normal_small_rmsd = (rmsd(df[!, :calculated_mass_large], df[!, :ground_truth_mass_large]) + rmsd(df[!, :calculated_mass_medium], df[!, :ground_truth_mass_medium]) + rmsd(df[!, :calculated_mass_small], df[!, :ground_truth_mass_small]))/ 3
end

# ╔═╡ 1f5eae73-2a34-4266-8e9d-3c1c448920f8
let
	df = df_i_normal_medium
	global int_normal_medium_rmsd
	int_normal_medium_rmsd = (rmsd(df[!, :calculated_mass_large], df[!, :ground_truth_mass_large]) + rmsd(df[!, :calculated_mass_medium], df[!, :ground_truth_mass_medium]) + rmsd(df[!, :calculated_mass_small], df[!, :ground_truth_mass_small]))/ 3
end

# ╔═╡ eb82b841-2942-4830-a3eb-4c339c0a3f81
let
	df = df_i_normal_large
	global int_normal_large_rmsd
	int_normal_large_rmsd = (rmsd(df[!, :calculated_mass_large], df[!, :ground_truth_mass_large]) + rmsd(df[!, :calculated_mass_medium], df[!, :ground_truth_mass_medium]) + rmsd(df[!, :calculated_mass_small], df[!, :ground_truth_mass_small])) / 3
end

# ╔═╡ 0b64b3cd-ce34-4842-90f6-7135bf9f0029
df_i_rmsd_normal = DataFrame(
	small_rmsd = int_normal_small_rmsd,
	medium_rmsd = int_normal_medium_rmsd,
	large_rmsd = int_normal_large_rmsd
)

# ╔═╡ 1dcd69da-1781-436d-9f54-be6d1599f073
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

# ╔═╡ cb049113-fd9c-4ce0-a326-44ba65dbf735
md"""
### Low Density
"""

# ╔═╡ 39328eb2-6fe9-4156-a51b-845b6621443a
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

# ╔═╡ 27e0fabf-6d24-4335-b04b-32e4096885a8
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

# ╔═╡ 0e2f7272-84d0-4f91-ae20-85814af65199
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

# ╔═╡ b29500f3-38e4-48c2-ac8f-bbe7650c364b
md"""
### -- Small Inserts
"""

# ╔═╡ 83c3d582-abaf-49cf-a7e5-7a06cc322cbf
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
	ylims!(ax, -1.5, 2.5)
	
	f[1, 2] = Legend(f, ax, framevisible = false)
	f
end

# ╔═╡ d13b9731-0a9e-4995-bde0-d369c2e8e3be
let
	df = df_i_low_small
	global int_low_small_rms
	int_low_small_rms = (rmsd(df[!, :calculated_mass_large], df[!, :ground_truth_mass_large]) + rmsd(df[!, :calculated_mass_medium], df[!, :ground_truth_mass_medium]) + rmsd(df[!, :calculated_mass_small], df[!, :ground_truth_mass_small]))/ 3
end

# ╔═╡ fd562c23-d505-4ec3-9d84-0dcf3c7c8d52
let
	df = df_i_low_medium
	global int_low_medium_rms
	int_low_medium_rms = (rmsd(df[!, :calculated_mass_large], df[!, :ground_truth_mass_large]) + rmsd(df[!, :calculated_mass_medium], df[!, :ground_truth_mass_medium]) + rmsd(df[!, :calculated_mass_small], df[!, :ground_truth_mass_small]))/ 3
end

# ╔═╡ 551528f7-2df3-46a6-be34-a959de286cb3
let
	df = df_i_low_large
	global int_low_large_rms
	int_low_large_rms = (rmsd(df[!, :calculated_mass_large], df[!, :ground_truth_mass_large]) + rmsd(df[!, :calculated_mass_medium], df[!, :ground_truth_mass_medium]) + rmsd(df[!, :calculated_mass_small], df[!, :ground_truth_mass_small])) / 3
end

# ╔═╡ b58b6d36-b170-4a2f-9937-c0f04ec0cd3e
df_i_rms_low = DataFrame(
	small_rms = int_low_small_rms,
	medium_rms = int_low_medium_rms,
	large_rms = int_low_large_rms
)

# ╔═╡ 7cdc22f2-1301-48ec-a62a-747dc0b549b1
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

# ╔═╡ de8c33f2-e75d-4753-a04c-c1f95f600382
df_a_low

# ╔═╡ ac2136d9-6100-4304-ad84-d7b6cd75e648
md"""
### Normal Density
"""

# ╔═╡ e5050b8a-2c93-4da7-a209-9b405566dfda
let
	f = Figure()
	ax = Axis(f[1, 1])

	df = df_a_normal_small

	scatter!(df[!, :ground_truth_mass_large], df[!, :calculated_mass_large], label="Large Inserts")
	scatter!(df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium], label="Medium Inserts")
	scatter!(df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], label="Small Inserts", color=:red)
	lines!([0, 200], [0, 200], label="Unity")
	
	ax.title = "Mass vs Known Mass (Small Patient)"
	ax.ylabel = "Calculated Mass (mg)"
	ax.xlabel = "Known Mass (mg)"
	# ax.xticks = [0, 25, 50, 75, 100]
	# ax.yticks = [0, 25, 50, 75, 100]

	# xlims!(ax, 0, 100)
	# ylims!(ax, 0, 100)
	
	f[1, 2] = Legend(f, ax, framevisible = false)
	f
end

# ╔═╡ 3602c899-8558-424c-91bb-c881ada941bf
let
	f = Figure()
	ax = Axis(f[1, 1])

	df = df_a_normal_medium

	scatter!(df[!, :ground_truth_mass_large], df[!, :calculated_mass_large], label="Large Inserts")
	scatter!(df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium], label="Medium Inserts")
	scatter!(df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], label="Small Inserts", color=:red)
	lines!([0, 200], [0, 200], label="Unity")
	
	ax.title = "Mass vs Known Mass (Medium Patient)"
	ax.ylabel = "Calculated Mass (mg)"
	ax.xlabel = "Known Mass (mg)"
	# ax.xticks = [0, 25, 50, 75, 100]
	# ax.yticks = [0, 25, 50, 75, 100]

	# xlims!(ax, 0, 100)
	# ylims!(ax, 0, 100)
	
	f[1, 2] = Legend(f, ax, framevisible = false)
	f
end

# ╔═╡ 73bc295c-47e6-42ae-96ed-abed9f4fcb0d
let
	f = Figure()
	ax = Axis(f[1, 1])

	df = df_a_normal_large

	scatter!(df[!, :ground_truth_mass_large], df[!, :calculated_mass_large], label="Large Inserts")
	scatter!(df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium], label="Medium Inserts")
	scatter!(df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], label="Small Inserts", color=:red)
	lines!([0, 200], [0, 200], label="Unity")
	
	ax.title = "Mass vs Known Mass (Large Patient)"
	ax.ylabel = "Calculated Mass (mg)"
	ax.xlabel = "Known Mass (mg)"
	# ax.xticks = [0, 25, 50, 75, 100]
	# ax.yticks = [0, 25, 50, 75, 100]

	# xlims!(ax, 0, 100)
	# ylims!(ax, 0, 100)
	
	f[1, 2] = Legend(f, ax, framevisible = false)
	f
end

# ╔═╡ 7ea77749-7cf4-430f-9ec3-c60c08003f19
let
	df = df_a_normal_small
	global agat_normal_small_rms
	agat_normal_small_rms = (rmsd(df[!, :calculated_mass_large], df[!, :ground_truth_mass_large]) + rmsd(df[!, :calculated_mass_medium], df[!, :ground_truth_mass_medium]) + rmsd(df[!, :calculated_mass_small], df[!, :ground_truth_mass_small]))/ 3
end

# ╔═╡ 3e00829f-9b69-41da-8c0d-8d230d00a647
let
	df = df_a_normal_medium
	global agat_normal_medium_rms
	agat_normal_medium_rms = (rmsd(df[!, :calculated_mass_large], df[!, :ground_truth_mass_large]) + rmsd(df[!, :calculated_mass_medium], df[!, :ground_truth_mass_medium]) + rmsd(df[!, :calculated_mass_small], df[!, :ground_truth_mass_small]))/ 3
end

# ╔═╡ fc587c42-2f3c-442f-a744-144b0150da87
let
	df = df_a_normal_large
	global agat_normal_large_rms
	agat_normal_large_rms = (rmsd(df[!, :calculated_mass_large], df[!, :ground_truth_mass_large]) + rmsd(df[!, :calculated_mass_medium], df[!, :ground_truth_mass_medium]) + rmsd(df[!, :calculated_mass_small], df[!, :ground_truth_mass_small])) / 3
end

# ╔═╡ b230c46b-f7b4-48e9-816a-d0b50d91f7a0
df_a_rms = DataFrame(
	normal_small_rms = agat_normal_small_rms,
	normal_medium_rms = agat_normal_medium_rms,
	normal_large_rms = agat_normal_large_rms
)

# ╔═╡ 787c27bf-250f-413d-98df-8a637001f522
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

# ╔═╡ c91266a8-ff76-43c9-abb9-548ae789b74e
md"""
### Low Density
"""

# ╔═╡ 431b5585-c6c6-4f67-a7ea-471a63892bba
let
	f = Figure()
	ax = Axis(f[1, 1])

	df = df_a_low_small

	scatter!(df[!, :ground_truth_mass_large], df[!, :calculated_mass_large], label="Large Inserts")
	scatter!(df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium], label="Medium Inserts")
	scatter!(df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], label="Small Inserts", color=:red)
	lines!([0, 30], [0, 30], label="Unity")
	
	ax.title = "Mass vs Known Mass (Small Patient)"
	ax.ylabel = "Calculated Mass (mg)"
	ax.xlabel = "Known Mass (mg)"
	# ax.xticks = [0, 5, 10, 15, 20]
	# ax.yticks = [0, 5, 10, 15, 20]

	# xlims!(ax, 0, 20)
	# ylims!(ax, 0, 20)
	
	f[1, 2] = Legend(f, ax, framevisible = false)
	f
end

# ╔═╡ 6c0744c3-9977-4094-be58-a47208c18356
let
	f = Figure()
	ax = Axis(f[1, 1])

	df = df_a_low_medium

	scatter!(df[!, :ground_truth_mass_large], df[!, :calculated_mass_large], label="Large Inserts")
	scatter!(df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium], label="Medium Inserts")
	scatter!(df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], label="Small Inserts", color=:red)
	lines!([0, 30], [0, 30], label="Unity")
	
	ax.title = "Mass vs Known Mass (Medium Patient)"
	ax.ylabel = "Calculated Mass (mg)"
	ax.xlabel = "Known Mass (mg)"
	# ax.xticks = [0, 5, 10, 15, 20]
	# ax.yticks = [0, 5, 10, 15, 20]

	# xlims!(ax, 0, 20)
	# ylims!(ax, 0, 20)
	
	f[1, 2] = Legend(f, ax, framevisible = false)
	f
end

# ╔═╡ 15161fba-5577-4ad3-bdfc-ef51a46a8728
let
	f = Figure()
	ax = Axis(f[1, 1])

	df = df_a_low_large

	scatter!(df[!, :ground_truth_mass_large], df[!, :calculated_mass_large], label="Large Inserts")
	scatter!(df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium], label="Medium Inserts")
	scatter!(df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], label="Small Inserts", color=:red)
	lines!([0, 80], [0, 80], label="Unity")
	
	ax.title = "Mass vs Known Mass (Large Patient)"
	ax.ylabel = "Calculated Mass (mg)"
	ax.xlabel = "Known Mass (mg)"
	# ax.xticks = [0, 5, 10, 15, 20]
	# ax.yticks = [0, 5, 10, 15, 20]

	# xlims!(ax, 0, 20)
	# ylims!(ax, 0, 20)
	
	f[1, 2] = Legend(f, ax, framevisible = false)
	f
end

# ╔═╡ 1a7e89a2-a8ca-48b3-ac8f-0ce4256a3bc3
let
	df = df_a_low_small
	global agat_low_small_rms
	agat_low_small_rms = (rmsd(df[!, :calculated_mass_large], df[!, :ground_truth_mass_large]) + rmsd(df[!, :calculated_mass_medium], df[!, :ground_truth_mass_medium]) + rmsd(df[!, :calculated_mass_small], df[!, :ground_truth_mass_small]))/ 3
end

# ╔═╡ 9b8d2206-0f33-4b2b-af6f-01547b1a8e23
let
	df = df_a_low_medium
	global agat_low_medium_rms
	agat_low_medium_rms = (rmsd(df[!, :calculated_mass_large], df[!, :ground_truth_mass_large]) + rmsd(df[!, :calculated_mass_medium], df[!, :ground_truth_mass_medium]) + rmsd(df[!, :calculated_mass_small], df[!, :ground_truth_mass_small]))/ 3
end

# ╔═╡ 21d86397-2984-4ad7-9886-81910749eae1
let
	df = df_a_low_large
	global agat_low_large_rms
	agat_low_large_rms = (rmsd(df[!, :calculated_mass_large], df[!, :ground_truth_mass_large]) + rmsd(df[!, :calculated_mass_medium], df[!, :ground_truth_mass_medium]) + rmsd(df[!, :calculated_mass_small], df[!, :ground_truth_mass_small])) / 3
end

# ╔═╡ 7ae079cd-333e-4143-9f31-116785181a67
df_a_rms_low = DataFrame(
	small_rms = agat_low_small_rms,
	medium_rms = agat_low_medium_rms,
	large_rms = agat_low_large_rms
)

# ╔═╡ c45ce6b0-1efd-4d72-8851-aff501131fd3
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

# ╔═╡ 9d88ebf5-300e-42ec-bf71-d50ae3e96868
df_s_low

# ╔═╡ aa7260ed-f91f-4b77-aa8f-05d082f4f28b
md"""
### Normal Density
"""

# ╔═╡ e9d4030b-f3de-4266-bff5-5497f10c67b6
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

# ╔═╡ dc945c2b-3799-4013-8740-758c7b478298
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

# ╔═╡ b33a3456-05ac-44bf-83e8-515af110bbf3
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

# ╔═╡ 5c89f8ce-7ce2-4738-b380-5c0a8a417bdc
md"""
### -- Small Inserts
"""

# ╔═╡ 5407c9ae-6300-4bbe-98e1-64091c7a9342
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

# ╔═╡ 1b1d7f6d-b251-47dd-bf83-74dde6c43614
md"""
### Low Density
"""

# ╔═╡ a894efec-d5b4-42a0-8542-f8e255b25180
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

	# xlims!(ax, 0, 30)
	# ylims!(ax, 0, 200)
	
	f[1, 2] = Legend(f, ax, framevisible = false)
	f
end

# ╔═╡ 8b00cc41-fcfd-422a-b4b9-4ffabcee5c10
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

	# xlims!(ax, 0, 30)
	# ylims!(ax, 0, 200)
	
	f[1, 2] = Legend(f, ax, framevisible = false)
	f
end

# ╔═╡ 32f8096d-672b-45a7-b81e-8ae8af4cd610
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
	
	# xlims!(ax, 0, 30)
	# ylims!(ax, 0, 200)
	
	f[1, 2] = Legend(f, ax, framevisible = false)
	f
end

# ╔═╡ 38d2e535-0eb4-4805-898a-49b597f7b00b
md"""
### -- Small Inserts
"""

# ╔═╡ 42bb0c22-db46-4d90-ab78-d4ddaf623560
let
	f = Figure()
	ax = Axis(f[1, 1])

	df = df_s_low_large
	scatter!(df[!, :ground_truth_mass_small], df[!, :calculated_swcs_small], label="Small Inserts", color=:red)
	
	ax.title = "Mass vs Known Mass (Large Patient)"
	ax.ylabel = "Calculated Mass (mg)"
	ax.xlabel = "Known Mass (mg)"

	# xlims!(ax, 0, 1)
	# ylims!(ax, 0, 45)
	
	f[1, 2] = Legend(f, ax, framevisible = false)
	f
end

# ╔═╡ 9b202cd7-0c47-43be-a703-30c239e6f3b3
md"""
## Comparisons
"""

# ╔═╡ 11b3703d-77b1-4a01-8080-4d826761691e
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

# ╔═╡ 16203ff4-e60d-4c59-88b2-f18d10b46a25
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

# ╔═╡ 6af2f1b2-452a-450b-a1ec-92c066f88abb
md"""
### Normal Density
"""

# ╔═╡ d444ca27-669a-4855-bb72-bdd466aded38
md"""
### --SWCS
"""

# ╔═╡ cd78f5a4-9cf8-4cea-ae91-f8bb70d04f73
let
	f = Figure()
	ax = Axis(f[1, 1])

	df1 = df_s_normal_large
	df2 = df_s_normal_medium
	df3 = df_s_normal_small
	scatter!(df1[!, :ground_truth_mass_small], df1[!, :calculated_swcs_small], label="Large Patient", color=:red)
	scatter!(df2[!, :ground_truth_mass_small], df2[!, :calculated_swcs_small], label="Medium Patient", color=:blue)
	scatter!(df3[!, :ground_truth_mass_small], df3[!, :calculated_swcs_small], label="Small Patient", color=:green)
	# lines!([0, 100], [0, 100], label="Unity")
	
	ax.title = "Small Inserts"
	ax.ylabel = "SWCS"
	ax.xlabel = "Known Mass (mg)"

	xlims!(ax, 0, 8)
	ylims!(ax, 0, 200)
	
	f[1, 2] = Legend(f, ax, framevisible = false)
	f
end

# ╔═╡ 30a5c203-dea7-4018-b2f8-b7662d5af43c
md"""
### --Integrated Score
"""

# ╔═╡ d173391f-441f-4d80-9b73-02111d3f6f7c
let
	f = Figure()
	ax = Axis(f[1, 1])

	df1 = df_i_normal_large
	df2 = df_i_normal_medium
	df3 = df_i_normal_small
	scatter!(df1[!, :ground_truth_mass_small], df1[!, :calculated_mass_small], label="Large Patient", color=:red)
	scatter!(df2[!, :ground_truth_mass_small], df2[!, :calculated_mass_small], label="Medium Patient", color=:blue)
	scatter!(df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small], label="Small Patient", color=:green)
	lines!([0, 100], [0, 100], label="Unity")
	
	ax.title = "Small Inserts"
	ax.ylabel = "Integrated Mass Score"
	ax.xlabel = "Known Mass (mg)"

	xlims!(ax, 0, 7)
	ylims!(ax, 0, 7)
	
	f[1, 2] = Legend(f, ax, framevisible = false)
	f
end

# ╔═╡ 31c907b5-0e94-4805-8684-521e07a6ecf5
md"""
### --Agatston Score
"""

# ╔═╡ 4798270f-2a98-4573-8f09-4cbbbff7e631
let
	f = Figure()
	ax = Axis(f[1, 1])

	df1 = df_a_normal_large
	df2 = df_a_normal_medium
	df3 = df_a_normal_small
	scatter!(df1[!, :ground_truth_mass_small], df1[!, :calculated_mass_small], label="Large Patient", color=:red)
	scatter!(df2[!, :ground_truth_mass_small], df2[!, :calculated_mass_small], label="Medium Patient", color=:blue)
	scatter!(df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small], label="Small Patient", color=:green)
	lines!([0, 100], [0, 100], label="Unity")
	
	ax.title = "Small Inserts"
	ax.ylabel = "Agatston Mass Score"
	ax.xlabel = "Known Mass (mg)"

	xlims!(ax, 0, 7)
	ylims!(ax, 0, 15)
	
	f[1, 2] = Legend(f, ax, framevisible = false)
	f
end

# ╔═╡ cf008537-8b66-48cb-b432-75c18c9b4af6
md"""
### Low Density
"""

# ╔═╡ 39899ec2-cde0-450f-92c3-2165f23dbe27
md"""
### --SWCS
"""

# ╔═╡ a526508c-5e0f-4211-99dc-17f951d0556d
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

# ╔═╡ 1f3f75e2-c680-4651-840e-adb6712ce9f9
md"""
### --Integrated Score
"""

# ╔═╡ e902c002-6181-4f26-b855-146616a34715
let
	f = Figure()
	ax = Axis(f[1, 1])

	df1 = df_i_low_large
	df2 = df_i_low_medium
	df3 = df_i_low_small
	scatter!(df1[!, :ground_truth_mass_small], df1[!, :calculated_mass_small], label="Large Patient", color=:red)
	scatter!(df2[!, :ground_truth_mass_small], df2[!, :calculated_mass_small], label="Medium Patient", color=:blue)
	scatter!(df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small], label="Small Patient", color=:green)
	lines!([0, 100], [0, 100], label="Unity")
	
	ax.title = "Small Inserts"
	ax.ylabel = "Integrated Mass"
	ax.xlabel = "Known Mass (mg)"

	xlims!(ax, 0, 1)
	ylims!(ax, -0.1, 1)
	
	f[1, 2] = Legend(f, ax, framevisible = false)
	f
end

# ╔═╡ 4208fa71-8441-4de8-95b3-1c075ab7c94c
md"""
### --Agatston Score
"""

# ╔═╡ 522a64dc-3c4e-4892-80a0-8e2383087e4b
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
total_zero_i

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
array_i = hcat(df_i[!, 7], df_i[!, 9], df_i[!, 11])

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

# ╔═╡ 564644cc-cfe7-40d5-a660-637e140a46cd
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

# ╔═╡ 2be9bf56-bd59-434d-bd99-3285d50771b8
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
df_i_r = CSV.read(string(path_integrated, "/full.csv"), DataFrame);

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
df_a_r = CSV.read(string(path_agat, "/full.csv"), DataFrame);

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

# ╔═╡ de038336-aa35-417c-a723-17d6745dca3d
df_s_r

# ╔═╡ a85f777a-76f2-4c64-9973-ea9dec245600
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
	# axbottom = Axis(gb[1, 1])
	
	# df2 = df_i_low
	# df2_r = df_i_low_r
	# scatter!(axbottom, df2_r[!, :calculated_mass_large], df2[!, :calculated_mass_large], label="Large Inserts")
	# scatter!(axbottom, df2_r[!, :calculated_mass_medium], df2[!, :calculated_mass_medium], label="Medium Inserts")
	# scatter!(axbottom, df2_r[!, :calculated_mass_small], df2[!, :calculated_mass_small], label="Small Inserts", color=:red)
	# lines!(axbottom, [-1000, 1000], [-1000, 1000], label="Unity")

	# xlims!(axbottom, low=0, high=25)
	# ylims!(axbottom, low=0, high=25)
	# axbottom.xticks = [0, 5, 10, 15, 20, 25]
	# axbottom.yticks = [0, 5, 10, 15, 20, 25]
	# axbottom.xlabel = "Known Mass (mg)"
	# axbottom.ylabel = "Calculated Mass (mg)"
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

	# ##-- D --##
	# axbottomright = Axis(gd[1, 1])
	
	# df4 = df_s
	# df4_r = df_s_r
	# scatter!(axbottomright, df4_r[!, :calculated_swcs_large], df4[!, :calculated_swcs_large], label="Large Inserts")
	# scatter!(axbottomright, df4_r[!, :calculated_swcs_medium], df4[!, :calculated_swcs_medium], label="Medium Inserts")
	# scatter!(axbottomright, df4_r[!, :calculated_swcs_small], df4[!, :calculated_swcs_small], label="Small Inserts", color=:red)
	# lines!(axbottomright, [-1000, 1000], [-1000, 1000], label="Unity")

	# xlims!(axbottomright, low=0, high=500)
	# ylims!(axbottomright, low=0, high=500)
	# # axbottomright.xticks = [0, 5, 10, 15, 20, 25]
	# # axbottomright.yticks = [0, 5, 10, 15, 20, 25]
	# axbottomright.xlabel = "SWCS 1"
	# axbottomright.ylabel = "SWCS 1"

	##-- LABELS --##

	f[1:2, 3] = Legend(f, axbottomright, framevisible = false)

	
	for (label, layout) in zip(["A", "B", "C"], [ga, gb, gc])
	    Label(layout[1, 1, TopLeft()], label,
	        textsize = 25,
	        padding = (0, 60, 25, 0),
	        halign = :right)
	end

	# colsize!(f.layout, 1, Auto(0.5))
	# rowsize!(gcd, 1, Auto(1.5))

	f
end

# ╔═╡ 40f443f8-6e6b-4de3-9c2e-b70599640c5d
md"""
## TEMPLATE FIGURE
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
# ╟─bab5cbcc-2193-49c1-8baa-f1e3d83ebaeb
# ╟─05074eb9-996b-40e9-8c66-c56562d54f84
# ╟─1874695d-4aba-443c-b71e-c071d49025ef
# ╟─a97dbbed-31d2-44a8-b7f2-1d1b59e3f493
# ╟─a4c255ee-1f96-466a-87a6-6963f03da0b3
# ╟─3a6a0850-62c0-48c9-a793-901267ea7912
# ╠═7f103d79-215d-4a5c-9ee3-1afc90eae60d
# ╠═1f5eae73-2a34-4266-8e9d-3c1c448920f8
# ╠═eb82b841-2942-4830-a3eb-4c339c0a3f81
# ╠═0b64b3cd-ce34-4842-90f6-7135bf9f0029
# ╟─1dcd69da-1781-436d-9f54-be6d1599f073
# ╟─cb049113-fd9c-4ce0-a326-44ba65dbf735
# ╟─39328eb2-6fe9-4156-a51b-845b6621443a
# ╟─27e0fabf-6d24-4335-b04b-32e4096885a8
# ╟─0e2f7272-84d0-4f91-ae20-85814af65199
# ╟─b29500f3-38e4-48c2-ac8f-bbe7650c364b
# ╟─83c3d582-abaf-49cf-a7e5-7a06cc322cbf
# ╠═d13b9731-0a9e-4995-bde0-d369c2e8e3be
# ╠═fd562c23-d505-4ec3-9d84-0dcf3c7c8d52
# ╠═551528f7-2df3-46a6-be34-a959de286cb3
# ╠═b58b6d36-b170-4a2f-9937-c0f04ec0cd3e
# ╟─7cdc22f2-1301-48ec-a62a-747dc0b549b1
# ╟─a58beea4-b186-44c6-bf75-98d368f7ebe4
# ╠═4fd5b649-cea1-4d1f-95ac-35aa12f8fbff
# ╠═3fa3748a-22d6-49d8-9888-a749dca99ce2
# ╠═c38b9994-7698-4c81-b65b-534b994d8ff0
# ╠═79952ee5-3ccd-4d85-95b4-f7610c90d63d
# ╠═7223ae32-7241-463f-be13-2d2ce6bc7de7
# ╠═de8c33f2-e75d-4753-a04c-c1f95f600382
# ╟─ac2136d9-6100-4304-ad84-d7b6cd75e648
# ╟─e5050b8a-2c93-4da7-a209-9b405566dfda
# ╟─3602c899-8558-424c-91bb-c881ada941bf
# ╟─73bc295c-47e6-42ae-96ed-abed9f4fcb0d
# ╠═7ea77749-7cf4-430f-9ec3-c60c08003f19
# ╠═3e00829f-9b69-41da-8c0d-8d230d00a647
# ╠═fc587c42-2f3c-442f-a744-144b0150da87
# ╠═b230c46b-f7b4-48e9-816a-d0b50d91f7a0
# ╟─787c27bf-250f-413d-98df-8a637001f522
# ╟─c91266a8-ff76-43c9-abb9-548ae789b74e
# ╟─431b5585-c6c6-4f67-a7ea-471a63892bba
# ╟─6c0744c3-9977-4094-be58-a47208c18356
# ╟─15161fba-5577-4ad3-bdfc-ef51a46a8728
# ╠═1a7e89a2-a8ca-48b3-ac8f-0ce4256a3bc3
# ╠═9b8d2206-0f33-4b2b-af6f-01547b1a8e23
# ╠═21d86397-2984-4ad7-9886-81910749eae1
# ╠═7ae079cd-333e-4143-9f31-116785181a67
# ╟─c45ce6b0-1efd-4d72-8851-aff501131fd3
# ╟─f39a9ed5-c897-4b90-8494-bf231312970e
# ╠═32623778-9a19-4318-81c6-1fc7aa0157fe
# ╠═317f8c26-7b39-45c3-b9d3-02925d4b0514
# ╠═06deefa7-9170-4e2f-ac49-d6dc65c0af76
# ╠═5a4cc121-fd02-49bf-bdf5-07b71b48ce19
# ╠═d5d8e9da-eba8-4788-b787-95b3a2b41329
# ╠═9d88ebf5-300e-42ec-bf71-d50ae3e96868
# ╟─aa7260ed-f91f-4b77-aa8f-05d082f4f28b
# ╟─e9d4030b-f3de-4266-bff5-5497f10c67b6
# ╟─dc945c2b-3799-4013-8740-758c7b478298
# ╟─b33a3456-05ac-44bf-83e8-515af110bbf3
# ╟─5c89f8ce-7ce2-4738-b380-5c0a8a417bdc
# ╟─5407c9ae-6300-4bbe-98e1-64091c7a9342
# ╟─1b1d7f6d-b251-47dd-bf83-74dde6c43614
# ╟─a894efec-d5b4-42a0-8542-f8e255b25180
# ╟─8b00cc41-fcfd-422a-b4b9-4ffabcee5c10
# ╟─32f8096d-672b-45a7-b81e-8ae8af4cd610
# ╟─38d2e535-0eb4-4805-898a-49b597f7b00b
# ╟─42bb0c22-db46-4d90-ab78-d4ddaf623560
# ╟─9b202cd7-0c47-43be-a703-30c239e6f3b3
# ╟─11b3703d-77b1-4a01-8080-4d826761691e
# ╟─16203ff4-e60d-4c59-88b2-f18d10b46a25
# ╟─6af2f1b2-452a-450b-a1ec-92c066f88abb
# ╟─d444ca27-669a-4855-bb72-bdd466aded38
# ╟─cd78f5a4-9cf8-4cea-ae91-f8bb70d04f73
# ╟─30a5c203-dea7-4018-b2f8-b7662d5af43c
# ╟─d173391f-441f-4d80-9b73-02111d3f6f7c
# ╟─31c907b5-0e94-4805-8684-521e07a6ecf5
# ╟─4798270f-2a98-4573-8f09-4cbbbff7e631
# ╟─cf008537-8b66-48cb-b432-75c18c9b4af6
# ╟─39899ec2-cde0-450f-92c3-2165f23dbe27
# ╟─a526508c-5e0f-4211-99dc-17f951d0556d
# ╟─1f3f75e2-c680-4651-840e-adb6712ce9f9
# ╟─e902c002-6181-4f26-b855-146616a34715
# ╟─4208fa71-8441-4de8-95b3-1c075ab7c94c
# ╟─522a64dc-3c4e-4892-80a0-8e2383087e4b
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
# ╠═d04ce8a0-8224-49c7-a057-f1c7d0d65ef8
# ╠═1ef8e352-7e89-461f-9ac8-8a39cf2c14d7
# ╠═c86c3d16-ace0-4de2-a7fb-267ea2925302
# ╠═d4c89efd-fb4b-4553-8920-847c393cb4bc
# ╠═564644cc-cfe7-40d5-a660-637e140a46cd
# ╠═2be9bf56-bd59-434d-bd99-3285d50771b8
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
# ╠═de038336-aa35-417c-a723-17d6745dca3d
# ╠═a85f777a-76f2-4c64-9973-ea9dec245600
# ╟─40f443f8-6e6b-4de3-9c2e-b70599640c5d
# ╠═f33bd1c6-009a-4071-9cb5-93f4f71c5415

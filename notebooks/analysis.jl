### A Pluto.jl notebook ###
# v0.19.2

using Markdown
using InteractiveUtils

# ╔═╡ c2496cdd-8ddb-4eeb-9747-f6b5397d4be6
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

# ╔═╡ 77c51949-937b-48fc-9621-5a071461f65d
TableOfContents()

# ╔═╡ 7bf86054-44de-456c-adc2-27617153faba
md"""
## Integrated Mass Scoring
"""

# ╔═╡ 41012935-cb85-4fc5-8f7c-7f4d27c2c7e4
canon_path = "/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/output/integrated_scoring";

# ╔═╡ be63d9b6-bf33-4e05-a3fa-51d77e5288c5
df_i_canon = CSV.read(string(canon_path, "/full.csv"), DataFrame);

# ╔═╡ b4a759df-05d7-47fb-a52a-e9848070da15
md"""
#### No Noise
"""

# ╔═╡ de94dfde-06a5-4991-9854-d286bd01c7bf
df_i_canon0 = df_i_canon[1:12, :];

# ╔═╡ 687297b7-6f31-41fb-9549-6262c27b5cdb
begin
	fInt1 = Figure()
	axInt1 = Axis(fInt1[1, 1])

	scatter!(df_i_canon0[!, :ground_truth_mass_large], df_i_canon0[!, :calculated_mass_large], label="Large Inserts")
	scatter!(df_i_canon0[!, :ground_truth_mass_medium], df_i_canon0[!, :calculated_mass_medium], label="Medium Inserts")
	scatter!(df_i_canon0[!, :ground_truth_mass_small], df_i_canon0[!, :calculated_mass_small], label="Small Inserts", color=:red)
	lines!([0, 100], [0, 100], label="Unity")
	
	axInt1.title = "IHU Mass vs Known Mass"
	axInt1.ylabel = "Calculated Mass (mg)"
	axInt1.xlabel = "Known Mass (mg)"
	axInt1.xticks = [0, 25, 50, 75, 100]
	axInt1.yticks = [0, 25, 50, 75, 100]

	xlims!(axInt1, 0, 100)
	ylims!(axInt1, 0, 100)
	
	fInt1[1, 2] = Legend(fInt1, axInt1, framevisible = false)
	fInt1
end

# ╔═╡ 5eeb159c-57b0-470b-81b4-493c2112314b
begin
	fInt1_small = Figure()
	axInt1_small = Axis(fInt1_small[1, 1])

	scatter!(df_i_canon0[!, :ground_truth_mass_small], df_i_canon0[!, :calculated_mass_small], label="Small Inserts", color=:red)
	lines!([0, 6], [0, 6], label="Unity")
	
	axInt1_small.title = "IHU Mass vs Known Mass (Small)"
	axInt1_small.ylabel = "Calculated Mass (mg)"
	axInt1_small.xlabel = "Known Mass (mg)"
	
	fInt1_small[1, 2] = Legend(fInt1_small, axInt1_small, framevisible = false)
	fInt1_small
end

# ╔═╡ 4977450f-882b-4913-b838-40fc69a66f4c
md"""
#### Medium Noise
"""

# ╔═╡ 0753c940-289c-4b91-8bfa-8544d553bed5
df_i_canon3 = df_i_canon[13:25, :];

# ╔═╡ f83b5d8f-7f03-4d1c-9ba0-8baf11e3d364
begin
	fInt3 = Figure()
	axInt3 = Axis(fInt3[1, 1])

	scatter!(df_i_canon3[!, :ground_truth_mass_large], df_i_canon3[!, :calculated_mass_large], label="Large Inserts")
	scatter!(df_i_canon3[!, :ground_truth_mass_medium], df_i_canon3[!, :calculated_mass_medium], label="Medium Inserts")
	scatter!(df_i_canon3[!, :ground_truth_mass_small], df_i_canon3[!, :calculated_mass_small], label="Small Inserts", color=:red)
	lines!([0, 100], [0, 100], label="Unity")
	
	axInt3.title = "IHU Mass vs Known Mass"
	axInt3.ylabel = "Calculated Mass (mg)"
	axInt3.xlabel = "Known Mass (mg)"
	axInt3.xticks = [0, 25, 50, 75, 100]
	axInt3.yticks = [0, 25, 50, 75, 100]

	xlims!(axInt3, 0, 100)
	ylims!(axInt3, 0, 100)
	
	fInt3[1, 2] = Legend(fInt3, axInt3, framevisible = false)
	fInt3
end

# ╔═╡ 84cb4868-b30e-496c-8a1f-753b28e73834
begin
	fInt3_small = Figure()
	axInt3_small = Axis(fInt3_small[1, 1])

	scatter!(df_i_canon3[!, :ground_truth_mass_small], df_i_canon3[!, :calculated_mass_small], label="Small Inserts", color=:red)
	lines!([0, 6], [0, 6], label="Unity")
	
	axInt3_small.title = "IHU Mass vs Known Mass (Small)"
	axInt3_small.ylabel = "Calculated Mass (mg)"
	axInt3_small.xlabel = "Known Mass (mg)"
	
	fInt3_small[1, 2] = Legend(fInt3_small, axInt3_small, framevisible = false)
	fInt3_small
end

# ╔═╡ f8c5964e-ac6b-4f78-bdc6-fe41e2e6365d
df_i_canon5 = df_i_canon[25:end, :];

# ╔═╡ 6e695941-f077-409c-9147-305c004028c3
md"""
#### High Noise
"""

# ╔═╡ 93b746a3-a4e6-4003-a01d-82ef2c86c10d
begin
	fInt5 = Figure()
	axInt5 = Axis(fInt5[1, 1])

	scatter!(df_i_canon5[!, :ground_truth_mass_large], df_i_canon5[!, :calculated_mass_large], label="Large Inserts")
	scatter!(df_i_canon5[!, :ground_truth_mass_medium], df_i_canon5[!, :calculated_mass_medium], label="Medium Inserts")
	scatter!(df_i_canon5[!, :ground_truth_mass_small], df_i_canon5[!, :calculated_mass_small], label="Small Inserts", color=:red)
	lines!([0, 100], [0, 100], label="Unity")
	
	axInt5.title = "IHU Mass vs Known Mass"
	axInt5.ylabel = "Calculated Mass (mg)"
	axInt5.xlabel = "Known Mass (mg)"
	axInt5.xticks = [0, 25, 50, 75, 100]
	axInt5.yticks = [0, 25, 50, 75, 100]

	xlims!(axInt5, 0, 100)
	ylims!(axInt5, 0, 100)
	
	fInt5[1, 2] = Legend(fInt5, axInt5, framevisible = false)
	fInt5
end

# ╔═╡ b76b5430-f9e7-4160-9046-72272393c810
begin
	fInt5_small = Figure()
	axInt5_small = Axis(fInt5_small[1, 1])

	scatter!(df_i_canon5[!, :ground_truth_mass_small], df_i_canon5[!, :calculated_mass_small], label="Small Inserts", color=:red)
	lines!([0, 8], [0, 8], label="Unity")
	
	axInt5_small.title = "IHU Mass vs Known Mass (Small)"
	axInt5_small.ylabel = "Calculated Mass (mg)"
	axInt5_small.xlabel = "Known Mass (mg)"
	
	fInt5_small[1, 2] = Legend(fInt5_small, axInt5_small, framevisible = false)
	fInt5_small
end

# ╔═╡ 8fc75b24-498a-469e-bf38-09e765ba723a
md"""
### RMSD measurements (large, medium, and small inserts)
"""

# ╔═╡ 130a8f9c-17bb-4aee-9e0a-4a394ec2546c
md"""
#### No Noise
"""

# ╔═╡ 971045b5-4ce1-499d-b86c-0616be49b18f
int1_rms_large_c = rmsd(df_i_canon0[!, :calculated_mass_large], df_i_canon0[!, :ground_truth_mass_large])

# ╔═╡ 3fbcc4ec-4a2c-476d-bea9-59eacd9d2993
int1_rms_medium_c = rmsd(df_i_canon0[!, :calculated_mass_medium], df_i_canon0[!, :ground_truth_mass_medium])

# ╔═╡ cb0de4f3-9918-472f-88aa-2b20297e3113
int1_rms_small_c = rmsd(df_i_canon0[!, :calculated_mass_small], df_i_canon0[!, :ground_truth_mass_small])

# ╔═╡ bcdba2c4-bbcb-4943-b90f-f2e6610accd7
md"""
#### Medium Noise
"""

# ╔═╡ f5ed82e9-4804-4adf-bf21-d647b2a495cd
int1_rms_large_c3 = rmsd(df_i_canon3[!, :calculated_mass_large], df_i_canon3[!, :ground_truth_mass_large])

# ╔═╡ 64e9f3ac-7ac3-4afd-aa41-9418b02696d1
int1_rms_medium_c3 = rmsd(df_i_canon3[!, :calculated_mass_medium], df_i_canon3[!, :ground_truth_mass_medium])

# ╔═╡ 5193ff4e-9ac4-4102-8f7c-e05ac7a192ee
int1_rms_small_c3 = rmsd(df_i_canon3[!, :calculated_mass_small], df_i_canon3[!, :ground_truth_mass_small])

# ╔═╡ 8cec0152-e592-48f4-9062-b9c18f28bf1b
md"""
#### High Noise
"""

# ╔═╡ 453cf000-9286-4b0e-ab6d-2799b0ddbbc0
int1_rms_large_c5 = rmsd(df_i_canon5[!, :calculated_mass_large], df_i_canon5[!, :ground_truth_mass_large])

# ╔═╡ 9bb69120-11e0-44a3-87a2-b50cc7d3456c
int1_rms_medium_c5 = rmsd(df_i_canon5[!, :calculated_mass_medium], df_i_canon5[!, :ground_truth_mass_medium])

# ╔═╡ e73629fa-73f3-4606-a00f-17f1520b7a3c
int1_rms_small_c5 = rmsd(df_i_canon5[!, :calculated_mass_small], df_i_canon5[!, :ground_truth_mass_small])

# ╔═╡ 87dd5b3b-204f-45c1-bf6b-e6f00ae0c3ec
md"""
## Agatston Mass Scoring
"""

# ╔═╡ 25dc45c9-3d20-4657-82e8-3a20dd2f2b9f
canon_path_agat = "/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/output/agatston";

# ╔═╡ fd5c8ed9-3f29-4af6-8e20-47a1d82bafbe
df_a_canon = CSV.read(string(canon_path_agat, "/full.csv"), DataFrame);

# ╔═╡ 87105910-8250-4332-b470-2650a3a425a7
md"""
#### No Noise
"""

# ╔═╡ 55ff7f68-0380-42a3-ba84-2a9fad909cd5
df_a_canon0 = df_a_canon[1:12, :];

# ╔═╡ d2cbd9fb-1f36-4096-86f2-dc9ebf3d8d85
begin
	fAgat = Figure()
	axAgat = Axis(fAgat[1, 1])

	scatter!(df_a_canon0[!, :ground_truth_mass_large], df_a_canon0[!, :calculated_mass_large], label="Mass: Large Inserts")
	scatter!(df_a_canon0[!, :ground_truth_mass_medium], df_a_canon0[!, :calculated_mass_medium], label="Mass: Medium Inserts")
	scatter!(df_a_canon0[!, :ground_truth_mass_small], df_a_canon0[!, :calculated_mass_small], label="Mass: Small Inserts", color=:red)
	lines!([0, 100], [0, 100], label="Unity Line")
	
	axAgat.title = "AS Mass vs Known Mass"
	axAgat.ylabel = "Calculated Mass (mg)"
	axAgat.xlabel = "Known Mass (mg)"
	axAgat.xticks = [0, 25, 50, 75, 100]
	axAgat.yticks = [0, 25, 50, 75, 100]

	xlims!(axAgat, 0, 100)
	ylims!(axAgat, 0, 100)
	
	fAgat[1, 2] = Legend(fAgat, axAgat, framevisible = false)

	fAgat
end

# ╔═╡ 82556fa3-84f2-4f3e-8345-a5d07465ebac
begin
	fAgat_small = Figure()
	axAgat_small = Axis(fAgat_small[1, 1])

	scatter!(df_a_canon0[!, :ground_truth_mass_small], df_a_canon0[!, :calculated_mass_small], label="Mass: Small Inserts", color=:red)
	lines!([0, 8], [0, 8], label="Unity Line")
	
	axAgat_small.title = "AS Mass vs Known Mass (Small)"
	axAgat_small.ylabel = "Calculated Mass (mg)"
	axAgat_small.xlabel = "Known Mass (mg)"
	
	fAgat_small[1, 2] = Legend(fAgat_small, axAgat_small, framevisible = false)

	fAgat_small
end

# ╔═╡ 62a683d2-68ba-40cf-be5f-089c61c64e04
md"""
#### Medium Noise
"""

# ╔═╡ ce7e5964-7ef5-49af-ac1b-d54954c45453
df_a_canon3 = df_a_canon[13:24, :];

# ╔═╡ d0071c6d-1ed7-4be1-8bca-4f6bd9459459
begin
	fAgat3 = Figure()
	axAgat3 = Axis(fAgat3[1, 1])

	scatter!(df_a_canon3[!, :ground_truth_mass_large], df_a_canon3[!, :calculated_mass_large], label="Mass: Large Inserts")
	scatter!(df_a_canon3[!, :ground_truth_mass_medium], df_a_canon3[!, :calculated_mass_medium], label="Mass: Medium Inserts")
	scatter!(df_a_canon3[!, :ground_truth_mass_small], df_a_canon3[!, :calculated_mass_small], label="Mass: Small Inserts", color=:red)
	lines!([0, 100], [0, 100], label="Unity Line")
	
	axAgat3.title = "AS Mass vs Known Mass"
	axAgat3.ylabel = "Calculated Mass (mg)"
	axAgat3.xlabel = "Known Mass (mg)"
	axAgat3.xticks = [0, 25, 50, 75, 100]
	axAgat3.yticks = [0, 25, 50, 75, 100]

	xlims!(axAgat3, 0, 100)
	ylims!(axAgat3, 0, 100)
	
	fAgat3[1, 2] = Legend(fAgat3, axAgat3, framevisible = false)

	fAgat3
end

# ╔═╡ d66a4a26-a24f-4e5a-989d-e86309ce2158
begin
	fAgat3_small = Figure()
	axAgat3_small = Axis(fAgat3_small[1, 1])

	scatter!(df_a_canon3[!, :ground_truth_mass_small], df_a_canon3[!, :calculated_mass_small], label="Mass: Small Inserts", color=:red)
	lines!([0, 10], [0, 10], label="Unity Line")
	
	axAgat3_small.title = "AS Mass vs Known Mass (Small)"
	axAgat3_small.ylabel = "Calculated Mass (mg)"
	axAgat3_small.xlabel = "Known Mass (mg)"
	
	fAgat3_small[1, 2] = Legend(fAgat3_small, axAgat3_small, framevisible = false)

	fAgat3_small
end

# ╔═╡ 570dda7e-5edc-4a5a-86bf-3af32855dbbb
md"""
#### High Noise
"""

# ╔═╡ b4f96b18-a85b-4772-b3a8-da1858655778
df_a_canon5 = df_a_canon[25:end, :];

# ╔═╡ 494d1b78-55e8-4b40-bb71-280d877cb9d7
begin
	fAgat5 = Figure()
	axAgat5 = Axis(fAgat5[1, 1])

	scatter!(df_a_canon5[!, :ground_truth_mass_large], df_a_canon5[!, :calculated_mass_large], label="Mass: Large Inserts")
	scatter!(df_a_canon5[!, :ground_truth_mass_medium], df_a_canon5[!, :calculated_mass_medium], label="Mass: Medium Inserts")
	scatter!(df_a_canon5[!, :ground_truth_mass_small], df_a_canon5[!, :calculated_mass_small], label="Mass: Small Inserts", color=:red)
	lines!([0, 100], [0, 100], label="Unity Line")
	
	axAgat5.title = "AS Mass vs Known Mass"
	axAgat5.ylabel = "Calculated Mass (mg)"
	axAgat5.xlabel = "Known Mass (mg)"
	axAgat5.xticks = [0, 25, 50, 75, 100]
	axAgat5.yticks = [0, 25, 50, 75, 100]

	xlims!(axAgat5, 0, 100)
	ylims!(axAgat5, 0, 100)
	
	fAgat5[1, 2] = Legend(fAgat5, axAgat5, framevisible = false)

	fAgat5
end

# ╔═╡ dc52c625-fc52-4073-8e18-e5e670fd1466
begin
	fAgat5_small = Figure()
	axAgat5_small = Axis(fAgat5_small[1, 1])

	scatter!(df_a_canon5[!, :ground_truth_mass_small], df_a_canon5[!, :calculated_mass_small], label="Mass: Small Inserts", color=:red)
	lines!([0, 12], [0, 12], label="Unity Line")
	
	axAgat5_small.title = "AS Mass vs Known Mass (Small)"
	axAgat5_small.ylabel = "Calculated Mass (mg)"
	axAgat5_small.xlabel = "Known Mass (mg)"
	
	fAgat5_small[1, 2] = Legend(fAgat5_small, axAgat5_small, framevisible = false)

	fAgat5_small
end

# ╔═╡ 1971e352-b0ea-43b5-81ef-d2b3962579aa
md"""
### RMSD measurements (large, medium, and small inserts)
"""

# ╔═╡ ec62d6e8-8bf5-4580-89f6-5148c3d921b6
md"""
#### No Noise
"""

# ╔═╡ 6be8cf2d-5e2a-4e78-b13e-5695e2929d28
agat_rms_large_c = rmsd(df_a_canon0[!, :calculated_mass_large], df_a_canon0[!, :ground_truth_mass_large])

# ╔═╡ 1401e167-732a-400f-9fff-ae2d8262bfff
agat_rms_medium_c = rmsd(df_a_canon0[!, :calculated_mass_medium], df_a_canon0[!, :ground_truth_mass_medium])

# ╔═╡ 757eddf5-ccfd-4b16-a55e-89bfc1807e3f
agat_rms_small_c = rmsd(df_a_canon0[!, :calculated_mass_small], df_a_canon0[!, :ground_truth_mass_small])

# ╔═╡ 5f87806b-a22a-4771-b302-05badfd30907
md"""
#### Medium Noise
"""

# ╔═╡ 983f2f6e-2e63-424b-aab0-7f1f629e13ef
agat_rms_large_c3 = rmsd(df_a_canon3[!, :calculated_mass_large], df_a_canon3[!, :ground_truth_mass_large])

# ╔═╡ 3fb10834-02b2-4110-80b1-2850a287c80f
agat_rms_medium_c3 = rmsd(df_a_canon3[!, :calculated_mass_medium], df_a_canon3[!, :ground_truth_mass_medium])

# ╔═╡ cc4111cc-a3aa-4b62-a3ec-5d7b82ff94d3
agat_rms_small_c3 = rmsd(df_a_canon3[!, :calculated_mass_small], df_a_canon3[!, :ground_truth_mass_small])

# ╔═╡ 901f2679-bebc-491c-b19c-26be3b6e7f68
md"""
#### High Noise
"""

# ╔═╡ 3cffbdec-a1b0-4791-a060-1af60a2316ec
agat_rms_large_c5 = rmsd(df_a_canon5[!, :calculated_mass_large], df_a_canon5[!, :ground_truth_mass_large])

# ╔═╡ 99221ba5-e14f-4908-9adc-8d88005c0a2b
agat_rms_medium_c5 = rmsd(df_a_canon5[!, :calculated_mass_medium], df_a_canon5[!, :ground_truth_mass_medium])

# ╔═╡ 22f86532-f7fa-4f4e-adac-52933abdb8c4
agat_rms_small_c5 = rmsd(df_a_canon5[!, :calculated_mass_small], df_a_canon5[!, :ground_truth_mass_small])

# ╔═╡ Cell order:
# ╠═c2496cdd-8ddb-4eeb-9747-f6b5397d4be6
# ╠═77c51949-937b-48fc-9621-5a071461f65d
# ╟─7bf86054-44de-456c-adc2-27617153faba
# ╠═41012935-cb85-4fc5-8f7c-7f4d27c2c7e4
# ╠═be63d9b6-bf33-4e05-a3fa-51d77e5288c5
# ╟─b4a759df-05d7-47fb-a52a-e9848070da15
# ╠═de94dfde-06a5-4991-9854-d286bd01c7bf
# ╟─687297b7-6f31-41fb-9549-6262c27b5cdb
# ╟─5eeb159c-57b0-470b-81b4-493c2112314b
# ╟─4977450f-882b-4913-b838-40fc69a66f4c
# ╠═0753c940-289c-4b91-8bfa-8544d553bed5
# ╟─f83b5d8f-7f03-4d1c-9ba0-8baf11e3d364
# ╟─84cb4868-b30e-496c-8a1f-753b28e73834
# ╠═f8c5964e-ac6b-4f78-bdc6-fe41e2e6365d
# ╟─6e695941-f077-409c-9147-305c004028c3
# ╟─93b746a3-a4e6-4003-a01d-82ef2c86c10d
# ╟─b76b5430-f9e7-4160-9046-72272393c810
# ╟─8fc75b24-498a-469e-bf38-09e765ba723a
# ╟─130a8f9c-17bb-4aee-9e0a-4a394ec2546c
# ╠═971045b5-4ce1-499d-b86c-0616be49b18f
# ╠═3fbcc4ec-4a2c-476d-bea9-59eacd9d2993
# ╠═cb0de4f3-9918-472f-88aa-2b20297e3113
# ╟─bcdba2c4-bbcb-4943-b90f-f2e6610accd7
# ╠═f5ed82e9-4804-4adf-bf21-d647b2a495cd
# ╠═64e9f3ac-7ac3-4afd-aa41-9418b02696d1
# ╠═5193ff4e-9ac4-4102-8f7c-e05ac7a192ee
# ╟─8cec0152-e592-48f4-9062-b9c18f28bf1b
# ╠═453cf000-9286-4b0e-ab6d-2799b0ddbbc0
# ╠═9bb69120-11e0-44a3-87a2-b50cc7d3456c
# ╠═e73629fa-73f3-4606-a00f-17f1520b7a3c
# ╟─87dd5b3b-204f-45c1-bf6b-e6f00ae0c3ec
# ╠═25dc45c9-3d20-4657-82e8-3a20dd2f2b9f
# ╠═fd5c8ed9-3f29-4af6-8e20-47a1d82bafbe
# ╟─87105910-8250-4332-b470-2650a3a425a7
# ╠═55ff7f68-0380-42a3-ba84-2a9fad909cd5
# ╟─d2cbd9fb-1f36-4096-86f2-dc9ebf3d8d85
# ╟─82556fa3-84f2-4f3e-8345-a5d07465ebac
# ╟─62a683d2-68ba-40cf-be5f-089c61c64e04
# ╠═ce7e5964-7ef5-49af-ac1b-d54954c45453
# ╟─d0071c6d-1ed7-4be1-8bca-4f6bd9459459
# ╟─d66a4a26-a24f-4e5a-989d-e86309ce2158
# ╟─570dda7e-5edc-4a5a-86bf-3af32855dbbb
# ╠═b4f96b18-a85b-4772-b3a8-da1858655778
# ╟─494d1b78-55e8-4b40-bb71-280d877cb9d7
# ╟─dc52c625-fc52-4073-8e18-e5e670fd1466
# ╟─1971e352-b0ea-43b5-81ef-d2b3962579aa
# ╟─ec62d6e8-8bf5-4580-89f6-5148c3d921b6
# ╠═6be8cf2d-5e2a-4e78-b13e-5695e2929d28
# ╠═1401e167-732a-400f-9fff-ae2d8262bfff
# ╠═757eddf5-ccfd-4b16-a55e-89bfc1807e3f
# ╟─5f87806b-a22a-4771-b302-05badfd30907
# ╠═983f2f6e-2e63-424b-aab0-7f1f629e13ef
# ╠═3fb10834-02b2-4110-80b1-2850a287c80f
# ╠═cc4111cc-a3aa-4b62-a3ec-5d7b82ff94d3
# ╟─901f2679-bebc-491c-b19c-26be3b6e7f68
# ╠═3cffbdec-a1b0-4791-a060-1af60a2316ec
# ╠═99221ba5-e14f-4908-9adc-8d88005c0a2b
# ╠═22f86532-f7fa-4f4e-adac-52933abdb8c4

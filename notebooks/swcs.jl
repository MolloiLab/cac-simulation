### A Pluto.jl notebook ###
# v0.19.2

using Markdown
using InteractiveUtils

# ╔═╡ 2a9b41e8-bab2-11ec-3294-ddf4409d19b6
# ╠═╡ show_logs = false
begin
	let
		using Pkg
		Pkg.activate(mktempdir())
		Pkg.Registry.update()
		Pkg.add("PlutoUI")
		Pkg.add("Statistics")
		Pkg.add("LinearAlgebra")
		Pkg.add("Images")
		Pkg.add("DSP")
		Pkg.add("CairoMakie")
		Pkg.add("Noise")
		Pkg.add("Distributions")
		Pkg.add(url="https://github.com/Dale-Black/CalciumScoring.jl")
	end
	
	using PlutoUI
	using Statistics
	using LinearAlgebra
	using Images
	using DSP
	using CairoMakie
	using Noise
	using Distributions
	using CalciumScoring
end

# ╔═╡ abc5db86-a240-46a4-8114-6aa7b7445317
TableOfContents()

# ╔═╡ b0d93309-7d9d-4d83-9477-2f534a05ddfe
begin
	abstract type CalciumScore end
	struct SpatiallyWeighted <: CalciumScore end
end

# ╔═╡ 9f4a38f7-9f7d-4bca-a7d3-5daefe798c14
function score_test(vol::AbstractMatrix, calibration, alg::SpatiallyWeighted)
    μ, σ = mean(calibration), std(calibration)
	d = Distributions.Normal(μ, σ)

	scaled_array = zeros(size(vol))
	for i in 1:size(vol, 1)
		for j in 1:size(vol, 2)
			scaled_array[i, j] = Distributions.cdf(d, vol[i, j])
		end
	end

	kern = [0 0.2 0; 0.2 0.2 0.2; 0 0.2 0]
	weighted_arr = conv(scaled_array, kern)[2:end-1, 2:end-1]
    
	return sum(weighted_arr)
end

# ╔═╡ 41c4242b-b025-497b-b5b8-f0e029ecbd8e
function score_test(vol::AbstractMatrix, μ, σ, alg::SpatiallyWeighted)
	d = Distributions.Normal(μ, σ)

	scaled_array = zeros(size(vol))
	for i in 1:size(vol, 1)
		for j in 1:size(vol, 2)
			scaled_array[i, j] = Distributions.cdf(d, vol[i, j])
		end
	end

	kern = [0 0.2 0; 0.2 0.2 0.2; 0 0.2 0]
	weighted_arr = conv(scaled_array, kern)[2:end-1, 2:end-1]
    
	return sum(weighted_arr)
end

# ╔═╡ 3cfc6aca-b1f0-41cf-9b7f-f1c84363ae34
function score_test(vol::AbstractArray, calibration, alg::SpatiallyWeighted)
    μ, σ = mean(calibration), std(calibration)
	d = Distributions.Normal(μ, σ)

	scaled_array = zeros(size(vol))
	for i in 1:size(vol, 1)
		for j in 1:size(vol, 2)
			for z in 1:size(vol, 3)
				scaled_array[i, j, z] = Distributions.cdf(d, vol[i, j, z])
			end
		end
	end

	weighted_arr = zeros(size(vol))
	for z in 1:size(scaled_array, 3)
		kern = [0 0.2 0; 0.2 0.2 0.2; 0 0.2 0]
		weighted_arr[:, :, z] = conv(scaled_array[:, :, z], kern)[2:end-1, 2:end-1]
	end
    
	return sum(weighted_arr)
end

# ╔═╡ 0e4a1d69-a342-4b0a-91b0-3667db8267b4
function score_test(vol::AbstractArray, μ, σ, alg::SpatiallyWeighted)
	d = Distributions.Normal(μ, σ)

	scaled_array = zeros(size(vol))
	for i in 1:size(vol, 1)
		for j in 1:size(vol, 2)
			for z in 1:size(vol, 3)
				scaled_array[i, j, z] = Distributions.cdf(d, vol[i, j, z])
			end
		end
	end

	weighted_arr = zeros(size(vol))
	for z in 1:size(scaled_array, 3)
		kern = [0 0.2 0; 0.2 0.2 0.2; 0 0.2 0]
		weighted_arr[:, :, z] = conv(scaled_array[:, :, z], kern)[2:end-1, 2:end-1]
	end
    
	return sum(weighted_arr)
end

# ╔═╡ 9929b3c2-6c07-4da0-9eb4-26d6c407d7e2
md"""
## Step 1
"""

# ╔═╡ 1b5d567f-4337-423b-8674-b4485c336ce8
vol = [
	0 46 58 123 133 104 65 7
	19 25 20 125 163 83 -22 -134
	127 99 65 104 139 57 -51 -102
	123 88 112 140 104 57 46 34
	122 31 71 121 93 83 101 72
	97 24 59 99 68 27 22 52
	0 22 29 72 80 54 28 72
	0 17 -20 15 42 61 80 83
	0 66 14 -3 11 54 110 148
	0 121 127 102 103 93 118 167
]

# ╔═╡ 1803d49e-287d-4779-ad34-157b6e977d0f
vol3D = cat(vol, vol, dims=3);

# ╔═╡ df800b73-6567-43c5-b7bc-df18234041dd
begin
	calibration = zeros(5, 5, 100)
	for i in 1:100
		calibration[:, :, i] = mult_gauss(ones(5, 5) * 170)
	end
end

# ╔═╡ 4baf3584-ddfe-4dd0-b63d-d9a43121b3cd
begin
	μ, σ = mean(calibration), std(calibration) + 30
	d = Distributions.Normal(μ, σ)
end

# ╔═╡ 6ddd523e-f43f-4410-bb50-d4b41989208c
scatter(d)

# ╔═╡ a076ec71-5a46-439d-938f-d7dcf4679248
begin
	# low, high = quantile.(d, [0.00001, 0.99999])
	xs = []
	ys = []
	for x in -100:300
		push!(xs, x)
		push!(ys, cdf(d, x))
	end
	xs = round.(xs, digits = 10)
	ys = round.(ys, digits = 10)
end

# ╔═╡ 3ef92c9c-59b7-42a4-abe9-38cbc8209358
scatter(xs, ys)

# ╔═╡ 85adf61f-f7a4-4576-921f-50f37631f134
begin
	scaled_array = zeros(size(vol))
	for i in 1:size(vol, 1)
		for j in 1:size(vol, 2)
			scaled_array[i, j] = Distributions.cdf(d, vol[i, j])
		end
	end
	round.(scaled_array, digits=2)
end

# ╔═╡ 91940958-2d95-472a-a7ec-7cfe88fe12cb
md"""
## Step 2
"""

# ╔═╡ c8b0da1f-5394-42d1-80e0-87b17417cb18
begin
	kern = [0 0.2 0; 0.2 0.2 0.2; 0 0.2 0]
	weighted_arr = DSP.conv(scaled_array, kern)[2:end-1, 2:end-1]
	round.(weighted_arr, digits=2)
end

# ╔═╡ c12d83f6-cfe0-4210-bb00-0ec8de046310
md"""
## Step 3
"""

# ╔═╡ 85179b38-a120-4d1f-8ce5-8cae6331e17a
sum(weighted_arr)

# ╔═╡ 41236dd9-fdcb-4c00-b89b-3ed742fe7f01
md"""
# Testing
"""

# ╔═╡ 3c5436d3-a925-4b63-ab0b-1adeac1146aa
md"""
### Agatston
"""

# ╔═╡ a5b06b4b-db94-42c5-a976-fd01f4da8661
struct Agatston <: CalciumScore end

# ╔═╡ 07ecf4e4-7e86-4ac6-be16-fe61d9ebed6b
function score_agat(vol, spacing, alg::Agatston; threshold=130, min_size_mm=1)
    area_mm = spacing[1] * spacing[2]
    min_size_pixels = Int(round(min_size_mm / area_mm))
    num_slices = size(vol, 3)
    score = 0
    for z in range(1, num_slices)
        slice = vol[:, :, z]
        thresholded_slice = slice .> threshold
        max_intensity = maximum(slice)
        if max_intensity < threshold
            continue
        end
        comp_connect = trues(min_size_pixels + 1, min_size_pixels + 1)
        lesion_map = label_components(thresholded_slice, comp_connect)
        num_non_zero = 0
        number_islands = 0
        slice_score = 0
        num_labels = length(unique(lesion_map))
        for label_idx in 0:num_labels
            if label_idx == 0
                continue
            end

            idxs = findall(x -> x == label_idx, lesion_map)
            num_label_idxs = length(idxs)
            if num_label_idxs < min_size_pixels
                continue
            end

            intensities = slice[idxs]
            max_intensity = maximum(intensities)
            weight = floor(max_intensity / 100)
            if weight > 4
                weight = 4.0
            end
            num_non_zero_in_island = num_label_idxs
            slice_score += num_non_zero_in_island * area_mm * weight
            num_non_zero += num_non_zero_in_island
            number_islands += 1
        end
        score += slice_score
    end
    return score
end

# ╔═╡ 9fff444c-d222-4fcc-8438-8c3217c336c4
function score_agat(vol, spacing, mass_cal_factor, alg::Agatston; threshold=130, min_size_mm=1)
    area_mm = spacing[1] * spacing[2]
    slice_thickness = spacing[3]
    min_size_pixels = Int(round(min_size_mm / area_mm))
    mass_score = 0
    score = 0
    for z in 1:size(vol, 3)
        slice = vol[:, :, z]
        thresholded_slice = slice .> threshold
        max_intensity = maximum(slice)
        if max_intensity < threshold
            continue
        end
        comp_connect = trues(min_size_pixels + 1, min_size_pixels + 1)
        lesion_map = label_components(thresholded_slice, comp_connect)
        num_non_zero = 0
        number_islands = 0
        mass_slice_score = 0
        slice_score = 0
        num_labels = length(unique(lesion_map))
        for label_idx in 1:(num_labels - 1)
            idxs = findall(x -> x == label_idx, lesion_map)
            num_label_idxs = length(idxs)
            if num_label_idxs < min_size_pixels
                continue
            end

            intensities = slice[idxs]
            max_intensity = maximum(intensities)
            weight = floor(max_intensity / 100)
            if weight > 4
                weight = 4.0
            end
            num_non_zero_in_island = num_label_idxs
            slice_score += num_non_zero_in_island * area_mm * weight
            num_non_zero += num_non_zero_in_island
            number_islands += 1
            mass_slice_score += mean(intensities)
        end
        plaque_vol = length(findall(x -> x != 0, lesion_map))
        plaque_vol = plaque_vol * area_mm * slice_thickness
        mass_score += mass_slice_score * plaque_vol * mass_cal_factor
        score += slice_score
    end
    return score, mass_score
end

# ╔═╡ 51a5a3b7-bab6-441d-92ca-4cec857e83a0
md"""
## 2D
"""

# ╔═╡ 94023c0d-7da5-4e6b-af0e-2730c5e3f84a
alg = SpatiallyWeighted()

# ╔═╡ 3a4d9b2c-a786-4c59-9f5c-05dfd38379ca
score_test(vol, μ, σ, alg)

# ╔═╡ 29e107be-4018-4b72-9489-84c8ec0ff685
alg2 = Agatston()

# ╔═╡ 0465cdaa-15eb-4a30-b218-78c888601f04
spacing = [1, 1, 3]

# ╔═╡ 2148710c-5ee2-4648-86c0-bf337af72b99
score(vol3D, spacing, alg2)

# ╔═╡ a3a9e388-fb66-4058-ba3b-3d6fe8d4dadd
md"""
## 3D
"""

# ╔═╡ 6f61f901-94b3-4475-bff6-be338c12ff90
# vol3D = cat(vol, vol, dims=3)

# ╔═╡ 71f3c1d5-82db-400d-a5d8-aee347521f5c
# calibration_ = calibration[:, :, 1]

# ╔═╡ b52f41a1-3f25-41e6-bbbb-2041dec891b9
# score(vol3D, calibration_, alg)

# ╔═╡ 49bbbd77-2bb1-489a-aaf2-8e43c86c680f
# pixel_size = [0.488, 0.488, 3.0]

# ╔═╡ 6ba1585e-3c91-475e-a301-132626b97c5b
# agat_l_hd = score(vol3D, pixel_size, alg)

# ╔═╡ Cell order:
# ╠═2a9b41e8-bab2-11ec-3294-ddf4409d19b6
# ╠═abc5db86-a240-46a4-8114-6aa7b7445317
# ╠═b0d93309-7d9d-4d83-9477-2f534a05ddfe
# ╠═9f4a38f7-9f7d-4bca-a7d3-5daefe798c14
# ╠═41c4242b-b025-497b-b5b8-f0e029ecbd8e
# ╠═3cfc6aca-b1f0-41cf-9b7f-f1c84363ae34
# ╠═0e4a1d69-a342-4b0a-91b0-3667db8267b4
# ╟─9929b3c2-6c07-4da0-9eb4-26d6c407d7e2
# ╠═1b5d567f-4337-423b-8674-b4485c336ce8
# ╠═1803d49e-287d-4779-ad34-157b6e977d0f
# ╠═df800b73-6567-43c5-b7bc-df18234041dd
# ╠═4baf3584-ddfe-4dd0-b63d-d9a43121b3cd
# ╠═6ddd523e-f43f-4410-bb50-d4b41989208c
# ╠═a076ec71-5a46-439d-938f-d7dcf4679248
# ╠═3ef92c9c-59b7-42a4-abe9-38cbc8209358
# ╠═85adf61f-f7a4-4576-921f-50f37631f134
# ╟─91940958-2d95-472a-a7ec-7cfe88fe12cb
# ╠═c8b0da1f-5394-42d1-80e0-87b17417cb18
# ╟─c12d83f6-cfe0-4210-bb00-0ec8de046310
# ╠═85179b38-a120-4d1f-8ce5-8cae6331e17a
# ╟─41236dd9-fdcb-4c00-b89b-3ed742fe7f01
# ╟─3c5436d3-a925-4b63-ab0b-1adeac1146aa
# ╠═a5b06b4b-db94-42c5-a976-fd01f4da8661
# ╠═07ecf4e4-7e86-4ac6-be16-fe61d9ebed6b
# ╠═9fff444c-d222-4fcc-8438-8c3217c336c4
# ╟─51a5a3b7-bab6-441d-92ca-4cec857e83a0
# ╠═94023c0d-7da5-4e6b-af0e-2730c5e3f84a
# ╠═3a4d9b2c-a786-4c59-9f5c-05dfd38379ca
# ╠═29e107be-4018-4b72-9489-84c8ec0ff685
# ╠═0465cdaa-15eb-4a30-b218-78c888601f04
# ╠═2148710c-5ee2-4648-86c0-bf337af72b99
# ╟─a3a9e388-fb66-4058-ba3b-3d6fe8d4dadd
# ╠═6f61f901-94b3-4475-bff6-be338c12ff90
# ╠═71f3c1d5-82db-400d-a5d8-aee347521f5c
# ╠═b52f41a1-3f25-41e6-bbbb-2041dec891b9
# ╠═49bbbd77-2bb1-489a-aaf2-8e43c86c680f
# ╠═6ba1585e-3c91-475e-a301-132626b97c5b

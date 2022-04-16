### A Pluto.jl notebook ###
# v0.18.1

using Markdown
using InteractiveUtils

# ╔═╡ 2a9b41e8-bab2-11ec-3294-ddf4409d19b6
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
	end
	
	using PlutoUI
	using Statistics
	using LinearAlgebra
	using Images
	using DSP
	using CairoMakie
	using Noise
	using Distributions
end

# ╔═╡ abc5db86-a240-46a4-8114-6aa7b7445317
TableOfContents()

# ╔═╡ c538773a-ed7d-42ee-bd2c-f2e792d369b7
md"""
## Create arrays
"""

# ╔═╡ 47d01a99-7aab-4810-90d7-160dc7081519
array_calcium = [
	10 30 135 145 5 150 100 80 9
	9 40 135 130 5 150 190 80 9
	10 100 100 130 5 9 0 80 18
	4 40 135 -10 5 150 100 0 9
	10 0 189 130 5 150 40 80 22
]

# ╔═╡ 1e6d27ee-943c-464f-912a-207edf61a3df
CairoMakie.heatmap(transpose(array_calcium))

# ╔═╡ 680787bc-797e-4454-8450-2cb0495efe0c
array_calibration = ones(10, 10) .* 200;

# ╔═╡ 60041821-dcbc-4cb6-b000-86379d01a815
noisy_array = mult_gauss(salt_pepper(array_calibration))

# ╔═╡ 5313555e-239e-4534-9fec-46f84e3973b1
CairoMakie.heatmap(noisy_array)

# ╔═╡ cb10a1e7-e292-466f-bd65-39a117e9a303
md"""
## Step 1
Create weighting function
"""

# ╔═╡ 0e6887f1-b68d-4584-89fe-a275caf85a55
vec_noisy = vec(noisy_array)

# ╔═╡ d9bb2391-907f-4673-8dba-70fc68b725b0
begin
	μ200 = mean(vec_noisy)
	σ200 = std(vec_noisy)
	d = Normal(μ200, σ200)
end

# ╔═╡ 1f5d5c9b-67dc-4ad1-bb9d-75b7d6d4b7fc
CairoMakie.scatter(d)

# ╔═╡ a024c1d0-9c61-49de-841b-ca2104d79bee
begin
	low, high = quantile.(d, [0.00001, 0.99999])
	cdf_array = cdf(d, range(low, high))
end

# ╔═╡ 0ef24b04-fc1f-4e88-8c8d-15f414542368
CairoMakie.scatter(cdf_array)

# ╔═╡ e7c7fe28-0abc-4b2c-8ad9-c5856928a5bd
md"""
Looks like `cdf(d, array[i])` is the weighting function we want. We can now apply that to the noisy array on all pixels
"""

# ╔═╡ 23b743a5-6438-4b65-a1a6-42935dd8c22b
begin
	scaled_array = zeros(size(noisy_array))
	for i in 1:size(noisy_array, 1)
		for j in 1:size(noisy_array, 2)
			scaled_array[i, j] =  cdf(d, noisy_array[i, j])
		end
	end
end

# ╔═╡ ceee29a6-9df1-4abd-ac86-d15cf71af780
scaled_array

# ╔═╡ 64919a4d-2833-4644-853a-65084e518f34
CairoMakie.heatmap(scaled_array)

# ╔═╡ 00777d97-194f-43da-b4f3-d15d432b0981
md"""
## Step 2
"""

# ╔═╡ 7d087894-e95e-47d9-84d3-013952240563
kern = [
		0   0.2 0
		0.2 0.2 0.2
		0   0.2 0
	]

# ╔═╡ 99c7f282-fc93-425b-aca6-472bdb352b8a
weighted_arr = conv(scaled_array, kern)[2:end-1, 2:end-1]

# ╔═╡ 4125e4e8-6eed-442d-9cec-ad072b2bd9c7
CairoMakie.heatmap(weighted_arr)

# ╔═╡ 99b81fa5-3a36-4a7a-8078-0154d0032443
md"""
## Step 3
"""

# ╔═╡ 7ba3a955-b6e1-498f-9b59-1698968347bb
swcs = sum(weighted_arr)

# ╔═╡ 41e00d27-1edf-4387-a777-7287578e8dc4
md"""
## Full function
"""

# ╔═╡ b0d93309-7d9d-4d83-9477-2f534a05ddfe
begin
	abstract type CalciumScore end
	struct SpatiallyWeighted <: CalciumScore end
end

# ╔═╡ 9f4a38f7-9f7d-4bca-a7d3-5daefe798c14
function score(vol, calibration, alg::SpatiallyWeighted)
    μ, σ= mean(calibration), std(calibration)
	d = Distributions.Normal(μ, σ)

	scaled_array = zeros(size(vol))
	for i in 1:size(vol, 1)
		for j in 1:size(vol, 2)
            if length(size(vol)) == 3
                for z in 1:size(vol, 3)
                    scaled_array[i, j, z] = Distributions.cdf(d, vol[i, j, z])
                end
            else
			    scaled_array[i, j] = Distributions.cdf(d, vol[i, j])
            end
		end
	end

    local weighted_arr
    if length(size(vol)) == 3
		kern1 = [0 0.2 0; 0.2 0.2 0.2; 0 0.2 0]
		kern2 = [0 0.2 0; 0.2 0.2 0.2; 0 0.2 0]
		kern3 = [0 0.2 0; 0.2 0.2 0.2; 0 0.2 0]
		kern = cat(kern1, kern2, kern3, dims=3)
        weighted_arr = conv(scaled_array, kern)[2:end-1, 2:end-1, 2:end-1]
    else
		kern = [0 0.2 0; 0.2 0.2 0.2; 0 0.2 0]
	    weighted_arr = conv(scaled_array, kern)[2:end-1, 2:end-1]
    end
    
	return sum(weighted_arr)
end

# ╔═╡ 94023c0d-7da5-4e6b-af0e-2730c5e3f84a
alg = SpatiallyWeighted()

# ╔═╡ 32bb7e29-0aa3-441e-94fb-d47a5bc46320
score(array_calcium, noisy_array, alg)

# ╔═╡ 68452530-776b-43eb-9ac9-fe68977ab1ed
md"""
### Test 3D
"""

# ╔═╡ 2f848513-54ac-4153-9f07-3eef0c7f6652
array_calcium3D = cat(array_calcium, array_calcium, dims=3);

# ╔═╡ e66e2f21-6fa6-4b47-bff3-e9c9a5f41320
noisy_array3D = cat(noisy_array, noisy_array, dims=3);

# ╔═╡ 3ae79fb3-18d6-4b18-91a7-7dbba80bc305
score(array_calcium3D, noisy_array3D, alg)

# ╔═╡ 41236dd9-fdcb-4c00-b89b-3ed742fe7f01
md"""
# Testing
"""

# ╔═╡ 51a5a3b7-bab6-441d-92ca-4cec857e83a0
md"""
## 2D
"""

# ╔═╡ 1b5d567f-4337-423b-8674-b4485c336ce8
vol = [
	10 30 135 145 5 150 100 80 9
	9 40 135 130 5 150 190 80 9
	10 100 100 130 5 9 0 80 18
	4 40 135 -10 5 150 100 0 9
	10 0 189 130 5 150 40 80 22
]

# ╔═╡ df800b73-6567-43c5-b7bc-df18234041dd
calibration = [
	 196.569  207.709  177.004  207.147  178.809
	 207.359  211.824    0.0    168.913  202.965
	 160.564    0.0    208.593  214.862  220.32
	 158.685  188.71   226.872  202.076  226.995
	 186.996  212.137  185.208  228.175  234.261
]

# ╔═╡ 3a4d9b2c-a786-4c59-9f5c-05dfd38379ca
score(vol, calibration, alg)

# ╔═╡ a3a9e388-fb66-4058-ba3b-3d6fe8d4dadd
md"""
## 3D
"""

# ╔═╡ 6f61f901-94b3-4475-bff6-be338c12ff90
vol3D = cat(vol, vol, dims=3)

# ╔═╡ 71f3c1d5-82db-400d-a5d8-aee347521f5c
calibration3D = cat(calibration, calibration, dims=3)

# ╔═╡ b52f41a1-3f25-41e6-bbbb-2041dec891b9
score(vol3D, calibration3D, alg)

# ╔═╡ Cell order:
# ╠═2a9b41e8-bab2-11ec-3294-ddf4409d19b6
# ╠═abc5db86-a240-46a4-8114-6aa7b7445317
# ╟─c538773a-ed7d-42ee-bd2c-f2e792d369b7
# ╠═47d01a99-7aab-4810-90d7-160dc7081519
# ╠═1e6d27ee-943c-464f-912a-207edf61a3df
# ╠═680787bc-797e-4454-8450-2cb0495efe0c
# ╠═60041821-dcbc-4cb6-b000-86379d01a815
# ╠═5313555e-239e-4534-9fec-46f84e3973b1
# ╟─cb10a1e7-e292-466f-bd65-39a117e9a303
# ╠═0e6887f1-b68d-4584-89fe-a275caf85a55
# ╠═d9bb2391-907f-4673-8dba-70fc68b725b0
# ╠═1f5d5c9b-67dc-4ad1-bb9d-75b7d6d4b7fc
# ╠═a024c1d0-9c61-49de-841b-ca2104d79bee
# ╠═0ef24b04-fc1f-4e88-8c8d-15f414542368
# ╟─e7c7fe28-0abc-4b2c-8ad9-c5856928a5bd
# ╠═23b743a5-6438-4b65-a1a6-42935dd8c22b
# ╠═ceee29a6-9df1-4abd-ac86-d15cf71af780
# ╠═64919a4d-2833-4644-853a-65084e518f34
# ╟─00777d97-194f-43da-b4f3-d15d432b0981
# ╠═7d087894-e95e-47d9-84d3-013952240563
# ╠═99c7f282-fc93-425b-aca6-472bdb352b8a
# ╠═4125e4e8-6eed-442d-9cec-ad072b2bd9c7
# ╟─99b81fa5-3a36-4a7a-8078-0154d0032443
# ╠═7ba3a955-b6e1-498f-9b59-1698968347bb
# ╟─41e00d27-1edf-4387-a777-7287578e8dc4
# ╠═b0d93309-7d9d-4d83-9477-2f534a05ddfe
# ╠═9f4a38f7-9f7d-4bca-a7d3-5daefe798c14
# ╠═94023c0d-7da5-4e6b-af0e-2730c5e3f84a
# ╠═32bb7e29-0aa3-441e-94fb-d47a5bc46320
# ╟─68452530-776b-43eb-9ac9-fe68977ab1ed
# ╠═2f848513-54ac-4153-9f07-3eef0c7f6652
# ╠═e66e2f21-6fa6-4b47-bff3-e9c9a5f41320
# ╠═3ae79fb3-18d6-4b18-91a7-7dbba80bc305
# ╟─41236dd9-fdcb-4c00-b89b-3ed742fe7f01
# ╟─51a5a3b7-bab6-441d-92ca-4cec857e83a0
# ╠═1b5d567f-4337-423b-8674-b4485c336ce8
# ╠═df800b73-6567-43c5-b7bc-df18234041dd
# ╠═3a4d9b2c-a786-4c59-9f5c-05dfd38379ca
# ╟─a3a9e388-fb66-4058-ba3b-3d6fe8d4dadd
# ╠═6f61f901-94b3-4475-bff6-be338c12ff90
# ╠═71f3c1d5-82db-400d-a5d8-aee347521f5c
# ╠═b52f41a1-3f25-41e6-bbbb-2041dec891b9

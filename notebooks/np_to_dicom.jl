### A Pluto.jl notebook ###
# v0.19.16

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 752b1e22-dc1d-4dcb-acf6-2ac6defe0781
# ╠═╡ show_logs = false
begin
    let
        using Pkg
        Pkg.activate(".")
    end

    using PlutoUI
    using CairoMakie
    using NPZ
    using DICOM
    using DICOMUtils
end

# ╔═╡ 63bcb593-6481-49a5-842c-e24819424910
TableOfContents()

# ╔═╡ 47af611d-5880-4c28-a5d1-e103b56d1d1e
densities = ["low", "normal"]

# ╔═╡ c3f731a7-dd5b-44b1-8172-d37496fd1514
sizes = ["small", "medium", "large"]

# ╔═╡ bafb59ed-c807-4078-8dc5-56e6151fbadf
energies_motion = ["80-motion-3.npy", "100-motion.npy", "120-motion.npy", "135-motion.npy"]

# ╔═╡ 22c41290-b952-4e45-87a8-b020fd6466fd
slices = [1, 2, 3, 4, 5, 6]

# ╔═╡ 09986775-f1ee-434a-a085-97de5d14618b
root = joinpath(dirname(pwd()), "images_new")

# ╔═╡ bc8c5e58-0094-476f-9f24-aee0f9a7ff92
for SIZE in sizes
	for DENSITY in densities
        for ENERGY_MOTION in energies_motion
			ENERGY = split(ENERGY_MOTION, "-")[1]

			
			np_path = joinpath(root, SIZE, DENSITY, ENERGY_MOTION)
			motion_path = joinpath(root, SIZE, DENSITY, ENERGY * "-motion")
			np = npzread(np_path)
			np = permutedims(np, (2, 1, 3))

			dcm_path = "/Users/daleblack/Google Drive/Datasets/Canon_Aquilion_One_Vision/Large_rep1/96E1EB4F"
			dcm = dcm_parse(dcm_path)

			@info motion_path
			if !isdir(motion_path)
				mkpath(motion_path)
			end
			# rm(motion_path; force=true)
			for z in slices
				dcm[tag"Pixel Data"] = Int16.(round.(np[:, :, z]))
				dcm[tag"Instance Number"] = z
				dcm[tag"Rows"] = size(np, 1)
				dcm[tag"Columns"] = size(np, 2)
				output_path = joinpath(motion_path, string(z) * ".dcm")
				# if !isdir(dirname(output_path))
				# 	mkpath(dirname(output_path))
				# end
				# @info output_path
            	dcm_write(output_path, dcm)
			end
		end
	end
end

# ╔═╡ ec512e95-7be0-4f27-98fd-635a5e4f759b
md"""
## Check DICOM image(s)
"""

# ╔═╡ 9ccad42f-6e1a-4a15-8230-426c47cf976d
path = "/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/images_new/large/normal/80"

# ╔═╡ caadc496-9fa7-49f8-a6dd-8df22e99e00b
dcmdir_combined = dcmdir_parse(path);

# ╔═╡ 9d3621de-9c59-4e64-a2c3-82669a8bc749
vol_combined = load_dcm_array(dcmdir_combined);

# ╔═╡ 9b16416f-e894-41b9-8494-b19fad598f08
@bind c PlutoUI.Slider(1:size(vol_combined, 3); default=1, show_value=true)

# ╔═╡ b4ada9d2-b955-4b53-92ce-805231fb781f
heatmap(transpose(vol_combined[:, :, c]); colormap=:grays)

# Cell order:
# ╠═4427b686-e049-11ec-3271-ab17a8f19bfc
# ╠═20bf9bda-4a51-4f64-a299-86d6cb6ea22f
# ╠═97d46d78-70ff-42b0-8e99-51721390e171
# ╠═726188c6-0938-4647-8aa1-fca66c983ee3
# ╠═ac21a7ea-a2c5-401e-affd-43152cdc85aa
# ╠═f02d5b27-9184-45e8-b23f-df6a86867fe9
# ╠═8d9777f7-22c6-4c61-b29d-85e3300b5c45
# ╟─258103d8-16ff-46ec-9a93-0fb526f17d3f
# ╠═b4a9e9b5-fcd2-4d62-a320-b8651674064a
# ╠═7eff2120-6f38-43c0-bf57-6ccdc49a9f7a
# ╠═9b7bb48e-6e79-4c06-824e-88f45bde68f0
# ╟─0b0590d0-8e45-48d4-bc70-9c7c33b3fed1
# ╠═a07b20bc-3891-4d43-b63b-340c50d3183f

# ╔═╡ 8e5e4d70-7a18-4016-822d-5aba6a4e17cf
# function to_array(dcm_data::Vector{DICOM.DICOMData})
#     return array = cat(
#         [dcm_data[i][(0x7fe0, 0x0010)] for i in 1:length(dcm_data)]...; dims=3
#     )
# end

# ╔═╡ fcd5f547-de61-4373-9920-e93bd4fce0e9
# begin
# 	sizes = ["large", "medium", "small"]
# 	densities = ["low", "normal"]
# 	energies = ["80", "100", "120", "135"]
# 	energies_motion = ["80-motion-3.npy", "100-motion.npy", "120-motion.npy", "135-motion.npy"]
# end

# ╔═╡ d9aa08c3-f5cd-4792-8c67-7d7be206a0cf
# np_path = joinpath(dirname(pwd()), "images_new", sizes[1], densities[1], energies_motion[1])

# ╔═╡ 2e81595a-6547-429c-9143-709b7c4b6a2d
# dcm_path = joinpath(dirname(pwd()), "images_new", sizes[1], densities[1], energies[1])

# ╔═╡ ce1c1613-6070-47fa-9950-50201d8245bf
# begin
# 	np = npzread(np_path);
# 	np_array = permutedims(np, (2, 1, 3))
# end;

# ╔═╡ c8813c40-cbbb-468a-8c86-b5fb928ffd7e
# begin
# 	dcm = dcmdir_parse(dcm_path)
# 	dcm_array = to_array(dcm)
# end;

# ╔═╡ 23c795db-a65a-4229-a594-26c94f69a167


# ╔═╡ 397a186f-76f9-4063-b12d-aca8bc0f2228
# @bind a PlutoUI.Slider(axes(np_array, 3), default=4, show_value=true)

# ╔═╡ 02fa1ae2-3e06-4c6e-bdb4-ad8bc95eccea
# heatmap(np_array[:, :, a]; colormap=:grays)

# ╔═╡ 152e5f0f-f020-40a2-914a-8400d1d49424
# @bind b PlutoUI.Slider(axes(np_array, 3), default=4, show_value=true)

# ╔═╡ ff7eebf3-3602-4a51-80a0-aebabe5b5721
# heatmap(dcm_array[:, :, b]; colormap=:grays)

# ╔═╡ Cell order:
# ╠═752b1e22-dc1d-4dcb-acf6-2ac6defe0781
# ╠═63bcb593-6481-49a5-842c-e24819424910
# ╠═47af611d-5880-4c28-a5d1-e103b56d1d1e
# ╠═c3f731a7-dd5b-44b1-8172-d37496fd1514
# ╠═bafb59ed-c807-4078-8dc5-56e6151fbadf
# ╠═22c41290-b952-4e45-87a8-b020fd6466fd
# ╠═09986775-f1ee-434a-a085-97de5d14618b
# ╠═bc8c5e58-0094-476f-9f24-aee0f9a7ff92
# ╟─ec512e95-7be0-4f27-98fd-635a5e4f759b
# ╠═9ccad42f-6e1a-4a15-8230-426c47cf976d
# ╠═caadc496-9fa7-49f8-a6dd-8df22e99e00b
# ╠═9d3621de-9c59-4e64-a2c3-82669a8bc749
# ╠═9b16416f-e894-41b9-8494-b19fad598f08
# ╠═b4ada9d2-b955-4b53-92ce-805231fb781f
# ╠═8e5e4d70-7a18-4016-822d-5aba6a4e17cf
# ╠═fcd5f547-de61-4373-9920-e93bd4fce0e9
# ╠═d9aa08c3-f5cd-4792-8c67-7d7be206a0cf
# ╠═2e81595a-6547-429c-9143-709b7c4b6a2d
# ╠═ce1c1613-6070-47fa-9950-50201d8245bf
# ╠═c8813c40-cbbb-468a-8c86-b5fb928ffd7e
# ╠═23c795db-a65a-4229-a594-26c94f69a167
# ╠═397a186f-76f9-4063-b12d-aca8bc0f2228
# ╠═02fa1ae2-3e06-4c6e-bdb4-ad8bc95eccea
# ╠═152e5f0f-f020-40a2-914a-8400d1d49424
# ╠═ff7eebf3-3602-4a51-80a0-aebabe5b5721

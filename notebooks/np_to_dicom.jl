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
energies_motion = ["80-motion-3.npy", "100-motion-3.npy", "120-motion-3.npy", "135-motion-3.npy"]

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
path = "/Users/daleblack/Google Drive/dev/MolloiLab/cac-simulation/images_new/large/normal/80-motion"

# ╔═╡ caadc496-9fa7-49f8-a6dd-8df22e99e00b
dcmdir_combined = dcmdir_parse(path);

# ╔═╡ 9d3621de-9c59-4e64-a2c3-82669a8bc749
vol_combined = load_dcm_array(dcmdir_combined);

# ╔═╡ 9b16416f-e894-41b9-8494-b19fad598f08
@bind c PlutoUI.Slider(1:size(vol_combined, 3); default=1, show_value=true)

# ╔═╡ b4ada9d2-b955-4b53-92ce-805231fb781f
heatmap(transpose(vol_combined[:, :, c]); colormap=:grays)

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
# ╟─9b16416f-e894-41b9-8494-b19fad598f08
# ╠═b4ada9d2-b955-4b53-92ce-805231fb781f

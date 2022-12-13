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

# ╔═╡ 5571cce0-224e-4584-881f-b136dad2d05b
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

# ╔═╡ 7fb18296-946b-4d51-9ed0-17ca134f5270
TableOfContents()

# ╔═╡ 3bd14493-23e8-4343-a370-57f53c4d2240
densities = ["low", "normal"]

# ╔═╡ 9301179c-4485-41b1-9971-62669707e312
sizes = ["small", "medium", "large"]

# ╔═╡ 94fc52bf-318d-426c-9b3c-bac271a1e20d
energies_motion = ["80-motion-3.npy", "100-motion-3.npy", "120-motion-3.npy", "135-motion-3.npy"]

# ╔═╡ d4251bb1-f785-4221-b5a1-63063f3d9b0e
slices = [1, 2, 3, 4, 5, 6]

# ╔═╡ 1e388da5-8aa5-447d-b6c4-a7d9776de3fd
root = joinpath(dirname(pwd()), "images_reproducibility1")

# ╔═╡ c7d71e90-ff42-43b2-a2dc-70c79de97d12
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

# ╔═╡ edab1de7-56bf-494f-8193-442011f5aef2
md"""
## Check DICOM image(s)
"""

# ╔═╡ 2decbb68-e600-47ba-adf2-38456f74580a
path = "/Users/daleblack/Google Drive/dev/MolloiLab/cac-simulation/images_reproducibility1/large/normal/80-motion"

# ╔═╡ c5dc9e05-0e80-4214-8d1b-f0553ff3f3cb
dcmdir_combined = dcmdir_parse(path);

# ╔═╡ 775e5949-3451-4b64-a6e2-15d5b1545bbc
vol_combined = load_dcm_array(dcmdir_combined);

# ╔═╡ ce993af3-8736-4e02-ab67-13fe91b72c96
@bind c PlutoUI.Slider(1:size(vol_combined, 3); default=1, show_value=true)

# ╔═╡ 8fbf59a8-f469-490e-8df1-ddc462019254
heatmap(transpose(vol_combined[:, :, c]); colormap=:grays)

# ╔═╡ 17e6c3aa-99a9-4f33-9201-20f40dace26b


# ╔═╡ Cell order:
# ╠═5571cce0-224e-4584-881f-b136dad2d05b
# ╠═7fb18296-946b-4d51-9ed0-17ca134f5270
# ╠═3bd14493-23e8-4343-a370-57f53c4d2240
# ╠═9301179c-4485-41b1-9971-62669707e312
# ╠═94fc52bf-318d-426c-9b3c-bac271a1e20d
# ╠═d4251bb1-f785-4221-b5a1-63063f3d9b0e
# ╠═1e388da5-8aa5-447d-b6c4-a7d9776de3fd
# ╠═c7d71e90-ff42-43b2-a2dc-70c79de97d12
# ╟─edab1de7-56bf-494f-8193-442011f5aef2
# ╠═2decbb68-e600-47ba-adf2-38456f74580a
# ╠═c5dc9e05-0e80-4214-8d1b-f0553ff3f3cb
# ╠═775e5949-3451-4b64-a6e2-15d5b1545bbc
# ╟─ce993af3-8736-4e02-ab67-13fe91b72c96
# ╠═8fbf59a8-f469-490e-8df1-ddc462019254
# ╠═17e6c3aa-99a9-4f33-9201-20f40dace26b

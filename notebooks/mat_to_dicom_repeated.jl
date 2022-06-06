### A Pluto.jl notebook ###
# v0.19.8

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

# ╔═╡ 6b158f74-4722-4be8-a998-5d36ae6e0a18
# ╠═╡ show_logs = false
begin
	let
		using Pkg
		Pkg.activate(mktempdir())
		Pkg.Registry.update()
		Pkg.add("PlutoUI")
		Pkg.add("CairoMakie")
		Pkg.add("MAT")
		Pkg.add(url="https://github.com/JuliaHealth/DICOM.jl")
		Pkg.add(url="https://github.com/Dale-Black/DICOMUtils.jl")
	end

	using PlutoUI
	using CairoMakie
	using MAT
	using DICOM
	using DICOMUtils
end

# ╔═╡ 8cfd76db-b78f-4388-bca4-5b414583cf93
TableOfContents()

# ╔═╡ ea12d90d-5cb5-4312-b9f0-1b2d42323c2f
densities = ["low", "normal"]

# ╔═╡ 00f7256b-28fd-423c-9b8d-3a392aac52b8
sizes = ["small", "medium", "large"]

# ╔═╡ 073390f4-e457-496a-a970-5232891fbcae
energies = [80, 100, 120, 135]

# ╔═╡ d590a6c0-2078-4f59-971b-e5b9097685d8
files = [1, 2, 3]

# ╔═╡ ddd861df-312e-4780-9fa0-982b3ac38605
for density in densities
	for _size in sizes
		for energy in energies
			for file in files
				file_inserts = file + 3
				ENERGY = energy
				ROD = string("QRM", ENERGY, "rod.mat")
				VESSEL = string("QRM", ENERGY, "vessels.mat")
				BASE_PATH = string("/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/mat_files/new_exposure_repro1/", _size, "/", density, "/")

				path = string(BASE_PATH, ROD)
				vars1 = matread(path)
				array1 = vars1[string("I")]
				array1 = Int16.(round.(array1))

				path2 = string(BASE_PATH, VESSEL)
				vars2 = matread(path2)
				array2 = vars2[string("I")]
				array2 = Int16.(round.(array2))

				dcm_path = "/Users/daleblack/Google Drive/Datasets/Canon_Aquilion_One_Vision/Large_rep1/96E1EB4F"
				
				dcm = dcm_parse(dcm_path)
				dcm[tag"Pixel Data"] = array1
				dcm[tag"Instance Number"] = file
				dcm[tag"Rows"] = size(array1, 1)
				dcm[tag"Columns"] = size(array1, 2)

				output_root1 = string("/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/images_reproducibility1/")
				if !isdir(output_root1)
					mkdir(output_root1)
				end
				output_root2 = string("/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/images_reproducibility1/", _size)
				if !isdir(output_root2)
					mkdir(output_root2)
				end
				output_root3 = string("/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/images_reproducibility1/", _size, "/", density)
				if !isdir(output_root3)
					mkdir(output_root3)
				end
				output_root = string("/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/images_reproducibility1/", _size, "/", density, "/", ENERGY)
				if !isdir(output_root)
					mkdir(output_root)
				end
				
				output_path_rod = string(output_root, "/", file, ".dcm")
				dcm_write(output_path_rod, dcm)

				dcm2 = dcm_parse(dcm_path)
				dcm2[tag"Pixel Data"] = array2
				dcm2[tag"Instance Number"] = file_inserts
				dcm2[tag"Rows"] = size(array2, 1)
				dcm2[tag"Columns"] = size(array2, 2)

				output_path_inserts = string(output_root, "/", file_inserts, ".dcm")
				dcm_write(output_path_inserts, dcm2)

			end
		end
	end
end

# ╔═╡ 93be84e2-d02a-4c39-b5c2-3850d6ab208e
md"""
## Check DICOM image(s)
"""

# ╔═╡ 23d5c197-a07f-4d5d-9227-2f080a95d325
path = "/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/images_reproducibility1/large/low/80"

# ╔═╡ 26be7a8d-84df-4e1d-a362-fbcadad80dbd
dcmdir_combined = dcmdir_parse(path);

# ╔═╡ 953e417e-351a-42fc-b4c0-75490179b44b
vol_combined = load_dcm_array(dcmdir_combined);

# ╔═╡ 2404e806-88ce-4941-a2e1-fc0d471477bb
@bind c PlutoUI.Slider(1:size(vol_combined, 3); default=1, show_value=true)

# ╔═╡ fbf72446-9b2a-45b6-94a7-7f53e2d6e817
heatmap(transpose(vol_combined[:, :, c]), colormap=:grays)

# ╔═╡ Cell order:
# ╠═6b158f74-4722-4be8-a998-5d36ae6e0a18
# ╠═8cfd76db-b78f-4388-bca4-5b414583cf93
# ╠═ea12d90d-5cb5-4312-b9f0-1b2d42323c2f
# ╠═00f7256b-28fd-423c-9b8d-3a392aac52b8
# ╠═073390f4-e457-496a-a970-5232891fbcae
# ╠═d590a6c0-2078-4f59-971b-e5b9097685d8
# ╠═ddd861df-312e-4780-9fa0-982b3ac38605
# ╟─93be84e2-d02a-4c39-b5c2-3850d6ab208e
# ╠═23d5c197-a07f-4d5d-9227-2f080a95d325
# ╠═26be7a8d-84df-4e1d-a362-fbcadad80dbd
# ╠═953e417e-351a-42fc-b4c0-75490179b44b
# ╟─2404e806-88ce-4941-a2e1-fc0d471477bb
# ╠═fbf72446-9b2a-45b6-94a7-7f53e2d6e817

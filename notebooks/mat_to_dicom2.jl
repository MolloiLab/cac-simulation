### A Pluto.jl notebook ###
# v0.18.1

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

# ╔═╡ 29fd80b0-3ba4-4e2a-9993-e19a74a237ec
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

# ╔═╡ eac4fea9-1967-4552-ac4f-b7f857bbc65e
TableOfContents()

# ╔═╡ ccf1d25d-7415-42b1-bf85-26335a7a9f3c
file = 3

# ╔═╡ 830f5666-6547-4fb9-9a4c-755564ea4bc1
file_inserts = file + 3

# ╔═╡ 324c688e-e63c-42a0-b6c0-75539554d7fa
begin
	ENERGY = 120
	ROD = string("QRM", ENERGY, "rod.mat")
	VESSEL = string("QRM", ENERGY, "vessels.mat")
	BASE_PATH = "/Users/daleblack/Desktop/QRM images new/"
end

# ╔═╡ 960f0da6-fd54-4456-8969-6ab440850b3a
md"""
## Convert .mat to Array
"""

# ╔═╡ 14fec7ab-4349-49be-80d8-aba636ce5510
begin
	path = string(BASE_PATH, ROD)
	vars1 = matread(path)
	array1 = vars1[string("test_rod")]
	array1 = Int16.(round.(array1))
end

# ╔═╡ e0744c9e-b260-4a32-90e5-c8da8ad417a0
heatmap(transpose(array1), colormap=:grays)

# ╔═╡ 32f16d8e-7049-4978-bbf2-f506ce07aa84
begin
	path2 = string(BASE_PATH, VESSEL)
	vars2 = matread(path2)
	array2 = vars2[string("test_phantom")]
	array2 = Int16.(round.(array2))
end

# ╔═╡ d6eaa14b-bf96-453b-a670-84127f3eb2fb
heatmap(transpose(array2), colormap=:grays)

# ╔═╡ a91105e0-bb2e-48de-a309-ef96ef347a62
md"""
## Create DICOM rod image
"""

# ╔═╡ 40ca726c-5ccc-4c20-bfe5-65a6c0ad4334
dcm_path = "/Users/daleblack/Google Drive/Datasets/Canon_Aquilion_One_Vision/Large_rep1/96E1EB4F";

# ╔═╡ 00c5c577-765d-4f5c-a0a9-7050bf293d12
begin
	dcm = dcm_parse(dcm_path)
	dcm[tag"Pixel Data"] = array1
	dcm[tag"Instance Number"] = file
	dcm[tag"Rows"] = size(array1, 1)
	dcm[tag"Columns"] = size(array1, 2)
	# dcm[tag"Pixel Spacing"] = [0.28, 0.28]
	# dcm[tag"Slice Thickness"] = [0.5]
end;

# ╔═╡ bff73dcc-ec55-48dc-ba4f-6fea6a85c69a
begin
	output_root = string("/Users/daleblack/Google Drive/Datasets/Simulated/", ENERGY, "new")
	if !isdir(output_root)
		mkdir(output_root)
	end
	output_path_rod = string(output_root, "/", file, ".dcm")
end

# ╔═╡ d1ad70c9-3dbc-4530-8bf2-0c6a5c5cfac5
dcm_write(output_path_rod, dcm)

# ╔═╡ 95ecbc6b-ec3d-4c37-88ad-39f08a3b5fa8
md"""
## Create DICOM inserts image
"""

# ╔═╡ 06af8de1-8e81-4506-b1f4-26ab02446993
begin
	dcm2 = dcm_parse(dcm_path)
	dcm2[tag"Pixel Data"] = array2
	dcm2[tag"Instance Number"] = file_inserts
	dcm2[tag"Rows"] = size(array2, 1)
	dcm2[tag"Columns"] = size(array2, 2)
	# dcm2[tag"Pixel Spacing"] = [0.28, 0.28]
	# dcm2[tag"Slice Thickness"] = [0.5]
end;

# ╔═╡ 236e3246-e16c-40e1-a0a1-f237fab4bd0b
output_path_inserts = string(output_root, "/", file_inserts, ".dcm");

# ╔═╡ f42f30bc-e830-4191-8ecc-75d12200f607
dcm_write(output_path_inserts, dcm2)

# ╔═╡ 6689d372-ad5d-48d5-92d1-d98cd7620e85
md"""
## Check DICOM image(s)
"""

# ╔═╡ 22cefe80-4ff6-4790-b06a-e777794a714c
# begin
# 	pth_root = "/Users/daleblack/Google Drive/Datasets/Simulated/"
# 	pth = string(pth_root, "combined")
# 	if ~isdir(pth)
# 		mkdir(pth)
# 	end
	
# 	pth_inserts = string(pth_root, inserts)
# 	files_inserts = readdir(pth_inserts)
# 	for i in files_inserts
# 		string(pth_inserts, )
# end

# ╔═╡ 02a0def9-910f-4d42-8a2a-248ad629ceb8
dcmdir_combined = dcmdir_parse(output_root);

# ╔═╡ 80efa582-15c8-4eda-a3f4-d292e017b770
vol_combined = load_dcm_array(dcmdir_combined);

# ╔═╡ dfe4b803-0b86-47e2-aa09-a2803f45aae5
@bind c PlutoUI.Slider(1:size(vol_combined, 3); default=1, show_value=true)

# ╔═╡ 7be0048c-21e2-4daf-84eb-3f2d35dc4c05
heatmap(transpose(vol_combined[:, :, c]), colormap=:grays)

# ╔═╡ Cell order:
# ╠═29fd80b0-3ba4-4e2a-9993-e19a74a237ec
# ╠═eac4fea9-1967-4552-ac4f-b7f857bbc65e
# ╠═ccf1d25d-7415-42b1-bf85-26335a7a9f3c
# ╠═830f5666-6547-4fb9-9a4c-755564ea4bc1
# ╠═324c688e-e63c-42a0-b6c0-75539554d7fa
# ╠═960f0da6-fd54-4456-8969-6ab440850b3a
# ╠═14fec7ab-4349-49be-80d8-aba636ce5510
# ╠═e0744c9e-b260-4a32-90e5-c8da8ad417a0
# ╠═32f16d8e-7049-4978-bbf2-f506ce07aa84
# ╠═d6eaa14b-bf96-453b-a670-84127f3eb2fb
# ╟─a91105e0-bb2e-48de-a309-ef96ef347a62
# ╠═40ca726c-5ccc-4c20-bfe5-65a6c0ad4334
# ╠═00c5c577-765d-4f5c-a0a9-7050bf293d12
# ╠═bff73dcc-ec55-48dc-ba4f-6fea6a85c69a
# ╠═d1ad70c9-3dbc-4530-8bf2-0c6a5c5cfac5
# ╠═95ecbc6b-ec3d-4c37-88ad-39f08a3b5fa8
# ╠═06af8de1-8e81-4506-b1f4-26ab02446993
# ╠═236e3246-e16c-40e1-a0a1-f237fab4bd0b
# ╠═f42f30bc-e830-4191-8ecc-75d12200f607
# ╠═6689d372-ad5d-48d5-92d1-d98cd7620e85
# ╠═22cefe80-4ff6-4790-b06a-e777794a714c
# ╠═02a0def9-910f-4d42-8a2a-248ad629ceb8
# ╠═80efa582-15c8-4eda-a3f4-d292e017b770
# ╠═dfe4b803-0b86-47e2-aa09-a2803f45aae5
# ╠═7be0048c-21e2-4daf-84eb-3f2d35dc4c05

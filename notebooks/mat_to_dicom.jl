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

# ╔═╡ 52c184d2-9a5e-11ec-0e10-976cdfa5f253
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

# ╔═╡ 5abcc458-7733-4ce4-a56f-bb7cb8882540
TableOfContents()

# ╔═╡ 1417f262-16f4-41f7-bed6-c895fb49fa80
file = 3

# ╔═╡ 3c5bfe99-3098-43b6-9d0e-0b3ed5884dea
file_inserts = file + 3

# ╔═╡ c1fde068-612a-45c2-a116-d6af918c721a
begin
	ENERGY = 135
	ROD = string("QRM", ENERGY, "rod.mat")
	VESSEL = string("QRM", ENERGY, "vessels.mat")
	BASE_PATH = "/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/mat_files/small/low/"
end

# ╔═╡ 2ae4e5ff-3cfd-4c80-a59a-e91c52480ecf
md"""
## Convert .mat to Array
"""

# ╔═╡ 0fab4688-9f2c-4df6-9cff-d3f9086e7a4b
begin
	path = string(BASE_PATH, ROD)
	vars1 = matread(path)
	array1 = vars1[string("r", ENERGY)]
	array1 = Int16.(round.(array1))
end;

# ╔═╡ 2cd3beb0-ce9f-461e-803c-59eb20b02dee
heatmap(transpose(array1), colormap=:grays)

# ╔═╡ 506a5e6b-94ed-4502-af57-0092eba7d817
begin
	path2 = string(BASE_PATH, VESSEL)
	vars2 = matread(path2)
	array2 = vars2[string("v", ENERGY, "_r")]
	array2 = Int16.(round.(array2))
end;

# ╔═╡ 65c27fd3-3959-4d87-a120-2be6a221ca87
heatmap(transpose(array2), colormap=:grays)

# ╔═╡ 42bdac01-ab21-43ea-928e-04231f0428ba
md"""
## Create DICOM rod image
"""

# ╔═╡ 7ef29728-2fd5-4c70-ae42-4d35ce6b45da
dcm_path = "/Users/daleblack/Google Drive/Datasets/Canon_Aquilion_One_Vision/Large_rep1/96E1EB4F";

# ╔═╡ a70027f1-5fbf-4e9d-8d09-4dd8718fc081
begin
	dcm = dcm_parse(dcm_path)
	dcm[tag"Pixel Data"] = array1
	dcm[tag"Instance Number"] = file
	dcm[tag"Rows"] = size(array1, 1)
	dcm[tag"Columns"] = size(array1, 2)
	# dcm[tag"Pixel Spacing"] = [0.28, 0.28]
	# dcm[tag"Slice Thickness"] = [0.5]
end;

# ╔═╡ 4062b5d4-9fe7-41f5-92e3-fb19cb6a1274
begin
	output_root = string("/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/images/small/low/", ENERGY)
	if !isdir(output_root)
		mkdir(output_root)
	end
	output_path_rod = string(output_root, "/", file, ".dcm")
end

# ╔═╡ 042bc6e8-edd2-4d56-b3e4-a0ec6ee94ffe
dcm_write(output_path_rod, dcm)

# ╔═╡ 68ab4f72-1d54-444d-b8a8-7d9b96108f36
md"""
## Create DICOM inserts image
"""

# ╔═╡ 6b558de5-751d-4fda-b6ab-c90ca4f5a6e6
begin
	dcm2 = dcm_parse(dcm_path)
	dcm2[tag"Pixel Data"] = array2
	dcm2[tag"Instance Number"] = file_inserts
	dcm2[tag"Rows"] = size(array2, 1)
	dcm2[tag"Columns"] = size(array2, 2)
	# dcm2[tag"Pixel Spacing"] = [0.28, 0.28]
	# dcm2[tag"Slice Thickness"] = [0.5]
end;

# ╔═╡ 3371dd16-6a49-430f-afad-7d94c9bc63ad
output_path_inserts = string(output_root, "/", file_inserts, ".dcm");

# ╔═╡ b822d070-167e-4ce4-a876-57ba142e2751
dcm_write(output_path_inserts, dcm2)

# ╔═╡ f78c8cb1-6f5a-439a-873f-68c71c95516f
md"""
## Check DICOM image(s)
"""

# ╔═╡ 719a8d5c-ef3a-4859-9fc9-7c02c911c71d
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

# ╔═╡ 93d48ee0-83e6-46cd-8567-16273bc269f9
dcmdir_combined = dcmdir_parse(output_root);

# ╔═╡ 81f0b564-52ff-4078-85a7-765375ad99f0
vol_combined = load_dcm_array(dcmdir_combined);

# ╔═╡ a79b50b9-7efb-49b8-8970-ecce09a284d4
@bind c PlutoUI.Slider(1:size(vol_combined, 3); default=1, show_value=true)

# ╔═╡ ec615641-95ee-43f8-9433-587ca30f3551
heatmap(transpose(vol_combined[:, :, c]), colormap=:grays)

# ╔═╡ Cell order:
# ╠═52c184d2-9a5e-11ec-0e10-976cdfa5f253
# ╠═5abcc458-7733-4ce4-a56f-bb7cb8882540
# ╠═1417f262-16f4-41f7-bed6-c895fb49fa80
# ╠═3c5bfe99-3098-43b6-9d0e-0b3ed5884dea
# ╠═c1fde068-612a-45c2-a116-d6af918c721a
# ╟─2ae4e5ff-3cfd-4c80-a59a-e91c52480ecf
# ╠═0fab4688-9f2c-4df6-9cff-d3f9086e7a4b
# ╠═2cd3beb0-ce9f-461e-803c-59eb20b02dee
# ╠═506a5e6b-94ed-4502-af57-0092eba7d817
# ╠═65c27fd3-3959-4d87-a120-2be6a221ca87
# ╟─42bdac01-ab21-43ea-928e-04231f0428ba
# ╠═7ef29728-2fd5-4c70-ae42-4d35ce6b45da
# ╠═a70027f1-5fbf-4e9d-8d09-4dd8718fc081
# ╠═4062b5d4-9fe7-41f5-92e3-fb19cb6a1274
# ╠═042bc6e8-edd2-4d56-b3e4-a0ec6ee94ffe
# ╟─68ab4f72-1d54-444d-b8a8-7d9b96108f36
# ╠═6b558de5-751d-4fda-b6ab-c90ca4f5a6e6
# ╠═3371dd16-6a49-430f-afad-7d94c9bc63ad
# ╠═b822d070-167e-4ce4-a876-57ba142e2751
# ╟─f78c8cb1-6f5a-439a-873f-68c71c95516f
# ╠═719a8d5c-ef3a-4859-9fc9-7c02c911c71d
# ╠═93d48ee0-83e6-46cd-8567-16273bc269f9
# ╠═81f0b564-52ff-4078-85a7-765375ad99f0
# ╟─a79b50b9-7efb-49b8-8970-ecce09a284d4
# ╠═ec615641-95ee-43f8-9433-587ca30f3551

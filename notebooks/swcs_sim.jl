### A Pluto.jl notebook ###
# v0.19.4

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

# ╔═╡ 98114689-d1cc-4aa9-b88f-c048059e95d3
# ╠═╡ show_logs = false
begin
	let
		using Pkg
		Pkg.activate(mktempdir())
		Pkg.Registry.update()
		Pkg.add("PlutoUI")
		Pkg.add("CairoMakie")
		Pkg.add("Statistics")
		Pkg.add("Images")
		Pkg.add("ImageMorphology")
		Pkg.add("ImageFiltering")
		Pkg.add("CSV")
		Pkg.add("DataFrames")
		Pkg.add(url="https://github.com/JuliaHealth/DICOM.jl")
		Pkg.add(url="https://github.com/Dale-Black/DICOMUtils.jl")
		Pkg.add(url="https://github.com/Dale-Black/PhantomSegmentation.jl")
		Pkg.add(url="https://github.com/Dale-Black/CalciumScoring.jl")
		Pkg.add("Distributions")
		Pkg.add("DSP")
	end
	
	using PlutoUI
	using CairoMakie
	using Statistics
	using Images
	using ImageMorphology
	using ImageFiltering
	using CSV
	using DataFrames
	using DICOM
	using DICOMUtils
	using PhantomSegmentation
	using CalciumScoring
	using Distributions
	using DSP
end

# ╔═╡ 33bd4068-abd8-45d4-babd-7b397310a645
TableOfContents()

# ╔═╡ 92ac2129-98f8-4dad-a446-6d3664aa783e
md"""
## Load DICOMS
"""

# ╔═╡ c02953cc-0020-452e-938b-694c7f7de44e
begin
	SCAN_NUMBER = 1
	VENDER = "80"
	# SIZE = "small"
	# SIZE = "medium"
	SIZE = "large"
	DENSITY = "low"
	# DENSITY = "normal"
	TYPE = "agatston"
	BASE_PATH = string("/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/images/", SIZE, "/", DENSITY, "/")
end

# ╔═╡ bdc233e2-a17f-42cf-a782-0b040c5b870d
md"""
**Everything below should be automatic, just scroll through to visually inspect that things make sense**
"""

# ╔═╡ e0bf7c37-4de9-4bd4-802d-c9df53dfa35c
root_path = string(BASE_PATH, VENDER)

# ╔═╡ d29cfbc0-aa1d-475f-8d5b-6572ae96bb10
dcm_path_list = dcm_list_builder(root_path)

# ╔═╡ 19220269-909e-4de7-a284-c6bb702922f6
pth = dcm_path_list[SCAN_NUMBER]

# ╔═╡ f243324f-4275-4466-b492-c51ca0986e49
scan = basename(pth)

# ╔═╡ 49d5afb9-3d52-418f-8c50-380660872184
header, dcm_array, slice_thick_ori1 = dcm_reader(pth);

# ╔═╡ b6f33cae-29c6-4d13-9b0c-a2551bd1e08e
md"""
## Helper Functions
"""

# ╔═╡ 818788f1-fb7b-4adc-9111-1f5cca892f86
function collect_tuple(tuple_array)
	row_num = size(tuple_array)
	col_num = length(tuple_array[1])
	container = zeros(Int64, row_num..., col_num)
	for i in 1:length(tuple_array)
		container[i,:] = collect(tuple_array[i])
	end
	return container
end

# ╔═╡ 78e52610-fae1-44fa-b931-275434a3ba73
function overlay_mask_bind(mask)
	indices = findall(x -> x == 1, mask)
	indices = Tuple.(indices)
	label_array = collect_tuple(indices)
	zs = unique(label_array[:,3])
	return PlutoUI.Slider(1:length(zs), default=3, show_value=true)
end

# ╔═╡ e1147d97-3f08-4a2f-ae30-ed3cb82df965
function overlay_mask_plot(array, mask, var, title::AbstractString)
	indices = findall(x -> x == 1, mask)
	indices = Tuple.(indices)
	label_array = collect_tuple(indices)
	zs = unique(label_array[:,3])
	indices_lbl = findall(x -> x == zs[var], label_array[:,3])
	
	fig = Figure()
	ax = Makie.Axis(fig[1, 1])
	ax.title = title
	heatmap!(array[:, :, zs[var]], colormap=:grays)
	scatter!(label_array[:, 1][indices_lbl], label_array[:, 2][indices_lbl], markersize=1, color=:red)
	fig
end

# ╔═╡ 9466ff87-9cef-4f4e-905c-d19070dd6890
function show_matrix(A::Matrix, red::Union{BitMatrix,Matrix{Bool}}=zeros(Bool, size(A)))
	base = RGB.(Gray.(A))

	base[red] .= RGB(1.0, 0.1, 0.1)

	# some tricks to show the pixels more clearly:
	s = max(size(A)...)
	if s >= 20
		min_size = 1200
		factor = min(5,min_size ÷ s)
		
		kron(base, ones(factor, factor))
	else
		base
	end
end

# ╔═╡ 1c2fa6f8-201d-40ae-8c03-7e5c5893fc1e
function create_mask(array, mask)
	@assert size(array) == size(mask)
	idxs = findall(x -> x == true, mask)
	overlayed_mask = zeros(size(array))
	for idx in idxs
		overlayed_mask[idx] = array[idx]
	end
	return overlayed_mask
end

# ╔═╡ a4ae678d-a7f1-4bfa-8d22-554583683089
md"""
## Segment Heart
"""

# ╔═╡ d8b8d96e-16ac-4412-b049-6287a933730c
masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3)÷2);

# ╔═╡ 043cffd4-1cc1-4a73-9fef-0a4519b45ee4
center_insert

# ╔═╡ 0e91e32b-aa3a-4b3b-a4b9-7b8fbaab132b
@bind a PlutoUI.Slider(1:size(masked_array, 3), default=10, show_value=true)

# ╔═╡ 1ec97def-1513-4309-a6e8-19d7ecfddc25
heatmap(masked_array[:, :, a], colormap=:grays)

# ╔═╡ fdad818f-1aff-44f7-8fd8-d625b7b873dc
begin
	fig = Figure()
	
	ax = Makie.Axis(fig[1, 1])
	ax.title = "Raw DICOM Array"
	heatmap!(transpose(dcm_array[:, :, 4]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2], center_insert[1]:center_insert[1], markersize=10, color=:red)
	fig
end

# ╔═╡ ce1cb4e1-9515-463b-b5c7-79c7bc295f49
begin
	fig2 = Figure()
	
	ax2 = Makie.Axis(fig2[1, 1])
	ax2.title = "Mask Array"
	heatmap!(transpose(mask), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2], center_insert[1]:center_insert[1], markersize=10, color=:red)
	fig2
end

# ╔═╡ 7d341908-e452-4975-bed6-d7f1d3e858af
begin
	fig3 = Figure()
	
	ax3 = Makie.Axis(fig3[1, 1])
	ax3.title = "Masked DICOM Array"
	heatmap!(transpose(masked_array[:, :, 5]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2], center_insert[1]:center_insert[1], markersize=10, color=:red)
	fig3
end

# ╔═╡ 1de1ae0b-5c08-4035-be82-53a9fdf2a7a9
md"""
## Segment Calcium Rod
"""

# ╔═╡ 47081781-ad0c-4ebb-af72-ee6ea6ccf335
begin
	global thresh
	if DENSITY == "low" && SIZE == "small"
		thresh = 60
	elseif DENSITY == "low" && SIZE == "large" && (VENDER == "120" || VENDER == "135")
		thresh = 75
	elseif DENSITY == "low" && SIZE == "medium"
		thresh = 75
	elseif DENSITY == "low"
		thresh = 60
	elseif DENSITY ==  "normal"
		thresh = 130
	end
end

# ╔═╡ daa65953-ac65-4689-bb28-344e822b2ba2
calcium_image, slice_CCI, quality_slice, cal_rod_slice = mask_rod(masked_array, header; calcium_threshold=thresh);

# ╔═╡ e1a9f2b8-be8d-4119-9a4a-e1c73e74b33e
@bind c PlutoUI.Slider(1:size(calcium_image, 3), default=cal_rod_slice, show_value=true)

# ╔═╡ 3ef8300d-56a6-4d72-9287-5e4eda72ee24
heatmap(transpose(calcium_image[:, :, c]), colormap=:grays)

# ╔═╡ 8351f18e-7fc7-4c3b-8688-cae0a1e6ec3b
array_filtered = abs.(mapwindow(median, calcium_image[:, :, 2], (3, 3)));

# ╔═╡ d2bf1ee1-fe3e-4ea1-b8c4-1218ac81ff91
bool_arr = array_filtered .> 0;

# ╔═╡ 08527140-2604-4685-a294-43ee44d4aa22
bool_arr_erode = ((erode(erode(erode(bool_arr)))));

# ╔═╡ 436ee3b5-a140-4a7a-a50b-1d90eeb04df6
heatmap(transpose(bool_arr), colormap=:grays)

# ╔═╡ 6ad451a3-6693-4355-8346-9469f7d94167
heatmap(transpose(bool_arr_erode), colormap=:grays)

# ╔═╡ 67cf1029-ebea-4b48-9168-b0d0b491a277
c_img = calcium_image[:, :, 1:3];

# ╔═╡ da901979-3472-48c7-81d5-ff36c1cb0366
begin
	mask_cal_3D = Array{Bool}(undef, size(c_img))
	for z in 1:size(c_img, 3)
		mask_cal_3D[:, :, z] = bool_arr_erode
	end
end;

# ╔═╡ 28445fd5-9798-4ace-bc96-ab57a402d56e
# hist(c_img[mask_cal_3D] .- offset)

# ╔═╡ 0bbac64e-c842-4dc4-8062-8c99cc3f029a
# calibration = c_img[mask_cal_3D] .- offset

# ╔═╡ 3a7a7501-30ef-4952-a2b7-1628f1af520e
# quantile!(c_img[mask_cal_3D] .- offset, 0.05:0.95)

# ╔═╡ 13948415-9b64-4ec0-b70a-8dd991bb4e9c
# mean(calibration), std(calibration)

# ╔═╡ b441160a-84fd-49d0-afae-1b63248cd90f
md"""
## Segment Calcium Inserts
"""

# ╔═╡ 325edc29-0752-432f-a1bc-ecb7976667ad
# mask_L_HD, mask_M_HD, mask_S_HD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_LD, mask_M_LD, mask_S_LD = mask_inserts_simulation(
#             dcm_array, masked_array, header, slice_CCI, center_insert
# );

# ╔═╡ 7208bc08-acad-44f7-854c-cea5d0190b24
begin 
	root_new = string("/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/julia_arrays/", SIZE, "/") 
	mask_L_HD = Array(CSV.read(string(root_new, "mask_L_HD.csv"), DataFrame; header=false))
	mask_M_HD = Array(CSV.read(string(root_new, "mask_M_HD.csv"), DataFrame; header=false))
	mask_S_HD = Array(CSV.read(string(root_new, "mask_S_HD.csv"), DataFrame; header=false))
	mask_L_MD = Array(CSV.read(string(root_new, "mask_L_MD.csv"), DataFrame;
	header=false))
	mask_M_MD = Array(CSV.read(string(root_new, "mask_M_MD.csv"), DataFrame; header=false))
	mask_S_MD = Array(CSV.read(string(root_new, "mask_S_MD.csv"), DataFrame; header=false))
	mask_L_LD = Array(CSV.read(string(root_new, "mask_L_LD.csv"), DataFrame; header=false))
	mask_M_LD = Array(CSV.read(string(root_new, "mask_M_LD.csv"), DataFrame; header=false))
	mask_S_LD = Array(CSV.read(string(root_new, "mask_S_LD.csv"), DataFrame; header=false))
end;

# ╔═╡ c3b659fe-b9bc-44d6-8cf1-5ad18c4fa346
md"""
## Mass cal factor
"""

# ╔═╡ 84df1941-72a6-4cf5-9e6e-7f79d2534dd6
output = calc_output(masked_array, header, 5, thresh, trues(3, 3));

# ╔═╡ 8b601177-be58-4f32-a32f-2bc675c6364a
heatmap(output[2])

# ╔═╡ bd6a4dee-1055-4c14-bafd-74b64daa0715
insert_centers = calc_centers(dcm_array, output, header, center_insert, slice_CCI)

# ╔═╡ d0194709-f4f5-4f7a-acba-76b600dcff92
center_large_LD = insert_centers[:Large_LD]

# ╔═╡ 72dcf6e7-8fb8-47f7-bfaa-f5529210e8d0
rows, cols = Int(header[tag"Rows"]), Int(header[tag"Columns"])

# ╔═╡ 65c9df3c-ea01-48fa-b192-52b3f3b76f9b
md"""
# Score Large Inserts
"""

# ╔═╡ 1f353ddd-4e9a-4437-8d2e-97a4c466ee28
md"""
## High Density
"""

# ╔═╡ 10899fba-616b-4bfc-9d75-49383512e747
arr = masked_array[:, :, 4:6];

# ╔═╡ a228cbf8-4d79-41a4-b419-afe26b981a1c
begin
	mask_L_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_HD_3D[:, :, z] = dilate(dilate(mask_L_HD))
	end
end;

# ╔═╡ 334e9473-fc4f-466c-a046-6c66bf65b738
md"""
#### Dilated mask
"""

# ╔═╡ a2c59aa5-6c86-4d4a-aa84-c6675446083d
dilated_mask_L_HD = dilate(dilate(mask_L_HD_3D));

# ╔═╡ aa851337-3dd3-4521-a402-52a429bec1e6
@bind g2 overlay_mask_bind(dilated_mask_L_HD)

# ╔═╡ e853a3d8-0a18-43b7-949b-c374203d4838
overlay_mask_plot(arr, dilated_mask_L_HD, g2, "dilated mask")

# ╔═╡ 0373a338-9a93-42b8-98ee-9768a50dda67
pixel_size = DICOMUtils.get_pixel_size(header)

# ╔═╡ 03e64fbb-49b0-4b68-ae1a-c5d155651146
mass_cal_factor, angle_0_200HA, water_rod_metrics = mass_calibration(masked_array, center_large_LD, center_insert, 2, cols, rows, pixel_size)

# ╔═╡ 268ad486-9907-4aaa-a271-a77664fea675
overlayed_mask_l_hd = create_mask(arr, dilated_mask_L_HD);

# ╔═╡ 405f2f7e-123a-4ccf-9fc1-9b6eee4dae12
alg = Agatston()

# ╔═╡ 19c04849-68a6-4ea3-8a1b-e9183c9b852e
alg2 = SpatiallyWeighted()

# ╔═╡ b5686627-bfa4-4802-9122-eed55fb7db96
agat_l_hd, mass_l_hd = score(overlayed_mask_l_hd, pixel_size, mass_cal_factor, alg)

# ╔═╡ 385df0f4-bd29-4702-a699-f70ea0093d06
begin
	global μ, σ
	if VENDER == "80"
		μ, σ = 170, 30
	elseif VENDER == "100"
		μ, σ = 165, 30
	elseif VENDER == "120"
		μ, σ = 160, 30
	else VENDER == "135"
		μ, σ = 155, 30
	end
end

# ╔═╡ 05cf3e18-950e-427c-91d9-be289da5c5c5
swcs_l_hd = score(overlayed_mask_l_hd, μ, σ, alg2)

# ╔═╡ 6fa0cf17-0c5a-4285-ab2c-09c55dd2e31b
md"""
### Mass cal 2
"""

# ╔═╡ f6a0b43f-a56a-41d2-a89f-c1f17f0c6154
begin
	dens_calicium_rod = 200 # mg/cc
	mean_calcium_rod = mean(c_img[mask_cal_3D])
	calibration_factor = dens_calicium_rod / mean_calcium_rod
end

# ╔═╡ 5a09303d-5360-48c9-9554-c769e63b1ee0
md"""
## Medium Density
"""

# ╔═╡ 20a53747-25a1-4642-9043-64109ff28711
begin
	mask_L_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_MD_3D[:, :, z] = mask_L_MD
	end
end;

# ╔═╡ 5f0f4233-ae8b-432f-93e9-f2b3a54dc890
md"""
#### Dilated mask
"""

# ╔═╡ bd5aa941-8ab9-4b43-b6f5-f2320f1f9e03
dilated_mask_L_MD = dilate(dilate(mask_L_MD_3D));

# ╔═╡ 1d87f244-f6c5-4e12-9d0d-5c2e7c4805e5
@bind h2 overlay_mask_bind(dilated_mask_L_MD)

# ╔═╡ 8f9e3da4-baab-431d-a4c7-67e753144aef
overlay_mask_plot(arr, dilated_mask_L_MD, h2, "dilated mask")

# ╔═╡ 31708cf8-1d56-4b16-ad4e-df356a8924fc
overlayed_mask_l_md = create_mask(arr, dilated_mask_L_MD);

# ╔═╡ cafed3ea-9f5a-4549-bc66-1181c1e7c6c8
agat_l_md, mass_l_md = score(overlayed_mask_l_md, pixel_size, mass_cal_factor, alg)

# ╔═╡ 2269beda-ae30-49ca-af6a-65a681dc7405
swcs_l_md = score(overlayed_mask_l_md, μ, σ, alg2)

# ╔═╡ 9d8e3d8d-a4e6-418c-a241-b4b1cd119a14
md"""
## Low Density
"""

# ╔═╡ 8f125f2e-ec57-40f0-96d6-e04a544db0c3
begin
	mask_L_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_LD_3D[:, :, z] = mask_L_LD
	end
end;

# ╔═╡ ad1c418f-a4b7-45fa-8798-534b5e474b99
md"""
#### Dilated mask
"""

# ╔═╡ f72f5020-6e3f-4b47-9340-55e5d687b56c
dilated_mask_L_LD = dilate(dilate(mask_L_LD_3D));

# ╔═╡ d3ce01ea-a5d6-41e3-8f98-e0042026da9b
@bind i2 overlay_mask_bind(dilated_mask_L_LD)

# ╔═╡ 2e85a77c-7bb0-47ba-949e-4fdccea6b8e8
overlay_mask_plot(arr, dilated_mask_L_LD, i2, "dilated mask")

# ╔═╡ db774a1b-dbda-42a3-abb2-89ae26a2a31c
overlayed_mask_l_ld = create_mask(arr, dilated_mask_L_LD);

# ╔═╡ a9ec77d3-27bc-4a0c-b74a-555b00e6f29a
agat_l_ld, mass_l_ld = score(overlayed_mask_l_ld, pixel_size, mass_cal_factor, alg)

# ╔═╡ 8b27b710-b790-441d-8541-17fae3a72b32
swcs_l_ld = score(overlayed_mask_l_ld, μ, σ, alg2)

# ╔═╡ 341e35e2-f280-47c6-9da5-0d2dc137bce0
md"""
# Score Medium Inserts
"""

# ╔═╡ 2b8e12f7-8504-4c28-8208-fd9afb52d1ae
md"""
## High Density
"""

# ╔═╡ 12240afd-a739-4c51-8570-23770faf48fd
begin
	mask_M_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_HD_3D[:, :, z] = mask_M_HD
	end
end;

# ╔═╡ 89eab716-678e-407e-acc5-3bd4a490aaf5
md"""
#### Dilated mask
"""

# ╔═╡ e12840cf-d45e-4f66-aba9-c178d46f744b
dilated_mask_M_HD = dilate(dilate(dilate(dilate(mask_M_HD_3D))));

# ╔═╡ ea52b48a-5b66-4d00-b12e-a41ed8a4bde1
@bind j2 overlay_mask_bind(dilated_mask_M_HD)

# ╔═╡ be4f9b9c-1f42-414e-b74b-f2df17811e8e
overlay_mask_plot(arr, dilated_mask_M_HD, j2, "dilated mask")

# ╔═╡ 84957e23-4d9c-409f-b1cf-ab0c8db1a1b6
overlayed_mask_m_hd = create_mask(arr, dilated_mask_M_HD);

# ╔═╡ 4c3d0757-f378-46ca-944a-cc19d56dd11b
agat_m_hd, mass_m_hd = score(overlayed_mask_m_hd, pixel_size, mass_cal_factor, alg)

# ╔═╡ 33b85892-7ce9-4392-a3dc-b0a69592dba9
swcs_m_hd = score(overlayed_mask_m_hd, μ, σ, alg2)

# ╔═╡ 3d5aa5e2-cfd9-4c5e-83fc-d855ea5c2f30
md"""
## Medium Density
"""

# ╔═╡ 602601ca-081a-4abc-8726-6007be5558cf
begin
	mask_M_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_MD_3D[:, :, z] = mask_M_MD
	end
end;

# ╔═╡ ac3acca8-868d-47a4-a14a-9696876fdc78
md"""
#### Dilated mask
"""

# ╔═╡ 28d97f1f-c865-4f15-a7a4-0224039a6d83
dilated_mask_M_MD = dilate(dilate(dilate(dilate(mask_M_MD_3D))));

# ╔═╡ b7d25c42-dc50-406c-93d6-fd5aec863c10
@bind k2 overlay_mask_bind(dilated_mask_M_MD)

# ╔═╡ 80c2adfa-f6b6-414e-bf68-74f807b992ff
overlay_mask_plot(arr, dilated_mask_M_MD, k2, "dilated mask")

# ╔═╡ 830871a6-3f35-4c97-8525-7967b0fbb427
overlayed_mask_m_md = create_mask(arr, dilated_mask_M_MD);

# ╔═╡ 3aaac942-afae-4ffc-b1ea-6c023d6e932b
agat_m_md, mass_m_md = score(overlayed_mask_m_md, pixel_size, mass_cal_factor, alg)

# ╔═╡ 0011a165-5fb0-4964-8065-de0663249ee2
swcs_m_md = score(overlayed_mask_m_md, μ, σ, alg2)

# ╔═╡ d179a807-dcf1-4524-81e4-2071c4933f13
md"""
## Low Density
"""

# ╔═╡ 4a769949-2e00-4825-b852-d419149b47fd
begin
	mask_M_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_LD_3D[:, :, z] = mask_M_LD
	end
end;

# ╔═╡ 81bcc301-c5a9-4ac6-92d9-0f7dcd8676b8
md"""
#### Dilated mask
"""

# ╔═╡ c13b7546-2094-42fc-bc78-86fe8e6e5a59
dilated_mask_M_LD = dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))));

# ╔═╡ 30fc4df4-6e6f-42ac-bc4e-222efb20291d
@bind l2 overlay_mask_bind(dilated_mask_M_LD)

# ╔═╡ 4f4bc464-c398-472f-861a-dff6cf96388d
overlay_mask_plot(arr, dilated_mask_M_LD, l2, "dilated mask")

# ╔═╡ b4aab606-8d2e-4110-a0a6-cca477375a45
overlayed_mask_m_ld = create_mask(arr, dilated_mask_M_LD);

# ╔═╡ b02751b0-c28c-4d1b-976b-b15617ddc1f1
agat_m_ld, mass_m_ld = score(overlayed_mask_m_ld, pixel_size, mass_cal_factor, alg)

# ╔═╡ e5eb2459-faad-43d9-a702-cae081f5ae8f
swcs_m_ld = score(overlayed_mask_m_ld, μ, σ, alg2)

# ╔═╡ aeb23eb0-6ac0-4791-95ce-563bbd9f551b
md"""
# Score Small Inserts
"""

# ╔═╡ a39be3e4-cd38-4bea-ac3e-4d29c847a5e7
md"""
## High Density
"""

# ╔═╡ 38aa73f8-9ec2-43b8-983f-688c7caefa15
begin
	mask_S_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_HD_3D[:, :, z] = mask_S_HD
	end
end;

# ╔═╡ b1dda8e1-2e73-4e09-a73a-e151a066b99b
md"""
#### Dilated mask
"""

# ╔═╡ 2cbc8541-9c85-473e-af24-a723d609b590
dilated_mask_S_HD = dilate(dilate((dilate(dilate((mask_S_HD_3D))))));

# ╔═╡ 6aaa21ae-afeb-4099-9bb6-096ed8e9c093
@bind m2 overlay_mask_bind(dilated_mask_S_HD)

# ╔═╡ 9c5f120e-f9bb-4737-8031-95da4ea26e23
overlay_mask_plot(arr, dilated_mask_S_HD, m2, "dilated mask")

# ╔═╡ 8a409f2e-7ab8-47b3-a855-fb1a84c5881a
overlayed_mask_s_hd = create_mask(arr, dilated_mask_S_HD);

# ╔═╡ ce3b04f9-1461-4007-825f-7481f5da9312
agat_s_hd, mass_s_hd = score(overlayed_mask_s_hd, pixel_size, mass_cal_factor, alg)

# ╔═╡ 6cfe8bf9-7869-41e6-9df2-d2c18a159255
swcs_s_hd = score(overlayed_mask_s_hd, μ, σ, alg2)

# ╔═╡ 26609afb-14a5-4a5f-a90e-55fa9f3b17dd
md"""
## Medium Density
"""

# ╔═╡ ff7b8b70-e0ca-4c83-aa3a-ed6a5b60fc55
begin
	mask_S_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_MD_3D[:, :, z] = mask_S_MD
	end
end;

# ╔═╡ 60c507e8-e885-491e-a8b3-db61d61ee4fa
md"""
#### Dilated mask
"""

# ╔═╡ 513dfbb8-9c34-40f4-93c9-8ae21d5dffb6
dilated_mask_S_MD = dilate(dilate((dilate(dilate(mask_S_MD_3D)))));

# ╔═╡ bf1a2681-60f3-425f-b938-98934aba4588
@bind n2 overlay_mask_bind(dilated_mask_S_MD)

# ╔═╡ 435321e8-c027-496d-8280-f9aad31f3632
overlay_mask_plot(arr, dilated_mask_S_MD, n2, "dilated mask")

# ╔═╡ 673d4151-ec74-434e-8d0c-bfe162f4b6a0
overlayed_mask_s_md = create_mask(arr, dilated_mask_S_MD);

# ╔═╡ 2ba692fb-7a93-4521-b1a8-d334f8a8d763
agat_s_md, mass_s_md = score(overlayed_mask_s_md, pixel_size, mass_cal_factor, alg)

# ╔═╡ 78ce4c60-51c7-4c64-a3ce-d71d1529a8c7
swcs_s_md = score(overlayed_mask_s_md, μ, σ, alg2)

# ╔═╡ 4238b68f-b4fc-40b3-b859-3a5130834017
md"""
## Low Density
"""

# ╔═╡ bacc56de-f082-436e-bffa-23379fdb9c61
begin
	mask_S_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_LD_3D[:, :, z] = mask_S_LD
	end
end;

# ╔═╡ 3e321593-3949-41c6-a734-cfa05fcc471e
md"""
#### Dilated mask
"""

# ╔═╡ 58cc58d7-4706-41c4-a2f7-04d9cd6eb31d
dilated_mask_S_LD = dilate(dilate((dilate(dilate(mask_S_LD_3D)))));

# ╔═╡ 5d3a539f-9ff3-4cbb-8957-88035698f69f
@bind o2 overlay_mask_bind(dilated_mask_S_LD)

# ╔═╡ 55ae3b8f-54cc-4d4a-a3c5-37a1afad2c0b
overlay_mask_plot(arr, dilated_mask_S_LD, o2, "dilated mask")

# ╔═╡ bd53b5f1-73a2-4f65-aa20-a7cdba456ce0
overlayed_mask_s_ld = create_mask(arr, dilated_mask_S_LD);

# ╔═╡ 9fafb79b-0706-470d-8a24-f2682621e39c
agat_s_ld, mass_s_ld = score(overlayed_mask_s_ld, pixel_size, mass_cal_factor, alg)

# ╔═╡ c2116577-2f05-40dc-9094-92c621955f93
swcs_s_ld = score(overlayed_mask_s_ld, μ, σ, alg2)

# ╔═╡ 258caf22-6f46-4341-80eb-d572aa26ca59
md"""
# Results
"""

# ╔═╡ 3d573812-8bc2-4dce-99fb-6883b86433f4
density_array = [0, 200, 400, 800]

# ╔═╡ 52f3d56d-9219-43f3-ba32-925faa2e0e07
inserts = [
	"Low Density",
	"Medium Density",
	"High Density"
]

# ╔═╡ 260b441d-12eb-44d5-adb2-dc5401abdb20
md"""
## Agatston
"""

# ╔═╡ c0b060b6-ef4d-415e-9f8c-a39ab4373541
calculated_agat_large = [
	agat_l_ld,
	agat_l_md,
	agat_l_hd
]

# ╔═╡ 3f79ca52-b120-43af-958e-6a0e5edf079f
calculated_agat_medium = [
	agat_m_ld,
	agat_m_md,
	agat_m_hd
]

# ╔═╡ 0de07ba8-f1d6-4c01-86fd-82e42e9a29bb
calculated_agat_small = [
	agat_s_ld,
	agat_s_md,
	agat_s_hd
]

# ╔═╡ 40f4a8e5-8a97-43d6-8d83-4993d5f59b83
md"""
## Spatially Weighted
"""

# ╔═╡ e89bd455-71ba-41c0-870c-2dad225e129e
calculated_swcs_large = [
	swcs_l_ld,
	swcs_l_md,
	swcs_l_hd
]

# ╔═╡ 36130f3d-8a03-49ce-ab36-f9ff0b452c8d
calculated_swcs_medium = [
	swcs_m_ld,
	swcs_m_md,
	swcs_m_hd
]

# ╔═╡ e53a38e7-ba34-45af-bc39-86026763241f
calculated_swcs_small = [
	swcs_s_ld,
	swcs_s_md,
	swcs_s_hd
]

# ╔═╡ bdf7c9e5-765d-40b2-b4ae-19e068348639
md"""
## SWCS vs Agatston
"""

# ╔═╡ 14b79f5f-87c5-4db1-855b-607cec0946b0
md"""
## All
"""

# ╔═╡ 5757d60a-6944-4286-ad55-95ac78acd903
md"""
## Mass
"""

# ╔═╡ d10d9943-948e-4d91-881d-93ac9c8fdf5b
volume_gt = [
	7.065,
	63.585,
	176.625
]

# ╔═╡ 9b39f9e2-2080-4ab7-ba84-7e58b0047805
volume_gt

# ╔═╡ b9ac10ef-1ea5-46ea-8f9a-471b2a25380d
mass_cal_factor * volume_gt[3] * calibration_factor

# ╔═╡ 5c3f138e-d46e-4a3d-9bc4-e2ade6784907
ground_truth_mass_large = [
	volume_gt[3] * density_array[2] * 1e-3,
	volume_gt[3] * density_array[3] * 1e-3,
	volume_gt[3] * density_array[4] * 1e-3
] # mg

# ╔═╡ c8066892-8ea6-4383-8771-699eb1394a5c
calculated_mass_large = [
	mass_l_ld,
	mass_l_md,
	mass_l_hd
]

# ╔═╡ ee2bd00c-e949-48f9-8e3c-1e05eb749075
ground_truth_mass_medium = [
	volume_gt[2] * density_array[2] * 1e-3,
	volume_gt[2] * density_array[3] * 1e-3,
	volume_gt[2] * density_array[4] * 1e-3
]

# ╔═╡ 455c1094-548f-4983-918e-2a7eaa47ff00
calculated_mass_medium = [
	mass_m_ld,
	mass_m_md,
	mass_m_hd
]

# ╔═╡ acf1c58f-24c3-4b6a-a003-67774960112e
ground_truth_mass_small = [
	volume_gt[1] * density_array[2] * 1e-3,
	volume_gt[1] * density_array[3] * 1e-3,
	volume_gt[1] * density_array[4] * 1e-3
]

# ╔═╡ 5c94eae2-9e45-49c9-a123-fbd46932c95f
calculated_mass_small = [
	mass_s_ld,
	mass_s_md,
	mass_s_hd
]

# ╔═╡ 7ee1efe9-d6ac-4076-a5f0-8e18b5d64702
df = DataFrame(
	scan = scan,
	inserts = inserts,
	calculated_agat_large = calculated_agat_large,
	calculated_agat_medium = calculated_agat_medium,
	calculated_agat_small = calculated_agat_small,
	calculated_swcs_large = calculated_swcs_large,
	calculated_swcs_medium = calculated_swcs_medium,
	calculated_swcs_small = calculated_swcs_small,
	ground_truth_mass_large = ground_truth_mass_large,
	calculated_mass_large = calculated_mass_large,
	ground_truth_mass_medium = ground_truth_mass_medium,
	calculated_mass_medium = calculated_mass_medium,
	ground_truth_mass_small = ground_truth_mass_small,
	calculated_mass_small = calculated_mass_small,
	mass_cal_factor = mass_cal_factor
)

# ╔═╡ 59a17ced-27c0-415f-a95a-7a8af8aba7e3
begin
	fswcs222 = Figure()
	axswcs222 = Makie.Axis(fswcs222[1, 1])
	
	scatter!(df[!, :calculated_agat_large], df[!, :calculated_swcs_large], label="label")
	
	axswcs222.title = "SWCS vs Agatston(Large)"
	axswcs222.ylabel = "SWCS"
	axswcs222.xlabel = "Agatston"

	# xlims!(axswcs222, 0, 400)
	# ylims!(axswcs222, 0, 400)
	
	fswcs222[1, 2] = Legend(fswcs222, axswcs222, framevisible = false)
	
	fswcs222
end

# ╔═╡ adacf893-01ef-4cd3-bcc5-66b27e77e3be
begin
	fswcs223 = Figure()
	axswcs223 = Makie.Axis(fswcs223[1, 1])
	
	scatter!(df[!, :calculated_agat_medium], df[!, :calculated_swcs_medium],  label="label")
	
	axswcs223.title = "SWCS vs Agatston (Medium)"
	axswcs223.ylabel = "SWCS"
	axswcs223.xlabel = "Agatston"

	# xlims!(axswcs223, 0, 200)
	# ylims!(axswcs223, 0, 200)
	
	fswcs223[1, 2] = Legend(fswcs223, axswcs223, framevisible = false)
	
	fswcs223
end

# ╔═╡ b24231d2-aa2e-45ad-80f4-e875fa282cd0
begin
	fswcs224 = Figure()
	axswcs224 = Makie.Axis(fswcs224[1, 1])
	
	scatter!(df[!, :calculated_agat_small], df[!, :calculated_swcs_small], label="label")
	
	axswcs224.title = "SWCS vs Agatston (Small)"
	axswcs224.ylabel = "SWCS"
	axswcs224.xlabel = "Agatston"

	# xlims!(axswcs224, 0, 1)
	# ylims!(axswcs224, 0, 1)
	
	fswcs224[1, 2] = Legend(fswcs224, axswcs224, framevisible = false)
	
	fswcs224
end

# ╔═╡ e9d6989d-93fb-43ee-b18a-582dd637eed4
begin
	fmass22 = Figure()
	axmass22 = Makie.Axis(fmass22[1, 1])
	
	scatter!(density_array[2:end], df[!, :ground_truth_mass_large], label="ground_truth_mass_large")
	scatter!(density_array[2:end], df[!, :calculated_mass_large], label="calculated_mass_large")
	
	axmass22.title = "Mass Measurements (Large)"
	axmass22.ylabel = "Mass (mg)"
	axmass22.xlabel = "Density (mg/cm^3)"

	xlims!(axmass22, 0, 850)
	ylims!(axmass22, 0, 200)
	
	fmass22[1, 2] = Legend(fmass22, axmass22, framevisible = false)
	
	fmass22
end

# ╔═╡ 9e533473-6fd4-4155-9c35-d6901580b45f
begin
	fmass32 = Figure()
	axmass32 = Makie.Axis(fmass32[1, 1])
	
	scatter!(density_array[2:end], df[!, :ground_truth_mass_medium], label="ground_truth_mass_medium")
	scatter!(density_array[2:end], df[!, :calculated_mass_medium], label="calculated_mass_medium")
	
	axmass32.title = "Mass Measurements (Medium)"
	axmass32.ylabel = "Mass (mg)"
	axmass32.xlabel = "Density (mg/cm^3)"

	xlims!(axmass32, 0, 850)
	ylims!(axmass32, 0, 85)
	
	fmass32[1, 2] = Legend(fmass32, axmass32, framevisible = false)
	
	fmass32
end

# ╔═╡ 27be579e-87fc-4b77-833c-01c3361375af
begin
	fmass42 = Figure()
	axmass42 = Makie.Axis(fmass42[1, 1])
	
	scatter!(density_array[2:end], df[!, :ground_truth_mass_small], label="ground_truth_mass_small")
	scatter!(density_array[2:end], df[!, :calculated_mass_small], label="calculated_mass_small")
	
	axmass42.title = "Mass Measurements (Small)"
	axmass42.ylabel = "Mass (mg)"
	axmass42.xlabel = "Density (mg/cm^3)"

	xlims!(axmass42, 0, 850)
	ylims!(axmass42, 0, 10)
	
	fmass42[1, 2] = Legend(fmass42, axmass42, framevisible = false)
	
	fmass42
end

# ╔═╡ ade95452-1c8d-4888-9313-b092409ca609
md"""
### Save Results
"""

# ╔═╡ 844754c6-a0ae-4cc7-9d6d-38204d7b24f5
# if ~isdir(string(cd(pwd, "..") , "/output/", TYPE))
# 	mkdir(string(cd(pwd, "..") , "/output/", TYPE))
# end

# ╔═╡ c7b49892-0079-4c79-92d7-3717460efd0f
# output_path = string(cd(pwd, "..") , "/output/", TYPE, "/", scan, ".csv")

# ╔═╡ 49275899-e5ec-47f4-b4f7-43401529642f
md"""
### Save full df
"""

# ╔═╡ 6471de28-e427-4152-b778-74c3d208b5b1
dfs = []

# ╔═╡ b99951e3-6292-42e4-801e-4824bca11720
push!(dfs, df)

# ╔═╡ 92be24c9-e0a5-4539-9501-93a8d6c52f10
if length(dfs) == 12
	global new_df = vcat(dfs[1:12]...)
	output_path_new = string(cd(pwd, "..") , "/output/", TYPE, "/", "full.csv")
	CSV.write(output_path_new, new_df)
end

# ╔═╡ 1723c48d-341f-4898-9b57-4202664a9415
# output_path_new = string(cd(pwd, "..") , "/output/", TYPE, "/", "full.csv")

# ╔═╡ 4d1825ea-b5cd-4445-99d5-a218845e629c


# ╔═╡ 3912fb58-7431-40f7-9823-9a4c91487bf6


# ╔═╡ Cell order:
# ╠═98114689-d1cc-4aa9-b88f-c048059e95d3
# ╠═33bd4068-abd8-45d4-babd-7b397310a645
# ╟─92ac2129-98f8-4dad-a446-6d3664aa783e
# ╠═c02953cc-0020-452e-938b-694c7f7de44e
# ╟─bdc233e2-a17f-42cf-a782-0b040c5b870d
# ╠═e0bf7c37-4de9-4bd4-802d-c9df53dfa35c
# ╠═d29cfbc0-aa1d-475f-8d5b-6572ae96bb10
# ╠═19220269-909e-4de7-a284-c6bb702922f6
# ╠═f243324f-4275-4466-b492-c51ca0986e49
# ╠═49d5afb9-3d52-418f-8c50-380660872184
# ╟─b6f33cae-29c6-4d13-9b0c-a2551bd1e08e
# ╟─818788f1-fb7b-4adc-9111-1f5cca892f86
# ╟─78e52610-fae1-44fa-b931-275434a3ba73
# ╟─e1147d97-3f08-4a2f-ae30-ed3cb82df965
# ╟─9466ff87-9cef-4f4e-905c-d19070dd6890
# ╟─1c2fa6f8-201d-40ae-8c03-7e5c5893fc1e
# ╟─a4ae678d-a7f1-4bfa-8d22-554583683089
# ╠═d8b8d96e-16ac-4412-b049-6287a933730c
# ╠═043cffd4-1cc1-4a73-9fef-0a4519b45ee4
# ╟─0e91e32b-aa3a-4b3b-a4b9-7b8fbaab132b
# ╠═1ec97def-1513-4309-a6e8-19d7ecfddc25
# ╟─fdad818f-1aff-44f7-8fd8-d625b7b873dc
# ╟─ce1cb4e1-9515-463b-b5c7-79c7bc295f49
# ╟─7d341908-e452-4975-bed6-d7f1d3e858af
# ╟─1de1ae0b-5c08-4035-be82-53a9fdf2a7a9
# ╠═47081781-ad0c-4ebb-af72-ee6ea6ccf335
# ╠═daa65953-ac65-4689-bb28-344e822b2ba2
# ╟─e1a9f2b8-be8d-4119-9a4a-e1c73e74b33e
# ╠═3ef8300d-56a6-4d72-9287-5e4eda72ee24
# ╠═8351f18e-7fc7-4c3b-8688-cae0a1e6ec3b
# ╠═d2bf1ee1-fe3e-4ea1-b8c4-1218ac81ff91
# ╠═08527140-2604-4685-a294-43ee44d4aa22
# ╠═436ee3b5-a140-4a7a-a50b-1d90eeb04df6
# ╠═6ad451a3-6693-4355-8346-9469f7d94167
# ╠═67cf1029-ebea-4b48-9168-b0d0b491a277
# ╠═da901979-3472-48c7-81d5-ff36c1cb0366
# ╠═28445fd5-9798-4ace-bc96-ab57a402d56e
# ╠═0bbac64e-c842-4dc4-8062-8c99cc3f029a
# ╠═3a7a7501-30ef-4952-a2b7-1628f1af520e
# ╠═13948415-9b64-4ec0-b70a-8dd991bb4e9c
# ╟─b441160a-84fd-49d0-afae-1b63248cd90f
# ╠═325edc29-0752-432f-a1bc-ecb7976667ad
# ╠═7208bc08-acad-44f7-854c-cea5d0190b24
# ╟─c3b659fe-b9bc-44d6-8cf1-5ad18c4fa346
# ╠═84df1941-72a6-4cf5-9e6e-7f79d2534dd6
# ╠═8b601177-be58-4f32-a32f-2bc675c6364a
# ╠═bd6a4dee-1055-4c14-bafd-74b64daa0715
# ╠═d0194709-f4f5-4f7a-acba-76b600dcff92
# ╠═72dcf6e7-8fb8-47f7-bfaa-f5529210e8d0
# ╟─65c9df3c-ea01-48fa-b192-52b3f3b76f9b
# ╟─1f353ddd-4e9a-4437-8d2e-97a4c466ee28
# ╠═10899fba-616b-4bfc-9d75-49383512e747
# ╠═a228cbf8-4d79-41a4-b419-afe26b981a1c
# ╟─334e9473-fc4f-466c-a046-6c66bf65b738
# ╠═a2c59aa5-6c86-4d4a-aa84-c6675446083d
# ╟─aa851337-3dd3-4521-a402-52a429bec1e6
# ╠═e853a3d8-0a18-43b7-949b-c374203d4838
# ╠═0373a338-9a93-42b8-98ee-9768a50dda67
# ╠═03e64fbb-49b0-4b68-ae1a-c5d155651146
# ╠═268ad486-9907-4aaa-a271-a77664fea675
# ╠═405f2f7e-123a-4ccf-9fc1-9b6eee4dae12
# ╠═19c04849-68a6-4ea3-8a1b-e9183c9b852e
# ╠═b5686627-bfa4-4802-9122-eed55fb7db96
# ╠═05cf3e18-950e-427c-91d9-be289da5c5c5
# ╠═385df0f4-bd29-4702-a699-f70ea0093d06
# ╟─6fa0cf17-0c5a-4285-ab2c-09c55dd2e31b
# ╠═f6a0b43f-a56a-41d2-a89f-c1f17f0c6154
# ╠═9b39f9e2-2080-4ab7-ba84-7e58b0047805
# ╠═b9ac10ef-1ea5-46ea-8f9a-471b2a25380d
# ╟─5a09303d-5360-48c9-9554-c769e63b1ee0
# ╠═20a53747-25a1-4642-9043-64109ff28711
# ╟─5f0f4233-ae8b-432f-93e9-f2b3a54dc890
# ╠═bd5aa941-8ab9-4b43-b6f5-f2320f1f9e03
# ╟─1d87f244-f6c5-4e12-9d0d-5c2e7c4805e5
# ╠═8f9e3da4-baab-431d-a4c7-67e753144aef
# ╠═31708cf8-1d56-4b16-ad4e-df356a8924fc
# ╠═cafed3ea-9f5a-4549-bc66-1181c1e7c6c8
# ╠═2269beda-ae30-49ca-af6a-65a681dc7405
# ╟─9d8e3d8d-a4e6-418c-a241-b4b1cd119a14
# ╠═8f125f2e-ec57-40f0-96d6-e04a544db0c3
# ╟─ad1c418f-a4b7-45fa-8798-534b5e474b99
# ╠═f72f5020-6e3f-4b47-9340-55e5d687b56c
# ╟─d3ce01ea-a5d6-41e3-8f98-e0042026da9b
# ╠═2e85a77c-7bb0-47ba-949e-4fdccea6b8e8
# ╠═db774a1b-dbda-42a3-abb2-89ae26a2a31c
# ╠═a9ec77d3-27bc-4a0c-b74a-555b00e6f29a
# ╠═8b27b710-b790-441d-8541-17fae3a72b32
# ╟─341e35e2-f280-47c6-9da5-0d2dc137bce0
# ╟─2b8e12f7-8504-4c28-8208-fd9afb52d1ae
# ╠═12240afd-a739-4c51-8570-23770faf48fd
# ╟─89eab716-678e-407e-acc5-3bd4a490aaf5
# ╠═e12840cf-d45e-4f66-aba9-c178d46f744b
# ╟─ea52b48a-5b66-4d00-b12e-a41ed8a4bde1
# ╠═be4f9b9c-1f42-414e-b74b-f2df17811e8e
# ╠═84957e23-4d9c-409f-b1cf-ab0c8db1a1b6
# ╠═4c3d0757-f378-46ca-944a-cc19d56dd11b
# ╠═33b85892-7ce9-4392-a3dc-b0a69592dba9
# ╟─3d5aa5e2-cfd9-4c5e-83fc-d855ea5c2f30
# ╠═602601ca-081a-4abc-8726-6007be5558cf
# ╟─ac3acca8-868d-47a4-a14a-9696876fdc78
# ╠═28d97f1f-c865-4f15-a7a4-0224039a6d83
# ╟─b7d25c42-dc50-406c-93d6-fd5aec863c10
# ╠═80c2adfa-f6b6-414e-bf68-74f807b992ff
# ╠═830871a6-3f35-4c97-8525-7967b0fbb427
# ╠═3aaac942-afae-4ffc-b1ea-6c023d6e932b
# ╠═0011a165-5fb0-4964-8065-de0663249ee2
# ╟─d179a807-dcf1-4524-81e4-2071c4933f13
# ╠═4a769949-2e00-4825-b852-d419149b47fd
# ╟─81bcc301-c5a9-4ac6-92d9-0f7dcd8676b8
# ╠═c13b7546-2094-42fc-bc78-86fe8e6e5a59
# ╟─30fc4df4-6e6f-42ac-bc4e-222efb20291d
# ╠═4f4bc464-c398-472f-861a-dff6cf96388d
# ╠═b4aab606-8d2e-4110-a0a6-cca477375a45
# ╠═b02751b0-c28c-4d1b-976b-b15617ddc1f1
# ╠═e5eb2459-faad-43d9-a702-cae081f5ae8f
# ╟─aeb23eb0-6ac0-4791-95ce-563bbd9f551b
# ╟─a39be3e4-cd38-4bea-ac3e-4d29c847a5e7
# ╠═38aa73f8-9ec2-43b8-983f-688c7caefa15
# ╟─b1dda8e1-2e73-4e09-a73a-e151a066b99b
# ╠═2cbc8541-9c85-473e-af24-a723d609b590
# ╟─6aaa21ae-afeb-4099-9bb6-096ed8e9c093
# ╠═9c5f120e-f9bb-4737-8031-95da4ea26e23
# ╠═8a409f2e-7ab8-47b3-a855-fb1a84c5881a
# ╠═ce3b04f9-1461-4007-825f-7481f5da9312
# ╠═6cfe8bf9-7869-41e6-9df2-d2c18a159255
# ╟─26609afb-14a5-4a5f-a90e-55fa9f3b17dd
# ╠═ff7b8b70-e0ca-4c83-aa3a-ed6a5b60fc55
# ╟─60c507e8-e885-491e-a8b3-db61d61ee4fa
# ╠═513dfbb8-9c34-40f4-93c9-8ae21d5dffb6
# ╟─bf1a2681-60f3-425f-b938-98934aba4588
# ╠═435321e8-c027-496d-8280-f9aad31f3632
# ╠═673d4151-ec74-434e-8d0c-bfe162f4b6a0
# ╠═2ba692fb-7a93-4521-b1a8-d334f8a8d763
# ╠═78ce4c60-51c7-4c64-a3ce-d71d1529a8c7
# ╟─4238b68f-b4fc-40b3-b859-3a5130834017
# ╠═bacc56de-f082-436e-bffa-23379fdb9c61
# ╟─3e321593-3949-41c6-a734-cfa05fcc471e
# ╠═58cc58d7-4706-41c4-a2f7-04d9cd6eb31d
# ╟─5d3a539f-9ff3-4cbb-8957-88035698f69f
# ╠═55ae3b8f-54cc-4d4a-a3c5-37a1afad2c0b
# ╠═bd53b5f1-73a2-4f65-aa20-a7cdba456ce0
# ╠═9fafb79b-0706-470d-8a24-f2682621e39c
# ╠═c2116577-2f05-40dc-9094-92c621955f93
# ╟─258caf22-6f46-4341-80eb-d572aa26ca59
# ╠═3d573812-8bc2-4dce-99fb-6883b86433f4
# ╠═52f3d56d-9219-43f3-ba32-925faa2e0e07
# ╟─260b441d-12eb-44d5-adb2-dc5401abdb20
# ╠═c0b060b6-ef4d-415e-9f8c-a39ab4373541
# ╠═3f79ca52-b120-43af-958e-6a0e5edf079f
# ╠═0de07ba8-f1d6-4c01-86fd-82e42e9a29bb
# ╟─40f4a8e5-8a97-43d6-8d83-4993d5f59b83
# ╠═e89bd455-71ba-41c0-870c-2dad225e129e
# ╠═36130f3d-8a03-49ce-ab36-f9ff0b452c8d
# ╠═e53a38e7-ba34-45af-bc39-86026763241f
# ╟─bdf7c9e5-765d-40b2-b4ae-19e068348639
# ╟─59a17ced-27c0-415f-a95a-7a8af8aba7e3
# ╟─adacf893-01ef-4cd3-bcc5-66b27e77e3be
# ╟─b24231d2-aa2e-45ad-80f4-e875fa282cd0
# ╠═14b79f5f-87c5-4db1-855b-607cec0946b0
# ╟─5757d60a-6944-4286-ad55-95ac78acd903
# ╠═d10d9943-948e-4d91-881d-93ac9c8fdf5b
# ╠═5c3f138e-d46e-4a3d-9bc4-e2ade6784907
# ╠═c8066892-8ea6-4383-8771-699eb1394a5c
# ╠═ee2bd00c-e949-48f9-8e3c-1e05eb749075
# ╠═455c1094-548f-4983-918e-2a7eaa47ff00
# ╠═acf1c58f-24c3-4b6a-a003-67774960112e
# ╠═5c94eae2-9e45-49c9-a123-fbd46932c95f
# ╠═7ee1efe9-d6ac-4076-a5f0-8e18b5d64702
# ╟─e9d6989d-93fb-43ee-b18a-582dd637eed4
# ╟─9e533473-6fd4-4155-9c35-d6901580b45f
# ╟─27be579e-87fc-4b77-833c-01c3361375af
# ╟─ade95452-1c8d-4888-9313-b092409ca609
# ╠═844754c6-a0ae-4cc7-9d6d-38204d7b24f5
# ╠═c7b49892-0079-4c79-92d7-3717460efd0f
# ╟─49275899-e5ec-47f4-b4f7-43401529642f
# ╠═6471de28-e427-4152-b778-74c3d208b5b1
# ╠═b99951e3-6292-42e4-801e-4824bca11720
# ╠═92be24c9-e0a5-4539-9501-93a8d6c52f10
# ╠═1723c48d-341f-4898-9b57-4202664a9415
# ╠═4d1825ea-b5cd-4445-99d5-a218845e629c
# ╠═3912fb58-7431-40f7-9823-9a4c91487bf6

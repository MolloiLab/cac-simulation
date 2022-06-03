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

# ╔═╡ 45f704d4-66e5-49db-aef5-05132f3853ee
# ╠═╡ show_logs = false
begin
	let
		using Pkg
		Pkg.activate(mktempdir())
		Pkg.Registry.update()
		Pkg.add("PlutoUI")
		Pkg.add("CairoMakie")
		Pkg.add("Statistics")
		Pkg.add("StatsBase")
		Pkg.add("ImageMorphology")
		Pkg.add("ImageFiltering")
		Pkg.add("CSV")
		Pkg.add("DataFrames")
		Pkg.add("GLM")
		Pkg.add("Noise")
		Pkg.add(url="https://github.com/JuliaHealth/DICOM.jl")
		Pkg.add(url="https://github.com/Dale-Black/DICOMUtils.jl")
		Pkg.add(url="https://github.com/Dale-Black/PhantomSegmentation.jl")
		Pkg.add(url="https://github.com/Dale-Black/CalciumScoring.jl")
		Pkg.add("ImageComponentAnalysis")
		Pkg.add("Tables")
	end
	
	using PlutoUI
	using CairoMakie
	using Statistics
	using StatsBase: quantile!
	using ImageMorphology
	using ImageFiltering
	using CSV
	using DataFrames
	using GLM
	using Noise
	using DICOM
	using DICOMUtils
	using PhantomSegmentation
	using CalciumScoring
	using ImageComponentAnalysis
	using Tables
end

# ╔═╡ e7c4aaab-83ef-4256-b392-f8ef7a899a05
TableOfContents()

# ╔═╡ 48e3097f-0767-4dad-9d7c-94e0899f790b
md"""
## Load DICOMS

All you need to do is set `base_path` once and leave it. After that, the only thing that should change is the `VENDER`, once for every set, and the `SCAN_NUMBER`, once for each scan.
"""

# ╔═╡ bc9383c0-e477-4e1a-a2fa-7f5c1d29f103
begin
	SCAN_NUMBER = 1
	VENDER = "135"
	# SIZE = "small"
	SIZE = "medium"
	# SIZE = "large"
	# DENSITY = "low"
	DENSITY = "normal"
	TYPE = "integrated_scoring"
	BASE_PATH = string("/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/images_new/", SIZE, "/", DENSITY, "/")
end

# ╔═╡ 1df1a491-6207-4b5c-a5e8-15a07323b6e7
md"""
**Scans**
- normal (consistent underestimation of HD inserts; probably solved by using all inserts for calibration line instead of calibration rod)
  - small
    - 80 (thresh: 130, ✅)
    - 100 (thresh: 130, ✅)
    - 120 (thresh: 130, ✅)
    - 135 (thresh: 130, ✅)
  - medium
    - 80 (thresh: 130, ✅)
    - 100 (thresh: 130, ✅)
    - 120 (thresh: 130, ✅)
    - 135 (thresh: 130, ✅)
  - large
    - 80 (thresh: 130, ✅)
    - 100 (thresh: 130, ✅)
    - 120 (thresh: 130, ✅)
    - 135 (thresh: 130, ✅)
- low
  - small
    - 80 (thresh: 60, ✅)
    - 100 (thresh: 60, ✅)
    - 120 (thresh: 60, ✅)
    - 135 (thresh: 60, ✅)
  - medium
    - 80 (thresh: 75, ✅)
    - 100 (thresh: 75, ✅)
    - 120 (thresh: 75, ✅)
    - 135 (thresh: 70, ✅)
  - large
    - 80
    - 100
    - 120
    - 135
"""

# ╔═╡ 6aa51429-981a-4dea-a0f6-2935867d5b2a
md"""
**Everything below should be automatic, just scroll through to visually inspect that things make sense**
"""

# ╔═╡ e849cf69-65c7-4e5e-9688-fc249d471f2c
root_path = string(BASE_PATH, VENDER)

# ╔═╡ 3e649fee-a2e9-4f0f-b705-3966f69d97ea
dcm_path_list = dcm_list_builder(root_path)

# ╔═╡ a327ef18-2941-4783-9266-7332826eaf58
pth = dcm_path_list[SCAN_NUMBER]

# ╔═╡ b44d18f8-1b86-4235-bb5b-7a77a1af55e0
scan = basename(pth)

# ╔═╡ b7e0b678-0627-44fe-b0cb-3ef2bccae6a7
header, dcm_array, slice_thick_ori1 = dcm_reader(pth);

# ╔═╡ 66f5fc66-d5ba-45ba-8624-55adc58085e4
md"""
## Helper Functions
"""

# ╔═╡ 6e812172-6371-4461-9365-22f68ef16e53
function collect_tuple(tuple_array)
	row_num = size(tuple_array)
	col_num = length(tuple_array[1])
	container = zeros(Int64, row_num..., col_num)
	for i in 1:length(tuple_array)
		container[i,:] = collect(tuple_array[i])
	end
	return container
end

# ╔═╡ 483a14dd-e798-41ed-9144-13678f8b8461
function overlay_mask_bind(mask)
	indices = findall(x -> x == 1, mask)
	indices = Tuple.(indices)
	label_array = collect_tuple(indices)
	zs = unique(label_array[:,3])
	return PlutoUI.Slider(1:length(zs), default=3, show_value=true)
end

# ╔═╡ 2cfdeb06-71a2-411b-accb-2b4a3f56b477
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

# ╔═╡ b3bae2ad-16ed-4381-b6a4-448cb5f0c6c1
md"""
## Segment Heart
"""

# ╔═╡ 4ff19af0-7f53-4f24-b315-08bd54a488e3
masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3)÷2);

# ╔═╡ b74941c4-12d4-4be1-81f7-ef1f4d288984
@bind a PlutoUI.Slider(1:size(masked_array, 3), default=10, show_value=true)

# ╔═╡ 84e56c91-14a4-45b7-81a1-77e778bae695
heatmap(transpose(masked_array[:, :, a]), colormap=:grays)

# ╔═╡ abfd965c-df00-4028-a0f3-8356a261b527
begin
	fig = Figure()
	
	ax = Makie.Axis(fig[1, 1])
	ax.title = "Raw DICOM Array"
	heatmap!(transpose(dcm_array[:, :, 5]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2], center_insert[1]:center_insert[1], markersize=10, color=:red)
	fig
end

# ╔═╡ 41541a09-609b-4f46-9a25-eb974422efc0
begin
	fig2 = Figure()
	
	ax2 = Makie.Axis(fig2[1, 1])
	ax2.title = "Mask Array"
	heatmap!(transpose(mask), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2], center_insert[1]:center_insert[1], markersize=10, color=:red)
	fig2
end

# ╔═╡ e6753c9a-d172-47fe-b066-78cb43df2c89
begin
	fig3 = Figure()
	
	ax3 = Makie.Axis(fig3[1, 1])
	ax3.title = "Masked DICOM Array"
	heatmap!(transpose(masked_array[:, :, 1]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2], center_insert[1]:center_insert[1], markersize=10, color=:red)
	fig3
end

# ╔═╡ f8a4eacd-2f69-4755-b017-5a0b0dda1004
md"""
## Segment Calcium Rod
"""

# ╔═╡ ac9d2652-6d69-4cf3-acbf-2e141fb633f0
begin
	global thresh
	if DENSITY == "low"
		thresh = 60
	elseif DENSITY ==  "normal"
		thresh = 130
	end
end

# ╔═╡ dd399c0e-f8e4-464c-8964-bdc0dd657202
calcium_image, slice_CCI, quality_slice, cal_rod_slice = mask_rod(masked_array, header; calcium_threshold=55);

# ╔═╡ 417d150e-0df9-4963-b59e-ebd9acf6d4a0
@bind c PlutoUI.Slider(1:size(calcium_image, 3), default=5, show_value=true)

# ╔═╡ 5f10075b-4d95-4d53-a22f-016749fb7583
heatmap(transpose(calcium_image[:, :, c]), colormap=:grays)

# ╔═╡ b921dcf8-54ea-420b-9285-23e38d5ce433
md"""
## Load (segment) Calcium Inserts
"""

# ╔═╡ ece7d76b-93c8-430c-8679-07cf92585949
# mask_L_HD, mask_M_HD, mask_S_HD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_LD, mask_M_LD, mask_S_LD = mask_inserts_simulation(
#             dcm_array, masked_array, header, slice_CCI, center_insert; calcium_threshold=thresh
# );

# ╔═╡ ca408955-54ff-40ff-9a0f-77ee2666d7a2
# begin
# 	root = string("/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/julia_arrays/", SIZE, "/") 
# 	CSV.write(string(root, "mask_L_HD.csv"),  Tables.table(mask_L_HD), writeheader=false)
# 	CSV.write(string(root, "mask_M_HD.csv"),  Tables.table(mask_M_HD), writeheader=false)
# 	CSV.write(string(root, "mask_S_HD.csv"),  Tables.table(mask_S_HD), writeheader=false)
# 	CSV.write(string(root, "mask_L_MD.csv"),  Tables.table(mask_L_MD), writeheader=false)
# 	CSV.write(string(root, "mask_M_MD.csv"),  Tables.table(mask_M_MD), writeheader=false)
# 	CSV.write(string(root, "mask_S_MD.csv"),  Tables.table(mask_S_MD), writeheader=false)
# 	CSV.write(string(root, "mask_L_LD.csv"),  Tables.table(mask_L_LD), writeheader=false)
# 	CSV.write(string(root, "mask_M_LD.csv"),  Tables.table(mask_M_LD), writeheader=false)
# 	CSV.write(string(root, "mask_S_LD.csv"),  Tables.table(mask_S_LD), writeheader=false)
# end

# ╔═╡ 05ed60db-b4d8-4cdc-9a54-b108ade22557
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

# ╔═╡ 8bb2d451-b575-40bf-b31e-51e80392ef51
masks = mask_L_HD + mask_M_HD + mask_S_HD + mask_L_MD + mask_M_MD + mask_S_MD + mask_L_LD + mask_M_LD + mask_S_LD;

# ╔═╡ 4862496e-464d-40c3-baef-e5c507028d66
heatmap(transpose(masks), colormap=:grays)

# ╔═╡ 997952cd-dc01-4528-a30a-fdc11a46dad8
md"""
## Calibration Prep
"""

# ╔═╡ 144d3db4-343b-44eb-bd13-160013b7a58c
array_filtered = abs.(mapwindow(median, calcium_image[:, :, 2], (3, 3)));

# ╔═╡ 26b4ff9f-1b55-4fcc-99c8-3e1079d8df10
bool_arr = array_filtered .> 0;

# ╔═╡ c3a87ddd-64ad-4aeb-95e9-7f7b1ba1bd28
bool_arr_erode = ((erode(erode(erode(bool_arr)))));

# ╔═╡ da3c96b9-efdb-40f4-a959-747c097e16b2
heatmap(bool_arr, colormap=:grays)

# ╔═╡ 5fec186d-82bb-4b37-b0c6-c3781f2a166d
heatmap(bool_arr_erode, colormap=:grays)

# ╔═╡ f221d7ac-7b27-4dc1-8081-af55d9e9b925
c_img = calcium_image[:, :, 1:3];

# ╔═╡ 0ce913ad-4844-4df0-92f2-79af4d85113e
begin
	mask_cal_3D = Array{Bool}(undef, size(c_img))
	for z in 1:size(c_img, 3)
		mask_cal_3D[:, :, z] = bool_arr_erode
	end
end;

# ╔═╡ c1d1b644-db3b-45b6-9050-342eca241a05
hist(c_img[mask_cal_3D])

# ╔═╡ d302f987-19aa-4401-9b79-20ca04103050
cal_insert_mean = mean(c_img[mask_cal_3D])

# ╔═╡ 37e25efb-abf5-4568-a4d4-ff5c8aff1f64
md"""
### Calibration Line
"""

# ╔═╡ 6c3f0dc5-9db2-4db8-abb0-64e4cda6d6fe
begin
	global density_array
	if DENSITY == "low"
		density_array = [0, 25, 50, 100]
	elseif DENSITY == "normal"
		density_array = [0, 200, 400, 800]
	end
end

# ╔═╡ 6ceb2a0a-8fcf-4f8c-8167-5b2ab770df48
density_array_cal = [0, 200]

# ╔═╡ 43cf1721-57db-4048-8528-b0430c5e8043
intensity_array_cal = [0, cal_insert_mean]

# ╔═╡ 444ec4ef-6140-4556-9305-599d819e58ca
# df_cal = DataFrame(:density => [density_array[1], density_array[4]], :intensity => [intensity_array[1], intensity_array[4]])

# ╔═╡ e50baf62-0a02-47ee-8440-551f68baee0d
md"""
We can see from above that the linear regression returns a best fit line with the formula:

```math
y = mx + b
```

Which can be solved for ``x`` and then used to calculate the density (``x``) given some measured intensity (``y``)

```math
x = \frac{y - b}{m}
```

"""

# ╔═╡ 3b30e707-02b0-49c2-a0b3-25b7fcbf79ee
md"""
# Score Large Inserts
"""

# ╔═╡ d27df860-7ca8-4e93-90e4-a58605aa1aaa
arr = masked_array[:, :, 4:6];

# ╔═╡ 5dfb50b9-7545-4948-8026-2074ce5d86cf
single_arr = masked_array[:, :, slice_CCI];

# ╔═╡ 177ab471-2783-4e83-95fd-fcd8d02cb104
md"""
## Background
"""

# ╔═╡ aa02f7d8-6b01-4da9-b0f5-d1df587c7d4d
md"""
#### Dilated background
"""

# ╔═╡ dec0f541-4f89-45e3-b084-ecdbc837b794
begin
	background_mask = zeros(size(arr)...)
	background_mask[center_insert[1]-5:center_insert[1]+5, center_insert[2]-5:center_insert[2]+5, 2] .= 1
	background_mask = Bool.(background_mask)
	background_mask_dil = dilate(dilate(background_mask))
end;

# ╔═╡ 94c66847-765d-4d58-8510-043285d6cdcd
md"""
#### Ring background
"""

# ╔═╡ 703c7556-3904-41d7-972b-48ea31890a7d
begin
	ring_background = background_mask_dil - background_mask
	single_ring_bkg = Bool.(ring_background[:, :, 2])
	s_bkg = mean(single_arr[single_ring_bkg])
end

# ╔═╡ 983fa5b0-45fe-4ba9-8eea-ad5c030ed32f
md"""
#### Noise
"""

# ╔═╡ 2a6961fc-f94f-4250-aeb8-63bbfd9e29a0
md"""
## High Density
"""

# ╔═╡ db4992af-e78b-4f84-b0b3-c13b6cd6cde5
begin
	mask_L_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_HD_3D[:, :, z] = mask_L_HD
	end
end;

# ╔═╡ 03427ccf-69b9-4a1c-b988-30ad9f0a214c
begin
	eroded_mask_L_HD = erode(erode(mask_L_HD_3D))
	high_density_cal = mean(arr[eroded_mask_L_HD])
end

# ╔═╡ 58172907-2026-4caa-a6ad-ab484b86b329
md"""
#### Dilated mask
"""

# ╔═╡ f8e646a2-8b1e-47f6-bc5d-760941d542e8
dilated_mask_L_HD = dilate(dilate(mask_L_HD_3D));

# ╔═╡ 080422ff-77a2-4b86-8a0f-c2bea1f0d722
@bind g2 overlay_mask_bind(dilated_mask_L_HD)

# ╔═╡ 1cba2d1a-19c2-49b7-8757-5843bb16a0f1
overlay_mask_plot(arr, dilated_mask_L_HD, g2, "dilated mask")

# ╔═╡ cf4bc01d-5603-4d47-8863-e413be42f223
md"""
#### Ring (background) mask
"""

# ╔═╡ ccf2c3d2-7243-412f-bbe1-84b4ac3f6272
ring_mask_L_HD = dilate(dilate(dilate(dilate(mask_L_HD_3D)))) - dilate(dilate(dilate(mask_L_HD_3D)));

# ╔═╡ b990d904-e359-45fe-990b-87ab625ba1b0
@bind g4 overlay_mask_bind(ring_mask_L_HD)

# ╔═╡ 0af76504-ffc9-44ac-aecd-a3b956166273
overlay_mask_plot(arr, ring_mask_L_HD, g4, "ring mask")

# ╔═╡ dc94ad4a-0692-4e0e-84e5-f706ab1ca316
begin
	single_ring_mask_L_HD = Bool.(ring_mask_L_HD[:, :, 3])
	s_bkg_L_HD = mean(single_arr[single_ring_mask_L_HD])
end

# ╔═╡ 28c869df-75eb-42dc-9033-867ca4782f41
pixel_size = DICOMUtils.get_pixel_size(header)

# ╔═╡ e93b26b5-ba84-4410-87b7-eefad3d8ab3e
md"""
## Medium Density
"""

# ╔═╡ 1b97a53d-6ed8-4414-befa-146225f0f95a
begin
	mask_L_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_MD_3D[:, :, z] = mask_L_MD
	end
end;

# ╔═╡ fa126b46-4dbb-4b8e-a3a7-7cf28a9ed88d
begin
	eroded_mask_L_MD = erode(erode(mask_L_MD_3D))
	med_density_cal = mean(arr[eroded_mask_L_MD])
end

# ╔═╡ 336cf5d0-ada0-4f9a-9788-bfd0ab72087a
md"""
#### Dilated mask
"""

# ╔═╡ 94779096-626f-4cb6-9937-0bc4bc8ab45b
dilated_mask_L_MD = dilate(dilate(mask_L_MD_3D));

# ╔═╡ 937a5188-e51e-4017-be9a-2d8051e79fd1
@bind h2 overlay_mask_bind(dilated_mask_L_MD)

# ╔═╡ f4ee897f-d6ad-4dfd-a91f-a759ecd170e5
overlay_mask_plot(arr, dilated_mask_L_MD, h2, "dilated mask")

# ╔═╡ afbc96dc-f5fa-4305-bd3d-119ef8d416a8
md"""
#### Ring (background) mask
"""

# ╔═╡ 078b6320-34d9-4318-863e-a2a6dcdeb500
ring_mask_L_MD = dilate(dilate(dilate(dilate(mask_L_MD_3D)))) - dilate(dilate(dilate(mask_L_MD_3D)));

# ╔═╡ 8e310e91-c557-48e5-a5d1-042c1ac8f7b3
@bind h4 overlay_mask_bind(ring_mask_L_MD)

# ╔═╡ ed5c5da5-21dc-4052-bc18-6b209edeabfd
overlay_mask_plot(arr, ring_mask_L_MD, h4, "ring mask")

# ╔═╡ 40c7c265-f558-4f74-8d96-40cd83c2d5bb
begin
	single_ring_mask_L_MD = Bool.(ring_mask_L_MD[:, :, 3])
	s_bkg_L_MD = mean(single_arr[single_ring_mask_L_MD])
end

# ╔═╡ 73105b5a-05b7-429f-9a7f-b8c12684df91
md"""
## Low Density
"""

# ╔═╡ 31e9a969-6514-404c-89f6-0f038cc115f2
begin
	mask_L_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_LD_3D[:, :, z] = mask_L_LD
	end
end;

# ╔═╡ d279e6b4-fba0-4fcb-a1ef-6c2c33ba8f63
begin
	eroded_mask_L_LD = erode(erode(mask_L_LD_3D))
	low_density_cal = mean(arr[eroded_mask_L_LD])
end

# ╔═╡ 2f06a3e2-5c65-4f7d-af65-829a8e657a92
intensity_array = [0, low_density_cal, med_density_cal, high_density_cal] # HU

# ╔═╡ b6d389f3-6ecf-4798-978c-c28ef4b8e5da
df_cal = DataFrame(:density => density_array, :intensity => intensity_array)

# ╔═╡ 6c3c7d22-78e6-428e-ac62-90a473b06e2d
linearRegressor = lm(@formula(intensity ~ density), df_cal);

# ╔═╡ f42a62d8-11e9-48b8-b6bc-ae25b6797707
linearFit = predict(linearRegressor)

# ╔═╡ b407d628-0db0-4d13-9a6c-8e4ff5a84cbd
m = linearRegressor.model.pp.beta0[2]

# ╔═╡ dbec08a2-96c5-488b-89b9-5df6fe1f5a12
b = linearRegressor.model.rr.mu[1]

# ╔═╡ 3fa80bee-20f6-4e9c-a14d-75a77482a2c3
density(intensity) = (intensity - b) / m

# ╔═╡ 04c77c9d-c498-4e09-abff-9c361e887c56
intensity(ρ) = m*ρ + b

# ╔═╡ 3895fafc-a999-4b06-974d-caef190f7830
begin
	single_bkg_center = Bool.(background_mask[:, :, 2])
	# S_Obj_bkg = mean(single_arr[single_bkg_center])
	S_Obj_bkg = intensity(200)
end

# ╔═╡ a7de43e0-795d-4267-b7dd-3be85a201485
begin
	alg_bkg = Integrated(arr[background_mask])
	ρ_bkg = 0.2 # mg/mm^3
	mass_bkg = score(s_bkg, S_Obj_bkg, pixel_size, ρ_bkg, alg_bkg)
end

# ╔═╡ fd089ad4-4e08-4e1c-b049-3378de42b2df
S_Obj_HD = intensity(800)

# ╔═╡ 7d61e536-fac5-43e1-9ce9-5952d7d20e8d
begin
	alg_L_HD = Integrated(arr[mask_L_HD_3D])
	ρ_hd = 0.8 # mg/mm^3
	mass_l_hd = score(s_bkg_L_HD, S_Obj_HD, pixel_size, ρ_hd, alg_L_HD)
end

# ╔═╡ 0bf239a7-233c-42f2-b7eb-29336b472ff9
S_Obj_MD = intensity(400)

# ╔═╡ 367750d8-a287-4313-b6df-d6533849623f
begin
	alg_L_MD = Integrated(arr[mask_L_MD_3D])
	ρ_md = 0.4
	mass_l_md = score(s_bkg_L_MD, S_Obj_MD, pixel_size, ρ_md, alg_L_MD)
end

# ╔═╡ b03c711f-9e0e-4089-add1-5abbde81cce1
begin
	f = Figure()
	ax1 = Axis(f[1, 1])
	
	scatter!(density_array, intensity_array)
	lines!(density_array, linearFit, color = :red)
	ax1.title = "Calibration Line (Intensity vs Density)"
	ax1.ylabel = "Intensity (HU)"
	ax1.xlabel = "Density (mg/cm^3)"
	
	f
end

# ╔═╡ d8d7a93f-85e7-43e7-a23d-11e7c6ce27be
md"""
#### Dilated mask
"""

# ╔═╡ 466cfc6f-2915-44b4-9231-4d0047accd6d
dilated_mask_L_LD = dilate(dilate(mask_L_LD_3D));

# ╔═╡ b849a64f-965a-40b7-b2f5-1a27fd8ceeb5
@bind i2 overlay_mask_bind(dilated_mask_L_LD)

# ╔═╡ 7f52b670-611d-4bd3-bc80-77ba59aa2d70
overlay_mask_plot(arr, dilated_mask_L_LD, i2, "dilated mask")

# ╔═╡ 4c425cf3-2ece-4df1-b3ab-4058c7d7aa13
md"""
#### Ring (background) mask
"""

# ╔═╡ 52b9efbd-80b3-4205-b4a7-b18809d6ba24
ring_mask_L_LD = dilate(dilate(dilate(dilate(mask_L_LD_3D)))) - dilate(dilate(dilate(mask_L_LD_3D)));

# ╔═╡ 7a618b0c-9429-472a-b076-774b9d229cab
@bind i4 overlay_mask_bind(ring_mask_L_LD)

# ╔═╡ 544151a5-3b47-4571-a6f0-9800d5b11085
overlay_mask_plot(arr, ring_mask_L_LD, i4, "ring mask")

# ╔═╡ 9da779a0-86bd-450d-9401-7d3507aaec7d
begin	
	single_ring_mask_L_LD = Bool.(ring_mask_L_LD[:, :, 3])
	s_bkg_L_LD = mean(single_arr[single_ring_mask_L_LD])
end

# ╔═╡ 187c110a-21a9-4b31-82b6-d9c239dbd899
S_Obj_LD = intensity(200)

# ╔═╡ 7560dc46-fa36-4706-925c-0675397c5922
begin
	alg_L_LD = Integrated(arr[mask_L_LD_3D])
	ρ_LD = 0.2
	mass_l_ld = score(s_bkg_L_LD, cal_insert_mean, pixel_size, ρ_LD, alg_L_LD)
end

# ╔═╡ ee7092e6-52fb-4a6a-ad21-90cd2b46c7fa
md"""
# Score Medium Inserts
"""

# ╔═╡ 62faf81f-eb3b-431a-8e2b-239735db3257
md"""
## High Density
"""

# ╔═╡ 6fba5f51-7940-45fe-8a0b-30ea12962a0c
begin
	mask_M_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_HD_3D[:, :, z] = mask_M_HD
	end
end;

# ╔═╡ 10765555-8507-413f-a582-1f47619f0181
md"""
#### Dilated mask
"""

# ╔═╡ f8519d1a-c5e7-41a3-90d3-5f6200f767cf
dilated_mask_M_HD = dilate(dilate(dilate(dilate(mask_M_HD_3D))));

# ╔═╡ 5aff6e7a-b635-44d7-a91e-1a30c5f4fd17
@bind j2 overlay_mask_bind(dilated_mask_M_HD)

# ╔═╡ 4aabe7d9-e730-4b4f-9dd2-e51e922ac08c
overlay_mask_plot(arr, dilated_mask_M_HD, j2, "dilated mask")

# ╔═╡ 706c30df-ce63-4c81-9c99-227d56ab2991
md"""
#### Ring (background) mask
"""

# ╔═╡ 6598b28c-5850-44bc-8bd0-995ff3d59f09
ring_mask_M_HD = dilate(dilate(dilate(dilate(dilate(mask_M_HD_3D))))) - dilate(dilate(dilate(dilate(mask_M_HD_3D))));

# ╔═╡ 47615958-a8fc-4fc2-852f-e11f69041a86
@bind j4 overlay_mask_bind(ring_mask_M_HD)

# ╔═╡ 85097de3-d336-4c85-a759-ad1b3e4083cc
overlay_mask_plot(arr, ring_mask_M_HD, j4, "ring mask")

# ╔═╡ 0f02e277-06b7-46f7-b6e8-0d5b6516dbe2
begin
	single_ring_mask_M_HD = Bool.(ring_mask_M_HD[:, :, 3])
	s_bkg_M_HD = mean(single_arr[single_ring_mask_M_HD])
end

# ╔═╡ bd4cec3c-446f-48e4-b06b-fbde4848daaf
begin
	alg_M_HD = Integrated(arr[mask_M_HD_3D])
	mass_m_hd = score(s_bkg_M_HD, S_Obj_HD, pixel_size, ρ_hd, alg_M_HD)
end

# ╔═╡ 0f71036e-b118-430d-bce7-420b2591f1a2
md"""
## Medium Density
"""

# ╔═╡ 539d9e60-1f76-47c5-ab4f-cb520c520965
begin
	mask_M_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_MD_3D[:, :, z] = mask_M_MD
	end
end;

# ╔═╡ c99672f4-9ee7-4d42-8a15-eee68be2a8da
md"""
#### Dilated mask
"""

# ╔═╡ c34e32c6-4f28-435b-b124-7d3727f107a7
dilated_mask_M_MD = dilate(dilate(dilate(dilate(mask_M_MD_3D))));

# ╔═╡ b1c10feb-3148-4de5-b305-b0721894b76e
@bind k2 overlay_mask_bind(dilated_mask_M_MD)

# ╔═╡ bd91bd0d-f368-4340-916a-ea197969c0e4
overlay_mask_plot(arr, dilated_mask_M_MD, k2, "dilated mask")

# ╔═╡ b5cb8919-df84-48a6-b00b-cf7dba301c76
md"""
#### Ring (background) mask
"""

# ╔═╡ 48eb8f20-1d00-4434-8dcd-4960b2392d28
ring_mask_M_MD = dilate(dilate(dilate(dilate(dilate(dilate(mask_M_MD_3D)))))) - dilate(dilate(dilate(dilate(dilate(mask_M_MD_3D)))));

# ╔═╡ 6c3c8eb5-c651-45e3-a912-8fd0dc330b51
@bind k4 overlay_mask_bind(ring_mask_M_MD)

# ╔═╡ f2bafa39-c1ac-42ca-89c0-df4d219c2c74
overlay_mask_plot(arr, ring_mask_M_MD, k4, "ring mask")

# ╔═╡ 906c678f-2006-4f66-a135-97e25182a3d6
begin
	single_ring_mask_M_MD = Bool.(ring_mask_M_MD[:, :, 3])
	s_bkg_M_MD = mean(single_arr[single_ring_mask_M_MD])
end

# ╔═╡ 43afbb30-b1b5-452b-8c82-25d173f24699
begin
	alg_M_MD = Integrated(arr[mask_M_MD_3D])
	mass_m_md = score(s_bkg_M_MD, S_Obj_MD, pixel_size, ρ_md, alg_M_MD)
end

# ╔═╡ 386ea37d-70ba-4530-be9e-549d5e0ab1ce
md"""
## Low Density
"""

# ╔═╡ 25434519-be6d-4616-9353-a0d342129ea4
begin
	mask_M_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_LD_3D[:, :, z] = mask_M_LD
	end
end;

# ╔═╡ b37af6a8-8f51-4263-ac85-bea25fd2aae2
md"""
#### Dilated mask
"""

# ╔═╡ b9dbf6c2-891f-490d-a89c-bd50891cfdee
dilated_mask_M_LD = dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))));

# ╔═╡ 7b8af8e8-e979-4b57-a65a-ad99a40f08b9
@bind l2 overlay_mask_bind(dilated_mask_M_LD)

# ╔═╡ 493ef94e-d6ac-4ccb-b342-0812792c2536
overlay_mask_plot(arr, dilated_mask_M_LD, l2, "dilated mask")

# ╔═╡ 43f5bf4e-c882-4596-a4b5-6063ce63e17c
md"""
#### Ring (background) mask
"""

# ╔═╡ c44aa6b4-bada-4ab0-a748-71d0e215c1c7
ring_mask_M_LD = dilate(dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))))) - dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))));

# ╔═╡ 48a0f2c4-0bea-4c7d-bf8a-dc3632a490fa
@bind l4 overlay_mask_bind(ring_mask_M_LD)

# ╔═╡ e03ccab0-b786-49ca-b7dd-b2b7f76e04a2
overlay_mask_plot(arr, ring_mask_M_LD, l4, "ring mask")

# ╔═╡ 5133060e-c1fb-492f-93d7-f161e8da73be
begin
	single_ring_mask_M_LD = Bool.(ring_mask_M_LD[:, :, 3])
	s_bkg_M_LD = mean(single_arr[single_ring_mask_M_LD])
end

# ╔═╡ 907b3f24-47b8-4519-b70a-de44ca764bb2
begin
	alg_M_LD = Integrated(arr[mask_M_LD_3D])
	mass_m_ld = score(s_bkg_M_LD, S_Obj_LD, pixel_size, ρ_LD, alg_M_LD)
end

# ╔═╡ c52478ca-2012-4690-bc39-6d32997427f7
md"""
# Score Small Inserts
"""

# ╔═╡ 5de67b93-72d7-4f86-b94e-245f4cfff576
md"""
## High Density
"""

# ╔═╡ 6323a01d-100a-4d09-b3ad-5821e88bed0a
begin
	mask_S_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_HD_3D[:, :, z] = mask_S_HD
	end
end;

# ╔═╡ d0ea5ea6-5c98-44b8-94d4-4104eb1f6080
md"""
#### Dilated mask
"""

# ╔═╡ 594ed641-0f62-461d-b031-2a5de1467d35
dilated_mask_S_HD = dilate(((dilate(dilate((mask_S_HD_3D))))));

# ╔═╡ 018d2392-7133-4da8-b4c9-268af6fa6f3a
@bind m2 overlay_mask_bind(dilated_mask_S_HD)

# ╔═╡ 9d72b647-73d9-45c7-b473-bbdea2817cb5
overlay_mask_plot(arr, dilated_mask_S_HD, m2, "dilated mask")

# ╔═╡ 27ffaff1-12c5-43f1-ac21-84fe25877f24
md"""
#### Ring (background) mask
"""

# ╔═╡ a4ddb30d-011c-4855-96f3-01ba441297f1
ring_mask_S_HD = dilate((dilate((dilate(mask_S_HD_3D))))) - dilate(dilate(((mask_S_HD_3D))));

# ╔═╡ 46110aad-07a8-4a56-9a6d-db71054a5883
@bind m4 overlay_mask_bind(ring_mask_S_HD)

# ╔═╡ 070aedea-d5e9-43f7-a2cf-14fdad760a72
overlay_mask_plot(arr, ring_mask_S_HD, m4, "ring mask")

# ╔═╡ 76b6bf9f-58fb-4a88-a21f-0695d2c76073
begin
	single_ring_mask_S_HD = Bool.(ring_mask_S_HD[:, :, 3])
	s_bkg_S_HD = mean(single_arr[single_ring_mask_S_HD])
end

# ╔═╡ 1590d3c5-d9f1-4c00-9730-639677635056
begin
	alg_S_HD = Integrated(arr[mask_S_HD_3D])
	mass_s_hd = score(s_bkg_S_HD, S_Obj_HD, pixel_size, ρ_hd, alg_S_HD)
	if mass_s_hd < 0
		mass_s_hd = 0
	end
	mass_s_hd
end

# ╔═╡ aa3f6b52-2c1a-462a-a89c-bcde6f0f59ba
md"""
## Medium Density
"""

# ╔═╡ e17a81a4-35a6-4d15-bf14-275114a49608
begin
	mask_S_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_MD_3D[:, :, z] = mask_S_MD
	end
end;

# ╔═╡ 01dce72f-6793-4360-80fc-bc2e97d579de
md"""
#### Dilated mask
"""

# ╔═╡ 5250fe6b-b94a-47b0-9142-525ca34f1d16
dilated_mask_S_MD = dilate(((dilate(dilate(mask_S_MD_3D)))));

# ╔═╡ ed6c5395-2780-47d2-a3ca-322468638828
@bind n2 overlay_mask_bind(dilated_mask_S_MD)

# ╔═╡ f5599896-297e-4cac-89b0-2dd1344f665c
overlay_mask_plot(arr, dilated_mask_S_MD, n2, "dilated mask")

# ╔═╡ 64247b3c-eef8-4d19-9bf9-71cd10f2fc02
md"""
#### Ring (background) mask
"""

# ╔═╡ 452331b3-50d0-4a7f-a930-f832ca9aeb6e
ring_mask_S_MD = dilate(((dilate(dilate(mask_S_MD_3D))))) - dilate(dilate(((mask_S_MD_3D))));

# ╔═╡ b7fba435-4496-4671-a377-8b73b03949d4
@bind n4 overlay_mask_bind(ring_mask_S_MD)

# ╔═╡ 386ec444-7e85-458e-a213-d2f5e87374d3
overlay_mask_plot(arr, ring_mask_S_MD, n4, "ring mask")

# ╔═╡ cdfab961-1588-4830-9014-f1540dabe5ee
begin
	single_ring_mask_S_MD = Bool.(ring_mask_S_MD[:, :, 3])
	s_bkg_S_MD = mean(single_arr[single_ring_mask_S_MD])
end

# ╔═╡ dcd5144e-408f-4b13-995e-f3e716f27c3b
begin
	alg_S_MD = Integrated(arr[mask_S_MD_3D])
	mass_s_md = score(s_bkg_S_MD, S_Obj_MD, pixel_size, ρ_md, alg_S_MD)
end

# ╔═╡ 5d625a7e-f027-443b-aa97-726350616c11
md"""
## Low Density
"""

# ╔═╡ 7e5ec60e-33de-4963-a343-3f063091005c
begin
	mask_S_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_LD_3D[:, :, z] = mask_S_LD
	end
end;

# ╔═╡ 776fa172-d964-4c86-8220-725c11f6768c
md"""
#### Dilated mask
"""

# ╔═╡ 31301174-a282-4620-b243-ee386cdf8408
dilated_mask_S_LD = dilate(((dilate(dilate(mask_S_LD_3D)))));

# ╔═╡ d2dbe7a4-fc39-406b-9ad2-313b191dcf63
@bind o2 overlay_mask_bind(dilated_mask_S_LD)

# ╔═╡ cdbd4ccf-1fc4-4ff8-af3a-847113b85ef2
overlay_mask_plot(arr, dilated_mask_S_LD, o2, "dilated mask")

# ╔═╡ 79cf02a6-98f4-4d35-85a4-323c5d160c8b
md"""
#### Ring (background) mask
"""

# ╔═╡ 71e867c5-ecde-4643-9ffa-c6959246deb8
ring_mask_S_LD = (dilate((dilate(dilate(mask_S_LD_3D))))) - dilate(dilate(((mask_S_LD_3D))));

# ╔═╡ 4969b058-f3f8-4916-9ed9-c3b99f23f972
@bind o4 overlay_mask_bind(ring_mask_S_LD)

# ╔═╡ 9c3c3eb7-17b9-440d-b694-987f65aa0513
overlay_mask_plot(arr, ring_mask_S_LD, o4, "ring mask")

# ╔═╡ 1a66befd-3f95-443a-b2a5-ce4d40b97a75
begin
	single_ring_mask_S_LD = Bool.(ring_mask_S_LD[:, :, 3])
	s_bkg_S_LD = mean(single_arr[single_ring_mask_S_LD])
end

# ╔═╡ b7f78d93-dd53-45b5-b82d-4ffe3f8ffc03
begin
	alg_S_LD = Integrated(arr[mask_S_LD_3D])
	mass_s_ld = score(s_bkg_S_LD, S_Obj_LD, pixel_size, ρ_LD, alg_S_LD)
end

# ╔═╡ 7ba16b69-a495-4c18-bb10-d61ac48cde4b
md"""
# Results
"""

# ╔═╡ 5d176221-bb23-4ca3-84d3-accf691cca59
PhantomSegmentation.get_pixel_size(header)

# ╔═╡ b4dc3d85-068c-40ca-b121-b3591becb667
inserts = [
	"Low Density",
	"Medium Density",
	"High Density"
]

# ╔═╡ bb84ab31-0cb4-4f72-b6d6-8147de327030
volume_gt = [
	7.065,
	63.585,
	176.625
]

# ╔═╡ 85b0a6b1-c80e-439b-8913-79cfb87e80da
ground_truth_mass_large = [
	volume_gt[3] * density_array[2] * 1e-3,
	volume_gt[3] * density_array[3] * 1e-3,
	volume_gt[3] * density_array[4] * 1e-3
] # mg

# ╔═╡ d64d8e9d-0d9b-4ffb-853d-aa1f4c173b7c
calculated_mass_large = [
	mass_l_ld,
	mass_l_md,
	mass_l_hd
]

# ╔═╡ 529f7ae4-804e-4156-898a-1eb26978eaf2
ground_truth_mass_medium = [
	volume_gt[2] * density_array[2] * 1e-3,
	volume_gt[2] * density_array[3] * 1e-3,
	volume_gt[2] * density_array[4] * 1e-3
]

# ╔═╡ 96f0f0bd-aca4-4471-b6fd-d2dd3b1f7a99
calculated_mass_medium = [
	mass_m_ld,
	mass_m_md,
	mass_m_hd
]

# ╔═╡ f3593c6b-03df-4e06-949f-172de68590e8
ground_truth_mass_small = [
	volume_gt[1] * density_array[2] * 1e-3,
	volume_gt[1] * density_array[3] * 1e-3,
	volume_gt[1] * density_array[4] * 1e-3
]

# ╔═╡ a959ffa3-1d93-413c-9f87-4a426b3665e5
calculated_mass_small = [
	mass_s_ld,
	mass_s_md,
	mass_s_hd
]

# ╔═╡ 82540ec2-d054-40db-a561-8d7ecb254756
df = DataFrame(
	DENSITY = DENSITY,
	scan = scan,
	inserts = inserts,
	ground_truth_mass_large = ground_truth_mass_large,
	calculated_mass_large = calculated_mass_large,
	ground_truth_mass_medium = ground_truth_mass_medium,
	calculated_mass_medium = calculated_mass_medium,
	ground_truth_mass_small = ground_truth_mass_small,
	calculated_mass_small = calculated_mass_small,
	background = mass_bkg
)

# ╔═╡ 053072ba-1eb5-46de-be89-8bdd78658fcb
begin
	fmass2 = Figure()
	axmass2 = Axis(fmass2[1, 1])
	
	scatter!(density_array[2:end], df[!, :ground_truth_mass_large], label="ground_truth_mass_large")
	scatter!(density_array[2:end], df[!, :calculated_mass_large], label="calculated_mass_large")
	
	axmass2.title = "Mass Measurements (Large)"
	axmass2.ylabel = "Mass (mg)"
	axmass2.xlabel = "Density (mg/cm^3)"
	local xlim
	local ylim
	if DENSITY == "low"
		xlim = 130
		ylim = 30
	elseif DENSITY == "normal"
		xlim = 850
		ylim = 200
	end
	xlims!(axmass2, 0, xlim)
	ylims!(axmass2, 0, ylim)
	
	fmass2[1, 2] = Legend(fmass2, axmass2, framevisible = false)
	
	fmass2
end

# ╔═╡ b2d10283-d9ed-4001-8d72-dbded8e8d845
begin
	fmass3 = Figure()
	axmass3 = Axis(fmass3[1, 1])
	
	scatter!(density_array[2:end], df[!, :ground_truth_mass_medium], label="ground_truth_mass_medium")
	scatter!(density_array[2:end], df[!, :calculated_mass_medium], label="calculated_mass_medium")
	
	axmass3.title = "Mass Measurements (Medium)"
	axmass3.ylabel = "Mass (mg)"
	axmass3.xlabel = "Density (mg/cm^3)"

	local xlim
	local ylim
	if DENSITY == "low"
		xlim = 130
		ylim = 10
	elseif DENSITY == "normal"
		xlim = 850
		ylim = 70
	end
	xlims!(axmass3, 0, xlim)
	ylims!(axmass3, 0, ylim)
	
	fmass3[1, 2] = Legend(fmass3, axmass3, framevisible = false)
	
	fmass3
end

# ╔═╡ c68ee643-ab29-4afc-be97-820e7019a3be
begin
	fmass4 = Figure()
	axmass4 = Axis(fmass4[1, 1])
	
	scatter!(density_array[2:end], df[!, :ground_truth_mass_small], label="ground_truth_mass_small")
	scatter!(density_array[2:end], df[!, :calculated_mass_small], label="calculated_mass_small")
	
	axmass4.title = "Mass Measurements (Small)"
	axmass4.ylabel = "Mass (mg)"
	axmass4.xlabel = "Density (mg/cm^3)"
	
	local xlim
	local ylim
	if DENSITY == "low"
		xlim = 130
		ylim = 1.5
	elseif DENSITY == "normal"
		xlim = 850
		ylim = 10
	end
	xlims!(axmass4, 0, xlim)
	ylims!(axmass4, 0, ylim)
	
	fmass4[1, 2] = Legend(fmass4, axmass4, framevisible = false)
	
	fmass4
end

# ╔═╡ 7fd00263-a6d3-477c-b258-11b63b7a4e6f
percent_error_large = (abs.(ground_truth_mass_large - calculated_mass_large) ./ ground_truth_mass_large) .* 100

# ╔═╡ dfb7c91a-e360-4ae9-afbd-f9a96a891990
percent_error_medium = (abs.(ground_truth_mass_medium - calculated_mass_medium) ./ ground_truth_mass_medium) .* 100

# ╔═╡ 67b575e5-d3ab-4ae6-842c-6e3d96725dd1
percent_error_small= (abs.(ground_truth_mass_small - calculated_mass_small) ./ ground_truth_mass_small) .* 100

# ╔═╡ b12639a9-cbf1-422c-8d30-cf68426ea9c2
md"""
### Save Results
"""

# ╔═╡ b1e90509-0020-467d-a4c4-ae5b527a135c
# if ~isdir(string(cd(pwd, "..") , "/output/", TYPE))
# 	mkdir(string(cd(pwd, "..") , "/output/", TYPE))
# end

# ╔═╡ de4eb9c7-dbcc-4188-8dd0-8ba7eb544fcb
# output_path = string(cd(pwd, "..") , "/output/", TYPE, "/", scan, ".csv")

# ╔═╡ 22ce98f6-53cb-494b-a140-460538e4359a
# CSV.write(output_path, df)

# ╔═╡ 66eba716-0a7a-49fc-a54c-da28ef2312c9
md"""
### Save full df
"""

# ╔═╡ ccfc2cbb-990d-49ca-94b3-6ca8cd665535
dfs = []

# ╔═╡ ca4c048e-1fe0-481f-b4f3-a85223bb5ab4
push!(dfs, df)

# ╔═╡ 5d7f1e9e-500d-4812-b966-6b42948da3f6
# if length(dfs) == 24
# 	global new_df = vcat(dfs[1:24]...)
# 	output_path_new = string(cd(pwd, "..") , "/output/", TYPE, "/", "full.csv")
# 	CSV.write(output_path_new, new_df)
# end

# ╔═╡ f53ba93d-e1c0-4614-9efc-25de15d17cf7


# ╔═╡ Cell order:
# ╠═45f704d4-66e5-49db-aef5-05132f3853ee
# ╠═e7c4aaab-83ef-4256-b392-f8ef7a899a05
# ╟─48e3097f-0767-4dad-9d7c-94e0899f790b
# ╠═bc9383c0-e477-4e1a-a2fa-7f5c1d29f103
# ╟─1df1a491-6207-4b5c-a5e8-15a07323b6e7
# ╟─6aa51429-981a-4dea-a0f6-2935867d5b2a
# ╠═e849cf69-65c7-4e5e-9688-fc249d471f2c
# ╠═3e649fee-a2e9-4f0f-b705-3966f69d97ea
# ╠═a327ef18-2941-4783-9266-7332826eaf58
# ╠═b44d18f8-1b86-4235-bb5b-7a77a1af55e0
# ╠═b7e0b678-0627-44fe-b0cb-3ef2bccae6a7
# ╟─66f5fc66-d5ba-45ba-8624-55adc58085e4
# ╟─6e812172-6371-4461-9365-22f68ef16e53
# ╟─483a14dd-e798-41ed-9144-13678f8b8461
# ╟─2cfdeb06-71a2-411b-accb-2b4a3f56b477
# ╟─b3bae2ad-16ed-4381-b6a4-448cb5f0c6c1
# ╠═4ff19af0-7f53-4f24-b315-08bd54a488e3
# ╟─b74941c4-12d4-4be1-81f7-ef1f4d288984
# ╠═84e56c91-14a4-45b7-81a1-77e778bae695
# ╟─abfd965c-df00-4028-a0f3-8356a261b527
# ╟─41541a09-609b-4f46-9a25-eb974422efc0
# ╟─e6753c9a-d172-47fe-b066-78cb43df2c89
# ╟─f8a4eacd-2f69-4755-b017-5a0b0dda1004
# ╠═ac9d2652-6d69-4cf3-acbf-2e141fb633f0
# ╠═dd399c0e-f8e4-464c-8964-bdc0dd657202
# ╟─417d150e-0df9-4963-b59e-ebd9acf6d4a0
# ╠═5f10075b-4d95-4d53-a22f-016749fb7583
# ╟─b921dcf8-54ea-420b-9285-23e38d5ce433
# ╠═ece7d76b-93c8-430c-8679-07cf92585949
# ╠═ca408955-54ff-40ff-9a0f-77ee2666d7a2
# ╠═05ed60db-b4d8-4cdc-9a54-b108ade22557
# ╠═8bb2d451-b575-40bf-b31e-51e80392ef51
# ╠═4862496e-464d-40c3-baef-e5c507028d66
# ╟─997952cd-dc01-4528-a30a-fdc11a46dad8
# ╠═144d3db4-343b-44eb-bd13-160013b7a58c
# ╠═26b4ff9f-1b55-4fcc-99c8-3e1079d8df10
# ╠═c3a87ddd-64ad-4aeb-95e9-7f7b1ba1bd28
# ╠═da3c96b9-efdb-40f4-a959-747c097e16b2
# ╠═5fec186d-82bb-4b37-b0c6-c3781f2a166d
# ╠═0ce913ad-4844-4df0-92f2-79af4d85113e
# ╠═f221d7ac-7b27-4dc1-8081-af55d9e9b925
# ╠═c1d1b644-db3b-45b6-9050-342eca241a05
# ╠═d302f987-19aa-4401-9b79-20ca04103050
# ╠═03427ccf-69b9-4a1c-b988-30ad9f0a214c
# ╠═fa126b46-4dbb-4b8e-a3a7-7cf28a9ed88d
# ╠═d279e6b4-fba0-4fcb-a1ef-6c2c33ba8f63
# ╟─37e25efb-abf5-4568-a4d4-ff5c8aff1f64
# ╠═6c3f0dc5-9db2-4db8-abb0-64e4cda6d6fe
# ╠═6ceb2a0a-8fcf-4f8c-8167-5b2ab770df48
# ╠═2f06a3e2-5c65-4f7d-af65-829a8e657a92
# ╠═43cf1721-57db-4048-8528-b0430c5e8043
# ╠═444ec4ef-6140-4556-9305-599d819e58ca
# ╠═b6d389f3-6ecf-4798-978c-c28ef4b8e5da
# ╠═6c3c7d22-78e6-428e-ac62-90a473b06e2d
# ╠═f42a62d8-11e9-48b8-b6bc-ae25b6797707
# ╠═b407d628-0db0-4d13-9a6c-8e4ff5a84cbd
# ╠═dbec08a2-96c5-488b-89b9-5df6fe1f5a12
# ╟─e50baf62-0a02-47ee-8440-551f68baee0d
# ╠═3fa80bee-20f6-4e9c-a14d-75a77482a2c3
# ╠═04c77c9d-c498-4e09-abff-9c361e887c56
# ╟─b03c711f-9e0e-4089-add1-5abbde81cce1
# ╟─3b30e707-02b0-49c2-a0b3-25b7fcbf79ee
# ╠═d27df860-7ca8-4e93-90e4-a58605aa1aaa
# ╠═5dfb50b9-7545-4948-8026-2074ce5d86cf
# ╟─177ab471-2783-4e83-95fd-fcd8d02cb104
# ╟─aa02f7d8-6b01-4da9-b0f5-d1df587c7d4d
# ╠═dec0f541-4f89-45e3-b084-ecdbc837b794
# ╠═3895fafc-a999-4b06-974d-caef190f7830
# ╟─94c66847-765d-4d58-8510-043285d6cdcd
# ╠═703c7556-3904-41d7-972b-48ea31890a7d
# ╟─983fa5b0-45fe-4ba9-8eea-ad5c030ed32f
# ╠═a7de43e0-795d-4267-b7dd-3be85a201485
# ╟─2a6961fc-f94f-4250-aeb8-63bbfd9e29a0
# ╠═db4992af-e78b-4f84-b0b3-c13b6cd6cde5
# ╟─58172907-2026-4caa-a6ad-ab484b86b329
# ╠═f8e646a2-8b1e-47f6-bc5d-760941d542e8
# ╟─080422ff-77a2-4b86-8a0f-c2bea1f0d722
# ╠═1cba2d1a-19c2-49b7-8757-5843bb16a0f1
# ╟─cf4bc01d-5603-4d47-8863-e413be42f223
# ╠═ccf2c3d2-7243-412f-bbe1-84b4ac3f6272
# ╟─b990d904-e359-45fe-990b-87ab625ba1b0
# ╠═0af76504-ffc9-44ac-aecd-a3b956166273
# ╠═dc94ad4a-0692-4e0e-84e5-f706ab1ca316
# ╠═fd089ad4-4e08-4e1c-b049-3378de42b2df
# ╠═28c869df-75eb-42dc-9033-867ca4782f41
# ╠═7d61e536-fac5-43e1-9ce9-5952d7d20e8d
# ╟─e93b26b5-ba84-4410-87b7-eefad3d8ab3e
# ╠═1b97a53d-6ed8-4414-befa-146225f0f95a
# ╟─336cf5d0-ada0-4f9a-9788-bfd0ab72087a
# ╠═94779096-626f-4cb6-9937-0bc4bc8ab45b
# ╟─937a5188-e51e-4017-be9a-2d8051e79fd1
# ╠═f4ee897f-d6ad-4dfd-a91f-a759ecd170e5
# ╟─afbc96dc-f5fa-4305-bd3d-119ef8d416a8
# ╠═078b6320-34d9-4318-863e-a2a6dcdeb500
# ╟─8e310e91-c557-48e5-a5d1-042c1ac8f7b3
# ╠═ed5c5da5-21dc-4052-bc18-6b209edeabfd
# ╠═40c7c265-f558-4f74-8d96-40cd83c2d5bb
# ╠═0bf239a7-233c-42f2-b7eb-29336b472ff9
# ╠═367750d8-a287-4313-b6df-d6533849623f
# ╟─73105b5a-05b7-429f-9a7f-b8c12684df91
# ╠═31e9a969-6514-404c-89f6-0f038cc115f2
# ╟─d8d7a93f-85e7-43e7-a23d-11e7c6ce27be
# ╠═466cfc6f-2915-44b4-9231-4d0047accd6d
# ╟─b849a64f-965a-40b7-b2f5-1a27fd8ceeb5
# ╠═7f52b670-611d-4bd3-bc80-77ba59aa2d70
# ╟─4c425cf3-2ece-4df1-b3ab-4058c7d7aa13
# ╠═52b9efbd-80b3-4205-b4a7-b18809d6ba24
# ╟─7a618b0c-9429-472a-b076-774b9d229cab
# ╠═544151a5-3b47-4571-a6f0-9800d5b11085
# ╠═9da779a0-86bd-450d-9401-7d3507aaec7d
# ╠═187c110a-21a9-4b31-82b6-d9c239dbd899
# ╠═7560dc46-fa36-4706-925c-0675397c5922
# ╟─ee7092e6-52fb-4a6a-ad21-90cd2b46c7fa
# ╟─62faf81f-eb3b-431a-8e2b-239735db3257
# ╠═6fba5f51-7940-45fe-8a0b-30ea12962a0c
# ╟─10765555-8507-413f-a582-1f47619f0181
# ╠═f8519d1a-c5e7-41a3-90d3-5f6200f767cf
# ╟─5aff6e7a-b635-44d7-a91e-1a30c5f4fd17
# ╠═4aabe7d9-e730-4b4f-9dd2-e51e922ac08c
# ╟─706c30df-ce63-4c81-9c99-227d56ab2991
# ╠═6598b28c-5850-44bc-8bd0-995ff3d59f09
# ╟─47615958-a8fc-4fc2-852f-e11f69041a86
# ╠═85097de3-d336-4c85-a759-ad1b3e4083cc
# ╠═0f02e277-06b7-46f7-b6e8-0d5b6516dbe2
# ╠═bd4cec3c-446f-48e4-b06b-fbde4848daaf
# ╟─0f71036e-b118-430d-bce7-420b2591f1a2
# ╠═539d9e60-1f76-47c5-ab4f-cb520c520965
# ╟─c99672f4-9ee7-4d42-8a15-eee68be2a8da
# ╠═c34e32c6-4f28-435b-b124-7d3727f107a7
# ╟─b1c10feb-3148-4de5-b305-b0721894b76e
# ╠═bd91bd0d-f368-4340-916a-ea197969c0e4
# ╟─b5cb8919-df84-48a6-b00b-cf7dba301c76
# ╠═48eb8f20-1d00-4434-8dcd-4960b2392d28
# ╟─6c3c8eb5-c651-45e3-a912-8fd0dc330b51
# ╠═f2bafa39-c1ac-42ca-89c0-df4d219c2c74
# ╠═906c678f-2006-4f66-a135-97e25182a3d6
# ╠═43afbb30-b1b5-452b-8c82-25d173f24699
# ╟─386ea37d-70ba-4530-be9e-549d5e0ab1ce
# ╠═25434519-be6d-4616-9353-a0d342129ea4
# ╟─b37af6a8-8f51-4263-ac85-bea25fd2aae2
# ╠═b9dbf6c2-891f-490d-a89c-bd50891cfdee
# ╟─7b8af8e8-e979-4b57-a65a-ad99a40f08b9
# ╠═493ef94e-d6ac-4ccb-b342-0812792c2536
# ╟─43f5bf4e-c882-4596-a4b5-6063ce63e17c
# ╠═c44aa6b4-bada-4ab0-a748-71d0e215c1c7
# ╟─48a0f2c4-0bea-4c7d-bf8a-dc3632a490fa
# ╠═e03ccab0-b786-49ca-b7dd-b2b7f76e04a2
# ╠═5133060e-c1fb-492f-93d7-f161e8da73be
# ╠═907b3f24-47b8-4519-b70a-de44ca764bb2
# ╟─c52478ca-2012-4690-bc39-6d32997427f7
# ╟─5de67b93-72d7-4f86-b94e-245f4cfff576
# ╠═6323a01d-100a-4d09-b3ad-5821e88bed0a
# ╟─d0ea5ea6-5c98-44b8-94d4-4104eb1f6080
# ╠═594ed641-0f62-461d-b031-2a5de1467d35
# ╟─018d2392-7133-4da8-b4c9-268af6fa6f3a
# ╠═9d72b647-73d9-45c7-b473-bbdea2817cb5
# ╟─27ffaff1-12c5-43f1-ac21-84fe25877f24
# ╠═a4ddb30d-011c-4855-96f3-01ba441297f1
# ╟─46110aad-07a8-4a56-9a6d-db71054a5883
# ╠═070aedea-d5e9-43f7-a2cf-14fdad760a72
# ╠═76b6bf9f-58fb-4a88-a21f-0695d2c76073
# ╠═1590d3c5-d9f1-4c00-9730-639677635056
# ╟─aa3f6b52-2c1a-462a-a89c-bcde6f0f59ba
# ╠═e17a81a4-35a6-4d15-bf14-275114a49608
# ╟─01dce72f-6793-4360-80fc-bc2e97d579de
# ╠═5250fe6b-b94a-47b0-9142-525ca34f1d16
# ╟─ed6c5395-2780-47d2-a3ca-322468638828
# ╠═f5599896-297e-4cac-89b0-2dd1344f665c
# ╟─64247b3c-eef8-4d19-9bf9-71cd10f2fc02
# ╠═452331b3-50d0-4a7f-a930-f832ca9aeb6e
# ╟─b7fba435-4496-4671-a377-8b73b03949d4
# ╠═386ec444-7e85-458e-a213-d2f5e87374d3
# ╠═cdfab961-1588-4830-9014-f1540dabe5ee
# ╠═dcd5144e-408f-4b13-995e-f3e716f27c3b
# ╟─5d625a7e-f027-443b-aa97-726350616c11
# ╠═7e5ec60e-33de-4963-a343-3f063091005c
# ╟─776fa172-d964-4c86-8220-725c11f6768c
# ╠═31301174-a282-4620-b243-ee386cdf8408
# ╟─d2dbe7a4-fc39-406b-9ad2-313b191dcf63
# ╠═cdbd4ccf-1fc4-4ff8-af3a-847113b85ef2
# ╟─79cf02a6-98f4-4d35-85a4-323c5d160c8b
# ╠═71e867c5-ecde-4643-9ffa-c6959246deb8
# ╟─4969b058-f3f8-4916-9ed9-c3b99f23f972
# ╠═9c3c3eb7-17b9-440d-b694-987f65aa0513
# ╠═1a66befd-3f95-443a-b2a5-ce4d40b97a75
# ╠═b7f78d93-dd53-45b5-b82d-4ffe3f8ffc03
# ╟─7ba16b69-a495-4c18-bb10-d61ac48cde4b
# ╠═5d176221-bb23-4ca3-84d3-accf691cca59
# ╠═b4dc3d85-068c-40ca-b121-b3591becb667
# ╠═bb84ab31-0cb4-4f72-b6d6-8147de327030
# ╠═85b0a6b1-c80e-439b-8913-79cfb87e80da
# ╠═d64d8e9d-0d9b-4ffb-853d-aa1f4c173b7c
# ╠═529f7ae4-804e-4156-898a-1eb26978eaf2
# ╠═96f0f0bd-aca4-4471-b6fd-d2dd3b1f7a99
# ╠═f3593c6b-03df-4e06-949f-172de68590e8
# ╠═a959ffa3-1d93-413c-9f87-4a426b3665e5
# ╠═82540ec2-d054-40db-a561-8d7ecb254756
# ╟─053072ba-1eb5-46de-be89-8bdd78658fcb
# ╟─b2d10283-d9ed-4001-8d72-dbded8e8d845
# ╟─c68ee643-ab29-4afc-be97-820e7019a3be
# ╠═7fd00263-a6d3-477c-b258-11b63b7a4e6f
# ╠═dfb7c91a-e360-4ae9-afbd-f9a96a891990
# ╠═67b575e5-d3ab-4ae6-842c-6e3d96725dd1
# ╟─b12639a9-cbf1-422c-8d30-cf68426ea9c2
# ╠═b1e90509-0020-467d-a4c4-ae5b527a135c
# ╠═de4eb9c7-dbcc-4188-8dd0-8ba7eb544fcb
# ╠═22ce98f6-53cb-494b-a140-460538e4359a
# ╟─66eba716-0a7a-49fc-a54c-da28ef2312c9
# ╠═ccfc2cbb-990d-49ca-94b3-6ca8cd665535
# ╠═ca4c048e-1fe0-481f-b4f3-a85223bb5ab4
# ╠═5d7f1e9e-500d-4812-b966-6b42948da3f6
# ╠═f53ba93d-e1c0-4614-9efc-25de15d17cf7

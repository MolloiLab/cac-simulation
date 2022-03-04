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

# ╔═╡ 11bc4d31-8569-4ddf-9833-017f6b4cca66
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
		Pkg.add(url="https://github.com/JuliaHealth/DICOM.jl")
		Pkg.add(url="https://github.com/Dale-Black/DICOMUtils.jl")
		Pkg.add(url="https://github.com/Dale-Black/Phantoms.jl")
		Pkg.add(url="https://github.com/Dale-Black/CalciumScoring.jl")
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
	using DICOM
	using DICOMUtils
	using Phantoms
	using CalciumScoring
end

# ╔═╡ 201127f2-ae5e-4542-9d9f-c134e16d2c29
TableOfContents()

# ╔═╡ bfda24da-a8cf-4833-893e-2e039d9bcaf1
md"""
## Load DICOMS

All you need to do is set `base_path` once and leave it. After that, the only thing that should change is the `VENDER`, once for every set, and the `SCAN_NUMBER`, once for each scan.
"""

# ╔═╡ bdd77f47-b5ff-4855-b3c4-1fcd5fe31777
begin
	SCAN_NUMBER = 1
	VENDER = "135"
	BASE_PATH = "/Users/daleblack/Google Drive/Datasets/Simulated/"
end

# ╔═╡ 773eff82-436b-4b36-8203-0764a669ce62
md"""
**Everything below should be automatic, just scroll through to visually inspect that things make sense**
"""

# ╔═╡ 377b58ac-28ff-4bfe-8ce2-c21481e8e973
root_path = string(BASE_PATH, VENDER)

# ╔═╡ e2fffe13-12b0-4645-aecd-b2e2359eb9d6
dcm_path_list = dcm_list_builder(root_path)

# ╔═╡ a023991a-6814-4375-a0b7-8b4c70c7f2b4
pth = dcm_path_list[SCAN_NUMBER]

# ╔═╡ 5b02e50f-0f8d-4987-aaf8-08cadf89562d
scan = basename(pth)

# ╔═╡ ebb232b1-a4ba-441f-a79f-23d8489aff2b
header, dcm_array, slice_thick_ori1 = dcm_reader(pth);

# ╔═╡ 604fd08b-6677-428f-a263-84bf9c2eedac
md"""
## Helper Functions
"""

# ╔═╡ 041a5099-8278-428f-9d00-db6295f32897
function collect_tuple(tuple_array)
	row_num = size(tuple_array)
	col_num = length(tuple_array[1])
	container = zeros(Int64, row_num..., col_num)
	for i in 1:length(tuple_array)
		container[i,:] = collect(tuple_array[i])
	end
	return container
end

# ╔═╡ c7868bf0-611f-4d02-8247-b87edc82313b
function overlay_mask_bind(mask)
	indices = findall(x -> x == 1, mask)
	indices = Tuple.(indices)
	label_array = collect_tuple(indices)
	zs = unique(label_array[:,3])
	return PlutoUI.Slider(1:length(zs), default=3, show_value=true)
end

# ╔═╡ 9ca9a957-17e3-4b9a-9f3b-e1d042ad0bcb
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

# ╔═╡ 8723bd5f-6526-43b3-9de0-95c42fd30bbb
md"""
## Segment Heart
"""

# ╔═╡ 07f993c4-43fd-4960-934f-abba867e8231
masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3)÷2);

# ╔═╡ de654d7a-ecba-4526-afa3-417cc5820b98
@bind a PlutoUI.Slider(1:size(masked_array, 3), default=10, show_value=true)

# ╔═╡ 3aec1071-d863-46b2-a66e-3d8500cf383c
heatmap(masked_array[:, :, a], colormap=:grays)

# ╔═╡ b01b496c-b80b-424f-a11f-d8bd5205ebbd
begin
	fig = Figure()
	
	ax = Makie.Axis(fig[1, 1])
	ax.title = "Raw DICOM Array"
	heatmap!(transpose(dcm_array[:, :, 4]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig
end

# ╔═╡ bb1df3e9-793b-4727-b806-5962c38c21c3
begin
	fig2 = Figure()
	
	ax2 = Makie.Axis(fig2[1, 1])
	ax2.title = "Mask Array"
	heatmap!(transpose(mask), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig2
end

# ╔═╡ 38c984b0-f132-42c6-bd99-27802bc5d191
begin
	fig3 = Figure()
	
	ax3 = Makie.Axis(fig3[1, 1])
	ax3.title = "Masked DICOM Array"
	heatmap!(transpose(masked_array[:, :, 1]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig3
end

# ╔═╡ e9e2cc57-61bc-479f-bddd-fd5f930ffa65
md"""
## Segment Calcium Rod
"""

# ╔═╡ c19dff18-eadb-43ee-bfb4-4fa3e573d1e4
calcium_image, slice_CCI, quality_slice, cal_rod_slice = mask_rod(masked_array, header);

# ╔═╡ c8f1c6d7-3717-4b14-9bb7-d72a37849ea8
@bind c PlutoUI.Slider(1:size(calcium_image, 3), default=cal_rod_slice, show_value=true)

# ╔═╡ 301e3bb2-3302-4251-931c-e32211136eea
heatmap(transpose(calcium_image[:, :, c]), colormap=:grays)

# ╔═╡ e60c281b-46cd-4a94-90a1-de0fdbb05491
md"""
## Segment Calcium Inserts
"""

# ╔═╡ 73ede908-a805-4d6b-93ae-3fc724b56423
mask_L_HD, mask_M_HD, mask_S_HD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_LD, mask_M_LD, mask_S_LD = mask_inserts_simulation(
            dcm_array, masked_array, header, slice_CCI, center_insert
);

# ╔═╡ 0e39214d-99c2-4af3-8f8e-78ab87f7fef6
slice_CCI

# ╔═╡ 64d526d4-b5a7-44cd-9db0-7b5265b37857
masks = mask_L_HD + mask_M_HD + mask_S_HD + mask_L_MD + mask_M_MD + mask_S_MD + mask_L_LD + mask_M_LD + mask_S_LD;

# ╔═╡ bf16de9b-3e2b-42ca-be6c-a40f2c2b87be
heatmap(masks, colormap=:grays)

# ╔═╡ 27cc2663-2080-48f6-a3a0-7fae64868ab7
md"""
## Calibration Prep
"""

# ╔═╡ 0a85351b-74a7-47ec-8baa-fe9f2b6e81f9
array_filtered = abs.(mapwindow(median, calcium_image[:, :, 2], (3, 3)));

# ╔═╡ ac5895f8-e629-4500-a895-b781ed76a203
bool_arr = array_filtered .> 0;

# ╔═╡ 5f8bc4e3-3bd3-456a-8b0f-55bd5912d81c
bool_arr_erode = (((erode(erode(bool_arr)))));

# ╔═╡ d1c93fea-1ced-4a68-b66a-9d8978ba7a2a
heatmap(bool_arr, colormap=:grays)

# ╔═╡ 6b357f4c-73a7-44ab-a554-76dfa0ebfb34
heatmap(bool_arr_erode, colormap=:grays)

# ╔═╡ f1355c86-cdc4-4f76-967a-5a6285ecbac6
c_img = calcium_image[:, :, 1:3];

# ╔═╡ 5ed7dba0-b78b-4484-94b3-0b34ae94a305
begin
	mask_cal_3D = Array{Bool}(undef, size(c_img))
	for z in 1:size(c_img, 3)
		mask_cal_3D[:, :, z] = bool_arr_erode
	end
end;

# ╔═╡ 4cc8d210-44fe-4524-aeb7-0688a383989d
hist(c_img[mask_cal_3D])

# ╔═╡ 9ce3d6a3-dbb5-480e-a34f-5084b5572d7c
cal_insert_mean = quantile!(c_img[mask_cal_3D], 0.7)

# ╔═╡ 66582769-d393-48e6-826e-98337ebbe219
md"""
### Calibration Line
"""

# ╔═╡ 035fd1c1-dadf-438f-b48b-7f6e3313e9dc
density_array = [0, 200, 400, 800]

# ╔═╡ 6a5c9b4d-c8ad-44e0-bcc5-e4509d266ba5
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

# ╔═╡ ca9554e3-147c-403a-8027-c286c1cdcd5f
md"""
# Score Large Inserts
"""

# ╔═╡ 2f9e82ff-ede9-4d0b-997e-2c814d997527
arr = masked_array[:, :, 4:6];

# ╔═╡ 6220799c-323e-49b3-83be-f828d75dfc4d
single_arr = masked_array[:, :, slice_CCI];

# ╔═╡ 3bedd4e1-696a-4a40-a18d-2e2575402399
md"""
## High Density
"""

# ╔═╡ d730dac7-ae9c-4102-ab32-d24e4d0c7709
begin
	mask_L_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_HD_3D[:, :, z] = mask_L_HD
	end
end;

# ╔═╡ 9368bf8f-54de-4619-86c5-cdd01e44a555
begin
	eroded_mask_L_HD = erode(erode(mask_L_HD_3D))
	high_density_cal = mean(arr[eroded_mask_L_HD])
end

# ╔═╡ b1679304-f66b-4199-9911-d9d02a757b62
md"""
#### Dilated mask
"""

# ╔═╡ 212554d6-71a2-47fd-8002-8edb87e1dce5
dilated_mask_L_HD = dilate(dilate(mask_L_HD_3D));

# ╔═╡ c4098846-e781-4a52-b280-f4a313c2cd2a
@bind g2 overlay_mask_bind(dilated_mask_L_HD)

# ╔═╡ e7efea19-74d6-4513-be91-f010c81b9fd2
overlay_mask_plot(arr, dilated_mask_L_HD, g2, "dilated mask")

# ╔═╡ 143c8c29-0063-4868-bf4e-8fd3c5c08dfb
md"""
#### Ring (background) mask
"""

# ╔═╡ e77a2392-1fbc-4fb0-8bef-d7a9283782be
ring_mask_L_HD = dilate(dilate(dilate(dilate(mask_L_HD_3D)))) - dilate(dilate(dilate(mask_L_HD_3D)));

# ╔═╡ 6784b04f-ece0-447c-8820-a0644eaf823e
@bind g4 overlay_mask_bind(ring_mask_L_HD)

# ╔═╡ 876eae18-f388-4f9b-bff5-4ec58657d3a0
overlay_mask_plot(arr, ring_mask_L_HD, g4, "ring mask")

# ╔═╡ 889052c7-d431-49f8-aa71-bd2c59aed37f
md"""
### Calculations
"""

# ╔═╡ ada32081-4eaf-4189-801a-f91d9c2f189a
begin
	single_ring_mask_L_HD = Bool.(ring_mask_L_HD[:, :, 3])
	s_bkg_L_HD = mean(single_arr[single_ring_mask_L_HD])
end

# ╔═╡ e5cdb203-85ba-4217-8b6c-4a9a782d54ad
pixel_size = DICOMUtils.get_pixel_size(header)

# ╔═╡ 696a45ad-aa73-4daa-b8f6-81584e1cfb9e
md"""
## Medium Density
"""

# ╔═╡ c3d0bb9e-99dc-49d2-b30d-a2534e300a0f
begin
	mask_L_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_MD_3D[:, :, z] = mask_L_MD
	end
end;

# ╔═╡ e20edcfd-1d64-44b5-9b08-27a77f44299f
begin
	eroded_mask_L_MD = erode(erode(mask_L_MD_3D))
	med_density_cal = mean(arr[eroded_mask_L_MD])
end

# ╔═╡ aac2fcd5-9ab1-42da-8464-76e1a3fbe443
md"""
#### Dilated mask
"""

# ╔═╡ b4f9c38e-dd42-49d2-9257-3d2c13653878
dilated_mask_L_MD = dilate(dilate(mask_L_MD_3D));

# ╔═╡ 53dcc1a6-1cc1-4d24-bd13-d84d9335b1fa
@bind h2 overlay_mask_bind(dilated_mask_L_MD)

# ╔═╡ b720e961-a907-45db-b86b-3330a809da4f
overlay_mask_plot(arr, dilated_mask_L_MD, h2, "dilated mask")

# ╔═╡ 75d2f8af-13e0-425d-a40b-c8cffb371e7b
md"""
#### Ring (background) mask
"""

# ╔═╡ edaceb23-784f-4149-8b74-df3089208a62
ring_mask_L_MD = dilate(dilate(dilate(dilate(mask_L_MD_3D)))) - dilate(dilate(dilate(mask_L_MD_3D)));

# ╔═╡ 1d3bab04-b000-4d71-a729-67620a4b4e73
@bind h4 overlay_mask_bind(ring_mask_L_MD)

# ╔═╡ 01333998-bf78-4b97-b45f-6a7b1ef22ea3
overlay_mask_plot(arr, ring_mask_L_MD, h4, "ring mask")

# ╔═╡ 533aa777-c1cc-4067-b6de-2019bb11ee77
md"""
### Calculations
"""

# ╔═╡ fe50b0e4-af31-441b-a513-cfc914584cb2
begin
	single_ring_mask_L_MD = Bool.(ring_mask_L_MD[:, :, 3])
	s_bkg_L_MD = mean(single_arr[single_ring_mask_L_MD])
end

# ╔═╡ 73047cbf-17df-430b-8a8f-69af936a89af
md"""
## Low Density
"""

# ╔═╡ a5767c29-300c-40d7-aeff-9f373f211b7a
begin
	mask_L_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_LD_3D[:, :, z] = mask_L_LD
	end
end;

# ╔═╡ e1d6755a-bb5f-44a9-abc8-240beac0f5ee
begin
	eroded_mask_L_LD = erode(erode(mask_L_LD_3D))
	low_density_cal = mean(arr[eroded_mask_L_LD])
end

# ╔═╡ 283d2191-06cd-40f3-bab0-d66f1ded8d27
intensity_array = [0, low_density_cal, med_density_cal, high_density_cal] # HU

# ╔═╡ b4d32db5-b5be-4fb4-9666-2bfe0aa75d65
df_cal = DataFrame(:density => density_array, :intensity => intensity_array)

# ╔═╡ 7803d8e7-a3fe-4043-af64-4e9b334c8f82
linearRegressor = lm(@formula(intensity ~ density), df_cal);

# ╔═╡ ee1bfe95-106c-4e43-846e-6293a1413830
linearFit = predict(linearRegressor)

# ╔═╡ 22540057-0179-49bb-8437-092839d4e7b5
m = linearRegressor.model.pp.beta0[2]

# ╔═╡ 8c3a1111-c5c6-4bd4-ad14-e1c578471354
b = linearRegressor.model.rr.mu[1]

# ╔═╡ 7b76786b-3cc7-447d-b78f-007d6fdc16dc
density(intensity) = (intensity - b) / m

# ╔═╡ 87fa7e66-945c-442a-8cd3-e4e302379074
intensity(ρ) = m*ρ + b

# ╔═╡ 0499d5a5-22a2-48b4-afc3-59a03eed10df
S_Obj_HD = intensity(800)

# ╔═╡ 60ccc705-6931-40b8-8f46-955398fe13c9
begin
	alg_L_HD = Integrated(arr[mask_L_HD_3D])
	ρ_hd = 0.8 # mg/mm^3
	mass_l_hd = score(s_bkg_L_HD, S_Obj_HD, pixel_size, ρ_hd, alg_L_HD)
end

# ╔═╡ 4c24aa2b-4c92-435a-9bac-1240741126fc
S_Obj_MD = intensity(400)

# ╔═╡ 4fd5d0ae-87b3-41a9-a4a8-7b4687592e2e
begin
	alg_L_MD = Integrated(arr[mask_L_MD_3D])
	ρ_md = 0.4
	mass_l_md = score(s_bkg_L_MD, S_Obj_MD, pixel_size, ρ_md, alg_L_MD)
end

# ╔═╡ 1ef692d4-91e0-44cc-a3bc-0002f145ca80
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

# ╔═╡ 64cd83f3-f9ed-4215-8bd2-be3453bcc4dc
md"""
#### Dilated mask
"""

# ╔═╡ 607aed18-2914-4b48-a2cf-e2f4e78524d2
dilated_mask_L_LD = dilate(dilate(mask_L_LD_3D));

# ╔═╡ f0fdd9ae-8b61-47b5-8826-923d1779bc8b
@bind i2 overlay_mask_bind(dilated_mask_L_LD)

# ╔═╡ e1c1f280-a844-4ccd-8237-0762b0ffc7c4
overlay_mask_plot(arr, dilated_mask_L_LD, i2, "dilated mask")

# ╔═╡ 44e8196c-6aa5-4797-b7a2-289325633c84
md"""
#### Ring (background) mask
"""

# ╔═╡ 422fdee4-6685-46e3-961b-64425f205ec8
ring_mask_L_LD = dilate(dilate(dilate(dilate(mask_L_LD_3D)))) - dilate(dilate(dilate(mask_L_LD_3D)));

# ╔═╡ e5e751bc-cee3-4593-9978-68578bc4021c
@bind i4 overlay_mask_bind(ring_mask_L_LD)

# ╔═╡ 161a8f24-e102-4ceb-b40d-77756139f089
overlay_mask_plot(arr, ring_mask_L_LD, i4, "ring mask")

# ╔═╡ 9b089a9d-6037-4bed-bf11-0c494aef3581
md"""
### Calculations
"""

# ╔═╡ b834eae9-a7c8-4ad3-8161-2d3463fa99d7
begin	
	single_ring_mask_L_LD = Bool.(ring_mask_L_LD[:, :, 3])
	s_bkg_L_LD = mean(single_arr[single_ring_mask_L_LD])
end

# ╔═╡ f6ca11d6-ff98-49f6-98e4-046eaa77931f
S_Obj_LD = intensity(200)

# ╔═╡ aec7f59f-b0d0-4a6e-96b2-94bc2192b2d2
begin
	alg_L_LD = Integrated(arr[mask_L_LD_3D])
	ρ_LD = 0.2
	mass_l_ld = score(s_bkg_L_LD, cal_insert_mean, pixel_size, ρ_LD, alg_L_LD)
end

# ╔═╡ 1a998d70-1867-4571-b534-5ef8b061e60d
md"""
# Score Medium Inserts
"""

# ╔═╡ 653c8222-07f7-4633-9bf1-b7b02e67b3c3
md"""
## High Density
"""

# ╔═╡ 4e5151cb-7acc-4331-b011-a4cd17f84f78
begin
	mask_M_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_HD_3D[:, :, z] = mask_M_HD
	end
end;

# ╔═╡ 6ed7f390-de7e-4489-bc82-d4a214baca2d
md"""
#### Dilated mask
"""

# ╔═╡ 6d092d23-46e0-4a47-a945-009a48471256
dilated_mask_M_HD = dilate(dilate(dilate(dilate(mask_M_HD_3D))));

# ╔═╡ 4caade91-81bf-4bfe-bbe0-cb428c57e735
@bind j2 overlay_mask_bind(dilated_mask_M_HD)

# ╔═╡ 84b41be9-f685-4a36-9b18-ed25a5c283d9
overlay_mask_plot(arr, dilated_mask_M_HD, j2, "dilated mask")

# ╔═╡ 8a5f08e3-0ef3-4b87-93b1-649328e07d96
md"""
#### Ring (background) mask
"""

# ╔═╡ 3a930ac5-8cdf-4677-9303-6e13eef7bbef
ring_mask_M_HD = dilate(dilate(dilate(dilate(dilate(mask_M_HD_3D))))) - dilate(dilate(dilate(dilate(mask_M_HD_3D))));

# ╔═╡ 7651bc64-161f-46e5-8947-a58840cae9c0
@bind j4 overlay_mask_bind(ring_mask_M_HD)

# ╔═╡ d80c7da7-3fb0-48ad-ac51-84f6005e5276
overlay_mask_plot(arr, ring_mask_M_HD, j4, "ring mask")

# ╔═╡ ded2cd96-0bb2-4718-8503-ca3361cdf168
md"""
### Calculations
"""

# ╔═╡ 7833403f-c17f-4346-bcd5-c4f66664e348
begin
	single_ring_mask_M_HD = Bool.(ring_mask_M_HD[:, :, 3])
	s_bkg_M_HD = mean(single_arr[single_ring_mask_M_HD])
end

# ╔═╡ ee3bfc69-15d0-4695-a1e3-3a72060f91de
begin
	alg_M_HD = Integrated(arr[mask_M_HD_3D])
	mass_m_hd = score(s_bkg_M_HD, S_Obj_HD, pixel_size, ρ_hd, alg_M_HD)
end

# ╔═╡ 4214559c-8c4e-41f9-a003-f926bedb1dda
md"""
## Medium Density
"""

# ╔═╡ 5044a80d-fa5d-413f-98fc-ba8208280803
begin
	mask_M_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_MD_3D[:, :, z] = mask_M_MD
	end
end;

# ╔═╡ fc119ac1-f6ad-4ef5-b641-a7961938b4cd
md"""
#### Dilated mask
"""

# ╔═╡ 843bbb31-f437-453c-ba61-864fa4234257
dilated_mask_M_MD = dilate(dilate(dilate(dilate(mask_M_MD_3D))));

# ╔═╡ 7ff5e49e-6f18-48d0-b97f-1abd9586867a
@bind k2 overlay_mask_bind(dilated_mask_M_MD)

# ╔═╡ aa70689c-adb7-45d8-935c-40ea9b65901e
overlay_mask_plot(arr, dilated_mask_M_MD, k2, "dilated mask")

# ╔═╡ b2974709-0f5f-4c14-b902-6021a624dc7c
md"""
#### Ring (background) mask
"""

# ╔═╡ 95d5c859-daf3-49b3-87e9-8b45b226d99b
ring_mask_M_MD = dilate(dilate(dilate(dilate(dilate(dilate(mask_M_MD_3D)))))) - dilate(dilate(dilate(dilate(dilate(mask_M_MD_3D)))));

# ╔═╡ 1ccae49d-9d4d-49ed-baeb-232eb31f1aae
@bind k4 overlay_mask_bind(ring_mask_M_MD)

# ╔═╡ 790998a8-542c-46f9-a850-0d07034bf22c
overlay_mask_plot(arr, ring_mask_M_MD, k4, "ring mask")

# ╔═╡ b4cc3468-9173-4168-8149-bd29127be392
md"""
### Calculations
"""

# ╔═╡ 1d16adec-06f5-4816-a277-cadef49ca558
begin
	single_ring_mask_M_MD = Bool.(ring_mask_M_MD[:, :, 3])
	s_bkg_M_MD = mean(single_arr[single_ring_mask_M_MD])
end

# ╔═╡ dcf74ca2-cf93-4f6d-a53c-b47b8f3435a3
begin
	alg_M_MD = Integrated(arr[mask_M_MD_3D])
	mass_m_md = score(s_bkg_M_MD, S_Obj_MD, pixel_size, ρ_md, alg_M_MD)
end

# ╔═╡ 08837a72-b43a-4b36-80ba-92cd64e59c79
md"""
## Low Density
"""

# ╔═╡ 67bd5d43-7d7e-4843-9288-61e67b22a5cf
begin
	mask_M_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_LD_3D[:, :, z] = mask_M_LD
	end
end;

# ╔═╡ d01bfc39-f257-482c-957f-f86ea7b5fbcc
md"""
#### Dilated mask
"""

# ╔═╡ 6b2c2f08-cdfd-45d4-ab2f-f23a5e1bb1a6
dilated_mask_M_LD = dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))));

# ╔═╡ 52eb3999-6403-43df-96de-0b2fa3930b1f
@bind l2 overlay_mask_bind(dilated_mask_M_LD)

# ╔═╡ 73ea58c2-77cc-451f-921b-d6db7c9ce56f
overlay_mask_plot(arr, dilated_mask_M_LD, l2, "dilated mask")

# ╔═╡ 35e5cc8e-bcc5-4365-a9ce-f75591af4520
md"""
#### Ring (background) mask
"""

# ╔═╡ 31adbea2-ccc3-4fa7-ae32-bfb92e8d0fbe
ring_mask_M_LD = dilate(dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))))) - dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))));

# ╔═╡ 37cd771a-75ac-42f1-a401-04ca838b199f
@bind l4 overlay_mask_bind(ring_mask_M_LD)

# ╔═╡ a0f644b2-5e8b-4ccd-a866-1a54e992061e
overlay_mask_plot(arr, ring_mask_M_LD, l4, "ring mask")

# ╔═╡ 04d0cc19-9ae8-43b2-812c-b1b4532526bf
md"""
### Calculations
"""

# ╔═╡ 710a6131-f5e9-4dd4-811d-c816af39987d
begin
	single_ring_mask_M_LD = Bool.(ring_mask_M_LD[:, :, 3])
	s_bkg_M_LD = mean(single_arr[single_ring_mask_M_LD])
end

# ╔═╡ 7ae39137-dd69-4083-85a2-3418e7ee3a00
begin
	alg_M_LD = Integrated(arr[mask_M_LD_3D])
	mass_m_ld = score(s_bkg_M_LD, S_Obj_LD, pixel_size, ρ_LD, alg_M_LD)
end

# ╔═╡ d4cddcea-bdac-4ff8-8ec3-d628d4c897b5
md"""
# Score Small Inserts
"""

# ╔═╡ 3003a11d-2fb1-4448-b5ef-aec6e75bb668
md"""
## High Density
"""

# ╔═╡ 816611f8-7d13-4ecf-9590-8462e4043e37
begin
	mask_S_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_HD_3D[:, :, z] = mask_S_HD
	end
end;

# ╔═╡ b5d99cf0-6179-4dfd-a95a-33c16fb93853
md"""
#### Dilated mask
"""

# ╔═╡ 0314504d-75e0-4166-af40-ca1a072f50ec
dilated_mask_S_HD = dilate(dilate(dilate(dilate(dilate((mask_S_HD_3D))))));

# ╔═╡ c758677b-585d-4fa9-a54a-04b4ae09d66a
@bind m2 overlay_mask_bind(dilated_mask_S_HD)

# ╔═╡ f361eec2-55c3-4dc7-bc35-1f0fe00b6432
overlay_mask_plot(arr, dilated_mask_S_HD, m2, "dilated mask")

# ╔═╡ b3422bcc-8887-4d2e-94ba-18252a683ca0
md"""
#### Ring (background) mask
"""

# ╔═╡ 3e67144b-d511-4ab2-b2bb-44840ac5137d
ring_mask_S_HD = dilate(dilate(dilate(dilate(dilate(mask_S_HD_3D))))) - dilate(dilate(dilate(dilate(mask_S_HD_3D))));

# ╔═╡ 02b55afe-4beb-44a6-8e25-aa595b5ae37d
@bind m4 overlay_mask_bind(ring_mask_S_HD)

# ╔═╡ ab238843-177c-4ee6-9115-ab2c22d45e4f
overlay_mask_plot(arr, ring_mask_S_HD, m4, "ring mask")

# ╔═╡ 5d21acbc-7de1-4a7c-8b58-9d9a639fc47c
md"""
### Calculations
"""

# ╔═╡ 2ccf896b-493d-4343-b8c9-511cf6070512
begin
	single_ring_mask_S_HD = Bool.(ring_mask_S_HD[:, :, 3])
	s_bkg_S_HD = mean(single_arr[single_ring_mask_S_HD])
end

# ╔═╡ 258a0a68-bb63-46a2-873a-9664bdce2a60
begin
	alg_S_HD = Integrated(arr[mask_S_HD_3D])
	mass_s_hd = score(s_bkg_S_HD, S_Obj_HD, pixel_size, ρ_hd, alg_S_HD)
	if mass_s_hd < 0
		mass_s_hd = 0
	end
	mass_s_hd
end

# ╔═╡ bf3ca79a-cddc-41ab-8be5-8c7a6d63f5fe
md"""
## Medium Density
"""

# ╔═╡ 686ea6ee-1b36-4ec2-989b-d5dd3bf5030b
begin
	mask_S_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_MD_3D[:, :, z] = mask_S_MD
	end
end;

# ╔═╡ e86618e2-40a1-4184-9d88-df288ed75d64
md"""
#### Dilated mask
"""

# ╔═╡ ae7fd1e1-70f1-44a7-80c7-ceb7f57705bf
dilated_mask_S_MD = dilate(dilate(dilate(dilate(dilate(mask_S_MD_3D)))));

# ╔═╡ 29dd6e86-5556-4a34-b913-2848790701e2
@bind n2 overlay_mask_bind(dilated_mask_S_MD)

# ╔═╡ 9725746f-9575-43d6-abd2-de8500d07784
overlay_mask_plot(arr, dilated_mask_S_MD, n2, "dilated mask")

# ╔═╡ cdba29c4-4893-4558-97cc-b727ca0df5fd
md"""
#### Ring (background) mask
"""

# ╔═╡ 54804999-657c-4a5a-8878-e36b9bb8d3c5
ring_mask_S_MD = dilate(dilate(dilate(dilate(dilate(mask_S_MD_3D))))) - dilate(dilate(dilate(dilate(mask_S_MD_3D))));

# ╔═╡ b776cf19-98b3-4a33-94ee-9cac1519c3b1
@bind n4 overlay_mask_bind(ring_mask_S_MD)

# ╔═╡ 86ecfdc8-2ded-45bc-8664-20ec66d46317
overlay_mask_plot(arr, ring_mask_S_MD, n4, "ring mask")

# ╔═╡ ffacd4d0-6650-414c-b4b8-d9329dcbdd00
md"""
### Calculations
"""

# ╔═╡ 8251e4e4-a504-4a50-9820-312b5e4d684b
begin
	single_ring_mask_S_MD = Bool.(ring_mask_S_MD[:, :, 3])
	s_bkg_S_MD = mean(single_arr[single_ring_mask_S_MD])
end

# ╔═╡ e8cccd5a-700e-4aa0-9574-872fc17083aa
begin
	alg_S_MD = Integrated(arr[mask_S_MD_3D])
	mass_s_md = score(s_bkg_S_MD, S_Obj_MD, pixel_size, ρ_md, alg_S_MD)
	if mass_s_md < 0
		mass_s_md = 0
	end
	mass_s_md
end

# ╔═╡ 7662b1cb-b735-4bd6-9d62-2fc22becb6c0
md"""
## Low Density
"""

# ╔═╡ 7ef55870-a9fc-4770-ae46-46dbe599fa8a
begin
	mask_S_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_LD_3D[:, :, z] = mask_S_LD
	end
end;

# ╔═╡ df17f71d-df4d-402b-8f65-979b76e86c36
md"""
#### Dilated mask
"""

# ╔═╡ e4b3fa55-9b0f-4077-8d49-cf76f9ae4ca9
dilated_mask_S_LD = dilate(dilate(dilate(dilate(dilate(mask_S_LD_3D)))));

# ╔═╡ a4136d62-8077-4629-92fd-0c92d650b700
@bind o2 overlay_mask_bind(dilated_mask_S_LD)

# ╔═╡ a3e6a26c-90c4-4686-97b0-c95f7dc2049e
overlay_mask_plot(arr, dilated_mask_S_LD, o2, "dilated mask")

# ╔═╡ 8a26d2c4-af51-433f-abc5-581fbfa93e3a
md"""
#### Ring (background) mask
"""

# ╔═╡ c1c3531d-1452-482f-837c-878fb895b4d8
ring_mask_S_LD = dilate(dilate(dilate(dilate(dilate(mask_S_LD_3D))))) - dilate(dilate(dilate(dilate(mask_S_LD_3D))));

# ╔═╡ 5b211a9b-1df5-4d0b-8faa-b64860ba3ad9
@bind o4 overlay_mask_bind(ring_mask_S_LD)

# ╔═╡ 1ed50733-7a24-49e4-88ae-6010d2e88c95
overlay_mask_plot(arr, ring_mask_S_LD, o4, "ring mask")

# ╔═╡ 81d4ab94-8541-492e-b0d6-b6d5c5877094
md"""
### Calculations
"""

# ╔═╡ 265156ba-50ae-4b0f-b7fb-e1437cf82c25
begin
	single_ring_mask_S_LD = Bool.(ring_mask_S_LD[:, :, 3])
	s_bkg_S_LD = mean(single_arr[single_ring_mask_S_LD])
end

# ╔═╡ bed1604a-cc5f-4ed1-8aa0-f4a394e56ff4
begin
	alg_S_LD = Integrated(arr[mask_S_LD_3D])
	mass_s_ld = score(s_bkg_S_LD, S_Obj_LD, pixel_size, ρ_LD, alg_S_LD)
	if mass_s_ld < 0
		mass_s_ld = 0
	end
	mass_s_ld
end

# ╔═╡ ccbc4e73-6d1e-46d6-b81d-6c5c3dd0bc26
md"""
# Results
"""

# ╔═╡ 5029bc99-2350-4ede-a7f2-8f008ae82212
Phantoms.get_pixel_size(header)

# ╔═╡ 4677af95-37d4-4a8d-9cfe-d8b3a7b994de
inserts = [
	"Low Density",
	"Medium Density",
	"High Density"
]

# ╔═╡ d969236a-c10b-4dc0-875b-51e6eb96b5a0
volume_gt = [
	7.065,
	63.585,
	176.625
]

# ╔═╡ 6ef6c38a-7f5b-4241-952e-2151f5eba962
ground_truth_mass_large = [
	volume_gt[3] * density_array[2] * 1e-3,
	volume_gt[3] * density_array[3] * 1e-3,
	volume_gt[3] * density_array[4] * 1e-3
] # mg

# ╔═╡ 33f61073-7a3a-469e-b9ea-c682fb2e0d9f
calculated_mass_large = [
	mass_l_ld,
	mass_l_md,
	mass_l_hd
]

# ╔═╡ 1feed0de-862d-4291-b5a2-be9c990d6d0d
ground_truth_mass_medium = [
	volume_gt[2] * density_array[2] * 1e-3,
	volume_gt[2] * density_array[3] * 1e-3,
	volume_gt[2] * density_array[4] * 1e-3
]

# ╔═╡ 8b53f99f-3ae0-4e12-8e2e-34d8e7e0b8ac
calculated_mass_medium = [
	mass_m_ld,
	mass_m_md,
	mass_m_hd
]

# ╔═╡ 5815285e-de49-479a-a499-675e9f788a61
ground_truth_mass_small = [
	volume_gt[1] * density_array[2] * 1e-3,
	volume_gt[1] * density_array[3] * 1e-3,
	volume_gt[1] * density_array[4] * 1e-3
]

# ╔═╡ 59ce3f8c-e56a-4488-a46b-8c597db1cdc4
calculated_mass_small = [
	mass_s_ld,
	mass_s_md,
	mass_s_hd
]

# ╔═╡ c354f281-d4ff-4605-90d1-4859f9ab3205
df = DataFrame(
	inserts = inserts,
	ground_truth_mass_large = ground_truth_mass_large,
	calculated_mass_large = calculated_mass_large,
	ground_truth_mass_medium = ground_truth_mass_medium,
	calculated_mass_medium = calculated_mass_medium,
	ground_truth_mass_small = ground_truth_mass_small,
	calculated_mass_small = calculated_mass_small
)

# ╔═╡ cc147016-f73a-40ce-ac2b-0a0ef5dfe312
begin
	fmass2 = Figure()
	axmass2 = Axis(fmass2[1, 1])
	
	scatter!(density_array[2:end], df[!, :ground_truth_mass_large], label="ground_truth_mass_large")
	scatter!(density_array[2:end], df[!, :calculated_mass_large], label="calculated_mass_large")
	
	axmass2.title = "Mass Measurements (Large)"
	axmass2.ylabel = "Mass (mg)"
	axmass2.xlabel = "Density (mg/cm^3)"

	xlims!(axmass2, 0, 850)
	ylims!(axmass2, 0, 200)
	
	fmass2[1, 2] = Legend(fmass2, axmass2, framevisible = false)
	
	fmass2
end

# ╔═╡ 641f48d4-96a8-4d71-bfaf-cfd41af78548
begin
	fmass3 = Figure()
	axmass3 = Axis(fmass3[1, 1])
	
	scatter!(density_array[2:end], df[!, :ground_truth_mass_medium], label="ground_truth_mass_medium")
	scatter!(density_array[2:end], df[!, :calculated_mass_medium], label="calculated_mass_medium")
	
	axmass3.title = "Mass Measurements (Medium)"
	axmass3.ylabel = "Mass (mg)"
	axmass3.xlabel = "Density (mg/cm^3)"

	xlims!(axmass3, 0, 850)
	ylims!(axmass3, 0, 70)
	
	fmass3[1, 2] = Legend(fmass3, axmass3, framevisible = false)
	
	fmass3
end

# ╔═╡ b283112c-dc6e-4cc7-9a66-78ff7861dbda
begin
	fmass4 = Figure()
	axmass4 = Axis(fmass4[1, 1])
	
	scatter!(density_array[2:end], df[!, :ground_truth_mass_small], label="ground_truth_mass_small")
	scatter!(density_array[2:end], df[!, :calculated_mass_small], label="calculated_mass_small")
	
	axmass4.title = "Mass Measurements (Small)"
	axmass4.ylabel = "Mass (mg)"
	axmass4.xlabel = "Density (mg/cm^3)"

	xlims!(axmass4, 0, 850)
	ylims!(axmass4, 0, 10)
	
	fmass4[1, 2] = Legend(fmass4, axmass4, framevisible = false)
	
	fmass4
end

# ╔═╡ af7c0cf1-1814-4ef1-8ee7-a0ec5778370f
percent_error_large = (abs.(ground_truth_mass_large - calculated_mass_large) ./ ground_truth_mass_large) .* 100

# ╔═╡ c709f97f-b593-4ed2-8a62-baaccdb6aacf
percent_error_medium = (abs.(ground_truth_mass_medium - calculated_mass_medium) ./ ground_truth_mass_medium) .* 100

# ╔═╡ 2f8b945c-a825-4bdc-b492-3da1e6d1800b
percent_error_small= (abs.(ground_truth_mass_small - calculated_mass_small) ./ ground_truth_mass_small) .* 100

# ╔═╡ 1c58f834-b95f-456b-b87f-2c7d2cb00dd3
md"""
### Save Results
"""

# ╔═╡ f0974424-139b-4db8-8a85-56716a193ad9
if ~isdir(string(cd(pwd, "..") , "/output/", VENDER))
	mkdir(string(cd(pwd, "..") , "/output/", VENDER))
end

# ╔═╡ 68135e69-97d9-455b-b4e7-f7b0287d0629
output_path = string(cd(pwd, "..") , "/output/", VENDER, "/", scan, ".csv")

# ╔═╡ 2ce24b1c-f9a5-4d0c-a18b-069a317592e7
CSV.write(output_path, df)

# ╔═╡ Cell order:
# ╠═11bc4d31-8569-4ddf-9833-017f6b4cca66
# ╠═201127f2-ae5e-4542-9d9f-c134e16d2c29
# ╟─bfda24da-a8cf-4833-893e-2e039d9bcaf1
# ╠═bdd77f47-b5ff-4855-b3c4-1fcd5fe31777
# ╟─773eff82-436b-4b36-8203-0764a669ce62
# ╠═377b58ac-28ff-4bfe-8ce2-c21481e8e973
# ╠═e2fffe13-12b0-4645-aecd-b2e2359eb9d6
# ╠═a023991a-6814-4375-a0b7-8b4c70c7f2b4
# ╠═5b02e50f-0f8d-4987-aaf8-08cadf89562d
# ╠═ebb232b1-a4ba-441f-a79f-23d8489aff2b
# ╟─604fd08b-6677-428f-a263-84bf9c2eedac
# ╟─041a5099-8278-428f-9d00-db6295f32897
# ╟─c7868bf0-611f-4d02-8247-b87edc82313b
# ╟─9ca9a957-17e3-4b9a-9f3b-e1d042ad0bcb
# ╟─8723bd5f-6526-43b3-9de0-95c42fd30bbb
# ╠═07f993c4-43fd-4960-934f-abba867e8231
# ╟─de654d7a-ecba-4526-afa3-417cc5820b98
# ╠═3aec1071-d863-46b2-a66e-3d8500cf383c
# ╠═b01b496c-b80b-424f-a11f-d8bd5205ebbd
# ╠═bb1df3e9-793b-4727-b806-5962c38c21c3
# ╠═38c984b0-f132-42c6-bd99-27802bc5d191
# ╟─e9e2cc57-61bc-479f-bddd-fd5f930ffa65
# ╠═c19dff18-eadb-43ee-bfb4-4fa3e573d1e4
# ╟─c8f1c6d7-3717-4b14-9bb7-d72a37849ea8
# ╠═301e3bb2-3302-4251-931c-e32211136eea
# ╟─e60c281b-46cd-4a94-90a1-de0fdbb05491
# ╠═73ede908-a805-4d6b-93ae-3fc724b56423
# ╠═0e39214d-99c2-4af3-8f8e-78ab87f7fef6
# ╠═64d526d4-b5a7-44cd-9db0-7b5265b37857
# ╠═bf16de9b-3e2b-42ca-be6c-a40f2c2b87be
# ╟─27cc2663-2080-48f6-a3a0-7fae64868ab7
# ╠═0a85351b-74a7-47ec-8baa-fe9f2b6e81f9
# ╠═ac5895f8-e629-4500-a895-b781ed76a203
# ╠═5f8bc4e3-3bd3-456a-8b0f-55bd5912d81c
# ╠═d1c93fea-1ced-4a68-b66a-9d8978ba7a2a
# ╠═6b357f4c-73a7-44ab-a554-76dfa0ebfb34
# ╠═5ed7dba0-b78b-4484-94b3-0b34ae94a305
# ╠═f1355c86-cdc4-4f76-967a-5a6285ecbac6
# ╠═4cc8d210-44fe-4524-aeb7-0688a383989d
# ╠═9ce3d6a3-dbb5-480e-a34f-5084b5572d7c
# ╠═9368bf8f-54de-4619-86c5-cdd01e44a555
# ╠═e20edcfd-1d64-44b5-9b08-27a77f44299f
# ╠═e1d6755a-bb5f-44a9-abc8-240beac0f5ee
# ╟─66582769-d393-48e6-826e-98337ebbe219
# ╠═035fd1c1-dadf-438f-b48b-7f6e3313e9dc
# ╠═283d2191-06cd-40f3-bab0-d66f1ded8d27
# ╠═b4d32db5-b5be-4fb4-9666-2bfe0aa75d65
# ╠═7803d8e7-a3fe-4043-af64-4e9b334c8f82
# ╠═ee1bfe95-106c-4e43-846e-6293a1413830
# ╠═22540057-0179-49bb-8437-092839d4e7b5
# ╠═8c3a1111-c5c6-4bd4-ad14-e1c578471354
# ╟─6a5c9b4d-c8ad-44e0-bcc5-e4509d266ba5
# ╠═7b76786b-3cc7-447d-b78f-007d6fdc16dc
# ╠═87fa7e66-945c-442a-8cd3-e4e302379074
# ╠═1ef692d4-91e0-44cc-a3bc-0002f145ca80
# ╟─ca9554e3-147c-403a-8027-c286c1cdcd5f
# ╠═2f9e82ff-ede9-4d0b-997e-2c814d997527
# ╠═6220799c-323e-49b3-83be-f828d75dfc4d
# ╟─3bedd4e1-696a-4a40-a18d-2e2575402399
# ╠═d730dac7-ae9c-4102-ab32-d24e4d0c7709
# ╟─b1679304-f66b-4199-9911-d9d02a757b62
# ╠═212554d6-71a2-47fd-8002-8edb87e1dce5
# ╟─c4098846-e781-4a52-b280-f4a313c2cd2a
# ╠═e7efea19-74d6-4513-be91-f010c81b9fd2
# ╟─143c8c29-0063-4868-bf4e-8fd3c5c08dfb
# ╠═e77a2392-1fbc-4fb0-8bef-d7a9283782be
# ╟─6784b04f-ece0-447c-8820-a0644eaf823e
# ╠═876eae18-f388-4f9b-bff5-4ec58657d3a0
# ╟─889052c7-d431-49f8-aa71-bd2c59aed37f
# ╠═ada32081-4eaf-4189-801a-f91d9c2f189a
# ╠═0499d5a5-22a2-48b4-afc3-59a03eed10df
# ╠═e5cdb203-85ba-4217-8b6c-4a9a782d54ad
# ╠═60ccc705-6931-40b8-8f46-955398fe13c9
# ╟─696a45ad-aa73-4daa-b8f6-81584e1cfb9e
# ╠═c3d0bb9e-99dc-49d2-b30d-a2534e300a0f
# ╟─aac2fcd5-9ab1-42da-8464-76e1a3fbe443
# ╠═b4f9c38e-dd42-49d2-9257-3d2c13653878
# ╟─53dcc1a6-1cc1-4d24-bd13-d84d9335b1fa
# ╠═b720e961-a907-45db-b86b-3330a809da4f
# ╟─75d2f8af-13e0-425d-a40b-c8cffb371e7b
# ╠═edaceb23-784f-4149-8b74-df3089208a62
# ╠═1d3bab04-b000-4d71-a729-67620a4b4e73
# ╠═01333998-bf78-4b97-b45f-6a7b1ef22ea3
# ╟─533aa777-c1cc-4067-b6de-2019bb11ee77
# ╠═fe50b0e4-af31-441b-a513-cfc914584cb2
# ╠═4c24aa2b-4c92-435a-9bac-1240741126fc
# ╠═4fd5d0ae-87b3-41a9-a4a8-7b4687592e2e
# ╟─73047cbf-17df-430b-8a8f-69af936a89af
# ╠═a5767c29-300c-40d7-aeff-9f373f211b7a
# ╟─64cd83f3-f9ed-4215-8bd2-be3453bcc4dc
# ╠═607aed18-2914-4b48-a2cf-e2f4e78524d2
# ╟─f0fdd9ae-8b61-47b5-8826-923d1779bc8b
# ╠═e1c1f280-a844-4ccd-8237-0762b0ffc7c4
# ╟─44e8196c-6aa5-4797-b7a2-289325633c84
# ╠═422fdee4-6685-46e3-961b-64425f205ec8
# ╟─e5e751bc-cee3-4593-9978-68578bc4021c
# ╠═161a8f24-e102-4ceb-b40d-77756139f089
# ╟─9b089a9d-6037-4bed-bf11-0c494aef3581
# ╠═b834eae9-a7c8-4ad3-8161-2d3463fa99d7
# ╠═f6ca11d6-ff98-49f6-98e4-046eaa77931f
# ╠═aec7f59f-b0d0-4a6e-96b2-94bc2192b2d2
# ╟─1a998d70-1867-4571-b534-5ef8b061e60d
# ╟─653c8222-07f7-4633-9bf1-b7b02e67b3c3
# ╠═4e5151cb-7acc-4331-b011-a4cd17f84f78
# ╟─6ed7f390-de7e-4489-bc82-d4a214baca2d
# ╠═6d092d23-46e0-4a47-a945-009a48471256
# ╟─4caade91-81bf-4bfe-bbe0-cb428c57e735
# ╠═84b41be9-f685-4a36-9b18-ed25a5c283d9
# ╟─8a5f08e3-0ef3-4b87-93b1-649328e07d96
# ╠═3a930ac5-8cdf-4677-9303-6e13eef7bbef
# ╟─7651bc64-161f-46e5-8947-a58840cae9c0
# ╠═d80c7da7-3fb0-48ad-ac51-84f6005e5276
# ╟─ded2cd96-0bb2-4718-8503-ca3361cdf168
# ╠═7833403f-c17f-4346-bcd5-c4f66664e348
# ╠═ee3bfc69-15d0-4695-a1e3-3a72060f91de
# ╟─4214559c-8c4e-41f9-a003-f926bedb1dda
# ╠═5044a80d-fa5d-413f-98fc-ba8208280803
# ╟─fc119ac1-f6ad-4ef5-b641-a7961938b4cd
# ╠═843bbb31-f437-453c-ba61-864fa4234257
# ╠═7ff5e49e-6f18-48d0-b97f-1abd9586867a
# ╠═aa70689c-adb7-45d8-935c-40ea9b65901e
# ╟─b2974709-0f5f-4c14-b902-6021a624dc7c
# ╠═95d5c859-daf3-49b3-87e9-8b45b226d99b
# ╠═1ccae49d-9d4d-49ed-baeb-232eb31f1aae
# ╠═790998a8-542c-46f9-a850-0d07034bf22c
# ╠═b4cc3468-9173-4168-8149-bd29127be392
# ╠═1d16adec-06f5-4816-a277-cadef49ca558
# ╠═dcf74ca2-cf93-4f6d-a53c-b47b8f3435a3
# ╟─08837a72-b43a-4b36-80ba-92cd64e59c79
# ╠═67bd5d43-7d7e-4843-9288-61e67b22a5cf
# ╟─d01bfc39-f257-482c-957f-f86ea7b5fbcc
# ╠═6b2c2f08-cdfd-45d4-ab2f-f23a5e1bb1a6
# ╠═52eb3999-6403-43df-96de-0b2fa3930b1f
# ╠═73ea58c2-77cc-451f-921b-d6db7c9ce56f
# ╟─35e5cc8e-bcc5-4365-a9ce-f75591af4520
# ╠═31adbea2-ccc3-4fa7-ae32-bfb92e8d0fbe
# ╟─37cd771a-75ac-42f1-a401-04ca838b199f
# ╠═a0f644b2-5e8b-4ccd-a866-1a54e992061e
# ╟─04d0cc19-9ae8-43b2-812c-b1b4532526bf
# ╠═710a6131-f5e9-4dd4-811d-c816af39987d
# ╠═7ae39137-dd69-4083-85a2-3418e7ee3a00
# ╟─d4cddcea-bdac-4ff8-8ec3-d628d4c897b5
# ╟─3003a11d-2fb1-4448-b5ef-aec6e75bb668
# ╠═816611f8-7d13-4ecf-9590-8462e4043e37
# ╟─b5d99cf0-6179-4dfd-a95a-33c16fb93853
# ╠═0314504d-75e0-4166-af40-ca1a072f50ec
# ╟─c758677b-585d-4fa9-a54a-04b4ae09d66a
# ╠═f361eec2-55c3-4dc7-bc35-1f0fe00b6432
# ╟─b3422bcc-8887-4d2e-94ba-18252a683ca0
# ╠═3e67144b-d511-4ab2-b2bb-44840ac5137d
# ╠═02b55afe-4beb-44a6-8e25-aa595b5ae37d
# ╠═ab238843-177c-4ee6-9115-ab2c22d45e4f
# ╠═5d21acbc-7de1-4a7c-8b58-9d9a639fc47c
# ╠═2ccf896b-493d-4343-b8c9-511cf6070512
# ╠═258a0a68-bb63-46a2-873a-9664bdce2a60
# ╟─bf3ca79a-cddc-41ab-8be5-8c7a6d63f5fe
# ╠═686ea6ee-1b36-4ec2-989b-d5dd3bf5030b
# ╟─e86618e2-40a1-4184-9d88-df288ed75d64
# ╠═ae7fd1e1-70f1-44a7-80c7-ceb7f57705bf
# ╟─29dd6e86-5556-4a34-b913-2848790701e2
# ╠═9725746f-9575-43d6-abd2-de8500d07784
# ╠═cdba29c4-4893-4558-97cc-b727ca0df5fd
# ╠═54804999-657c-4a5a-8878-e36b9bb8d3c5
# ╠═b776cf19-98b3-4a33-94ee-9cac1519c3b1
# ╠═86ecfdc8-2ded-45bc-8664-20ec66d46317
# ╠═ffacd4d0-6650-414c-b4b8-d9329dcbdd00
# ╠═8251e4e4-a504-4a50-9820-312b5e4d684b
# ╠═e8cccd5a-700e-4aa0-9574-872fc17083aa
# ╟─7662b1cb-b735-4bd6-9d62-2fc22becb6c0
# ╠═7ef55870-a9fc-4770-ae46-46dbe599fa8a
# ╟─df17f71d-df4d-402b-8f65-979b76e86c36
# ╠═e4b3fa55-9b0f-4077-8d49-cf76f9ae4ca9
# ╟─a4136d62-8077-4629-92fd-0c92d650b700
# ╠═a3e6a26c-90c4-4686-97b0-c95f7dc2049e
# ╟─8a26d2c4-af51-433f-abc5-581fbfa93e3a
# ╠═c1c3531d-1452-482f-837c-878fb895b4d8
# ╟─5b211a9b-1df5-4d0b-8faa-b64860ba3ad9
# ╠═1ed50733-7a24-49e4-88ae-6010d2e88c95
# ╟─81d4ab94-8541-492e-b0d6-b6d5c5877094
# ╠═265156ba-50ae-4b0f-b7fb-e1437cf82c25
# ╠═bed1604a-cc5f-4ed1-8aa0-f4a394e56ff4
# ╟─ccbc4e73-6d1e-46d6-b81d-6c5c3dd0bc26
# ╠═5029bc99-2350-4ede-a7f2-8f008ae82212
# ╠═4677af95-37d4-4a8d-9cfe-d8b3a7b994de
# ╠═d969236a-c10b-4dc0-875b-51e6eb96b5a0
# ╠═6ef6c38a-7f5b-4241-952e-2151f5eba962
# ╠═33f61073-7a3a-469e-b9ea-c682fb2e0d9f
# ╠═1feed0de-862d-4291-b5a2-be9c990d6d0d
# ╠═8b53f99f-3ae0-4e12-8e2e-34d8e7e0b8ac
# ╠═5815285e-de49-479a-a499-675e9f788a61
# ╠═59ce3f8c-e56a-4488-a46b-8c597db1cdc4
# ╠═c354f281-d4ff-4605-90d1-4859f9ab3205
# ╟─cc147016-f73a-40ce-ac2b-0a0ef5dfe312
# ╟─641f48d4-96a8-4d71-bfaf-cfd41af78548
# ╟─b283112c-dc6e-4cc7-9a66-78ff7861dbda
# ╠═af7c0cf1-1814-4ef1-8ee7-a0ec5778370f
# ╠═c709f97f-b593-4ed2-8a62-baaccdb6aacf
# ╠═2f8b945c-a825-4bdc-b492-3da1e6d1800b
# ╟─1c58f834-b95f-456b-b87f-2c7d2cb00dd3
# ╠═f0974424-139b-4db8-8a85-56716a193ad9
# ╠═68135e69-97d9-455b-b4e7-f7b0287d0629
# ╠═2ce24b1c-f9a5-4d0c-a18b-069a317592e7

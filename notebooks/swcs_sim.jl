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

# ╔═╡ 2c064fd4-3719-4503-baa3-61a453e3642d
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
end

# ╔═╡ 355b6a05-065e-4b2c-b65f-af9c324756d0
TableOfContents()

# ╔═╡ 739ccb4d-1378-4241-9e15-fffaa4522080
md"""
## Load DICOMS

All you need to do is set `base_path` once and leave it. After that, the only thing that should change is the `VENDER`, once for every set, and the `SCAN_NUMBER`, once for each scan.
"""

# ╔═╡ 29290cd1-945c-4ec2-9cc1-fe65c7d21dae
begin
	SCAN_NUMBER = 1
	VENDER = "80"
	kern = 0
	TYPE = "swcs"
	BASE_PATH = "/Users/daleblack/Google Drive/Datasets/Simulated/"
end

# ╔═╡ a5b8ce68-4ad4-499c-8caf-c0c477e3d53d
md"""
**Everything below should be automatic, just scroll through to visually inspect that things make sense**
"""

# ╔═╡ 6aaf8df0-a546-432c-ad17-d8d735fa0919
root_path = string(BASE_PATH, VENDER)

# ╔═╡ dd590806-ee35-47a8-898c-d158d6133eee
dcm_path_list = dcm_list_builder(root_path)

# ╔═╡ 2ab91f08-fac4-48d6-9776-ac47002159ab
pth = dcm_path_list[SCAN_NUMBER]

# ╔═╡ ce9c455d-ad21-4c27-8468-84be6f14f5ed
scan = basename(pth)

# ╔═╡ 56c069d9-d1d6-4f6f-b2d1-cef68a970d9b
# header, dcm_array, slice_thick_ori1 = dcm_reader(pth);

# ╔═╡ 621eb1c3-a294-4620-b79c-1c98e501a21f
begin
	header, dcm_array, slice_thick_ori1 = dcm_reader(pth)
	if kern != 0
		for z in size(dcm_array, 3)
			dcm_array[:, :, z] = mult_gauss(dcm_array[:, :, z], kern)
		end
	end
end;

# ╔═╡ d4cea938-0455-4f3d-b9f9-9a5d6c9f6967
md"""
## Helper Functions
"""

# ╔═╡ 80fff6f2-c369-4efc-84a2-da15e7b718f7
function collect_tuple(tuple_array)
	row_num = size(tuple_array)
	col_num = length(tuple_array[1])
	container = zeros(Int64, row_num..., col_num)
	for i in 1:length(tuple_array)
		container[i,:] = collect(tuple_array[i])
	end
	return container
end

# ╔═╡ 4f9e6c87-1d73-441d-a0e6-b0c2934dd035
function overlay_mask_bind(mask)
	indices = findall(x -> x == 1, mask)
	indices = Tuple.(indices)
	label_array = collect_tuple(indices)
	zs = unique(label_array[:,3])
	return PlutoUI.Slider(1:length(zs), default=3, show_value=true)
end

# ╔═╡ 2cdc2b8c-7b26-421e-a008-62cc56c026fd
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

# ╔═╡ 4256495e-319a-4587-89a5-0185a2cb22e4
md"""
## Segment Heart
"""

# ╔═╡ 6f5312d8-0148-4b5b-ab62-8a02e1e66b5a
masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3)÷2);

# ╔═╡ 14c64d65-06d9-4919-b8fd-fc5d7d5eac63
@bind a PlutoUI.Slider(1:size(masked_array, 3), default=10, show_value=true)

# ╔═╡ 91eef4a3-80df-4764-97d0-1d23cb34de7b
heatmap(masked_array[:, :, a], colormap=:grays)

# ╔═╡ 9df0b48d-e053-4d40-bad8-c12c0813beb8
begin
	fig = Figure()
	
	ax = Makie.Axis(fig[1, 1])
	ax.title = "Raw DICOM Array"
	heatmap!(transpose(dcm_array[:, :, 5]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig
end

# ╔═╡ 080e44cd-7918-4a79-bfa5-720830738388
begin
	fig2 = Figure()
	
	ax2 = Makie.Axis(fig2[1, 1])
	ax2.title = "Mask Array"
	heatmap!(transpose(mask), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig2
end

# ╔═╡ 312dbc4e-2f66-48da-8801-92d6f9e9301f
begin
	fig3 = Figure()
	
	ax3 = Makie.Axis(fig3[1, 1])
	ax3.title = "Masked DICOM Array"
	heatmap!(transpose(masked_array[:, :, 1]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig3
end

# ╔═╡ b6ae6e8f-b1b4-492c-977f-aa0e72d76119
md"""
## Segment Calcium Rod
"""

# ╔═╡ 773573d7-97c6-4b93-9e3d-ce544918bb4f
calcium_image, slice_CCI, quality_slice, cal_rod_slice = mask_rod(masked_array, header);

# ╔═╡ e7c9f252-ca91-4a7a-84c3-7bcd647af9b2
@bind c PlutoUI.Slider(1:size(calcium_image, 3), default=cal_rod_slice, show_value=true)

# ╔═╡ adfe2070-3ae3-412f-b7f3-99622d6003c5
heatmap(transpose(calcium_image[:, :, c]), colormap=:grays)

# ╔═╡ e14c4255-b10d-4502-8b7c-94a5daed1b09
md"""
## Segment Calcium Inserts
"""

# ╔═╡ 58ff4866-983d-4481-9442-d9085e97f855
mask_L_HD, mask_M_HD, mask_S_HD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_LD, mask_M_LD, mask_S_LD = mask_inserts_simulation(
            dcm_array, masked_array, header, slice_CCI, center_insert
);

# ╔═╡ 47982924-7917-41f0-989b-3233c5e5f402
slice_CCI

# ╔═╡ 1d800cab-b66a-4497-9e7c-1367bfdd3a87
masks = mask_L_HD + mask_M_HD + mask_S_HD + mask_L_MD + mask_M_MD + mask_S_MD + mask_L_LD + mask_M_LD + mask_S_LD;

# ╔═╡ c4c6c326-9b15-478b-b462-f4f20b4601fb
heatmap(masks, colormap=:grays)

# ╔═╡ 2630ef4c-1870-4d83-b191-027a1db1bd47
md"""
## Calibration Prep
"""

# ╔═╡ b7f9e7cd-c626-4df5-833e-e4491baaaf87
array_filtered = abs.(mapwindow(median, calcium_image[:, :, 2], (3, 3)));

# ╔═╡ 742b3c55-5131-4554-8980-1284fa6d49a2
bool_arr = array_filtered .> 0;

# ╔═╡ fc069260-4b7f-4ee5-8bc3-29d481e0d6cc
bool_arr_erode = (((erode(erode(bool_arr)))));

# ╔═╡ 5ad725f9-2997-4af7-b324-1224fd0e3017
heatmap(bool_arr, colormap=:grays)

# ╔═╡ bed91830-f655-4b9f-845e-0353aecf7734
heatmap(bool_arr_erode, colormap=:grays)

# ╔═╡ 8ab94c76-47f2-48d0-888c-f4ea960cf276
c_img = calcium_image[:, :, 1:3];

# ╔═╡ d721be4b-5599-49bc-abd9-160175455e22
begin
	mask_cal_3D = Array{Bool}(undef, size(c_img))
	for z in 1:size(c_img, 3)
		mask_cal_3D[:, :, z] = bool_arr_erode
	end
end;

# ╔═╡ 0f40aa11-8ca6-4574-905b-eb6ee9772b90
hist(c_img[mask_cal_3D])

# ╔═╡ d815512f-3ef6-4481-8b12-33b1288edc87
cal_insert_mean = quantile!(c_img[mask_cal_3D], 0.7)

# ╔═╡ 65f63514-38ac-4662-b468-4883604276ee
md"""
### Calibration Line
"""

# ╔═╡ 15fe3d10-a572-405d-b1f6-09f3637d58f8
density_array = [0, 200, 400, 800]

# ╔═╡ 8fe16dbf-14aa-4f4a-8f22-e5bffac397d5
density_array_cal = [0, 200]

# ╔═╡ 76d9bbbe-db91-4bb6-b271-1aea6e37e9c6
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

# ╔═╡ 0cd096ca-2b26-4dda-9cb6-33ad2b296db0
md"""
# Score Large Inserts
"""

# ╔═╡ d070d007-72a0-460f-b1d3-c1798e25a7af
arr = masked_array[:, :, 4:6];

# ╔═╡ dca26385-2c02-4085-9cdf-2e401de044d9
single_arr = masked_array[:, :, slice_CCI];

# ╔═╡ bdbe21c9-87c1-45ea-bf8a-4cc7bd65491d
md"""
## High Density
"""

# ╔═╡ afc1d81f-1e02-47da-8773-bedcbf309166
begin
	mask_L_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_HD_3D[:, :, z] = mask_L_HD
	end
end;

# ╔═╡ 8aeca764-a3f1-482a-b9b0-975d5d08fb6a
begin
	eroded_mask_L_HD = erode(erode(mask_L_HD_3D))
	high_density_cal = mean(arr[eroded_mask_L_HD])
end

# ╔═╡ 92e85638-3f71-41db-8d1e-6fcf46ceaf69
md"""
#### Dilated mask
"""

# ╔═╡ 5e98dad3-5bf1-497f-93ea-69d3acdf3411
dilated_mask_L_HD = dilate(dilate(mask_L_HD_3D));

# ╔═╡ 8494fff0-91a3-413f-8a0b-cb9c79223fd0
@bind g2 overlay_mask_bind(dilated_mask_L_HD)

# ╔═╡ 53bf9d17-537b-4d70-a67c-5c88ceef08a9
overlay_mask_plot(arr, dilated_mask_L_HD, g2, "dilated mask")

# ╔═╡ a02e06fe-ff2e-4a46-93f5-e661fae20e4f
md"""
#### Ring (background) mask
"""

# ╔═╡ 9aee1a0a-2204-45bf-a6e1-b42be746d351
ring_mask_L_HD = dilate(dilate(dilate(dilate(mask_L_HD_3D)))) - dilate(dilate(dilate(mask_L_HD_3D)));

# ╔═╡ 8ea82534-97c5-4494-ad06-a7979f4ff8b3
@bind g4 overlay_mask_bind(ring_mask_L_HD)

# ╔═╡ 4a2a49a9-c435-4caa-92f2-74b5f92a782b
overlay_mask_plot(arr, ring_mask_L_HD, g4, "ring mask")

# ╔═╡ 118efcb2-f712-49f2-9589-7453e4e0e342
begin
	single_ring_mask_L_HD = Bool.(ring_mask_L_HD[:, :, 3])
	s_bkg_L_HD = mean(single_arr[single_ring_mask_L_HD])
end

# ╔═╡ 0b505164-0978-4c40-b044-b0211b2bab7c
pixel_size = DICOMUtils.get_pixel_size(header)

# ╔═╡ f4585dda-4493-4caf-93f0-36a924b50bea
md"""
## Medium Density
"""

# ╔═╡ 96323aa7-b24f-4942-9197-08cb49e96427
begin
	mask_L_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_MD_3D[:, :, z] = mask_L_MD
	end
end;

# ╔═╡ 5270b57b-5cb0-48b6-bfd8-ffac3463d9e6
begin
	eroded_mask_L_MD = erode(erode(mask_L_MD_3D))
	med_density_cal = mean(arr[eroded_mask_L_MD])
end

# ╔═╡ 4a045f83-bc81-42c9-a048-7c20855fbaec
md"""
#### Dilated mask
"""

# ╔═╡ b7e539c3-56a6-47e0-ae4a-8a15679331c2
dilated_mask_L_MD = dilate(dilate(mask_L_MD_3D));

# ╔═╡ e539b82b-945e-463e-a2c9-8240dbb7dd12
@bind h2 overlay_mask_bind(dilated_mask_L_MD)

# ╔═╡ e976610d-5c02-414e-b8c6-db989b5b098a
overlay_mask_plot(arr, dilated_mask_L_MD, h2, "dilated mask")

# ╔═╡ b0fb5549-e0ec-47a4-bb56-122a1f76c6a8
md"""
#### Ring (background) mask
"""

# ╔═╡ 029bd311-f8f8-4a06-bde3-24be773b6970
ring_mask_L_MD = dilate(dilate(dilate(dilate(mask_L_MD_3D)))) - dilate(dilate(dilate(mask_L_MD_3D)));

# ╔═╡ ac775f8c-77e7-4a60-a844-f57b9be48c93
@bind h4 overlay_mask_bind(ring_mask_L_MD)

# ╔═╡ 4bfdd6ee-2bfe-4e6f-8fd4-acdb90775cbb
overlay_mask_plot(arr, ring_mask_L_MD, h4, "ring mask")

# ╔═╡ 7894faf9-83c4-431f-a6b4-9baa1f97ffcb
begin
	single_ring_mask_L_MD = Bool.(ring_mask_L_MD[:, :, 3])
	s_bkg_L_MD = mean(single_arr[single_ring_mask_L_MD])
end

# ╔═╡ 2f3a1e98-1c33-47cc-b90a-b80fcf552356
md"""
## Low Density
"""

# ╔═╡ 7894ca8c-e6d8-402e-8f0d-99e936de513e
begin
	mask_L_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_LD_3D[:, :, z] = mask_L_LD
	end
end;

# ╔═╡ bc7f7c97-4fbd-44a5-9ee0-f6390b2e2674
begin
	eroded_mask_L_LD = erode(erode(mask_L_LD_3D))
	low_density_cal = mean(arr[eroded_mask_L_LD])
end

# ╔═╡ 0d2c6476-cf89-445f-ac9f-db54972aec02
intensity_array = [0, low_density_cal, med_density_cal, high_density_cal] # HU

# ╔═╡ 3dcb027b-87a7-411d-a51f-23b7606cd672
intensity_array_cal = [0, low_density_cal]

# ╔═╡ bcaafe33-b0dd-4ec6-96cf-ff8fe04d9ad5
df_cal = DataFrame(:density => density_array_cal, :intensity => intensity_array_cal)

# ╔═╡ 05c4ff2b-7d4f-45a4-9a8d-ce1afd429f4b
linearRegressor = lm(@formula(intensity ~ density), df_cal);

# ╔═╡ 7332bcca-cd2a-4809-84e1-19c0829b0055
linearFit = predict(linearRegressor)

# ╔═╡ 5a754fcc-3dd8-4ecb-b762-ae913a2c3ad5
m = linearRegressor.model.pp.beta0[2]

# ╔═╡ b6f7ed3a-20f1-4136-afc4-0c5ba9632343
b = linearRegressor.model.rr.mu[1]

# ╔═╡ 450f264d-cd87-476c-9036-c64ba9de3c63
density(intensity) = (intensity - b) / m

# ╔═╡ d4c5d18a-c7a1-48e1-876c-d4c4b94eab8f
intensity(ρ) = m*ρ + b

# ╔═╡ b06deeb6-7090-4d17-bdba-1d90dad6b202
S_Obj_HD = intensity(800)

# ╔═╡ 145e5acb-9e35-40c8-b326-fc71c4a07dc5
begin
	alg_L_HD = Integrated(arr[mask_L_HD_3D])
	ρ_hd = 0.8 # mg/mm^3
	mass_l_hd = score(s_bkg_L_HD, S_Obj_HD, pixel_size, ρ_hd, alg_L_HD)
end

# ╔═╡ 5a3aadde-5120-46a1-b861-072d27ce46b7
S_Obj_MD = intensity(400)

# ╔═╡ 95cac882-2c04-4728-8c69-46f263dbbd1d
begin
	alg_L_MD = Integrated(arr[mask_L_MD_3D])
	ρ_md = 0.4
	mass_l_md = score(s_bkg_L_MD, S_Obj_MD, pixel_size, ρ_md, alg_L_MD)
end

# ╔═╡ c3cc88b1-76b4-4a61-82bf-f01c403124a3
begin
	f = Figure()
	ax1 = Axis(f[1, 1])
	
	scatter!(density_array_cal, intensity_array_cal)
	lines!(density_array_cal, linearFit, color = :red)
	ax1.title = "Calibration Line (Intensity vs Density)"
	ax1.ylabel = "Intensity (HU)"
	ax1.xlabel = "Density (mg/cm^3)"
	
	f
end

# ╔═╡ c25633c7-12ee-4497-b733-fa6d785bd77c
md"""
#### Dilated mask
"""

# ╔═╡ d9b0b828-f910-4f10-ba56-b4e6819a9004
dilated_mask_L_LD = dilate(dilate(mask_L_LD_3D));

# ╔═╡ 2cc31f1c-8a04-4d00-8ec9-c21c6289a587
@bind i2 overlay_mask_bind(dilated_mask_L_LD)

# ╔═╡ 44145515-28f3-4f57-8c29-01681249c9dd
overlay_mask_plot(arr, dilated_mask_L_LD, i2, "dilated mask")

# ╔═╡ 5ad20c65-58c1-4fa5-9ce0-15336843895b
md"""
#### Ring (background) mask
"""

# ╔═╡ f496aef2-a7cc-4dc9-bcfc-ab67ef00e8a3
ring_mask_L_LD = dilate(dilate(dilate(dilate(mask_L_LD_3D)))) - dilate(dilate(dilate(mask_L_LD_3D)));

# ╔═╡ 14d13b95-252f-4bf1-9f39-94aa0a004959
@bind i4 overlay_mask_bind(ring_mask_L_LD)

# ╔═╡ f64803d5-610a-4e01-a38e-48d03eddd1da
overlay_mask_plot(arr, ring_mask_L_LD, i4, "ring mask")

# ╔═╡ fb586eda-bfbe-4427-bd85-476fec08ec09
begin	
	single_ring_mask_L_LD = Bool.(ring_mask_L_LD[:, :, 3])
	s_bkg_L_LD = mean(single_arr[single_ring_mask_L_LD])
end

# ╔═╡ b16ad6b0-12d8-4a28-8278-f286f7f0f11f
S_Obj_LD = intensity(200)

# ╔═╡ 14fa6374-1d4c-4194-9b7f-47a35232a328
begin
	alg_L_LD = Integrated(arr[mask_L_LD_3D])
	ρ_LD = 0.2
	mass_l_ld = score(s_bkg_L_LD, cal_insert_mean, pixel_size, ρ_LD, alg_L_LD)
end

# ╔═╡ 060a220d-b6be-4cf2-a8b5-ed3f99ce3ca8
md"""
# Score Medium Inserts
"""

# ╔═╡ 9ee074e8-5627-4ece-9150-057e4b3cf8a4
md"""
## High Density
"""

# ╔═╡ ef0fa44f-75fe-45bd-b61d-34c34dbee506
begin
	mask_M_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_HD_3D[:, :, z] = mask_M_HD
	end
end;

# ╔═╡ 13509cd8-557a-4883-8152-2ab1b29830cf
md"""
#### Dilated mask
"""

# ╔═╡ 1f8b288a-25bc-41d3-9d24-7b325543dca2
dilated_mask_M_HD = dilate(dilate(dilate(dilate(mask_M_HD_3D))));

# ╔═╡ 371e86b7-5a7e-44fd-b8e6-d2cf85cb4783
@bind j2 overlay_mask_bind(dilated_mask_M_HD)

# ╔═╡ e5bed8d6-e403-4ed3-b21e-2eb513a2c112
overlay_mask_plot(arr, dilated_mask_M_HD, j2, "dilated mask")

# ╔═╡ 52499aa5-93f2-4ecf-ad28-6454196de7b5
md"""
#### Ring (background) mask
"""

# ╔═╡ a3f66bd2-13d9-42db-8660-4a135e608d03
ring_mask_M_HD = dilate(dilate(dilate(dilate(dilate(mask_M_HD_3D))))) - dilate(dilate(dilate(dilate(mask_M_HD_3D))));

# ╔═╡ 58645445-5011-4139-97c0-1891ae651a5c
@bind j4 overlay_mask_bind(ring_mask_M_HD)

# ╔═╡ c635a53e-ff58-4270-a4a2-76c8358567d4
overlay_mask_plot(arr, ring_mask_M_HD, j4, "ring mask")

# ╔═╡ 0f5e6427-2480-4958-bbad-c54499288f61
begin
	single_ring_mask_M_HD = Bool.(ring_mask_M_HD[:, :, 3])
	s_bkg_M_HD = mean(single_arr[single_ring_mask_M_HD])
end

# ╔═╡ 839873c2-2927-44c7-86f6-cf729a42d6ea
begin
	alg_M_HD = Integrated(arr[mask_M_HD_3D])
	mass_m_hd = score(s_bkg_M_HD, S_Obj_HD, pixel_size, ρ_hd, alg_M_HD)
end

# ╔═╡ 63ba13f6-c17e-47ad-90e2-9c20a3c606b9
md"""
## Medium Density
"""

# ╔═╡ 8b868d1c-06b0-4169-bd52-f66c6409b360
begin
	mask_M_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_MD_3D[:, :, z] = mask_M_MD
	end
end;

# ╔═╡ ccae0b6e-bd74-43e0-ae77-bf26f4a12a1e
md"""
#### Dilated mask
"""

# ╔═╡ 8c34dc16-ce3e-4ddd-bab8-fdf8fb633596
dilated_mask_M_MD = dilate(dilate(dilate(dilate(mask_M_MD_3D))));

# ╔═╡ 30c313d4-5851-4ed4-9d0e-d64b36dcbacf
@bind k2 overlay_mask_bind(dilated_mask_M_MD)

# ╔═╡ 1511cbf4-cdf4-4bb0-a678-fd1c327b23d7
overlay_mask_plot(arr, dilated_mask_M_MD, k2, "dilated mask")

# ╔═╡ ff09ac24-388b-4f94-82ce-1c8f9ce22684
md"""
#### Ring (background) mask
"""

# ╔═╡ 35d35d55-397d-45eb-93c2-b846acef8553
ring_mask_M_MD = dilate(dilate(dilate(dilate(dilate(dilate(mask_M_MD_3D)))))) - dilate(dilate(dilate(dilate(dilate(mask_M_MD_3D)))));

# ╔═╡ f459d9c5-705a-4e54-a97e-a1f5710f4757
@bind k4 overlay_mask_bind(ring_mask_M_MD)

# ╔═╡ 3e4ab08d-77c7-4304-acec-13857ac5c4bb
overlay_mask_plot(arr, ring_mask_M_MD, k4, "ring mask")

# ╔═╡ d3f08c32-1245-4021-a55a-02c5b4e39936
begin
	single_ring_mask_M_MD = Bool.(ring_mask_M_MD[:, :, 3])
	s_bkg_M_MD = mean(single_arr[single_ring_mask_M_MD])
end

# ╔═╡ 3a44a556-9c2a-4ffe-9e63-de5aef808b4e
begin
	alg_M_MD = Integrated(arr[mask_M_MD_3D])
	mass_m_md = score(s_bkg_M_MD, S_Obj_MD, pixel_size, ρ_md, alg_M_MD)
end

# ╔═╡ 7442a4c4-611a-42bb-976e-2fd081fe9f5c
md"""
## Low Density
"""

# ╔═╡ f77b6b75-8cf3-4d25-84c8-87b836127d54
begin
	mask_M_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_LD_3D[:, :, z] = mask_M_LD
	end
end;

# ╔═╡ a7ad848f-7650-41a3-b59a-a03cb9e5bb72
md"""
#### Dilated mask
"""

# ╔═╡ 161915c6-755c-43ab-b420-209fd1ecc3ff
dilated_mask_M_LD = dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))));

# ╔═╡ f3e673c1-9894-4977-9bdc-27f171ad996a
@bind l2 overlay_mask_bind(dilated_mask_M_LD)

# ╔═╡ af4040fa-bbb9-49bb-815f-696bd591c172
overlay_mask_plot(arr, dilated_mask_M_LD, l2, "dilated mask")

# ╔═╡ a34cf899-0ca4-4f7e-9c7a-f5e9ab4a6e44
md"""
#### Ring (background) mask
"""

# ╔═╡ a1ac80fe-ff14-49b4-9955-89f779bb2f51
ring_mask_M_LD = dilate(dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))))) - dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))));

# ╔═╡ ab7c72ad-8b8a-414c-a535-2e580056fc81
@bind l4 overlay_mask_bind(ring_mask_M_LD)

# ╔═╡ b1f93bc7-378a-4fdf-8bd5-30e5bcf6a233
overlay_mask_plot(arr, ring_mask_M_LD, l4, "ring mask")

# ╔═╡ 91c3c1a7-4f5c-40fa-a232-26071d036cea
begin
	single_ring_mask_M_LD = Bool.(ring_mask_M_LD[:, :, 3])
	s_bkg_M_LD = mean(single_arr[single_ring_mask_M_LD])
end

# ╔═╡ 7144c4c0-2a66-4c48-ac7a-60329492202a
begin
	alg_M_LD = Integrated(arr[mask_M_LD_3D])
	mass_m_ld = score(s_bkg_M_LD, S_Obj_LD, pixel_size, ρ_LD, alg_M_LD)
end

# ╔═╡ 3f3bd27b-956c-4fe7-8dd2-a1990ab9112b
md"""
# Score Small Inserts
"""

# ╔═╡ 60612e8a-9be7-49f2-bb8a-8780bfd3cb8e
md"""
## High Density
"""

# ╔═╡ 729ddba3-2000-4d83-a738-4185d475cd0a
begin
	mask_S_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_HD_3D[:, :, z] = mask_S_HD
	end
end;

# ╔═╡ 9bdddaf9-5679-4c9b-849c-44b1aef37aa5
md"""
#### Dilated mask
"""

# ╔═╡ b976bd96-eed4-463f-9387-27d2a88b1320
dilated_mask_S_HD = dilate(((dilate(dilate((mask_S_HD_3D))))));

# ╔═╡ 38047c33-debf-41f9-b19e-d76fc9c55cdc
@bind m2 overlay_mask_bind(dilated_mask_S_HD)

# ╔═╡ 4ee3ea5b-b2e1-47be-a456-96a6c8d47914
overlay_mask_plot(arr, dilated_mask_S_HD, m2, "dilated mask")

# ╔═╡ aab1058e-5a3e-4f42-a2cf-804c16da25f8
md"""
#### Ring (background) mask
"""

# ╔═╡ 880204fc-8ace-476d-b839-e45b3032fb85
ring_mask_S_HD = dilate((dilate((dilate(mask_S_HD_3D))))) - dilate(dilate(((mask_S_HD_3D))));

# ╔═╡ d7541fbf-4c61-43ba-8822-8abf2f818538
@bind m4 overlay_mask_bind(ring_mask_S_HD)

# ╔═╡ c1d86ff7-4eb2-41cf-8e1d-25acb9d544df
overlay_mask_plot(arr, ring_mask_S_HD, m4, "ring mask")

# ╔═╡ 983b8e00-4ec8-490b-8549-15ef8940a8b6
begin
	single_ring_mask_S_HD = Bool.(ring_mask_S_HD[:, :, 3])
	s_bkg_S_HD = mean(single_arr[single_ring_mask_S_HD])
end

# ╔═╡ d27e9e59-7743-4df4-b093-22d474c50390
begin
	alg_S_HD = Integrated(arr[mask_S_HD_3D])
	mass_s_hd = score(s_bkg_S_HD, S_Obj_HD, pixel_size, ρ_hd, alg_S_HD)
	if mass_s_hd < 0
		mass_s_hd = 0
	end
	mass_s_hd
end

# ╔═╡ a3ef8573-0a5a-4a67-a7f0-51dfae9258a4
md"""
## Medium Density
"""

# ╔═╡ 52c69f0b-925f-474d-bd58-1fe9ad8295df
begin
	mask_S_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_MD_3D[:, :, z] = mask_S_MD
	end
end;

# ╔═╡ e4c07fd9-9f9d-433e-8f20-864fbd5bbc50
md"""
#### Dilated mask
"""

# ╔═╡ dbc3b24e-75ce-493e-9f51-9cf53730c8b8
dilated_mask_S_MD = dilate(((dilate(dilate(mask_S_MD_3D)))));

# ╔═╡ 14cefd3b-bbc0-466c-af9c-280ae47edb0a
@bind n2 overlay_mask_bind(dilated_mask_S_MD)

# ╔═╡ 789727b8-c575-480e-bacf-5384994ce992
overlay_mask_plot(arr, dilated_mask_S_MD, n2, "dilated mask")

# ╔═╡ f35e9a1c-78bb-40d3-963f-0505f1e9c2d2
md"""
#### Ring (background) mask
"""

# ╔═╡ 1dba9f10-0f08-49bb-8f3c-b3ff444e6b92
ring_mask_S_MD = dilate(((dilate(dilate(mask_S_MD_3D))))) - dilate(dilate(((mask_S_MD_3D))));

# ╔═╡ 989b53f7-e9c5-4786-bf5f-de69d0cf8208
@bind n4 overlay_mask_bind(ring_mask_S_MD)

# ╔═╡ 28299916-5ee5-487f-b912-49077e47a4c4
overlay_mask_plot(arr, ring_mask_S_MD, n4, "ring mask")

# ╔═╡ 2299fda8-be4f-456f-ad3b-d7a335902063
begin
	single_ring_mask_S_MD = Bool.(ring_mask_S_MD[:, :, 3])
	s_bkg_S_MD = mean(single_arr[single_ring_mask_S_MD])
end

# ╔═╡ 5c9918ea-9b4c-416b-ba54-5e347a5abc5d
begin
	alg_S_MD = Integrated(arr[mask_S_MD_3D])
	mass_s_md = score(s_bkg_S_MD, S_Obj_MD, pixel_size, ρ_md, alg_S_MD)
end

# ╔═╡ 520f4533-e478-452e-accb-b6e879c8973a
md"""
## Low Density
"""

# ╔═╡ 7862cb1b-7e47-42a3-85be-7bd00dbd4b2e
begin
	mask_S_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_LD_3D[:, :, z] = mask_S_LD
	end
end;

# ╔═╡ 2d33dbe5-d13f-4ae7-ae1e-70a281b5dcc7
md"""
#### Dilated mask
"""

# ╔═╡ 7267153d-4d66-40d2-b0ed-f21e086a118a
dilated_mask_S_LD = dilate(((dilate(dilate(mask_S_LD_3D)))));

# ╔═╡ d1bbc7a3-5595-4991-8f69-cb8073612194
@bind o2 overlay_mask_bind(dilated_mask_S_LD)

# ╔═╡ 413a75be-42f9-4c26-91c1-38e68b1b32d3
overlay_mask_plot(arr, dilated_mask_S_LD, o2, "dilated mask")

# ╔═╡ 46eeb5c1-25c6-4b24-b81c-7a94d5221b2a
md"""
#### Ring (background) mask
"""

# ╔═╡ 6ec91801-5e22-4e4f-8e34-fde6fa6956ac
ring_mask_S_LD = (dilate((dilate(dilate(mask_S_LD_3D))))) - dilate(dilate(((mask_S_LD_3D))));

# ╔═╡ 891b8711-5531-4623-9234-57e1e492d28c
@bind o4 overlay_mask_bind(ring_mask_S_LD)

# ╔═╡ 31e85f55-682c-456c-807e-46abeb281a3f
overlay_mask_plot(arr, ring_mask_S_LD, o4, "ring mask")

# ╔═╡ f38c7244-3ba7-46fb-b918-6e35f227a284
begin
	single_ring_mask_S_LD = Bool.(ring_mask_S_LD[:, :, 3])
	s_bkg_S_LD = mean(single_arr[single_ring_mask_S_LD])
end

# ╔═╡ b9263d5a-dcc6-4326-83d4-c12d259deaf5
begin
	alg_S_LD = Integrated(arr[mask_S_LD_3D])
	mass_s_ld = score(s_bkg_S_LD, S_Obj_LD, pixel_size, ρ_LD, alg_S_LD)
end

# ╔═╡ b7e97710-502e-4896-8e33-b1dccb7d091e
md"""
# Results
"""

# ╔═╡ f42fed63-fa26-4f1c-9598-e5dd3ec16f9c
PhantomSegmentation.get_pixel_size(header)

# ╔═╡ c58e23d2-10bf-4b25-b9e1-0c76389c82e5
inserts = [
	"Low Density",
	"Medium Density",
	"High Density"
]

# ╔═╡ 227f9f48-a8fb-46cf-9324-4b9daefcfbce
volume_gt = [
	7.065,
	63.585,
	176.625
]

# ╔═╡ bc923671-e629-4d73-b386-1b296455af25
ground_truth_mass_large = [
	volume_gt[3] * density_array[2] * 1e-3,
	volume_gt[3] * density_array[3] * 1e-3,
	volume_gt[3] * density_array[4] * 1e-3
] # mg

# ╔═╡ 70613704-fd90-4c81-a00e-958520c3b00f
calculated_mass_large = [
	mass_l_ld,
	mass_l_md,
	mass_l_hd
]

# ╔═╡ 7557b68a-baa3-4404-9af5-40c4512eea65
ground_truth_mass_medium = [
	volume_gt[2] * density_array[2] * 1e-3,
	volume_gt[2] * density_array[3] * 1e-3,
	volume_gt[2] * density_array[4] * 1e-3
]

# ╔═╡ 46fa5daa-b558-48c6-aa78-8b6f8be1382f
calculated_mass_medium = [
	mass_m_ld,
	mass_m_md,
	mass_m_hd
]

# ╔═╡ 0ed73fdc-f00d-45c3-a99c-c5ca72b16525
ground_truth_mass_small = [
	volume_gt[1] * density_array[2] * 1e-3,
	volume_gt[1] * density_array[3] * 1e-3,
	volume_gt[1] * density_array[4] * 1e-3
]

# ╔═╡ 95e78cac-2114-4e45-841f-93bfc87a4786
calculated_mass_small = [
	mass_s_ld,
	mass_s_md,
	mass_s_hd
]

# ╔═╡ e625e6ac-a430-4c20-b91b-d979d4d7d1c8
df = DataFrame(
	kern = kern,
	scan = scan,
	inserts = inserts,
	ground_truth_mass_large = ground_truth_mass_large,
	calculated_mass_large = calculated_mass_large,
	ground_truth_mass_medium = ground_truth_mass_medium,
	calculated_mass_medium = calculated_mass_medium,
	ground_truth_mass_small = ground_truth_mass_small,
	calculated_mass_small = calculated_mass_small
)

# ╔═╡ 3fd1f454-6aeb-464f-9a29-8c4c498a4762
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

# ╔═╡ 9ed2bd38-f96f-4163-a759-518be99604ca
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

# ╔═╡ 4ca0c660-4e53-4ef8-9bab-183bb2d88001
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

# ╔═╡ 59b6082c-24c2-4aac-b1bd-45f799dfcb5d
percent_error_large = (abs.(ground_truth_mass_large - calculated_mass_large) ./ ground_truth_mass_large) .* 100

# ╔═╡ 8e21b221-200a-4ca0-ad7e-6fdadc94a969
percent_error_medium = (abs.(ground_truth_mass_medium - calculated_mass_medium) ./ ground_truth_mass_medium) .* 100

# ╔═╡ 66e71346-a642-41af-bde2-d42fceb5749a
percent_error_small= (abs.(ground_truth_mass_small - calculated_mass_small) ./ ground_truth_mass_small) .* 100

# ╔═╡ 7aa406f4-31a5-42b9-8198-074ca3a2363e
md"""
### Save Results
"""

# ╔═╡ 30a3ee3b-8e82-4282-9c66-f959b4e3df14
# if ~isdir(string(cd(pwd, "..") , "/output/", TYPE))
# 	mkdir(string(cd(pwd, "..") , "/output/", TYPE))
# end

# ╔═╡ 3447c350-d4da-4e32-80c4-b2cd36e9c730
# output_path = string(cd(pwd, "..") , "/output/", TYPE, "/", scan, ".csv")

# ╔═╡ 6646b4f9-9cba-4ec4-b526-0f847c543473
# CSV.write(output_path, df)

# ╔═╡ ddb689fa-a417-4719-b1fa-d0acb0c8d51a
md"""
### Save full df
"""

# ╔═╡ b740e1ce-98bb-4a54-906e-8fabb10f6ef8
dfs = []

# ╔═╡ 33153ff9-2010-4cae-a355-9448ea6fbb95
push!(dfs, df)

# ╔═╡ 8d12c98a-3ad4-411d-81e9-8a63fdb18081
if length(dfs) == 12
	global new_df = vcat(dfs[1:12]...)
	output_path_new = string(cd(pwd, "..") , "/output/", TYPE, "/", "full.csv")
	CSV.write(output_path_new, new_df)
end

# Cell order:
# ╠═45f704d4-66e5-49db-aef5-05132f3853ee
# ╠═e7c4aaab-83ef-4256-b392-f8ef7a899a05
# ╟─48e3097f-0767-4dad-9d7c-94e0899f790b
# ╠═bc9383c0-e477-4e1a-a2fa-7f5c1d29f103
# ╟─6aa51429-981a-4dea-a0f6-2935867d5b2a
# ╠═e849cf69-65c7-4e5e-9688-fc249d471f2c
# ╠═3e649fee-a2e9-4f0f-b705-3966f69d97ea
# ╠═a327ef18-2941-4783-9266-7332826eaf58
# ╠═b44d18f8-1b86-4235-bb5b-7a77a1af55e0
# ╠═b7e0b678-0627-44fe-b0cb-3ef2bccae6a7
# ╠═3813d24f-eef7-4d96-8adc-69a9e22fe9a1
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
# ╠═dd399c0e-f8e4-464c-8964-bdc0dd657202
# ╟─417d150e-0df9-4963-b59e-ebd9acf6d4a0
# ╠═5f10075b-4d95-4d53-a22f-016749fb7583
# ╟─b921dcf8-54ea-420b-9285-23e38d5ce433
# ╠═ece7d76b-93c8-430c-8679-07cf92585949
# ╠═fa58a5ba-7aaf-4be6-8796-ecef516d8d53
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
# ╠═6c3c7d22-78e6-428e-ac62-90a473b06e2d
# ╠═f42a62d8-11e9-48b8-b6bc-ae25b6797707
# ╠═b407d628-0db0-4d13-9a6c-8e4ff5a84cbd
# ╠═dbec08a2-96c5-488b-89b9-5df6fe1f5a12
# ╟─e50baf62-0a02-47ee-8440-551f68baee0d
# ╠═3fa80bee-20f6-4e9c-a14d-75a77482a2c3
# ╠═04c77c9d-c498-4e09-abff-9c361e887c56
# ╠═b03c711f-9e0e-4089-add1-5abbde81cce1
# ╟─3b30e707-02b0-49c2-a0b3-25b7fcbf79ee
# ╠═d27df860-7ca8-4e93-90e4-a58605aa1aaa
# ╠═5dfb50b9-7545-4948-8026-2074ce5d86cf
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

# ╔═╡ Cell order:
# ╠═2c064fd4-3719-4503-baa3-61a453e3642d
# ╠═355b6a05-065e-4b2c-b65f-af9c324756d0
# ╠═739ccb4d-1378-4241-9e15-fffaa4522080
# ╠═29290cd1-945c-4ec2-9cc1-fe65c7d21dae
# ╠═a5b8ce68-4ad4-499c-8caf-c0c477e3d53d
# ╠═6aaf8df0-a546-432c-ad17-d8d735fa0919
# ╠═dd590806-ee35-47a8-898c-d158d6133eee
# ╠═2ab91f08-fac4-48d6-9776-ac47002159ab
# ╠═ce9c455d-ad21-4c27-8468-84be6f14f5ed
# ╠═56c069d9-d1d6-4f6f-b2d1-cef68a970d9b
# ╠═621eb1c3-a294-4620-b79c-1c98e501a21f
# ╠═d4cea938-0455-4f3d-b9f9-9a5d6c9f6967
# ╠═80fff6f2-c369-4efc-84a2-da15e7b718f7
# ╠═4f9e6c87-1d73-441d-a0e6-b0c2934dd035
# ╠═2cdc2b8c-7b26-421e-a008-62cc56c026fd
# ╠═4256495e-319a-4587-89a5-0185a2cb22e4
# ╠═6f5312d8-0148-4b5b-ab62-8a02e1e66b5a
# ╠═14c64d65-06d9-4919-b8fd-fc5d7d5eac63
# ╠═91eef4a3-80df-4764-97d0-1d23cb34de7b
# ╠═9df0b48d-e053-4d40-bad8-c12c0813beb8
# ╠═080e44cd-7918-4a79-bfa5-720830738388
# ╠═312dbc4e-2f66-48da-8801-92d6f9e9301f
# ╠═b6ae6e8f-b1b4-492c-977f-aa0e72d76119
# ╠═773573d7-97c6-4b93-9e3d-ce544918bb4f
# ╠═e7c9f252-ca91-4a7a-84c3-7bcd647af9b2
# ╠═adfe2070-3ae3-412f-b7f3-99622d6003c5
# ╠═e14c4255-b10d-4502-8b7c-94a5daed1b09
# ╠═58ff4866-983d-4481-9442-d9085e97f855
# ╠═47982924-7917-41f0-989b-3233c5e5f402
# ╠═1d800cab-b66a-4497-9e7c-1367bfdd3a87
# ╠═c4c6c326-9b15-478b-b462-f4f20b4601fb
# ╠═2630ef4c-1870-4d83-b191-027a1db1bd47
# ╠═b7f9e7cd-c626-4df5-833e-e4491baaaf87
# ╠═742b3c55-5131-4554-8980-1284fa6d49a2
# ╠═fc069260-4b7f-4ee5-8bc3-29d481e0d6cc
# ╠═5ad725f9-2997-4af7-b324-1224fd0e3017
# ╠═bed91830-f655-4b9f-845e-0353aecf7734
# ╠═8ab94c76-47f2-48d0-888c-f4ea960cf276
# ╠═d721be4b-5599-49bc-abd9-160175455e22
# ╠═0f40aa11-8ca6-4574-905b-eb6ee9772b90
# ╠═d815512f-3ef6-4481-8b12-33b1288edc87
# ╠═65f63514-38ac-4662-b468-4883604276ee
# ╠═15fe3d10-a572-405d-b1f6-09f3637d58f8
# ╠═8fe16dbf-14aa-4f4a-8f22-e5bffac397d5
# ╠═76d9bbbe-db91-4bb6-b271-1aea6e37e9c6
# ╠═0cd096ca-2b26-4dda-9cb6-33ad2b296db0
# ╠═d070d007-72a0-460f-b1d3-c1798e25a7af
# ╠═dca26385-2c02-4085-9cdf-2e401de044d9
# ╠═bdbe21c9-87c1-45ea-bf8a-4cc7bd65491d
# ╠═afc1d81f-1e02-47da-8773-bedcbf309166
# ╠═8aeca764-a3f1-482a-b9b0-975d5d08fb6a
# ╠═92e85638-3f71-41db-8d1e-6fcf46ceaf69
# ╠═5e98dad3-5bf1-497f-93ea-69d3acdf3411
# ╠═8494fff0-91a3-413f-8a0b-cb9c79223fd0
# ╠═53bf9d17-537b-4d70-a67c-5c88ceef08a9
# ╠═a02e06fe-ff2e-4a46-93f5-e661fae20e4f
# ╠═9aee1a0a-2204-45bf-a6e1-b42be746d351
# ╠═8ea82534-97c5-4494-ad06-a7979f4ff8b3
# ╠═4a2a49a9-c435-4caa-92f2-74b5f92a782b
# ╠═118efcb2-f712-49f2-9589-7453e4e0e342
# ╠═0b505164-0978-4c40-b044-b0211b2bab7c
# ╠═f4585dda-4493-4caf-93f0-36a924b50bea
# ╠═96323aa7-b24f-4942-9197-08cb49e96427
# ╠═5270b57b-5cb0-48b6-bfd8-ffac3463d9e6
# ╠═4a045f83-bc81-42c9-a048-7c20855fbaec
# ╠═b7e539c3-56a6-47e0-ae4a-8a15679331c2
# ╠═e539b82b-945e-463e-a2c9-8240dbb7dd12
# ╠═e976610d-5c02-414e-b8c6-db989b5b098a
# ╠═b0fb5549-e0ec-47a4-bb56-122a1f76c6a8
# ╠═029bd311-f8f8-4a06-bde3-24be773b6970
# ╠═ac775f8c-77e7-4a60-a844-f57b9be48c93
# ╠═4bfdd6ee-2bfe-4e6f-8fd4-acdb90775cbb
# ╠═7894faf9-83c4-431f-a6b4-9baa1f97ffcb
# ╠═2f3a1e98-1c33-47cc-b90a-b80fcf552356
# ╠═7894ca8c-e6d8-402e-8f0d-99e936de513e
# ╠═bc7f7c97-4fbd-44a5-9ee0-f6390b2e2674
# ╠═0d2c6476-cf89-445f-ac9f-db54972aec02
# ╠═3dcb027b-87a7-411d-a51f-23b7606cd672
# ╠═bcaafe33-b0dd-4ec6-96cf-ff8fe04d9ad5
# ╠═05c4ff2b-7d4f-45a4-9a8d-ce1afd429f4b
# ╠═7332bcca-cd2a-4809-84e1-19c0829b0055
# ╠═5a754fcc-3dd8-4ecb-b762-ae913a2c3ad5
# ╠═b6f7ed3a-20f1-4136-afc4-0c5ba9632343
# ╠═450f264d-cd87-476c-9036-c64ba9de3c63
# ╠═d4c5d18a-c7a1-48e1-876c-d4c4b94eab8f
# ╠═b06deeb6-7090-4d17-bdba-1d90dad6b202
# ╠═145e5acb-9e35-40c8-b326-fc71c4a07dc5
# ╠═5a3aadde-5120-46a1-b861-072d27ce46b7
# ╠═95cac882-2c04-4728-8c69-46f263dbbd1d
# ╠═c3cc88b1-76b4-4a61-82bf-f01c403124a3
# ╠═c25633c7-12ee-4497-b733-fa6d785bd77c
# ╠═d9b0b828-f910-4f10-ba56-b4e6819a9004
# ╠═2cc31f1c-8a04-4d00-8ec9-c21c6289a587
# ╠═44145515-28f3-4f57-8c29-01681249c9dd
# ╠═5ad20c65-58c1-4fa5-9ce0-15336843895b
# ╠═f496aef2-a7cc-4dc9-bcfc-ab67ef00e8a3
# ╠═14d13b95-252f-4bf1-9f39-94aa0a004959
# ╠═f64803d5-610a-4e01-a38e-48d03eddd1da
# ╠═fb586eda-bfbe-4427-bd85-476fec08ec09
# ╠═b16ad6b0-12d8-4a28-8278-f286f7f0f11f
# ╠═14fa6374-1d4c-4194-9b7f-47a35232a328
# ╠═060a220d-b6be-4cf2-a8b5-ed3f99ce3ca8
# ╠═9ee074e8-5627-4ece-9150-057e4b3cf8a4
# ╠═ef0fa44f-75fe-45bd-b61d-34c34dbee506
# ╠═13509cd8-557a-4883-8152-2ab1b29830cf
# ╠═1f8b288a-25bc-41d3-9d24-7b325543dca2
# ╠═371e86b7-5a7e-44fd-b8e6-d2cf85cb4783
# ╠═e5bed8d6-e403-4ed3-b21e-2eb513a2c112
# ╠═52499aa5-93f2-4ecf-ad28-6454196de7b5
# ╠═a3f66bd2-13d9-42db-8660-4a135e608d03
# ╠═58645445-5011-4139-97c0-1891ae651a5c
# ╠═c635a53e-ff58-4270-a4a2-76c8358567d4
# ╠═0f5e6427-2480-4958-bbad-c54499288f61
# ╠═839873c2-2927-44c7-86f6-cf729a42d6ea
# ╠═63ba13f6-c17e-47ad-90e2-9c20a3c606b9
# ╠═8b868d1c-06b0-4169-bd52-f66c6409b360
# ╠═ccae0b6e-bd74-43e0-ae77-bf26f4a12a1e
# ╠═8c34dc16-ce3e-4ddd-bab8-fdf8fb633596
# ╠═30c313d4-5851-4ed4-9d0e-d64b36dcbacf
# ╠═1511cbf4-cdf4-4bb0-a678-fd1c327b23d7
# ╠═ff09ac24-388b-4f94-82ce-1c8f9ce22684
# ╠═35d35d55-397d-45eb-93c2-b846acef8553
# ╠═f459d9c5-705a-4e54-a97e-a1f5710f4757
# ╠═3e4ab08d-77c7-4304-acec-13857ac5c4bb
# ╠═d3f08c32-1245-4021-a55a-02c5b4e39936
# ╠═3a44a556-9c2a-4ffe-9e63-de5aef808b4e
# ╠═7442a4c4-611a-42bb-976e-2fd081fe9f5c
# ╠═f77b6b75-8cf3-4d25-84c8-87b836127d54
# ╠═a7ad848f-7650-41a3-b59a-a03cb9e5bb72
# ╠═161915c6-755c-43ab-b420-209fd1ecc3ff
# ╠═f3e673c1-9894-4977-9bdc-27f171ad996a
# ╠═af4040fa-bbb9-49bb-815f-696bd591c172
# ╠═a34cf899-0ca4-4f7e-9c7a-f5e9ab4a6e44
# ╠═a1ac80fe-ff14-49b4-9955-89f779bb2f51
# ╠═ab7c72ad-8b8a-414c-a535-2e580056fc81
# ╠═b1f93bc7-378a-4fdf-8bd5-30e5bcf6a233
# ╠═91c3c1a7-4f5c-40fa-a232-26071d036cea
# ╠═7144c4c0-2a66-4c48-ac7a-60329492202a
# ╠═3f3bd27b-956c-4fe7-8dd2-a1990ab9112b
# ╠═60612e8a-9be7-49f2-bb8a-8780bfd3cb8e
# ╠═729ddba3-2000-4d83-a738-4185d475cd0a
# ╠═9bdddaf9-5679-4c9b-849c-44b1aef37aa5
# ╠═b976bd96-eed4-463f-9387-27d2a88b1320
# ╠═38047c33-debf-41f9-b19e-d76fc9c55cdc
# ╠═4ee3ea5b-b2e1-47be-a456-96a6c8d47914
# ╠═aab1058e-5a3e-4f42-a2cf-804c16da25f8
# ╠═880204fc-8ace-476d-b839-e45b3032fb85
# ╠═d7541fbf-4c61-43ba-8822-8abf2f818538
# ╠═c1d86ff7-4eb2-41cf-8e1d-25acb9d544df
# ╠═983b8e00-4ec8-490b-8549-15ef8940a8b6
# ╠═d27e9e59-7743-4df4-b093-22d474c50390
# ╠═a3ef8573-0a5a-4a67-a7f0-51dfae9258a4
# ╠═52c69f0b-925f-474d-bd58-1fe9ad8295df
# ╠═e4c07fd9-9f9d-433e-8f20-864fbd5bbc50
# ╠═dbc3b24e-75ce-493e-9f51-9cf53730c8b8
# ╠═14cefd3b-bbc0-466c-af9c-280ae47edb0a
# ╠═789727b8-c575-480e-bacf-5384994ce992
# ╠═f35e9a1c-78bb-40d3-963f-0505f1e9c2d2
# ╠═1dba9f10-0f08-49bb-8f3c-b3ff444e6b92
# ╠═989b53f7-e9c5-4786-bf5f-de69d0cf8208
# ╠═28299916-5ee5-487f-b912-49077e47a4c4
# ╠═2299fda8-be4f-456f-ad3b-d7a335902063
# ╠═5c9918ea-9b4c-416b-ba54-5e347a5abc5d
# ╠═520f4533-e478-452e-accb-b6e879c8973a
# ╠═7862cb1b-7e47-42a3-85be-7bd00dbd4b2e
# ╠═2d33dbe5-d13f-4ae7-ae1e-70a281b5dcc7
# ╠═7267153d-4d66-40d2-b0ed-f21e086a118a
# ╠═d1bbc7a3-5595-4991-8f69-cb8073612194
# ╠═413a75be-42f9-4c26-91c1-38e68b1b32d3
# ╠═46eeb5c1-25c6-4b24-b81c-7a94d5221b2a
# ╠═6ec91801-5e22-4e4f-8e34-fde6fa6956ac
# ╠═891b8711-5531-4623-9234-57e1e492d28c
# ╠═31e85f55-682c-456c-807e-46abeb281a3f
# ╠═f38c7244-3ba7-46fb-b918-6e35f227a284
# ╠═b9263d5a-dcc6-4326-83d4-c12d259deaf5
# ╠═b7e97710-502e-4896-8e33-b1dccb7d091e
# ╠═f42fed63-fa26-4f1c-9598-e5dd3ec16f9c
# ╠═c58e23d2-10bf-4b25-b9e1-0c76389c82e5
# ╠═227f9f48-a8fb-46cf-9324-4b9daefcfbce
# ╠═bc923671-e629-4d73-b386-1b296455af25
# ╠═70613704-fd90-4c81-a00e-958520c3b00f
# ╠═7557b68a-baa3-4404-9af5-40c4512eea65
# ╠═46fa5daa-b558-48c6-aa78-8b6f8be1382f
# ╠═0ed73fdc-f00d-45c3-a99c-c5ca72b16525
# ╠═95e78cac-2114-4e45-841f-93bfc87a4786
# ╠═e625e6ac-a430-4c20-b91b-d979d4d7d1c8
# ╠═3fd1f454-6aeb-464f-9a29-8c4c498a4762
# ╠═9ed2bd38-f96f-4163-a759-518be99604ca
# ╠═4ca0c660-4e53-4ef8-9bab-183bb2d88001
# ╠═59b6082c-24c2-4aac-b1bd-45f799dfcb5d
# ╠═8e21b221-200a-4ca0-ad7e-6fdadc94a969
# ╠═66e71346-a642-41af-bde2-d42fceb5749a
# ╠═7aa406f4-31a5-42b9-8198-074ca3a2363e
# ╠═30a3ee3b-8e82-4282-9c66-f959b4e3df14
# ╠═3447c350-d4da-4e32-80c4-b2cd36e9c730
# ╠═6646b4f9-9cba-4ec4-b526-0f847c543473
# ╠═ddb689fa-a417-4719-b1fa-d0acb0c8d51a
# ╠═b740e1ce-98bb-4a54-906e-8fabb10f6ef8
# ╠═33153ff9-2010-4cae-a355-9448ea6fbb95
# ╠═8d12c98a-3ad4-411d-81e9-8a63fdb18081

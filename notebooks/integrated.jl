### A Pluto.jl notebook ###
# v0.19.18

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
	using Pkg
	Pkg.activate(".")

    using PlutoUI, Statistics, CSV, DataFrames, GLM, CairoMakie, HypothesisTests, Colors, MLJBase, DICOM, DICOMUtils, PhantomSegmentation, CalciumScoring, ImageMorphology, ImageFiltering, Noise
    using StatsBase: quantile!, rmsd
end

# ╔═╡ e7c4aaab-83ef-4256-b392-f8ef7a899a05
TableOfContents()

# ╔═╡ 48e3097f-0767-4dad-9d7c-94e0899f790b
md"""
## Load DICOMS

All you need to do is set `base_path` once and leave it. After that, the only thing that should change is the `VENDOR`, once for every set, and the `SCAN_NUMBER`, once for each scan.
"""

# ╔═╡ bc9383c0-e477-4e1a-a2fa-7f5c1d29f103
begin
	VENDORS = ["80", "100", "120", "135"]
	SIZES = ["small", "medium", "large"]
	DENSITIES = ["low", "normal"]

	IMAGES = "images_new"

	VENDOR = VENDORS[1]
	SIZE = SIZES[1]
	DENSITY = DENSITIES[2]
	
    BASE_PATH = joinpath("/Users/daleblack/Google Drive/dev/MolloiLab/cac-simulation", IMAGES, SIZE, DENSITY)
	root_path = joinpath(BASE_PATH, VENDOR)
	dcm_path_list = dcm_list_builder(root_path)
	pth = dcm_path_list[1]
	scan = basename(pth)
	header, dcm_array, slice_thick_ori1 = dcm_reader(pth)
	kV = parse(Int64, VENDOR)
end

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
        container[i, :] = collect(tuple_array[i])
    end
    return container
end

# ╔═╡ 483a14dd-e798-41ed-9144-13678f8b8461
function overlay_mask_bind(mask)
    indices = findall(x -> x == 1, mask)
    indices = Tuple.(indices)
    label_array = collect_tuple(indices)
    zs = unique(label_array[:, 3])
    return PlutoUI.Slider(1:length(zs); default=3, show_value=true)
end

# ╔═╡ 2cfdeb06-71a2-411b-accb-2b4a3f56b477
function overlay_mask_plot(array, mask, var, title::AbstractString)
    indices = findall(x -> x == 1, mask)
    indices = Tuple.(indices)
    label_array = collect_tuple(indices)
    zs = unique(label_array[:, 3])
    indices_lbl = findall(x -> x == zs[var], label_array[:, 3])

    fig = Figure()
    ax = Makie.Axis(fig[1, 1])
    ax.title = title
    heatmap!(array[:, :, zs[var]]; colormap=:grays)
    scatter!(
        label_array[:, 1][indices_lbl],
        label_array[:, 2][indices_lbl];
        markersize=1,
        color=:red,
    )
    return fig
end

# ╔═╡ ab294a53-0ac8-4a40-8a98-c5271c66d958
function create_mask(array, mask)
    @assert size(array) == size(mask)
    idxs = findall(x -> x == true, mask)
    overlayed_mask = zeros(size(array))
    for idx in idxs
        overlayed_mask[idx] = array[idx]
    end
    return overlayed_mask
end

# ╔═╡ 3f12b12c-90a9-48cc-99db-3545e4d0a617
function dilate_mask_large(mask)
    return dilate(mask)
end

# ╔═╡ 29bb75c2-ff32-4a01-b04c-614647c66c67
function ring_mask_large(dilated_mask)
    return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask)
end

# ╔═╡ eacd75fc-1a1f-4f55-81c6-f897154adc38
function dilate_mask_medium(mask)
    return (mask)
end

# ╔═╡ 6dfd0583-21c2-400a-9ab5-3a09a95e5adc
function ring_mask_medium(dilated_mask)
    return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask)
end

# ╔═╡ 6d306954-6c37-4138-924d-151f0b2dd068
function dilate_mask_small(mask)
    return (mask)
end

# ╔═╡ 8678df23-ef1c-42b0-8102-84ca44f4f7d8
function ring_mask_small(dilated_mask)
    return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask)
end

# ╔═╡ e4d10b8b-93fb-415e-bd35-aca8fee7510e
function dilate_mask_large_bkg(mask)
    return dilate(dilate(mask))
end

# ╔═╡ d6a34d2d-937f-4f59-9e77-c8d484bbe75b
function dilate_mask_medium_bkg(mask)
    return dilate(mask)
end

# ╔═╡ 2985cf19-4a56-478f-9eda-5387b9fbec08
function dilate_mask_small_bkg(mask)
    return (mask)
end

# ╔═╡ b3bae2ad-16ed-4381-b6a4-448cb5f0c6c1
md"""
## Segment Heart
"""

# ╔═╡ 4ff19af0-7f53-4f24-b315-08bd54a488e3
masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3) ÷ 2);

# ╔═╡ b74941c4-12d4-4be1-81f7-ef1f4d288984
@bind a PlutoUI.Slider(1:size(masked_array, 3), default=10, show_value=true)

# ╔═╡ 84e56c91-14a4-45b7-81a1-77e778bae695
heatmap(transpose(masked_array[:, :, a]); colormap=:grays)

# ╔═╡ 4a4cf730-1307-451e-be5c-f9c1183cdd5f
let
    fig = Figure()

    ax = Makie.Axis(fig[1, 1])
    heatmap!(transpose(dcm_array[:, :, 5]); colormap=:grays)

	save("/Users/daleblack/Google Drive/Research/Papers/My Papers/cac-simulation/figures-review/simulated_phantom.png", fig)
    fig
end

# ╔═╡ abfd965c-df00-4028-a0f3-8356a261b527
let
    fig = Figure()

    ax = Makie.Axis(fig[1, 1])
    ax.title = "Raw DICOM Array"
    heatmap!(transpose(dcm_array[:, :, 5]); colormap=:grays)
    scatter!(
        center_insert[2]:center_insert[2],
        center_insert[1]:center_insert[1];
        markersize=10,
        color=:red,
    )
    fig
end

# ╔═╡ 41541a09-609b-4f46-9a25-eb974422efc0
begin
    fig2 = Figure()

    ax2 = Makie.Axis(fig2[1, 1])
    ax2.title = "Mask Array"
    heatmap!(transpose(mask); colormap=:grays)
    scatter!(
        center_insert[2]:center_insert[2],
        center_insert[1]:center_insert[1];
        markersize=10,
        color=:red,
    )
    fig2
end

# ╔═╡ e6753c9a-d172-47fe-b066-78cb43df2c89
begin
    fig3 = Figure()

    ax3 = Makie.Axis(fig3[1, 1])
    ax3.title = "Masked DICOM Array"
    heatmap!(transpose(masked_array[:, :, 1]); colormap=:grays)
    scatter!(
        center_insert[2]:center_insert[2],
        center_insert[1]:center_insert[1];
        markersize=10,
        color=:red,
    )
    fig3
end

# ╔═╡ f8a4eacd-2f69-4755-b017-5a0b0dda1004
md"""
## Segment Calcium Rod
"""

# ╔═╡ ac9d2652-6d69-4cf3-acbf-2e141fb633f0
begin
	if DENSITY == "low"
		thresh = 55
	elseif DENSITY == "normal"
		thresh = 130
	end

	calcium_image, slice_CCI, quality_slice, cal_rod_slice = mask_rod(
		masked_array, header; calcium_threshold=thresh
	)
end;

# ╔═╡ 417d150e-0df9-4963-b59e-ebd9acf6d4a0
@bind c PlutoUI.Slider(1:size(calcium_image, 3), default=5, show_value=true)

# ╔═╡ 5f10075b-4d95-4d53-a22f-016749fb7583
heatmap(transpose(calcium_image[:, :, c]); colormap=:grays)

# ╔═╡ b921dcf8-54ea-420b-9285-23e38d5ce433
md"""
## Load Masks
"""

# ╔═╡ 05ed60db-b4d8-4cdc-9a54-b108ade22557
begin
    root_new = joinpath("/Users/daleblack/Google Drive/dev/MolloiLab/cac-simulation/julia_arrays", SIZE)
	mask_L_HD = Array(CSV.read(joinpath(root_new, "mask_L_HD.csv"), DataFrame; header=false))
	mask_M_HD = Array(CSV.read(joinpath(root_new, "mask_M_HD.csv"), DataFrame; header=false))
	mask_S_HD = Array(CSV.read(joinpath(root_new, "mask_S_HD.csv"), DataFrame; header=false))
	mask_L_MD = Array(CSV.read(joinpath(root_new, "mask_L_MD.csv"), DataFrame; header=false))
	mask_M_MD = Array(CSV.read(joinpath(root_new, "mask_M_MD.csv"), DataFrame; header=false))
	mask_S_MD = Array(CSV.read(joinpath(root_new, "mask_S_MD.csv"), DataFrame; header=false))
	mask_L_LD = Array(CSV.read(joinpath(root_new, "mask_L_LD.csv"), DataFrame; header=false))
	mask_M_LD = Array(CSV.read(joinpath(root_new, "mask_M_LD.csv"), DataFrame; header=false))
	mask_S_LD = Array(CSV.read(joinpath(root_new, "mask_S_LD.csv"), DataFrame; header=false))
end;

# ╔═╡ 8bb2d451-b575-40bf-b31e-51e80392ef51
masks =
    mask_L_HD +
    mask_M_HD +
    mask_S_HD +
    mask_L_MD +
    mask_M_MD +
    mask_S_MD +
    mask_L_LD +
    mask_M_LD +
    mask_S_LD;

# ╔═╡ 4862496e-464d-40c3-baef-e5c507028d66
heatmap(transpose(masks); colormap=:grays)

# ╔═╡ 1545510f-1753-4c60-88d0-75606e4903f1
md"""
## Prepare Array Masks
"""

# ╔═╡ 16fb4a53-116a-404c-95ae-bedcadbc8f9f
begin
	arr = masked_array[:, :, 4:6]
	pixel_size = DICOMUtils.get_pixel_size(header)
	voxel_size = pixel_size[1] * pixel_size[2] * pixel_size[3]
end;

# ╔═╡ ccec799c-e328-4c41-8279-c71f27bd92ec
begin
	## Large
	### High Density
	mask_L_HD_3D = Array{Bool}(undef, size(arr))
	for z in axes(arr, 3)
		mask_L_HD_3D[:, :, z] = mask_L_HD
	end
	dilated_mask_L_HD = dilate_mask_large(mask_L_HD_3D)
	ring_mask_L_HD = ring_mask_large(dilated_mask_L_HD)

	### Medium Density
	mask_L_MD_3D = Array{Bool}(undef, size(arr))
	for z in axes(arr, 3)
		mask_L_MD_3D[:, :, z] = mask_L_MD
	end
	dilated_mask_L_MD = dilate_mask_large(mask_L_MD_3D)
	ring_mask_L_MD = ring_mask_large(dilated_mask_L_MD)

	### Low Density
	mask_L_LD_3D = Array{Bool}(undef, size(arr))
	for z in axes(arr, 3)
		mask_L_LD_3D[:, :, z] = mask_L_LD
	end
	dilated_mask_L_LD = dilate_mask_large(mask_L_LD_3D)
	ring_mask_L_LD = ring_mask_large(dilated_mask_L_LD)


	## Medium 
	### High Density
	mask_M_HD_3D = Array{Bool}(undef, size(arr))
	for z in axes(arr, 3)
		mask_M_HD_3D[:, :, z] = mask_M_HD
	end
	dilated_mask_M_HD = dilate_mask_medium(mask_M_HD_3D)
	ring_mask_M_HD = ring_mask_medium(dilated_mask_M_HD)

	### Medium Density
	mask_M_MD_3D = Array{Bool}(undef, size(arr))
	for z in axes(arr, 3)
		mask_M_MD_3D[:, :, z] = mask_M_MD
	end
	dilated_mask_M_MD = dilate_mask_medium(mask_M_MD_3D)
	ring_mask_M_MD = ring_mask_medium(dilated_mask_M_MD)

	### Low Density
	mask_M_LD_3D = Array{Bool}(undef, size(arr))
	for z in axes(arr, 3)
		mask_M_LD_3D[:, :, z] = mask_M_LD
	end
	dilated_mask_M_LD = dilate_mask_medium(mask_M_LD_3D)
	ring_mask_M_LD = ring_mask_medium(dilated_mask_M_LD)

	## Small
	### High Density
	mask_S_HD_3D = Array{Bool}(undef, size(arr))
	for z in axes(arr, 3)
		mask_S_HD_3D[:, :, z] = mask_S_HD
	end
	dilated_mask_S_HD = dilate_mask_small(mask_S_HD_3D)
	ring_mask_S_HD = ring_mask_small(dilated_mask_S_HD)

	### Medium Density
	mask_S_MD_3D = Array{Bool}(undef, size(arr))
	for z in axes(arr, 3)
		mask_S_MD_3D[:, :, z] = mask_S_MD
	end
	dilated_mask_S_MD = dilate_mask_small(mask_S_MD_3D)
	ring_mask_S_MD = ring_mask_small(dilated_mask_S_MD)

	### Low Density
	mask_S_LD_3D = Array{Bool}(undef, size(arr))
	for z in axes(arr, 3)
		mask_S_LD_3D[:, :, z] = mask_S_LD
	end
	dilated_mask_S_LD = dilate_mask_small(mask_S_LD_3D)
	ring_mask_S_LD = ring_mask_small(dilated_mask_S_LD)
end;

# ╔═╡ 997952cd-dc01-4528-a30a-fdc11a46dad8
md"""
## Calibration Prep
"""

# ╔═╡ f3f73d41-0d98-461b-8673-7207857253ad
begin
	if DENSITY == "low"
		density_array = [0.025, 0.050, 0.100]
	elseif DENSITY == "normal"
		density_array = [0.200, 0.400, 0.800]
	end
	array_filtered = abs.(mapwindow(median, calcium_image[:, :, 2], (3, 3)))
	bool_arr = array_filtered .> 0
	bool_arr_erode = (((erode(erode(bool_arr)))))
	c_img = calcium_image[:, :, 1:3]
	mask_cal_3D = Array{Bool}(undef, size(c_img))
	for z in 1:size(c_img, 3)
		mask_cal_3D[:, :, z] = bool_arr_erode
	end

	hu_calcium = mean(c_img[mask_cal_3D])
	ρ_calcium = 0.2
end;

# ╔═╡ da3c96b9-efdb-40f4-a959-747c097e16b2
heatmap(bool_arr; colormap=:grays)

# ╔═╡ 5fec186d-82bb-4b37-b0c6-c3781f2a166d
heatmap(bool_arr_erode; colormap=:grays)

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

# ╔═╡ b03c711f-9e0e-4089-add1-5abbde81cce1
# let
#     f = Figure()
#     ax1 = Axis(f[1, 1])

#     scatter!([0, density_array...], intensity_array)
#     lines!([0, density_array...], linearFit; color=:red)
#     ax1.title = "Calibration Line (Intensity vs Density)"
#     ax1.ylabel = "Intensity (HU)"
#     ax1.xlabel = "Density (mg/cm^3)"

#     f
# end

# ╔═╡ 1e68e0e8-a605-4557-925a-501ba8b9a1e1
md"""
# Score Background
"""

# ╔═╡ cd5770d8-6e31-4842-a350-e506d23ee77d
begin
	background_mask = zeros(size(arr)...)
	background_mask[
		(center_insert[1]-5):(center_insert[1]+5),
		(center_insert[2]-5):(center_insert[2]+5),
		2,
	] .= 1
	
	dilated_mask_L_bkg = dilate_mask_large_bkg(Bool.(background_mask))
	ring_mask_L_bkg = ring_mask_large(dilated_mask_L_bkg)

	dilated_mask_M_bkg = dilate_mask_medium_bkg(Bool.(background_mask))
	ring_mask_M_bkg = ring_mask_medium(dilated_mask_M_bkg)

	dilated_mask_S_bkg = dilate_mask_small_bkg(Bool.(background_mask))
	ring_mask_S_bkg = ring_mask_small(dilated_mask_S_bkg)
end;

# ╔═╡ 9bafed48-5d63-4512-9ea7-bd17dc2e0737
hu_heart_tissue_large_bkg = mean(arr[ring_mask_L_bkg])

# ╔═╡ 6c3f0dc5-9db2-4db8-abb0-64e4cda6d6fe
begin
	# if VENDOR == "80"
	# 	intensity_array = [0, 69.7, 377.3, 1375.3]
	# elseif VENDOR == "100"
	# 	intensity_array = [0, 63.3, 326.0, 1179.4]
	# elseif VENDOR == "120"
	# 	intensity_array = [0, 58.1, 296.9, 1072.1]
	# else
	# 	intensity_array = [0, 55.5, 282.7, 1019.7]
	# end

	if VENDOR == "80"
		intensity_array = [23, 69.7, 377.3, 1375.3]
	elseif VENDOR == "100"
		intensity_array = [23, 63.3, 326.0, 1179.4]
	elseif VENDOR == "120"
		intensity_array = [20, 58.1, 296.9, 1072.1]
	else
		intensity_array = [20, 55.5, 282.7, 1019.7]
	end

	# df_cal = DataFrame(:density => [0, 0.025, 0.200, 0.800], :intensity => intensity_array)
	df_cal = DataFrame(:density => [0, 0.2], :intensity => [hu_heart_tissue_large_bkg, hu_calcium])

	linearRegressor = lm(@formula(intensity ~ density), df_cal)
	linearFit = GLM.predict(linearRegressor)
	m = linearRegressor.model.pp.beta0[2]
	b = linearRegressor.model.rr.mu[1]
	density(intensity) = (intensity - b) / m
	intensity(ρ) = m * ρ + b
end

# ╔═╡ fc20bc4c-beb6-4979-be77-00bc6f481021
begin
	# Score background
	S_Obj = intensity(ρ_calcium)
	ρ = ρ_calcium # mg/mm^3

	S_bkg_large = mean(arr[ring_mask_L_bkg])
	alg_bkg_large = Integrated(arr[dilated_mask_L_bkg])
	mass_large_bkg = score(S_bkg_large, S_Obj, pixel_size, ρ, alg_bkg_large)

	S_bkg_medium = mean(arr[ring_mask_M_bkg])
	alg_bkg_medium = Integrated(arr[dilated_mask_M_bkg])
	mass_medium_bkg = score(S_bkg_medium, S_Obj, pixel_size, ρ, alg_bkg_medium)

	S_bkg_small = mean(arr[ring_mask_S_bkg])
	alg_bkg_small = Integrated(arr[dilated_mask_S_bkg])
	mass_small_bkg = score(S_bkg_small, S_Obj, pixel_size, ρ, alg_bkg_small)

	mass_bkg = [mass_large_bkg, mass_medium_bkg, mass_small_bkg]
end

# ╔═╡ 3b30e707-02b0-49c2-a0b3-25b7fcbf79ee
md"""
# Score Large Inserts
"""

# ╔═╡ 2a6961fc-f94f-4250-aeb8-63bbfd9e29a0
md"""
## High Density
"""

# ╔═╡ 58172907-2026-4caa-a6ad-ab484b86b329
md"""
#### Dilated mask
"""

# ╔═╡ 080422ff-77a2-4b86-8a0f-c2bea1f0d722
@bind g2 overlay_mask_bind(dilated_mask_L_HD)

# ╔═╡ 1cba2d1a-19c2-49b7-8757-5843bb16a0f1
overlay_mask_plot(arr, dilated_mask_L_HD, g2, "dilated mask")

# ╔═╡ cf4bc01d-5603-4d47-8863-e413be42f223
md"""
#### Ring (background) mask
"""

# ╔═╡ b990d904-e359-45fe-990b-87ab625ba1b0
@bind g4 overlay_mask_bind(ring_mask_L_HD)

# ╔═╡ 0af76504-ffc9-44ac-aecd-a3b956166273
overlay_mask_plot(arr, ring_mask_L_HD, g4, "ring mask")

# ╔═╡ 27a3c119-cd91-4172-94ea-f92400ab4c83
function score_test(vol, S_Bkg, S_Obj, size, ρ)
    I = Float64.(sum(vol))
    N = Float64.(length(vol))

    N_Obj = (I - (N * S_Bkg)) / (S_Obj - S_Bkg)
    Vol_Obj = N_Obj * size[1] * size[2] * size[3]
    Mass_Obj = Vol_Obj * ρ
    return Mass_Obj
end

# ╔═╡ dc94ad4a-0692-4e0e-84e5-f706ab1ca316
begin
	s_bkg_L_HD = mean(arr[ring_mask_L_HD])
	alg_L_HD = Integrated(arr[mask_L_HD_3D])
	mass_l_hd = score(s_bkg_L_HD, S_Obj, pixel_size, ρ, alg_L_HD)
end

# ╔═╡ d73e2a4c-a4fa-480a-82e1-2e9519490045
s_bkg_L_HD

# ╔═╡ 9f69e44f-c927-4664-9cd8-0a0a0cf70c7c
mass_l_hd_test = score_test(arr[mask_L_HD_3D], s_bkg_L_HD, S_Obj, pixel_size, ρ)

# ╔═╡ e93b26b5-ba84-4410-87b7-eefad3d8ab3e
md"""
## Medium Density
"""

# ╔═╡ 336cf5d0-ada0-4f9a-9788-bfd0ab72087a
md"""
#### Dilated mask
"""

# ╔═╡ 937a5188-e51e-4017-be9a-2d8051e79fd1
@bind h2 overlay_mask_bind(dilated_mask_L_MD)

# ╔═╡ f4ee897f-d6ad-4dfd-a91f-a759ecd170e5
overlay_mask_plot(arr, dilated_mask_L_MD, h2, "dilated mask")

# ╔═╡ afbc96dc-f5fa-4305-bd3d-119ef8d416a8
md"""
#### Ring (background) mask
"""

# ╔═╡ 8e310e91-c557-48e5-a5d1-042c1ac8f7b3
@bind h4 overlay_mask_bind(ring_mask_L_MD)

# ╔═╡ ed5c5da5-21dc-4052-bc18-6b209edeabfd
overlay_mask_plot(arr, ring_mask_L_MD, h4, "ring mask")

# ╔═╡ 367750d8-a287-4313-b6df-d6533849623f
begin
	s_bkg_L_MD = mean(arr[ring_mask_L_MD])
	alg_L_MD = Integrated(arr[mask_L_MD_3D])
	mass_l_md = score(s_bkg_L_MD, S_Obj, pixel_size, ρ, alg_L_MD)
end

# ╔═╡ 73105b5a-05b7-429f-9a7f-b8c12684df91
md"""
## Low Density
"""

# ╔═╡ d8d7a93f-85e7-43e7-a23d-11e7c6ce27be
md"""
#### Dilated mask
"""

# ╔═╡ b849a64f-965a-40b7-b2f5-1a27fd8ceeb5
@bind i2 overlay_mask_bind(dilated_mask_L_LD)

# ╔═╡ 7f52b670-611d-4bd3-bc80-77ba59aa2d70
overlay_mask_plot(arr, dilated_mask_L_LD, i2, "dilated mask")

# ╔═╡ 4c425cf3-2ece-4df1-b3ab-4058c7d7aa13
md"""
#### Ring (background) mask
"""

# ╔═╡ 7a618b0c-9429-472a-b076-774b9d229cab
@bind i4 overlay_mask_bind(ring_mask_L_LD)

# ╔═╡ 544151a5-3b47-4571-a6f0-9800d5b11085
overlay_mask_plot(arr, ring_mask_L_LD, i4, "ring mask")

# ╔═╡ 7560dc46-fa36-4706-925c-0675397c5922
begin
	s_bkg_L_LD = mean(arr[ring_mask_L_LD])
	alg_L_LD = Integrated(arr[mask_L_LD_3D])
	mass_l_ld = score(s_bkg_L_LD, S_Obj, pixel_size, ρ, alg_L_LD)
end

# ╔═╡ ee7092e6-52fb-4a6a-ad21-90cd2b46c7fa
md"""
# Score Medium Inserts
"""

# ╔═╡ 62faf81f-eb3b-431a-8e2b-239735db3257
md"""
## High Density
"""

# ╔═╡ 10765555-8507-413f-a582-1f47619f0181
md"""
#### Dilated mask
"""

# ╔═╡ 5aff6e7a-b635-44d7-a91e-1a30c5f4fd17
@bind j2 overlay_mask_bind(dilated_mask_M_HD)

# ╔═╡ 4aabe7d9-e730-4b4f-9dd2-e51e922ac08c
overlay_mask_plot(arr, dilated_mask_M_HD, j2, "dilated mask")

# ╔═╡ 706c30df-ce63-4c81-9c99-227d56ab2991
md"""
#### Ring (background) mask
"""

# ╔═╡ 47615958-a8fc-4fc2-852f-e11f69041a86
@bind j4 overlay_mask_bind(ring_mask_M_HD)

# ╔═╡ 85097de3-d336-4c85-a759-ad1b3e4083cc
overlay_mask_plot(arr, ring_mask_M_HD, j4, "ring mask")

# ╔═╡ bd4cec3c-446f-48e4-b06b-fbde4848daaf
begin
	s_bkg_M_HD = mean(arr[ring_mask_M_HD])
	alg_M_HD = Integrated(arr[mask_M_HD_3D])
	mass_m_hd = score(s_bkg_M_HD, S_Obj, pixel_size, ρ, alg_M_HD)
end

# ╔═╡ 0f71036e-b118-430d-bce7-420b2591f1a2
md"""
## Medium Density
"""

# ╔═╡ c99672f4-9ee7-4d42-8a15-eee68be2a8da
md"""
#### Dilated mask
"""

# ╔═╡ b1c10feb-3148-4de5-b305-b0721894b76e
@bind k2 overlay_mask_bind(dilated_mask_M_MD)

# ╔═╡ bd91bd0d-f368-4340-916a-ea197969c0e4
overlay_mask_plot(arr, dilated_mask_M_MD, k2, "dilated mask")

# ╔═╡ b5cb8919-df84-48a6-b00b-cf7dba301c76
md"""
#### Ring (background) mask
"""

# ╔═╡ 6c3c8eb5-c651-45e3-a912-8fd0dc330b51
@bind k4 overlay_mask_bind(ring_mask_M_MD)

# ╔═╡ f2bafa39-c1ac-42ca-89c0-df4d219c2c74
overlay_mask_plot(arr, ring_mask_M_MD, k4, "ring mask")

# ╔═╡ 43afbb30-b1b5-452b-8c82-25d173f24699
begin
	s_bkg_M_MD = mean(arr[ring_mask_M_MD])
	alg_M_MD = Integrated(arr[mask_M_MD_3D])
	mass_m_md = score(s_bkg_M_MD, S_Obj, pixel_size, ρ, alg_M_MD)
end

# ╔═╡ 386ea37d-70ba-4530-be9e-549d5e0ab1ce
md"""
## Low Density
"""

# ╔═╡ b37af6a8-8f51-4263-ac85-bea25fd2aae2
md"""
#### Dilated mask
"""

# ╔═╡ 7b8af8e8-e979-4b57-a65a-ad99a40f08b9
@bind l2 overlay_mask_bind(dilated_mask_M_LD)

# ╔═╡ 493ef94e-d6ac-4ccb-b342-0812792c2536
overlay_mask_plot(arr, dilated_mask_M_LD, l2, "dilated mask")

# ╔═╡ 43f5bf4e-c882-4596-a4b5-6063ce63e17c
md"""
#### Ring (background) mask
"""

# ╔═╡ 48a0f2c4-0bea-4c7d-bf8a-dc3632a490fa
@bind l4 overlay_mask_bind(ring_mask_M_LD)

# ╔═╡ e03ccab0-b786-49ca-b7dd-b2b7f76e04a2
overlay_mask_plot(arr, ring_mask_M_LD, l4, "ring mask")

# ╔═╡ 907b3f24-47b8-4519-b70a-de44ca764bb2
begin
	s_bkg_M_LD = mean(arr[ring_mask_M_LD])
	alg_M_LD = Integrated(arr[mask_M_LD_3D])
	mass_m_ld = score(s_bkg_M_LD, S_Obj, pixel_size, ρ, alg_M_LD)
end

# ╔═╡ c52478ca-2012-4690-bc39-6d32997427f7
md"""
# Score Small Inserts
"""

# ╔═╡ 5de67b93-72d7-4f86-b94e-245f4cfff576
md"""
## High Density
"""

# ╔═╡ d0ea5ea6-5c98-44b8-94d4-4104eb1f6080
md"""
#### Dilated mask
"""

# ╔═╡ 018d2392-7133-4da8-b4c9-268af6fa6f3a
@bind m2 overlay_mask_bind(dilated_mask_S_HD)

# ╔═╡ 9d72b647-73d9-45c7-b473-bbdea2817cb5
overlay_mask_plot(arr, dilated_mask_S_HD, m2, "dilated mask")

# ╔═╡ 27ffaff1-12c5-43f1-ac21-84fe25877f24
md"""
#### Ring (background) mask
"""

# ╔═╡ 46110aad-07a8-4a56-9a6d-db71054a5883
@bind m4 overlay_mask_bind(ring_mask_S_HD)

# ╔═╡ 070aedea-d5e9-43f7-a2cf-14fdad760a72
overlay_mask_plot(arr, ring_mask_S_HD, m4, "ring mask")

# ╔═╡ 1590d3c5-d9f1-4c00-9730-639677635056
begin
	s_bkg_S_HD = mean(arr[ring_mask_S_HD])
	alg_S_HD = Integrated(arr[mask_S_HD_3D])
	mass_s_hd = score(s_bkg_S_HD, S_Obj, pixel_size, ρ, alg_S_HD)
end

# ╔═╡ aa3f6b52-2c1a-462a-a89c-bcde6f0f59ba
md"""
## Medium Density
"""

# ╔═╡ 01dce72f-6793-4360-80fc-bc2e97d579de
md"""
#### Dilated mask
"""

# ╔═╡ ed6c5395-2780-47d2-a3ca-322468638828
@bind n2 overlay_mask_bind(dilated_mask_S_MD)

# ╔═╡ f5599896-297e-4cac-89b0-2dd1344f665c
overlay_mask_plot(arr, dilated_mask_S_MD, n2, "dilated mask")

# ╔═╡ 64247b3c-eef8-4d19-9bf9-71cd10f2fc02
md"""
#### Ring (background) mask
"""

# ╔═╡ b7fba435-4496-4671-a377-8b73b03949d4
@bind n4 overlay_mask_bind(ring_mask_S_MD)

# ╔═╡ 386ec444-7e85-458e-a213-d2f5e87374d3
overlay_mask_plot(arr, ring_mask_S_MD, n4, "ring mask")

# ╔═╡ dcd5144e-408f-4b13-995e-f3e716f27c3b
begin
	s_bkg_S_MD = mean(arr[ring_mask_S_MD])
	alg_S_MD = Integrated(arr[mask_S_MD_3D])
	mass_s_md = score(s_bkg_S_MD, S_Obj, pixel_size, ρ, alg_S_MD)
end

# ╔═╡ 5d625a7e-f027-443b-aa97-726350616c11
md"""
## Low Density
"""

# ╔═╡ 776fa172-d964-4c86-8220-725c11f6768c
md"""
#### Dilated mask
"""

# ╔═╡ d2dbe7a4-fc39-406b-9ad2-313b191dcf63
@bind o2 overlay_mask_bind(dilated_mask_S_LD)

# ╔═╡ cdbd4ccf-1fc4-4ff8-af3a-847113b85ef2
overlay_mask_plot(arr, dilated_mask_S_LD, o2, "dilated mask")

# ╔═╡ 79cf02a6-98f4-4d35-85a4-323c5d160c8b
md"""
#### Ring (background) mask
"""

# ╔═╡ 4969b058-f3f8-4916-9ed9-c3b99f23f972
@bind o4 overlay_mask_bind(ring_mask_S_LD)

# ╔═╡ 9c3c3eb7-17b9-440d-b694-987f65aa0513
overlay_mask_plot(arr, ring_mask_S_LD, o4, "ring mask")

# ╔═╡ b7f78d93-dd53-45b5-b82d-4ffe3f8ffc03
begin
	s_bkg_S_LD = mean(arr[ring_mask_S_LD])
	alg_S_LD = Integrated(arr[mask_S_LD_3D])
	mass_s_ld = score(s_bkg_S_LD, S_Obj, pixel_size, ρ, alg_S_LD)
end

# ╔═╡ 7ba16b69-a495-4c18-bb10-d61ac48cde4b
md"""
# Results
"""

# ╔═╡ b43b3658-fa22-4380-b218-be26f0dace70
volume = (pi * 0.5^2) * 9

# ╔═╡ b10419ba-4172-497e-b42a-4937a86c9ea2
begin
	inserts = ["Low Density", "Medium Density", "High Density"]
	volume_gt = [7.065, 63.585, 176.625]
	ground_truth_mass_large = [
		volume_gt[3] * density_array[1],
		volume_gt[3] * density_array[2],
		volume_gt[3] * density_array[3],
	] # mg
	ground_truth_mass_medium = [
		volume_gt[2] * density_array[1],
		volume_gt[2] * density_array[2],
		volume_gt[2] * density_array[3],
	] # mg
	ground_truth_mass_small = [
		volume_gt[1] * density_array[1],
		volume_gt[1] * density_array[2],
		volume_gt[1] * density_array[3],
	] # mg
end

# ╔═╡ ef3aca58-e313-443b-9c6c-c4f6d2cfd86f
begin
	calculated_mass_large = [mass_l_ld, mass_l_md, mass_l_hd]
	calculated_mass_medium = [mass_m_ld, mass_m_md, mass_m_hd]
	calculated_mass_small = [mass_s_ld, mass_s_md, mass_s_hd]
end

# ╔═╡ 82540ec2-d054-40db-a561-8d7ecb254756
df = DataFrame(;
	VENDOR=VENDOR,
	SIZE=SIZE,
	DENSITY=DENSITY,
	scan=scan,
	inserts=inserts,
	ground_truth_mass_large=ground_truth_mass_large,
	calculated_mass_large=calculated_mass_large,
	ground_truth_mass_medium=ground_truth_mass_medium,
	calculated_mass_medium=calculated_mass_medium,
	ground_truth_mass_small=ground_truth_mass_small,
	calculated_mass_small=calculated_mass_small,
	mass_bkg=mass_bkg
)

# ╔═╡ Cell order:
# ╠═45f704d4-66e5-49db-aef5-05132f3853ee
# ╠═e7c4aaab-83ef-4256-b392-f8ef7a899a05
# ╟─48e3097f-0767-4dad-9d7c-94e0899f790b
# ╠═bc9383c0-e477-4e1a-a2fa-7f5c1d29f103
# ╟─66f5fc66-d5ba-45ba-8624-55adc58085e4
# ╟─6e812172-6371-4461-9365-22f68ef16e53
# ╟─483a14dd-e798-41ed-9144-13678f8b8461
# ╟─2cfdeb06-71a2-411b-accb-2b4a3f56b477
# ╟─ab294a53-0ac8-4a40-8a98-c5271c66d958
# ╟─3f12b12c-90a9-48cc-99db-3545e4d0a617
# ╟─29bb75c2-ff32-4a01-b04c-614647c66c67
# ╟─eacd75fc-1a1f-4f55-81c6-f897154adc38
# ╟─6dfd0583-21c2-400a-9ab5-3a09a95e5adc
# ╟─6d306954-6c37-4138-924d-151f0b2dd068
# ╟─8678df23-ef1c-42b0-8102-84ca44f4f7d8
# ╟─e4d10b8b-93fb-415e-bd35-aca8fee7510e
# ╟─d6a34d2d-937f-4f59-9e77-c8d484bbe75b
# ╟─2985cf19-4a56-478f-9eda-5387b9fbec08
# ╟─b3bae2ad-16ed-4381-b6a4-448cb5f0c6c1
# ╠═4ff19af0-7f53-4f24-b315-08bd54a488e3
# ╟─b74941c4-12d4-4be1-81f7-ef1f4d288984
# ╟─84e56c91-14a4-45b7-81a1-77e778bae695
# ╟─4a4cf730-1307-451e-be5c-f9c1183cdd5f
# ╟─abfd965c-df00-4028-a0f3-8356a261b527
# ╟─41541a09-609b-4f46-9a25-eb974422efc0
# ╟─e6753c9a-d172-47fe-b066-78cb43df2c89
# ╟─f8a4eacd-2f69-4755-b017-5a0b0dda1004
# ╠═ac9d2652-6d69-4cf3-acbf-2e141fb633f0
# ╟─417d150e-0df9-4963-b59e-ebd9acf6d4a0
# ╠═5f10075b-4d95-4d53-a22f-016749fb7583
# ╟─b921dcf8-54ea-420b-9285-23e38d5ce433
# ╠═05ed60db-b4d8-4cdc-9a54-b108ade22557
# ╠═8bb2d451-b575-40bf-b31e-51e80392ef51
# ╠═4862496e-464d-40c3-baef-e5c507028d66
# ╟─1545510f-1753-4c60-88d0-75606e4903f1
# ╠═16fb4a53-116a-404c-95ae-bedcadbc8f9f
# ╠═ccec799c-e328-4c41-8279-c71f27bd92ec
# ╟─997952cd-dc01-4528-a30a-fdc11a46dad8
# ╠═f3f73d41-0d98-461b-8673-7207857253ad
# ╠═da3c96b9-efdb-40f4-a959-747c097e16b2
# ╠═5fec186d-82bb-4b37-b0c6-c3781f2a166d
# ╠═9bafed48-5d63-4512-9ea7-bd17dc2e0737
# ╠═6c3f0dc5-9db2-4db8-abb0-64e4cda6d6fe
# ╟─e50baf62-0a02-47ee-8440-551f68baee0d
# ╟─b03c711f-9e0e-4089-add1-5abbde81cce1
# ╟─1e68e0e8-a605-4557-925a-501ba8b9a1e1
# ╠═cd5770d8-6e31-4842-a350-e506d23ee77d
# ╠═fc20bc4c-beb6-4979-be77-00bc6f481021
# ╟─3b30e707-02b0-49c2-a0b3-25b7fcbf79ee
# ╟─2a6961fc-f94f-4250-aeb8-63bbfd9e29a0
# ╟─58172907-2026-4caa-a6ad-ab484b86b329
# ╟─080422ff-77a2-4b86-8a0f-c2bea1f0d722
# ╠═1cba2d1a-19c2-49b7-8757-5843bb16a0f1
# ╟─cf4bc01d-5603-4d47-8863-e413be42f223
# ╟─b990d904-e359-45fe-990b-87ab625ba1b0
# ╠═0af76504-ffc9-44ac-aecd-a3b956166273
# ╠═27a3c119-cd91-4172-94ea-f92400ab4c83
# ╠═d73e2a4c-a4fa-480a-82e1-2e9519490045
# ╠═dc94ad4a-0692-4e0e-84e5-f706ab1ca316
# ╠═9f69e44f-c927-4664-9cd8-0a0a0cf70c7c
# ╟─e93b26b5-ba84-4410-87b7-eefad3d8ab3e
# ╟─336cf5d0-ada0-4f9a-9788-bfd0ab72087a
# ╟─937a5188-e51e-4017-be9a-2d8051e79fd1
# ╠═f4ee897f-d6ad-4dfd-a91f-a759ecd170e5
# ╟─afbc96dc-f5fa-4305-bd3d-119ef8d416a8
# ╟─8e310e91-c557-48e5-a5d1-042c1ac8f7b3
# ╠═ed5c5da5-21dc-4052-bc18-6b209edeabfd
# ╠═367750d8-a287-4313-b6df-d6533849623f
# ╟─73105b5a-05b7-429f-9a7f-b8c12684df91
# ╟─d8d7a93f-85e7-43e7-a23d-11e7c6ce27be
# ╟─b849a64f-965a-40b7-b2f5-1a27fd8ceeb5
# ╠═7f52b670-611d-4bd3-bc80-77ba59aa2d70
# ╟─4c425cf3-2ece-4df1-b3ab-4058c7d7aa13
# ╟─7a618b0c-9429-472a-b076-774b9d229cab
# ╠═544151a5-3b47-4571-a6f0-9800d5b11085
# ╠═7560dc46-fa36-4706-925c-0675397c5922
# ╟─ee7092e6-52fb-4a6a-ad21-90cd2b46c7fa
# ╟─62faf81f-eb3b-431a-8e2b-239735db3257
# ╟─10765555-8507-413f-a582-1f47619f0181
# ╟─5aff6e7a-b635-44d7-a91e-1a30c5f4fd17
# ╠═4aabe7d9-e730-4b4f-9dd2-e51e922ac08c
# ╟─706c30df-ce63-4c81-9c99-227d56ab2991
# ╟─47615958-a8fc-4fc2-852f-e11f69041a86
# ╠═85097de3-d336-4c85-a759-ad1b3e4083cc
# ╠═bd4cec3c-446f-48e4-b06b-fbde4848daaf
# ╟─0f71036e-b118-430d-bce7-420b2591f1a2
# ╟─c99672f4-9ee7-4d42-8a15-eee68be2a8da
# ╟─b1c10feb-3148-4de5-b305-b0721894b76e
# ╠═bd91bd0d-f368-4340-916a-ea197969c0e4
# ╟─b5cb8919-df84-48a6-b00b-cf7dba301c76
# ╟─6c3c8eb5-c651-45e3-a912-8fd0dc330b51
# ╠═f2bafa39-c1ac-42ca-89c0-df4d219c2c74
# ╠═43afbb30-b1b5-452b-8c82-25d173f24699
# ╟─386ea37d-70ba-4530-be9e-549d5e0ab1ce
# ╟─b37af6a8-8f51-4263-ac85-bea25fd2aae2
# ╟─7b8af8e8-e979-4b57-a65a-ad99a40f08b9
# ╠═493ef94e-d6ac-4ccb-b342-0812792c2536
# ╟─43f5bf4e-c882-4596-a4b5-6063ce63e17c
# ╟─48a0f2c4-0bea-4c7d-bf8a-dc3632a490fa
# ╠═e03ccab0-b786-49ca-b7dd-b2b7f76e04a2
# ╠═907b3f24-47b8-4519-b70a-de44ca764bb2
# ╟─c52478ca-2012-4690-bc39-6d32997427f7
# ╟─5de67b93-72d7-4f86-b94e-245f4cfff576
# ╟─d0ea5ea6-5c98-44b8-94d4-4104eb1f6080
# ╟─018d2392-7133-4da8-b4c9-268af6fa6f3a
# ╠═9d72b647-73d9-45c7-b473-bbdea2817cb5
# ╟─27ffaff1-12c5-43f1-ac21-84fe25877f24
# ╟─46110aad-07a8-4a56-9a6d-db71054a5883
# ╠═070aedea-d5e9-43f7-a2cf-14fdad760a72
# ╠═1590d3c5-d9f1-4c00-9730-639677635056
# ╟─aa3f6b52-2c1a-462a-a89c-bcde6f0f59ba
# ╟─01dce72f-6793-4360-80fc-bc2e97d579de
# ╟─ed6c5395-2780-47d2-a3ca-322468638828
# ╠═f5599896-297e-4cac-89b0-2dd1344f665c
# ╟─64247b3c-eef8-4d19-9bf9-71cd10f2fc02
# ╟─b7fba435-4496-4671-a377-8b73b03949d4
# ╠═386ec444-7e85-458e-a213-d2f5e87374d3
# ╠═dcd5144e-408f-4b13-995e-f3e716f27c3b
# ╟─5d625a7e-f027-443b-aa97-726350616c11
# ╟─776fa172-d964-4c86-8220-725c11f6768c
# ╟─d2dbe7a4-fc39-406b-9ad2-313b191dcf63
# ╠═cdbd4ccf-1fc4-4ff8-af3a-847113b85ef2
# ╟─79cf02a6-98f4-4d35-85a4-323c5d160c8b
# ╟─4969b058-f3f8-4916-9ed9-c3b99f23f972
# ╠═9c3c3eb7-17b9-440d-b694-987f65aa0513
# ╠═b7f78d93-dd53-45b5-b82d-4ffe3f8ffc03
# ╟─7ba16b69-a495-4c18-bb10-d61ac48cde4b
# ╠═b43b3658-fa22-4380-b218-be26f0dace70
# ╠═b10419ba-4172-497e-b42a-4937a86c9ea2
# ╠═ef3aca58-e313-443b-9c6c-c4f6d2cfd86f
# ╠═82540ec2-d054-40db-a561-8d7ecb254756

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

# ╔═╡ 9860656e-409f-4dc7-bd95-e96781a104f0
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.activate(".")

    using PlutoUI, Statistics, CSV, DataFrames, GLM, CairoMakie, HypothesisTests, Colors, MLJBase, DICOM, DICOMUtils, PhantomSegmentation, CalciumScoring, ImageMorphology, ImageFiltering, Noise
    using StatsBase: quantile!, rmsd
end

# ╔═╡ 163ffa5f-0509-44db-a5e5-1f93d625079b
TableOfContents()

# ╔═╡ 01526785-a645-4b1a-afc2-0bc828c2db3a
md"""
## Load DICOMS

All you need to do is set `base_path` once and leave it. After that, the only thing that should change is the `VENDOR`, once for every set, and the `SCAN_NUMBER`, once for each scan.
"""

# ╔═╡ 0a2c9749-e2bd-4d5e-9ef5-a07834f6684b
begin
	SAVE_DF = "physical.csv"
	VENDORS = ["Canon_Aquilion_One_Vision"]
	# SCANS = [1, 2, 6, 8]
	SCANS = [1]

	VENDOR = VENDORS[1]
	SCAN = SCANS[1]
	
	BASE_PATH = "/Users/daleblack/Google Drive/Datasets/"
	root_path = joinpath(BASE_PATH, VENDOR)
	dcm_path_list = dcm_list_builder(root_path)
	pth = dcm_path_list[SCAN]
	scan = basename(pth)
	header, dcm_array, slice_thick_ori1 = dcm_reader(pth)
	kV = header[tag"KVP"]
end

# ╔═╡ 28a061b7-5a67-4454-80d8-d8db807f2d33
md"""
## Helper Functions
"""

# ╔═╡ f07b68cb-97fe-4073-8363-55bbf4e49a3f
function collect_tuple(tuple_array)
    row_num = size(tuple_array)
    col_num = length(tuple_array[1])
    container = zeros(Int64, row_num..., col_num)
    for i in 1:length(tuple_array)
        container[i, :] = collect(tuple_array[i])
    end
    return container
end

# ╔═╡ 5ae40e25-6c2b-4a69-ae75-71e022efc50f
function overlay_mask_bind(mask)
    indices = findall(x -> x == 1, mask)
    indices = Tuple.(indices)
    label_array = collect_tuple(indices)
    zs = unique(label_array[:, 3])
    return PlutoUI.Slider(1:length(zs); default=3, show_value=true)
end

# ╔═╡ 8ed6a20b-e2d9-4596-a296-8e0dccbd7e5b
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

# ╔═╡ 37b0192b-a590-421f-877a-339d96493d08
function create_mask(array, mask)
    @assert size(array) == size(mask)
    idxs = findall(x -> x == true, mask)
    overlayed_mask = zeros(size(array))
    for idx in idxs
        overlayed_mask[idx] = array[idx]
    end
    return overlayed_mask
end

# ╔═╡ 0fdfc1be-f7ee-4d93-8186-2e0a3e121366
function dilate_mask_large(mask)
    return dilate(mask)
end

# ╔═╡ 4fea86c3-2c9c-4549-b4c6-b67ffb6bad04
function ring_mask_large(dilated_mask)
    return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask)
end

# ╔═╡ 1e773785-94a6-490a-99d3-5da4478fa9a0
function dilate_mask_medium(mask)
    return dilate(dilate(dilate(mask)))
end

# ╔═╡ 8a68d139-8712-40ef-a2c5-ee33ba278174
function ring_mask_medium(dilated_mask)
    return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask)
end

# ╔═╡ 593b9d85-dcf7-4db4-9e78-ca284007d30e
function dilate_mask_small(mask)
    return dilate(dilate(dilate(mask)))
end

# ╔═╡ 94200cb0-0c4e-4fa8-8ed9-418a9aa25469
function ring_mask_small(dilated_mask)
    return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask)
end

# ╔═╡ 509af748-95e8-4122-9f81-aecc42c33901
md"""
## Segment Heart
"""

# ╔═╡ 9a93a6b5-9ec6-4667-b9b0-ad0f6a58709f
masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3) ÷ 2);

# ╔═╡ cd4c2e05-c305-4da4-a411-054bc64a810c
@bind a PlutoUI.Slider(1:size(masked_array, 3), default=10, show_value=true)

# ╔═╡ 09a6b6e2-9aee-4d77-b101-8e143a5d7fb3
heatmap(transpose(masked_array[:, :, a]); colormap=:grays)

# ╔═╡ 97d7830f-befe-4f4e-b5d7-a37c9cf690ff
let
    fig = Figure()

    ax = Makie.Axis(fig[1, 1])
    heatmap!(transpose(dcm_array[:, :, 15]); colormap=:grays)
	
    fig
end

# ╔═╡ c52cf513-a741-4fe6-b8a2-117a1059c4fa
let
    fig = Figure()

    ax = Makie.Axis(fig[1, 1])
    ax.title = "Raw DICOM Array"
    heatmap!(transpose(dcm_array[:, :, 15]); colormap=:grays)
    scatter!(
        center_insert[2]:center_insert[2],
        center_insert[1]:center_insert[1];
        markersize=10,
        color=:red,
    )
    fig
end

# ╔═╡ 6d0b5f0c-5a9b-470d-acbe-07bf5bafe27d
let
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

# ╔═╡ 60d95a36-f3d1-448f-a0cb-d92eaf9808a0
let
    fig3 = Figure()

    ax3 = Makie.Axis(fig3[1, 1])
    ax3.title = "Masked DICOM Array"
    heatmap!(transpose(masked_array[:, :, 15]); colormap=:grays)
    scatter!(
        center_insert[2]:center_insert[2],
        center_insert[1]:center_insert[1];
        markersize=10,
        color=:red,
    )
    fig3
end

# ╔═╡ 63ed1025-518a-4ac0-b6be-c87fca2a53ca
md"""
## Segment Calcium Rod
"""

# ╔═╡ 0db64444-a5c6-4c42-b129-e6f24f2ce426
begin
	thresh = 130
	calcium_image, slice_CCI, quality_slice, cal_rod_slice = mask_rod(
		masked_array, header; calcium_threshold=thresh
	)
end;

# ╔═╡ 9cad1e1c-9e66-4e5a-86c6-9b02b0a0e30f
@bind c PlutoUI.Slider(axes(calcium_image, 3), default=5, show_value=true)

# ╔═╡ e4dccc7e-d47c-41a5-9ffd-d535b1615f07
heatmap(transpose(calcium_image[:, :, c]); colormap=:grays)

# ╔═╡ 6af669a9-77cd-480c-a7b6-1d9532ce62e1
md"""
## Load Masks
"""

# ╔═╡ 0fbf3a06-0850-444e-aa06-a3508661d0a4
mask_L_HD, mask_M_HD, mask_S_HD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_LD, mask_M_LD, mask_S_LD = mask_inserts(
	dcm_array, masked_array, header, slice_CCI, center_insert;
	calcium_threshold=thresh
);

# ╔═╡ 807b06fc-5246-4c8e-b7b1-389cdbda395c
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

# ╔═╡ 5aa9dd7b-82a5-40e0-a8ce-4c4079369c8f
heatmap(transpose(masks); colormap=:grays)

# ╔═╡ bdf43911-b1f3-4567-bb84-6781dbc320cc
md"""
## Prepare Array Masks
"""

# ╔═╡ 8a978945-b378-4ed4-8aab-3e6e5451a298
begin
	arr = masked_array[:, :, slice_CCI-2:slice_CCI+2]
	pixel_size = DICOMUtils.get_pixel_size(header)
	voxel_size = pixel_size[1] * pixel_size[2] * pixel_size[3]
end;

# ╔═╡ 0a78b66a-a545-409e-bc0b-5f575e6dea36
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

# ╔═╡ 2346e495-e1d4-4955-8ff0-cdff7fb07032
md"""
## Calibration Prep
"""

# ╔═╡ 6173c719-886d-4142-b2e3-cee8f9babdc6
begin
	density_array = [0.200, 0.400, 0.800]
	array_filtered = abs.(mapwindow(median, calcium_image[:, :, cal_rod_slice], (3, 3)))
	bool_arr = array_filtered .> 0
	bool_arr_erode = (erode(erode(erode(erode(bool_arr)))))
	c_img = calcium_image[:, :, cal_rod_slice-1:cal_rod_slice+2]
	mask_cal_3D = zeros(size(c_img))
	for z in 1:size(c_img, 3)
		mask_cal_3D[:, :, z] = Bool.(erode(bool_arr_erode))
	end

	hu_calcium = mean(c_img[Bool.(mask_cal_3D)])
	# hu_calcium = 296.9
	ρ_calcium = 0.2
end

# ╔═╡ 1cd4003d-02bb-4c6b-8bf0-0b5e1299d148
md"""
# Background
"""

# ╔═╡ 072167d0-65e9-427f-991c-6e3aa23c172c
begin
	background_mask = zeros(size(arr)...)
	background_mask[
		(center_insert[1]-5):(center_insert[1]+5),
		(center_insert[2]-5):(center_insert[2]+5),
		2,
	] .= 1
	background_mask = Bool.(dilate(dilate(background_mask)))
	background_ring = ring_mask_large(background_mask)
	hu_heart_tissue_bkg = mean(arr[background_ring])
	mass_bkg = score(arr[background_mask], hu_calcium, hu_heart_tissue_bkg, voxel_size, ρ_calcium, VolumeFraction())
end

# ╔═╡ 97dbd4e9-ac4a-4f4e-aa2b-de123102bec1
hu_heart_tissue_bkg

# ╔═╡ 6a01dbc6-d8ca-4cad-a83c-995442c0b6e5
size(background_mask)

# ╔═╡ e80a1e90-a89c-4530-88c0-dd2b80c7bd3f
@bind q1 overlay_mask_bind(background_mask)

# ╔═╡ c095470d-c1a8-46eb-99c9-65cb3dfffd34
overlay_mask_plot(arr, background_mask, q1, "dilated mask")

# ╔═╡ ff0b1b99-ec2d-4610-97f2-ccb68a854460
overlay_mask_plot(arr, background_ring, q1, "dilated mask")

# ╔═╡ 5aaabb67-d87b-4761-a927-2f53448ee697
md"""
# Score Large Inserts
"""

# ╔═╡ 7d06e368-c46c-451a-84c8-424fbf12f06b
md"""
## High Density
"""

# ╔═╡ 24d9ad56-8c07-4af9-b216-de9bdcd16fb3
md"""
#### Dilated mask
"""

# ╔═╡ 208f81ed-06e6-4eb0-9c50-9c487547edb9
@bind g2 overlay_mask_bind(dilated_mask_L_HD)

# ╔═╡ 6983213e-1c80-4738-9695-b83bab339207
overlay_mask_plot(arr, dilated_mask_L_HD, g2, "dilated mask")

# ╔═╡ bf1474f6-9c56-4b0a-8309-42430dd5ef7d
md"""
#### Ring (background) mask
"""

# ╔═╡ 04df393b-7a3f-48be-944c-902146303088
@bind g4 overlay_mask_bind(ring_mask_L_HD)

# ╔═╡ a1391510-889f-46f8-bf29-48a4d9266c3c
overlay_mask_plot(arr, ring_mask_L_HD, g4, "ring mask")

# ╔═╡ 6869fc6c-ca3d-4c44-9a47-9e964a5dbc02
hu_heart_tissue_large_hd = mean(arr[ring_mask_L_HD])

# ╔═╡ cecb93a2-d6fc-4b0d-acb7-c3e816d52498
mass_large_hd = score(arr[dilated_mask_L_HD], hu_calcium, hu_heart_tissue_bkg, voxel_size, ρ_calcium, VolumeFraction())

# ╔═╡ 6e93a14b-7c0a-4555-bd81-bbd002aaef50
md"""
## Medium Density
"""

# ╔═╡ 4d675576-16e4-4241-8548-e40593e04893
md"""
#### Dilated mask
"""

# ╔═╡ e9296b3e-3b28-4c8b-b44d-7ce65d6c3de7
@bind h2 overlay_mask_bind(dilated_mask_L_MD)

# ╔═╡ a144f068-6818-4a4d-8b61-87966d89c789
overlay_mask_plot(arr, dilated_mask_L_MD, h2, "dilated mask")

# ╔═╡ 30f68359-e661-4b74-8046-03cf8123a1a0
md"""
#### Ring (background) mask
"""

# ╔═╡ 94bfca4c-3068-4683-b8a2-239c7afc4d7d
@bind h4 overlay_mask_bind(ring_mask_L_MD)

# ╔═╡ ab58e65d-940c-4058-8a02-56d23ae3bca5
overlay_mask_plot(arr, ring_mask_L_MD, h4, "ring mask")

# ╔═╡ 911ddf88-0149-4146-8e08-25995de7d0d4
hu_heart_tissue_large_md = mean(arr[ring_mask_L_MD])

# ╔═╡ 97c24681-62d5-410e-a2ff-8cad19d13147
mass_large_md = score(arr[dilated_mask_L_MD], hu_calcium, hu_heart_tissue_large_md, voxel_size, ρ_calcium, VolumeFraction())

# ╔═╡ 7a133076-83c2-4543-8c2e-0cbd9d234b27
md"""
## Low Density
"""

# ╔═╡ c20e99b7-2a62-4afb-ad8a-cd2f1815a7ca
md"""
#### Dilated mask
"""

# ╔═╡ 89181d61-763a-47fc-910e-7bd5c52a2db6
@bind j2 overlay_mask_bind(dilated_mask_L_LD)

# ╔═╡ f41e74e1-01a7-4581-bc1d-b9592005adcb
overlay_mask_plot(arr, dilated_mask_L_LD, j2, "dilated mask")

# ╔═╡ 2b6a668b-ad5b-4677-b1b5-3ae0e309d2cf
md"""
#### Ring (background) mask
"""

# ╔═╡ 40035644-21a6-4b1e-ace9-845af80fb9e6
@bind j4 overlay_mask_bind(ring_mask_L_LD)

# ╔═╡ 39be734a-5f4d-43a2-89c7-c49fcd2ab49d
overlay_mask_plot(arr, ring_mask_L_LD, j4, "ring mask")

# ╔═╡ 0b81a43b-6f34-4955-966b-3198ee522304
hu_heart_tissue_large_ld = mean(arr[ring_mask_L_LD])

# ╔═╡ 9cd46d40-d35d-4fdc-8d3e-fc11b6072358
mass_large_ld = score(arr[dilated_mask_L_LD], hu_calcium, hu_heart_tissue_large_ld, voxel_size, ρ_calcium, VolumeFraction())

# ╔═╡ 0537e975-ad8c-475a-ac2e-dd38dc9c1f8f
md"""
# Score Medium Inserts
"""

# ╔═╡ 75e2d259-c54b-4f05-a2ee-cea87e26866a
md"""
## High Density
"""

# ╔═╡ 2f53e813-b084-4fa4-9db1-edda60737890
md"""
#### Dilated mask
"""

# ╔═╡ f17b1b5f-2cfa-4850-91e1-001dd55662df
@bind z2 overlay_mask_bind(dilated_mask_M_HD)

# ╔═╡ 2cd125c6-c2c0-43a1-8128-532af278b121
overlay_mask_plot(arr, dilated_mask_M_HD, z2, "dilated mask")

# ╔═╡ 4d06ad97-0ec0-4c70-8b02-d689e591411f
md"""
#### Ring (background) mask
"""

# ╔═╡ 9a535e6a-7e3c-413a-a02d-6f41e26da67f
@bind z4 overlay_mask_bind(ring_mask_M_HD)

# ╔═╡ c6e690d8-af7b-418b-aec6-67ce457eae54
overlay_mask_plot(arr, ring_mask_M_HD, z4, "ring mask")

# ╔═╡ 92c74076-20a3-4e81-a31b-57c7d3198a07
hu_heart_tissue_medium_hd = mean(arr[ring_mask_M_HD])

# ╔═╡ 0d4f00cb-a956-4bf3-9b24-d60b1c55df3a
mass_medium_hd = score(arr[dilated_mask_M_HD], hu_calcium, hu_heart_tissue_medium_hd, voxel_size, ρ_calcium, VolumeFraction())

# ╔═╡ 23e4ec31-263a-4be1-817b-94b84b89d3b4
md"""
## Medium Density
"""

# ╔═╡ 35dbd2cc-90bb-4cc8-a274-742a2b61fd47
md"""
#### Dilated mask
"""

# ╔═╡ aace1599-d1a0-43a7-8c13-f94281b1e6a7
@bind x2 overlay_mask_bind(dilated_mask_M_MD)

# ╔═╡ be540eea-76a9-4734-94e7-7ff34f3372fe
overlay_mask_plot(arr, dilated_mask_M_MD, x2, "dilated mask")

# ╔═╡ f9c1fdc6-c19e-47d5-a50a-2e47398b2c6d
md"""
#### Ring (background) mask
"""

# ╔═╡ 6a245113-0cbf-43d1-80b4-8cab748d6055
@bind x4 overlay_mask_bind(ring_mask_M_MD)

# ╔═╡ 9b6ac02f-be21-4158-9e48-bc01a3e688f9
overlay_mask_plot(arr, ring_mask_M_MD, x4, "ring mask")

# ╔═╡ 503ac757-dc03-417f-a5d7-6f4f7b90c2a7
hu_heart_tissue_medium_md = mean(arr[ring_mask_M_MD])

# ╔═╡ ce9820de-a672-43e7-bfc0-46d104c5d55a
mass_medium_md = score(arr[dilated_mask_M_MD], hu_calcium, hu_heart_tissue_medium_md, voxel_size, ρ_calcium, VolumeFraction())

# ╔═╡ 74d7fc3a-58e5-4c39-91ee-8239d11bd32c
md"""
## Low Density
"""

# ╔═╡ 52ea8424-87c5-45e2-a51f-4b48eac0a4b4
md"""
#### Dilated mask
"""

# ╔═╡ e7f43d90-4122-4bd1-a172-f8a541aff7dd
@bind y2 overlay_mask_bind(dilated_mask_M_LD)

# ╔═╡ 13f8c58e-3cc4-4229-9ff8-bca00263226b
overlay_mask_plot(arr, dilated_mask_M_LD, y2, "dilated mask")

# ╔═╡ bfdbe59b-f2ec-4b76-b7af-a77ff299e086
md"""
#### Ring (background) mask
"""

# ╔═╡ 68ef024a-aed7-4c9f-829a-475424e97ec5
@bind y4 overlay_mask_bind(ring_mask_M_LD)

# ╔═╡ 4d66583d-594b-4639-a073-430b49d4dd26
overlay_mask_plot(arr, ring_mask_M_LD, y4, "ring mask")

# ╔═╡ afd52af3-02db-4bde-99a1-836f1edbf1df
hu_heart_tissue_medium_ld = mean(arr[ring_mask_M_LD])

# ╔═╡ bcaed557-7114-4e44-acaf-7b768489d488
mass_medium_ld = score(arr[dilated_mask_M_LD], hu_calcium, hu_heart_tissue_medium_ld, voxel_size, ρ_calcium, VolumeFraction())

# ╔═╡ 45347bdc-59b2-4258-9413-80db932c13c7
md"""
# Score Small Inserts
"""

# ╔═╡ c8a3e0c0-a785-4446-bd61-042a659ecc47
md"""
## High Density
"""

# ╔═╡ 02ec8856-6c12-4a85-a78a-b9aa64efe8c4
md"""
#### Dilated mask
"""

# ╔═╡ 68cab80d-4fb9-47d9-aee1-b885205a1fe6
@bind l2 overlay_mask_bind(dilated_mask_S_HD)

# ╔═╡ efc90700-7f38-4048-a952-615046914668
overlay_mask_plot(arr, dilated_mask_S_HD, l2, "dilated mask")

# ╔═╡ 993064dd-7c27-42e8-a0e3-b71a037d5f50
md"""
#### Ring (background) mask
"""

# ╔═╡ 61a081c6-8650-4ba5-91af-9e604199dc3b
@bind l4 overlay_mask_bind(ring_mask_S_HD)

# ╔═╡ 0c00b98a-7439-41f2-b7e4-a13d9d4b9f9c
overlay_mask_plot(arr, ring_mask_S_HD, l4, "ring mask")

# ╔═╡ 98bf0863-29ce-436d-bc2f-279be0009d47
hu_heart_tissue_small_hd = mean(arr[ring_mask_S_HD])

# ╔═╡ 9a3fce26-7fe4-407b-943c-258f809aee1d
mass_small_hd = score(arr[dilated_mask_S_HD], hu_calcium, hu_heart_tissue_large_hd, voxel_size, ρ_calcium, VolumeFraction())

# ╔═╡ 6113ae63-0775-45c8-9c53-cc14fc99bbd6
md"""
## Medium Density
"""

# ╔═╡ f4566aae-9fc8-456a-a9bd-59f6e263f12f
md"""
#### Dilated mask
"""

# ╔═╡ 26f75883-dceb-4d00-a378-1f993e8505af
@bind k2 overlay_mask_bind(dilated_mask_S_MD)

# ╔═╡ a9a3fd12-1b83-437a-80f6-967b3dfbaa31
overlay_mask_plot(arr, dilated_mask_S_MD, k2, "dilated mask")

# ╔═╡ a1db64a8-e873-4c25-8eb5-aa59c67e26b2
md"""
#### Ring (background) mask
"""

# ╔═╡ e318c9fc-92aa-46fb-bf82-c723468a5761
@bind k4 overlay_mask_bind(ring_mask_S_MD)

# ╔═╡ 28b410d5-9143-4911-a92f-4f26fdae6900
overlay_mask_plot(arr, ring_mask_S_MD, k4, "ring mask")

# ╔═╡ deea4a88-0f90-4f76-8133-6a448a363d99
hu_heart_tissue_small_md = mean(arr[ring_mask_S_MD])

# ╔═╡ 0092662e-570a-4be2-a901-761b5da4bc22
mass_small_md = score(arr[dilated_mask_S_MD], hu_calcium, hu_heart_tissue_large_md, voxel_size, ρ_calcium, VolumeFraction())

# ╔═╡ db97b437-e959-4991-a0a0-14436bc07775
md"""
## Low Density
"""

# ╔═╡ e778b1e9-8a43-4052-926b-b5ae4a83904c
md"""
#### Dilated mask
"""

# ╔═╡ 4524347e-b729-4bd3-a4a6-3ed7bde51f9b
@bind m2 overlay_mask_bind(dilated_mask_S_LD)

# ╔═╡ 734b71d4-dbec-44fb-b3f7-1e92349d184f
overlay_mask_plot(arr, dilated_mask_S_LD, m2, "dilated mask")

# ╔═╡ 8cc5642c-a9cf-43c7-a4cd-1307a0c6ec3d
md"""
#### Ring (background) mask
"""

# ╔═╡ 731ded23-4289-4076-bf88-94f9b9500446
@bind m4 overlay_mask_bind(ring_mask_S_LD)

# ╔═╡ 381602f5-e8e0-49f4-9acd-b93edf179fce
overlay_mask_plot(arr, ring_mask_S_LD, m4, "ring mask")

# ╔═╡ 757e08b1-5ebe-409e-a003-7ba8ff31009a
hu_heart_tissue_small_ld = mean(arr[ring_mask_S_LD])

# ╔═╡ 9ee300e1-288e-4cb8-9b48-8e3b8afb9b49
mass_small_ld = score(arr[dilated_mask_S_LD], hu_calcium, hu_heart_tissue_large_ld, voxel_size, ρ_calcium, VolumeFraction())

# ╔═╡ 6f691ada-c9fc-4c12-8c05-f7c7c0f43f83
md"""
# Results
"""

# ╔═╡ f94bf7f8-6c95-4414-adef-2955f2218f93
PhantomSegmentation.get_pixel_size(header)

# ╔═╡ 3d2a470d-1b18-4efb-b2ef-f36b50c7c57b
inserts = ["Low Density", "Medium Density", "High Density"]

# ╔═╡ 62aa0caa-5e69-496c-9982-5d989c91afa4
volume_gt = [7.065, 63.585, 176.625]

# ╔═╡ 94c98db0-aa3d-4c0b-9f26-76e1ef08946a
ground_truth_mass_large = [
    volume_gt[3] * density_array[1],
    volume_gt[3] * density_array[2],
    volume_gt[3] * density_array[3],
] # mg

# ╔═╡ 31a4b17a-c819-4eaa-a8da-571a35a838ee
calculated_mass_large = [mass_large_ld, mass_large_md, mass_large_hd]

# ╔═╡ ce7a828f-fc28-4ae4-a400-99a593c2955f
ground_truth_mass_medium = [
    volume_gt[2] * density_array[1],
    volume_gt[2] * density_array[2],
    volume_gt[2] * density_array[3],
] # mg

# ╔═╡ 8c0c6c21-f8ca-43ce-a0e7-6644384d587b
calculated_mass_medium = [mass_medium_ld, mass_medium_md, mass_medium_hd]

# ╔═╡ b184e938-2e5c-40e5-8e52-ffdfff7eca5f
ground_truth_mass_small = [
    volume_gt[1] * density_array[1],
    volume_gt[1] * density_array[2],
    volume_gt[1] * density_array[3],
] # mg

# ╔═╡ 03c98f17-fae8-47ff-bddf-5acd08c0be1b
calculated_mass_small = [mass_small_ld, mass_small_md, mass_small_hd]

# ╔═╡ 31ac7d78-1780-4a2e-8c78-f3f7554300e5
df = DataFrame(;
    scan=scan,
    inserts=inserts,
    ground_truth_mass_large=ground_truth_mass_large,
    calculated_mass_large=calculated_mass_large,
    ground_truth_mass_medium=ground_truth_mass_medium,
    calculated_mass_medium=calculated_mass_medium,
    ground_truth_mass_small=ground_truth_mass_small,
    calculated_mass_small=calculated_mass_small,
)

# ╔═╡ 1e470c8e-52fa-49e1-83d2-4bcbbb7f7c15
md"""
### Save Results
"""

# ╔═╡ 8408801a-f00a-4458-8bbd-7ca56aaa8140
# if ~isdir(string(cd(pwd, "..") , "/output/", TYPE))
# 	mkdir(string(cd(pwd, "..") , "/output/", TYPE))
# end

# ╔═╡ 9814530b-c46e-42f2-8582-c4b2e556a179
# output_path = string(cd(pwd, "..") , "/output/", TYPE, "/", scan, ".csv")

# ╔═╡ ccdc9543-e5a3-4e2a-ac81-8f97bab08eec
# CSV.write(output_path, df)

# ╔═╡ 1b378b95-df3a-4f7c-b570-6e41743b9013
md"""
### Save full df
"""

# ╔═╡ f5bd0678-6e8e-4d42-9cf3-11e9ec8689a5
dfs = []

# ╔═╡ 7bd99563-f83f-4489-815b-62d136d1a1cd
push!(dfs, df)

# ╔═╡ 90c19e6e-5bc2-4360-86de-dd9e3bfd1e6c
# if length(dfs) == 24
# 	global new_df = vcat(dfs[1:24]...)
# 	output_path_new = string(cd(pwd, "..") , "/output/", TYPE, "/", "full.csv")
# 	CSV.write(output_path_new, new_df)
# end

# ╔═╡ Cell order:
# ╠═9860656e-409f-4dc7-bd95-e96781a104f0
# ╠═163ffa5f-0509-44db-a5e5-1f93d625079b
# ╟─01526785-a645-4b1a-afc2-0bc828c2db3a
# ╠═0a2c9749-e2bd-4d5e-9ef5-a07834f6684b
# ╟─28a061b7-5a67-4454-80d8-d8db807f2d33
# ╟─f07b68cb-97fe-4073-8363-55bbf4e49a3f
# ╟─5ae40e25-6c2b-4a69-ae75-71e022efc50f
# ╟─8ed6a20b-e2d9-4596-a296-8e0dccbd7e5b
# ╟─37b0192b-a590-421f-877a-339d96493d08
# ╟─0fdfc1be-f7ee-4d93-8186-2e0a3e121366
# ╟─4fea86c3-2c9c-4549-b4c6-b67ffb6bad04
# ╟─1e773785-94a6-490a-99d3-5da4478fa9a0
# ╟─8a68d139-8712-40ef-a2c5-ee33ba278174
# ╟─593b9d85-dcf7-4db4-9e78-ca284007d30e
# ╟─94200cb0-0c4e-4fa8-8ed9-418a9aa25469
# ╟─509af748-95e8-4122-9f81-aecc42c33901
# ╠═9a93a6b5-9ec6-4667-b9b0-ad0f6a58709f
# ╟─cd4c2e05-c305-4da4-a411-054bc64a810c
# ╠═09a6b6e2-9aee-4d77-b101-8e143a5d7fb3
# ╟─97d7830f-befe-4f4e-b5d7-a37c9cf690ff
# ╟─c52cf513-a741-4fe6-b8a2-117a1059c4fa
# ╟─6d0b5f0c-5a9b-470d-acbe-07bf5bafe27d
# ╟─60d95a36-f3d1-448f-a0cb-d92eaf9808a0
# ╟─63ed1025-518a-4ac0-b6be-c87fca2a53ca
# ╠═0db64444-a5c6-4c42-b129-e6f24f2ce426
# ╟─9cad1e1c-9e66-4e5a-86c6-9b02b0a0e30f
# ╟─e4dccc7e-d47c-41a5-9ffd-d535b1615f07
# ╟─6af669a9-77cd-480c-a7b6-1d9532ce62e1
# ╠═0fbf3a06-0850-444e-aa06-a3508661d0a4
# ╠═807b06fc-5246-4c8e-b7b1-389cdbda395c
# ╠═5aa9dd7b-82a5-40e0-a8ce-4c4079369c8f
# ╟─bdf43911-b1f3-4567-bb84-6781dbc320cc
# ╠═8a978945-b378-4ed4-8aab-3e6e5451a298
# ╠═0a78b66a-a545-409e-bc0b-5f575e6dea36
# ╟─2346e495-e1d4-4955-8ff0-cdff7fb07032
# ╠═6173c719-886d-4142-b2e3-cee8f9babdc6
# ╟─1cd4003d-02bb-4c6b-8bf0-0b5e1299d148
# ╠═072167d0-65e9-427f-991c-6e3aa23c172c
# ╠═97dbd4e9-ac4a-4f4e-aa2b-de123102bec1
# ╠═6a01dbc6-d8ca-4cad-a83c-995442c0b6e5
# ╟─e80a1e90-a89c-4530-88c0-dd2b80c7bd3f
# ╟─c095470d-c1a8-46eb-99c9-65cb3dfffd34
# ╟─ff0b1b99-ec2d-4610-97f2-ccb68a854460
# ╟─5aaabb67-d87b-4761-a927-2f53448ee697
# ╟─7d06e368-c46c-451a-84c8-424fbf12f06b
# ╟─24d9ad56-8c07-4af9-b216-de9bdcd16fb3
# ╟─208f81ed-06e6-4eb0-9c50-9c487547edb9
# ╠═6983213e-1c80-4738-9695-b83bab339207
# ╟─bf1474f6-9c56-4b0a-8309-42430dd5ef7d
# ╟─04df393b-7a3f-48be-944c-902146303088
# ╠═a1391510-889f-46f8-bf29-48a4d9266c3c
# ╠═6869fc6c-ca3d-4c44-9a47-9e964a5dbc02
# ╠═cecb93a2-d6fc-4b0d-acb7-c3e816d52498
# ╟─6e93a14b-7c0a-4555-bd81-bbd002aaef50
# ╟─4d675576-16e4-4241-8548-e40593e04893
# ╟─e9296b3e-3b28-4c8b-b44d-7ce65d6c3de7
# ╟─a144f068-6818-4a4d-8b61-87966d89c789
# ╟─30f68359-e661-4b74-8046-03cf8123a1a0
# ╟─94bfca4c-3068-4683-b8a2-239c7afc4d7d
# ╠═ab58e65d-940c-4058-8a02-56d23ae3bca5
# ╠═911ddf88-0149-4146-8e08-25995de7d0d4
# ╠═97c24681-62d5-410e-a2ff-8cad19d13147
# ╟─7a133076-83c2-4543-8c2e-0cbd9d234b27
# ╟─c20e99b7-2a62-4afb-ad8a-cd2f1815a7ca
# ╟─89181d61-763a-47fc-910e-7bd5c52a2db6
# ╠═f41e74e1-01a7-4581-bc1d-b9592005adcb
# ╟─2b6a668b-ad5b-4677-b1b5-3ae0e309d2cf
# ╟─40035644-21a6-4b1e-ace9-845af80fb9e6
# ╠═39be734a-5f4d-43a2-89c7-c49fcd2ab49d
# ╠═0b81a43b-6f34-4955-966b-3198ee522304
# ╠═9cd46d40-d35d-4fdc-8d3e-fc11b6072358
# ╟─0537e975-ad8c-475a-ac2e-dd38dc9c1f8f
# ╟─75e2d259-c54b-4f05-a2ee-cea87e26866a
# ╟─2f53e813-b084-4fa4-9db1-edda60737890
# ╟─f17b1b5f-2cfa-4850-91e1-001dd55662df
# ╠═2cd125c6-c2c0-43a1-8128-532af278b121
# ╟─4d06ad97-0ec0-4c70-8b02-d689e591411f
# ╟─9a535e6a-7e3c-413a-a02d-6f41e26da67f
# ╠═c6e690d8-af7b-418b-aec6-67ce457eae54
# ╠═92c74076-20a3-4e81-a31b-57c7d3198a07
# ╠═0d4f00cb-a956-4bf3-9b24-d60b1c55df3a
# ╟─23e4ec31-263a-4be1-817b-94b84b89d3b4
# ╟─35dbd2cc-90bb-4cc8-a274-742a2b61fd47
# ╟─aace1599-d1a0-43a7-8c13-f94281b1e6a7
# ╟─be540eea-76a9-4734-94e7-7ff34f3372fe
# ╟─f9c1fdc6-c19e-47d5-a50a-2e47398b2c6d
# ╟─6a245113-0cbf-43d1-80b4-8cab748d6055
# ╠═9b6ac02f-be21-4158-9e48-bc01a3e688f9
# ╠═503ac757-dc03-417f-a5d7-6f4f7b90c2a7
# ╠═ce9820de-a672-43e7-bfc0-46d104c5d55a
# ╟─74d7fc3a-58e5-4c39-91ee-8239d11bd32c
# ╟─52ea8424-87c5-45e2-a51f-4b48eac0a4b4
# ╟─e7f43d90-4122-4bd1-a172-f8a541aff7dd
# ╟─13f8c58e-3cc4-4229-9ff8-bca00263226b
# ╟─bfdbe59b-f2ec-4b76-b7af-a77ff299e086
# ╟─68ef024a-aed7-4c9f-829a-475424e97ec5
# ╠═4d66583d-594b-4639-a073-430b49d4dd26
# ╠═afd52af3-02db-4bde-99a1-836f1edbf1df
# ╠═bcaed557-7114-4e44-acaf-7b768489d488
# ╟─45347bdc-59b2-4258-9413-80db932c13c7
# ╟─c8a3e0c0-a785-4446-bd61-042a659ecc47
# ╟─02ec8856-6c12-4a85-a78a-b9aa64efe8c4
# ╟─68cab80d-4fb9-47d9-aee1-b885205a1fe6
# ╠═efc90700-7f38-4048-a952-615046914668
# ╟─993064dd-7c27-42e8-a0e3-b71a037d5f50
# ╟─61a081c6-8650-4ba5-91af-9e604199dc3b
# ╠═0c00b98a-7439-41f2-b7e4-a13d9d4b9f9c
# ╠═98bf0863-29ce-436d-bc2f-279be0009d47
# ╠═9a3fce26-7fe4-407b-943c-258f809aee1d
# ╟─6113ae63-0775-45c8-9c53-cc14fc99bbd6
# ╟─f4566aae-9fc8-456a-a9bd-59f6e263f12f
# ╟─26f75883-dceb-4d00-a378-1f993e8505af
# ╠═a9a3fd12-1b83-437a-80f6-967b3dfbaa31
# ╟─a1db64a8-e873-4c25-8eb5-aa59c67e26b2
# ╟─e318c9fc-92aa-46fb-bf82-c723468a5761
# ╠═28b410d5-9143-4911-a92f-4f26fdae6900
# ╠═deea4a88-0f90-4f76-8133-6a448a363d99
# ╠═0092662e-570a-4be2-a901-761b5da4bc22
# ╟─db97b437-e959-4991-a0a0-14436bc07775
# ╟─e778b1e9-8a43-4052-926b-b5ae4a83904c
# ╟─4524347e-b729-4bd3-a4a6-3ed7bde51f9b
# ╠═734b71d4-dbec-44fb-b3f7-1e92349d184f
# ╟─8cc5642c-a9cf-43c7-a4cd-1307a0c6ec3d
# ╟─731ded23-4289-4076-bf88-94f9b9500446
# ╠═381602f5-e8e0-49f4-9acd-b93edf179fce
# ╠═757e08b1-5ebe-409e-a003-7ba8ff31009a
# ╠═9ee300e1-288e-4cb8-9b48-8e3b8afb9b49
# ╟─6f691ada-c9fc-4c12-8c05-f7c7c0f43f83
# ╠═f94bf7f8-6c95-4414-adef-2955f2218f93
# ╠═3d2a470d-1b18-4efb-b2ef-f36b50c7c57b
# ╠═62aa0caa-5e69-496c-9982-5d989c91afa4
# ╠═94c98db0-aa3d-4c0b-9f26-76e1ef08946a
# ╠═31a4b17a-c819-4eaa-a8da-571a35a838ee
# ╠═ce7a828f-fc28-4ae4-a400-99a593c2955f
# ╠═8c0c6c21-f8ca-43ce-a0e7-6644384d587b
# ╠═b184e938-2e5c-40e5-8e52-ffdfff7eca5f
# ╠═03c98f17-fae8-47ff-bddf-5acd08c0be1b
# ╠═31ac7d78-1780-4a2e-8c78-f3f7554300e5
# ╟─1e470c8e-52fa-49e1-83d2-4bcbbb7f7c15
# ╠═8408801a-f00a-4458-8bbd-7ca56aaa8140
# ╠═9814530b-c46e-42f2-8582-c4b2e556a179
# ╠═ccdc9543-e5a3-4e2a-ac81-8f97bab08eec
# ╟─1b378b95-df3a-4f7c-b570-6e41743b9013
# ╠═f5bd0678-6e8e-4d42-9cf3-11e9ec8689a5
# ╠═7bd99563-f83f-4489-815b-62d136d1a1cd
# ╠═90c19e6e-5bc2-4360-86de-dd9e3bfd1e6c

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

# ╔═╡ 5e9fdac1-d4e6-4b87-ace9-002fb609b599
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.activate(".")

    using PlutoUI, Statistics, CSV, DataFrames, GLM, CairoMakie, HypothesisTests, Colors, MLJBase, DICOM, DICOMUtils, PhantomSegmentation, CalciumScoring, ImageMorphology, ImageFiltering, Noise
    using StatsBase: quantile!, rmsd
end

# ╔═╡ c08b5164-86e3-464e-9162-f9b40f86325a
TableOfContents()

# ╔═╡ 8eb50b0a-4502-473d-91cf-4c143b5c8774
md"""
## Load DICOMS

All you need to do is set `base_path` once and leave it. After that, the only thing that should change is the `VENDOR`, once for every set, and the `SCAN_NUMBER`, once for each scan.
"""

# ╔═╡ 4afeee33-e030-430f-899a-d1de13cade18
begin
	VENDORS = ["80", "100", "120", "135"]
	SIZES = ["small", "medium", "large"]
	DENSITIES = ["low", "normal"]

	IMAGES = "images_new"

	VENDOR = VENDORS[4]
	SIZE = SIZES[1]
	DENSITY = DENSITIES[2]
	
    BASE_PATH = joinpath("/Users/daleblack/Google Drive/dev/MolloiLab/cac-simulation", IMAGES, SIZE, DENSITY)
	root_path = joinpath(BASE_PATH, VENDOR * "-motion")
	dcm_path_list = dcm_list_builder(root_path)
	pth = dcm_path_list[1]
	scan = basename(pth)
	header, dcm_array, slice_thick_ori1 = dcm_reader(pth)
	kV = parse(Int64, VENDOR)
end

# ╔═╡ ab78ec80-d30d-4e58-8371-a82fc9c2f1d7
md"""
## Helper Functions
"""

# ╔═╡ aa8d4145-e221-4755-b7da-1643668f2309
function collect_tuple(tuple_array)
    row_num = size(tuple_array)
    col_num = length(tuple_array[1])
    container = zeros(Int64, row_num..., col_num)
    for i in 1:length(tuple_array)
        container[i, :] = collect(tuple_array[i])
    end
    return container
end

# ╔═╡ 569f7b08-2b2d-430e-9336-40301d75b8c8
function overlay_mask_bind(mask)
    indices = findall(x -> x == 1, mask)
    indices = Tuple.(indices)
    label_array = collect_tuple(indices)
    zs = unique(label_array[:, 3])
    return PlutoUI.Slider(1:length(zs); default=3, show_value=true)
end

# ╔═╡ 841b3032-9ae7-4584-bea9-c3e54e800d61
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

# ╔═╡ 8fe3bc3f-09a6-41ce-999a-e27a3f6092c8
function create_mask(array, mask)
    @assert size(array) == size(mask)
    idxs = findall(x -> x == true, mask)
    overlayed_mask = zeros(size(array))
    for idx in idxs
        overlayed_mask[idx] = array[idx]
    end
    return overlayed_mask
end

# ╔═╡ ef121fd3-4390-4fa1-8e06-23837d9dacba
function dilate_mask_large(mask)
    return dilate(mask)
end

# ╔═╡ 5bd36e23-7f1c-4bcf-aae7-5fd5f3d23776
function ring_mask_large(dilated_mask)
    return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask)
end

# ╔═╡ 299c48b4-fb04-4147-850c-b5781e9f75da
function dilate_mask_medium(mask)
    return (mask)
end

# ╔═╡ 11bc564a-0b41-42c3-9fd3-27ffc411b87c
function ring_mask_medium(dilated_mask)
    return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask)
end

# ╔═╡ 75555168-c295-40ed-9b28-1548cb18a062
function dilate_mask_small(mask)
    return (mask)
end

# ╔═╡ 31828197-d025-49c1-8949-e5ba6531f449
function ring_mask_small(dilated_mask)
    return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask)
end

# ╔═╡ c06447b6-54b2-403d-874e-35b6b87b2b12
md"""
## Segment Heart
"""

# ╔═╡ b65c8925-d60b-47a8-aa55-b4fa96cf3c9f
masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3) ÷ 2);

# ╔═╡ 16c0c5a4-5c50-43d6-89e7-f7b804ac08ef
@bind a PlutoUI.Slider(1:size(masked_array, 3), default=10, show_value=true)

# ╔═╡ 8910bea9-b4dc-4f99-875a-4de329caac12
heatmap(transpose(masked_array[:, :, a]); colormap=:grays)

# ╔═╡ 05a79599-4aa1-43fa-9339-004617dbba9c
let
    fig = Figure()

    ax = Makie.Axis(fig[1, 1])
    heatmap!(transpose(dcm_array[:, :, 5]); colormap=:grays)

	save("/Users/daleblack/Google Drive/Research/Papers/My Papers/cac-simulation/figures-review/simulated_phantom.png", fig)
    fig
end

# ╔═╡ 4b536ab0-c98a-4f77-8cbc-57188c3feef1
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

# ╔═╡ 7af38d16-3186-4818-9158-934afd9da116
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

# ╔═╡ 729ba9a0-bb10-4692-a552-78f6d147ce21
let
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

# ╔═╡ a727b372-e880-401d-ab3c-20c7412707c4
md"""
## Segment Calcium Rod
"""

# ╔═╡ bddf5fb7-fa61-4529-90d8-5aaedac3d79e
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

# ╔═╡ f2b29838-ec4a-4f28-9f97-36510a73a600
@bind c PlutoUI.Slider(1:size(calcium_image, 3), default=5, show_value=true)

# ╔═╡ d9465f89-36b8-447b-939e-a5495f1e14d3
heatmap(transpose(calcium_image[:, :, c]); colormap=:grays)

# ╔═╡ 715aa9ec-8c8b-44e0-850f-19dccd4b79f7
md"""
## Load Masks
"""

# ╔═╡ 47bab2d4-d8b6-44ca-8b04-615906fb1e96
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

# ╔═╡ 56928438-f9d9-47b5-bfef-5549aa183614
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

# ╔═╡ 920c5619-9943-484a-99b6-b11e9e5a4a89
heatmap(transpose(masks); colormap=:grays)

# ╔═╡ 341c5fbb-321d-4702-a443-c6a11a11f22b
md"""
## Prepare Array Masks
"""

# ╔═╡ 908780f5-8280-48ce-bb0d-98d5e0ee3c6b
begin
	arr = masked_array[:, :, 4:6]
	pixel_size = DICOMUtils.get_pixel_size(header)
	voxel_size = pixel_size[1] * pixel_size[2] * pixel_size[3]
end;

# ╔═╡ 8a1b53d3-ec47-47ec-ad23-30bf6ed4ec0c
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

# ╔═╡ 852eda80-17c6-42c6-a015-032afaadc44e
md"""
## Calibration Prep
"""

# ╔═╡ 0b55ef6a-913d-4c9d-8d5f-1bfa7e555b61
begin
	if DENSITY == "low"
		density_array = [0.025, 0.050, 0.100]
	elseif DENSITY == "normal"
		density_array = [0.200, 0.400, 0.800]
	end
	array_filtered = abs.(mapwindow(median, calcium_image[:, :, 2], (3, 3)))
	bool_arr = array_filtered .> 0
	bool_arr_erode = (erode(erode(erode(erode(bool_arr)))))
	c_img = calcium_image[:, :, 1:3]
	mask_cal_3D = Array{Bool}(undef, size(c_img))
	for z in 1:size(c_img, 3)
		mask_cal_3D[:, :, z] = bool_arr_erode
	end

	# hu_calcium = mean(c_img[mask_cal_3D])
	if VENDOR == "80"
		hu_calcium = 377.3
	elseif VENDOR == "100"
		hu_calcium = 326.0
	elseif VENDOR == "120"
		hu_calcium = 296.9
	else
		hu_calcium = 282.7
	end
	ρ_calcium = 0.2
end

# ╔═╡ 295ff81c-cc6d-44fb-b10c-86f24fb3bf9e
hu_calcium

# ╔═╡ 15a30f21-f31f-4c56-b955-b8f119834c32
md"""
# Background
"""

# ╔═╡ 21e1ae40-dfd0-44c5-8c1d-ef5a68d13f79
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

# ╔═╡ 0a062ed7-661e-478b-80a3-d0d54de035fb
hu_heart_tissue_bkg

# ╔═╡ f9c9f86a-1e90-448a-b0e2-270e82068f27
size(background_mask)

# ╔═╡ 6dc6ba21-30c9-48ec-aec1-9ea817bcae9c
@bind q1 overlay_mask_bind(background_mask)

# ╔═╡ 0853e72a-0e38-4ff8-adad-40f740dd1ddd
overlay_mask_plot(arr, background_mask, q1, "dilated mask")

# ╔═╡ 6b475b90-7f3e-4f44-8929-7056b7182046
overlay_mask_plot(arr, background_ring, q1, "dilated mask")

# ╔═╡ 1dbad5ed-2493-4e39-8ec9-1e5e3720ddf7
md"""
# Score Large Inserts
"""

# ╔═╡ 2b9cc53e-bcd3-407b-bf90-c65947b28581
md"""
## High Density
"""

# ╔═╡ ddb290bf-b754-4b41-8430-963bf835ae48
md"""
#### Dilated mask
"""

# ╔═╡ 7e95a1de-6021-48ea-bd77-3480abf7b92d
@bind g2 overlay_mask_bind(dilated_mask_L_HD)

# ╔═╡ 00998177-1039-40ce-ada9-1dcdf611ac98
overlay_mask_plot(arr, dilated_mask_L_HD, g2, "dilated mask")

# ╔═╡ 3dc58caa-68d3-4af5-96b4-0d02c7116a00
md"""
#### Ring (background) mask
"""

# ╔═╡ 4833bb01-1b20-41a5-a5b9-023b67f41c54
@bind g4 overlay_mask_bind(ring_mask_L_HD)

# ╔═╡ 62b6a7ef-59e9-4ec3-adc6-f9992d9ebfce
overlay_mask_plot(arr, ring_mask_L_HD, g4, "ring mask")

# ╔═╡ 54de1790-5194-490c-b0fd-e18692dd933a
hu_heart_tissue_large_hd = mean(arr[ring_mask_L_HD])

# ╔═╡ bf788e3f-538e-425e-a06a-6e4499259fee
mass_large_hd = score(arr[dilated_mask_L_HD], hu_calcium, hu_heart_tissue_bkg, voxel_size, ρ_calcium, VolumeFraction())

# ╔═╡ 1b0d7eeb-d6dc-4c38-b5bb-624b54235b4b
md"""
## Medium Density
"""

# ╔═╡ c1fea4b2-7740-4d62-b908-d89261d25064
md"""
#### Dilated mask
"""

# ╔═╡ 45900f95-7501-4655-aa3b-2c4979ae140a
@bind h2 overlay_mask_bind(dilated_mask_L_MD)

# ╔═╡ 43a80610-f2ad-44a2-a2e2-33717211b29e
overlay_mask_plot(arr, dilated_mask_L_MD, h2, "dilated mask")

# ╔═╡ 57275d3f-1707-408e-96d3-bc1b31e61e96
md"""
#### Ring (background) mask
"""

# ╔═╡ d667363c-ee29-47e6-8737-23a7ec0b3de3
@bind h4 overlay_mask_bind(ring_mask_L_MD)

# ╔═╡ 319c8a0e-4365-41b7-8a4a-95a436f62ff7
overlay_mask_plot(arr, ring_mask_L_MD, h4, "ring mask")

# ╔═╡ 2269784e-540c-4dcd-94f4-8a53285b65f0
hu_heart_tissue_large_md = mean(arr[ring_mask_L_MD])

# ╔═╡ 42a5aa69-8d08-4eba-a952-00650bad6905
mass_large_md = score(arr[dilated_mask_L_MD], hu_calcium, hu_heart_tissue_large_md, voxel_size, ρ_calcium, VolumeFraction())

# ╔═╡ 39215a15-ca20-4eb4-ac7c-07bb022e410b
md"""
## Low Density
"""

# ╔═╡ bb4670bf-0c80-4039-b610-2db7103efd1e
md"""
#### Dilated mask
"""

# ╔═╡ 50214c11-1a09-4137-b867-525327ce45d5
@bind j2 overlay_mask_bind(dilated_mask_L_LD)

# ╔═╡ 1f1a5df3-8e8c-49f1-b9e4-ab966340725f
overlay_mask_plot(arr, dilated_mask_L_LD, j2, "dilated mask")

# ╔═╡ 49ef1452-3219-4f4d-aa45-42bd291ecc31
md"""
#### Ring (background) mask
"""

# ╔═╡ 29b88c85-4255-4080-995c-30bd7be7ddec
@bind j4 overlay_mask_bind(ring_mask_L_LD)

# ╔═╡ 2e7566eb-e208-472d-9cbf-4b528fe74a25
overlay_mask_plot(arr, ring_mask_L_LD, j4, "ring mask")

# ╔═╡ c1ff9bcc-4485-4e41-998a-29e5772add07
hu_heart_tissue_large_ld = mean(arr[ring_mask_L_LD])

# ╔═╡ c70bc1a1-7ffe-43bb-8955-874d6fc3550d
mass_large_ld = score(arr[dilated_mask_L_LD], hu_calcium, hu_heart_tissue_large_ld, voxel_size, ρ_calcium, VolumeFraction())

# ╔═╡ a2577753-a937-47f8-b14f-8c2309940bae
md"""
# Score Medium Inserts
"""

# ╔═╡ e9f53bce-0d59-456f-8add-ec3e7e0aa47e
md"""
## High Density
"""

# ╔═╡ 4a65c23f-942c-45e0-ae58-0712240e9668
md"""
#### Dilated mask
"""

# ╔═╡ 3a970363-9789-440b-bc47-3d4dbf9d3fd3
@bind z2 overlay_mask_bind(dilated_mask_M_HD)

# ╔═╡ 6b560d89-2ee0-4b26-88c3-a2554f01bc83
overlay_mask_plot(arr, dilated_mask_M_HD, z2, "dilated mask")

# ╔═╡ 6cc2e7d6-ca60-4f44-b16b-08c8d6a80bb8
md"""
#### Ring (background) mask
"""

# ╔═╡ 0d532334-ddc4-4cf3-9453-d1a0e65e7a60
@bind z4 overlay_mask_bind(ring_mask_M_HD)

# ╔═╡ a77d1304-b887-4d77-8b3a-350527549c23
overlay_mask_plot(arr, ring_mask_M_HD, z4, "ring mask")

# ╔═╡ d0fd3ef5-8e1c-4c75-bcef-94bfcb52f69b
hu_heart_tissue_medium_hd = mean(arr[ring_mask_M_HD])

# ╔═╡ 93c6c475-cbae-442c-aceb-0b7127b268dd
mass_medium_hd = score(arr[dilated_mask_M_HD], hu_calcium, hu_heart_tissue_medium_hd, voxel_size, ρ_calcium, VolumeFraction())

# ╔═╡ 945a529f-c24e-4122-886f-1402672ed625
md"""
## Medium Density
"""

# ╔═╡ 040a62b8-7c41-4a39-8ecd-0a6242d3ff5b
md"""
#### Dilated mask
"""

# ╔═╡ e57598a1-35c6-42fa-82d4-2c976fccb629
@bind x2 overlay_mask_bind(dilated_mask_M_MD)

# ╔═╡ f2652dbd-28a3-406a-9bb3-8ba226780aa0
overlay_mask_plot(arr, dilated_mask_M_MD, x2, "dilated mask")

# ╔═╡ b8be5890-a0fc-49ca-85e3-915ff81979ec
md"""
#### Ring (background) mask
"""

# ╔═╡ a20af150-7341-47f0-9d1a-d5db8e8abc09
@bind x4 overlay_mask_bind(ring_mask_M_MD)

# ╔═╡ 81f4f874-f309-4764-b492-b2912108fcce
overlay_mask_plot(arr, ring_mask_M_MD, x4, "ring mask")

# ╔═╡ 352ed759-eb3b-4abb-8886-2f20d53bbfb9
hu_heart_tissue_medium_md = mean(arr[ring_mask_M_MD])

# ╔═╡ 2a152532-109a-461d-8eec-3a37ed236161
mass_medium_md = score(arr[dilated_mask_M_MD], hu_calcium, hu_heart_tissue_medium_md, voxel_size, ρ_calcium, VolumeFraction())

# ╔═╡ 3632354d-4850-45d6-a9d1-15949f8c0eb4
md"""
## Low Density
"""

# ╔═╡ 62b7d043-47fe-4aa2-8936-2739113be27d
md"""
#### Dilated mask
"""

# ╔═╡ 74722206-da16-4772-8ad3-b351bf0dd174
@bind y2 overlay_mask_bind(dilated_mask_M_LD)

# ╔═╡ 49d8dc4d-acac-42be-b8c5-44160d42cec6
overlay_mask_plot(arr, dilated_mask_M_LD, y2, "dilated mask")

# ╔═╡ 2e8ad4d6-4d1d-41c4-9e22-9049a03d374f
md"""
#### Ring (background) mask
"""

# ╔═╡ 4391a50d-0642-4804-bf0a-0695ff103b28
@bind y4 overlay_mask_bind(ring_mask_M_LD)

# ╔═╡ dad32de5-4c26-4d93-996a-9a2111d39ada
overlay_mask_plot(arr, ring_mask_M_LD, y4, "ring mask")

# ╔═╡ e70b9be1-711b-41b0-9888-c31b3a2dcba4
hu_heart_tissue_medium_ld = mean(arr[ring_mask_M_LD])

# ╔═╡ ec9dfcd9-13eb-41a7-be8c-405483d84a92
mass_medium_ld = score(arr[dilated_mask_M_LD], hu_calcium, hu_heart_tissue_medium_ld, voxel_size, ρ_calcium, VolumeFraction())

# ╔═╡ 63fc39c5-16c1-4b42-99de-feb449776beb
md"""
# Score Small Inserts
"""

# ╔═╡ 7cccce6c-954f-4720-9dd8-31dcf442208d
md"""
## High Density
"""

# ╔═╡ b524f8b8-64ed-4e9d-b1f5-2e0bd49bdaa4
md"""
#### Dilated mask
"""

# ╔═╡ a44c104b-577e-4746-ad94-53e1405a442f
@bind l2 overlay_mask_bind(dilated_mask_S_HD)

# ╔═╡ 8aa91548-24bb-460f-a1eb-072938b96641
overlay_mask_plot(arr, dilated_mask_S_HD, l2, "dilated mask")

# ╔═╡ 09f0a38f-64c5-4795-96ba-34e96883f595
md"""
#### Ring (background) mask
"""

# ╔═╡ 73d06ecc-baa2-4b5e-9411-a2f8c253efaa
@bind l4 overlay_mask_bind(ring_mask_S_HD)

# ╔═╡ be7ece4d-bf25-4c9d-baaa-1aebb696c4e2
overlay_mask_plot(arr, ring_mask_S_HD, l4, "ring mask")

# ╔═╡ 1aabeba7-5619-4774-96dc-eaab08b5b14a
hu_heart_tissue_small_hd = mean(arr[ring_mask_S_HD])

# ╔═╡ 2df030b1-87fd-4eaa-923d-56035a9a825f
mass_small_hd = score(arr[dilated_mask_S_HD], hu_calcium, hu_heart_tissue_large_hd, voxel_size, ρ_calcium, VolumeFraction())

# ╔═╡ 6bdf5633-319f-4e94-b69e-469cccdd017a
md"""
## Medium Density
"""

# ╔═╡ c7c1b57e-9f03-4cb6-b141-d4bf15d08905
md"""
#### Dilated mask
"""

# ╔═╡ c4d7421e-3046-431f-aabc-cab0f88de453
@bind k2 overlay_mask_bind(dilated_mask_S_MD)

# ╔═╡ 5ebb2e49-79e4-462d-92a6-6c3224c64e36
overlay_mask_plot(arr, dilated_mask_S_MD, k2, "dilated mask")

# ╔═╡ e2f2a55f-9934-45fd-9541-4d30de449cfe
md"""
#### Ring (background) mask
"""

# ╔═╡ 2fa9406a-0333-45a6-a674-5a9942aa85d7
@bind k4 overlay_mask_bind(ring_mask_S_MD)

# ╔═╡ 60293d04-c27a-4c17-8b96-59f9d42ba781
overlay_mask_plot(arr, ring_mask_S_MD, k4, "ring mask")

# ╔═╡ 7622bbd6-b851-4f59-836d-0f53ca8c610f
hu_heart_tissue_small_md = mean(arr[ring_mask_S_MD])

# ╔═╡ f5924025-4673-4218-a222-630584a52dc6
mass_small_md = score(arr[dilated_mask_S_MD], hu_calcium, hu_heart_tissue_large_md, voxel_size, ρ_calcium, VolumeFraction())

# ╔═╡ 3afec334-cc16-4f3d-a068-926a31e21ee5
md"""
## Low Density
"""

# ╔═╡ d50286ee-6dca-468e-9000-9a32620703ec
md"""
#### Dilated mask
"""

# ╔═╡ 8c7021cd-d89f-422d-9377-a5fdba5759ab
@bind m2 overlay_mask_bind(dilated_mask_S_LD)

# ╔═╡ 122be89d-bffa-418b-8449-b273bf73ec67
overlay_mask_plot(arr, dilated_mask_S_LD, m2, "dilated mask")

# ╔═╡ 486207ce-fd95-410f-9286-75d7b2f21d32
md"""
#### Ring (background) mask
"""

# ╔═╡ 45be9e86-390c-4ab7-8d01-c2a2ebb1e245
@bind m4 overlay_mask_bind(ring_mask_S_LD)

# ╔═╡ 355fd701-c652-4421-a906-1eb1f7506cc4
overlay_mask_plot(arr, ring_mask_S_LD, m4, "ring mask")

# ╔═╡ 6e8e836b-5069-4de7-abee-34f7bf882a13
hu_heart_tissue_small_ld = mean(arr[ring_mask_S_LD])

# ╔═╡ 9e3b40f2-9a8a-43e3-854a-b1b8b3e7be00
mass_small_ld = score(arr[dilated_mask_S_LD], hu_calcium, hu_heart_tissue_large_ld, voxel_size, ρ_calcium, VolumeFraction())

# ╔═╡ 269dde91-76cf-499b-8a51-a449c7ee7c40
md"""
# Results
"""

# ╔═╡ d631b3c8-65f5-439a-acb6-b15efe8e4425
PhantomSegmentation.get_pixel_size(header)

# ╔═╡ fd8677c2-0c45-42fc-b7c3-c8c618ce4a1e
inserts = ["Low Density", "Medium Density", "High Density"]

# ╔═╡ e3670840-3582-4256-8c3d-a253712c9101
volume_gt = [7.065, 63.585, 176.625]

# ╔═╡ ebcb1fbe-6c28-4c58-a856-3c6ac9466b23
ground_truth_mass_large = [
    volume_gt[3] * density_array[1],
    volume_gt[3] * density_array[2],
    volume_gt[3] * density_array[3],
] # mg

# ╔═╡ 5a8c91e9-1d1b-4bc7-b1df-96729792a920
calculated_mass_large = [mass_large_ld, mass_large_md, mass_large_hd]

# ╔═╡ 786e8ba3-1f22-459e-9273-f623262f136a
ground_truth_mass_medium = [
    volume_gt[2] * density_array[1],
    volume_gt[2] * density_array[2],
    volume_gt[2] * density_array[3],
] # mg

# ╔═╡ 6f8faab9-d447-4e6d-a8cd-19bfeb023dec
calculated_mass_medium = [mass_medium_ld, mass_medium_md, mass_medium_hd]

# ╔═╡ 6cf9cb98-7ee1-4762-81e3-9996f25ceff7
ground_truth_mass_small = [
    volume_gt[1] * density_array[1],
    volume_gt[1] * density_array[2],
    volume_gt[1] * density_array[3],
] # mg

# ╔═╡ 8f909b69-739c-40a2-b922-0589a5aac8cf
calculated_mass_small = [mass_small_ld, mass_small_md, mass_small_hd]

# ╔═╡ 9086dcbc-8033-4e79-a311-eb14438dac4c
df = DataFrame(;
    DENSITY=DENSITY,
    scan=scan,
    inserts=inserts,
    ground_truth_mass_large=ground_truth_mass_large,
    calculated_mass_large=calculated_mass_large,
    ground_truth_mass_medium=ground_truth_mass_medium,
    calculated_mass_medium=calculated_mass_medium,
    ground_truth_mass_small=ground_truth_mass_small,
    calculated_mass_small=calculated_mass_small,
)

# ╔═╡ 9754b9de-0a56-42eb-a1df-cd5a46d4f40c
df

# ╔═╡ 3ff725ed-be27-496a-ba11-3734d47e32cd
df

# ╔═╡ 2b51aaa8-7363-435d-b1f8-5252e85c21f3


# ╔═╡ 13191733-be50-4428-a5d1-57145cea221a


# ╔═╡ 2991a411-c270-4425-9cf6-1a3b4f759668


# ╔═╡ 082cf4fb-382f-4d5c-9ade-08018a66ac87
md"""
### Save Results
"""

# ╔═╡ 8d31c67f-db73-41b3-9393-16c86c4b05d8
# if ~isdir(string(cd(pwd, "..") , "/output/", TYPE))
# 	mkdir(string(cd(pwd, "..") , "/output/", TYPE))
# end

# ╔═╡ 890a10df-5beb-41bb-8c90-62de5148bd4c
# output_path = string(cd(pwd, "..") , "/output/", TYPE, "/", scan, ".csv")

# ╔═╡ deaad3fd-2eeb-42ab-b1ff-96eec9c6c2fd
# CSV.write(output_path, df)

# ╔═╡ d6f38052-f243-4d48-8ecb-492758c6cd28
md"""
### Save full df
"""

# ╔═╡ 23effaf5-b7a3-495b-aad5-209c8c640306
dfs = []

# ╔═╡ f4615674-173f-4245-9e09-1ba49ff2dedb
push!(dfs, df)

# ╔═╡ fa324943-c6d9-41c5-9122-d268783bed15
# if length(dfs) == 24
# 	global new_df = vcat(dfs[1:24]...)
# 	output_path_new = string(cd(pwd, "..") , "/output/", TYPE, "/", "full.csv")
# 	CSV.write(output_path_new, new_df)
# end

# ╔═╡ Cell order:
# ╠═5e9fdac1-d4e6-4b87-ace9-002fb609b599
# ╠═c08b5164-86e3-464e-9162-f9b40f86325a
# ╟─8eb50b0a-4502-473d-91cf-4c143b5c8774
# ╠═4afeee33-e030-430f-899a-d1de13cade18
# ╟─ab78ec80-d30d-4e58-8371-a82fc9c2f1d7
# ╟─aa8d4145-e221-4755-b7da-1643668f2309
# ╟─569f7b08-2b2d-430e-9336-40301d75b8c8
# ╟─841b3032-9ae7-4584-bea9-c3e54e800d61
# ╟─8fe3bc3f-09a6-41ce-999a-e27a3f6092c8
# ╟─ef121fd3-4390-4fa1-8e06-23837d9dacba
# ╟─5bd36e23-7f1c-4bcf-aae7-5fd5f3d23776
# ╟─299c48b4-fb04-4147-850c-b5781e9f75da
# ╟─11bc564a-0b41-42c3-9fd3-27ffc411b87c
# ╟─75555168-c295-40ed-9b28-1548cb18a062
# ╟─31828197-d025-49c1-8949-e5ba6531f449
# ╟─c06447b6-54b2-403d-874e-35b6b87b2b12
# ╠═b65c8925-d60b-47a8-aa55-b4fa96cf3c9f
# ╟─16c0c5a4-5c50-43d6-89e7-f7b804ac08ef
# ╟─8910bea9-b4dc-4f99-875a-4de329caac12
# ╟─05a79599-4aa1-43fa-9339-004617dbba9c
# ╟─4b536ab0-c98a-4f77-8cbc-57188c3feef1
# ╟─7af38d16-3186-4818-9158-934afd9da116
# ╟─729ba9a0-bb10-4692-a552-78f6d147ce21
# ╟─a727b372-e880-401d-ab3c-20c7412707c4
# ╠═bddf5fb7-fa61-4529-90d8-5aaedac3d79e
# ╟─f2b29838-ec4a-4f28-9f97-36510a73a600
# ╟─d9465f89-36b8-447b-939e-a5495f1e14d3
# ╟─715aa9ec-8c8b-44e0-850f-19dccd4b79f7
# ╠═47bab2d4-d8b6-44ca-8b04-615906fb1e96
# ╠═56928438-f9d9-47b5-bfef-5549aa183614
# ╠═920c5619-9943-484a-99b6-b11e9e5a4a89
# ╟─341c5fbb-321d-4702-a443-c6a11a11f22b
# ╠═908780f5-8280-48ce-bb0d-98d5e0ee3c6b
# ╠═8a1b53d3-ec47-47ec-ad23-30bf6ed4ec0c
# ╟─852eda80-17c6-42c6-a015-032afaadc44e
# ╠═0b55ef6a-913d-4c9d-8d5f-1bfa7e555b61
# ╠═295ff81c-cc6d-44fb-b10c-86f24fb3bf9e
# ╟─15a30f21-f31f-4c56-b955-b8f119834c32
# ╠═21e1ae40-dfd0-44c5-8c1d-ef5a68d13f79
# ╠═0a062ed7-661e-478b-80a3-d0d54de035fb
# ╠═f9c9f86a-1e90-448a-b0e2-270e82068f27
# ╟─6dc6ba21-30c9-48ec-aec1-9ea817bcae9c
# ╟─0853e72a-0e38-4ff8-adad-40f740dd1ddd
# ╟─6b475b90-7f3e-4f44-8929-7056b7182046
# ╟─1dbad5ed-2493-4e39-8ec9-1e5e3720ddf7
# ╟─2b9cc53e-bcd3-407b-bf90-c65947b28581
# ╟─ddb290bf-b754-4b41-8430-963bf835ae48
# ╟─7e95a1de-6021-48ea-bd77-3480abf7b92d
# ╠═00998177-1039-40ce-ada9-1dcdf611ac98
# ╟─3dc58caa-68d3-4af5-96b4-0d02c7116a00
# ╟─4833bb01-1b20-41a5-a5b9-023b67f41c54
# ╠═62b6a7ef-59e9-4ec3-adc6-f9992d9ebfce
# ╠═54de1790-5194-490c-b0fd-e18692dd933a
# ╠═bf788e3f-538e-425e-a06a-6e4499259fee
# ╟─1b0d7eeb-d6dc-4c38-b5bb-624b54235b4b
# ╟─c1fea4b2-7740-4d62-b908-d89261d25064
# ╟─45900f95-7501-4655-aa3b-2c4979ae140a
# ╠═43a80610-f2ad-44a2-a2e2-33717211b29e
# ╟─57275d3f-1707-408e-96d3-bc1b31e61e96
# ╟─d667363c-ee29-47e6-8737-23a7ec0b3de3
# ╠═319c8a0e-4365-41b7-8a4a-95a436f62ff7
# ╠═2269784e-540c-4dcd-94f4-8a53285b65f0
# ╠═42a5aa69-8d08-4eba-a952-00650bad6905
# ╟─39215a15-ca20-4eb4-ac7c-07bb022e410b
# ╟─bb4670bf-0c80-4039-b610-2db7103efd1e
# ╟─50214c11-1a09-4137-b867-525327ce45d5
# ╠═1f1a5df3-8e8c-49f1-b9e4-ab966340725f
# ╟─49ef1452-3219-4f4d-aa45-42bd291ecc31
# ╟─29b88c85-4255-4080-995c-30bd7be7ddec
# ╠═2e7566eb-e208-472d-9cbf-4b528fe74a25
# ╠═c1ff9bcc-4485-4e41-998a-29e5772add07
# ╠═c70bc1a1-7ffe-43bb-8955-874d6fc3550d
# ╟─a2577753-a937-47f8-b14f-8c2309940bae
# ╠═9754b9de-0a56-42eb-a1df-cd5a46d4f40c
# ╟─e9f53bce-0d59-456f-8add-ec3e7e0aa47e
# ╟─4a65c23f-942c-45e0-ae58-0712240e9668
# ╟─3a970363-9789-440b-bc47-3d4dbf9d3fd3
# ╠═6b560d89-2ee0-4b26-88c3-a2554f01bc83
# ╟─6cc2e7d6-ca60-4f44-b16b-08c8d6a80bb8
# ╟─0d532334-ddc4-4cf3-9453-d1a0e65e7a60
# ╠═a77d1304-b887-4d77-8b3a-350527549c23
# ╠═d0fd3ef5-8e1c-4c75-bcef-94bfcb52f69b
# ╠═93c6c475-cbae-442c-aceb-0b7127b268dd
# ╟─945a529f-c24e-4122-886f-1402672ed625
# ╟─040a62b8-7c41-4a39-8ecd-0a6242d3ff5b
# ╟─e57598a1-35c6-42fa-82d4-2c976fccb629
# ╠═f2652dbd-28a3-406a-9bb3-8ba226780aa0
# ╟─b8be5890-a0fc-49ca-85e3-915ff81979ec
# ╟─a20af150-7341-47f0-9d1a-d5db8e8abc09
# ╠═81f4f874-f309-4764-b492-b2912108fcce
# ╠═352ed759-eb3b-4abb-8886-2f20d53bbfb9
# ╠═2a152532-109a-461d-8eec-3a37ed236161
# ╟─3632354d-4850-45d6-a9d1-15949f8c0eb4
# ╟─62b7d043-47fe-4aa2-8936-2739113be27d
# ╟─74722206-da16-4772-8ad3-b351bf0dd174
# ╠═49d8dc4d-acac-42be-b8c5-44160d42cec6
# ╟─2e8ad4d6-4d1d-41c4-9e22-9049a03d374f
# ╟─4391a50d-0642-4804-bf0a-0695ff103b28
# ╠═dad32de5-4c26-4d93-996a-9a2111d39ada
# ╠═e70b9be1-711b-41b0-9888-c31b3a2dcba4
# ╠═ec9dfcd9-13eb-41a7-be8c-405483d84a92
# ╟─63fc39c5-16c1-4b42-99de-feb449776beb
# ╠═3ff725ed-be27-496a-ba11-3734d47e32cd
# ╟─7cccce6c-954f-4720-9dd8-31dcf442208d
# ╟─b524f8b8-64ed-4e9d-b1f5-2e0bd49bdaa4
# ╟─a44c104b-577e-4746-ad94-53e1405a442f
# ╠═8aa91548-24bb-460f-a1eb-072938b96641
# ╟─09f0a38f-64c5-4795-96ba-34e96883f595
# ╟─73d06ecc-baa2-4b5e-9411-a2f8c253efaa
# ╠═be7ece4d-bf25-4c9d-baaa-1aebb696c4e2
# ╠═1aabeba7-5619-4774-96dc-eaab08b5b14a
# ╠═2df030b1-87fd-4eaa-923d-56035a9a825f
# ╟─6bdf5633-319f-4e94-b69e-469cccdd017a
# ╟─c7c1b57e-9f03-4cb6-b141-d4bf15d08905
# ╟─c4d7421e-3046-431f-aabc-cab0f88de453
# ╠═5ebb2e49-79e4-462d-92a6-6c3224c64e36
# ╟─e2f2a55f-9934-45fd-9541-4d30de449cfe
# ╟─2fa9406a-0333-45a6-a674-5a9942aa85d7
# ╠═60293d04-c27a-4c17-8b96-59f9d42ba781
# ╠═7622bbd6-b851-4f59-836d-0f53ca8c610f
# ╠═f5924025-4673-4218-a222-630584a52dc6
# ╟─3afec334-cc16-4f3d-a068-926a31e21ee5
# ╟─d50286ee-6dca-468e-9000-9a32620703ec
# ╟─8c7021cd-d89f-422d-9377-a5fdba5759ab
# ╠═122be89d-bffa-418b-8449-b273bf73ec67
# ╟─486207ce-fd95-410f-9286-75d7b2f21d32
# ╟─45be9e86-390c-4ab7-8d01-c2a2ebb1e245
# ╠═355fd701-c652-4421-a906-1eb1f7506cc4
# ╠═6e8e836b-5069-4de7-abee-34f7bf882a13
# ╠═9e3b40f2-9a8a-43e3-854a-b1b8b3e7be00
# ╟─269dde91-76cf-499b-8a51-a449c7ee7c40
# ╠═d631b3c8-65f5-439a-acb6-b15efe8e4425
# ╠═fd8677c2-0c45-42fc-b7c3-c8c618ce4a1e
# ╠═e3670840-3582-4256-8c3d-a253712c9101
# ╠═ebcb1fbe-6c28-4c58-a856-3c6ac9466b23
# ╠═5a8c91e9-1d1b-4bc7-b1df-96729792a920
# ╠═786e8ba3-1f22-459e-9273-f623262f136a
# ╠═6f8faab9-d447-4e6d-a8cd-19bfeb023dec
# ╠═6cf9cb98-7ee1-4762-81e3-9996f25ceff7
# ╠═8f909b69-739c-40a2-b922-0589a5aac8cf
# ╠═9086dcbc-8033-4e79-a311-eb14438dac4c
# ╠═2b51aaa8-7363-435d-b1f8-5252e85c21f3
# ╠═13191733-be50-4428-a5d1-57145cea221a
# ╠═2991a411-c270-4425-9cf6-1a3b4f759668
# ╟─082cf4fb-382f-4d5c-9ade-08018a66ac87
# ╠═8d31c67f-db73-41b3-9393-16c86c4b05d8
# ╠═890a10df-5beb-41bb-8c90-62de5148bd4c
# ╠═deaad3fd-2eeb-42ab-b1ff-96eec9c6c2fd
# ╟─d6f38052-f243-4d48-8ecb-492758c6cd28
# ╠═23effaf5-b7a3-495b-aad5-209c8c640306
# ╠═f4615674-173f-4245-9e09-1ba49ff2dedb
# ╠═fa324943-c6d9-41c5-9122-d268783bed15

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

# ╔═╡ 27d888e3-7e41-4f0f-b540-6d0c3d479375
begin
	using Pkg
	Pkg.activate(".")

    using PlutoUI, Statistics, CSV, DataFrames, GLM, CairoMakie, HypothesisTests, Colors, MLJBase, DICOM, DICOMUtils, PhantomSegmentation, CalciumScoring, ImageMorphology, ImageFiltering, Noise
    using StatsBase: quantile!, rmsd
end

# ╔═╡ 98367079-bc14-4e8e-a089-5de6a29e39d1
TableOfContents()

# ╔═╡ 551c30ea-1c8c-48b5-98f0-9a75077c2b3b
md"""
## Load DICOMS

All you need to do is set `base_path` once and leave it. After that, the only thing that should change is the `VENDOR`, once for every set, and the `SCAN_NUMBER`, once for each scan.
"""

# ╔═╡ 4b719c89-732b-44ff-9139-8a4e93b49183
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

# ╔═╡ dce447e0-9839-4de1-a413-45d9122bd9d3
md"""
## Helper Functions
"""

# ╔═╡ 067a433e-cad2-4521-aaeb-5f75b0c87ca5
function collect_tuple(tuple_array)
    row_num = size(tuple_array)
    col_num = length(tuple_array[1])
    container = zeros(Int64, row_num..., col_num)
    for i in 1:length(tuple_array)
        container[i, :] = collect(tuple_array[i])
    end
    return container
end

# ╔═╡ 155fc21f-0ff9-4bf0-bda4-3b94f92415bc
function overlay_mask_bind(mask)
    indices = findall(x -> x == 1, mask)
    indices = Tuple.(indices)
    label_array = collect_tuple(indices)
    zs = unique(label_array[:, 3])
    return PlutoUI.Slider(1:length(zs); default=3, show_value=true)
end

# ╔═╡ 0bc01e84-bfaa-43a8-9e4b-3425ed15688e
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

# ╔═╡ b0ceca81-100c-47ba-996b-1f035b71da22
function create_mask(array, mask)
    @assert size(array) == size(mask)
    idxs = findall(x -> x == true, mask)
    overlayed_mask = zeros(size(array))
    for idx in idxs
        overlayed_mask[idx] = array[idx]
    end
    return overlayed_mask
end

# ╔═╡ 516e7416-4dc0-45dc-a47f-6e479a1e17aa
function dilate_mask_large(mask)
    return dilate(mask)
end

# ╔═╡ bf879358-352d-4f41-9264-0f22ae2c6b07
function ring_mask_large(dilated_mask)
    return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask)
end

# ╔═╡ 7fd50487-a6e8-48ad-a60e-f6634bc1a2a5
function dilate_mask_medium(mask)
    return dilate(dilate(dilate(mask)))
end

# ╔═╡ 3db782bf-2a2e-48bf-91ec-e5db53a1dd62
function ring_mask_medium(dilated_mask)
    return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask)
end

# ╔═╡ 7a8032b6-bdd0-48dd-b5bc-d6aa78fcc7dd
function dilate_mask_small(mask)
    return dilate(dilate(dilate(mask)))
end

# ╔═╡ 36da69fa-69ad-4719-be95-c0a5dd918e67
function ring_mask_small(dilated_mask)
    return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask)
end

# ╔═╡ 3fb1c51a-5b80-4b30-9389-541f75d90210
md"""
## Segment Heart
"""

# ╔═╡ 587d3654-1e04-4712-8041-eb5cc3bbaa01
masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3) ÷ 2);

# ╔═╡ afe92437-c623-4b97-8d44-ecd1f9b49462
@bind a PlutoUI.Slider(1:size(masked_array, 3), default=10, show_value=true)

# ╔═╡ bad6f13d-22b2-4727-873f-c0eb4e913742
heatmap(transpose(masked_array[:, :, a]); colormap=:grays)

# ╔═╡ c368787c-75ad-449a-86f8-4cea771cc551
let
    fig = Figure()

    ax = Makie.Axis(fig[1, 1])
    heatmap!(transpose(dcm_array[:, :, 15]); colormap=:grays)
	
    fig
end

# ╔═╡ 82490d92-bef6-4f75-923d-0febc03d032a
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

# ╔═╡ f79369ae-27b6-4216-bddd-0ccb172df6c4
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

# ╔═╡ 7d37e53d-a52d-4e2f-bdc7-ca7bd4827bd7
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

# ╔═╡ 13515784-6b01-46d5-be0a-10ad05b5f7a3
md"""
## Segment Calcium Rod
"""

# ╔═╡ 9064ffe0-1127-4b43-a98a-bd229569b3c3
begin
	thresh = 130
	calcium_image, slice_CCI, quality_slice, cal_rod_slice = mask_rod(
		masked_array, header; calcium_threshold=thresh
	)
end;

# ╔═╡ a60e1ac1-9fc9-4e9e-b09a-527f76c19e3c
@bind c PlutoUI.Slider(axes(calcium_image, 3), default=5, show_value=true)

# ╔═╡ 283975c3-b16a-44f6-b3de-d05af2c4a339
heatmap(transpose(calcium_image[:, :, c]); colormap=:grays)

# ╔═╡ dee99976-34e0-4cfe-b3ff-79ca5a802cbd
md"""
## Load Masks
"""

# ╔═╡ f1f40ec6-5706-485d-9047-7a42abd011b6
mask_L_HD, mask_M_HD, mask_S_HD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_LD, mask_M_LD, mask_S_LD = mask_inserts(
	dcm_array, masked_array, header, slice_CCI, center_insert;
	calcium_threshold=thresh
);

# ╔═╡ 4fb47f99-d0c1-4849-a56f-6c1f92d0c1b7
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

# ╔═╡ cdadb4ce-1423-49df-8c80-d881223557bb
heatmap(transpose(masks); colormap=:grays)

# ╔═╡ ddabff78-b0af-41fe-af5d-ca6fbb1ceadb
md"""
## Prepare Array Masks
"""

# ╔═╡ 6d3c0888-ee9f-4cd4-b7d1-ffa8aae91736
begin
	arr = masked_array[:, :, slice_CCI-2:slice_CCI+2]
	pixel_size = DICOMUtils.get_pixel_size(header)
	voxel_size = pixel_size[1] * pixel_size[2] * pixel_size[3]
end;

# ╔═╡ cf59f276-5799-4cea-9df8-a17e0001ab88
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

# ╔═╡ ce76bf28-5a60-4cd5-a18d-ff0e0da5894e
md"""
## Calibration Prep
"""

# ╔═╡ c1a12d5e-0616-4116-a369-091188f6b95a
begin
	μ, σ = 160, 40
	output = calc_output(masked_array, header, slice_CCI, thresh, trues(3, 3))
	insert_centers = calc_centers(dcm_array, output, header, center_insert, slice_CCI)
	rows, cols = Int(header[tag"Rows"]), Int(header[tag"Columns"])
	mass_cal_factor, angle_0_200HA, water_rod_metrics = mass_calibration(masked_array, insert_centers[:Large_LD], center_insert, cal_rod_slice, rows, cols, pixel_size)
end

# ╔═╡ 3ed8f1b4-ca63-4010-a13f-28e7921a6db1
md"""
# Background
"""

# ╔═╡ 64ddf074-0152-44d6-bc83-e428e5c88cbc
begin
	background_mask = zeros(size(arr)...)
	background_mask[
		(center_insert[1]-5):(center_insert[1]+5),
		(center_insert[2]-5):(center_insert[2]+5),
		2,
	] .= 1
	background_mask = Bool.(background_mask)
	overlayed_bkg_mask = create_mask(arr, background_mask)
	alg2 = SpatiallyWeighted()
	swcs_bkg = score(overlayed_bkg_mask, μ, σ, alg2)
end

# ╔═╡ dccff266-4289-4a59-a627-2316425e346b
@bind q1 overlay_mask_bind(background_mask)

# ╔═╡ c4d68bf6-5704-43e5-bbc8-f551808ac8af
overlay_mask_plot(arr, background_mask, q1, "dilated mask")

# ╔═╡ d9b028e0-3217-4725-a172-1a8ba7ef3d2c
md"""
# Score Large Inserts
"""

# ╔═╡ 61e8a95d-cc58-4ee4-8384-777b26a8faa6
md"""
## High Density
"""

# ╔═╡ 3062d461-bc1e-459c-ae8e-136949f1b45e
md"""
#### Dilated mask
"""

# ╔═╡ 9eb2ae1f-7802-43cd-9df6-6ed13869fbb6
@bind g2 overlay_mask_bind(dilated_mask_L_HD)

# ╔═╡ dcb2146d-7cbe-4f6c-aea3-9bc1e894ed08
overlay_mask_plot(arr, dilated_mask_L_HD, g2, "dilated mask")

# ╔═╡ 028e8039-059b-480e-a606-a9a06fb01d6a
begin
	alg = Agatston()
	overlayed_mask_l_hd = create_mask(arr, dilated_mask_L_HD)
	agat_l_hd, mass_l_hd = score(
		overlayed_mask_l_hd,
		pixel_size,
		mass_cal_factor,
		alg;
		kV=kV
	)
end

# ╔═╡ 039c94dd-2d13-41a9-b729-f686d9eb7b91
md"""
## Medium Density
"""

# ╔═╡ e5adc836-ed1e-42ce-946a-260c65876436
md"""
#### Dilated mask
"""

# ╔═╡ e760a64c-f314-4643-a9c3-d7dda51af515
@bind h2 overlay_mask_bind(dilated_mask_L_MD)

# ╔═╡ 4125229f-dbde-400e-836d-967c43909dec
overlay_mask_plot(arr, dilated_mask_L_MD, h2, "dilated mask")

# ╔═╡ dc23071b-292c-4b23-9098-2658de2f8774
begin
	overlayed_mask_l_md = create_mask(arr, dilated_mask_L_MD)
	agat_l_md, mass_l_md = score(
		overlayed_mask_l_md,
		pixel_size,
		mass_cal_factor,
		alg;
		kV=kV
	)
end

# ╔═╡ ca085092-fb1a-4eeb-8231-a87f70c6938f
md"""
## Low Density
"""

# ╔═╡ 46570e16-be18-4d05-9353-02cfb839b05d
md"""
#### Dilated mask
"""

# ╔═╡ a73cf5b4-428b-498c-8fe4-a14ada38b4e8
@bind j2 overlay_mask_bind(dilated_mask_L_LD)

# ╔═╡ 65b3cf51-bd5e-4579-9d9a-d6b350b28e3d
overlay_mask_plot(arr, dilated_mask_L_LD, j2, "dilated mask")

# ╔═╡ 8ac0b3da-f98b-4e92-8e85-6916e1e4cdfc
begin
	overlayed_mask_l_ld = create_mask(arr, dilated_mask_L_LD)
	agat_l_ld, mass_l_ld = score(
		overlayed_mask_l_ld,
		pixel_size,
		mass_cal_factor,
		alg;
		kV=kV
	)
end

# ╔═╡ b55fadec-8289-4efe-8967-228a07673978
mass_l_ld

# ╔═╡ 696796b8-4390-46b3-8558-cae127cfc505
md"""
# Score Medium Inserts
"""

# ╔═╡ 1a559f19-ec63-4e9d-a020-25ec2e542f7a
md"""
## High Density
"""

# ╔═╡ 89ca7ecd-fccf-444c-9c0f-1c6ec38ab473
md"""
#### Dilated mask
"""

# ╔═╡ 3c1984cb-0e03-4625-aa03-37b96e0d23d8
@bind z2 overlay_mask_bind(dilated_mask_M_HD)

# ╔═╡ a6a809fd-6784-4e99-ab06-7c3529b1be8e
overlay_mask_plot(arr, dilated_mask_M_HD, z2, "dilated mask")

# ╔═╡ 64abde3a-09b2-44d3-be9b-53f16792ed65
begin
	overlayed_mask_m_hd = create_mask(arr, dilated_mask_M_HD)
	agat_m_hd, mass_m_hd = score(
		overlayed_mask_m_hd,
		pixel_size,
		mass_cal_factor,
		alg;
		kV=kV
	)
end

# ╔═╡ c03d3856-ca48-4d25-bbc0-9d685d9ae30f
md"""
## Medium Density
"""

# ╔═╡ b8522dce-d3a0-43b6-8440-83f0330430c3
md"""
#### Dilated mask
"""

# ╔═╡ 49ef8688-3151-4e05-97ed-3fe7a478fac6
@bind x2 overlay_mask_bind(dilated_mask_M_MD)

# ╔═╡ fefb7123-80f3-4551-9a1d-5152976273d9
overlay_mask_plot(arr, dilated_mask_M_MD, x2, "dilated mask")

# ╔═╡ 4be15540-666f-4671-97a4-f3a774378c62
begin
	overlayed_mask_m_md = create_mask(arr, dilated_mask_M_MD)
	agat_m_md, mass_m_md = score(
		overlayed_mask_m_md,
		pixel_size,
		mass_cal_factor,
		alg;
		kV=kV
	)
end

# ╔═╡ 35e26834-9ad6-40a2-8d84-e4b3d03f6d00
md"""
## Low Density
"""

# ╔═╡ 0954e8b3-d9b9-43a8-8ba4-c71d02d30592
md"""
#### Dilated mask
"""

# ╔═╡ 0471936a-df05-445c-9784-bb8b7cc0bbe4
@bind y2 overlay_mask_bind(dilated_mask_M_LD)

# ╔═╡ 3cb639a9-088a-4474-895d-46c2eff5a742
overlay_mask_plot(arr, dilated_mask_M_LD, y2, "dilated mask")

# ╔═╡ f31754b2-1d96-479d-acba-1ae3a52315a2
begin
	overlayed_mask_m_ld = create_mask(arr, dilated_mask_M_LD)
	agat_m_ld, mass_m_ld = score(
		overlayed_mask_m_ld,
		pixel_size,
		mass_cal_factor,
		alg;
		kV=kV
	)
end

# ╔═╡ 0b6a7aa9-c1ed-4c05-9b98-2f3ea6b962b5
md"""
# Score Small Inserts
"""

# ╔═╡ 8350eb00-54ab-4265-81e6-97f3283ac07e
md"""
## High Density
"""

# ╔═╡ 48cb3c65-b3a1-482f-8e99-2db0d7d50b00
md"""
#### Dilated mask
"""

# ╔═╡ 4d30f03e-0e16-46c6-8c9f-5cf67a16804d
@bind l2 overlay_mask_bind(dilated_mask_S_HD)

# ╔═╡ 30c71d60-28a4-49ee-b68d-600b15801e62
overlay_mask_plot(arr, dilated_mask_S_HD, l2, "dilated mask")

# ╔═╡ d4548c3c-b3fd-4637-80b2-29998ad5a48d
begin
	overlayed_mask_s_hd = create_mask(arr, dilated_mask_S_HD)
	agat_s_hd, mass_s_hd = score(
		overlayed_mask_s_hd,
		pixel_size,
		mass_cal_factor,
		alg;
		kV=kV
	)
end

# ╔═╡ 03eb1cb3-c1fd-4b88-997e-01aa622c4c6a
md"""
## Medium Density
"""

# ╔═╡ d76226d0-d56c-4ee6-841f-cbf08311c1b2
md"""
#### Dilated mask
"""

# ╔═╡ aaf29e1d-81cf-4156-955c-07e98325988a
@bind k2 overlay_mask_bind(dilated_mask_S_MD)

# ╔═╡ 3c56dad7-c08f-44ef-b0d8-c7e3f30c5595
overlay_mask_plot(arr, dilated_mask_S_MD, k2, "dilated mask")

# ╔═╡ 9a36d070-222d-4a5f-aeb0-07eae3a3fd75
begin
	overlayed_mask_s_md = create_mask(arr, dilated_mask_S_MD)
	agat_s_md, mass_s_md = score(
		overlayed_mask_s_md,
		pixel_size,
		mass_cal_factor,
		alg;
		kV=kV
	)
end

# ╔═╡ 7d7e1f64-2bbc-48e5-bbb1-b96767a16ff5
md"""
## Low Density
"""

# ╔═╡ 71a8f052-56d6-43be-b926-b60f42cd430f
md"""
#### Dilated mask
"""

# ╔═╡ a7dbdade-9a8e-4282-b296-4f2bfede1d57
@bind m2 overlay_mask_bind(dilated_mask_S_LD)

# ╔═╡ 86907265-3a6d-4d85-9480-1cc2a62041d8
overlay_mask_plot(arr, dilated_mask_S_LD, m2, "dilated mask")

# ╔═╡ 41028c8e-b49b-43cd-9b97-532df3521674
begin
	overlayed_mask_s_ld = create_mask(arr, dilated_mask_S_LD)
	agat_s_ld, mass_s_ld = score(
		overlayed_mask_s_ld,
		pixel_size,
		mass_cal_factor,
		alg;
		kV=kV
	)
end

# ╔═╡ b5c4eb96-3531-492e-8c0d-741a30d174e0
md"""
# Results
"""

# ╔═╡ 02570396-5153-4460-827f-884573001ea3
PhantomSegmentation.get_pixel_size(header)

# ╔═╡ 8668c62c-4a56-4a52-b7da-6e5ae2cafac3
inserts = ["Low Density", "Medium Density", "High Density"]

# ╔═╡ ffd9ba0c-70ed-43bc-b350-eb1d7c8e8818
volume_gt = [7.065, 63.585, 176.625]

# ╔═╡ 08c9c09b-dacf-4682-bc6e-138546dab56b
density_array = [0.200, 0.400, 0.800]

# ╔═╡ 11751ee9-a4d3-48be-9dd0-e85e25a9e76e
ground_truth_mass_large = [
    volume_gt[3] * density_array[1],
    volume_gt[3] * density_array[2],
    volume_gt[3] * density_array[3],
] # mg

# ╔═╡ 815687b4-5d83-47b8-b534-03a261cee32f
calculated_mass_large = [mass_l_ld, mass_l_md, mass_l_hd]

# ╔═╡ dc918132-fd8c-4120-a01c-f8baecf2c6a2
ground_truth_mass_medium = [
    volume_gt[2] * density_array[1],
    volume_gt[2] * density_array[2],
    volume_gt[2] * density_array[3],
] # mg

# ╔═╡ 906410b4-64a3-4873-899b-7031bf05a501
calculated_mass_medium = [mass_m_ld, mass_m_md, mass_m_hd]

# ╔═╡ e164e898-4e68-4c10-9916-f42520c1a132
ground_truth_mass_small = [
    volume_gt[1] * density_array[1],
    volume_gt[1] * density_array[2],
    volume_gt[1] * density_array[3],
] # mg

# ╔═╡ eb2d5ad2-19de-494e-bdc2-bb9a38037672
calculated_mass_small = [mass_s_ld, mass_s_md, mass_s_hd]

# ╔═╡ 56f599d6-1cd3-45fc-a252-de2bdd087fdc
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

# ╔═╡ 8216bfb0-0416-4c3b-b0da-8cdce7abb88a
md"""
### Save Results
"""

# ╔═╡ f5d3b0c2-bbcb-4087-8958-b37d05ccf112
# if ~isdir(string(cd(pwd, "..") , "/output/", TYPE))
# 	mkdir(string(cd(pwd, "..") , "/output/", TYPE))
# end

# ╔═╡ d270a8a1-4d0d-4a16-869a-94e451e54248
# output_path = string(cd(pwd, "..") , "/output/", TYPE, "/", scan, ".csv")

# ╔═╡ eb7b01e6-71eb-45af-96ac-7740523da5d3
# CSV.write(output_path, df)

# ╔═╡ e763c1e7-0cf6-4aa5-9085-3ce8260d841b
md"""
### Save full df
"""

# ╔═╡ 2b633e1b-4bac-4ace-ad35-ae5b276f53e3
dfs = []

# ╔═╡ 9c1c7cbc-35b3-4c42-868a-dddbe978d99f
push!(dfs, df)

# ╔═╡ 84bd8bf6-a748-43de-a009-6a7009eaf2dc
# if length(dfs) == 24
# 	global new_df = vcat(dfs[1:24]...)
# 	output_path_new = string(cd(pwd, "..") , "/output/", TYPE, "/", "full.csv")
# 	CSV.write(output_path_new, new_df)
# end

# ╔═╡ Cell order:
# ╠═27d888e3-7e41-4f0f-b540-6d0c3d479375
# ╠═98367079-bc14-4e8e-a089-5de6a29e39d1
# ╟─551c30ea-1c8c-48b5-98f0-9a75077c2b3b
# ╠═4b719c89-732b-44ff-9139-8a4e93b49183
# ╟─dce447e0-9839-4de1-a413-45d9122bd9d3
# ╟─067a433e-cad2-4521-aaeb-5f75b0c87ca5
# ╟─155fc21f-0ff9-4bf0-bda4-3b94f92415bc
# ╟─0bc01e84-bfaa-43a8-9e4b-3425ed15688e
# ╟─b0ceca81-100c-47ba-996b-1f035b71da22
# ╟─516e7416-4dc0-45dc-a47f-6e479a1e17aa
# ╟─bf879358-352d-4f41-9264-0f22ae2c6b07
# ╟─7fd50487-a6e8-48ad-a60e-f6634bc1a2a5
# ╟─3db782bf-2a2e-48bf-91ec-e5db53a1dd62
# ╟─7a8032b6-bdd0-48dd-b5bc-d6aa78fcc7dd
# ╟─36da69fa-69ad-4719-be95-c0a5dd918e67
# ╟─3fb1c51a-5b80-4b30-9389-541f75d90210
# ╠═587d3654-1e04-4712-8041-eb5cc3bbaa01
# ╟─afe92437-c623-4b97-8d44-ecd1f9b49462
# ╟─bad6f13d-22b2-4727-873f-c0eb4e913742
# ╟─c368787c-75ad-449a-86f8-4cea771cc551
# ╟─82490d92-bef6-4f75-923d-0febc03d032a
# ╟─f79369ae-27b6-4216-bddd-0ccb172df6c4
# ╟─7d37e53d-a52d-4e2f-bdc7-ca7bd4827bd7
# ╟─13515784-6b01-46d5-be0a-10ad05b5f7a3
# ╠═9064ffe0-1127-4b43-a98a-bd229569b3c3
# ╟─a60e1ac1-9fc9-4e9e-b09a-527f76c19e3c
# ╟─283975c3-b16a-44f6-b3de-d05af2c4a339
# ╟─dee99976-34e0-4cfe-b3ff-79ca5a802cbd
# ╠═f1f40ec6-5706-485d-9047-7a42abd011b6
# ╠═4fb47f99-d0c1-4849-a56f-6c1f92d0c1b7
# ╟─cdadb4ce-1423-49df-8c80-d881223557bb
# ╟─ddabff78-b0af-41fe-af5d-ca6fbb1ceadb
# ╠═6d3c0888-ee9f-4cd4-b7d1-ffa8aae91736
# ╠═cf59f276-5799-4cea-9df8-a17e0001ab88
# ╟─ce76bf28-5a60-4cd5-a18d-ff0e0da5894e
# ╠═c1a12d5e-0616-4116-a369-091188f6b95a
# ╟─3ed8f1b4-ca63-4010-a13f-28e7921a6db1
# ╠═64ddf074-0152-44d6-bc83-e428e5c88cbc
# ╠═dccff266-4289-4a59-a627-2316425e346b
# ╠═c4d68bf6-5704-43e5-bbc8-f551808ac8af
# ╠═d9b028e0-3217-4725-a172-1a8ba7ef3d2c
# ╠═61e8a95d-cc58-4ee4-8384-777b26a8faa6
# ╠═3062d461-bc1e-459c-ae8e-136949f1b45e
# ╠═9eb2ae1f-7802-43cd-9df6-6ed13869fbb6
# ╠═dcb2146d-7cbe-4f6c-aea3-9bc1e894ed08
# ╠═028e8039-059b-480e-a606-a9a06fb01d6a
# ╠═039c94dd-2d13-41a9-b729-f686d9eb7b91
# ╠═e5adc836-ed1e-42ce-946a-260c65876436
# ╠═e760a64c-f314-4643-a9c3-d7dda51af515
# ╠═4125229f-dbde-400e-836d-967c43909dec
# ╠═dc23071b-292c-4b23-9098-2658de2f8774
# ╠═ca085092-fb1a-4eeb-8231-a87f70c6938f
# ╠═46570e16-be18-4d05-9353-02cfb839b05d
# ╠═a73cf5b4-428b-498c-8fe4-a14ada38b4e8
# ╠═65b3cf51-bd5e-4579-9d9a-d6b350b28e3d
# ╠═8ac0b3da-f98b-4e92-8e85-6916e1e4cdfc
# ╠═b55fadec-8289-4efe-8967-228a07673978
# ╠═696796b8-4390-46b3-8558-cae127cfc505
# ╠═1a559f19-ec63-4e9d-a020-25ec2e542f7a
# ╠═89ca7ecd-fccf-444c-9c0f-1c6ec38ab473
# ╠═3c1984cb-0e03-4625-aa03-37b96e0d23d8
# ╠═a6a809fd-6784-4e99-ab06-7c3529b1be8e
# ╠═64abde3a-09b2-44d3-be9b-53f16792ed65
# ╠═c03d3856-ca48-4d25-bbc0-9d685d9ae30f
# ╠═b8522dce-d3a0-43b6-8440-83f0330430c3
# ╠═49ef8688-3151-4e05-97ed-3fe7a478fac6
# ╠═fefb7123-80f3-4551-9a1d-5152976273d9
# ╠═4be15540-666f-4671-97a4-f3a774378c62
# ╠═35e26834-9ad6-40a2-8d84-e4b3d03f6d00
# ╠═0954e8b3-d9b9-43a8-8ba4-c71d02d30592
# ╠═0471936a-df05-445c-9784-bb8b7cc0bbe4
# ╠═3cb639a9-088a-4474-895d-46c2eff5a742
# ╠═f31754b2-1d96-479d-acba-1ae3a52315a2
# ╠═0b6a7aa9-c1ed-4c05-9b98-2f3ea6b962b5
# ╠═8350eb00-54ab-4265-81e6-97f3283ac07e
# ╠═48cb3c65-b3a1-482f-8e99-2db0d7d50b00
# ╠═4d30f03e-0e16-46c6-8c9f-5cf67a16804d
# ╠═30c71d60-28a4-49ee-b68d-600b15801e62
# ╠═d4548c3c-b3fd-4637-80b2-29998ad5a48d
# ╠═03eb1cb3-c1fd-4b88-997e-01aa622c4c6a
# ╠═d76226d0-d56c-4ee6-841f-cbf08311c1b2
# ╠═aaf29e1d-81cf-4156-955c-07e98325988a
# ╠═3c56dad7-c08f-44ef-b0d8-c7e3f30c5595
# ╠═9a36d070-222d-4a5f-aeb0-07eae3a3fd75
# ╠═7d7e1f64-2bbc-48e5-bbb1-b96767a16ff5
# ╠═71a8f052-56d6-43be-b926-b60f42cd430f
# ╠═a7dbdade-9a8e-4282-b296-4f2bfede1d57
# ╠═86907265-3a6d-4d85-9480-1cc2a62041d8
# ╠═41028c8e-b49b-43cd-9b97-532df3521674
# ╠═b5c4eb96-3531-492e-8c0d-741a30d174e0
# ╠═02570396-5153-4460-827f-884573001ea3
# ╠═8668c62c-4a56-4a52-b7da-6e5ae2cafac3
# ╠═ffd9ba0c-70ed-43bc-b350-eb1d7c8e8818
# ╠═08c9c09b-dacf-4682-bc6e-138546dab56b
# ╠═11751ee9-a4d3-48be-9dd0-e85e25a9e76e
# ╠═815687b4-5d83-47b8-b534-03a261cee32f
# ╠═dc918132-fd8c-4120-a01c-f8baecf2c6a2
# ╠═906410b4-64a3-4873-899b-7031bf05a501
# ╠═e164e898-4e68-4c10-9916-f42520c1a132
# ╠═eb2d5ad2-19de-494e-bdc2-bb9a38037672
# ╠═56f599d6-1cd3-45fc-a252-de2bdd087fdc
# ╠═8216bfb0-0416-4c3b-b0da-8cdce7abb88a
# ╠═f5d3b0c2-bbcb-4087-8958-b37d05ccf112
# ╠═d270a8a1-4d0d-4a16-869a-94e451e54248
# ╠═eb7b01e6-71eb-45af-96ac-7740523da5d3
# ╠═e763c1e7-0cf6-4aa5-9085-3ce8260d841b
# ╠═2b633e1b-4bac-4ace-ad35-ae5b276f53e3
# ╠═9c1c7cbc-35b3-4c42-868a-dddbe978d99f
# ╠═84bd8bf6-a748-43de-a009-6a7009eaf2dc

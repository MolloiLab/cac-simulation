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

# ╔═╡ 742070ad-ad57-4153-abbb-7ef7d65c6ac3
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.activate(".")

    using PlutoUI, Statistics, CSV, DataFrames, GLM, CairoMakie, HypothesisTests, Colors, MLJBase, DICOM, DICOMUtils, PhantomSegmentation, CalciumScoring, ImageMorphology, ImageFiltering, Noise
    using StatsBase: quantile!, rmsd
end

# ╔═╡ 2bf6a8ad-8203-4543-86bf-668dc437e25f
TableOfContents()

# ╔═╡ 4b417377-dfa2-4bdc-b64b-c0543f487c24
md"""
## Load DICOMS

All you need to do is set `base_path` once and leave it. After that, the only thing that should change is the `VENDOR`, once for every set, and the `SCAN_NUMBER`, once for each scan.
"""

# ╔═╡ ac86dfe5-30c7-4d22-adae-152763354408
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

# ╔═╡ 9c998e5a-0f36-4b61-ad70-069ed600808c
md"""
## Helper Functions
"""

# ╔═╡ 9bdedef5-2b33-46b7-bed1-f47637b08eef
function collect_tuple(tuple_array)
    row_num = size(tuple_array)
    col_num = length(tuple_array[1])
    container = zeros(Int64, row_num..., col_num)
    for i in 1:length(tuple_array)
        container[i, :] = collect(tuple_array[i])
    end
    return container
end

# ╔═╡ 6f79a0d3-7a22-4047-a679-b215713f4bfa
function overlay_mask_bind(mask)
    indices = findall(x -> x == 1, mask)
    indices = Tuple.(indices)
    label_array = collect_tuple(indices)
    zs = unique(label_array[:, 3])
    return PlutoUI.Slider(1:length(zs); default=3, show_value=true)
end

# ╔═╡ c43e6414-e329-4a27-9ba4-9a9f3b1d94ae
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

# ╔═╡ ff73c57d-f1cf-4f7e-80d1-b9b0e1a0a90f
function create_mask(array, mask)
    @assert size(array) == size(mask)
    idxs = findall(x -> x == true, mask)
    overlayed_mask = zeros(size(array))
    for idx in idxs
        overlayed_mask[idx] = array[idx]
    end
    return overlayed_mask
end

# ╔═╡ 7cb015df-3b20-4168-8a68-0956dc9256ff
function dilate_mask_large(mask)
    return dilate(mask)
end

# ╔═╡ 50ab516f-f988-45c8-a780-30b678b7cfe6
function ring_mask_large(dilated_mask)
    return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask)
end

# ╔═╡ 74ddfe10-c39e-45bc-b1e7-3335c293c123
function dilate_mask_medium(mask)
    return dilate(dilate(dilate(mask)))
end

# ╔═╡ 4812bf96-7699-4dc7-8c7d-b408705e7f31
function ring_mask_medium(dilated_mask)
    return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask)
end

# ╔═╡ f9b901e4-6fd6-453b-9eda-a31a1601c196
function dilate_mask_small(mask)
    return dilate(dilate(dilate(mask)))
end

# ╔═╡ d439d790-8c2f-40dd-bef7-72c1baf0848e
function ring_mask_small(dilated_mask)
    return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask)
end

# ╔═╡ 4e5c0dd3-f9fe-461d-b526-8314a048c1f0
md"""
## Segment Heart
"""

# ╔═╡ b11ec843-d599-4c40-a258-28d87699d7c5
masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3) ÷ 2);

# ╔═╡ 05543739-3299-4556-b9eb-ba7366cc5c1d
@bind a PlutoUI.Slider(1:size(masked_array, 3), default=10, show_value=true)

# ╔═╡ 42ad7c3c-521c-4d91-a236-2b56dd3ab674
heatmap(transpose(masked_array[:, :, a]); colormap=:grays)

# ╔═╡ e40c8ac9-60e7-4caa-9ed4-3647c302b47d
let
    fig = Figure()

    ax = Makie.Axis(fig[1, 1])
    heatmap!(transpose(dcm_array[:, :, 15]); colormap=:grays)
	
    fig
end

# ╔═╡ ccfb5f8c-0d47-4d95-b84e-2afa9fd0e14a
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

# ╔═╡ 764351f3-841e-4d88-b537-f08b711f5f87
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

# ╔═╡ a00d562f-a0aa-45c9-aaf0-7cedc1c02953
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

# ╔═╡ 4f0cd797-150a-40fe-ae70-d631ef68c627
md"""
## Segment Calcium Rod
"""

# ╔═╡ f05ae4e1-a6e7-4351-bf15-426cab19cd10
begin
	thresh = 130
	calcium_image, slice_CCI, quality_slice, cal_rod_slice = mask_rod(
		masked_array, header; calcium_threshold=thresh
	)
end;

# ╔═╡ 4e6aaf97-58d9-4d8c-9e3f-1cb39acfb06f
@bind c PlutoUI.Slider(axes(calcium_image, 3), default=5, show_value=true)

# ╔═╡ 98774ca0-e8a1-4d4d-b333-fa9236aa8f96
heatmap(transpose(calcium_image[:, :, c]); colormap=:grays)

# ╔═╡ 0086df28-73ae-4b69-ac2a-48fb946ae3b4
md"""
## Load Masks
"""

# ╔═╡ 6a8b205c-a41f-4e6c-aaf2-b27137dd041b
mask_L_HD, mask_M_HD, mask_S_HD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_LD, mask_M_LD, mask_S_LD = mask_inserts(
	dcm_array, masked_array, header, slice_CCI, center_insert;
	calcium_threshold=thresh
);

# ╔═╡ ca7dbfa8-bc4d-410e-8276-e9afbdfdbcf6
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

# ╔═╡ 763b44dc-7323-46e0-b80c-a3b939986483
heatmap(transpose(masks); colormap=:grays)

# ╔═╡ 87b784cd-1dfa-43e3-a0f1-8965330d2af3
md"""
## Prepare Array Masks
"""

# ╔═╡ e816098b-a584-4301-8fcf-5849a4b6916c
begin
	arr = masked_array[:, :, slice_CCI-2:slice_CCI+2]
	pixel_size = DICOMUtils.get_pixel_size(header)
	voxel_size = pixel_size[1] * pixel_size[2] * pixel_size[3]
end;

# ╔═╡ 78d1716c-af9a-47da-a320-70042116d752
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

# ╔═╡ 4cce344a-693d-457c-8c4b-7c479e72f2ab
md"""
## Calibration Prep
"""

# ╔═╡ 9eab7433-c807-4d61-8178-a8400a8f7e24
begin
	output = calc_output(masked_array, header, slice_CCI, thresh, trues(3, 3))
	insert_centers = calc_centers(dcm_array, output, header, center_insert, slice_CCI)
	rows, cols = Int(header[tag"Rows"]), Int(header[tag"Columns"])
	mass_cal_factor, angle_0_200HA, water_rod_metrics = mass_calibration(masked_array, insert_centers[:Large_LD], center_insert, cal_rod_slice, rows, cols, pixel_size)
end

# ╔═╡ 583cc301-29f8-48e7-b208-281a2c552745
md"""
# Background
"""

# ╔═╡ aa6a1783-ea13-455c-9118-663b8db76820
md"""
# Score Large Inserts
"""

# ╔═╡ 2fd2f313-5651-429c-a6e0-25c5da0968d7
md"""
## High Density
"""

# ╔═╡ 8705619d-7acf-4246-9eb7-658cf7f71547
md"""
#### Dilated mask
"""

# ╔═╡ 3c78c1c0-8fe6-49d3-81a4-84100b461a14
@bind g2 overlay_mask_bind(dilated_mask_L_HD)

# ╔═╡ 2d426f83-904e-4fac-a0a6-37c286a33a7a
overlay_mask_plot(arr, dilated_mask_L_HD, g2, "dilated mask")

# ╔═╡ 76c6c43e-cc4d-43aa-b557-4c0c6860a071
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

# ╔═╡ 92710a4d-7f1d-4145-b56c-d18836ddb438
begin
	background_mask = zeros(size(arr)...)
	background_mask[
		(center_insert[1]-5):(center_insert[1]+5),
		(center_insert[2]-5):(center_insert[2]+5),
		2,
	] .= 1
	background_mask = Bool.(background_mask)
	overlayed_bkg_mask = create_mask(arr, background_mask)
	agat_bkg, mass_bkg = score(
		overlayed_bkg_mask,
		pixel_size,
		mass_cal_factor,
		alg;
		kV=kV
	)
end

# ╔═╡ f4fd148b-9e95-4826-baa7-e736e9b883a0
@bind q1 overlay_mask_bind(background_mask)

# ╔═╡ b99c68ae-86e3-40c2-8b7e-847d9f3fc402
overlay_mask_plot(arr, background_mask, q1, "dilated mask")

# ╔═╡ e04673e6-5e70-4944-8c40-9162df11e7e9
md"""
## Medium Density
"""

# ╔═╡ 012a1ca7-24d3-447f-a0a4-c8550150936a
md"""
#### Dilated mask
"""

# ╔═╡ 9fd4245a-0098-41af-8ae4-6c7fc502c4e5
@bind h2 overlay_mask_bind(dilated_mask_L_MD)

# ╔═╡ 5a77d2b9-259d-4f95-bfdf-a2f09ed42cbb
overlay_mask_plot(arr, dilated_mask_L_MD, h2, "dilated mask")

# ╔═╡ 6f46e044-50b0-4d78-b8ff-895e1b575331
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

# ╔═╡ 2363391c-0cef-4744-a182-5a1553d25d06
md"""
## Low Density
"""

# ╔═╡ 8ec96204-8a2a-4143-b999-6a0a5efc7199
md"""
#### Dilated mask
"""

# ╔═╡ 837da46d-c1e5-43cf-ab7b-642a17a2f764
@bind j2 overlay_mask_bind(dilated_mask_L_LD)

# ╔═╡ bbf4b093-455d-4097-9d7a-0aeaa21cdeca
overlay_mask_plot(arr, dilated_mask_L_LD, j2, "dilated mask")

# ╔═╡ 0c37fdc1-25f9-4c33-a7d1-706210b8ef95
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

# ╔═╡ 6651c12d-2c0f-4699-906f-1d4089422bf0
mass_l_ld

# ╔═╡ 400ffd0c-084b-470e-807f-1f8ad045d631
md"""
# Score Medium Inserts
"""

# ╔═╡ 9c1401a2-5f58-41f3-a869-c17ebbf98a64
md"""
## High Density
"""

# ╔═╡ f8781db1-8898-484f-acd8-bd84c0aaac59
md"""
#### Dilated mask
"""

# ╔═╡ 2a2dced3-055c-44a7-a620-f3e7a8522c68
@bind z2 overlay_mask_bind(dilated_mask_M_HD)

# ╔═╡ e199a541-c83e-45ba-9650-97d6831aba55
overlay_mask_plot(arr, dilated_mask_M_HD, z2, "dilated mask")

# ╔═╡ a1df38e8-f8bb-47bb-98c5-47be9cd5c8a6
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

# ╔═╡ 48114b90-5c6c-4d24-8d6b-c6ba5e4ae38c
md"""
## Medium Density
"""

# ╔═╡ 317d5fa2-3aac-480b-abda-dbbf824b8b75
md"""
#### Dilated mask
"""

# ╔═╡ 4d91ed1d-2caf-463c-b80a-893614cc3c44
@bind x2 overlay_mask_bind(dilated_mask_M_MD)

# ╔═╡ c469b6fa-e11d-424b-81ea-2512f6bc1b9e
overlay_mask_plot(arr, dilated_mask_M_MD, x2, "dilated mask")

# ╔═╡ 78cbca0a-78eb-4101-9e63-565300173df6
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

# ╔═╡ b83ddade-0e5b-447f-801a-dfd81b2d998f
md"""
## Low Density
"""

# ╔═╡ c27a5300-0be0-479c-a368-cc81c8991291
md"""
#### Dilated mask
"""

# ╔═╡ a963b7eb-d5b9-461d-9f4d-bf94596eb96e
@bind y2 overlay_mask_bind(dilated_mask_M_LD)

# ╔═╡ a170c6b0-5673-4ec5-a4ff-4c697135f94c
overlay_mask_plot(arr, dilated_mask_M_LD, y2, "dilated mask")

# ╔═╡ c63cf606-c3f4-41ad-b877-0a87b15759d7
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

# ╔═╡ 8beb225c-4285-4e8c-9fdc-14d8e93011c2
md"""
# Score Small Inserts
"""

# ╔═╡ e2e9d93d-0d34-4236-9e3f-74ed77df4f38
md"""
## High Density
"""

# ╔═╡ 2b87600b-1c67-4987-a37e-c9860a52fb80
md"""
#### Dilated mask
"""

# ╔═╡ 8e395a11-5572-4a96-b2f4-5954c92b58a6
@bind l2 overlay_mask_bind(dilated_mask_S_HD)

# ╔═╡ 6a600411-fd67-4cd8-b38f-7a89a5914749
overlay_mask_plot(arr, dilated_mask_S_HD, l2, "dilated mask")

# ╔═╡ 54aad26c-51c7-42cb-972a-ed61c37ea8d3
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

# ╔═╡ 660ec03c-8f6b-439d-8600-6f6248ad6f54
md"""
## Medium Density
"""

# ╔═╡ 9d07d8a3-6052-4992-a34b-773eebce7241
md"""
#### Dilated mask
"""

# ╔═╡ e19546f2-ea12-47d8-9398-72daa619b8c0
@bind k2 overlay_mask_bind(dilated_mask_S_MD)

# ╔═╡ b2199735-aa64-4a54-9422-6e94ab59cc38
overlay_mask_plot(arr, dilated_mask_S_MD, k2, "dilated mask")

# ╔═╡ 43fb1435-b69f-487d-bcaf-1f8416e3dcd6
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

# ╔═╡ 15e994bc-c4ee-4b89-b3c7-23e5874291b0
md"""
## Low Density
"""

# ╔═╡ 32c98b0c-e55e-4cf9-a393-bf5ef777799e
md"""
#### Dilated mask
"""

# ╔═╡ 672b519f-55bd-4456-9d71-8db8861e9710
@bind m2 overlay_mask_bind(dilated_mask_S_LD)

# ╔═╡ 22f4f737-ad6d-492d-86c2-7155ec3f07d5
overlay_mask_plot(arr, dilated_mask_S_LD, m2, "dilated mask")

# ╔═╡ 4d37ec3f-ed5a-48df-8dfd-9bd7ad8b5714
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

# ╔═╡ 66d72b92-3975-4d01-8a7e-deaa58720044
md"""
# Results
"""

# ╔═╡ 3bccc7f7-8230-4495-82fb-23113e69a91a
PhantomSegmentation.get_pixel_size(header)

# ╔═╡ 6c547b31-6a36-4ffa-a8d8-bbdcc8b9bf83
inserts = ["Low Density", "Medium Density", "High Density"]

# ╔═╡ c320a097-5d86-4af7-83a0-4eb8407a314c
volume_gt = [7.065, 63.585, 176.625]

# ╔═╡ d39490dc-cc2c-46dd-87ab-150d420cced6
density_array = [0.200, 0.400, 0.800]

# ╔═╡ 413400a1-de67-4f52-85e0-fb8b1fa9f1aa
ground_truth_mass_large = [
    volume_gt[3] * density_array[1],
    volume_gt[3] * density_array[2],
    volume_gt[3] * density_array[3],
] # mg

# ╔═╡ 832dde61-330c-41e8-afa4-4586d22def72
calculated_mass_large = [mass_l_ld, mass_l_md, mass_l_hd]

# ╔═╡ 407e21cf-7cc0-454e-8e55-f97abc3e6c3c
ground_truth_mass_medium = [
    volume_gt[2] * density_array[1],
    volume_gt[2] * density_array[2],
    volume_gt[2] * density_array[3],
] # mg

# ╔═╡ f8e648d6-c42c-472d-b771-22bcfad39f87
calculated_mass_medium = [mass_m_ld, mass_m_md, mass_m_hd]

# ╔═╡ 8a1602e8-44fe-45cf-be29-10accca25519
ground_truth_mass_small = [
    volume_gt[1] * density_array[1],
    volume_gt[1] * density_array[2],
    volume_gt[1] * density_array[3],
] # mg

# ╔═╡ 5499d8b7-35b2-41c6-8a4e-321654897c24
calculated_mass_small = [mass_s_ld, mass_s_md, mass_s_hd]

# ╔═╡ 63cb2095-d00e-4e50-8c4f-c7c50f047bdd
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

# ╔═╡ b7cea572-e416-4ee4-9342-73ded1ad133e
md"""
### Save Results
"""

# ╔═╡ 6a2d949f-4083-4605-9455-5d7b9974ac4a
# if ~isdir(string(cd(pwd, "..") , "/output/", TYPE))
# 	mkdir(string(cd(pwd, "..") , "/output/", TYPE))
# end

# ╔═╡ 0e00f908-458f-4cf3-aa79-de46aed9fb3f
# output_path = string(cd(pwd, "..") , "/output/", TYPE, "/", scan, ".csv")

# ╔═╡ d3dd64da-dfa3-4a72-83bd-59690ffa0d73
# CSV.write(output_path, df)

# ╔═╡ a61b162e-26c9-45bf-8376-3c4708cb319f
md"""
### Save full df
"""

# ╔═╡ 4a5dd05c-67ff-493e-b4ec-d86711cb28f6
dfs = []

# ╔═╡ 76cfc410-c728-4bf5-a2dd-916136b8ba35
push!(dfs, df)

# ╔═╡ c219f04d-efca-4d43-a0c8-52a57b8949fd
# if length(dfs) == 24
# 	global new_df = vcat(dfs[1:24]...)
# 	output_path_new = string(cd(pwd, "..") , "/output/", TYPE, "/", "full.csv")
# 	CSV.write(output_path_new, new_df)
# end

# ╔═╡ Cell order:
# ╠═742070ad-ad57-4153-abbb-7ef7d65c6ac3
# ╠═2bf6a8ad-8203-4543-86bf-668dc437e25f
# ╟─4b417377-dfa2-4bdc-b64b-c0543f487c24
# ╠═ac86dfe5-30c7-4d22-adae-152763354408
# ╟─9c998e5a-0f36-4b61-ad70-069ed600808c
# ╟─9bdedef5-2b33-46b7-bed1-f47637b08eef
# ╟─6f79a0d3-7a22-4047-a679-b215713f4bfa
# ╟─c43e6414-e329-4a27-9ba4-9a9f3b1d94ae
# ╟─ff73c57d-f1cf-4f7e-80d1-b9b0e1a0a90f
# ╟─7cb015df-3b20-4168-8a68-0956dc9256ff
# ╟─50ab516f-f988-45c8-a780-30b678b7cfe6
# ╟─74ddfe10-c39e-45bc-b1e7-3335c293c123
# ╟─4812bf96-7699-4dc7-8c7d-b408705e7f31
# ╟─f9b901e4-6fd6-453b-9eda-a31a1601c196
# ╟─d439d790-8c2f-40dd-bef7-72c1baf0848e
# ╟─4e5c0dd3-f9fe-461d-b526-8314a048c1f0
# ╠═b11ec843-d599-4c40-a258-28d87699d7c5
# ╟─05543739-3299-4556-b9eb-ba7366cc5c1d
# ╟─42ad7c3c-521c-4d91-a236-2b56dd3ab674
# ╟─e40c8ac9-60e7-4caa-9ed4-3647c302b47d
# ╟─ccfb5f8c-0d47-4d95-b84e-2afa9fd0e14a
# ╟─764351f3-841e-4d88-b537-f08b711f5f87
# ╟─a00d562f-a0aa-45c9-aaf0-7cedc1c02953
# ╟─4f0cd797-150a-40fe-ae70-d631ef68c627
# ╠═f05ae4e1-a6e7-4351-bf15-426cab19cd10
# ╟─4e6aaf97-58d9-4d8c-9e3f-1cb39acfb06f
# ╟─98774ca0-e8a1-4d4d-b333-fa9236aa8f96
# ╟─0086df28-73ae-4b69-ac2a-48fb946ae3b4
# ╠═6a8b205c-a41f-4e6c-aaf2-b27137dd041b
# ╠═ca7dbfa8-bc4d-410e-8276-e9afbdfdbcf6
# ╠═763b44dc-7323-46e0-b80c-a3b939986483
# ╟─87b784cd-1dfa-43e3-a0f1-8965330d2af3
# ╠═e816098b-a584-4301-8fcf-5849a4b6916c
# ╠═78d1716c-af9a-47da-a320-70042116d752
# ╟─4cce344a-693d-457c-8c4b-7c479e72f2ab
# ╠═9eab7433-c807-4d61-8178-a8400a8f7e24
# ╟─583cc301-29f8-48e7-b208-281a2c552745
# ╠═92710a4d-7f1d-4145-b56c-d18836ddb438
# ╟─f4fd148b-9e95-4826-baa7-e736e9b883a0
# ╠═b99c68ae-86e3-40c2-8b7e-847d9f3fc402
# ╟─aa6a1783-ea13-455c-9118-663b8db76820
# ╟─2fd2f313-5651-429c-a6e0-25c5da0968d7
# ╟─8705619d-7acf-4246-9eb7-658cf7f71547
# ╟─3c78c1c0-8fe6-49d3-81a4-84100b461a14
# ╠═2d426f83-904e-4fac-a0a6-37c286a33a7a
# ╠═76c6c43e-cc4d-43aa-b557-4c0c6860a071
# ╟─e04673e6-5e70-4944-8c40-9162df11e7e9
# ╟─012a1ca7-24d3-447f-a0a4-c8550150936a
# ╟─9fd4245a-0098-41af-8ae4-6c7fc502c4e5
# ╠═5a77d2b9-259d-4f95-bfdf-a2f09ed42cbb
# ╠═6f46e044-50b0-4d78-b8ff-895e1b575331
# ╟─2363391c-0cef-4744-a182-5a1553d25d06
# ╟─8ec96204-8a2a-4143-b999-6a0a5efc7199
# ╟─837da46d-c1e5-43cf-ab7b-642a17a2f764
# ╠═bbf4b093-455d-4097-9d7a-0aeaa21cdeca
# ╠═0c37fdc1-25f9-4c33-a7d1-706210b8ef95
# ╠═6651c12d-2c0f-4699-906f-1d4089422bf0
# ╟─400ffd0c-084b-470e-807f-1f8ad045d631
# ╟─9c1401a2-5f58-41f3-a869-c17ebbf98a64
# ╟─f8781db1-8898-484f-acd8-bd84c0aaac59
# ╟─2a2dced3-055c-44a7-a620-f3e7a8522c68
# ╠═e199a541-c83e-45ba-9650-97d6831aba55
# ╠═a1df38e8-f8bb-47bb-98c5-47be9cd5c8a6
# ╟─48114b90-5c6c-4d24-8d6b-c6ba5e4ae38c
# ╟─317d5fa2-3aac-480b-abda-dbbf824b8b75
# ╟─4d91ed1d-2caf-463c-b80a-893614cc3c44
# ╠═c469b6fa-e11d-424b-81ea-2512f6bc1b9e
# ╠═78cbca0a-78eb-4101-9e63-565300173df6
# ╟─b83ddade-0e5b-447f-801a-dfd81b2d998f
# ╟─c27a5300-0be0-479c-a368-cc81c8991291
# ╟─a963b7eb-d5b9-461d-9f4d-bf94596eb96e
# ╠═a170c6b0-5673-4ec5-a4ff-4c697135f94c
# ╠═c63cf606-c3f4-41ad-b877-0a87b15759d7
# ╟─8beb225c-4285-4e8c-9fdc-14d8e93011c2
# ╟─e2e9d93d-0d34-4236-9e3f-74ed77df4f38
# ╟─2b87600b-1c67-4987-a37e-c9860a52fb80
# ╟─8e395a11-5572-4a96-b2f4-5954c92b58a6
# ╠═6a600411-fd67-4cd8-b38f-7a89a5914749
# ╠═54aad26c-51c7-42cb-972a-ed61c37ea8d3
# ╟─660ec03c-8f6b-439d-8600-6f6248ad6f54
# ╟─9d07d8a3-6052-4992-a34b-773eebce7241
# ╟─e19546f2-ea12-47d8-9398-72daa619b8c0
# ╠═b2199735-aa64-4a54-9422-6e94ab59cc38
# ╠═43fb1435-b69f-487d-bcaf-1f8416e3dcd6
# ╟─15e994bc-c4ee-4b89-b3c7-23e5874291b0
# ╟─32c98b0c-e55e-4cf9-a393-bf5ef777799e
# ╟─672b519f-55bd-4456-9d71-8db8861e9710
# ╠═22f4f737-ad6d-492d-86c2-7155ec3f07d5
# ╠═4d37ec3f-ed5a-48df-8dfd-9bd7ad8b5714
# ╟─66d72b92-3975-4d01-8a7e-deaa58720044
# ╠═3bccc7f7-8230-4495-82fb-23113e69a91a
# ╠═6c547b31-6a36-4ffa-a8d8-bbdcc8b9bf83
# ╠═c320a097-5d86-4af7-83a0-4eb8407a314c
# ╠═d39490dc-cc2c-46dd-87ab-150d420cced6
# ╠═413400a1-de67-4f52-85e0-fb8b1fa9f1aa
# ╠═832dde61-330c-41e8-afa4-4586d22def72
# ╠═407e21cf-7cc0-454e-8e55-f97abc3e6c3c
# ╠═f8e648d6-c42c-472d-b771-22bcfad39f87
# ╠═8a1602e8-44fe-45cf-be29-10accca25519
# ╠═5499d8b7-35b2-41c6-8a4e-321654897c24
# ╠═63cb2095-d00e-4e50-8c4f-c7c50f047bdd
# ╠═b7cea572-e416-4ee4-9342-73ded1ad133e
# ╠═6a2d949f-4083-4605-9455-5d7b9974ac4a
# ╠═0e00f908-458f-4cf3-aa79-de46aed9fb3f
# ╠═d3dd64da-dfa3-4a72-83bd-59690ffa0d73
# ╠═a61b162e-26c9-45bf-8376-3c4708cb319f
# ╠═4a5dd05c-67ff-493e-b4ec-d86711cb28f6
# ╠═76cfc410-c728-4bf5-a2dd-916136b8ba35
# ╠═c219f04d-efca-4d43-a0c8-52a57b8949fd

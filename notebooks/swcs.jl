### A Pluto.jl notebook ###
# v0.19.22

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
	using Pkg
	Pkg.activate(".")

    using PlutoUI, Statistics, CSV, DataFrames, GLM, CairoMakie, HypothesisTests, Colors, MLJBase, DICOM, DICOMUtils, PhantomSegmentation, CalciumScoring, ImageMorphology, ImageFiltering, Noise, Distributions, DSP
    using StatsBase: quantile!, rmsd
end

# ╔═╡ 33bd4068-abd8-45d4-babd-7b397310a645
TableOfContents()

# ╔═╡ 92ac2129-98f8-4dad-a446-6d3664aa783e
md"""
## Load DICOMS
"""

# ╔═╡ c02953cc-0020-452e-938b-694c7f7de44e
begin
	VENDORS = ["80", "100", "120", "135"]
	SIZES = ["small", "medium", "large"]
	DENSITIES = ["low", "normal"]

	IMAGES = "images_new"

	VENDOR = VENDORS[4]
	SIZE = SIZES[1]
	DENSITY = DENSITIES[2]
	
    BASE_PATH = joinpath(dirname(pwd()), IMAGES, SIZE, DENSITY)
	root_path = joinpath(BASE_PATH, VENDOR)
	dcm_path_list = dcm_list_builder(root_path)
	pth = dcm_path_list[1]
	scan = basename(pth)
	header, dcm_array, slice_thick_ori1 = dcm_reader(pth)
	kV = parse(Int64, VENDOR)
end

# ╔═╡ bdc233e2-a17f-42cf-a782-0b040c5b870d
md"""
**Everything below should be automatic, just scroll through to visually inspect that things make sense**
"""

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
        container[i, :] = collect(tuple_array[i])
    end
    return container
end

# ╔═╡ 78e52610-fae1-44fa-b931-275434a3ba73
function overlay_mask_bind(mask)
    indices = findall(x -> x == 1, mask)
    indices = Tuple.(indices)
    label_array = collect_tuple(indices)
    zs = unique(label_array[:, 3])
    return PlutoUI.Slider(1:length(zs); default=3, show_value=true)
end

# ╔═╡ e1147d97-3f08-4a2f-ae30-ed3cb82df965
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

# ╔═╡ 9466ff87-9cef-4f4e-905c-d19070dd6890
function show_matrix(A::Matrix, red::Union{BitMatrix,Matrix{Bool}}=zeros(Bool, size(A)))
    base = RGB.(Gray.(A))

    base[red] .= RGB(1.0, 0.1, 0.1)

    # some tricks to show the pixels more clearly:
    s = max(size(A)...)
    if s >= 20
        min_size = 1200
        factor = min(5, min_size ÷ s)

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

# ╔═╡ 228e339f-ddf0-47b3-821e-5f0e506905ac
function dilate_mask_large(mask)
    return dilate(mask)
end

# ╔═╡ 1a56361a-c5e8-4248-9074-516bf7904466
function ring_mask_large(dilated_mask)
    return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask)
end

# ╔═╡ a7f24627-a821-4396-b2aa-9a2472ff6353
function dilate_mask_medium(mask)
    return (mask)
end

# ╔═╡ c5819fb3-1556-404b-a7b1-9088123227f5
function ring_mask_medium(dilated_mask)
    return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask)
end

# ╔═╡ 01b7ed1f-7d9f-474d-8631-f62bb68eee1a
function dilate_mask_small(mask)
    return (mask)
end

# ╔═╡ 2f421532-c045-4af9-be30-f3a12f66f24b
function ring_mask_small(dilated_mask)
    return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask)
end

# ╔═╡ 24d372ec-3d14-46c4-8319-91c23fd60435
function dilate_mask_large_bkg(mask)
    return dilate(dilate(mask))
end

# ╔═╡ ed442d7f-d151-4b31-8369-0ade86aad81b
function dilate_mask_medium_bkg(mask)
    return dilate(mask)
end

# ╔═╡ 6592a519-4066-4a75-a4e9-4524352a6ec9
function dilate_mask_small_bkg(mask)
    return (mask)
end

# ╔═╡ a4ae678d-a7f1-4bfa-8d22-554583683089
md"""
## Segment Heart
"""

# ╔═╡ d8b8d96e-16ac-4412-b049-6287a933730c
masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3) ÷ 2);

# ╔═╡ 7cde12ee-db55-4030-8d1d-8031cd679674
center_insert

# ╔═╡ 0e91e32b-aa3a-4b3b-a4b9-7b8fbaab132b
@bind a PlutoUI.Slider(1:size(masked_array, 3), default=10, show_value=true)

# ╔═╡ 1ec97def-1513-4309-a6e8-19d7ecfddc25
heatmap(masked_array[:, :, a]; colormap=:grays)

# ╔═╡ fdad818f-1aff-44f7-8fd8-d625b7b873dc
begin
    fig = Figure()

    ax = Makie.Axis(fig[1, 1])
    ax.title = "Raw DICOM Array"
    heatmap!(transpose(dcm_array[:, :, 4]); colormap=:grays)
    scatter!(
        center_insert[2]:center_insert[2],
        center_insert[1]:center_insert[1];
        markersize=10,
        color=:red,
    )
    fig
end

# ╔═╡ ce1cb4e1-9515-463b-b5c7-79c7bc295f49
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

# ╔═╡ 7d341908-e452-4975-bed6-d7f1d3e858af
begin
    fig3 = Figure()

    ax3 = Makie.Axis(fig3[1, 1])
    ax3.title = "Masked DICOM Array"
    heatmap!(transpose(masked_array[:, :, 5]); colormap=:grays)
    scatter!(
        center_insert[2]:center_insert[2],
        center_insert[1]:center_insert[1];
        markersize=10,
        color=:red,
    )
    fig3
end

# ╔═╡ 1de1ae0b-5c08-4035-be82-53a9fdf2a7a9
md"""
## Segment Calcium Rod
"""

# ╔═╡ 47081781-ad0c-4ebb-af72-ee6ea6ccf335
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

# ╔═╡ e1a9f2b8-be8d-4119-9a4a-e1c73e74b33e
@bind c PlutoUI.Slider(1:size(calcium_image, 3), default=cal_rod_slice, show_value=true)

# ╔═╡ 3ef8300d-56a6-4d72-9287-5e4eda72ee24
heatmap(transpose(calcium_image[:, :, c]); colormap=:grays)

# ╔═╡ 8351f18e-7fc7-4c3b-8688-cae0a1e6ec3b
array_filtered = abs.(mapwindow(median, calcium_image[:, :, 2], (3, 3)));

# ╔═╡ d2bf1ee1-fe3e-4ea1-b8c4-1218ac81ff91
bool_arr = array_filtered .> 0;

# ╔═╡ 08527140-2604-4685-a294-43ee44d4aa22
bool_arr_erode = ((erode(erode(erode(bool_arr)))));

# ╔═╡ 436ee3b5-a140-4a7a-a50b-1d90eeb04df6
heatmap(transpose(bool_arr); colormap=:grays)

# ╔═╡ 6ad451a3-6693-4355-8346-9469f7d94167
heatmap(transpose(bool_arr_erode); colormap=:grays)

# ╔═╡ 67cf1029-ebea-4b48-9168-b0d0b491a277
c_img = calcium_image[:, :, 1:3];

# ╔═╡ da901979-3472-48c7-81d5-ff36c1cb0366
begin
    mask_cal_3D = Array{Bool}(undef, size(c_img))
    for z in 1:size(c_img, 3)
        mask_cal_3D[:, :, z] = bool_arr_erode
    end
end;

# ╔═╡ b441160a-84fd-49d0-afae-1b63248cd90f
md"""
## Segment Calcium Inserts
"""

# ╔═╡ 7208bc08-acad-44f7-854c-cea5d0190b24
begin
    root_new = joinpath(dirname(pwd()), "julia_arrays", SIZE)
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

# ╔═╡ 22723c3e-0d5d-4676-b33a-11faa8476d2f
md"""
# Background
"""

# ╔═╡ 6bb03789-2fdb-4bee-a75a-bb4eb30d9ba9


# ╔═╡ 65c9df3c-ea01-48fa-b192-52b3f3b76f9b
md"""
# Score Large Inserts
"""

# ╔═╡ 46c2bdd7-f333-41de-9cd5-2ec5d3a21f2c
arr = masked_array[:, :, 4:6];

# ╔═╡ 61e7055c-01ce-4982-8bf7-6a12c0fad6a5
begin
	## Background
	background_mask1 = zeros(size(arr)...)
	background_mask2 = zeros(size(arr)...)
	background_mask3 = zeros(size(arr)...)

	offs = 10
	rnge = 10
	background_mask1[
		(center_insert[1]-rnge+offs):(center_insert[1]+rnge+offs),
		(center_insert[2]-rnge):(center_insert[2]+rnge),
		2,
	] .= 1
	background_mask2[
		(center_insert[1]-rnge-offs):(center_insert[1]+rnge-offs),
		(center_insert[2]-rnge):(center_insert[2]+rnge),
		2,
	] .= 1
	background_mask3[
		(center_insert[1]-rnge):(center_insert[1]+rnge),
		(center_insert[2]-rnge+offs):(center_insert[2]+rnge+offs),
		2,
	] .= 1

	dilated_mask_L_bkg = dilate_mask_small_bkg(Bool.(background_mask1))
	ring_mask_L_bkg = ring_mask_small(dilated_mask_L_bkg)

	dilated_mask_M_bkg = dilate_mask_small_bkg(Bool.(background_mask2))
	ring_mask_M_bkg = ring_mask_small(dilated_mask_M_bkg)

	dilated_mask_S_bkg = dilate_mask_small_bkg(Bool.(background_mask3))
	ring_mask_S_bkg = ring_mask_small(dilated_mask_S_bkg)
end;

# ╔═╡ 9d7a0818-c0de-47f0-9c53-381b0ddee674
let
	f = Figure()
	ax = CairoMakie.Axis(f[1, 1])
	heatmap!(dilated_mask_L_bkg[:, :, 2], colormap=:grays)

	ax = CairoMakie.Axis(f[1, 2])
	heatmap!(dilated_mask_M_bkg[:, :, 2], colormap=:grays)

	ax = CairoMakie.Axis(f[2, 1])
	heatmap!(dilated_mask_S_bkg[:, :, 2], colormap=:grays)
	f
end

# ╔═╡ 00c43b6b-69ff-48ed-8d6d-fc3e44a6a4f2
single_arr = masked_array[:, :, slice_CCI];

# ╔═╡ 1f353ddd-4e9a-4437-8d2e-97a4c466ee28
md"""
## High Density
"""

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

# ╔═╡ 64fcb0ab-7106-4bf2-b268-cc5e814dc942
μ, σ = mean(c_img[Bool.(mask_cal_3D)]) / 2, std(c_img[Bool.(mask_cal_3D)])

# ╔═╡ aa851337-3dd3-4521-a402-52a429bec1e6
@bind g2 overlay_mask_bind(dilated_mask_L_HD)

# ╔═╡ e853a3d8-0a18-43b7-949b-c374203d4838
overlay_mask_plot(arr, dilated_mask_L_HD, g2, "dilated mask")

# ╔═╡ 0373a338-9a93-42b8-98ee-9768a50dda67
pixel_size = DICOMUtils.get_pixel_size(header)

# ╔═╡ 268ad486-9907-4aaa-a271-a77664fea675
overlayed_mask_l_hd = create_mask(arr, dilated_mask_L_HD);

# ╔═╡ 405f2f7e-123a-4ccf-9fc1-9b6eee4dae12
alg = Agatston()

# ╔═╡ 19c04849-68a6-4ea3-8a1b-e9183c9b852e
alg2 = SpatiallyWeighted()

# ╔═╡ c604fc41-00ef-498b-ac8b-b5e5cd438306
begin
	overlayed_mask_l_bkg = create_mask(arr, dilated_mask_L_bkg)
	overlayed_mask_m_bkg = create_mask(arr, dilated_mask_M_bkg)
	overlayed_mask_s_bkg = create_mask(arr, dilated_mask_S_bkg)
	
	swcs_bkg_large = score(overlayed_mask_l_bkg, μ, σ, alg2)
	swcs_bkg_medium = score(overlayed_mask_m_bkg, μ, σ, alg2)
	swcs_bkg_small = score(overlayed_mask_s_bkg, μ, σ, alg2)

	swcs_bkg = [swcs_bkg_large, swcs_bkg_medium, swcs_bkg_small]
end;

# ╔═╡ 38437b4c-e4e4-478d-be0b-3f19d1e4db97
swcs_bkg[1]

# ╔═╡ 310822fe-5664-4177-9033-6945717fa5b6
swcs_bkg[2]

# ╔═╡ 3007cc2f-9073-4ba1-8b93-b6005d1daba6
swcs_bkg[3]

# ╔═╡ 9354aeaa-053b-4601-b83c-0c4270d6dc5d
agat_l_hd, vol_l_hd = score(overlayed_mask_l_hd, pixel_size, alg)

# ╔═╡ 05cf3e18-950e-427c-91d9-be289da5c5c5
swcs_l_hd = score(overlayed_mask_l_hd, μ, σ, alg2)

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

# ╔═╡ 69e06cfd-0cd2-485b-9338-014f40f5d4f9
agat_l_md, vol_l_md = score(overlayed_mask_l_md, pixel_size, alg)

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

# ╔═╡ 755379df-d43a-42a3-8d2f-56fdd2d85f7e
agat_l_ld, vol_l_ld = score(overlayed_mask_l_ld, pixel_size, alg)

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
agat_m_hd, vol_m_hd = score(overlayed_mask_m_hd, pixel_size, alg)

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
agat_m_md, vol_m_md = score(overlayed_mask_m_md, pixel_size, alg)

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
agat_m_ld, vol_m_ld = score(overlayed_mask_m_ld, pixel_size, alg)

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
agat_s_hd, vol_s_hd = score(overlayed_mask_s_hd, pixel_size, alg)

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
agat_s_md, vol_s_md = score(overlayed_mask_s_md, pixel_size, alg)

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
agat_s_ld, vol_s_ld = score(overlayed_mask_s_ld, pixel_size, alg)

# ╔═╡ c2116577-2f05-40dc-9094-92c621955f93
swcs_s_ld = score(overlayed_mask_s_ld, μ, σ, alg2)

# ╔═╡ 258caf22-6f46-4341-80eb-d572aa26ca59
md"""
# Results
"""

# ╔═╡ 3d573812-8bc2-4dce-99fb-6883b86433f4
density_array = [0, 200, 400, 800]

# ╔═╡ 52f3d56d-9219-43f3-ba32-925faa2e0e07
inserts = ["Low Density", "Medium Density", "High Density"]

# ╔═╡ 260b441d-12eb-44d5-adb2-dc5401abdb20
md"""
## Agatston
"""

# ╔═╡ c0b060b6-ef4d-415e-9f8c-a39ab4373541
calculated_agat_large = [agat_l_ld, agat_l_md, agat_l_hd]

# ╔═╡ 3f79ca52-b120-43af-958e-6a0e5edf079f
calculated_agat_medium = [agat_m_ld, agat_m_md, agat_m_hd]

# ╔═╡ 0de07ba8-f1d6-4c01-86fd-82e42e9a29bb
calculated_agat_small = [agat_s_ld, agat_s_md, agat_s_hd]

# ╔═╡ 40f4a8e5-8a97-43d6-8d83-4993d5f59b83
md"""
## Spatially Weighted
"""

# ╔═╡ e89bd455-71ba-41c0-870c-2dad225e129e
calculated_swcs_large = [swcs_l_ld, swcs_l_md, swcs_l_hd]

# ╔═╡ 36130f3d-8a03-49ce-ab36-f9ff0b452c8d
calculated_swcs_medium = [swcs_m_ld, swcs_m_md, swcs_m_hd]

# ╔═╡ e53a38e7-ba34-45af-bc39-86026763241f
calculated_swcs_small = [swcs_s_ld, swcs_s_md, swcs_s_hd]

# ╔═╡ Cell order:
# ╠═98114689-d1cc-4aa9-b88f-c048059e95d3
# ╠═33bd4068-abd8-45d4-babd-7b397310a645
# ╟─92ac2129-98f8-4dad-a446-6d3664aa783e
# ╠═c02953cc-0020-452e-938b-694c7f7de44e
# ╟─bdc233e2-a17f-42cf-a782-0b040c5b870d
# ╟─b6f33cae-29c6-4d13-9b0c-a2551bd1e08e
# ╟─818788f1-fb7b-4adc-9111-1f5cca892f86
# ╟─78e52610-fae1-44fa-b931-275434a3ba73
# ╟─e1147d97-3f08-4a2f-ae30-ed3cb82df965
# ╟─9466ff87-9cef-4f4e-905c-d19070dd6890
# ╟─1c2fa6f8-201d-40ae-8c03-7e5c5893fc1e
# ╟─228e339f-ddf0-47b3-821e-5f0e506905ac
# ╟─1a56361a-c5e8-4248-9074-516bf7904466
# ╟─a7f24627-a821-4396-b2aa-9a2472ff6353
# ╟─c5819fb3-1556-404b-a7b1-9088123227f5
# ╟─01b7ed1f-7d9f-474d-8631-f62bb68eee1a
# ╟─2f421532-c045-4af9-be30-f3a12f66f24b
# ╟─24d372ec-3d14-46c4-8319-91c23fd60435
# ╟─ed442d7f-d151-4b31-8369-0ade86aad81b
# ╟─6592a519-4066-4a75-a4e9-4524352a6ec9
# ╟─a4ae678d-a7f1-4bfa-8d22-554583683089
# ╠═7cde12ee-db55-4030-8d1d-8031cd679674
# ╠═d8b8d96e-16ac-4412-b049-6287a933730c
# ╟─0e91e32b-aa3a-4b3b-a4b9-7b8fbaab132b
# ╠═1ec97def-1513-4309-a6e8-19d7ecfddc25
# ╟─fdad818f-1aff-44f7-8fd8-d625b7b873dc
# ╟─ce1cb4e1-9515-463b-b5c7-79c7bc295f49
# ╟─7d341908-e452-4975-bed6-d7f1d3e858af
# ╟─1de1ae0b-5c08-4035-be82-53a9fdf2a7a9
# ╠═47081781-ad0c-4ebb-af72-ee6ea6ccf335
# ╟─e1a9f2b8-be8d-4119-9a4a-e1c73e74b33e
# ╠═3ef8300d-56a6-4d72-9287-5e4eda72ee24
# ╠═8351f18e-7fc7-4c3b-8688-cae0a1e6ec3b
# ╠═d2bf1ee1-fe3e-4ea1-b8c4-1218ac81ff91
# ╠═08527140-2604-4685-a294-43ee44d4aa22
# ╠═436ee3b5-a140-4a7a-a50b-1d90eeb04df6
# ╠═6ad451a3-6693-4355-8346-9469f7d94167
# ╠═67cf1029-ebea-4b48-9168-b0d0b491a277
# ╠═da901979-3472-48c7-81d5-ff36c1cb0366
# ╟─b441160a-84fd-49d0-afae-1b63248cd90f
# ╠═7208bc08-acad-44f7-854c-cea5d0190b24
# ╟─22723c3e-0d5d-4676-b33a-11faa8476d2f
# ╠═61e7055c-01ce-4982-8bf7-6a12c0fad6a5
# ╟─9d7a0818-c0de-47f0-9c53-381b0ddee674
# ╠═c604fc41-00ef-498b-ac8b-b5e5cd438306
# ╠═38437b4c-e4e4-478d-be0b-3f19d1e4db97
# ╠═310822fe-5664-4177-9033-6945717fa5b6
# ╠═3007cc2f-9073-4ba1-8b93-b6005d1daba6
# ╠═6bb03789-2fdb-4bee-a75a-bb4eb30d9ba9
# ╟─65c9df3c-ea01-48fa-b192-52b3f3b76f9b
# ╠═46c2bdd7-f333-41de-9cd5-2ec5d3a21f2c
# ╠═00c43b6b-69ff-48ed-8d6d-fc3e44a6a4f2
# ╟─1f353ddd-4e9a-4437-8d2e-97a4c466ee28
# ╠═a228cbf8-4d79-41a4-b419-afe26b981a1c
# ╟─334e9473-fc4f-466c-a046-6c66bf65b738
# ╠═a2c59aa5-6c86-4d4a-aa84-c6675446083d
# ╠═64fcb0ab-7106-4bf2-b268-cc5e814dc942
# ╟─aa851337-3dd3-4521-a402-52a429bec1e6
# ╠═e853a3d8-0a18-43b7-949b-c374203d4838
# ╠═0373a338-9a93-42b8-98ee-9768a50dda67
# ╠═268ad486-9907-4aaa-a271-a77664fea675
# ╠═405f2f7e-123a-4ccf-9fc1-9b6eee4dae12
# ╠═19c04849-68a6-4ea3-8a1b-e9183c9b852e
# ╠═9354aeaa-053b-4601-b83c-0c4270d6dc5d
# ╠═05cf3e18-950e-427c-91d9-be289da5c5c5
# ╟─5a09303d-5360-48c9-9554-c769e63b1ee0
# ╠═20a53747-25a1-4642-9043-64109ff28711
# ╟─5f0f4233-ae8b-432f-93e9-f2b3a54dc890
# ╠═bd5aa941-8ab9-4b43-b6f5-f2320f1f9e03
# ╟─1d87f244-f6c5-4e12-9d0d-5c2e7c4805e5
# ╠═8f9e3da4-baab-431d-a4c7-67e753144aef
# ╠═31708cf8-1d56-4b16-ad4e-df356a8924fc
# ╠═69e06cfd-0cd2-485b-9338-014f40f5d4f9
# ╠═2269beda-ae30-49ca-af6a-65a681dc7405
# ╟─9d8e3d8d-a4e6-418c-a241-b4b1cd119a14
# ╠═8f125f2e-ec57-40f0-96d6-e04a544db0c3
# ╟─ad1c418f-a4b7-45fa-8798-534b5e474b99
# ╠═f72f5020-6e3f-4b47-9340-55e5d687b56c
# ╟─d3ce01ea-a5d6-41e3-8f98-e0042026da9b
# ╠═2e85a77c-7bb0-47ba-949e-4fdccea6b8e8
# ╠═db774a1b-dbda-42a3-abb2-89ae26a2a31c
# ╠═755379df-d43a-42a3-8d2f-56fdd2d85f7e
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

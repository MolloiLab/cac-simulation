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

# ╔═╡ a2647e62-e189-498a-981c-9980485fceca
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.activate(".")

    using PlutoUI, Statistics, CSV, DataFrames, GLM, CairoMakie, HypothesisTests, Colors, MLJBase, DICOM, DICOMUtils, PhantomSegmentation, CalciumScoring, ImageMorphology, ImageFiltering, Noise, NPZ
    using StatsBase: quantile!, rmsd
end

# ╔═╡ f4e4decf-d8bf-4203-b44a-0d2cb9830a3a
TableOfContents()

# ╔═╡ d161adaf-c02a-49d1-a5dd-6929ea0614c8
pwd()

# ╔═╡ bb5f1496-9b8f-4d26-8f8d-03773fa546a9
md"""
## Load Arrays
"""

# ╔═╡ f1777984-c7e3-4fd2-92ba-7649d32caab8
begin
	sizes_arr = ["large", "medium", "small"]
	densities_arr = ["low", "normal"]
	energies_arr = ["80", "100", "120", "135"]
	energies_motion = ["80-motion.npy", "100-motion-3.npy", "120-motion-3.npy", "135-motion-3.npy"]
end

# ╔═╡ 7b6e134b-2249-479d-a3f6-f0a2791a25f7
begin
	SIZE = sizes_arr[1]
	DENSITY = densities_arr[2]
	ENERGY = energies_arr[1]
	ENERGY_MOTION = energies_motion[1]
end

# ╔═╡ 44a31da5-9914-44aa-b3ef-551ab3ea30ae
md"""
**Everything below should be automatic, just scroll through to visually inspect that things make sense**
"""

# ╔═╡ 3de3c36e-6350-4aab-a139-92d052641136
np_path = joinpath(dirname(pwd()), "images_new", SIZE, DENSITY, ENERGY_MOTION)

# ╔═╡ 9155bc22-344b-4ff1-b629-ce7b935a3370
dcm_path = joinpath(dirname(pwd()), "images_new", SIZE, DENSITY, ENERGY)

# ╔═╡ 88ab1259-baa7-422a-9321-0f8f3d40b051
header, _, _ = dcm_reader(dcm_path);

# ╔═╡ f45b6ec1-330d-4e6e-a36f-f10002489f25
begin
	dcm = npzread(np_path);
	dcm_array = permutedims(dcm, (2, 1, 3))
end;

# ╔═╡ 38c4b7f8-1c7f-4776-922b-2a74dd4ffb3d
md"""
## Helper Functions
"""

# ╔═╡ 9f43b7c7-e99c-4893-aa4b-2252dc3ba91a
function collect_tuple(tuple_array)
    row_num = size(tuple_array)
    col_num = length(tuple_array[1])
    container = zeros(Int64, row_num..., col_num)
    for i in 1:length(tuple_array)
        container[i, :] = collect(tuple_array[i])
    end
    return container
end

# ╔═╡ a352f0bf-f2a3-420f-ad78-0e3a8fab20ae
function overlay_mask_bind(mask)
    indices = findall(x -> x == 1, mask)
    indices = Tuple.(indices)
    label_array = collect_tuple(indices)
    zs = unique(label_array[:, 3])
    return PlutoUI.Slider(1:length(zs); default=3, show_value=true)
end

# ╔═╡ 003e03c1-3796-4f58-853c-055bbc2dfb37
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

# ╔═╡ 66dc4873-d94c-4d6e-b08f-5dc31ae3f19b
function dilate_mask_large(mask)
	return dilate(mask)
end

# ╔═╡ 64717903-e79e-4e4d-9c16-f780770b786a
function ring_mask_large(dilated_mask)
	return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask);
end

# ╔═╡ 64eca4c3-d079-4c4e-9070-f39095d97163
md"""
## Segment Heart
"""

# ╔═╡ 12205cce-c9cf-47df-b929-92eb4f72309b
masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3) ÷ 2);

# ╔═╡ 512ae1c7-0184-4fe4-a492-899e407711ee
@bind a PlutoUI.Slider(1:size(masked_array, 3), default=10, show_value=true)

# ╔═╡ a159ae46-dfce-4a0a-80ab-d1695a9acefd
heatmap(transpose(dcm_array[:, :, a]); colormap=:grays)

# ╔═╡ b097a0aa-f5e3-48e5-ad32-275723015e1c
heatmap(transpose(masked_array[:, :, a]); colormap=:grays)

# ╔═╡ d0bc44a9-cf63-419f-b291-a6f40fd0d8fb
xs, yx = 200:250, 375:425

# ╔═╡ 7012f3c3-4a49-4ff0-8c9a-248b5e9fb8a3
heatmap(transpose(masked_array[xs, yx, a]); colormap=:grays)

# ╔═╡ 1d20bd27-1c22-403b-b02b-bc30f282e819
# let
#     fig = Figure()

#     ax = Makie.Axis(fig[1, 1])
#     heatmap!(transpose(dcm_array[:, :, 5]); colormap=:grays)

# 	save("/Users/daleblack/Google Drive/Research/Papers/My Papers/cac-simulation/figures-review/simulated_phantom.png", fig)
#     fig
# end

# ╔═╡ b7e9d017-da52-479a-b804-7cd74fd1af5e
# let
#     fig = Figure()

#     ax = Makie.Axis(fig[1, 1])
#     ax.title = "Raw DICOM Array"
#     heatmap!(transpose(dcm_array[:, :, 5]); colormap=:grays)
#     scatter!(
#         center_insert[2]:center_insert[2],
#         center_insert[1]:center_insert[1];
#         markersize=10,
#         color=:red,
#     )
#     fig
# end

# ╔═╡ 6722dbc9-0d9a-4570-89e2-d2fba04b7b5b
# let
#     fig2 = Figure()

#     ax2 = Makie.Axis(fig2[1, 1])
#     ax2.title = "Mask Array"
#     heatmap!(transpose(mask); colormap=:grays)
#     scatter!(
#         center_insert[2]:center_insert[2],
#         center_insert[1]:center_insert[1];
#         markersize=10,
#         color=:red,
#     )
#     fig2
# end

# ╔═╡ 8e4c21d5-50d0-49a1-9342-8cf027c40d6b
# let
#     fig3 = Figure()

#     ax3 = Makie.Axis(fig3[1, 1])
#     ax3.title = "Masked DICOM Array"
#     heatmap!(transpose(masked_array[:, :, 1]); colormap=:grays)
#     scatter!(
#         center_insert[2]:center_insert[2],
#         center_insert[1]:center_insert[1];
#         markersize=10,
#         color=:red,
#     )
#     fig3
# end

# ╔═╡ b2a82830-7240-49d4-9503-7a817f138d28
md"""
## Segment Calcium Rod
"""

# ╔═╡ 370115da-6071-4abe-afce-314bca318ed6
if DENSITY == "low"
	thresh = 55
elseif DENSITY == "normal"
	thresh = 130
end

# ╔═╡ 2e2953cb-cea8-4cf9-afdd-872ead7da567
calcium_image, slice_CCI, quality_slice, cal_rod_slice = mask_rod(
    masked_array, header; calcium_threshold=thresh
);

# ╔═╡ 0b88e378-cda1-4400-8005-51e25fcef33e
array_filtered = abs.(mapwindow(median, calcium_image[:, :, 2], (3, 3)));

# ╔═╡ 33408dd8-e97a-45c6-a99e-1709d30c924c
c_img = calcium_image[:, :, 1:3];

# ╔═╡ 24e19f9a-c820-4cdc-b3bb-cc081d1ca0e7
@bind c PlutoUI.Slider(1:size(calcium_image, 3), default=5, show_value=true)

# ╔═╡ d0941946-4338-4b8b-b4cc-e83f8780a253
heatmap(transpose(calcium_image[:, :, c]); colormap=:grays)

# ╔═╡ 4db2de22-7cf4-4f45-a895-76547816e1dd
md"""
## Load (segment) Calcium Inserts
"""

# ╔═╡ 9ae8c021-f465-4e55-ae5a-43cfb08a7e01
begin
    root_new = string(
        "/Users/daleblack/Google Drive/dev/MolloiLab/cac-simulation/julia_arrays/",
        SIZE,
        "/",
    )
    mask_L_HD = Array(CSV.read(string(root_new, "mask_L_HD.csv"), DataFrame; header=false))
    mask_M_HD = Array(CSV.read(string(root_new, "mask_M_HD.csv"), DataFrame; header=false))
    mask_S_HD = Array(CSV.read(string(root_new, "mask_S_HD.csv"), DataFrame; header=false))
    mask_L_MD = Array(CSV.read(string(root_new, "mask_L_MD.csv"), DataFrame; header=false))
    mask_M_MD = Array(CSV.read(string(root_new, "mask_M_MD.csv"), DataFrame; header=false))
    mask_S_MD = Array(CSV.read(string(root_new, "mask_S_MD.csv"), DataFrame; header=false))
    mask_L_LD = Array(CSV.read(string(root_new, "mask_L_LD.csv"), DataFrame; header=false))
    mask_M_LD = Array(CSV.read(string(root_new, "mask_M_LD.csv"), DataFrame; header=false))
    mask_S_LD = Array(CSV.read(string(root_new, "mask_S_LD.csv"), DataFrame; header=false))
end;

# ╔═╡ e3dbd551-e3f5-461b-931e-90c3de086270
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

# ╔═╡ 4a8fe98a-d2d1-4fe2-a4d4-5a1898c49f84
# heatmap(transpose(masks); colormap=:grays)

# ╔═╡ 0cbc4ec2-790f-45a0-801d-f04c25657bfb
md"""
## Calibration Prep
"""

# ╔═╡ 04e06029-e368-4ee4-8394-ad4333e68232
bool_arr = array_filtered .> 0;

# ╔═╡ 02010cf7-7bf8-44dc-a01b-8f5d57df5733
bool_arr_erode = (erode(erode(erode(erode(bool_arr)))));

# ╔═╡ d20b5d55-e8fb-43cd-b447-4b6cf558db6d
# heatmap(bool_arr; colormap=:grays)

# ╔═╡ 84419df2-9512-4dcb-b4af-e5d194fd65db
# heatmap(bool_arr_erode; colormap=:grays)

# ╔═╡ 60e001ba-cabd-4df8-bb48-827be2332513
begin
    mask_cal_3D = Array{Bool}(undef, size(c_img))
    for z in 1:size(c_img, 3)
        mask_cal_3D[:, :, z] = bool_arr_erode
    end
end;

# ╔═╡ d274581a-ae43-4398-b5ee-d85b980b56cc
hu_calcium = mean(c_img[mask_cal_3D])

# ╔═╡ 16e0bc1a-7b56-4f91-8b72-4677b782d3fd
# hist(c_img[mask_cal_3D])

# ╔═╡ 07d2a022-e61d-4737-b586-9816bcbd72d1
ρ_calcium = 0.2

# ╔═╡ 9995fadd-7291-46e4-be27-b7e1e8dbce0e
# heatmap(single_arr)

# ╔═╡ 2ed40040-bfec-4ad4-b87f-29caccc88712
md"""
# Score Large Inserts
"""

# ╔═╡ 78587458-4711-4c2b-8755-ed4030a0a18e
arr = masked_array[:, :, 4:6];

# ╔═╡ aef27e7c-8164-469b-b5ed-712bdf2cd645
single_arr = masked_array[:, :, slice_CCI];

# ╔═╡ e2f0e863-8f4b-46ab-af9c-9d59538a9eec
hu_heart_tissue = mean(single_arr[center_insert[1]-10:center_insert[1]+10, center_insert[2]-10:center_insert[2]+10])

# ╔═╡ 27eb5242-330a-4f7d-aad9-9067a3e869ea
begin
	pixel_size = DICOMUtils.get_pixel_size(header)
	voxel_size = pixel_size[1] * pixel_size[2] * pixel_size[3]
end

# ╔═╡ bab453c4-0057-42e3-a009-02c5490f847f
md"""
## High Density
"""

# ╔═╡ 7ede2237-aaec-4fc3-86aa-bade846dad15
begin
    mask_L_HD_3D = Array{Bool}(undef, size(arr))
    for z in 1:size(arr, 3)
        mask_L_HD_3D[:, :, z] = mask_L_HD
    end
end;

# ╔═╡ e8cbc4c1-bc45-4065-afb0-014eaf4b3638
md"""
#### Dilated mask
"""

# ╔═╡ f1e04fa6-f5bd-42c9-a5d0-93960c3e6b92
begin
	dilated_mask_L_HD = dilate_mask_large(mask_L_HD_3D)
	ring_mask_L_HD = ring_mask_large(dilated_mask_L_HD)
end;

# ╔═╡ d4f4cd88-38f1-4801-8973-32cdc6e172b6
@bind g2 overlay_mask_bind(dilated_mask_L_HD)

# ╔═╡ 2d0d5ff5-b2f3-4d44-be7b-50772222f277
overlay_mask_plot(arr, dilated_mask_L_HD, g2, "dilated mask")

# ╔═╡ 28324392-eecf-4865-a788-a871decc255d
md"""
#### Ring (background) mask
"""

# ╔═╡ 1608d803-cf85-4f99-84d4-0b89b90dea34
@bind g4 overlay_mask_bind(ring_mask_L_HD)

# ╔═╡ 9594574b-87aa-49da-822c-1d9a835f21c1
overlay_mask_plot(arr, ring_mask_L_HD, g4, "ring mask")

# ╔═╡ 52c0720b-6025-477d-9bb9-eba8e9457bb3
begin
    single_ring_mask_L_HD = Bool.(ring_mask_L_HD[:, :, 3])
    hu_heart_tissue_large_hd = mean(arr[ring_mask_L_HD])
end

# ╔═╡ f0eb66fd-95ab-466d-b201-67308e6f53db
mass_large_hd = score(arr[dilated_mask_L_HD], hu_calcium, hu_heart_tissue_large_hd, voxel_size, ρ_calcium, VolumeFraction())

# ╔═╡ 5377fc95-b66a-42f6-b1a7-5d7c83a0714a
md"""
## Medium Density
"""

# ╔═╡ fb6f86d5-531c-4a28-b225-5d4c1b134d82
begin
    mask_L_MD_3D = Array{Bool}(undef, size(arr))
    for z in 1:size(arr, 3)
        mask_L_MD_3D[:, :, z] = mask_L_MD
    end
end;

# ╔═╡ fa0acc1f-c333-4ead-bc03-0a748466d79d
md"""
#### Dilated mask
"""

# ╔═╡ f4f4c29e-6838-405b-b65b-1e0bed3d9854
begin
	dilated_mask_L_MD = dilate_mask_large(mask_L_MD_3D)
	ring_mask_L_MD = ring_mask_large(dilated_mask_L_MD)
end;

# ╔═╡ 58a299a7-e138-4c76-a016-a27cf71697a7
# dilated_mask_L_MD = dilate(dilate(mask_L_MD_3D));

# ╔═╡ 525199de-753c-44d6-967d-f5d3cbc30fef
# ring_mask_L_MD =
#     Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask_L_MD)))))) - dilated_mask_L_MD);

# ╔═╡ 28a1deed-248e-4b15-896c-ff75dba4ddd2
@bind h2 overlay_mask_bind(dilated_mask_L_MD)

# ╔═╡ b91850d3-aa07-461d-ae78-a805b0d4a23e
overlay_mask_plot(arr, dilated_mask_L_MD, h2, "dilated mask")

# ╔═╡ 6b8af2c6-fd50-44f2-b1b8-93d083e1bb21
md"""
#### Ring (background) mask
"""

# ╔═╡ 6adf7bde-eaba-4480-8229-7d56fd765f68
@bind h4 overlay_mask_bind(ring_mask_L_MD)

# ╔═╡ eadf72ff-6d0c-4e5f-a32e-dc42cb83f93e
overlay_mask_plot(arr, ring_mask_L_MD, h4, "ring mask")

# ╔═╡ 2c966c01-56db-4661-a0e5-276ba44b02e8
begin
    single_ring_mask_L_MD = Bool.(ring_mask_L_MD[:, :, 3])
    hu_heart_tissue_large_md = mean(arr[ring_mask_L_MD])
end

# ╔═╡ 4ab5580a-d5f2-42df-ade5-53b0cb0c9293
mass_large_md = score(arr[dilated_mask_L_MD], hu_calcium, hu_heart_tissue_large_md, voxel_size, ρ_calcium, VolumeFraction())

# ╔═╡ e32e2ddb-cc73-4474-8711-d7fd105b359d
md"""
## Low Density
"""

# ╔═╡ b5fd1d9c-de0a-4e74-8808-ce57c777c675
begin
    mask_L_LD_3D = Array{Bool}(undef, size(arr))
    for z in 1:size(arr, 3)
        mask_L_LD_3D[:, :, z] = mask_L_LD
    end
end;

# ╔═╡ 14d00a87-f964-450e-b631-7a36c40d6834
md"""
#### Dilated mask
"""

# ╔═╡ b209a40f-2434-4940-973e-e6e3debbbbe4
# dilated_mask_L_LD = dilate(dilate(mask_L_LD_3D));

# ╔═╡ ceac8866-7932-4bf3-a2ee-eb0393b4a036
# ring_mask_L_LD =
#     Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask_L_LD)))))) - dilated_mask_L_LD);

# ╔═╡ ebe4a819-0466-40fc-b29d-94fe2dad15a3
begin
	dilated_mask_L_LD = dilate_mask_large(mask_L_LD_3D)
	ring_mask_L_LD = ring_mask_large(dilated_mask_L_LD)
end;

# ╔═╡ 257efcfc-0c72-4b68-9a32-8b04b39da835
@bind j2 overlay_mask_bind(dilated_mask_L_LD)

# ╔═╡ 69ff66e2-afa3-44df-989b-7497e390b252
overlay_mask_plot(arr, dilated_mask_L_LD, j2, "dilated mask")

# ╔═╡ 69513ba5-9e1a-4932-b056-eac183aa3ec5
md"""
#### Ring (background) mask
"""

# ╔═╡ 2d9fb8ed-3b2a-4e6c-880b-9d2d1627bf60
@bind j4 overlay_mask_bind(ring_mask_L_LD)

# ╔═╡ 7fa3423e-6ad0-470f-b1b3-9e62f9533ae6
overlay_mask_plot(arr, ring_mask_L_LD, j4, "ring mask")

# ╔═╡ 260258de-360f-410f-b43c-bdfb72be5322
begin
    single_ring_mask_L_LD = Bool.(ring_mask_L_LD[:, :, 3])
    hu_heart_tissue_large_ld = mean(arr[ring_mask_L_LD])
end

# ╔═╡ 2d43321b-2836-4887-a4d6-0e6a99fbc0f2
mass_large_ld = score(arr[dilated_mask_L_LD], hu_calcium, hu_heart_tissue_large_ld, voxel_size, ρ_calcium, VolumeFraction())

# ╔═╡ 5f475234-cd14-47dc-b5a2-f4c8b0917849
md"""
# Score Medium Inserts
"""

# ╔═╡ 38b812bd-7291-4dc6-9838-9e0e65e59a0f
function dilate_mask_medium(mask)
	return (mask)
end

# ╔═╡ 300e41a3-bd66-4b1f-a225-e1b2b8cf39d7
function ring_mask_medium(dilated_mask)
	return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask);
end

# ╔═╡ 046a03f1-a285-4a80-bd64-e821252d07f0
md"""
## High Density
"""

# ╔═╡ 9a201bfe-beb9-4729-b8af-c4186bcd0d3b
begin
    mask_M_HD_3D = Array{Bool}(undef, size(arr))
    for z in 1:size(arr, 3)
        mask_M_HD_3D[:, :, z] = mask_M_HD
    end
end;

# ╔═╡ c76ad7d6-76a1-42b8-b5ef-3ebd342177a9
md"""
#### Dilated mask
"""

# ╔═╡ 57c2b58f-6bef-4088-a734-9fc21ddb6217
# dilated_mask_M_HD = dilate(dilate(mask_M_HD_3D));

# ╔═╡ 6e833dc6-412f-4273-b658-0a1e27850d55
# ring_mask_M_HD =
#     Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask_M_HD)))))) - dilated_mask_M_HD);

# ╔═╡ 34b427ef-d767-4216-b570-928ef05c7258
begin
	dilated_mask_M_HD = dilate_mask_medium(mask_M_HD_3D)
	ring_mask_M_HD = ring_mask_medium(dilated_mask_M_HD)
end;

# ╔═╡ b992a000-8138-4c70-9b38-342bfb4460a7
@bind z2 overlay_mask_bind(dilated_mask_M_HD)

# ╔═╡ 5452a3a7-ab92-4117-b357-cc7c806b6a83
overlay_mask_plot(arr, dilated_mask_M_HD, z2, "dilated mask")

# ╔═╡ adab3fe3-8ff7-49af-9aa0-c0e6e631fed3
md"""
#### Ring (background) mask
"""

# ╔═╡ ddbf8772-16e0-4ad1-a950-779a5d1c8a44
@bind z4 overlay_mask_bind(ring_mask_M_HD)

# ╔═╡ d3b28cdf-909d-45ad-bcd2-46524bc527f6
overlay_mask_plot(arr, ring_mask_M_HD, z4, "ring mask")

# ╔═╡ fb0a7476-684c-4c6c-a522-5bf1ea5c5247
begin
    single_ring_mask_M_HD = Bool.(ring_mask_M_HD[:, :, 3])
    hu_heart_tissue_medium_hd = mean(arr[ring_mask_M_HD])
end

# ╔═╡ e51f9107-94a5-490c-a03b-14ea5aadbb47
mass_medium_hd = score(arr[dilated_mask_M_HD], hu_calcium, hu_heart_tissue_medium_hd, voxel_size, ρ_calcium, VolumeFraction())

# ╔═╡ dcf534a7-8545-4a9f-87a8-5ca031541f43
md"""
## Medium Density
"""

# ╔═╡ 08ebf130-56fa-4b5f-8cc7-3d25e2326a00
begin
    mask_M_MD_3D = Array{Bool}(undef, size(arr))
    for z in 1:size(arr, 3)
        mask_M_MD_3D[:, :, z] = mask_M_MD
    end
end;

# ╔═╡ dcf7e7b4-0302-48f9-818f-f059c2873b2f
md"""
#### Dilated mask
"""

# ╔═╡ 60c3a5ac-0a81-4352-b2e1-3ea89513b0ee
# dilated_mask_M_MD = dilate(dilate(mask_M_MD_3D));

# ╔═╡ c77d924d-4496-4734-9e08-bb19d8a01a33
# ring_mask_M_MD =
#     Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask_M_MD)))))) - dilated_mask_M_MD);

# ╔═╡ 2bef5e72-b48d-41e8-b6c3-846cefc2281d
begin
	dilated_mask_M_MD = dilate_mask_medium(mask_M_MD_3D)
	ring_mask_M_MD = ring_mask_medium(dilated_mask_M_MD)
end;

# ╔═╡ 18cbce65-2543-4210-896a-b10e80ba7d95
@bind x2 overlay_mask_bind(dilated_mask_M_MD)

# ╔═╡ 56ede52c-706f-41f7-bae1-059ab0b3bb82
overlay_mask_plot(arr, dilated_mask_M_MD, x2, "dilated mask")

# ╔═╡ 6b8797f1-4d2d-4b64-ae3e-fb436b2ab8a2
md"""
#### Ring (background) mask
"""

# ╔═╡ f0d2a0b9-74be-4ca7-8efd-208a5010cac2
@bind x4 overlay_mask_bind(ring_mask_M_MD)

# ╔═╡ 4c3db4e5-eacd-4e5e-a6e8-a35504d9c01a
overlay_mask_plot(arr, ring_mask_M_MD, x4, "ring mask")

# ╔═╡ 2f02bdcb-1d2e-4313-b338-2cac8143662b
begin
    single_ring_mask_M_MD = Bool.(ring_mask_M_MD[:, :, 3])
    hu_heart_tissue_medium_md = mean(arr[ring_mask_M_MD])
end

# ╔═╡ 69c26c48-e86b-41c1-98ee-f19338fdc9f9
mass_medium_md = score(arr[dilated_mask_M_MD], hu_calcium, hu_heart_tissue_medium_md, voxel_size, ρ_calcium, VolumeFraction())

# ╔═╡ 9c7e1b01-8d8e-466e-bbb1-c216d70415cc
md"""
## Low Density
"""

# ╔═╡ cf6888ee-04ed-42f7-b646-08363a8ca5a2
begin
    mask_M_LD_3D = Array{Bool}(undef, size(arr))
    for z in 1:size(arr, 3)
        mask_M_LD_3D[:, :, z] = mask_M_LD
    end
end;

# ╔═╡ 4b819875-6f6f-4f4e-8d8f-e3b9951818cb
md"""
#### Dilated mask
"""

# ╔═╡ 6f3668e0-d0e0-424f-bc57-44bc9174e931
# dilated_mask_M_LD = dilate(dilate(mask_M_LD_3D));

# ╔═╡ b1f718c8-13fc-4e9a-9a7e-50fd35eb532c
# ring_mask_M_LD =
#     Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask_M_LD)))))) - dilated_mask_M_LD);

# ╔═╡ a8f7f1b0-7b1c-4551-b3c3-0d0ac76e548e
begin
	dilated_mask_M_LD = dilate_mask_medium(mask_M_LD_3D)
	ring_mask_M_LD = ring_mask_medium(dilated_mask_M_LD)
end;

# ╔═╡ c6dcaea6-0c11-4245-b5af-63286f6ec860
@bind y2 overlay_mask_bind(dilated_mask_M_LD)

# ╔═╡ 3a269c36-c68d-461d-9e13-d758c6a10b8f
overlay_mask_plot(arr, dilated_mask_M_LD, y2, "dilated mask")

# ╔═╡ 9069d26b-540e-498a-b028-dcece011ceab
md"""
#### Ring (background) mask
"""

# ╔═╡ 7951fa92-01d9-4a82-9ca2-3779ac63c5e6
@bind y4 overlay_mask_bind(ring_mask_M_LD)

# ╔═╡ 95ca91c1-2c35-4a0e-8121-edfe6e6c62d4
overlay_mask_plot(arr, ring_mask_M_LD, y4, "ring mask")

# ╔═╡ bccf2b86-8938-4f36-982e-47d7a64ffb0a
begin
    single_ring_mask_M_LD = Bool.(ring_mask_M_LD[:, :, 3])
    hu_heart_tissue_medium_ld = mean(arr[ring_mask_M_LD])
end

# ╔═╡ e31f6f94-dcd9-495d-9d93-4fba7144e98f
mass_medium_ld = score(arr[dilated_mask_M_LD], hu_calcium, hu_heart_tissue_medium_ld, voxel_size, ρ_calcium, VolumeFraction())

# ╔═╡ ebbe8346-6dd6-40ca-a2af-cbd540810435
md"""
# Score Small Inserts
"""

# ╔═╡ d062507c-e4d1-4f40-8061-04cbf07cfe5b
function dilate_mask_small(mask)
	return (mask)
end

# ╔═╡ 04e1cb00-89bb-4328-b020-40d03604b5f9
function ring_mask_small(dilated_mask)
	return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask);
end

# ╔═╡ 9b9858c8-0c97-4d86-9ddc-60fa93d4f6c5
md"""
## High Density
"""

# ╔═╡ 6f9e78c9-ac2a-4674-9c4e-81579337115d
begin
    mask_S_HD_3D = Array{Bool}(undef, size(arr))
    for z in 1:size(arr, 3)
        mask_S_HD_3D[:, :, z] = mask_S_HD
    end
end;

# ╔═╡ d9c14524-049e-48d1-8a7c-fed24ba6c564
md"""
#### Dilated mask
"""

# ╔═╡ 2515aa75-0ce8-4eb5-a4d8-cc9d4803b688
# dilated_mask_S_HD = dilate(dilate(mask_S_HD_3D));

# ╔═╡ 73f9fa72-9dcd-4aaa-ac5c-11bdae66e69f
# ring_mask_S_HD =
#     Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask_S_HD)))))) - dilated_mask_S_HD);

# ╔═╡ a42fb7d6-1b18-4876-80c2-24de6cc6325c
begin
	dilated_mask_S_HD = dilate_mask_small(mask_S_HD_3D)
	ring_mask_S_HD = ring_mask_small(dilated_mask_S_HD)
end;

# ╔═╡ 0dfe8e4a-7664-4363-957b-54c1a7c8ef72
@bind l2 overlay_mask_bind(dilated_mask_S_HD)

# ╔═╡ 815c7f53-650c-49c1-893f-1865c9efaea6
overlay_mask_plot(arr, dilated_mask_S_HD, l2, "dilated mask")

# ╔═╡ c8672e42-6821-44aa-95b6-f084b5ce9e08
md"""
#### Ring (background) mask
"""

# ╔═╡ b25d71d7-bb3a-4f02-8692-85c998e02c1f
@bind l4 overlay_mask_bind(ring_mask_S_HD)

# ╔═╡ 9c56d5c0-479d-4e8e-bf1f-3673272ff1bb
overlay_mask_plot(arr, ring_mask_S_HD, l4, "ring mask")

# ╔═╡ c5ec8c8a-e5c0-4f29-a694-d06e04380aca
begin
    single_ring_mask_S_HD = Bool.(ring_mask_S_HD[:, :, 3])
    hu_heart_tissue_small_hd = mean(arr[ring_mask_S_HD])
end

# ╔═╡ ebb59f35-6962-4903-bd6d-93137f682c0e
mass_small_hd = score(arr[dilated_mask_S_HD], hu_calcium, hu_heart_tissue_large_hd, voxel_size, ρ_calcium, VolumeFraction())

# ╔═╡ 90c1c4ab-d199-4b8e-a5f0-d1a1ccb6e2e9
md"""
## Medium Density
"""

# ╔═╡ b2474626-da7e-4c94-b9c1-c912be1ca7c0
begin
    mask_S_MD_3D = Array{Bool}(undef, size(arr))
    for z in 1:size(arr, 3)
        mask_S_MD_3D[:, :, z] = mask_S_MD
    end
end;

# ╔═╡ 9d4395fb-2b7c-4b8f-835e-b65558e8ca4b
md"""
#### Dilated mask
"""

# ╔═╡ d76ca9f3-5804-4dc3-b285-37260d51a7e3
# dilated_mask_S_MD = dilate(dilate(mask_S_MD_3D));

# ╔═╡ 755da448-5f99-4e34-8340-886e7d8b58a5
# ring_mask_S_MD =
#     Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask_S_MD)))))) - dilated_mask_S_MD);

# ╔═╡ 079c8994-2f63-40b5-a0cf-0450e1ba6e71
begin
	dilated_mask_S_MD = dilate_mask_small(mask_S_MD_3D)
	ring_mask_S_MD = ring_mask_small(dilated_mask_S_MD)
end;

# ╔═╡ 6899c448-0d1d-438c-9c29-15b8f24f2573
@bind k2 overlay_mask_bind(dilated_mask_S_MD)

# ╔═╡ d6ba5b97-0c89-457e-ad52-27c45286e7f4
overlay_mask_plot(arr, dilated_mask_S_MD, k2, "dilated mask")

# ╔═╡ 65d2898e-b425-458b-b2c8-afb6c2d117be
md"""
#### Ring (background) mask
"""

# ╔═╡ db0edfc5-5215-426e-acc9-bb84c23979a9
@bind k4 overlay_mask_bind(ring_mask_S_MD)

# ╔═╡ ff9a4a3a-3d8a-444e-8b28-8a62e861341d
overlay_mask_plot(arr, ring_mask_S_MD, k4, "ring mask")

# ╔═╡ bc7547a4-a33d-4f5f-9dfd-db19ecc0f386
begin
    single_ring_mask_S_MD = Bool.(ring_mask_S_MD[:, :, 3])
    hu_heart_tissue_small_md = mean(arr[ring_mask_S_MD])
end

# ╔═╡ 966b91fa-3379-40a6-9178-07b729510021
mass_small_md = score(arr[dilated_mask_S_MD], hu_calcium, hu_heart_tissue_large_md, voxel_size, ρ_calcium, VolumeFraction())

# ╔═╡ a95d6831-ae36-488b-88b3-0678be493309
md"""
## Low Density
"""

# ╔═╡ 648a4846-3a84-44fc-8a34-11dd1873aecc
begin
    mask_S_LD_3D = Array{Bool}(undef, size(arr))
    for z in 1:size(arr, 3)
        mask_S_LD_3D[:, :, z] = mask_S_LD
    end
end;

# ╔═╡ 1dc477b9-5e7a-4f57-9c9b-0c568da32403
md"""
#### Dilated mask
"""

# ╔═╡ ce7e346b-9d05-4572-95b1-78d1109fafcd
# dilated_mask_S_LD = dilate(dilate(mask_S_LD_3D));

# ╔═╡ 248d0f15-e9e8-4000-8312-dfa2a06061c8
# ring_mask_S_LD =
#     Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask_S_LD)))))) - dilated_mask_S_LD);

# ╔═╡ d99d5975-c579-4d6d-ace1-7e7756c28872
begin
	dilated_mask_S_LD = dilate_mask_small(mask_S_LD_3D)
	ring_mask_S_LD = ring_mask_small(dilated_mask_S_LD)
end;

# ╔═╡ 12ca48fa-d240-4b20-ad3d-3b7e90e576af
@bind m2 overlay_mask_bind(dilated_mask_S_LD)

# ╔═╡ 8376358b-25ee-4d58-ab46-29a3551bce88
overlay_mask_plot(arr, dilated_mask_S_LD, m2, "dilated mask")

# ╔═╡ 0d0ccfed-0052-4181-b618-410dd907e9c1
md"""
#### Ring (background) mask
"""

# ╔═╡ b71170fb-ed36-4145-9230-1b4b9ff8eb47
@bind m4 overlay_mask_bind(ring_mask_S_LD)

# ╔═╡ 1c4f5f08-4395-4350-a24d-2e7c1a41689e
overlay_mask_plot(arr, ring_mask_S_LD, m4, "ring mask")

# ╔═╡ 3dafa251-4619-49cc-b2dd-5417943292ab
begin
    single_ring_mask_S_LD = Bool.(ring_mask_S_LD[:, :, 3])
    hu_heart_tissue_small_ld = mean(arr[ring_mask_S_LD])
end

# ╔═╡ d3dfc3f3-b003-425f-be16-96e2223c7c02
mass_small_ld = score(arr[dilated_mask_S_LD], hu_calcium, hu_heart_tissue_large_ld, voxel_size, ρ_calcium, VolumeFraction())

# ╔═╡ 45737edc-607c-47e1-a593-a0e21b1066ac
md"""
# Results
"""

# ╔═╡ 046d2b0b-1b65-463b-8467-14a288b0935b
if DENSITY == "low"
	density_array = [0.025, 0.050, 0.100]
elseif DENSITY == "normal"
	density_array = [0.200, 0.400, 0.800]
end

# ╔═╡ d39f3dec-0266-4325-bab1-56b1458009aa
PhantomSegmentation.get_pixel_size(header)

# ╔═╡ d0cca38e-3e0a-479a-b575-81026246b60b
inserts = ["Low Density", "Medium Density", "High Density"]

# ╔═╡ 23dd052e-3a84-4541-acba-15f6a1952402
volume_gt = [7.065, 63.585, 176.625]

# ╔═╡ 4c251a39-0570-4bc5-acfa-be495f8dff7c
ground_truth_mass_large = [
    volume_gt[3] * density_array[1],
    volume_gt[3] * density_array[2],
    volume_gt[3] * density_array[3],
] # mg

# ╔═╡ 0dff0a81-d695-4917-bb36-70bdb8ec5add
calculated_mass_large = [mass_large_ld, mass_large_md, mass_large_hd]

# ╔═╡ 86d787cd-7e72-48cc-bb1b-1c89eb09ee35
ground_truth_mass_medium = [
    volume_gt[2] * density_array[1],
    volume_gt[2] * density_array[2],
    volume_gt[2] * density_array[3],
] # mg

# ╔═╡ e3300105-8a84-415c-88e3-c4c7807d1e33
calculated_mass_medium = [mass_medium_ld, mass_medium_md, mass_medium_hd]

# ╔═╡ 73bf315d-c04e-4402-82fe-360c26516a72
ground_truth_mass_small = [
    volume_gt[1] * density_array[1],
    volume_gt[1] * density_array[2],
    volume_gt[1] * density_array[3],
] # mg

# ╔═╡ f8dc8bbe-6cde-4ee3-90e4-96f213c64274
calculated_mass_small = [mass_small_ld, mass_small_md, mass_small_hd]

# ╔═╡ e066e985-1e63-4e5b-af48-1f8275952a64
df = DataFrame(;
    DENSITY=DENSITY,
	ENERGY=ENERGY,
    inserts=inserts,
    ground_truth_mass_large=ground_truth_mass_large,
    calculated_mass_large=calculated_mass_large,
    ground_truth_mass_medium=ground_truth_mass_medium,
    calculated_mass_medium=calculated_mass_medium,
    ground_truth_mass_small=ground_truth_mass_small,
    calculated_mass_small=calculated_mass_small,
)

# ╔═╡ 429a62af-e281-401f-98c0-c29958a6bcf9
percent_error_large =
    (abs.(ground_truth_mass_large - calculated_mass_large) ./ ground_truth_mass_large) .*
    100

# ╔═╡ 37caaa56-d9f6-4bee-8cef-1d38842c29ef
percent_error_medium =
    (abs.(ground_truth_mass_medium - calculated_mass_medium) ./ ground_truth_mass_medium) .*
    100

# ╔═╡ 9a373872-2cdc-4c74-8059-101db6ccd0db
percent_error_small =
    (abs.(ground_truth_mass_small - calculated_mass_small) ./ ground_truth_mass_small) .*
    100

# ╔═╡ 62510b39-6331-44ce-b240-e1df940dac43
md"""
### Save Results
"""

# ╔═╡ 1b67c78b-375c-47aa-8b45-5ac88c242bec
# if ~isdir(string(cd(pwd, "..") , "/output/", TYPE))
# 	mkdir(string(cd(pwd, "..") , "/output/", TYPE))
# end

# ╔═╡ a886c3d2-51c0-4df5-b82b-a61689e46422
# output_path = string(cd(pwd, "..") , "/output/", TYPE, "/", scan, ".csv")

# ╔═╡ 73a3746c-8ee6-457b-b582-297f7bf95dae
# CSV.write(output_path, df)

# ╔═╡ 9ffd9bf1-5965-452e-9c83-e6c60f045a38
md"""
### Save full df
"""

# ╔═╡ 53c38bb2-e166-45b3-880b-0762e2d98683
dfs = []

# ╔═╡ 54fde0dc-62a6-4064-bd05-501b349c2a25
push!(dfs, df)

# ╔═╡ e8fa3728-bc90-4b4f-9a8f-80b89cd7002b
# if length(dfs) == 24
# 	global new_df = vcat(dfs[1:24]...)
# 	output_path_new = string(cd(pwd, "..") , "/output/", TYPE, "/", "full.csv")
# 	CSV.write(output_path_new, new_df)
# end



# ╔═╡ 3f86acd8-8f4d-46f4-8856-574b65c83386


# ╔═╡ Cell order:
# ╠═a2647e62-e189-498a-981c-9980485fceca
# ╠═f4e4decf-d8bf-4203-b44a-0d2cb9830a3a
# ╠═d161adaf-c02a-49d1-a5dd-6929ea0614c8
# ╟─bb5f1496-9b8f-4d26-8f8d-03773fa546a9
# ╠═f1777984-c7e3-4fd2-92ba-7649d32caab8
# ╠═7b6e134b-2249-479d-a3f6-f0a2791a25f7
# ╟─44a31da5-9914-44aa-b3ef-551ab3ea30ae
# ╠═3de3c36e-6350-4aab-a139-92d052641136
# ╠═9155bc22-344b-4ff1-b629-ce7b935a3370
# ╠═88ab1259-baa7-422a-9321-0f8f3d40b051
# ╠═f45b6ec1-330d-4e6e-a36f-f10002489f25
# ╟─38c4b7f8-1c7f-4776-922b-2a74dd4ffb3d
# ╟─9f43b7c7-e99c-4893-aa4b-2252dc3ba91a
# ╟─a352f0bf-f2a3-420f-ad78-0e3a8fab20ae
# ╟─003e03c1-3796-4f58-853c-055bbc2dfb37
# ╠═66dc4873-d94c-4d6e-b08f-5dc31ae3f19b
# ╠═64717903-e79e-4e4d-9c16-f780770b786a
# ╟─64eca4c3-d079-4c4e-9070-f39095d97163
# ╠═12205cce-c9cf-47df-b929-92eb4f72309b
# ╟─512ae1c7-0184-4fe4-a492-899e407711ee
# ╟─a159ae46-dfce-4a0a-80ab-d1695a9acefd
# ╠═b097a0aa-f5e3-48e5-ad32-275723015e1c
# ╠═d0bc44a9-cf63-419f-b291-a6f40fd0d8fb
# ╠═7012f3c3-4a49-4ff0-8c9a-248b5e9fb8a3
# ╟─1d20bd27-1c22-403b-b02b-bc30f282e819
# ╟─b7e9d017-da52-479a-b804-7cd74fd1af5e
# ╟─6722dbc9-0d9a-4570-89e2-d2fba04b7b5b
# ╟─8e4c21d5-50d0-49a1-9342-8cf027c40d6b
# ╟─b2a82830-7240-49d4-9503-7a817f138d28
# ╠═370115da-6071-4abe-afce-314bca318ed6
# ╠═2e2953cb-cea8-4cf9-afdd-872ead7da567
# ╠═0b88e378-cda1-4400-8005-51e25fcef33e
# ╠═33408dd8-e97a-45c6-a99e-1709d30c924c
# ╟─24e19f9a-c820-4cdc-b3bb-cc081d1ca0e7
# ╠═d0941946-4338-4b8b-b4cc-e83f8780a253
# ╟─4db2de22-7cf4-4f45-a895-76547816e1dd
# ╠═9ae8c021-f465-4e55-ae5a-43cfb08a7e01
# ╠═e3dbd551-e3f5-461b-931e-90c3de086270
# ╠═4a8fe98a-d2d1-4fe2-a4d4-5a1898c49f84
# ╟─0cbc4ec2-790f-45a0-801d-f04c25657bfb
# ╠═04e06029-e368-4ee4-8394-ad4333e68232
# ╠═02010cf7-7bf8-44dc-a01b-8f5d57df5733
# ╠═d20b5d55-e8fb-43cd-b447-4b6cf558db6d
# ╠═84419df2-9512-4dcb-b4af-e5d194fd65db
# ╠═60e001ba-cabd-4df8-bb48-827be2332513
# ╠═d274581a-ae43-4398-b5ee-d85b980b56cc
# ╠═16e0bc1a-7b56-4f91-8b72-4677b782d3fd
# ╠═07d2a022-e61d-4737-b586-9816bcbd72d1
# ╠═9995fadd-7291-46e4-be27-b7e1e8dbce0e
# ╟─2ed40040-bfec-4ad4-b87f-29caccc88712
# ╠═78587458-4711-4c2b-8755-ed4030a0a18e
# ╠═aef27e7c-8164-469b-b5ed-712bdf2cd645
# ╠═e2f0e863-8f4b-46ab-af9c-9d59538a9eec
# ╠═27eb5242-330a-4f7d-aad9-9067a3e869ea
# ╟─bab453c4-0057-42e3-a009-02c5490f847f
# ╠═7ede2237-aaec-4fc3-86aa-bade846dad15
# ╟─e8cbc4c1-bc45-4065-afb0-014eaf4b3638
# ╠═f1e04fa6-f5bd-42c9-a5d0-93960c3e6b92
# ╟─d4f4cd88-38f1-4801-8973-32cdc6e172b6
# ╠═2d0d5ff5-b2f3-4d44-be7b-50772222f277
# ╟─28324392-eecf-4865-a788-a871decc255d
# ╠═1608d803-cf85-4f99-84d4-0b89b90dea34
# ╠═9594574b-87aa-49da-822c-1d9a835f21c1
# ╠═52c0720b-6025-477d-9bb9-eba8e9457bb3
# ╠═f0eb66fd-95ab-466d-b201-67308e6f53db
# ╟─5377fc95-b66a-42f6-b1a7-5d7c83a0714a
# ╠═fb6f86d5-531c-4a28-b225-5d4c1b134d82
# ╟─fa0acc1f-c333-4ead-bc03-0a748466d79d
# ╠═f4f4c29e-6838-405b-b65b-1e0bed3d9854
# ╠═58a299a7-e138-4c76-a016-a27cf71697a7
# ╠═525199de-753c-44d6-967d-f5d3cbc30fef
# ╟─28a1deed-248e-4b15-896c-ff75dba4ddd2
# ╠═b91850d3-aa07-461d-ae78-a805b0d4a23e
# ╟─6b8af2c6-fd50-44f2-b1b8-93d083e1bb21
# ╟─6adf7bde-eaba-4480-8229-7d56fd765f68
# ╠═eadf72ff-6d0c-4e5f-a32e-dc42cb83f93e
# ╠═2c966c01-56db-4661-a0e5-276ba44b02e8
# ╠═4ab5580a-d5f2-42df-ade5-53b0cb0c9293
# ╟─e32e2ddb-cc73-4474-8711-d7fd105b359d
# ╠═b5fd1d9c-de0a-4e74-8808-ce57c777c675
# ╟─14d00a87-f964-450e-b631-7a36c40d6834
# ╠═b209a40f-2434-4940-973e-e6e3debbbbe4
# ╠═ceac8866-7932-4bf3-a2ee-eb0393b4a036
# ╠═ebe4a819-0466-40fc-b29d-94fe2dad15a3
# ╟─257efcfc-0c72-4b68-9a32-8b04b39da835
# ╠═69ff66e2-afa3-44df-989b-7497e390b252
# ╟─69513ba5-9e1a-4932-b056-eac183aa3ec5
# ╟─2d9fb8ed-3b2a-4e6c-880b-9d2d1627bf60
# ╠═7fa3423e-6ad0-470f-b1b3-9e62f9533ae6
# ╠═260258de-360f-410f-b43c-bdfb72be5322
# ╠═2d43321b-2836-4887-a4d6-0e6a99fbc0f2
# ╟─5f475234-cd14-47dc-b5a2-f4c8b0917849
# ╠═38b812bd-7291-4dc6-9838-9e0e65e59a0f
# ╠═300e41a3-bd66-4b1f-a225-e1b2b8cf39d7
# ╟─046a03f1-a285-4a80-bd64-e821252d07f0
# ╠═9a201bfe-beb9-4729-b8af-c4186bcd0d3b
# ╟─c76ad7d6-76a1-42b8-b5ef-3ebd342177a9
# ╠═57c2b58f-6bef-4088-a734-9fc21ddb6217
# ╠═6e833dc6-412f-4273-b658-0a1e27850d55
# ╠═34b427ef-d767-4216-b570-928ef05c7258
# ╟─b992a000-8138-4c70-9b38-342bfb4460a7
# ╠═5452a3a7-ab92-4117-b357-cc7c806b6a83
# ╟─adab3fe3-8ff7-49af-9aa0-c0e6e631fed3
# ╟─ddbf8772-16e0-4ad1-a950-779a5d1c8a44
# ╠═d3b28cdf-909d-45ad-bcd2-46524bc527f6
# ╠═fb0a7476-684c-4c6c-a522-5bf1ea5c5247
# ╠═e51f9107-94a5-490c-a03b-14ea5aadbb47
# ╟─dcf534a7-8545-4a9f-87a8-5ca031541f43
# ╠═08ebf130-56fa-4b5f-8cc7-3d25e2326a00
# ╟─dcf7e7b4-0302-48f9-818f-f059c2873b2f
# ╠═60c3a5ac-0a81-4352-b2e1-3ea89513b0ee
# ╠═c77d924d-4496-4734-9e08-bb19d8a01a33
# ╠═2bef5e72-b48d-41e8-b6c3-846cefc2281d
# ╟─18cbce65-2543-4210-896a-b10e80ba7d95
# ╠═56ede52c-706f-41f7-bae1-059ab0b3bb82
# ╟─6b8797f1-4d2d-4b64-ae3e-fb436b2ab8a2
# ╟─f0d2a0b9-74be-4ca7-8efd-208a5010cac2
# ╠═4c3db4e5-eacd-4e5e-a6e8-a35504d9c01a
# ╠═2f02bdcb-1d2e-4313-b338-2cac8143662b
# ╠═69c26c48-e86b-41c1-98ee-f19338fdc9f9
# ╟─9c7e1b01-8d8e-466e-bbb1-c216d70415cc
# ╠═cf6888ee-04ed-42f7-b646-08363a8ca5a2
# ╟─4b819875-6f6f-4f4e-8d8f-e3b9951818cb
# ╠═6f3668e0-d0e0-424f-bc57-44bc9174e931
# ╠═b1f718c8-13fc-4e9a-9a7e-50fd35eb532c
# ╟─a8f7f1b0-7b1c-4551-b3c3-0d0ac76e548e
# ╟─c6dcaea6-0c11-4245-b5af-63286f6ec860
# ╠═3a269c36-c68d-461d-9e13-d758c6a10b8f
# ╟─9069d26b-540e-498a-b028-dcece011ceab
# ╟─7951fa92-01d9-4a82-9ca2-3779ac63c5e6
# ╠═95ca91c1-2c35-4a0e-8121-edfe6e6c62d4
# ╠═bccf2b86-8938-4f36-982e-47d7a64ffb0a
# ╠═e31f6f94-dcd9-495d-9d93-4fba7144e98f
# ╟─ebbe8346-6dd6-40ca-a2af-cbd540810435
# ╠═d062507c-e4d1-4f40-8061-04cbf07cfe5b
# ╠═04e1cb00-89bb-4328-b020-40d03604b5f9
# ╟─9b9858c8-0c97-4d86-9ddc-60fa93d4f6c5
# ╠═6f9e78c9-ac2a-4674-9c4e-81579337115d
# ╟─d9c14524-049e-48d1-8a7c-fed24ba6c564
# ╠═2515aa75-0ce8-4eb5-a4d8-cc9d4803b688
# ╠═73f9fa72-9dcd-4aaa-ac5c-11bdae66e69f
# ╠═a42fb7d6-1b18-4876-80c2-24de6cc6325c
# ╟─0dfe8e4a-7664-4363-957b-54c1a7c8ef72
# ╠═815c7f53-650c-49c1-893f-1865c9efaea6
# ╟─c8672e42-6821-44aa-95b6-f084b5ce9e08
# ╟─b25d71d7-bb3a-4f02-8692-85c998e02c1f
# ╠═9c56d5c0-479d-4e8e-bf1f-3673272ff1bb
# ╠═c5ec8c8a-e5c0-4f29-a694-d06e04380aca
# ╠═ebb59f35-6962-4903-bd6d-93137f682c0e
# ╟─90c1c4ab-d199-4b8e-a5f0-d1a1ccb6e2e9
# ╠═b2474626-da7e-4c94-b9c1-c912be1ca7c0
# ╟─9d4395fb-2b7c-4b8f-835e-b65558e8ca4b
# ╠═d76ca9f3-5804-4dc3-b285-37260d51a7e3
# ╠═755da448-5f99-4e34-8340-886e7d8b58a5
# ╠═079c8994-2f63-40b5-a0cf-0450e1ba6e71
# ╟─6899c448-0d1d-438c-9c29-15b8f24f2573
# ╠═d6ba5b97-0c89-457e-ad52-27c45286e7f4
# ╟─65d2898e-b425-458b-b2c8-afb6c2d117be
# ╟─db0edfc5-5215-426e-acc9-bb84c23979a9
# ╠═ff9a4a3a-3d8a-444e-8b28-8a62e861341d
# ╠═bc7547a4-a33d-4f5f-9dfd-db19ecc0f386
# ╠═966b91fa-3379-40a6-9178-07b729510021
# ╟─a95d6831-ae36-488b-88b3-0678be493309
# ╠═648a4846-3a84-44fc-8a34-11dd1873aecc
# ╟─1dc477b9-5e7a-4f57-9c9b-0c568da32403
# ╠═ce7e346b-9d05-4572-95b1-78d1109fafcd
# ╠═248d0f15-e9e8-4000-8312-dfa2a06061c8
# ╠═d99d5975-c579-4d6d-ace1-7e7756c28872
# ╟─12ca48fa-d240-4b20-ad3d-3b7e90e576af
# ╠═8376358b-25ee-4d58-ab46-29a3551bce88
# ╟─0d0ccfed-0052-4181-b618-410dd907e9c1
# ╠═b71170fb-ed36-4145-9230-1b4b9ff8eb47
# ╠═1c4f5f08-4395-4350-a24d-2e7c1a41689e
# ╠═3dafa251-4619-49cc-b2dd-5417943292ab
# ╠═d3dfc3f3-b003-425f-be16-96e2223c7c02
# ╟─45737edc-607c-47e1-a593-a0e21b1066ac
# ╠═046d2b0b-1b65-463b-8467-14a288b0935b
# ╠═d39f3dec-0266-4325-bab1-56b1458009aa
# ╠═d0cca38e-3e0a-479a-b575-81026246b60b
# ╠═23dd052e-3a84-4541-acba-15f6a1952402
# ╠═4c251a39-0570-4bc5-acfa-be495f8dff7c
# ╠═0dff0a81-d695-4917-bb36-70bdb8ec5add
# ╠═86d787cd-7e72-48cc-bb1b-1c89eb09ee35
# ╠═e3300105-8a84-415c-88e3-c4c7807d1e33
# ╠═73bf315d-c04e-4402-82fe-360c26516a72
# ╠═f8dc8bbe-6cde-4ee3-90e4-96f213c64274
# ╠═e066e985-1e63-4e5b-af48-1f8275952a64
# ╠═429a62af-e281-401f-98c0-c29958a6bcf9
# ╠═37caaa56-d9f6-4bee-8cef-1d38842c29ef
# ╠═9a373872-2cdc-4c74-8059-101db6ccd0db
# ╟─62510b39-6331-44ce-b240-e1df940dac43
# ╠═1b67c78b-375c-47aa-8b45-5ac88c242bec
# ╠═a886c3d2-51c0-4df5-b82b-a61689e46422
# ╠═73a3746c-8ee6-457b-b582-297f7bf95dae
# ╟─9ffd9bf1-5965-452e-9c83-e6c60f045a38
# ╠═53c38bb2-e166-45b3-880b-0762e2d98683
# ╠═54fde0dc-62a6-4064-bd05-501b349c2a25
# ╠═e8fa3728-bc90-4b4f-9a8f-80b89cd7002b
# ╠═3f86acd8-8f4d-46f4-8856-574b65c83386

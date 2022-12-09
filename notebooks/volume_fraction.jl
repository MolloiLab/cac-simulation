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

# ╔═╡ 3f9e21d6-a663-49ce-85a3-e5943b0be245
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.activate(".")

    using PlutoUI, Statistics, CSV, DataFrames, GLM, CairoMakie, HypothesisTests, Colors, MLJBase, DICOM, DICOMUtils, PhantomSegmentation, CalciumScoring, ImageMorphology, ImageFiltering, Noise
    using StatsBase: quantile!, rmsd
end

# ╔═╡ 70fea9a9-8d77-40c4-84ea-b084d72b1b4c
TableOfContents()

# ╔═╡ e3a4316d-481c-4e2c-a592-03727ae494e8
pwd()

# ╔═╡ fea757e2-e0c6-48eb-b997-cbabd91399bf
md"""
## Load DICOMS

All you need to do is set `base_path` once and leave it. After that, the only thing that should change is the `VENDER`, once for every set, and the `SCAN_NUMBER`, once for each scan.
"""

# ╔═╡ 3e211eda-2d67-45f9-a7ee-1825bd746cbb
begin
    SCAN_NUMBER = 1
	venders = ["80", "100", "120", "135"]
    VENDER = venders[1]
    SIZE = "small"
    # SIZE = "medium"
    # SIZE = "large"
    # DENSITY = "low"
    DENSITY = "normal"
    TYPE = "volume_fraction"
    BASE_PATH = string(
        "/Users/daleblack/Google Drive/dev/MolloiLab/cac-simulation/images_new/",
        SIZE,
        "/",
        DENSITY,
        "/",
    )
end

# ╔═╡ c62f5af3-a3d6-458b-bc72-3ac0f3372323
md"""
**Everything below should be automatic, just scroll through to visually inspect that things make sense**
"""

# ╔═╡ e35a7c50-0098-47ec-bad0-84e4c36011e8
begin
	root_path = string(BASE_PATH, VENDER)
	dcm_path_list = dcm_list_builder(root_path)
	pth = dcm_path_list[SCAN_NUMBER]
	scan = basename(pth)
	header, dcm_array, slice_thick_ori1 = dcm_reader(pth)
end;

# ╔═╡ 6c5111f0-de83-4e64-800c-2f9aefb5ee4f
md"""
## Helper Functions
"""

# ╔═╡ 3b8ee7a8-917c-4836-967b-ee7e51ad0331
function collect_tuple(tuple_array)
    row_num = size(tuple_array)
    col_num = length(tuple_array[1])
    container = zeros(Int64, row_num..., col_num)
    for i in 1:length(tuple_array)
        container[i, :] = collect(tuple_array[i])
    end
    return container
end

# ╔═╡ 67fd63c4-bac1-4ebc-9299-d17850e35d6d
function overlay_mask_bind(mask)
    indices = findall(x -> x == 1, mask)
    indices = Tuple.(indices)
    label_array = collect_tuple(indices)
    zs = unique(label_array[:, 3])
    return PlutoUI.Slider(1:length(zs); default=3, show_value=true)
end

# ╔═╡ 899a8652-c961-4e65-93ec-d75f0a082e87
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

# ╔═╡ 5cbb0bbb-96b9-4a2f-80d8-f582d4f822e0
md"""
## Segment Heart
"""

# ╔═╡ 5b663e10-4498-4f11-aaf9-7c49b25b9fe1
masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3) ÷ 2);

# ╔═╡ 05ae002a-7a61-473e-aa06-9ebfcc51f623
@bind a PlutoUI.Slider(1:size(masked_array, 3), default=10, show_value=true)

# ╔═╡ 8018fb0b-5d1a-4f50-904f-b2686c75993f
heatmap(transpose(masked_array[:, :, a]); colormap=:grays)

# ╔═╡ f20119be-1939-48a7-aa8d-507003b4ebb3
xs, yx = 185:235, 350:400

# ╔═╡ 28f2e565-fde8-40d7-a5fd-9b944e77efb2
heatmap(transpose(masked_array[xs, yx, a]); colormap=:grays)

# ╔═╡ fcf12109-023a-4c6c-82a2-64471bf8ac1b
# let
#     fig = Figure()

#     ax = Makie.Axis(fig[1, 1])
#     heatmap!(transpose(dcm_array[:, :, 5]); colormap=:grays)

# 	save("/Users/daleblack/Google Drive/Research/Papers/My Papers/cac-simulation/figures-review/simulated_phantom.png", fig)
#     fig
# end

# ╔═╡ 5c3e4fe5-1ad8-4cff-8382-9fda0a05d168
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

# ╔═╡ b925bdde-25e5-4560-bfef-0d822011f781
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

# ╔═╡ e836390f-1647-4681-87cc-242b0c48932e
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

# ╔═╡ 12973897-c65a-496c-bb8a-af665c14f2a7
md"""
## Segment Calcium Rod
"""

# ╔═╡ 099c6fdd-0f98-4d6c-9121-33c0ef0e84c5
if DENSITY == "low"
	thresh = 55
elseif DENSITY == "normal"
	thresh = 130
end

# ╔═╡ a85ae81e-daae-4b60-8443-fe5d53309adb
calcium_image, slice_CCI, quality_slice, cal_rod_slice = mask_rod(
    masked_array, header; calcium_threshold=thresh
);

# ╔═╡ f858d9aa-5282-401c-a5ec-d8c0a0b43238
array_filtered = abs.(mapwindow(median, calcium_image[:, :, 2], (3, 3)));

# ╔═╡ 5ab916de-71b9-4b73-9c3e-57ea89b608bb
c_img = calcium_image[:, :, 1:3];

# ╔═╡ 97e09a38-569c-4cab-85dc-15f05cc3e0ad
@bind c PlutoUI.Slider(1:size(calcium_image, 3), default=5, show_value=true)

# ╔═╡ 7fad421c-6a53-41a4-bfaa-59c25c4d4d09
# heatmap(transpose(calcium_image[:, :, c]); colormap=:grays)

# ╔═╡ 59e532c9-1550-4c20-82ea-a9df5d5dc4a5
md"""
## Load (segment) Calcium Inserts
"""

# ╔═╡ f5cfbc09-2455-4c5f-84bc-d704f54fee69
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

# ╔═╡ 542d221b-d765-4bf7-b6ca-5d3a2f5c2e52
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

# ╔═╡ d7d2fadd-6bca-4132-9093-27ef4502722b
# heatmap(transpose(masks); colormap=:grays)

# ╔═╡ 852eda80-17c6-42c6-a015-032afaadc44e
md"""
## Calibration Prep
"""

# ╔═╡ 9e3eff5b-6bc3-4799-aea8-311709caa170
bool_arr = array_filtered .> 0;

# ╔═╡ 60b42dec-ccfe-49e8-9cb0-c73ce6455b00
bool_arr_erode = (erode(erode(erode(erode(bool_arr)))));

# ╔═╡ 88fa07e8-430b-42ce-aa09-3bdc9a79a872
# heatmap(bool_arr; colormap=:grays)

# ╔═╡ b3a43484-00db-4837-845d-620ed4c96ff1
# heatmap(bool_arr_erode; colormap=:grays)

# ╔═╡ b4564806-2dd4-48fe-9cbf-39010fecb559
begin
    mask_cal_3D = Array{Bool}(undef, size(c_img))
    for z in 1:size(c_img, 3)
        mask_cal_3D[:, :, z] = bool_arr_erode
    end
end;

# ╔═╡ c2920879-e886-45bd-814b-96aba0af7f97
hu_calcium = mean(c_img[mask_cal_3D])

# ╔═╡ 5c981f22-89d9-436d-a94a-581250ea1ec2
# hist(c_img[mask_cal_3D])

# ╔═╡ 1ec199c4-6b17-4c58-b959-4904c320244d
ρ_calcium = 0.2

# ╔═╡ a99e68e4-9c0d-4f0d-a164-8ab67d5839ef
# heatmap(single_arr)

# ╔═╡ 1dbad5ed-2493-4e39-8ec9-1e5e3720ddf7
md"""
# Score Large Inserts
"""

# ╔═╡ 839bf6ad-3a75-4134-b15d-b612cf7e3ff8
arr = masked_array[:, :, 4:6];

# ╔═╡ 2bec75b1-c15b-42ce-a1c8-42af149d491e
single_arr = masked_array[:, :, slice_CCI];

# ╔═╡ da3bdc26-aeb3-4a92-86dc-5c3653bc41cd
hu_heart_tissue = mean(single_arr[center_insert[1]-10:center_insert[1]+10, center_insert[2]-10:center_insert[2]+10])

# ╔═╡ ba48ec9f-8970-4104-85d5-6f50e2a4c16c
begin
	pixel_size = DICOMUtils.get_pixel_size(header)
	voxel_size = pixel_size[1] * pixel_size[2] * pixel_size[3]
end

# ╔═╡ 8dedc517-d985-4658-b6b7-cb8322b79d03
function dilate_mask_large(mask)
	return dilate(mask)
end

# ╔═╡ c393ffaa-7e87-48af-8e2d-274dcef25447
function ring_mask_large(dilated_mask)
	return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask);
end

# ╔═╡ 2b9cc53e-bcd3-407b-bf90-c65947b28581
md"""
## High Density
"""

# ╔═╡ 77e907a4-164a-45a0-aed7-a2d599b688a4
begin
    mask_L_HD_3D = Array{Bool}(undef, size(arr))
    for z in 1:size(arr, 3)
        mask_L_HD_3D[:, :, z] = mask_L_HD
    end
end;

# ╔═╡ ddb290bf-b754-4b41-8430-963bf835ae48
md"""
#### Dilated mask
"""

# ╔═╡ 0d3b15df-f115-47e7-a15a-754025163025
begin
	dilated_mask_L_HD = dilate_mask_large(mask_L_HD_3D)
	ring_mask_L_HD = ring_mask_large(dilated_mask_L_HD)
end;

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
begin
    single_ring_mask_L_HD = Bool.(ring_mask_L_HD[:, :, 3])
    hu_heart_tissue_large_hd = mean(arr[ring_mask_L_HD])
end

# ╔═╡ bf788e3f-538e-425e-a06a-6e4499259fee
mass_large_hd = score(arr[dilated_mask_L_HD], hu_calcium, hu_heart_tissue_large_hd, voxel_size, ρ_calcium, VolumeFraction())

# ╔═╡ 1b0d7eeb-d6dc-4c38-b5bb-624b54235b4b
md"""
## Medium Density
"""

# ╔═╡ 188742b7-dc79-496c-9c2c-948937a7ced0
begin
    mask_L_MD_3D = Array{Bool}(undef, size(arr))
    for z in 1:size(arr, 3)
        mask_L_MD_3D[:, :, z] = mask_L_MD
    end
end;

# ╔═╡ c1fea4b2-7740-4d62-b908-d89261d25064
md"""
#### Dilated mask
"""

# ╔═╡ 1a91c447-08cc-49af-ade1-b52b8a8258f0
begin
	dilated_mask_L_MD = dilate_mask_large(mask_L_MD_3D)
	ring_mask_L_MD = ring_mask_large(dilated_mask_L_MD)
end;

# ╔═╡ 4c4ce42b-4a64-4f81-96e2-6da91186b62f
# dilated_mask_L_MD = dilate(dilate(mask_L_MD_3D));

# ╔═╡ 968ff813-9218-4c43-946e-543301c0278e
# ring_mask_L_MD =
#     Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask_L_MD)))))) - dilated_mask_L_MD);

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
begin
    single_ring_mask_L_MD = Bool.(ring_mask_L_MD[:, :, 3])
    hu_heart_tissue_large_md = mean(arr[ring_mask_L_MD])
end

# ╔═╡ 42a5aa69-8d08-4eba-a952-00650bad6905
mass_large_md = score(arr[dilated_mask_L_MD], hu_calcium, hu_heart_tissue_large_md, voxel_size, ρ_calcium, VolumeFraction())

# ╔═╡ 39215a15-ca20-4eb4-ac7c-07bb022e410b
md"""
## Low Density
"""

# ╔═╡ 442282b4-4be9-481e-a341-f94637e71c9a
begin
    mask_L_LD_3D = Array{Bool}(undef, size(arr))
    for z in 1:size(arr, 3)
        mask_L_LD_3D[:, :, z] = mask_L_LD
    end
end;

# ╔═╡ bb4670bf-0c80-4039-b610-2db7103efd1e
md"""
#### Dilated mask
"""

# ╔═╡ 50c43516-f0e6-469b-9dd6-c9be70e12a7b
# dilated_mask_L_LD = dilate(dilate(mask_L_LD_3D));

# ╔═╡ 42ba3ea2-c34c-4b5a-8b13-5dea4981c830
# ring_mask_L_LD =
#     Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask_L_LD)))))) - dilated_mask_L_LD);

# ╔═╡ bfd37f8f-987e-445f-adbb-e036bf9f2839
begin
	dilated_mask_L_LD = dilate_mask_large(mask_L_LD_3D)
	ring_mask_L_LD = ring_mask_large(dilated_mask_L_LD)
end;

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
begin
    single_ring_mask_L_LD = Bool.(ring_mask_L_LD[:, :, 3])
    hu_heart_tissue_large_ld = mean(arr[ring_mask_L_LD])
end

# ╔═╡ c70bc1a1-7ffe-43bb-8955-874d6fc3550d
mass_large_ld = score(arr[dilated_mask_L_LD], hu_calcium, hu_heart_tissue_large_ld, voxel_size, ρ_calcium, VolumeFraction())

# ╔═╡ a2577753-a937-47f8-b14f-8c2309940bae
md"""
# Score Medium Inserts
"""

# ╔═╡ 387eeef3-a7fb-4a0c-a0f7-9fc2c30c170e
function dilate_mask_medium(mask)
	return (mask)
end

# ╔═╡ 718ec3cb-3fed-49d2-bf49-0303a4514ee6
function ring_mask_medium(dilated_mask)
	return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask);
end

# ╔═╡ e9f53bce-0d59-456f-8add-ec3e7e0aa47e
md"""
## High Density
"""

# ╔═╡ 179cb850-1365-4b5f-87a6-c0b3432d4708
begin
    mask_M_HD_3D = Array{Bool}(undef, size(arr))
    for z in 1:size(arr, 3)
        mask_M_HD_3D[:, :, z] = mask_M_HD
    end
end;

# ╔═╡ 4a65c23f-942c-45e0-ae58-0712240e9668
md"""
#### Dilated mask
"""

# ╔═╡ b24be6e5-d41a-465c-be01-0c274c7e2124
# dilated_mask_M_HD = dilate(dilate(mask_M_HD_3D));

# ╔═╡ 206e570f-5a65-43b1-b999-86e75b9cdf9d
# ring_mask_M_HD =
#     Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask_M_HD)))))) - dilated_mask_M_HD);

# ╔═╡ 1d4fb6c8-225a-467c-ae32-8c8180796f2b
begin
	dilated_mask_M_HD = dilate_mask_medium(mask_M_HD_3D)
	ring_mask_M_HD = ring_mask_medium(dilated_mask_M_HD)
end;

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
begin
    single_ring_mask_M_HD = Bool.(ring_mask_M_HD[:, :, 3])
    hu_heart_tissue_medium_hd = mean(arr[ring_mask_M_HD])
end

# ╔═╡ 93c6c475-cbae-442c-aceb-0b7127b268dd
mass_medium_hd = score(arr[dilated_mask_M_HD], hu_calcium, hu_heart_tissue_medium_hd, voxel_size, ρ_calcium, VolumeFraction())

# ╔═╡ 945a529f-c24e-4122-886f-1402672ed625
md"""
## Medium Density
"""

# ╔═╡ 37261e95-caef-4b6e-b514-783fd0c72dfa
begin
    mask_M_MD_3D = Array{Bool}(undef, size(arr))
    for z in 1:size(arr, 3)
        mask_M_MD_3D[:, :, z] = mask_M_MD
    end
end;

# ╔═╡ 040a62b8-7c41-4a39-8ecd-0a6242d3ff5b
md"""
#### Dilated mask
"""

# ╔═╡ 10dfc8bf-cb1c-4907-a613-dd6d1c1782ab
# dilated_mask_M_MD = dilate(dilate(mask_M_MD_3D));

# ╔═╡ 6785b8f9-e348-4cbb-9329-f8bf3d3e89cd
# ring_mask_M_MD =
#     Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask_M_MD)))))) - dilated_mask_M_MD);

# ╔═╡ 1e7995fd-29b3-4bf9-9abf-e534e95d8375
begin
	dilated_mask_M_MD = dilate_mask_medium(mask_M_MD_3D)
	ring_mask_M_MD = ring_mask_medium(dilated_mask_M_MD)
end;

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
begin
    single_ring_mask_M_MD = Bool.(ring_mask_M_MD[:, :, 3])
    hu_heart_tissue_medium_md = mean(arr[ring_mask_M_MD])
end

# ╔═╡ 2a152532-109a-461d-8eec-3a37ed236161
mass_medium_md = score(arr[dilated_mask_M_MD], hu_calcium, hu_heart_tissue_medium_md, voxel_size, ρ_calcium, VolumeFraction())

# ╔═╡ 3632354d-4850-45d6-a9d1-15949f8c0eb4
md"""
## Low Density
"""

# ╔═╡ 2b3d1930-1251-4435-ac36-09fe19b266d4
begin
    mask_M_LD_3D = Array{Bool}(undef, size(arr))
    for z in 1:size(arr, 3)
        mask_M_LD_3D[:, :, z] = mask_M_LD
    end
end;

# ╔═╡ 62b7d043-47fe-4aa2-8936-2739113be27d
md"""
#### Dilated mask
"""

# ╔═╡ 04d1cc9d-9a57-47d4-902f-ed3787ef1cd9
# dilated_mask_M_LD = dilate(dilate(mask_M_LD_3D));

# ╔═╡ 6f77adfd-6577-4760-ab87-7078bc4cff99
# ring_mask_M_LD =
#     Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask_M_LD)))))) - dilated_mask_M_LD);

# ╔═╡ d27c155a-1abd-4fb2-a229-a7110827ffda
begin
	dilated_mask_M_LD = dilate_mask_medium(mask_M_LD_3D)
	ring_mask_M_LD = ring_mask_medium(dilated_mask_M_LD)
end;

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
begin
    single_ring_mask_M_LD = Bool.(ring_mask_M_LD[:, :, 3])
    hu_heart_tissue_medium_ld = mean(arr[ring_mask_M_LD])
end

# ╔═╡ ec9dfcd9-13eb-41a7-be8c-405483d84a92
mass_medium_ld = score(arr[dilated_mask_M_LD], hu_calcium, hu_heart_tissue_medium_ld, voxel_size, ρ_calcium, VolumeFraction())

# ╔═╡ 63fc39c5-16c1-4b42-99de-feb449776beb
md"""
# Score Small Inserts
"""

# ╔═╡ 3a360b3e-a7fb-4e49-addb-0ae4e88059e6
function dilate_mask_small(mask)
	return (mask)
end

# ╔═╡ ad0a8800-79b4-42c8-9bad-960d47d438b7
function ring_mask_small(dilated_mask)
	return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask);
end

# ╔═╡ 7cccce6c-954f-4720-9dd8-31dcf442208d
md"""
## High Density
"""

# ╔═╡ b770a618-b9df-4e22-9845-defdfe725f4f
begin
    mask_S_HD_3D = Array{Bool}(undef, size(arr))
    for z in 1:size(arr, 3)
        mask_S_HD_3D[:, :, z] = mask_S_HD
    end
end;

# ╔═╡ b524f8b8-64ed-4e9d-b1f5-2e0bd49bdaa4
md"""
#### Dilated mask
"""

# ╔═╡ 34589e91-e18e-4085-95b3-56a2c9f7c54b
# dilated_mask_S_HD = dilate(dilate(mask_S_HD_3D));

# ╔═╡ b30f7b10-aba4-44e7-aaa1-61aeae3dcf48
# ring_mask_S_HD =
#     Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask_S_HD)))))) - dilated_mask_S_HD);

# ╔═╡ 7fd762ed-84f0-4a0b-9dca-6f46459469e2
begin
	dilated_mask_S_HD = dilate_mask_small(mask_S_HD_3D)
	ring_mask_S_HD = ring_mask_small(dilated_mask_S_HD)
end;

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
begin
    single_ring_mask_S_HD = Bool.(ring_mask_S_HD[:, :, 3])
    hu_heart_tissue_small_hd = mean(arr[ring_mask_S_HD])
end

# ╔═╡ 2df030b1-87fd-4eaa-923d-56035a9a825f
mass_small_hd = score(arr[dilated_mask_S_HD], hu_calcium, hu_heart_tissue_large_hd, voxel_size, ρ_calcium, VolumeFraction())

# ╔═╡ 6bdf5633-319f-4e94-b69e-469cccdd017a
md"""
## Medium Density
"""

# ╔═╡ 215a621d-bc9a-4a16-bccf-c21d4d79bc5e
begin
    mask_S_MD_3D = Array{Bool}(undef, size(arr))
    for z in 1:size(arr, 3)
        mask_S_MD_3D[:, :, z] = mask_S_MD
    end
end;

# ╔═╡ c7c1b57e-9f03-4cb6-b141-d4bf15d08905
md"""
#### Dilated mask
"""

# ╔═╡ a41f1c13-45f0-4829-a368-0b9c3853cfff
# dilated_mask_S_MD = dilate(dilate(mask_S_MD_3D));

# ╔═╡ c0a6c308-39dd-419a-9d83-0eea2121aedc
# ring_mask_S_MD =
#     Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask_S_MD)))))) - dilated_mask_S_MD);

# ╔═╡ a267bdc5-a5be-4a3d-94de-26178f4dcc17
begin
	dilated_mask_S_MD = dilate_mask_small(mask_S_MD_3D)
	ring_mask_S_MD = ring_mask_small(dilated_mask_S_MD)
end;

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
begin
    single_ring_mask_S_MD = Bool.(ring_mask_S_MD[:, :, 3])
    hu_heart_tissue_small_md = mean(arr[ring_mask_S_MD])
end

# ╔═╡ f5924025-4673-4218-a222-630584a52dc6
mass_small_md = score(arr[dilated_mask_S_MD], hu_calcium, hu_heart_tissue_large_md, voxel_size, ρ_calcium, VolumeFraction())

# ╔═╡ 3afec334-cc16-4f3d-a068-926a31e21ee5
md"""
## Low Density
"""

# ╔═╡ 6fae9a93-b88c-49c6-af36-0ee19c770748
begin
    mask_S_LD_3D = Array{Bool}(undef, size(arr))
    for z in 1:size(arr, 3)
        mask_S_LD_3D[:, :, z] = mask_S_LD
    end
end;

# ╔═╡ d50286ee-6dca-468e-9000-9a32620703ec
md"""
#### Dilated mask
"""

# ╔═╡ ead8bbbd-847e-4411-b532-7c44acdfb389
# dilated_mask_S_LD = dilate(dilate(mask_S_LD_3D));

# ╔═╡ 8a6fccee-fd88-4631-ada6-01bbcf48a506
# ring_mask_S_LD =
#     Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask_S_LD)))))) - dilated_mask_S_LD);

# ╔═╡ ac01653c-c236-44bc-9fe9-a60fc6c806b3
begin
	dilated_mask_S_LD = dilate_mask_small(mask_S_LD_3D)
	ring_mask_S_LD = ring_mask_small(dilated_mask_S_LD)
end;

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
begin
    single_ring_mask_S_LD = Bool.(ring_mask_S_LD[:, :, 3])
    hu_heart_tissue_small_ld = mean(arr[ring_mask_S_LD])
end

# ╔═╡ 9e3b40f2-9a8a-43e3-854a-b1b8b3e7be00
mass_small_ld = score(arr[dilated_mask_S_LD], hu_calcium, hu_heart_tissue_large_ld, voxel_size, ρ_calcium, VolumeFraction())

# ╔═╡ 269dde91-76cf-499b-8a51-a449c7ee7c40
md"""
# Results
"""

# ╔═╡ d4b5f299-8086-4c65-b302-9cce86ba3888
if DENSITY == "low"
	density_array = [0.025, 0.050, 0.100]
elseif DENSITY == "normal"
	density_array = [0.200, 0.400, 0.800]
end

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

# ╔═╡ 761cab03-39fa-4a1a-ba47-5ed85b97afbb
df

# ╔═╡ 9754b9de-0a56-42eb-a1df-cd5a46d4f40c
df

# ╔═╡ 3ff725ed-be27-496a-ba11-3734d47e32cd
df

# ╔═╡ 2b51aaa8-7363-435d-b1f8-5252e85c21f3
percent_error_large =
    (abs.(ground_truth_mass_large - calculated_mass_large) ./ ground_truth_mass_large) .*
    100

# ╔═╡ 13191733-be50-4428-a5d1-57145cea221a
percent_error_medium =
    (abs.(ground_truth_mass_medium - calculated_mass_medium) ./ ground_truth_mass_medium) .*
    100

# ╔═╡ 2991a411-c270-4425-9cf6-1a3b4f759668
percent_error_small =
    (abs.(ground_truth_mass_small - calculated_mass_small) ./ ground_truth_mass_small) .*
    100

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
# ╠═3f9e21d6-a663-49ce-85a3-e5943b0be245
# ╠═70fea9a9-8d77-40c4-84ea-b084d72b1b4c
# ╠═e3a4316d-481c-4e2c-a592-03727ae494e8
# ╟─fea757e2-e0c6-48eb-b997-cbabd91399bf
# ╠═3e211eda-2d67-45f9-a7ee-1825bd746cbb
# ╟─c62f5af3-a3d6-458b-bc72-3ac0f3372323
# ╠═e35a7c50-0098-47ec-bad0-84e4c36011e8
# ╠═f858d9aa-5282-401c-a5ec-d8c0a0b43238
# ╠═5ab916de-71b9-4b73-9c3e-57ea89b608bb
# ╠═c2920879-e886-45bd-814b-96aba0af7f97
# ╟─6c5111f0-de83-4e64-800c-2f9aefb5ee4f
# ╟─3b8ee7a8-917c-4836-967b-ee7e51ad0331
# ╟─67fd63c4-bac1-4ebc-9299-d17850e35d6d
# ╟─899a8652-c961-4e65-93ec-d75f0a082e87
# ╟─5cbb0bbb-96b9-4a2f-80d8-f582d4f822e0
# ╠═5b663e10-4498-4f11-aaf9-7c49b25b9fe1
# ╟─05ae002a-7a61-473e-aa06-9ebfcc51f623
# ╠═8018fb0b-5d1a-4f50-904f-b2686c75993f
# ╠═f20119be-1939-48a7-aa8d-507003b4ebb3
# ╠═28f2e565-fde8-40d7-a5fd-9b944e77efb2
# ╠═fcf12109-023a-4c6c-82a2-64471bf8ac1b
# ╠═5c3e4fe5-1ad8-4cff-8382-9fda0a05d168
# ╠═b925bdde-25e5-4560-bfef-0d822011f781
# ╠═e836390f-1647-4681-87cc-242b0c48932e
# ╟─12973897-c65a-496c-bb8a-af665c14f2a7
# ╠═099c6fdd-0f98-4d6c-9121-33c0ef0e84c5
# ╠═a85ae81e-daae-4b60-8443-fe5d53309adb
# ╟─97e09a38-569c-4cab-85dc-15f05cc3e0ad
# ╠═7fad421c-6a53-41a4-bfaa-59c25c4d4d09
# ╟─59e532c9-1550-4c20-82ea-a9df5d5dc4a5
# ╠═f5cfbc09-2455-4c5f-84bc-d704f54fee69
# ╠═542d221b-d765-4bf7-b6ca-5d3a2f5c2e52
# ╠═d7d2fadd-6bca-4132-9093-27ef4502722b
# ╟─852eda80-17c6-42c6-a015-032afaadc44e
# ╠═9e3eff5b-6bc3-4799-aea8-311709caa170
# ╠═60b42dec-ccfe-49e8-9cb0-c73ce6455b00
# ╠═88fa07e8-430b-42ce-aa09-3bdc9a79a872
# ╠═b3a43484-00db-4837-845d-620ed4c96ff1
# ╠═b4564806-2dd4-48fe-9cbf-39010fecb559
# ╠═5c981f22-89d9-436d-a94a-581250ea1ec2
# ╠═1ec199c4-6b17-4c58-b959-4904c320244d
# ╠═da3bdc26-aeb3-4a92-86dc-5c3653bc41cd
# ╠═a99e68e4-9c0d-4f0d-a164-8ab67d5839ef
# ╟─1dbad5ed-2493-4e39-8ec9-1e5e3720ddf7
# ╠═839bf6ad-3a75-4134-b15d-b612cf7e3ff8
# ╠═2bec75b1-c15b-42ce-a1c8-42af149d491e
# ╠═ba48ec9f-8970-4104-85d5-6f50e2a4c16c
# ╠═8dedc517-d985-4658-b6b7-cb8322b79d03
# ╠═c393ffaa-7e87-48af-8e2d-274dcef25447
# ╠═761cab03-39fa-4a1a-ba47-5ed85b97afbb
# ╟─2b9cc53e-bcd3-407b-bf90-c65947b28581
# ╠═77e907a4-164a-45a0-aed7-a2d599b688a4
# ╟─ddb290bf-b754-4b41-8430-963bf835ae48
# ╠═0d3b15df-f115-47e7-a15a-754025163025
# ╟─7e95a1de-6021-48ea-bd77-3480abf7b92d
# ╠═00998177-1039-40ce-ada9-1dcdf611ac98
# ╟─3dc58caa-68d3-4af5-96b4-0d02c7116a00
# ╟─4833bb01-1b20-41a5-a5b9-023b67f41c54
# ╠═62b6a7ef-59e9-4ec3-adc6-f9992d9ebfce
# ╠═54de1790-5194-490c-b0fd-e18692dd933a
# ╠═bf788e3f-538e-425e-a06a-6e4499259fee
# ╟─1b0d7eeb-d6dc-4c38-b5bb-624b54235b4b
# ╠═188742b7-dc79-496c-9c2c-948937a7ced0
# ╟─c1fea4b2-7740-4d62-b908-d89261d25064
# ╠═1a91c447-08cc-49af-ade1-b52b8a8258f0
# ╠═4c4ce42b-4a64-4f81-96e2-6da91186b62f
# ╠═968ff813-9218-4c43-946e-543301c0278e
# ╟─45900f95-7501-4655-aa3b-2c4979ae140a
# ╠═43a80610-f2ad-44a2-a2e2-33717211b29e
# ╟─57275d3f-1707-408e-96d3-bc1b31e61e96
# ╟─d667363c-ee29-47e6-8737-23a7ec0b3de3
# ╠═319c8a0e-4365-41b7-8a4a-95a436f62ff7
# ╠═2269784e-540c-4dcd-94f4-8a53285b65f0
# ╠═42a5aa69-8d08-4eba-a952-00650bad6905
# ╟─39215a15-ca20-4eb4-ac7c-07bb022e410b
# ╠═442282b4-4be9-481e-a341-f94637e71c9a
# ╟─bb4670bf-0c80-4039-b610-2db7103efd1e
# ╠═50c43516-f0e6-469b-9dd6-c9be70e12a7b
# ╠═42ba3ea2-c34c-4b5a-8b13-5dea4981c830
# ╠═bfd37f8f-987e-445f-adbb-e036bf9f2839
# ╟─50214c11-1a09-4137-b867-525327ce45d5
# ╠═1f1a5df3-8e8c-49f1-b9e4-ab966340725f
# ╟─49ef1452-3219-4f4d-aa45-42bd291ecc31
# ╟─29b88c85-4255-4080-995c-30bd7be7ddec
# ╠═2e7566eb-e208-472d-9cbf-4b528fe74a25
# ╠═c1ff9bcc-4485-4e41-998a-29e5772add07
# ╠═c70bc1a1-7ffe-43bb-8955-874d6fc3550d
# ╟─a2577753-a937-47f8-b14f-8c2309940bae
# ╠═387eeef3-a7fb-4a0c-a0f7-9fc2c30c170e
# ╠═718ec3cb-3fed-49d2-bf49-0303a4514ee6
# ╠═9754b9de-0a56-42eb-a1df-cd5a46d4f40c
# ╟─e9f53bce-0d59-456f-8add-ec3e7e0aa47e
# ╠═179cb850-1365-4b5f-87a6-c0b3432d4708
# ╟─4a65c23f-942c-45e0-ae58-0712240e9668
# ╠═b24be6e5-d41a-465c-be01-0c274c7e2124
# ╠═206e570f-5a65-43b1-b999-86e75b9cdf9d
# ╠═1d4fb6c8-225a-467c-ae32-8c8180796f2b
# ╟─3a970363-9789-440b-bc47-3d4dbf9d3fd3
# ╠═6b560d89-2ee0-4b26-88c3-a2554f01bc83
# ╟─6cc2e7d6-ca60-4f44-b16b-08c8d6a80bb8
# ╟─0d532334-ddc4-4cf3-9453-d1a0e65e7a60
# ╠═a77d1304-b887-4d77-8b3a-350527549c23
# ╠═d0fd3ef5-8e1c-4c75-bcef-94bfcb52f69b
# ╠═93c6c475-cbae-442c-aceb-0b7127b268dd
# ╟─945a529f-c24e-4122-886f-1402672ed625
# ╠═37261e95-caef-4b6e-b514-783fd0c72dfa
# ╟─040a62b8-7c41-4a39-8ecd-0a6242d3ff5b
# ╠═10dfc8bf-cb1c-4907-a613-dd6d1c1782ab
# ╠═6785b8f9-e348-4cbb-9329-f8bf3d3e89cd
# ╠═1e7995fd-29b3-4bf9-9abf-e534e95d8375
# ╟─e57598a1-35c6-42fa-82d4-2c976fccb629
# ╠═f2652dbd-28a3-406a-9bb3-8ba226780aa0
# ╟─b8be5890-a0fc-49ca-85e3-915ff81979ec
# ╟─a20af150-7341-47f0-9d1a-d5db8e8abc09
# ╠═81f4f874-f309-4764-b492-b2912108fcce
# ╠═352ed759-eb3b-4abb-8886-2f20d53bbfb9
# ╠═2a152532-109a-461d-8eec-3a37ed236161
# ╟─3632354d-4850-45d6-a9d1-15949f8c0eb4
# ╠═2b3d1930-1251-4435-ac36-09fe19b266d4
# ╟─62b7d043-47fe-4aa2-8936-2739113be27d
# ╠═04d1cc9d-9a57-47d4-902f-ed3787ef1cd9
# ╠═6f77adfd-6577-4760-ab87-7078bc4cff99
# ╠═d27c155a-1abd-4fb2-a229-a7110827ffda
# ╟─74722206-da16-4772-8ad3-b351bf0dd174
# ╠═49d8dc4d-acac-42be-b8c5-44160d42cec6
# ╟─2e8ad4d6-4d1d-41c4-9e22-9049a03d374f
# ╟─4391a50d-0642-4804-bf0a-0695ff103b28
# ╠═dad32de5-4c26-4d93-996a-9a2111d39ada
# ╠═e70b9be1-711b-41b0-9888-c31b3a2dcba4
# ╠═ec9dfcd9-13eb-41a7-be8c-405483d84a92
# ╟─63fc39c5-16c1-4b42-99de-feb449776beb
# ╠═3a360b3e-a7fb-4e49-addb-0ae4e88059e6
# ╠═ad0a8800-79b4-42c8-9bad-960d47d438b7
# ╠═3ff725ed-be27-496a-ba11-3734d47e32cd
# ╟─7cccce6c-954f-4720-9dd8-31dcf442208d
# ╠═b770a618-b9df-4e22-9845-defdfe725f4f
# ╟─b524f8b8-64ed-4e9d-b1f5-2e0bd49bdaa4
# ╠═34589e91-e18e-4085-95b3-56a2c9f7c54b
# ╠═b30f7b10-aba4-44e7-aaa1-61aeae3dcf48
# ╠═7fd762ed-84f0-4a0b-9dca-6f46459469e2
# ╟─a44c104b-577e-4746-ad94-53e1405a442f
# ╠═8aa91548-24bb-460f-a1eb-072938b96641
# ╟─09f0a38f-64c5-4795-96ba-34e96883f595
# ╟─73d06ecc-baa2-4b5e-9411-a2f8c253efaa
# ╠═be7ece4d-bf25-4c9d-baaa-1aebb696c4e2
# ╠═1aabeba7-5619-4774-96dc-eaab08b5b14a
# ╠═2df030b1-87fd-4eaa-923d-56035a9a825f
# ╟─6bdf5633-319f-4e94-b69e-469cccdd017a
# ╠═215a621d-bc9a-4a16-bccf-c21d4d79bc5e
# ╟─c7c1b57e-9f03-4cb6-b141-d4bf15d08905
# ╠═a41f1c13-45f0-4829-a368-0b9c3853cfff
# ╠═c0a6c308-39dd-419a-9d83-0eea2121aedc
# ╠═a267bdc5-a5be-4a3d-94de-26178f4dcc17
# ╟─c4d7421e-3046-431f-aabc-cab0f88de453
# ╠═5ebb2e49-79e4-462d-92a6-6c3224c64e36
# ╟─e2f2a55f-9934-45fd-9541-4d30de449cfe
# ╟─2fa9406a-0333-45a6-a674-5a9942aa85d7
# ╠═60293d04-c27a-4c17-8b96-59f9d42ba781
# ╠═7622bbd6-b851-4f59-836d-0f53ca8c610f
# ╠═f5924025-4673-4218-a222-630584a52dc6
# ╟─3afec334-cc16-4f3d-a068-926a31e21ee5
# ╠═6fae9a93-b88c-49c6-af36-0ee19c770748
# ╟─d50286ee-6dca-468e-9000-9a32620703ec
# ╠═ead8bbbd-847e-4411-b532-7c44acdfb389
# ╠═8a6fccee-fd88-4631-ada6-01bbcf48a506
# ╠═ac01653c-c236-44bc-9fe9-a60fc6c806b3
# ╟─8c7021cd-d89f-422d-9377-a5fdba5759ab
# ╠═122be89d-bffa-418b-8449-b273bf73ec67
# ╟─486207ce-fd95-410f-9286-75d7b2f21d32
# ╟─45be9e86-390c-4ab7-8d01-c2a2ebb1e245
# ╠═355fd701-c652-4421-a906-1eb1f7506cc4
# ╠═6e8e836b-5069-4de7-abee-34f7bf882a13
# ╠═9e3b40f2-9a8a-43e3-854a-b1b8b3e7be00
# ╟─269dde91-76cf-499b-8a51-a449c7ee7c40
# ╠═d4b5f299-8086-4c65-b302-9cce86ba3888
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

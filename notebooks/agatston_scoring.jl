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

# ╔═╡ f5acf58b-8436-4b0b-a170-746ba5475c10
begin
	let
		using Pkg
		Pkg.activate(mktempdir())
		Pkg.Registry.update()
		Pkg.add("PlutoUI")
		Pkg.add("CairoMakie")
		Pkg.add("Statistics")
		Pkg.add("Images")
		Pkg.add("ImageMorphology")
		Pkg.add("ImageFiltering")
		Pkg.add("CSV")
		Pkg.add("DataFrames")
		Pkg.add(url="https://github.com/JuliaHealth/DICOM.jl")
		Pkg.add(url="https://github.com/Dale-Black/DICOMUtils.jl")
		Pkg.add(url="https://github.com/Dale-Black/Phantoms.jl")
		Pkg.add(url="https://github.com/Dale-Black/CalciumScoring.jl")
	end
	
	using PlutoUI
	using CairoMakie
	using Statistics
	using Images
	using ImageMorphology
	using ImageFiltering
	using CSV
	using DataFrames
	using DICOM
	using DICOMUtils
	using Phantoms
	using CalciumScoring
end

# ╔═╡ 8219f406-7175-4573-ae37-642ef9c45b1b
TableOfContents()

# ╔═╡ ef30dc5c-0f83-485b-9718-669f292a87ec
md"""
## Load DICOMS
"""

# ╔═╡ 422ba302-e559-48c0-90a5-73a0544ccbd5
begin
	SCAN_NUMBER = 1
	VENDER = "135"
	TYPE = "agatston"
	BASE_PATH = "/Users/daleblack/Google Drive/Datasets/Simulated/"
end

# ╔═╡ 80592b38-6672-419e-9315-0e33046c42a8
md"""
**Everything below should be automatic, just scroll through to visually inspect that things make sense**
"""

# ╔═╡ d2f9de75-8b2a-40f2-98c0-5df90c32c8d6
root_path = string(BASE_PATH, VENDER)

# ╔═╡ 813bf029-e475-4e65-a2ca-fdb6c23de115
dcm_path_list = dcm_list_builder(root_path)

# ╔═╡ ca10a4df-7e26-49c3-bdee-de9f0ffdc39c
pth = dcm_path_list[SCAN_NUMBER]

# ╔═╡ 670cfd3e-c32f-4357-a2b4-9b9d31d49c6a
scan = basename(pth)

# ╔═╡ 9dd34ad4-3576-4b85-b811-cafcd05b60a9
header, dcm_array, slice_thick_ori1 = dcm_reader(pth);

# ╔═╡ 29ef15d4-754a-4619-9abc-51182709bc5e
md"""
## Helper Functions
"""

# ╔═╡ e500136f-34f4-484c-8296-088ee69b0f48
function collect_tuple(tuple_array)
	row_num = size(tuple_array)
	col_num = length(tuple_array[1])
	container = zeros(Int64, row_num..., col_num)
	for i in 1:length(tuple_array)
		container[i,:] = collect(tuple_array[i])
	end
	return container
end

# ╔═╡ f6c6aa68-5e87-495e-b45b-ca96e61e9a69
function overlay_mask_bind(mask)
	indices = findall(x -> x == 1, mask)
	indices = Tuple.(indices)
	label_array = collect_tuple(indices)
	zs = unique(label_array[:,3])
	return PlutoUI.Slider(1:length(zs), default=3, show_value=true)
end

# ╔═╡ 8a8e8a42-a823-47dd-b6eb-ef3dba4182d5
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

# ╔═╡ 22ceda3c-3c30-4632-9052-4dd90cfabb24
function show_matrix(A::Matrix, red::Union{BitMatrix,Matrix{Bool}}=zeros(Bool, size(A)))
	base = RGB.(Gray.(A))

	base[red] .= RGB(1.0, 0.1, 0.1)

	# some tricks to show the pixels more clearly:
	s = max(size(A)...)
	if s >= 20
		min_size = 1200
		factor = min(5,min_size ÷ s)
		
		kron(base, ones(factor, factor))
	else
		base
	end
end

# ╔═╡ 4f486b1a-b272-4b7d-914c-baa9530ff4e3
function create_mask(array, mask)
	@assert size(array) == size(mask)
	idxs = findall(x -> x == true, mask)
	overlayed_mask = zeros(size(array))
	for idx in idxs
		overlayed_mask[idx] = array[idx]
	end
	return overlayed_mask
end

# ╔═╡ 8151a4eb-4139-48d2-8a25-14c350b30d6d
md"""
## Segment Heart
"""

# ╔═╡ 3c5e62f5-a435-43c4-a71e-256e3e0caacb
masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3)÷2);

# ╔═╡ 14f7ed28-00ad-4b4e-be1f-b85f50c39b23
center_insert

# ╔═╡ 88f63f47-ac6b-4670-938b-9eda91550c8c
@bind a PlutoUI.Slider(1:size(masked_array, 3), default=10, show_value=true)

# ╔═╡ 06e6ec7d-e2e7-4348-b38c-8fde60a24807
heatmap(masked_array[:, :, a], colormap=:grays)

# ╔═╡ 62e0db89-5a32-40d8-9264-684f7df58476
begin
	fig = Figure()
	
	ax = Makie.Axis(fig[1, 1])
	ax.title = "Raw DICOM Array"
	heatmap!(transpose(dcm_array[:, :, 4]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig
end

# ╔═╡ 09bd19ff-0db5-4a79-adb8-c3e202a4d7c7
begin
	fig2 = Figure()
	
	ax2 = Makie.Axis(fig2[1, 1])
	ax2.title = "Mask Array"
	heatmap!(transpose(mask), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig2
end

# ╔═╡ 85a5667d-6606-4cce-8ca2-b3682422648e
begin
	fig3 = Figure()
	
	ax3 = Makie.Axis(fig3[1, 1])
	ax3.title = "Masked DICOM Array"
	heatmap!(transpose(masked_array[:, :, 1]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig3
end

# ╔═╡ 935767d5-046d-46b0-b1cc-01f85b0698f5
md"""
## Segment Calcium Rod
"""

# ╔═╡ d1cc731f-7a3d-412e-ae6c-a4cc9c4d1e9d
calcium_image, slice_CCI, quality_slice, cal_rod_slice = mask_rod(masked_array, header);

# ╔═╡ a15bd533-776e-44a7-bac4-6c37478d5969
@bind c PlutoUI.Slider(1:size(calcium_image, 3), default=cal_rod_slice, show_value=true)

# ╔═╡ 1665243d-4c39-4c69-bf35-bf30c4662c48
heatmap(transpose(calcium_image[:, :, c]), colormap=:grays)

# ╔═╡ 638286bd-35c6-4d18-84a2-15e92d2e815f
md"""
## Segment Calcium Inserts
"""

# ╔═╡ 77ff2fb9-d20c-4ec8-95a2-e95a677eadda
mask_L_HD, mask_M_HD, mask_S_HD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_LD, mask_M_LD, mask_S_LD = mask_inserts_simulation(
            dcm_array, masked_array, header, slice_CCI, center_insert
);

# ╔═╡ bfd0e2a3-ef9a-4086-af43-1a1ebe552258
slice_CCI

# ╔═╡ 84938868-94a9-4fd0-8164-209122612f88
masks = mask_L_HD + mask_M_HD + mask_S_HD + mask_L_MD + mask_M_MD + mask_S_MD + mask_L_LD + mask_M_LD + mask_S_LD;

# ╔═╡ 229478be-eaac-4750-a3fd-790bfd9d541c
heatmap(masks, colormap=:grays)

# ╔═╡ c73cd156-dde4-446c-a4e2-b3c04283af4f
md"""
# Agatston Scoring
"""

# ╔═╡ a56d63f4-9f6b-48c4-8bf2-35f76acbf1d2
md"""
## High Density
"""

# ╔═╡ e12e4bb7-4ecd-4a7f-9ebd-70a80a18ac05
arr = masked_array[:, :, 4:6];

# ╔═╡ 03a30d17-76ae-4551-9d9f-5c5c22aaec73
begin
	mask_L_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_HD_3D[:, :, z] = dilate(dilate(mask_L_HD))
	end
end;

# ╔═╡ 9c68fb4e-d92d-4cb3-b04a-c3586212af4d
md"""
#### Dilated mask
"""

# ╔═╡ 12c89924-fc11-4696-b367-45457c5f255a
dilated_mask_L_HD = dilate(dilate(mask_L_HD_3D));

# ╔═╡ 262e4339-e041-4235-bc94-e6aae01d74a3
@bind g2 overlay_mask_bind(dilated_mask_L_HD)

# ╔═╡ d00f4649-e221-43a9-a504-5f61722eea91
overlay_mask_plot(arr, dilated_mask_L_HD, g2, "dilated mask")

# ╔═╡ a945ece4-2bd8-40db-bbe7-9747efa087e6
pixel_size = DICOMUtils.get_pixel_size(header)

# ╔═╡ 2d2d2600-dd70-4325-9032-4467882aff73
overlayed_mask_l_hd = create_mask(arr, dilated_mask_L_HD);

# ╔═╡ 73121968-a17a-4dcd-bbe5-96c0a04b5cb9
alg = Agatston()

# ╔═╡ 59dfd023-2fb8-4366-b223-1153e198efb9
agat_l_hd = score(overlayed_mask_l_hd, pixel_size, alg)

# ╔═╡ d390e652-23a9-455d-87a0-3f8c54e2bb09
md"""
## Medium Density
"""

# ╔═╡ fcc751cd-120c-4c6f-828d-cf070b3466fa
begin
	mask_L_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_MD_3D[:, :, z] = mask_L_MD
	end
end;

# ╔═╡ 96d09c06-e499-4390-b101-d663af7e5d75
md"""
#### Dilated mask
"""

# ╔═╡ fd344c0d-eb41-488c-9c6a-ef68cbd68637
dilated_mask_L_MD = dilate(dilate(mask_L_MD_3D));

# ╔═╡ 24fe0716-4e2f-47fb-b438-42d264bce384
@bind h2 overlay_mask_bind(dilated_mask_L_MD)

# ╔═╡ 37d44541-d30e-42be-b2d3-3423ac34c7dc
overlay_mask_plot(arr, dilated_mask_L_MD, h2, "dilated mask")

# ╔═╡ 8270b725-8748-4f48-8ac2-20ea81b80306
overlayed_mask_l_md = create_mask(arr, dilated_mask_L_MD);

# ╔═╡ acdf2d47-adfa-4ee0-9c9b-4d2bdb7e7aca
agat_l_md = score(overlayed_mask_l_md, pixel_size, alg)

# ╔═╡ d5f66a85-966a-4cfc-8e8c-c2211e103f74
md"""
## Low Density
"""

# ╔═╡ 5a84502c-39ba-4234-a1b2-f852c7955193
begin
	mask_L_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_LD_3D[:, :, z] = mask_L_LD
	end
end;

# ╔═╡ 53a04c71-9bcf-4af7-96de-71322df085b4
md"""
#### Dilated mask
"""

# ╔═╡ 473b7f74-904b-48b8-9a9f-d60eae5bbedd
dilated_mask_L_LD = dilate(dilate(mask_L_LD_3D));

# ╔═╡ 7237e02e-730e-4b6b-b7ef-60c8f19d5b3a
@bind i2 overlay_mask_bind(dilated_mask_L_LD)

# ╔═╡ bbd8b972-a2b7-4e2e-ae1f-e2d7a7c99a7b
overlay_mask_plot(arr, dilated_mask_L_LD, i2, "dilated mask")

# ╔═╡ 7c26c32f-d77d-41c4-b461-a2497549ecfb
overlayed_mask_l_ld = create_mask(arr, dilated_mask_L_LD);

# ╔═╡ ed7d7a23-3574-4cf2-9295-88386c013270
agat_l_ld = score(overlayed_mask_l_ld, pixel_size, alg)

# ╔═╡ 73ac04c0-c8c8-4d55-aa6c-7e1747b951d2
md"""
# Score Medium Inserts
"""

# ╔═╡ 569182e7-e7fa-4063-9ba4-0e84ac7d5202
md"""
## High Density
"""

# ╔═╡ 4d0872f8-b728-4bfe-a96c-dfb53710a571
begin
	mask_M_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_HD_3D[:, :, z] = mask_M_HD
	end
end;

# ╔═╡ d4e77ed5-d2c7-4288-a35d-34812617d089
md"""
#### Dilated mask
"""

# ╔═╡ 10546afe-656c-4837-9f21-f144023d56d4
dilated_mask_M_HD = dilate(dilate(dilate(dilate(mask_M_HD_3D))));

# ╔═╡ d1c3cd8a-258b-426f-bea0-aad31ad825f3
@bind j2 overlay_mask_bind(dilated_mask_M_HD)

# ╔═╡ 744681a4-70e5-49d9-bec1-29c4ecc80aea
overlay_mask_plot(arr, dilated_mask_M_HD, j2, "dilated mask")

# ╔═╡ fcffc8e9-2dca-425c-aad7-3a0be219b67c
overlayed_mask_m_hd = create_mask(arr, dilated_mask_M_HD);

# ╔═╡ 08141335-4566-4ba5-a29b-a39aa428b995
agat_m_hd = score(overlayed_mask_m_hd, pixel_size, alg)

# ╔═╡ c34ea971-4638-4348-8c7e-0342490421b1
md"""
## Medium Density
"""

# ╔═╡ 4b1a3793-fbc1-4d66-83c9-0dd8424c9334
begin
	mask_M_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_MD_3D[:, :, z] = mask_M_MD
	end
end;

# ╔═╡ f4d289cf-25bf-4024-9be6-4397e03af0d3
md"""
#### Dilated mask
"""

# ╔═╡ 6dc948d0-d637-41bb-8c8e-6721d06f38e9
dilated_mask_M_MD = dilate(dilate(dilate(dilate(mask_M_MD_3D))));

# ╔═╡ a0f4c8df-2e3b-4de8-87f0-a180ee9a4ce1
@bind k2 overlay_mask_bind(dilated_mask_M_MD)

# ╔═╡ 4389d5d2-1a5c-4f96-8f03-d8824cd4381f
overlay_mask_plot(arr, dilated_mask_M_MD, k2, "dilated mask")

# ╔═╡ d656924f-72a8-4d91-8119-aa8c63ded3b8
overlayed_mask_m_md = create_mask(arr, dilated_mask_M_MD);

# ╔═╡ 732f75b0-4967-46e2-9777-21b6b4fd7990
agat_m_md = score(overlayed_mask_m_md, pixel_size, alg)

# ╔═╡ 25743bee-e433-4b4d-91b9-27c9a24a3958
md"""
## Low Density
"""

# ╔═╡ f30257ed-4335-4ba6-9b91-93db7fd4ed9f
begin
	mask_M_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_LD_3D[:, :, z] = mask_M_LD
	end
end;

# ╔═╡ 451821dc-8ff2-4042-af14-e52b94c298cb
md"""
#### Dilated mask
"""

# ╔═╡ 036386c8-591c-4478-bde1-33d973c458b9
dilated_mask_M_LD = dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))));

# ╔═╡ 3fa6ba73-74c8-473a-8cd2-1cfe77381882
@bind l2 overlay_mask_bind(dilated_mask_M_LD)

# ╔═╡ 14cc065e-5150-4cad-a9bb-a562c0f4a9ee
overlay_mask_plot(arr, dilated_mask_M_LD, l2, "dilated mask")

# ╔═╡ fe2502fe-3ed6-4c16-94ad-6aafd87d1af9
overlayed_mask_m_ld = create_mask(arr, dilated_mask_M_LD);

# ╔═╡ 40d01237-f9fc-458e-adb1-4d3476e3a55c
agat_m_ld = score(overlayed_mask_m_ld, pixel_size, alg)

# ╔═╡ a2b6908e-06a8-465c-931c-2c7895eb7fdd
md"""
# Score Small Inserts
"""

# ╔═╡ f30a694c-a9a5-4096-b03c-20e018f344b4
md"""
## High Density
"""

# ╔═╡ 324f5030-2bb5-4c2d-87bf-21069f86f898
begin
	mask_S_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_HD_3D[:, :, z] = mask_S_HD
	end
end;

# ╔═╡ e2e5d32e-beab-470f-b753-c63bf58586ac
md"""
#### Dilated mask
"""

# ╔═╡ e0ae05d4-c910-453b-9851-445d36c715ae
dilated_mask_S_HD = dilate(dilate(dilate(dilate(dilate((mask_S_HD_3D))))));

# ╔═╡ 27d2277e-a4f5-46d5-b093-c41462622304
@bind m2 overlay_mask_bind(dilated_mask_S_HD)

# ╔═╡ 8637a05f-8db4-424e-a453-331bb564f181
overlay_mask_plot(arr, dilated_mask_S_HD, m2, "dilated mask")

# ╔═╡ 93158aaa-6622-4fa8-babc-959e7349ccba
overlayed_mask_s_hd = create_mask(arr, dilated_mask_S_HD);

# ╔═╡ 787547f0-cd60-43dd-8639-f9909d9197f3
agat_s_hd = score(overlayed_mask_s_hd, pixel_size, alg)

# ╔═╡ ad7c62c8-52a8-4877-8b68-5be216f709cd
md"""
## Medium Density
"""

# ╔═╡ b72035ba-89ff-4405-aefc-b6af912d7af2
begin
	mask_S_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_MD_3D[:, :, z] = mask_S_MD
	end
end;

# ╔═╡ 754d327d-2148-407f-aa8d-6cdf78bc55c3
md"""
#### Dilated mask
"""

# ╔═╡ bbb2d7ca-a25e-48bc-af10-af571d06009f
dilated_mask_S_MD = dilate(dilate(dilate(dilate(dilate(mask_S_MD_3D)))));

# ╔═╡ c311b2bc-0dd2-49a6-a295-cbecb44dbcd6
@bind n2 overlay_mask_bind(dilated_mask_S_MD)

# ╔═╡ 312bb63d-a933-462c-a513-5d20f905e950
overlay_mask_plot(arr, dilated_mask_S_MD, n2, "dilated mask")

# ╔═╡ 352a5ab8-b337-4e8a-bb44-8f977ec0675b
overlayed_mask_s_md = create_mask(arr, dilated_mask_S_MD);

# ╔═╡ f87fa716-e714-47fb-8e32-f7ed10c31999
agat_s_md = score(overlayed_mask_s_md, pixel_size, alg)

# ╔═╡ e6d4b895-03c7-4913-8aaf-1af09260a603
md"""
## Low Density
"""

# ╔═╡ 15747d6b-9b96-4300-97cf-b5ec4550be90
begin
	mask_S_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_LD_3D[:, :, z] = mask_S_LD
	end
end;

# ╔═╡ 0bd2712f-f02d-4287-a4b5-b7e626154234
md"""
#### Dilated mask
"""

# ╔═╡ 07e119cf-5443-44b4-8b24-4acf59803cae
dilated_mask_S_LD = dilate(dilate(dilate(dilate(dilate(mask_S_LD_3D)))));

# ╔═╡ c532d1b4-f852-4f49-882c-b9eab4ac721a
@bind o2 overlay_mask_bind(dilated_mask_S_LD)

# ╔═╡ 9698c301-59e0-4ca2-b3c3-3b989149558e
overlay_mask_plot(arr, dilated_mask_S_LD, o2, "dilated mask")

# ╔═╡ b48fbcd8-9c0b-4197-a6fd-2302ed0f644e
overlayed_mask_s_ld = create_mask(arr, dilated_mask_S_LD);

# ╔═╡ 5c9e783f-46d6-416d-a6ef-bb8b019b13ac
agat_s_ld = score(overlayed_mask_s_ld, pixel_size, alg)

# ╔═╡ 8f10515b-eb71-4557-a97b-cdb38878626a
md"""
# Results
"""

# ╔═╡ 2aaa3421-55a9-4424-b89d-ce627b303a27
density_array = [0, 200, 400, 800]

# ╔═╡ ffc3809f-798d-4fa7-ab14-1e80a949e8aa
inserts = [
	"Low Density",
	"Medium Density",
	"High Density"
]

# ╔═╡ c3540d63-ef20-410e-b5bc-930f5567b534
stanford_agat_large = [
	138.60,
	238.38,
	402.94
] # mg

# ╔═╡ 6ae8d140-06f5-4925-a1ca-ab994a4f6229
calculated_agat_large = [
	agat_l_ld,
	agat_l_md,
	agat_l_hd
]

# ╔═╡ 39e1efd3-3b8c-4c2b-a6cf-3f3a135303bb
stanford_agat_medium = [
	22.86,
	69.54,
	88.11
]

# ╔═╡ 45a79998-9dce-4b99-a54c-ff02ddea5420
calculated_agat_medium = [
	agat_m_ld,
	agat_m_md,
	agat_m_hd
]

# ╔═╡ 6fd9d10b-09b5-46ad-aaa2-0e9dcc446438
stanford_agat_small = [
	0.0,
	0.0,
	6.19
]

# ╔═╡ 0e1e11aa-9e6f-4aad-82f6-30402dd21603
calculated_agat_small = [
	agat_s_ld,
	agat_s_md,
	agat_s_hd
]

# ╔═╡ 76ebbe09-7075-4b4d-a76c-02982ac219d4
df = DataFrame(
	inserts = inserts,
	stanford_agat_large = stanford_agat_large,
	calculated_agat_large = calculated_agat_large,
	stanford_agat_medium = stanford_agat_medium,
	calculated_agat_medium = calculated_agat_medium,
	stanford_agat_small = stanford_agat_small,
	calculated_agat_small = calculated_agat_small
)

# ╔═╡ bbd00738-7632-4b0d-914a-d25aec8bb05f
begin
	fmass2 = Figure()
	axmass2 = Makie.Axis(fmass2[1, 1])
	
	scatter!(density_array[2:end], df[!, :stanford_agat_large], label="stanford_agat_large")
	scatter!(density_array[2:end], df[!, :calculated_agat_large], label="calculated_agat_large")
	
	axmass2.title = "Agatstom Scores (Large)"
	axmass2.ylabel = "Agatston Score (mg)"
	axmass2.xlabel = "Density (mg/cm^3)"

	xlims!(axmass2, 0, 850)
	ylims!(axmass2, 0, 450)
	
	fmass2[1, 2] = Legend(fmass2, axmass2, framevisible = false)
	
	fmass2
end

# ╔═╡ ad6aad8f-ce74-4ba6-90ed-48d50fdefd5b
begin
	fmass3 = Figure()
	axmass3 = Makie. Axis(fmass3[1, 1])
	
	scatter!(density_array[2:end], df[!, :stanford_agat_medium], label="stanford_agat_medium")
	scatter!(density_array[2:end], df[!, :calculated_agat_medium], label="calculated_agat_medium")
	
	axmass3.title = "Agatston Scores (Medium)"
	axmass3.ylabel = "Agatston Score"
	axmass3.xlabel = "Density (mg/cm^3)"

	xlims!(axmass3, 0, 850)
	ylims!(axmass3, 0, 110)
	
	fmass3[1, 2] = Legend(fmass3, axmass3, framevisible = false)
	
	fmass3
end

# ╔═╡ ccfc2e1c-3729-41d1-b50f-d99bd44bf625
begin
	fmass4 = Figure()
	axmass4 = Makie.Axis(fmass4[1, 1])
	
	scatter!(density_array[2:end], df[!, :stanford_agat_small], label="stanford_agat_small")
	scatter!(density_array[2:end], df[!, :calculated_agat_small], label="calculated_agat_small")
	
	axmass4.title = "Agatston Scores (Small)"
	axmass4.ylabel = "Atatston Score (mg)"
	axmass4.xlabel = "Density (mg/cm^3)"

	xlims!(axmass4, 0, 850)
	ylims!(axmass4, 0, 7.0)
	
	fmass4[1, 2] = Legend(fmass4, axmass4, framevisible = false)
	
	fmass4
end

# ╔═╡ 722fb5f8-7443-41a3-acd8-01359d26e03c
md"""
### Save Results
"""

# ╔═╡ 3aa596a8-141d-4ce5-9047-c54ba76bc957
if ~isdir(string(cd(pwd, "..") , "/output/", TYPE))
	mkdir(string(cd(pwd, "..") , "/output/", TYPE))
end

# ╔═╡ 963377a8-4c68-4516-be15-699a9fefcb5a
output_path = string(cd(pwd, "..") , "/output/", TYPE, "/", scan, ".csv")

# ╔═╡ fba127c6-2ffc-41ad-8450-ee7588f9e681
CSV.write(output_path, df)

# ╔═╡ Cell order:
# ╠═f5acf58b-8436-4b0b-a170-746ba5475c10
# ╠═8219f406-7175-4573-ae37-642ef9c45b1b
# ╟─ef30dc5c-0f83-485b-9718-669f292a87ec
# ╠═422ba302-e559-48c0-90a5-73a0544ccbd5
# ╟─80592b38-6672-419e-9315-0e33046c42a8
# ╠═d2f9de75-8b2a-40f2-98c0-5df90c32c8d6
# ╠═813bf029-e475-4e65-a2ca-fdb6c23de115
# ╠═ca10a4df-7e26-49c3-bdee-de9f0ffdc39c
# ╠═670cfd3e-c32f-4357-a2b4-9b9d31d49c6a
# ╠═9dd34ad4-3576-4b85-b811-cafcd05b60a9
# ╟─29ef15d4-754a-4619-9abc-51182709bc5e
# ╟─e500136f-34f4-484c-8296-088ee69b0f48
# ╟─f6c6aa68-5e87-495e-b45b-ca96e61e9a69
# ╟─8a8e8a42-a823-47dd-b6eb-ef3dba4182d5
# ╟─22ceda3c-3c30-4632-9052-4dd90cfabb24
# ╟─4f486b1a-b272-4b7d-914c-baa9530ff4e3
# ╟─8151a4eb-4139-48d2-8a25-14c350b30d6d
# ╠═3c5e62f5-a435-43c4-a71e-256e3e0caacb
# ╠═14f7ed28-00ad-4b4e-be1f-b85f50c39b23
# ╟─88f63f47-ac6b-4670-938b-9eda91550c8c
# ╠═06e6ec7d-e2e7-4348-b38c-8fde60a24807
# ╟─62e0db89-5a32-40d8-9264-684f7df58476
# ╟─09bd19ff-0db5-4a79-adb8-c3e202a4d7c7
# ╟─85a5667d-6606-4cce-8ca2-b3682422648e
# ╟─935767d5-046d-46b0-b1cc-01f85b0698f5
# ╠═d1cc731f-7a3d-412e-ae6c-a4cc9c4d1e9d
# ╠═1665243d-4c39-4c69-bf35-bf30c4662c48
# ╟─a15bd533-776e-44a7-bac4-6c37478d5969
# ╟─638286bd-35c6-4d18-84a2-15e92d2e815f
# ╠═77ff2fb9-d20c-4ec8-95a2-e95a677eadda
# ╠═bfd0e2a3-ef9a-4086-af43-1a1ebe552258
# ╠═84938868-94a9-4fd0-8164-209122612f88
# ╠═229478be-eaac-4750-a3fd-790bfd9d541c
# ╟─c73cd156-dde4-446c-a4e2-b3c04283af4f
# ╟─a56d63f4-9f6b-48c4-8bf2-35f76acbf1d2
# ╠═e12e4bb7-4ecd-4a7f-9ebd-70a80a18ac05
# ╠═03a30d17-76ae-4551-9d9f-5c5c22aaec73
# ╟─9c68fb4e-d92d-4cb3-b04a-c3586212af4d
# ╠═12c89924-fc11-4696-b367-45457c5f255a
# ╟─262e4339-e041-4235-bc94-e6aae01d74a3
# ╠═d00f4649-e221-43a9-a504-5f61722eea91
# ╠═a945ece4-2bd8-40db-bbe7-9747efa087e6
# ╠═2d2d2600-dd70-4325-9032-4467882aff73
# ╠═73121968-a17a-4dcd-bbe5-96c0a04b5cb9
# ╠═59dfd023-2fb8-4366-b223-1153e198efb9
# ╟─d390e652-23a9-455d-87a0-3f8c54e2bb09
# ╠═fcc751cd-120c-4c6f-828d-cf070b3466fa
# ╟─96d09c06-e499-4390-b101-d663af7e5d75
# ╠═fd344c0d-eb41-488c-9c6a-ef68cbd68637
# ╟─24fe0716-4e2f-47fb-b438-42d264bce384
# ╠═37d44541-d30e-42be-b2d3-3423ac34c7dc
# ╠═8270b725-8748-4f48-8ac2-20ea81b80306
# ╠═acdf2d47-adfa-4ee0-9c9b-4d2bdb7e7aca
# ╟─d5f66a85-966a-4cfc-8e8c-c2211e103f74
# ╠═5a84502c-39ba-4234-a1b2-f852c7955193
# ╟─53a04c71-9bcf-4af7-96de-71322df085b4
# ╠═473b7f74-904b-48b8-9a9f-d60eae5bbedd
# ╟─7237e02e-730e-4b6b-b7ef-60c8f19d5b3a
# ╠═bbd8b972-a2b7-4e2e-ae1f-e2d7a7c99a7b
# ╠═7c26c32f-d77d-41c4-b461-a2497549ecfb
# ╠═ed7d7a23-3574-4cf2-9295-88386c013270
# ╟─73ac04c0-c8c8-4d55-aa6c-7e1747b951d2
# ╟─569182e7-e7fa-4063-9ba4-0e84ac7d5202
# ╠═4d0872f8-b728-4bfe-a96c-dfb53710a571
# ╟─d4e77ed5-d2c7-4288-a35d-34812617d089
# ╠═10546afe-656c-4837-9f21-f144023d56d4
# ╟─d1c3cd8a-258b-426f-bea0-aad31ad825f3
# ╠═744681a4-70e5-49d9-bec1-29c4ecc80aea
# ╠═fcffc8e9-2dca-425c-aad7-3a0be219b67c
# ╠═08141335-4566-4ba5-a29b-a39aa428b995
# ╟─c34ea971-4638-4348-8c7e-0342490421b1
# ╠═4b1a3793-fbc1-4d66-83c9-0dd8424c9334
# ╟─f4d289cf-25bf-4024-9be6-4397e03af0d3
# ╠═6dc948d0-d637-41bb-8c8e-6721d06f38e9
# ╟─a0f4c8df-2e3b-4de8-87f0-a180ee9a4ce1
# ╠═4389d5d2-1a5c-4f96-8f03-d8824cd4381f
# ╠═d656924f-72a8-4d91-8119-aa8c63ded3b8
# ╠═732f75b0-4967-46e2-9777-21b6b4fd7990
# ╟─25743bee-e433-4b4d-91b9-27c9a24a3958
# ╠═f30257ed-4335-4ba6-9b91-93db7fd4ed9f
# ╟─451821dc-8ff2-4042-af14-e52b94c298cb
# ╠═036386c8-591c-4478-bde1-33d973c458b9
# ╟─3fa6ba73-74c8-473a-8cd2-1cfe77381882
# ╠═14cc065e-5150-4cad-a9bb-a562c0f4a9ee
# ╠═fe2502fe-3ed6-4c16-94ad-6aafd87d1af9
# ╠═40d01237-f9fc-458e-adb1-4d3476e3a55c
# ╟─a2b6908e-06a8-465c-931c-2c7895eb7fdd
# ╟─f30a694c-a9a5-4096-b03c-20e018f344b4
# ╠═324f5030-2bb5-4c2d-87bf-21069f86f898
# ╟─e2e5d32e-beab-470f-b753-c63bf58586ac
# ╠═e0ae05d4-c910-453b-9851-445d36c715ae
# ╟─27d2277e-a4f5-46d5-b093-c41462622304
# ╠═8637a05f-8db4-424e-a453-331bb564f181
# ╠═93158aaa-6622-4fa8-babc-959e7349ccba
# ╟─787547f0-cd60-43dd-8639-f9909d9197f3
# ╟─ad7c62c8-52a8-4877-8b68-5be216f709cd
# ╠═b72035ba-89ff-4405-aefc-b6af912d7af2
# ╟─754d327d-2148-407f-aa8d-6cdf78bc55c3
# ╠═bbb2d7ca-a25e-48bc-af10-af571d06009f
# ╟─c311b2bc-0dd2-49a6-a295-cbecb44dbcd6
# ╠═312bb63d-a933-462c-a513-5d20f905e950
# ╠═352a5ab8-b337-4e8a-bb44-8f977ec0675b
# ╠═f87fa716-e714-47fb-8e32-f7ed10c31999
# ╟─e6d4b895-03c7-4913-8aaf-1af09260a603
# ╠═15747d6b-9b96-4300-97cf-b5ec4550be90
# ╟─0bd2712f-f02d-4287-a4b5-b7e626154234
# ╠═07e119cf-5443-44b4-8b24-4acf59803cae
# ╟─c532d1b4-f852-4f49-882c-b9eab4ac721a
# ╠═9698c301-59e0-4ca2-b3c3-3b989149558e
# ╠═b48fbcd8-9c0b-4197-a6fd-2302ed0f644e
# ╠═5c9e783f-46d6-416d-a6ef-bb8b019b13ac
# ╟─8f10515b-eb71-4557-a97b-cdb38878626a
# ╠═2aaa3421-55a9-4424-b89d-ce627b303a27
# ╠═ffc3809f-798d-4fa7-ab14-1e80a949e8aa
# ╠═c3540d63-ef20-410e-b5bc-930f5567b534
# ╠═6ae8d140-06f5-4925-a1ca-ab994a4f6229
# ╠═39e1efd3-3b8c-4c2b-a6cf-3f3a135303bb
# ╠═45a79998-9dce-4b99-a54c-ff02ddea5420
# ╠═6fd9d10b-09b5-46ad-aaa2-0e9dcc446438
# ╠═0e1e11aa-9e6f-4aad-82f6-30402dd21603
# ╠═76ebbe09-7075-4b4d-a76c-02982ac219d4
# ╟─bbd00738-7632-4b0d-914a-d25aec8bb05f
# ╟─ad6aad8f-ce74-4ba6-90ed-48d50fdefd5b
# ╟─ccfc2e1c-3729-41d1-b50f-d99bd44bf625
# ╟─722fb5f8-7443-41a3-acd8-01359d26e03c
# ╠═3aa596a8-141d-4ce5-9047-c54ba76bc957
# ╠═963377a8-4c68-4516-be15-699a9fefcb5a
# ╠═fba127c6-2ffc-41ad-8450-ee7588f9e681

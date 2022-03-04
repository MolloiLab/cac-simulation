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

# ╔═╡ 763d240d-eac6-4b56-98de-ae16b32d516a
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
		Pkg.add(url="https://github.com/JuliaHealth/DICOM.jl")
		Pkg.add(url="https://github.com/Dale-Black/DICOMUtils.jl")
		Pkg.add(url="https://github.com/Dale-Black/Phantoms.jl")
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
	using DICOM
	using DICOMUtils
	using Phantoms
	using CalciumScoring
end

# ╔═╡ 3ca65cc1-c88b-4ddc-8f90-a0726aeb362c
TableOfContents()

# ╔═╡ d2e68b0d-d60f-4296-ad4e-dfd7a7ec6bb8
md"""
## Load DICOMS

All you need to do is set `base_path` once and leave it. After that, the only thing that should change is the `VENDER`, once for every set, and the `SCAN_NUMBER`, once for each scan.
"""

# ╔═╡ 514daca0-ecef-47b1-ac4c-d5c7d9042c73
begin
	SCAN_NUMBER = 1
	VENDER = "combined"
	BASE_PATH = "/Users/daleblack/Google Drive/dev/MolloiLab/CAC-stanford-data/data/simulated/"
end

# ╔═╡ 84bfd706-4cf4-4ea6-af92-1c2c3fcf9944
md"""
**Everything below should be automatic, just scroll through to visually inspect that things make sense**
"""

# ╔═╡ 1a9cccfb-70c9-4417-b4a6-cc795053676c
root_path = string(BASE_PATH, VENDER)

# ╔═╡ ef42a652-cd7e-41f4-9a26-6dc413a2e04c
dcm_path_list = dcm_list_builder(root_path)

# ╔═╡ 71e79062-183f-4823-856e-ba241cda6f56
pth = dcm_path_list[SCAN_NUMBER]

# ╔═╡ 27db1864-db41-4b20-97f7-306c97d363bc
pth

# ╔═╡ 7023a6ef-ea42-43f6-9f84-a1bbdf24ebc6
scan = basename(pth)

# ╔═╡ b300d7ef-81ce-48d3-9bec-45479a3709ce
header, dcm_array, slice_thick_ori1 = dcm_reader(pth);

# ╔═╡ 8168e2ab-ab93-4886-aed7-3d604fd075a3
md"""
## Helper Functions
"""

# ╔═╡ d6a67da3-b9aa-458e-ae31-5e5c7a4c720d
function collect_tuple(tuple_array)
	row_num = size(tuple_array)
	col_num = length(tuple_array[1])
	container = zeros(Int64, row_num..., col_num)
	for i in 1:length(tuple_array)
		container[i,:] = collect(tuple_array[i])
	end
	return container
end

# ╔═╡ 3fc241f6-bfe4-44f1-b38d-de446360cb33
function overlay_mask_bind(mask)
	indices = findall(x -> x == 1, mask)
	indices = Tuple.(indices)
	label_array = collect_tuple(indices)
	zs = unique(label_array[:,3])
	return PlutoUI.Slider(1:length(zs), default=3, show_value=true)
end

# ╔═╡ 030488fa-08b9-46ce-a959-71661eb96b38
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

# ╔═╡ d56684dd-71da-471d-8e74-79974b3bce38
md"""
## Segment Heart
"""

# ╔═╡ 61c2d1bf-2375-4889-8721-0dc041defdc1
masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3)÷2);

# ╔═╡ df505864-04f8-410a-b04b-efbe5d77be8d
@bind a PlutoUI.Slider(1:size(masked_array, 3), default=10, show_value=true)

# ╔═╡ 122d3232-ca6e-49e8-826a-b3f248185b4a
heatmap(masked_array[:, :, a], colormap=:grays)

# ╔═╡ 0cea0c6f-6034-42ab-987c-5d5431806207
begin
	fig = Figure()
	
	ax = Makie.Axis(fig[1, 1])
	ax.title = "Raw DICOM Array"
	heatmap!(transpose(dcm_array[:, :, 4]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig
end

# ╔═╡ 72bdf950-8f7f-409a-8a84-a8a6d94cdf24
begin
	fig2 = Figure()
	
	ax2 = Makie.Axis(fig2[1, 1])
	ax2.title = "Mask Array"
	heatmap!(transpose(mask), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig2
end

# ╔═╡ b4278464-e3a4-44d8-8c58-dc2b408e491a
begin
	fig3 = Figure()
	
	ax3 = Makie.Axis(fig3[1, 1])
	ax3.title = "Masked DICOM Array"
	heatmap!(transpose(masked_array[:, :, 1]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig3
end

# ╔═╡ 57ba893d-2698-4c7b-a530-1e9da11342ff
md"""
## Segment Calcium Rod
"""

# ╔═╡ 9955eff8-48e0-49da-af33-27b631655b0e
calcium_image, slice_CCI, quality_slice, cal_rod_slice = mask_rod(masked_array, header);

# ╔═╡ a04b7256-f90f-43bb-9e04-376923bb429b
@bind c PlutoUI.Slider(1:size(calcium_image, 3), default=cal_rod_slice, show_value=true)

# ╔═╡ 50ccf844-8ffb-43d5-bf59-900200bd407b
heatmap(transpose(calcium_image[:, :, c]), colormap=:grays)

# ╔═╡ 32427336-d528-4ae6-9add-2c0bdd13d26b
md"""
## Segment Calcium Inserts
"""

# ╔═╡ c7c0ad86-2d72-4ffc-9bfa-4bb9c1eada97
function calc_centers(dcm_array, output, header, tmp_center, CCI_slice, sim::Bool)
    PixelSpacing = Phantoms.get_pixel_size(header)
	center, center1, center2, center3 = center_points(dcm_array, output, header, tmp_center, CCI_slice)
    centers = Dict()
    for center_index in (center1, center2, center3)
        side_x = abs(center[1]-center_index[1])
        side_y = abs(center[2]-center_index[2])
        angle = angle_calc(side_x, side_y)
        if (center_index[1] < center[1] && center_index[2] < center[2])
			medium_calc = [center_index[1] + (10.5 / PixelSpacing[1]) * sin(angle), (center_index[2] + (10.5 / PixelSpacing[2]) * cos(angle))]
			low_calc = [center_index[1] + (17 / PixelSpacing[1]) * sin(angle), (center_index[2] + (17 / PixelSpacing[2]) * cos(angle))]
			
		elseif (center_index[1] < center[1] && center_index[2] > center[2])
			medium_calc = [center_index[1] + (10.5 / PixelSpacing[1]) * sin(angle), (center_index[2] - (10.5 / PixelSpacing[2]) * cos(angle))]
			low_calc = [center_index[1] + (17 / PixelSpacing[1]) * sin(angle), (center_index[2] - (17 / PixelSpacing[2]) * cos(angle))] 
			
		elseif (center_index[1] > center[1] && center_index[2] < center[2])
			medium_calc = [center_index[1] - (10.5 / PixelSpacing[1]) * sin(angle), (center_index[2] + (10.5 / PixelSpacing[2]) * cos(angle))]
			low_calc = [center_index[1] - (17 / PixelSpacing[1]) * sin(angle), (center_index[2] + (17 / PixelSpacing[2]) * cos(angle))]
			
		elseif (center_index[1] > center[1] && center_index[2] > center[2])
			medium_calc = [center_index[1] - (10.5 / PixelSpacing[1]) * sin(angle), (center_index[2] - (10.5 / PixelSpacing[2]) * cos(angle))]
			low_calc = [center_index[1] - (17 / PixelSpacing[1]) * sin(angle), (center_index[2] - (17 / PixelSpacing[2]) * cos(angle))]
			
		elseif (side_x == 0 && center_index[2] < center[2])
			medium_calc = [center_index[1], center_index[2] + (10.5 / PixelSpacing[2])]
			low_calc = [center_index[1], center_index[2] + (17 / PixelSpacing[2])]
			
		elseif (side_x == 0 && center_index[2] > center[2])
			medium_calc = [center_index[1], center_index[2] - (10.5 / PixelSpacing[2])]
			low_calc = [center_index[1], center_index[2] - (17 / PixelSpacing[2])]
			
		elseif (center_index[1] > center[1] && side_y == 0)
            medium_calc = [center_index[1] - (10.5 / PixelSpacing[1]), center_index[2]]
			low_calc = [center_index[1] - (17 / PixelSpacing[1]), center_index[2]]
			
		elseif (center_index[1] > center[1] && side_y == 0)
			medium_calc = [center_index[1] + (10.5 / PixelSpacing[1]), center_index[2]]
            low_calc = [(center_index[1] + (17 / PixelSpacing[1])), center_index[1]]
			
        else
			error("unknown angle")
		end
                
        if center_index == center1
            centers[:Large_HD] = Int.(round.(center_index))
            centers[:Medium_HD] = Int.(round.(medium_calc))
            centers[:Small_HD] = Int.(round.(low_calc))
        
		elseif center_index == center2
            centers[:Large_MD] = Int.(round.(center_index))
            centers[:Medium_MD] = Int.(round.(medium_calc))
            centers[:Small_MD] = Int.(round.(low_calc))
        
		elseif center_index == center3
            centers[:Large_LD] = Int.(round.(center_index))
            centers[:Medium_LD] = Int.(round.(medium_calc))
            centers[:Small_LD] = Int.(round.(low_calc))
        
        else
            nothing
		end
	end
    return centers
end

# ╔═╡ 7dbded2e-7f9b-42dc-be1f-83bb558b647f
function mask_inserts(
    dcm_array,
    masked_array,
    header,
    CCI_slice,
    center_insert,
    sim::Bool;
    calcium_threshold=130,
    comp_connect=trues(3, 3),
)
    output = calc_output(
        masked_array, header, CCI_slice, calcium_threshold, comp_connect
    )
    insert_centers = calc_centers(dcm_array, output, header, center_insert, CCI_slice, sim)

    PixelSpacing = Phantoms.get_pixel_size(header)
    rows, cols = Int(header[(0x0028, 0x0010)]), Int(header[(0x0028, 0x0011)])

    mask_L_HD = create_circular_mask(
        cols, rows, insert_centers[:Large_HD], (round(5 / PixelSpacing[1], RoundUp) / 2) + 1
    )
    mask_L_MD = create_circular_mask(
        cols, rows, insert_centers[:Large_MD], (round(5 / PixelSpacing[1], RoundUp) / 2) + 1
    )
    mask_L_LD = create_circular_mask(
        cols, rows, insert_centers[:Large_LD], (round(5 / PixelSpacing[1], RoundUp) / 2) + 1
    )
    mask_M_HD = create_circular_mask(
        cols,
        rows,
        insert_centers[:Medium_HD],
        (round(3 / PixelSpacing[1], RoundUp) / 2) + 1,
    )
    mask_M_MD = create_circular_mask(
        cols,
        rows,
        insert_centers[:Medium_MD],
        (round(3 / PixelSpacing[1], RoundUp) / 2) + 1,
    )
    mask_M_LD = create_circular_mask(
        cols,
        rows,
        insert_centers[:Medium_LD],
        (round(3 / PixelSpacing[1], RoundUp) / 2) + 1,
    )
    mask_S_HD = create_circular_mask(
        cols, rows, insert_centers[:Small_HD], (round(1 / PixelSpacing[1], RoundUp) / 2) + 1
    )
    mask_S_MD = create_circular_mask(
        cols, rows, insert_centers[:Small_MD], (round(1 / PixelSpacing[1], RoundUp) / 2) + 1
    )
    mask_S_LD = create_circular_mask(
        cols, rows, insert_centers[:Small_LD], (round(1 / PixelSpacing[1], RoundUp) / 2) + 1
    )

    return mask_L_HD,
    mask_M_HD, mask_S_HD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_LD, mask_M_LD,
    mask_S_LD
end

# ╔═╡ 70f704b2-2890-4922-bf1d-a7c4c0b24de9
mask_L_HD, mask_M_HD, mask_S_HD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_LD, mask_M_LD, mask_S_LD = mask_inserts(
            dcm_array, masked_array, header, slice_CCI, center_insert, true
);

# ╔═╡ 22277a4c-2076-4057-920f-0e8bc19c559b
masks = mask_L_HD + mask_M_HD + mask_S_HD + mask_L_MD + mask_M_MD + mask_S_MD + mask_L_LD + mask_M_LD + mask_S_LD;

# ╔═╡ 789bf4d5-80d3-421b-9530-f2522ff93f06
heatmap(masks, colormap=:grays)

# ╔═╡ 90b3343d-bdd3-4578-b7b5-1070774b7c04
md"""
## Calibration Prep
"""

# ╔═╡ b12ded79-c500-4757-be7d-f2e8bbb2bf06
array_filtered = abs.(mapwindow(median, calcium_image[:, :, 2], (3, 3)));

# ╔═╡ dce43cbe-ad7b-4fe5-b641-c1e7968dd8a8
bool_arr = array_filtered .> 0;

# ╔═╡ 1828d1f9-699b-488a-9247-6240182d0b10
bool_arr_erode = (((erode(erode(bool_arr)))));

# ╔═╡ c17e672b-a491-417e-8566-7b82d3e55fa9
heatmap(bool_arr, colormap=:grays)

# ╔═╡ b6d1fc06-7a63-4590-b1e7-6c09f8f9ce1b
heatmap(bool_arr_erode, colormap=:grays)

# ╔═╡ 3c757fa7-9e29-4be0-855d-f21998541453
c_img = calcium_image[:, :, 1:3];

# ╔═╡ e39cc0cb-5b26-4d1c-b5fe-e19556947673
begin
	mask_cal_3D = Array{Bool}(undef, size(c_img))
	for z in 1:size(c_img, 3)
		mask_cal_3D[:, :, z] = bool_arr_erode
	end
end;

# ╔═╡ 6343e16c-2659-4df0-b33c-6bfc540480df
hist(c_img[mask_cal_3D])

# ╔═╡ be75951b-a7a7-4e4b-a940-771050063bd5
# cal_insert_mean2 = mean(c_img[mask_cal_3D])

# ╔═╡ 439c32a1-a305-4f36-b041-f18d47e93fec
cal_insert_mean = quantile!(c_img[mask_cal_3D], 0.7)

# ╔═╡ 6cdf20c5-f60b-4ee9-b693-fcfcc0cb7a6d
md"""
# Score Large Inserts
"""

# ╔═╡ a0355895-f3c7-4d32-a6a0-e2cb8e195911
arr = masked_array[:, :, 4:6];

# ╔═╡ 9d8a066d-4b16-4309-888a-83dbf00ab182
single_arr = masked_array[:, :, slice_CCI];

# ╔═╡ be710eb7-1f52-4dae-9990-c42e49e3249a
md"""
## High Density
"""

# ╔═╡ cdfcc3a9-8566-4b2a-a253-5963ac47921e
begin
	mask_L_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_HD_3D[:, :, z] = transpose(mask_L_HD)
	end
end;

# ╔═╡ bf147c28-211d-4b85-92cd-f6ec98e88b2b
md"""
#### Dilated mask
"""

# ╔═╡ 76b89f07-2092-4f62-8f51-4ec13ad10c09
dilated_mask_L_HD = dilate(dilate(mask_L_HD_3D));

# ╔═╡ 7773df9a-c586-43ee-9908-653659c3f781
@bind g2 overlay_mask_bind(dilated_mask_L_HD)

# ╔═╡ ed27795e-dcad-4189-9a13-1e3b71a3b6ef
overlay_mask_plot(arr, dilated_mask_L_HD, g2, "dilated mask")

# ╔═╡ 9ee143bc-0354-43a8-aef8-5f43b89383cc
md"""
#### Ring (background) mask
"""

# ╔═╡ f18c1912-4539-4005-82ea-c8f7671603be
ring_mask_L_HD = dilate(dilate(dilate(dilate(mask_L_HD_3D)))) - dilate(dilate(dilate(mask_L_HD_3D)));

# ╔═╡ 7fbf6140-eb84-4a4b-ba9d-f7022021bdd9
@bind g4 overlay_mask_bind(ring_mask_L_HD)

# ╔═╡ 5977e54f-fcd2-4b89-b1b5-8864741f4bcd
overlay_mask_plot(arr, ring_mask_L_HD, g4, "ring mask")

# ╔═╡ 5ad33c0d-efbf-43c8-856b-3e52341dc610
md"""
### Calculations
"""

# ╔═╡ 84d926c3-c852-4b6d-8584-bdf376b07225
begin
	single_ring_mask_L_HD = Bool.(ring_mask_L_HD[:, :, 3])
	s_bkg_L_HD = mean(single_arr[single_ring_mask_L_HD])
end

# ╔═╡ 1271a003-1ac2-4700-bc16-4efb2ce9e4c2
pixel_size = DICOMUtils.get_pixel_size(header)

# ╔═╡ 8241e9de-9424-4944-8fae-a6d86d2792d9
begin
	alg_L_HD = Integrated(arr[mask_L_HD_3D])
	ρ = 0.2 # mg/mm^3
	mass_l_hd = score(s_bkg_L_HD, cal_insert_mean, pixel_size, ρ, alg_L_HD)
end

# ╔═╡ cc98187a-ab76-44b2-957e-715da75c1327
md"""
## Medium Density
"""

# ╔═╡ 1df1ad9a-63c5-4d3f-a269-007b30db63b2
begin
	mask_L_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_MD_3D[:, :, z] = transpose(mask_L_MD)
	end
end;

# ╔═╡ 28586ed7-5222-496a-9246-d611438ae626
md"""
#### Dilated mask
"""

# ╔═╡ 5d61a58c-0443-4014-baa9-e70fb86b143c
dilated_mask_L_MD = dilate(dilate(mask_L_MD_3D));

# ╔═╡ e8847548-2b76-47aa-9940-99f4b7262979
@bind h2 overlay_mask_bind(dilated_mask_L_MD)

# ╔═╡ 56434b2d-c78b-4e41-9013-00ba3c1dc616
overlay_mask_plot(arr, dilated_mask_L_MD, h2, "dilated mask")

# ╔═╡ 141be839-42c1-4393-a9dd-c7433e3b3090
md"""
#### Ring (background) mask
"""

# ╔═╡ 74aee07f-d487-40b2-a456-b69b8f1f4d10
ring_mask_L_MD = dilate(dilate(dilate(dilate(mask_L_MD_3D)))) - dilate(dilate(dilate(mask_L_MD_3D)));

# ╔═╡ fb0a75d5-6090-49e6-889c-fca5498285ec
@bind h4 overlay_mask_bind(ring_mask_L_MD)

# ╔═╡ ebd3c1e7-db9f-4d0c-af53-749dba8a645f
overlay_mask_plot(arr, ring_mask_L_MD, h4, "ring mask")

# ╔═╡ a8eeff97-93d5-43e9-99fb-44cdd6ba55b7
md"""
### Calculations
"""

# ╔═╡ 6702778e-2624-4b73-9955-22d7387c4377
begin
	single_ring_mask_L_MD = Bool.(ring_mask_L_MD[:, :, 3])
	s_bkg_L_MD = mean(single_arr[single_ring_mask_L_MD])
end

# ╔═╡ d8246022-8316-4470-9ac8-ab84db51f9a2
begin
	alg_L_MD = Integrated(arr[mask_L_MD_3D])
	mass_l_md = score(s_bkg_L_MD, cal_insert_mean, pixel_size, ρ, alg_L_MD)
end

# ╔═╡ 29b0ffbe-bfcb-491d-8586-cd349ce494b9
md"""
## Low Density
"""

# ╔═╡ 1b70c3b7-f1a1-4a5d-a30e-47e220309016
begin
	mask_L_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_LD_3D[:, :, z] = transpose(mask_L_LD)
	end
end;

# ╔═╡ b66deff6-0ca7-4474-b096-8559f2b12df2
md"""
#### Dilated mask
"""

# ╔═╡ 1f76d736-e7cf-4e1d-aac0-4447cb96cf45
dilated_mask_L_LD = dilate(dilate(mask_L_LD_3D));

# ╔═╡ 88e3f3da-965c-4d15-bbc5-6028574579fb
@bind i2 overlay_mask_bind(dilated_mask_L_LD)

# ╔═╡ 03abf090-704a-4faa-8478-9b405f77947c
overlay_mask_plot(arr, dilated_mask_L_LD, i2, "dilated mask")

# ╔═╡ 1d5a15e0-fc89-4f17-bea0-c7389b27e4d1
md"""
#### Ring (background) mask
"""

# ╔═╡ 7fbd33ac-e597-42c0-87c9-b83761a5a0bc
ring_mask_L_LD = dilate(dilate(dilate(dilate(mask_L_LD_3D)))) - dilate(dilate(dilate(mask_L_LD_3D)));

# ╔═╡ de711839-6210-4c5b-a97c-62de22585525
@bind i4 overlay_mask_bind(ring_mask_L_LD)

# ╔═╡ 12070026-7e8b-46a8-baee-916a6b356520
overlay_mask_plot(arr, ring_mask_L_LD, i4, "ring mask")

# ╔═╡ ac12ec14-c580-4eb3-b41c-a4ed59b297c3
md"""
### Calculations
"""

# ╔═╡ 545a3628-f217-433c-b415-cf9773c1149a
begin	
	single_ring_mask_L_LD = Bool.(ring_mask_L_LD[:, :, 3])
	s_bkg_L_LD = mean(single_arr[single_ring_mask_L_LD])
end

# ╔═╡ 1b902062-70df-4549-b5f2-4fa25462883b
begin
	alg_L_LD = Integrated(arr[mask_L_LD_3D])
	mass_l_ld = score(s_bkg_L_LD, cal_insert_mean, pixel_size, ρ, alg_L_LD)
end

# ╔═╡ 05cd6ce0-3995-474c-94aa-e77a06e2d535
md"""
# Score Medium Inserts
"""

# ╔═╡ 6fa1af38-0769-46d9-aa38-2682e7de4d81
md"""
## High Density
"""

# ╔═╡ 815ac5eb-71ea-4db6-b122-fa9469f45fdb
begin
	mask_M_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_HD_3D[:, :, z] = transpose(mask_M_HD)
	end
end;

# ╔═╡ 457b0987-8dbf-4917-97b8-f3bbd6ff815c
md"""
#### Dilated mask
"""

# ╔═╡ baadfcb3-900e-48f7-a66e-0397ddb4e092
dilated_mask_M_HD = dilate(dilate(dilate(dilate(mask_M_HD_3D))));

# ╔═╡ 52fb81ba-52ff-4caa-b6f2-d03a909193ee
@bind j2 overlay_mask_bind(dilated_mask_M_HD)

# ╔═╡ af6cce1f-07a0-4bc1-a8bd-048ed90bb0b5
overlay_mask_plot(arr, dilated_mask_M_HD, j2, "dilated mask")

# ╔═╡ 4a1ac4e2-17bf-4a99-9cb3-5e21073c2613
md"""
#### Ring (background) mask
"""

# ╔═╡ 2d8075bc-a8d5-4f22-ad4c-a61d6b5c36a0
ring_mask_M_HD = dilate(dilate(dilate(dilate(dilate(mask_M_HD_3D))))) - dilate(dilate(dilate(dilate(mask_M_HD_3D))));

# ╔═╡ 05f2e46b-9a0a-4241-8617-5ec0315be7c5
@bind j4 overlay_mask_bind(ring_mask_M_HD)

# ╔═╡ c3e4a14f-7dd0-4af7-a0cb-add9d2d813aa
overlay_mask_plot(arr, ring_mask_M_HD, j4, "ring mask")

# ╔═╡ f79e6080-b025-46d8-a9a9-2bcc41d68e2e
md"""
### Calculations
"""

# ╔═╡ 4bdc7df3-9ba6-4855-ae80-ca48468e0fd7
begin
	single_ring_mask_M_HD = Bool.(ring_mask_M_HD[:, :, 3])
	s_bkg_M_HD = mean(single_arr[single_ring_mask_M_HD])
end

# ╔═╡ 219dc72c-4959-4a41-a89b-f18915551dac
begin
	alg_M_HD = Integrated(arr[mask_M_HD_3D])
	mass_m_hd = score(s_bkg_M_HD, cal_insert_mean, pixel_size, ρ, alg_M_HD)
end

# ╔═╡ 5ce2f02b-e72a-417b-952b-c7f9d67c3ec1
md"""
## Medium Density
"""

# ╔═╡ 319a043d-4294-45b8-9503-9a690a049b50
begin
	mask_M_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_MD_3D[:, :, z] = transpose(mask_M_MD)
	end
end;

# ╔═╡ 65d9878f-7e81-4456-b7dc-79fc16cc22ae
md"""
#### Dilated mask
"""

# ╔═╡ 3d530e00-f7b7-42d4-b929-e34c73d8c9ec
dilated_mask_M_MD = dilate(dilate(dilate(dilate(mask_M_MD_3D))));

# ╔═╡ 69e3b369-33da-481b-a827-94c3c44b536d
@bind k2 overlay_mask_bind(dilated_mask_M_MD)

# ╔═╡ 816cac04-105e-4a7e-bf57-17433da47dfd
overlay_mask_plot(arr, dilated_mask_M_MD, k2, "dilated mask")

# ╔═╡ 11f1c325-d074-423f-95c5-65ed15afacd1
md"""
#### Ring (background) mask
"""

# ╔═╡ 04b7a3e1-1720-4d38-96fe-1ba5c3f97703
ring_mask_M_MD = dilate(dilate(dilate(dilate(dilate(dilate(mask_M_MD_3D)))))) - dilate(dilate(dilate(dilate(dilate(mask_M_MD_3D)))));

# ╔═╡ c1b276da-54c2-4f6e-b763-63b6b59ae7af
@bind k4 overlay_mask_bind(ring_mask_M_MD)

# ╔═╡ dc68b795-c991-4b67-9f8d-5361424b9a57
overlay_mask_plot(arr, ring_mask_M_MD, k4, "ring mask")

# ╔═╡ 61ebd4a0-d6eb-433d-b0a0-67071db1e74c
md"""
### Calculations
"""

# ╔═╡ 27e1a948-7b3f-4e6c-8b73-12505b0a5e3b
begin
	single_ring_mask_M_MD = Bool.(ring_mask_M_MD[:, :, 3])
	s_bkg_M_MD = mean(single_arr[single_ring_mask_M_MD])
end

# ╔═╡ faace9a2-533f-4bd1-bc66-00c6777895ac
begin
	alg_M_MD = Integrated(arr[mask_M_MD_3D])
	mass_m_md = score(s_bkg_M_MD, cal_insert_mean, pixel_size, ρ, alg_M_MD)
end

# ╔═╡ a57d43c3-4e8b-4061-96c2-00fab1b7144b
md"""
## Low Density
"""

# ╔═╡ a7db6579-1d4d-4d4b-93bd-f13e738bcf66
begin
	mask_M_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_LD_3D[:, :, z] = transpose(mask_M_LD)
	end
end;

# ╔═╡ 65bbff3d-0ad5-48d1-94b9-7566c9fe3b2e
md"""
#### Dilated mask
"""

# ╔═╡ 659b76a1-4295-42dc-b013-73d491c83009
dilated_mask_M_LD = dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))));

# ╔═╡ d78af735-4a8d-4c70-b3be-c24bc26d5ba1
@bind l2 overlay_mask_bind(dilated_mask_M_LD)

# ╔═╡ dcc82a1e-d33a-485b-9eba-12ae58993a3a
overlay_mask_plot(arr, dilated_mask_M_LD, l2, "dilated mask")

# ╔═╡ fb720343-55e3-4c98-8415-a91db70e384c
md"""
#### Ring (background) mask
"""

# ╔═╡ a3c90692-c021-4e01-9841-1ea930331b44
ring_mask_M_LD = dilate(dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))))) - dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))));

# ╔═╡ a90f193c-a81b-4feb-9a22-b494adfcc1ac
@bind l4 overlay_mask_bind(ring_mask_M_LD)

# ╔═╡ 283a555b-a506-4ffb-b589-995c56f6e86e
overlay_mask_plot(arr, ring_mask_M_LD, l4, "ring mask")

# ╔═╡ 207ae3c0-f87a-40f0-83f6-aa15b1266937
md"""
### Calculations
"""

# ╔═╡ a5c6a4fb-9217-4fc7-b896-62fe4256ada3
begin
	single_ring_mask_M_LD = Bool.(ring_mask_M_LD[:, :, 3])
	s_bkg_M_LD = mean(single_arr[single_ring_mask_M_LD])
end

# ╔═╡ 98bbc8d8-4542-475b-8e0d-c4ce00dcefa0
begin
	alg_M_LD = Integrated(arr[mask_M_LD_3D])
	mass_m_ld = score(s_bkg_M_LD, cal_insert_mean, pixel_size, ρ, alg_M_LD)
end

# ╔═╡ e938627e-c9f8-4c15-b14b-c1f5f0c3a72c
md"""
# Score Small Inserts
"""

# ╔═╡ 808ccc93-1154-4336-8bab-08172e0c1490
md"""
## High Density
"""

# ╔═╡ 6aa5d095-1565-4391-a46b-ada2ac17a87c
begin
	mask_S_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_HD_3D[:, :, z] = transpose(mask_S_HD)
	end
end;

# ╔═╡ 140ce3bb-e40f-4c80-9c50-2d4da8e7daaa
md"""
#### Dilated mask
"""

# ╔═╡ 67d46d38-c52a-4fbf-992d-635aca3fb174
dilated_mask_S_HD = dilate(dilate(dilate(dilate(dilate((mask_S_HD_3D))))));

# ╔═╡ 2dbb3af4-e97a-4f49-b4a2-651e14664ba9
@bind m2 overlay_mask_bind(dilated_mask_S_HD)

# ╔═╡ a0dac4d8-d9f7-44e2-8e2b-13711cb91d88
overlay_mask_plot(arr, dilated_mask_S_HD, m2, "dilated mask")

# ╔═╡ 27688ed8-e585-40cd-89ed-47a7b34845db
md"""
#### Ring (background) mask
"""

# ╔═╡ a70fceb8-26eb-42cc-95d2-e531ade92e5d
ring_mask_S_HD = dilate(dilate(dilate(dilate(dilate(mask_S_HD_3D))))) - dilate(dilate(dilate(dilate(mask_S_HD_3D))));

# ╔═╡ d6ab445d-fdea-4c96-84a6-27d270138bd9
@bind m4 overlay_mask_bind(ring_mask_S_HD)

# ╔═╡ 52c8bd34-d0ee-4518-92ee-6ee580377dc6
overlay_mask_plot(arr, ring_mask_S_HD, m4, "ring mask")

# ╔═╡ ac5b3982-3e17-4388-821f-4eae89022862
md"""
### Calculations
"""

# ╔═╡ 6a3fb1a6-db96-433b-add5-061b070aad92
begin
	single_ring_mask_S_HD = Bool.(ring_mask_S_HD[:, :, 3])
	s_bkg_S_HD = mean(single_arr[single_ring_mask_S_HD])
end

# ╔═╡ cfd12f3f-4208-4fe4-897c-7653d19ec255
begin
	alg_S_HD = Integrated(arr[mask_S_HD_3D])
	mass_s_hd = score(s_bkg_S_HD, cal_insert_mean, pixel_size, ρ, alg_S_HD)
	if mass_s_hd < 0
		mass_s_hd = 0
	end
	mass_s_hd
end

# ╔═╡ 61edc2be-eed0-43dd-bbae-df9551924829
md"""
## Medium Density
"""

# ╔═╡ c7704ff4-5382-48fa-bd6f-443974b6c4e5
begin
	mask_S_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_MD_3D[:, :, z] = transpose(mask_S_MD)
	end
end;

# ╔═╡ 64fb36f3-07b9-4d9f-90c0-c0903e6dd95e
md"""
#### Dilated mask
"""

# ╔═╡ 67ebae07-06c9-40b7-940d-c1e47c99c3bb
dilated_mask_S_MD = dilate(dilate(dilate(dilate(dilate(mask_S_MD_3D)))));

# ╔═╡ 6ae57425-9b79-45af-998c-58403b19d7f4
@bind n2 overlay_mask_bind(dilated_mask_S_MD)

# ╔═╡ cffdca3a-b6ef-493b-8289-88f5988ce493
overlay_mask_plot(arr, dilated_mask_S_MD, n2, "dilated mask")

# ╔═╡ cb5fda0b-5a7c-44fb-bf44-7dac4af3b0d7
md"""
#### Ring (background) mask
"""

# ╔═╡ eff04af2-658c-478b-9daa-0a1fe4c25437
ring_mask_S_MD = dilate(dilate(dilate(dilate(dilate(mask_S_MD_3D))))) - dilate(dilate(dilate(dilate(mask_S_MD_3D))));

# ╔═╡ a94b9885-9c0e-44ea-89eb-ca96e415b3b7
@bind n4 overlay_mask_bind(ring_mask_S_MD)

# ╔═╡ 81e2ecdd-e47d-4b2a-a5f8-7d65b66ea625
overlay_mask_plot(arr, ring_mask_S_MD, n4, "ring mask")

# ╔═╡ cdc8427c-6818-4711-a01c-2ae1c0424f9b
md"""
### Calculations
"""

# ╔═╡ 9a8aa8ec-1b3e-4325-8f3a-5368d80f7107
begin
	single_ring_mask_S_MD = Bool.(ring_mask_S_MD[:, :, 3])
	s_bkg_S_MD = mean(single_arr[single_ring_mask_S_MD])
end

# ╔═╡ 023f5749-937b-4e66-a076-c8ce912140b9
begin
	alg_S_MD = Integrated(arr[mask_S_MD_3D])
	mass_s_md = score(s_bkg_S_MD, cal_insert_mean, pixel_size, ρ, alg_S_MD)
	if mass_s_md < 0
		mass_s_md = 0
	end
	mass_s_md
end

# ╔═╡ 3885fe24-ee1d-4a8a-8828-4399d1cbdc9c
md"""
## Low Density
"""

# ╔═╡ 2df14c0d-3e07-47e3-80cc-a756145063e4
begin
	mask_S_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_LD_3D[:, :, z] = transpose(mask_S_LD)
	end
end;

# ╔═╡ 9eb94138-c66c-4201-bfd7-704628510248
md"""
#### Dilated mask
"""

# ╔═╡ f468c6f5-04bd-4c3a-9fb7-c74174b0c73c
dilated_mask_S_LD = dilate(dilate(dilate(dilate(dilate(mask_S_LD_3D)))));

# ╔═╡ c9aba8a8-4546-42a5-bc37-f34290c31c86
@bind o2 overlay_mask_bind(dilated_mask_S_LD)

# ╔═╡ 8bbb947c-1337-4965-98d1-983b40dca864
overlay_mask_plot(arr, dilated_mask_S_LD, o2, "dilated mask")

# ╔═╡ 4beec33f-f6b3-4f83-aa6b-7ff5a57379c8
md"""
#### Ring (background) mask
"""

# ╔═╡ 42f0b7dd-01f9-43ea-b7a6-f17402062651
ring_mask_S_LD = dilate(dilate(dilate(dilate(dilate(mask_S_LD_3D))))) - dilate(dilate(dilate(dilate(mask_S_LD_3D))));

# ╔═╡ e0f3a0bd-2e76-4083-8ee2-edf91fb08b6f
@bind o4 overlay_mask_bind(ring_mask_S_LD)

# ╔═╡ fb99738a-b45e-472d-a88c-2ae42a502eb3
overlay_mask_plot(arr, ring_mask_S_LD, o4, "ring mask")

# ╔═╡ c75c389c-fe68-4dae-8c03-2006785d9421
md"""
### Calculations
"""

# ╔═╡ 5222088c-3bd8-4af0-a6d2-5aea75265ae1
begin
	single_ring_mask_S_LD = Bool.(ring_mask_S_LD[:, :, 3])
	s_bkg_S_LD = mean(single_arr[single_ring_mask_S_LD])
end

# ╔═╡ d3f8aebe-4138-4476-b2a8-60e4c15d4e15
begin
	alg_S_LD = Integrated(arr[mask_S_LD_3D])
	mass_s_ld = score(s_bkg_S_LD, cal_insert_mean, pixel_size, ρ, alg_S_LD)
	if mass_s_ld < 0
		mass_s_ld = 0
	end
	mass_s_ld
end

# ╔═╡ d58dc0f5-4a25-4651-99b7-5ff13b18e9e0
md"""
# Results
"""

# ╔═╡ d0139063-3ead-4dc6-b89f-7051bd994db4
Phantoms.get_pixel_size(header)

# ╔═╡ 91ad40c6-d8c3-42c3-a76a-2a6a03946e82
density_array = [0, 200, 400, 800]

# ╔═╡ 74eefe0d-e5e3-4013-bad2-572fe87015c1
inserts = [
	"Low Density",
	"Medium Density",
	"High Density"
]

# ╔═╡ 5e52f3f3-af70-4463-9272-8e7cfa182078
volume_gt = [
	7.065,
	63.585,
	176.625
]

# ╔═╡ 7afca438-16e7-40b9-b169-b4325492fb3d
ground_truth_mass_large = [
	volume_gt[3] * density_array[2] * 1e-3,
	volume_gt[3] * density_array[3] * 1e-3,
	volume_gt[3] * density_array[4] * 1e-3
] # mg

# ╔═╡ c0f36d80-b667-4617-9694-363c40b61a52
calculated_mass_large = [
	mass_l_ld,
	mass_l_md,
	mass_l_hd
]

# ╔═╡ d2666acd-c4af-4e3f-ac02-0922085214ec
ground_truth_mass_medium = [
	volume_gt[2] * density_array[2] * 1e-3,
	volume_gt[2] * density_array[3] * 1e-3,
	volume_gt[2] * density_array[4] * 1e-3
]

# ╔═╡ 7b174c1d-a2f4-4162-8cf1-98e7e7594869
calculated_mass_medium = [
	mass_m_ld,
	mass_m_md,
	mass_m_hd
]

# ╔═╡ f0261cf3-120d-4205-b993-2365d1455a67
ground_truth_mass_small = [
	volume_gt[1] * density_array[2] * 1e-3,
	volume_gt[1] * density_array[3] * 1e-3,
	volume_gt[1] * density_array[4] * 1e-3
]

# ╔═╡ 3a164f16-970e-4a7a-93b4-e214165e4a42
calculated_mass_small = [
	mass_s_ld,
	mass_s_md,
	mass_s_hd
]

# ╔═╡ f8c23d6d-c1e3-412b-93db-8cdeae3376ba
df = DataFrame(
	inserts = inserts,
	ground_truth_mass_large = ground_truth_mass_large,
	calculated_mass_large = calculated_mass_large,
	ground_truth_mass_medium = ground_truth_mass_medium,
	calculated_mass_medium = calculated_mass_medium,
	ground_truth_mass_small = ground_truth_mass_small,
	calculated_mass_small = calculated_mass_small
)

# ╔═╡ d65fe9f9-adb3-455b-b30a-dea5cb58619d
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

# ╔═╡ f00c4ba8-beee-4221-bfd0-ea1bde6db445
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

# ╔═╡ a27cf678-2fb8-478b-a0e5-d3c8084a2045
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

# ╔═╡ d0641685-3d23-4f33-8dd4-42d95e92bf37
percent_error_large = (abs.(ground_truth_mass_large - calculated_mass_large) ./ ground_truth_mass_large) .* 100

# ╔═╡ f53f8b23-3d62-427f-9455-aee221ed232b
percent_error_medium = (abs.(ground_truth_mass_medium - calculated_mass_medium) ./ ground_truth_mass_medium) .* 100

# ╔═╡ 0fda4495-297a-4596-b4ea-fecaf7173ff7
percent_error_small= (abs.(ground_truth_mass_small - calculated_mass_small) ./ ground_truth_mass_small) .* 100

# ╔═╡ 7871ea81-5683-486e-988d-783ad85b96ec
md"""
### Save Results
"""

# ╔═╡ a4536396-c73f-4fc3-a72a-00884f6f40e8
if ~isdir(string(cd(pwd, "..") , "/data/output/", VENDER))
	mkdir(string(cd(pwd, "..") , "/data/output/", VENDER))
end

# ╔═╡ 267eb596-b901-47cb-840c-638c695f61d4
output_path = string(cd(pwd, "..") , "/data/output/", VENDER, "/", scan, ".csv")

# ╔═╡ a134478a-ac93-4119-b3c2-49192105336c
CSV.write(output_path, df)

# ╔═╡ Cell order:
# ╠═763d240d-eac6-4b56-98de-ae16b32d516a
# ╠═3ca65cc1-c88b-4ddc-8f90-a0726aeb362c
# ╟─d2e68b0d-d60f-4296-ad4e-dfd7a7ec6bb8
# ╠═514daca0-ecef-47b1-ac4c-d5c7d9042c73
# ╟─84bfd706-4cf4-4ea6-af92-1c2c3fcf9944
# ╠═1a9cccfb-70c9-4417-b4a6-cc795053676c
# ╠═ef42a652-cd7e-41f4-9a26-6dc413a2e04c
# ╠═71e79062-183f-4823-856e-ba241cda6f56
# ╠═27db1864-db41-4b20-97f7-306c97d363bc
# ╠═7023a6ef-ea42-43f6-9f84-a1bbdf24ebc6
# ╠═b300d7ef-81ce-48d3-9bec-45479a3709ce
# ╟─8168e2ab-ab93-4886-aed7-3d604fd075a3
# ╟─d6a67da3-b9aa-458e-ae31-5e5c7a4c720d
# ╟─3fc241f6-bfe4-44f1-b38d-de446360cb33
# ╟─030488fa-08b9-46ce-a959-71661eb96b38
# ╟─d56684dd-71da-471d-8e74-79974b3bce38
# ╠═61c2d1bf-2375-4889-8721-0dc041defdc1
# ╟─df505864-04f8-410a-b04b-efbe5d77be8d
# ╠═122d3232-ca6e-49e8-826a-b3f248185b4a
# ╟─0cea0c6f-6034-42ab-987c-5d5431806207
# ╟─72bdf950-8f7f-409a-8a84-a8a6d94cdf24
# ╟─b4278464-e3a4-44d8-8c58-dc2b408e491a
# ╟─57ba893d-2698-4c7b-a530-1e9da11342ff
# ╠═9955eff8-48e0-49da-af33-27b631655b0e
# ╟─a04b7256-f90f-43bb-9e04-376923bb429b
# ╠═50ccf844-8ffb-43d5-bf59-900200bd407b
# ╟─32427336-d528-4ae6-9add-2c0bdd13d26b
# ╟─c7c0ad86-2d72-4ffc-9bfa-4bb9c1eada97
# ╟─7dbded2e-7f9b-42dc-be1f-83bb558b647f
# ╠═70f704b2-2890-4922-bf1d-a7c4c0b24de9
# ╠═22277a4c-2076-4057-920f-0e8bc19c559b
# ╠═789bf4d5-80d3-421b-9530-f2522ff93f06
# ╟─90b3343d-bdd3-4578-b7b5-1070774b7c04
# ╠═b12ded79-c500-4757-be7d-f2e8bbb2bf06
# ╠═dce43cbe-ad7b-4fe5-b641-c1e7968dd8a8
# ╠═1828d1f9-699b-488a-9247-6240182d0b10
# ╠═c17e672b-a491-417e-8566-7b82d3e55fa9
# ╠═b6d1fc06-7a63-4590-b1e7-6c09f8f9ce1b
# ╠═e39cc0cb-5b26-4d1c-b5fe-e19556947673
# ╠═3c757fa7-9e29-4be0-855d-f21998541453
# ╠═6343e16c-2659-4df0-b33c-6bfc540480df
# ╠═be75951b-a7a7-4e4b-a940-771050063bd5
# ╠═439c32a1-a305-4f36-b041-f18d47e93fec
# ╟─6cdf20c5-f60b-4ee9-b693-fcfcc0cb7a6d
# ╠═a0355895-f3c7-4d32-a6a0-e2cb8e195911
# ╠═9d8a066d-4b16-4309-888a-83dbf00ab182
# ╟─be710eb7-1f52-4dae-9990-c42e49e3249a
# ╠═cdfcc3a9-8566-4b2a-a253-5963ac47921e
# ╟─bf147c28-211d-4b85-92cd-f6ec98e88b2b
# ╠═76b89f07-2092-4f62-8f51-4ec13ad10c09
# ╟─7773df9a-c586-43ee-9908-653659c3f781
# ╠═ed27795e-dcad-4189-9a13-1e3b71a3b6ef
# ╟─9ee143bc-0354-43a8-aef8-5f43b89383cc
# ╠═f18c1912-4539-4005-82ea-c8f7671603be
# ╟─7fbf6140-eb84-4a4b-ba9d-f7022021bdd9
# ╠═5977e54f-fcd2-4b89-b1b5-8864741f4bcd
# ╟─5ad33c0d-efbf-43c8-856b-3e52341dc610
# ╠═84d926c3-c852-4b6d-8584-bdf376b07225
# ╠═1271a003-1ac2-4700-bc16-4efb2ce9e4c2
# ╠═8241e9de-9424-4944-8fae-a6d86d2792d9
# ╟─cc98187a-ab76-44b2-957e-715da75c1327
# ╠═1df1ad9a-63c5-4d3f-a269-007b30db63b2
# ╟─28586ed7-5222-496a-9246-d611438ae626
# ╠═5d61a58c-0443-4014-baa9-e70fb86b143c
# ╠═e8847548-2b76-47aa-9940-99f4b7262979
# ╠═56434b2d-c78b-4e41-9013-00ba3c1dc616
# ╟─141be839-42c1-4393-a9dd-c7433e3b3090
# ╠═74aee07f-d487-40b2-a456-b69b8f1f4d10
# ╠═fb0a75d5-6090-49e6-889c-fca5498285ec
# ╠═ebd3c1e7-db9f-4d0c-af53-749dba8a645f
# ╟─a8eeff97-93d5-43e9-99fb-44cdd6ba55b7
# ╠═6702778e-2624-4b73-9955-22d7387c4377
# ╠═d8246022-8316-4470-9ac8-ab84db51f9a2
# ╟─29b0ffbe-bfcb-491d-8586-cd349ce494b9
# ╠═1b70c3b7-f1a1-4a5d-a30e-47e220309016
# ╟─b66deff6-0ca7-4474-b096-8559f2b12df2
# ╠═1f76d736-e7cf-4e1d-aac0-4447cb96cf45
# ╠═88e3f3da-965c-4d15-bbc5-6028574579fb
# ╠═03abf090-704a-4faa-8478-9b405f77947c
# ╟─1d5a15e0-fc89-4f17-bea0-c7389b27e4d1
# ╠═7fbd33ac-e597-42c0-87c9-b83761a5a0bc
# ╟─de711839-6210-4c5b-a97c-62de22585525
# ╠═12070026-7e8b-46a8-baee-916a6b356520
# ╟─ac12ec14-c580-4eb3-b41c-a4ed59b297c3
# ╠═545a3628-f217-433c-b415-cf9773c1149a
# ╠═1b902062-70df-4549-b5f2-4fa25462883b
# ╟─05cd6ce0-3995-474c-94aa-e77a06e2d535
# ╟─6fa1af38-0769-46d9-aa38-2682e7de4d81
# ╠═815ac5eb-71ea-4db6-b122-fa9469f45fdb
# ╟─457b0987-8dbf-4917-97b8-f3bbd6ff815c
# ╠═baadfcb3-900e-48f7-a66e-0397ddb4e092
# ╟─52fb81ba-52ff-4caa-b6f2-d03a909193ee
# ╠═af6cce1f-07a0-4bc1-a8bd-048ed90bb0b5
# ╟─4a1ac4e2-17bf-4a99-9cb3-5e21073c2613
# ╠═2d8075bc-a8d5-4f22-ad4c-a61d6b5c36a0
# ╠═05f2e46b-9a0a-4241-8617-5ec0315be7c5
# ╠═c3e4a14f-7dd0-4af7-a0cb-add9d2d813aa
# ╟─f79e6080-b025-46d8-a9a9-2bcc41d68e2e
# ╠═4bdc7df3-9ba6-4855-ae80-ca48468e0fd7
# ╠═219dc72c-4959-4a41-a89b-f18915551dac
# ╟─5ce2f02b-e72a-417b-952b-c7f9d67c3ec1
# ╠═319a043d-4294-45b8-9503-9a690a049b50
# ╟─65d9878f-7e81-4456-b7dc-79fc16cc22ae
# ╠═3d530e00-f7b7-42d4-b929-e34c73d8c9ec
# ╠═69e3b369-33da-481b-a827-94c3c44b536d
# ╠═816cac04-105e-4a7e-bf57-17433da47dfd
# ╟─11f1c325-d074-423f-95c5-65ed15afacd1
# ╠═04b7a3e1-1720-4d38-96fe-1ba5c3f97703
# ╠═c1b276da-54c2-4f6e-b763-63b6b59ae7af
# ╠═dc68b795-c991-4b67-9f8d-5361424b9a57
# ╟─61ebd4a0-d6eb-433d-b0a0-67071db1e74c
# ╠═27e1a948-7b3f-4e6c-8b73-12505b0a5e3b
# ╠═faace9a2-533f-4bd1-bc66-00c6777895ac
# ╟─a57d43c3-4e8b-4061-96c2-00fab1b7144b
# ╠═a7db6579-1d4d-4d4b-93bd-f13e738bcf66
# ╟─65bbff3d-0ad5-48d1-94b9-7566c9fe3b2e
# ╠═659b76a1-4295-42dc-b013-73d491c83009
# ╟─d78af735-4a8d-4c70-b3be-c24bc26d5ba1
# ╠═dcc82a1e-d33a-485b-9eba-12ae58993a3a
# ╟─fb720343-55e3-4c98-8415-a91db70e384c
# ╠═a3c90692-c021-4e01-9841-1ea930331b44
# ╟─a90f193c-a81b-4feb-9a22-b494adfcc1ac
# ╠═283a555b-a506-4ffb-b589-995c56f6e86e
# ╟─207ae3c0-f87a-40f0-83f6-aa15b1266937
# ╠═a5c6a4fb-9217-4fc7-b896-62fe4256ada3
# ╠═98bbc8d8-4542-475b-8e0d-c4ce00dcefa0
# ╟─e938627e-c9f8-4c15-b14b-c1f5f0c3a72c
# ╟─808ccc93-1154-4336-8bab-08172e0c1490
# ╠═6aa5d095-1565-4391-a46b-ada2ac17a87c
# ╟─140ce3bb-e40f-4c80-9c50-2d4da8e7daaa
# ╠═67d46d38-c52a-4fbf-992d-635aca3fb174
# ╠═2dbb3af4-e97a-4f49-b4a2-651e14664ba9
# ╠═a0dac4d8-d9f7-44e2-8e2b-13711cb91d88
# ╟─27688ed8-e585-40cd-89ed-47a7b34845db
# ╠═a70fceb8-26eb-42cc-95d2-e531ade92e5d
# ╟─d6ab445d-fdea-4c96-84a6-27d270138bd9
# ╠═52c8bd34-d0ee-4518-92ee-6ee580377dc6
# ╟─ac5b3982-3e17-4388-821f-4eae89022862
# ╠═6a3fb1a6-db96-433b-add5-061b070aad92
# ╠═cfd12f3f-4208-4fe4-897c-7653d19ec255
# ╟─61edc2be-eed0-43dd-bbae-df9551924829
# ╠═c7704ff4-5382-48fa-bd6f-443974b6c4e5
# ╟─64fb36f3-07b9-4d9f-90c0-c0903e6dd95e
# ╠═67ebae07-06c9-40b7-940d-c1e47c99c3bb
# ╟─6ae57425-9b79-45af-998c-58403b19d7f4
# ╠═cffdca3a-b6ef-493b-8289-88f5988ce493
# ╟─cb5fda0b-5a7c-44fb-bf44-7dac4af3b0d7
# ╠═eff04af2-658c-478b-9daa-0a1fe4c25437
# ╟─a94b9885-9c0e-44ea-89eb-ca96e415b3b7
# ╠═81e2ecdd-e47d-4b2a-a5f8-7d65b66ea625
# ╟─cdc8427c-6818-4711-a01c-2ae1c0424f9b
# ╠═9a8aa8ec-1b3e-4325-8f3a-5368d80f7107
# ╠═023f5749-937b-4e66-a076-c8ce912140b9
# ╟─3885fe24-ee1d-4a8a-8828-4399d1cbdc9c
# ╠═2df14c0d-3e07-47e3-80cc-a756145063e4
# ╟─9eb94138-c66c-4201-bfd7-704628510248
# ╠═f468c6f5-04bd-4c3a-9fb7-c74174b0c73c
# ╟─c9aba8a8-4546-42a5-bc37-f34290c31c86
# ╠═8bbb947c-1337-4965-98d1-983b40dca864
# ╟─4beec33f-f6b3-4f83-aa6b-7ff5a57379c8
# ╠═42f0b7dd-01f9-43ea-b7a6-f17402062651
# ╟─e0f3a0bd-2e76-4083-8ee2-edf91fb08b6f
# ╠═fb99738a-b45e-472d-a88c-2ae42a502eb3
# ╟─c75c389c-fe68-4dae-8c03-2006785d9421
# ╠═5222088c-3bd8-4af0-a6d2-5aea75265ae1
# ╠═d3f8aebe-4138-4476-b2a8-60e4c15d4e15
# ╟─d58dc0f5-4a25-4651-99b7-5ff13b18e9e0
# ╠═d0139063-3ead-4dc6-b89f-7051bd994db4
# ╠═91ad40c6-d8c3-42c3-a76a-2a6a03946e82
# ╠═74eefe0d-e5e3-4013-bad2-572fe87015c1
# ╠═5e52f3f3-af70-4463-9272-8e7cfa182078
# ╠═7afca438-16e7-40b9-b169-b4325492fb3d
# ╠═c0f36d80-b667-4617-9694-363c40b61a52
# ╠═d2666acd-c4af-4e3f-ac02-0922085214ec
# ╠═7b174c1d-a2f4-4162-8cf1-98e7e7594869
# ╠═f0261cf3-120d-4205-b993-2365d1455a67
# ╠═3a164f16-970e-4a7a-93b4-e214165e4a42
# ╠═f8c23d6d-c1e3-412b-93db-8cdeae3376ba
# ╟─d65fe9f9-adb3-455b-b30a-dea5cb58619d
# ╟─f00c4ba8-beee-4221-bfd0-ea1bde6db445
# ╟─a27cf678-2fb8-478b-a0e5-d3c8084a2045
# ╟─d0641685-3d23-4f33-8dd4-42d95e92bf37
# ╟─f53f8b23-3d62-427f-9455-aee221ed232b
# ╟─0fda4495-297a-4596-b4ea-fecaf7173ff7
# ╟─7871ea81-5683-486e-988d-783ad85b96ec
# ╠═a4536396-c73f-4fc3-a72a-00884f6f40e8
# ╠═267eb596-b901-47cb-840c-638c695f61d4
# ╠═a134478a-ac93-4119-b3c2-49192105336c

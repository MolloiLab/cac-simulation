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

# ╔═╡ e8179970-0653-4a05-8987-55889721ce70
begin
	let
		using Pkg
		Pkg.activate(mktempdir())
		Pkg.Registry.update()
		Pkg.add("PlutoUI")
		Pkg.add("CairoMakie")
		Pkg.add("Statistics")
		Pkg.add("StatsBase")
		Pkg.add("Images")
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
	using Images
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

# ╔═╡ 630218bd-af3e-44a4-9a25-3a808f3c1f2f
TableOfContents()

# ╔═╡ 6dfa4c77-cf36-4855-b6ac-da5789b3f301
md"""
## Load DICOMS
"""

# ╔═╡ a3038e01-851a-4189-bfcf-cac8567c027f
begin
	SCAN_NUMBER = 1
	VENDER = "120"
	TYPE = "integrated_scoring"
	BASE_PATH = "/Users/daleblack/Google Drive/Datasets/Simulated/"
end

# ╔═╡ cffc78f3-e6a3-4c2c-b10d-93115ae80088
md"""
**Everything below should be automatic, just scroll through to visually inspect that things make sense**
"""

# ╔═╡ f473edea-11f9-467d-9c6e-2fe2d93b3cc7
root_path = string(BASE_PATH, VENDER)

# ╔═╡ e1cbaedd-d941-45ea-8d2d-6c8cecb19bb0
dcm_path_list = dcm_list_builder(root_path)

# ╔═╡ 37f95e54-43c1-4eba-b79a-4517f0b65694
pth = dcm_path_list[SCAN_NUMBER]

# ╔═╡ a60b510c-fe72-4cad-9dd6-5836af8a0bb2
scan = basename(pth)

# ╔═╡ 1a0ab1eb-2377-4585-afc8-f35c9a550f76
header, dcm_array, slice_thick_ori1 = dcm_reader(pth);

# ╔═╡ e81d56f6-5490-48dd-b214-84a6740363a5
md"""
## Helper Functions
"""

# ╔═╡ 3208e7f3-9cb6-4d76-a850-209437e921c7
function collect_tuple(tuple_array)
	row_num = size(tuple_array)
	col_num = length(tuple_array[1])
	container = zeros(Int64, row_num..., col_num)
	for i in 1:length(tuple_array)
		container[i,:] = collect(tuple_array[i])
	end
	return container
end

# ╔═╡ b858b3c6-9d0e-4bf4-99b5-6091ee54c311
function overlay_mask_bind(mask)
	indices = findall(x -> x == 1, mask)
	indices = Tuple.(indices)
	label_array = collect_tuple(indices)
	zs = unique(label_array[:,3])
	return PlutoUI.Slider(1:length(zs), default=3, show_value=true)
end

# ╔═╡ 9984a458-1583-4034-82c6-98ecee0487de
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

# ╔═╡ 9cf82f15-c90d-4ec4-84a9-399c6c8159bd
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

# ╔═╡ 2d76a219-f9fc-4c1d-bb41-6d25406a675f
md"""
## Segment Heart
"""

# ╔═╡ 92f3ffe5-2f92-4add-89a8-5a575843af1d
masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3)÷2);

# ╔═╡ 410bb0d5-57a6-4cf4-92ff-8e4cdb45e23b
center_insert

# ╔═╡ 3c19320c-7a8b-4195-8e10-46bf1f11ad00
@bind a PlutoUI.Slider(1:size(masked_array, 3), default=10, show_value=true)

# ╔═╡ 81f5dd90-6985-4cdc-baf4-8644f97bd9ce
heatmap(masked_array[:, :, a], colormap=:grays)

# ╔═╡ 41bcecdf-5590-44d2-b8b8-1f2012889ca1
begin
	fig = Figure()
	
	ax = Makie.Axis(fig[1, 1])
	ax.title = "Raw DICOM Array"
	heatmap!(transpose(dcm_array[:, :, 4]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig
end

# ╔═╡ b4b818ba-18f1-4735-9979-989f4f72c8cb
begin
	fig2 = Figure()
	
	ax2 = Makie.Axis(fig2[1, 1])
	ax2.title = "Mask Array"
	heatmap!(transpose(mask), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig2
end

# ╔═╡ 047ba02b-a8bc-4396-91ad-001b8915be8c
begin
	fig3 = Figure()
	
	ax3 = Makie.Axis(fig3[1, 1])
	ax3.title = "Masked DICOM Array"
	heatmap!(transpose(masked_array[:, :, 1]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig3
end

# ╔═╡ a1ae1c1b-50bc-43d0-a500-81fd74efe9be
md"""
## Segment Calcium Rod
"""

# ╔═╡ 8e2d8b74-c548-4bee-ac6b-9e13abf5c5f0
calcium_image, slice_CCI, quality_slice, cal_rod_slice = mask_rod(masked_array, header);

# ╔═╡ 56cc97eb-3004-47cc-9e63-d67036789f14
@bind c PlutoUI.Slider(1:size(calcium_image, 3), default=cal_rod_slice, show_value=true)

# ╔═╡ 635affbc-b268-4b70-b9b8-3707cf104683
heatmap(transpose(calcium_image[:, :, c]), colormap=:grays)

# ╔═╡ f7d89fa9-b62f-47c7-9cf0-839a951dc953
md"""
## Segment Calcium Inserts
"""

# ╔═╡ 47bac03a-f389-4e85-b0a3-c07d7da0c00b
mask_L_HD, mask_M_HD, mask_S_HD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_LD, mask_M_LD, mask_S_LD = mask_inserts_simulation(
            dcm_array, masked_array, header, slice_CCI, center_insert
);

# ╔═╡ 95df6ad2-cced-4188-bd7f-869580981bf8
masks_tuple = mask_L_HD, mask_M_HD, mask_S_HD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_LD, mask_M_LD, mask_S_LD

# ╔═╡ ce5f4908-1fbb-4220-b194-b395a373c7e1
slice_CCI

# ╔═╡ 794790c7-75bb-4a56-83b6-a3994ffb6551
masks = mask_L_HD + mask_M_HD + mask_S_HD + mask_L_MD + mask_M_MD + mask_S_MD + mask_L_LD + mask_M_LD + mask_S_LD;

# ╔═╡ 4afcd6b7-e407-4b6b-843f-bc30f8c54312
heatmap(masks, colormap=:grays)

# ╔═╡ a2869b97-13fa-415c-8209-afb4449480ab
md"""
# Agatston Scoring
"""

# ╔═╡ 3d857c01-19ea-4709-af8a-217eae1b02b0
md"""
## High Density
"""

# ╔═╡ 371f02b1-89e7-42ac-9396-780027d6f451
arr = masked_array[:, :, 4:6];

# ╔═╡ e7a96263-9f04-488f-a095-3a0313c33e83
begin
	mask_L_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_HD_3D[:, :, z] = dilate(dilate(mask_L_HD))
	end
end;

# ╔═╡ df5aaed8-bfc5-4589-aea3-70efd06de9c2
md"""
### `mass_mean_def`
"""

# ╔═╡ 03cf3b8a-aa85-4a3b-9041-4d2d3fa1817d
function mass_mean_def(array, calcium_threshold)
    # Calculation of mean value within masked image
    if maximum(array) == 0
        mass_mean = 0
    else
        array = Int.(array .> calcium_threshold)
		mass_mean = sum(array) / length(array .!= 0)	
	end
	return mass_mean
end

# ╔═╡ 84dcd7c0-34a9-4c50-8356-2814d159eca5
CCI_array = dcm_array[:, :, 4:6];

# ╔═╡ 2c01f691-0c8c-49b6-99ea-df0801674986
output2 = Phantoms.calc_output(dcm_array, header, 5, 130, trues(3,3));

# ╔═╡ c8dab715-1825-48af-b5ae-353669cafc83
mask_mean_mass = copy(output2[2]);

# ╔═╡ d8960cda-bf64-4f41-a1f5-af6483af313f
inp = abs.(mask_mean_mass .* CCI_array[:, :, 2]);

# ╔═╡ 727a4e34-fb27-4e92-b4c9-56d1c1a08dfb
array = Int.(inp .> 130);

# ╔═╡ 2f81d4c7-cc1e-458c-8c8e-db62d942c830
mass_mean_def(inp, 130)

# ╔═╡ aa6a5e6e-a804-4ee2-a2cf-8efc1cdf54a2
md"""
### `AS_weight_calc`
"""

# ╔═╡ 25eef3b9-056f-4109-bc97-a46f182de38e
function AS_weight_calc(mask, calcium_image_slice, header)
	ImageType = header[(0x0008, 0x0008)]
	maximum_voxel = maximum(mask .* calcium_image_slice)
	local AS_weight
    try
        if calcium_threshold_declaration == "monoE"
            if "ME40KEV" in ImageType
                if maximum_voxel >= 1029
                    AS_weight = 4
				elseif maximum_voxel >= 773
                    AS_weight = 3
                elseif maximum_voxel >= 515
                    AS_weight = 2
                else
                    AS_weight = 1
				end
			elseif "ME50KEV" in ImageType
                if maximum_voxel >= 694
                    AS_weight = 4
                elseif maximum_voxel >= 522
                    AS_weight = 3
                elseif maximum_voxel >= 348
                    AS_weight = 2
                else
                    AS_weight = 1
				end
			elseif "ME60KEV" in ImageType
                if maximum_voxel >= 508
                    AS_weight = 4
                elseif maximum_voxel >= 381
                    AS_weight = 3
                elseif maximum_voxel >= 254
                    AS_weight = 2
                else
                    AS_weight = 1
				end
			elseif "ME70KEV" in ImageType
                if maximum_voxel >= 400
                    AS_weight = 4
                elseif maximum_voxel >= 300
                    AS_weight = 3
                elseif maximum_voxel >= 200
                    AS_weight = 2
                else
                    AS_weight = 1
				end
			elseif "ME80KEV" in ImageType
                if maximum_voxel >= 335
                    AS_weight = 4
                elseif maximum_voxel >= 251
                    AS_weight = 3
                elseif maximum_voxel >= 168
                    AS_weight = 2
                else
                    AS_weight = 1
				end
			elseif "ME90KEV" in ImageType
                if maximum_voxel >= 293
                    AS_weight = 4
                elseif maximum_voxel >= 220
                    AS_weight = 3
                elseif maximum_voxel >= 147
                    AS_weight = 2
                else
                    AS_weight = 1
				end
			elseif "ME100KEV" in ImageType
                if maximum_voxel >= 265
                    AS_weight = 4
                elseif maximum_voxel >= 199
                    AS_weight = 3
                elseif maximum_voxel >= 133
                    AS_weight = 2
                else
                    AS_weight = 1
				end
			elseif "ME110KEV" in ImageType
                if maximum_voxel >= 246
                    AS_weight = 4
                elseif maximum_voxel >= 185
                    AS_weight = 3
                elseif maximum_voxel >= 123
                    AS_weight = 2
                else
                    AS_weight = 1
				end
			elseif "ME120KEV" in ImageType
                if maximum_voxel >= 246
                    AS_weight = 4
                elseif maximum_voxel >= 185
                    AS_weight = 3
                elseif maximum_voxel >= 123
                    AS_weight = 2
                else
                    AS_weight = 1
				end
			elseif "ME130KEV" in ImageType
                if maximum_voxel >= 246
                    AS_weight = 4
                elseif maximum_voxel >= 185
                    AS_weight = 3
                elseif maximum_voxel >= 123
                    AS_weight = 2
                else
                    AS_weight = 1
				end
			elseif "ME140KEV" in ImageType
                if maximum_voxel >= 246
                    AS_weight = 4
                elseif maximum_voxel >= 185
                    AS_weight = 3
                elseif maximum_voxel >= 123
                    AS_weight = 2
                else
                    AS_weight = 1
				end
			elseif "ME150KEV" in ImageType
                if maximum_voxel >= 246
                    AS_weight = 4
                elseif maximum_voxel >= 185
                    AS_weight = 3
                elseif maximum_voxel >= 123
                    AS_weight = 2
                else
                    AS_weight = 1
				end
			elseif "ME160KEV" in ImageType
                if maximum_voxel >= 246
                    AS_weight = 4
                elseif maximum_voxel >= 185
                    AS_weight = 3
                elseif maximum_voxel >= 123
                    AS_weight = 2
                else
                    AS_weight = 1
				end
			elseif "ME170KEV" in ImageType
                if maximum_voxel >= 246
                    AS_weight = 4
                elseif maximum_voxel >= 185
                    AS_weight = 3
                elseif maximum_voxel >= 123
                    AS_weight = 2
                else
                    AS_weight = 1
				end
			elseif "ME180KEV" in ImageType
                if maximum_voxel >= 246
                    AS_weight = 4
                elseif maximum_voxel >= 185
                    AS_weight = 3
                elseif maximum_voxel >= 123
                    AS_weight = 2
                else
                    AS_weight = 1
				end
			elseif "ME190KEV" in ImageType
                if maximum_voxel >= 246
                    AS_weight = 4
                elseif maximum_voxel >= 185
                    AS_weight = 3
                elseif maximum_voxel >= 123
                    AS_weight = 2
                else
                    AS_weight = 1
				end
			end
        else
            if maximum_voxel >= 400
                AS_weight = 4
			elseif maximum_voxel >= 300
                AS_weight = 3
			elseif maximum_voxel >= 200
                AS_weight = 2
            else
                AS_weight = 1
			end
		end
	catch
        if maximum_voxel >= 400
            AS_weight = 4
		elseif maximum_voxel >= 300
            AS_weight = 3
		elseif maximum_voxel >= 200
            AS_weight = 2
        else
            AS_weight = 1
		end
	end
    return AS_weight
end

# ╔═╡ f728b397-79fe-4ba5-b18b-f1860cf275ce
AS_weight_calc(mask_L_HD, CCI_array[:, :, 2], header)

# ╔═╡ e7385151-db55-4a58-9cfb-93843ca3ed2e
md"""
### `calc_mean_calculator`
"""

# ╔═╡ eaead5ad-3af7-4f43-a546-064f92d58332
function calc_mean_calculator(calc_dict_tmp, key_calc_dict, arr_mean, ML_array, slice_ind)
	arr_mean = Int.(arr_mean .== slice_ind)
    tmp_mean_arr = arr_mean .* ML_array
   
    mean_tmp = sum(tmp_mean_arr)
    mean_count = length(tmp_mean_arr .!= 0)
    
    try
        calc_dict_tmp[key_calc_dict][1] += mean_tmp
        calc_dict_tmp[key_calc_dict][2] += mean_count
            
	catch
        calc_dict_tmp[key_calc_dict] = [mean_tmp, mean_count]
	end
    return calc_dict_tmp
end

# ╔═╡ 35a3a6e3-046e-4992-a099-2038fac63880
calc_dict = Dict()

# ╔═╡ b3d2ce33-46e3-4abe-885a-1220a2801959
key_calc_dict = :Large_HD

# ╔═╡ 190a5904-a6ac-40dc-84b6-9f840202c76b
arr_mean = copy(output2[2]);

# ╔═╡ d285f567-eef3-4bb3-bdd8-1d85cea8aee5
ML_array = copy(CCI_array[:, :, 2]);

# ╔═╡ f36e71e8-c3b8-4fa6-8f5b-f542856f5c10
slice_ind = 2;

# ╔═╡ 2acbe48a-9594-4807-b27a-75317ea39ace
md"""
### `min_area_det`
"""

# ╔═╡ ddc01dea-f488-4c55-b459-d728269f8aeb
function min_area_det(header)
	manufacturer = header[(0x0008, 0x0070)]
	PixelSpacing = Phantoms.get_pixel_size(header)
    if "SIEMENS" == manufacturer
        min_area = 0
	elseif "PHILIPS" == manufacturer
        min_area = 0.5
	elseif "CANON" == manufacturer
        min_area = 3 * PixelSpacing[1]^2 #min 3 connected voxels
	elseif "TOSHIBA" == manufacturer
        min_area = 3 * PixelSpacing[1]^2 #min 3 connected voxels
	elseif "GE" == manufacturer
        min_area = 1
	elseif "Literature" == manufacturer
        min_area = 1
    else
        min_area = 1
	end
    return min_area
end

# ╔═╡ 16af2fb8-ae0d-40c9-a2cf-754f196e9f07
min_area_det(header)

# ╔═╡ 25ee1a2c-8e4a-42b5-87ee-29bbb1a15e53
md"""
### `agat`
"""

# ╔═╡ 0f004cc9-f097-4f43-a39b-a82d08863d79
# Agatston and Mass score calculation
function agat(dcm_array, header, masks, CCI_slice, center_insert; calcium_threshold=130, comp_connect=trues(3, 3), min_area=4)
	
	mask_L_HD, mask_M_HD, mask_S_HD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_LD, mask_M_LD, mask_S_LD = masks[1], masks[2], masks[3], masks[4], masks[5], masks[6], masks[7], masks[8], masks[9]
	
	PixelSpacing = Phantoms.get_pixel_size(header)
    SliceThickness = header[(0x0018, 0x0050)]
    # CCI_min = Int((CCI_slice - round(5 / SliceThickness, RoundUp)) - 1)
    # CCI_max = Int((CCI_slice + round(5 / SliceThickness, RoundUp)) + 1)
	CCI_min = 4
	CCI_max = 6
	CCI_array = dcm_array[:,:,CCI_min:CCI_max]
	global outputs = []
	calc_mean_dict = Dict()
    for slice_index in 1:size(CCI_array, 3)
        output = Phantoms.calc_output(dcm_array, header, CCI_slice, calcium_threshold, comp_connect)
		centers = calc_centers(dcm_array, output, header, center_insert, CCI_slice)
		push!(outputs, output)
		centroids = output[4]
		for slice_index2 in 1:size(output, 1)
			coordinates = centroids[slice_index2][2], centroids[slice_index2][1]
			df_area = output[3]
			area = df_area[!, :area][slice_index2] * PixelSpacing[1]^2
			mask_mean_mass = copy(output[2])
			mask_mean_mass = Int.(mask_mean_mass .== slice_index2)
			min_area = min_area_det(header)
					
			if area > min_area
				MS = area * SliceThickness .* mass_mean_def(mask_mean_mass .* CCI_array[:, :, slice_index], calcium_threshold)
				
				if mask_L_HD[coordinates] != false
					AS_weight = AS_weight_calc(mask_L_HD, CCI_array[:,:,slice_index])
					centers[:Large_HD][2] += AS_weight * area
					
					centers[:Large_HD][3] += MS                        
					calc_mean_dict = calc_mean_calculator(calc_mean_dict, :Large_HD, copy(output[2]), copy(CCI_array[:, :, slice_index]), slice_index2)
				end
						
	
		# 	# 	elif mask_L_MD[coordinates] != False:
		# 	# 		AS_weight = AS_weight_calc(mask_L_MD, CCI_array[:,:,slice_index])
		# 	# 		calc_size_density_VS_AS_MS['Large_MD'][2] += AS_weight * area
					
		# 	# 		calc_size_density_VS_AS_MS['Large_MD'][3] += MS                    
		# 	# 		calc_mean_dict = calc_mean_calculator(calc_mean_dict, 'Large_MD', output[1].copy(),\
		# 	# 						   CCI_array[:,:,slice_index].copy(), slice_index2)
				
		# 	# 	elif mask_L_LD[coordinates] != False:
		# 	# 		AS_weight = AS_weight_calc(mask_L_LD, CCI_array[:,:,slice_index])
		# 	# 		calc_size_density_VS_AS_MS['Large_LD'][2] += AS_weight * area
	
		# 	# 		calc_size_density_VS_AS_MS['Large_LD'][3] += MS                     
		# 	# 		calc_mean_dict = calc_mean_calculator(calc_mean_dict, 'Large_LD', output[1].copy(),\
		# 	# 						   CCI_array[:,:,slice_index].copy(), slice_index2)
				
		# 	# 	elif mask_M_HD[coordinates] != False:
		# 	# 		AS_weight = AS_weight_calc(mask_M_HD, CCI_array[:,:,slice_index])
		# 	# 		calc_size_density_VS_AS_MS['Medium_HD'][2] += AS_weight * area
	
		# 	# 		calc_size_density_VS_AS_MS['Medium_HD'][3] += MS                    
		# 	# 		calc_mean_dict = calc_mean_calculator(calc_mean_dict, 'Medium_HD', output[1].copy(),\
		# 	# 						   CCI_array[:,:,slice_index].copy(), slice_index2)
				
		# 	# 	elif mask_M_MD[coordinates] != False:
		# 	# 		AS_weight = AS_weight_calc(mask_M_MD, CCI_array[:,:,slice_index])
		# 	# 		calc_size_density_VS_AS_MS['Medium_MD'][2] += AS_weight * area
					
		# 	# 		calc_size_density_VS_AS_MS['Medium_MD'][3] += MS                     
		# 	# 		calc_mean_dict = calc_mean_calculator(calc_mean_dict, 'Medium_MD', output[1].copy(),\
		# 	# 						   CCI_array[:,:,slice_index].copy(), slice_index2)
			   
		# 	# 	elif mask_M_LD[coordinates] != False:
		# 	# 		AS_weight = AS_weight_calc(mask_M_LD, CCI_array[:,:,slice_index])
		# 	# 		calc_size_density_VS_AS_MS['Medium_LD'][2] += AS_weight * area
					
		# 	# 		calc_size_density_VS_AS_MS['Medium_LD'][3] += MS  
					
		# 	# 		calc_mean_dict = calc_mean_calculator(calc_mean_dict, 'Medium_LD', output[1].copy(),\
		# 	# 						   CCI_array[:,:,slice_index].copy(), slice_index2)
				
		# 	# 	elif mask_S_HD[coordinates] != False:
		# 	# 		AS_weight = AS_weight_calc(mask_S_HD, CCI_array[:,:,slice_index])
		# 	# 		calc_size_density_VS_AS_MS['Small_HD'][2] += AS_weight * area
	
		# 	# 		calc_size_density_VS_AS_MS['Small_HD'][3] += MS 
		# 	# 		calc_mean_dict = calc_mean_calculator(calc_mean_dict, 'Small_HD', output[1].copy(),\
		# 	# 						   CCI_array[:,:,slice_index].copy(), slice_index2)
				
		# 	# 	elif mask_S_MD[coordinates] != False:
		# 	# 		AS_weight = AS_weight_calc(mask_S_MD, CCI_array[:,:,slice_index])
		# 	# 		calc_size_density_VS_AS_MS['Small_MD'][2] += AS_weight * area
					
		# 	# 		calc_size_density_VS_AS_MS['Small_MD'][3] += MS
		# 	# 		calc_mean_dict = calc_mean_calculator(calc_mean_dict, 'Small_MD', output[1].copy(),\
		# 	# 						   CCI_array[:,:,slice_index].copy(), slice_index2)
				
		# 	# 	elif mask_S_LD[coordinates] != False:
		# 	# 		AS_weight = AS_weight_calc(mask_S_LD, CCI_array[:,:,slice_index])
		# 	# 		calc_size_density_VS_AS_MS['Small_LD'][2] += AS_weight * area
	
		# 	# 		calc_size_density_VS_AS_MS['Small_LD'][3] += MS
		# 	# 		calc_mean_dict = calc_mean_calculator(calc_mean_dict, 'Small_LD', output[1].copy(),\
		# 	# 						   CCI_array[:,:,slice_index].copy(), slice_index2)
		# 	# 	else:
		# 	# 		pass
			end
		# 	# else:
		# 	# 	pass
		end
	end
end

# ╔═╡ b2f9b70e-daed-4269-9f43-0a0afec5c76f
agat(dcm_array, header, masks_tuple, slice_CCI, center_insert)

# ╔═╡ 9eb7ca8b-1b6d-4335-bebf-3aae56f8621e
mask_mean_mass

# ╔═╡ 4346b89b-be0b-4a96-be6f-ea7a5fadd440
show_matrix(mask_mean_mass)

# ╔═╡ 175a38e5-1980-4ab9-9234-8b65fad5639f
output = outputs[1]

# ╔═╡ 6086d9b7-9313-451f-9b12-75d93215b900
outputs

# ╔═╡ 45e2d42c-55e8-44e2-8ec0-fc5eab1f68e3


# ╔═╡ 1c6cbddc-0c53-4984-aae4-d8f902fad26a


# ╔═╡ 3b4b920c-9bf1-11ec-3eb9-c71a55cef1ac
function agatston_vol(
	image, spacing, mask;
	comp_connect=comp_connect=trues(3, 3), min_vol=nothing, max_vol=nothing
		)
	
    voxel_volume = spacing[1] * spacing[2] * spacing[3]
    agatston_score = 0
    calcium_volume = 0

    # Find individual lesions (in 3D) so that we can discard too small or too large lesions
	for z in size(mask, 3)
		lesion_map = label_components(mask[:, :, z], comp_connect)
	    n_lesions = length(unique(lesion_map))

	    for lesion in 1:(n_lesions - 1)
	        lesion_mask = map(x -> x == lesion, lesion_map)
	
	        # Ignore too small or too large lesions
	        lesion_volume = count(lesion_mask) * voxel_volume
	        if ((min_vol != nothing) && (lesion_volume < min_vol))
	            continue
	        end
	        if ((max_vol != nothing) && (lesion_volume > max_vol))
	            continue
	        end
	
	        calcium_volume += lesion_volume
	
	        # Calculate Agatston score for this lesion
	        slices = sort(unique(lesion_mask)) .+ 1
	        for z in 1:slices[end]
	            fragment_mask = lesion_mask
	            n_pixels = count(fragment_mask)
	            maximum_intensity = maximum(image[:, :, z])
	            if maximum_intensity < 200
	                coefficient = 1
	            elseif maximum_intensity < 300
	                coefficient = 2
	            elseif maximum_intensity < 400
	                coefficient = 3
	            else
	                coefficient = 4
	            end
	            agatston_score += coefficient * n_pixels
	        end
	    end
	end
    agatston_score *= spacing[1] / 3.0 * spacing[2] * spacing[3]
    return agatston_score, calcium_volume
end

# ╔═╡ 0e4902fe-3b31-4a56-ae20-b898e42b2eba
spacing = Phantoms.get_pixel_size(header)

# ╔═╡ 3b0119e3-1d46-46e7-9b93-db0b8a194e9f
append!(spacing, header[tag"Slice Thickness"])

# ╔═╡ a8ed22ff-ab66-4a78-a627-0d73af033a31
agatston_score, calcium_volume = agatston_vol(arr, spacing, mask_L_HD_3D)

# ╔═╡ Cell order:
# ╠═e8179970-0653-4a05-8987-55889721ce70
# ╠═630218bd-af3e-44a4-9a25-3a808f3c1f2f
# ╟─6dfa4c77-cf36-4855-b6ac-da5789b3f301
# ╠═a3038e01-851a-4189-bfcf-cac8567c027f
# ╟─cffc78f3-e6a3-4c2c-b10d-93115ae80088
# ╠═f473edea-11f9-467d-9c6e-2fe2d93b3cc7
# ╠═e1cbaedd-d941-45ea-8d2d-6c8cecb19bb0
# ╠═37f95e54-43c1-4eba-b79a-4517f0b65694
# ╠═a60b510c-fe72-4cad-9dd6-5836af8a0bb2
# ╠═1a0ab1eb-2377-4585-afc8-f35c9a550f76
# ╟─e81d56f6-5490-48dd-b214-84a6740363a5
# ╟─3208e7f3-9cb6-4d76-a850-209437e921c7
# ╟─b858b3c6-9d0e-4bf4-99b5-6091ee54c311
# ╟─9984a458-1583-4034-82c6-98ecee0487de
# ╟─9cf82f15-c90d-4ec4-84a9-399c6c8159bd
# ╟─2d76a219-f9fc-4c1d-bb41-6d25406a675f
# ╠═92f3ffe5-2f92-4add-89a8-5a575843af1d
# ╠═410bb0d5-57a6-4cf4-92ff-8e4cdb45e23b
# ╟─3c19320c-7a8b-4195-8e10-46bf1f11ad00
# ╠═81f5dd90-6985-4cdc-baf4-8644f97bd9ce
# ╟─41bcecdf-5590-44d2-b8b8-1f2012889ca1
# ╟─b4b818ba-18f1-4735-9979-989f4f72c8cb
# ╟─047ba02b-a8bc-4396-91ad-001b8915be8c
# ╟─a1ae1c1b-50bc-43d0-a500-81fd74efe9be
# ╠═8e2d8b74-c548-4bee-ac6b-9e13abf5c5f0
# ╠═635affbc-b268-4b70-b9b8-3707cf104683
# ╟─56cc97eb-3004-47cc-9e63-d67036789f14
# ╟─f7d89fa9-b62f-47c7-9cf0-839a951dc953
# ╠═47bac03a-f389-4e85-b0a3-c07d7da0c00b
# ╠═95df6ad2-cced-4188-bd7f-869580981bf8
# ╠═ce5f4908-1fbb-4220-b194-b395a373c7e1
# ╠═794790c7-75bb-4a56-83b6-a3994ffb6551
# ╠═4afcd6b7-e407-4b6b-843f-bc30f8c54312
# ╟─a2869b97-13fa-415c-8209-afb4449480ab
# ╟─3d857c01-19ea-4709-af8a-217eae1b02b0
# ╠═371f02b1-89e7-42ac-9396-780027d6f451
# ╠═e7a96263-9f04-488f-a095-3a0313c33e83
# ╟─df5aaed8-bfc5-4589-aea3-70efd06de9c2
# ╠═03cf3b8a-aa85-4a3b-9041-4d2d3fa1817d
# ╠═84dcd7c0-34a9-4c50-8356-2814d159eca5
# ╠═2c01f691-0c8c-49b6-99ea-df0801674986
# ╠═c8dab715-1825-48af-b5ae-353669cafc83
# ╠═d8960cda-bf64-4f41-a1f5-af6483af313f
# ╠═727a4e34-fb27-4e92-b4c9-56d1c1a08dfb
# ╠═2f81d4c7-cc1e-458c-8c8e-db62d942c830
# ╟─aa6a5e6e-a804-4ee2-a2cf-8efc1cdf54a2
# ╠═25eef3b9-056f-4109-bc97-a46f182de38e
# ╠═f728b397-79fe-4ba5-b18b-f1860cf275ce
# ╟─e7385151-db55-4a58-9cfb-93843ca3ed2e
# ╠═eaead5ad-3af7-4f43-a546-064f92d58332
# ╠═35a3a6e3-046e-4992-a099-2038fac63880
# ╠═b3d2ce33-46e3-4abe-885a-1220a2801959
# ╠═190a5904-a6ac-40dc-84b6-9f840202c76b
# ╠═d285f567-eef3-4bb3-bdd8-1d85cea8aee5
# ╠═f36e71e8-c3b8-4fa6-8f5b-f542856f5c10
# ╟─2acbe48a-9594-4807-b27a-75317ea39ace
# ╠═ddc01dea-f488-4c55-b459-d728269f8aeb
# ╠═16af2fb8-ae0d-40c9-a2cf-754f196e9f07
# ╟─25ee1a2c-8e4a-42b5-87ee-29bbb1a15e53
# ╠═0f004cc9-f097-4f43-a39b-a82d08863d79
# ╠═b2f9b70e-daed-4269-9f43-0a0afec5c76f
# ╠═9eb7ca8b-1b6d-4335-bebf-3aae56f8621e
# ╠═4346b89b-be0b-4a96-be6f-ea7a5fadd440
# ╠═175a38e5-1980-4ab9-9234-8b65fad5639f
# ╠═6086d9b7-9313-451f-9b12-75d93215b900
# ╠═45e2d42c-55e8-44e2-8ec0-fc5eab1f68e3
# ╠═1c6cbddc-0c53-4984-aae4-d8f902fad26a
# ╠═3b4b920c-9bf1-11ec-3eb9-c71a55cef1ac
# ╠═0e4902fe-3b31-4a56-ae20-b898e42b2eba
# ╠═3b0119e3-1d46-46e7-9b93-db0b8a194e9f
# ╠═a8ed22ff-ab66-4a78-a627-0d73af033a31

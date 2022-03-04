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

# ╔═╡ 946debe7-078f-4db6-a383-d0b5a973bd4c
begin
	let
		using Pkg
		Pkg.activate(mktempdir())
		Pkg.Registry.update()
		Pkg.add("PlutoUI")
		Pkg.add("ImageFiltering")
		Pkg.add("Images")
		Pkg.add("ImageMorphology")
		Pkg.add("ImageSegmentation")
		Pkg.add("ImageComponentAnalysis")
		Pkg.add("DataFrames")
		Pkg.add("Statistics")
		Pkg.add("CairoMakie")
		Pkg.add("DataStructures")
		Pkg.add("MAT")
		Pkg.add(url="https://github.com/JuliaHealth/DICOM.jl")
		Pkg.add(url="https://github.com/Dale-Black/DICOMUtils.jl")
		Pkg.add(url="https://github.com/Dale-Black/Phantoms.jl")
	end

	using PlutoUI
	using ImageFiltering
	using Images
	using ImageMorphology
	using ImageSegmentation
	using ImageComponentAnalysis
	using DataFrames
	using Statistics
	using CairoMakie
	using DataStructures
	using MAT
	using DICOM
	using DICOMUtils
	using Phantoms
end

# ╔═╡ 601c4b6d-ca82-4fb0-b66b-5b50d507fa83
TableOfContents()

# ╔═╡ 3a43e3e5-adc3-4a0d-ad98-1a7e39dc6a69
begin
	root_dir = dirname(pwd())
	root = string(root_dir, "/data/simulated/combined")
end

# ╔═╡ 20e6e9ba-2d70-4d47-9187-95b903d95814
# root = "/Users/daleblack/Google Drive/Datasets/Simulated/combined";

# ╔═╡ 65b807a2-8e49-4065-baf7-dca0534a7f20
md"""
## `dcm_list_builder`
"""

# ╔═╡ cfb8b5b9-cfaa-4b4f-81bd-3bc7b2eb30ba
function dcm_list_builder(path)
    dcm_path_list = []
    for (dirpath, dirnames, filenames) in walkdir(path, topdown=true)
        if (dirpath in dcm_path_list) == false
            for filename in filenames
                try
                    tmp_str = string(dirpath, "/", filename)
                    ds = dcm_parse(tmp_str)
                    if (dirpath in dcm_path_list) == false
                        push!(dcm_path_list, dirpath)
					end
				catch
                    nothing
				end
			end
        else
                nothing
		end
	end
    return dcm_path_list
end

# ╔═╡ 34221566-cdd5-480a-85c1-8262dda92a38
md"""
## `dcm_reader`
"""

# ╔═╡ defe2c1c-576f-451f-b4b9-ac88411bf45e
function dcm_reader(dcm_path)
    dcm_files = []
    for (dirpath, dirnames, filenames) in walkdir(dcm_path, topdown=false)
        for filename in filenames
            try
                if (filename == "DIRFILE") == false   
                    dcm_file = string(dirpath, "/", filename)
                    dcm_parse(dcm_file)
                    push!(dcm_files, dcm_file)
				end
			catch
				nothing
			end
		end
	end

    read_RefDs = true
	local RefDs
    while read_RefDs
        for index in range(1, length(dcm_files))
            try
                RefDs = dcm_parse(dcm_files[index])
                read_RefDs = false
                break
			catch
                nothing
			end
		end
	end

	header = RefDs.meta
	slice_thick_ori = header[(0x0018, 0x0050)]
	rows, cols = Int(header[(0x0028, 0x0010)]), Int(header[(0x0028, 0x0011)])
    
    ConstPixelDims = (rows, cols, length(dcm_files))
    dcm_array = zeros(ConstPixelDims...)

    instances = []    
    for filenameDCM in dcm_files
        try
            ds = dcm_parse(filenameDCM)
			head = ds.meta
			InstanceNumber = head[(0x0020, 0x0013)]
            push!(instances, InstanceNumber)
		catch
            nothing
		end
	end
    
    sort!(instances)

    index = 0
    for filenameDCM in dcm_files
        try
            ds = dcm_parse(filenameDCM)
			head = ds.meta
			InstanceNumber = head[(0x0020, 0x0013)]
			index = findall(x -> x==InstanceNumber, instances)
			pixel_array = head[(0x7fe0, 0x0010)]
            dcm_array[:, :, index] = pixel_array
            index += 1
		catch
            nothing
		end
	end
	
    RescaleSlope = header[(0x0028, 0x1053)]
	RescaleIntercept = header[(0x0028, 0x1052)]
    dcm_array = dcm_array .* RescaleSlope .+ RescaleIntercept
    return RefDs.meta, dcm_array, slice_thick_ori
end

# ╔═╡ cf853da4-c254-4a44-abc2-27b871bba172
md"""
## Load DICOMs
"""

# ╔═╡ 2f8f212d-d897-4034-a63c-e635de77a6fc
dcm_path_list = dcm_list_builder(root)

# ╔═╡ f8c7a77e-b37f-4c6f-801c-5e5636cb46d5
dcm_reader(dcm_path_list[1]);

# ╔═╡ 6bbd91e2-1832-4a01-87dc-0c253452e360
header, dcm_array, slice_thick_ori = dcm_reader(dcm_path_list[1]);

# ╔═╡ d39cd4ac-8fa4-46d6-b0ae-0d99c912b8d7
heatmap(transpose(dcm_array[:, :, 1]), colormap=:grays)

# ╔═╡ 936c2be7-072f-4774-8261-3780e7e0897e
maximum(dcm_array)

# ╔═╡ fd9bba1e-27a7-48d9-b44d-14c213582e07
minimum(dcm_array)

# ╔═╡ d1f83e63-e42a-4a1c-a5eb-c4dfe93679ac
md"""
## `find_circle`
"""

# ╔═╡ 5c18f738-2e97-4e05-9a4c-f688895f5d83
function find_circle(point_1, point_2, point_3)
    x1, y1 = point_1
    x2, y2 = point_2
    x3, y3 = point_3
    
    x12 = x1 - x2 
    x13 = x1 - x3  
    y12 = y1 - y2  
    y13 = y1 - y3 
    y31 = y3 - y1  
    y21 = y2 - y1
    x31 = x3 - x1  
    x21 = x2 - x1 
 
    sx13 = x1^2 - x3^2  
    sy13 = y1^2 - y3^2
    sx21 = x2^2 - x1^2  
    sy21 = y2^2 - y1^2  
  
    f = (((sx13) * (x12) + (sy13) * (x12) + (sx21) * (x13) + (sy21) * (x13)) ÷ (2 * ((y31) * (x12) - (y21) * (x13)))) 
              
    g = (((sx13) * (y12) + (sy13) * (y12) + (sx21) * (y13) + (sy21) * (y13)) ÷ (2 * ((x31) * (y12) - (x21) * (y13))))  
  
    # eqn of circle be x^2 + y^2 + 2*g*x + 2*f*y + c = 0 where center is (h = -g, k = -f)  
    center_insert = [-g, -f]

    return center_insert
end

# ╔═╡ 21ad5a93-7121-4fff-95d7-c5a6b8bff336
find_circle([309, 309], [312, 200], [155, 155])

# ╔═╡ 9ed39437-e066-46da-94bd-0cbacbbe9981
md"""
## `mask_heart`
"""

# ╔═╡ 5f389fbc-655f-4c32-97b8-1584f5990181
function mask_heart(
	header, array_used, slice_used_center; 
	radius_val=95, 
	)
	
	pixel_size = Phantoms.get_pixel_size(header)
    
    radius = (radius_val / 2) / pixel_size[1]
    central_image = copy(array_used[:, :, slice_used_center])
    central_image = Int.(central_image .< -200)
    kern = Int.(round(5 / pixel_size[1]))
    if kern % 2 == 0
        kern += 1
	end
    central_image = mapwindow(median, central_image, (kern, kern))
    center = [size(central_image, 1) ÷ 2, size(central_image, 2) ÷ 2]
    a = copy(central_image)
    local point_1
	for index in 1:size(central_image, 2) ÷ 2
		if (central_image[center[1] + index, center[2] + index] == 1 && central_image[center[1] + index, center[2] + index + 5] == 1) 
			point_1 = [center[1] + index, center[2] + index]
			break
		else
			a[center[1] + index, center[2] + index] = 2
		end
	end
    
	local point_2
	for index in 1:size(central_image, 2) ÷ 2
		if (central_image[center[1] + index, center[2] - index] == 1 && central_image[center[1] + index, center[2] - index - 5] == 1) 
			point_2 = [center[1] + index, center[2] - index]
			break
		else
			a[center[1] + index, center[2] - index] = 2
		end
	end
	
	local point_3
	for index in 1:size(central_image, 2) ÷ 2
		if (central_image[center[1] - index, center[2] - index] == 1 && central_image[center[1] - index, center[2] - index - 5] == 1)
			point_3 = [center[1] - index, center[2] - index]
			break
		else
			a[center[1] - index, center[2] - index] = 2
		end
	end
	
    center_insert = find_circle(point_1, point_2, point_3)
	rows, cols = Int(header[(0x0028, 0x0010)]), Int(header[(0x0028, 0x0011)])
    Y, X = collect(1:rows), collect(1:cols)'
    dist_from_center = @. sqrt((X - center_insert[2])^2 + (Y-center_insert[1])^2)

    mask = dist_from_center .<= radius[1]  
    masked_array = zeros(size(array_used))
    for index in 1:size(array_used, 3)
        masked_array[:, :, index] = array_used[:, :, index] .* mask
	end

    return masked_array, center_insert, mask
end

# ╔═╡ 26b57188-0d11-4f4d-bed3-37f7d8b37a96
masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3)÷2);

# ╔═╡ 544734b7-a7de-41ac-b2ea-9587d1409bab
center_insert

# ╔═╡ f1f2f519-2ee3-4392-94f0-6114da0c98f4
begin
	fig = Figure()
	
	ax = Makie.Axis(fig[1, 1])
	ax.title = "Raw DICOM Array"
	heatmap!(transpose(dcm_array[:, :, 2]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig
end

# ╔═╡ 0b40f43c-78d1-4714-a840-3788ba8867b6
begin
	fig2 = Figure()
	
	ax2 = Makie.Axis(fig2[1, 1])
	ax2.title = "Mask Array"
	heatmap!(transpose(mask), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig2
end

# ╔═╡ 6ce0efe4-6d4f-4ff4-983b-3932d75b5dc9
begin
	fig3 = Figure()
	
	ax3 = Makie.Axis(fig3[1, 1])
	ax3.title = "Masked DICOM Array"
	heatmap!(transpose(masked_array[:, :, 2]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig3
end

# ╔═╡ 3b0263c8-280d-4081-a1ef-b24bf1a841fe
@bind a2 PlutoUI.Slider(1:size(masked_array, 3), default=10, show_value=true)

# ╔═╡ 278f03d7-6f86-432f-9d48-6574e596fad4
heatmap(transpose(masked_array[:, :, a2]), colormap=:grays)

# ╔═╡ 01e277ce-dfa7-4b8f-bd7e-2bd5c52c8def
md"""
# Calcium rod mask
"""

# ╔═╡ 20c7b163-c2c2-4dd3-9cf2-1aeca69df7c0
md"""
## `get_calcium_slices`

Returns the slices that contain the calcium vessel inserts `slice_dict` and the slices that contain the calcium calibration rod insert `large_index`.

The calcium rod slice `large_index` usually omits the first and last slice of the rod. So to account for all of the slices containing the calcium rod, one would want a range like so: `(large_index - 1):(large_index + 1)`
"""

# ╔═╡ 8363142b-2909-4133-9a50-a43232f3e82d
function get_calcium_slices(dcm_array, header; calcium_threshold=130)
    array = copy(dcm_array)
    array = Int.(array .> (1.1 * calcium_threshold))

	pixel_size = Phantoms.get_pixel_size(header)
    CCI_5mm_num_pixels = Int(round(π * (5/2)^2 / pixel_size[1]^2))
    cal_rod_num_pixels = Int(round(π * (20/2)^2 / pixel_size[1]^2))
    
	kern = Int.(round(5 / pixel_size[1]))
    if kern % 2 == 0
        kern += 1
	end
    
    slice_dict = Dict()
    large_index = []
    cal_rod_dict = Dict()
    for idx in 1:size(array, 3)
		array_filtered = mapwindow(median, array[:,:,idx], (kern, kern))
        components = ImageComponentAnalysis.label_components(array_filtered)
		a1 = analyze_components(
			components, BasicMeasurement(area=true, perimeter=true)
		)
		a2 = analyze_components(components, BoundingBox(box_area = true))
		df = leftjoin(a1, a2, on = :l)

		count_5mm = 0
		count = 0
		for row in eachrow(df)
			count += 1
			df_area = Int(round(row[:area]))
			
			r1_1 = Int(round(CCI_5mm_num_pixels * 0.6))
			r1_2 = Int(round(CCI_5mm_num_pixels * 1.5))
			r2_1 = Int(round(cal_rod_num_pixels * 0.7))
			r2_2 = Int(round(cal_rod_num_pixels * 1.3))
			
			if df_area in r1_1:r1_2
				count_5mm += 1
			elseif df_area in r2_1:r2_2
				indices = row[:box_indices]
				x_point = ((indices[1][end] - indices[1][1]) ÷ 2) + indices[1][1]
				y_point = ((indices[2][end] - indices[2][1]) ÷ 2) + indices[2][1]
				cal_rod_dict[count] = [x_point, y_point]
			end
		end
		
        if count_5mm > 0 && count_5mm < 4
            slice_dict[idx] = count_5mm
		end
    
        poppable_keys = []
        for key in cal_rod_dict
            start_coordinate = [key[2][1], key[2][2]]
			
            x_right = 0
            while array_filtered[start_coordinate[1], start_coordinate[2] + x_right] == 1
                x_right += 1
			end
            
            x_left = 0
            while array_filtered[start_coordinate[1], start_coordinate[2] - x_left] == 1
                x_left += 1
			end
            
            y_top = 0
            while array_filtered[start_coordinate[1] + y_top, start_coordinate[2]] == 1
                y_top += 1
			end
            
            y_bottom = 0
            while array_filtered[start_coordinate[1] - y_bottom, start_coordinate[2]] == 1
                y_bottom += 1
			end
                
            x_dist = x_right + x_left
            y_dist = y_top + y_bottom

			range1 = round(0.7 * y_dist):round(1.2 * y_dist)
            if ((x_dist in range1) == false) || ((round(0.7 * y_dist) == 0) && (round(1.2 * y_dist) == 0))
                push!(poppable_keys, key)
            else
				nothing
			end
		end
		
        for key in poppable_keys
			pop!(cal_rod_dict)
		end
                
        if length(cal_rod_dict) == 0
			nothing
        else
            append!(large_index, idx)
		end
	end
	return slice_dict, large_index
end

# ╔═╡ a548b6b9-3777-437b-8575-1f91be77e4be
slice_dict, large_index = get_calcium_slices(masked_array, header);

# ╔═╡ 6c87123f-769e-468c-b8c8-032b0a121e56
slice_dict, large_index

# ╔═╡ 16a9ac09-132a-43cb-a415-5c98493eba65
md"""
## `get_calcium_center_slices`

Returns the slices that contain the calcium vessel inserts `slice_dict` and the center of the calcium calibration rod slice `flipped_index`
"""

# ╔═╡ 3e25a510-1344-4858-b6b3-8787f0886cd2
function get_calcium_center_slices(dcm_array, slice_dict, large_index)
	flipped_index = Int(round(median(large_index)))
    edge_index = []
    if flipped_index < (size(dcm_array, 3) / 2)
        flipped = -1
        for element in large_index
            if element > (size(dcm_array, 3) / 2)
                append!(edge_index, element)
			end
		end
        if length(edge_index) == 0
            nothing
        else
            for index_edge in minimum(edge_index):size(dcm_array, 3)
                try
                    delete!(slice_dict, index_edge)
				catch
                    nothing
				end
			end
            for element2 in edge_index
				deleteat!(large_index, findall(x->x==element2, large_index))
			end
		end
                
        for element in 1:maximum(large_index)
            try
                delete!(slice_dict, element)
			catch
                nothing
			end
		end
    else
        flipped = 1
        for element in large_index
            if element < (size(dcm_array, 3) / 2)
                append!(edge_index, element)
			end
		end
        if length(edge_index) == 0
            nothing
        else
            for index_edge in 1:maximum(edge_index)
                try
                    delete!(slice_dict, index_edge)
				catch
                    nothing
				end
			end
            for element2 in edge_index
				deleteat!(large_index, findall(x->x==element2, large_index))
			end
		end
        for element in minimum(large_index):size(dcm_array, 3)
            try
                delete!(slice_dict, element)
			catch
                nothing
			end
		end
	end
	return slice_dict, flipped, flipped_index
end

# ╔═╡ da0d4727-7726-4ae5-8e61-8ed7a48e5ee7
slice_dict2, flipped, flipped_index = get_calcium_center_slices(masked_array, slice_dict, large_index)

# ╔═╡ afbe2b03-df48-45d6-92f0-f35a12614b02
md"""
## `poppable_keys`
"""

# ╔═╡ d77fd459-2207-4a61-825c-527469e5c623
function poppable_keys(flipped, flipped_index, header, slice_dict)
	SliceThickness = header[(0x0018,0x0050)]
	poppable_keys = []        
    if flipped == -1
        for key in slice_dict
            if key[1] > (flipped_index + (55 / SliceThickness))
                append!(poppable_keys, key)
			elseif flipped == 1
		        for key in slice_dict
		            if key[1] < (flipped_index - (55 / SliceThickness))
		                append!(poppable_keys, key)
					end
				end
			end
		end
	end
    for key in poppable_keys
		pop!(slice_dict)           
	end
	return slice_dict
end

# ╔═╡ 2aa8ab66-d297-4514-a8ae-56c19dfbb27c
pop_keys = poppable_keys(flipped, flipped_index, header, slice_dict)

# ╔═╡ 0802301f-864e-45a7-b39d-648b8fda1f03
md"""
## `compute_CCI`
"""

# ╔═╡ 29b34c30-a280-48bd-893b-4a3fe788fd20
function compute_CCI(dcm_array, header, slice_dict, flipped; calcium_threshold=130)
	SliceThickness = header[(0x0018,0x0050)]
	max_key, _ = maximum(zip(values(slice_dict), keys(slice_dict)))
    max_keys = []
    for key in slice_dict
        if key[2] == max_key
            append!(max_keys, key[1])
		end
	end
    slice_CCI = Int(floor(median(max_keys)))
    
    array = copy(dcm_array)
    array = Int.(array .> calcium_threshold)
    
    calcium_image = array .* dcm_array
    quality_slice = Int.(round(slice_CCI - flipped * (20 / SliceThickness)))

    cal_rod_slice = slice_CCI + (flipped * Int(30 / SliceThickness))
    
    return calcium_image, slice_CCI, quality_slice, cal_rod_slice
end

# ╔═╡ df4c740d-cdd8-468a-8699-3e37669f7094
calcium_image1, slice_CCI1, quality_slice1, cal_rod_slice1 = compute_CCI(masked_array, header, slice_dict2, flipped);

# ╔═╡ 8f91a763-3458-443a-a634-ee392ce8d130
slice_CCI1, quality_slice1, cal_rod_slice1

# ╔═╡ 7dd55970-fca1-4cd6-87ed-aa5578be12b5
@bind a PlutoUI.Slider(1:size(calcium_image1, 3), default=10, show_value=true)

# ╔═╡ 57962d59-66ac-4e88-999d-d3ec825bb8ef
begin
	fig4 = Figure()
	
	ax4 = Makie.Axis(fig4[1, 1])
	heatmap!(transpose(calcium_image1[:, :, a]), colormap=:grays)
	ax5 = Makie.Axis(fig4[1,2])
	heatmap!(transpose(dcm_array[:, :, a]), colormap=:grays)
	fig4
end

# ╔═╡ 85a3591e-60c2-4522-9e73-1be0d99845d6
md"""
## `mask_rod`
"""

# ╔═╡ d46f59c3-a25c-4363-8563-853748a3ee84
function mask_rod(dcm_array, header; calcium_threshold=130)
    slice_dict, large_index = get_calcium_slices(
		dcm_array, header; 
		calcium_threshold=calcium_threshold
	)
    slice_dict, flipped, flipped_index = get_calcium_center_slices(
		dcm_array, slice_dict, large_index
	)
    slice_dict = poppable_keys(flipped, flipped_index, header, slice_dict)
    calcium_image, slice_CCI, quality_slice, cal_rod_slice = compute_CCI(
		dcm_array, header, slice_dict, flipped; calcium_threshold=calcium_threshold
	)
	return calcium_image, slice_CCI, quality_slice, cal_rod_slice
end

# ╔═╡ 212ed8c0-6aee-448a-aef8-6db18ea94d18
calcium_image, slice_CCI, quality_slice, cal_rod_slice = mask_rod(masked_array, header);

# ╔═╡ 9e1eb1a4-5a86-4aca-9a5e-c7a75ebcc75b
slice_CCI

# ╔═╡ 7b53db36-4bb5-4130-894b-b56be069ae6d
@bind b PlutoUI.Slider(1:size(calcium_image, 3), default=10, show_value=true)

# ╔═╡ fed039bf-2ace-4471-9d44-b55b2e6b5465
heatmap(transpose(calcium_image[:, :, b]), colormap=:grays)

# ╔═╡ e52eb19a-a8cb-4667-82b8-6fcf0db418ee
heatmap(transpose(calcium_image[:, :, slice_CCI]), colormap=:grays)

# ╔═╡ b7e7d1d2-8213-4f0d-8fc0-4359c8268537
md"""
# Calcium inserts mask
"""

# ╔═╡ 6a4ee530-dc95-4732-b3b2-5171585bcab4
md"""
## `angle_calc`
"""

# ╔═╡ 005659e3-80e9-4308-8230-508bff92264f
function angle_calc(side1, side2)
    #Calculate angle between two sides of rectangular triangle
    if side1 == 0
        angle = 0
	elseif side2 == 0
        angle = π / 2
    else
        angle = atan(side1 / side2)
	end
    
    return angle
end

# ╔═╡ 99299d76-9c48-4435-8d38-0abca676144a
angle_calc(4, 3)

# ╔═╡ 7ed5df88-34fd-4638-b25e-9f805e89e22e
md"""
## `create_circular_mask`
"""

# ╔═╡ 7c590aa2-cec9-4860-8cbb-ba44ddda81c1
function create_circular_mask(h, w, center_circle, radius_circle)
	Y, X = collect(1:h), collect(1:w)'
    dist_from_center = sqrt.((X .- center_circle[1]).^2 .+ (Y .- center_circle[2]).^2)

    mask = dist_from_center .<= radius_circle
    
    return mask
end

# ╔═╡ 3cebbb7a-724e-4f2f-ac8e-3faeab27d76f
mask1 = create_circular_mask(40, 40, [20, 20], 1);

# ╔═╡ f1b6675a-9755-48ff-9e49-c32b5e8c220b
heatmap(mask1, colormap=:grays)

# ╔═╡ 9255a372-52bd-4ea2-a9ce-171ab134e1e9
md"""
## `calc_output`
"""

# ╔═╡ 72ce8586-1e3d-4da4-8e3f-d70f711e56ba
function calc_output_sim(dcm_array, header, CCI_slice, calcium_threshold=130, comp_connect=trues(3, 3))
	# Actual scoring for CCI insert
    # First step is to remove slices without calcium from arrays
	PixelSpacing = Phantoms.get_pixel_size(header)
	SliceThickness = header[(0x0018, 0x0050)]
    CCI_min = Int((CCI_slice - round(5 / SliceThickness, RoundUp)))
    CCI_max = Int((CCI_slice + round(5 / SliceThickness, RoundUp)) + 1)
    central_CCI = Int(round((CCI_max - CCI_min) / 2))
    
    if CCI_min < 0
        CCI_min = 0
	end
    if CCI_max > size(dcm_array, 3)
        CCI_max = size(dcm_array, 3)
	end
    
    CCI_array = copy(dcm_array[:, :, CCI_min+1:CCI_max])

	image_kernel = Int(round(3 / PixelSpacing[1]))
    if image_kernel % 2 == 0
        image_kernel += 1
	end
    
    CCI_array_binary = copy(CCI_array)
	CCI_array_binary = Int.(CCI_array_binary .> 1.0*calcium_threshold)
	inp = CCI_array_binary[:, :, central_CCI - 1] + CCI_array_binary[:, :, central_CCI] + CCI_array_binary[:, :, central_CCI + 1]
	components = ImageComponentAnalysis.label_components(inp, comp_connect)
	a1 = analyze_components(components, BasicMeasurement(area=true, perimeter=true))
	a2 = analyze_components(components, BoundingBox(box_area = true))
	df = leftjoin(a1, a2, on = :l)
	centroids = []
	for row in eachrow(df)
		indices = row[:box_indices]
		x_point = ((indices[1][end] - indices[1][1]) ÷ 2) + indices[1][1]
		y_point = ((indices[2][end] - indices[2][1]) ÷ 2) + indices[2][1]
		push!(centroids, (x_point, y_point))
	end
    
    centroids = deleteat!(centroids, 1)
        
	i1 = mapwindow(median, CCI_array_binary[:,:,central_CCI - 1], (image_kernel, image_kernel))
	i2 = mapwindow(median, CCI_array_binary[:,:,central_CCI], (image_kernel, image_kernel))
	i3 = mapwindow(median, CCI_array_binary[:,:,central_CCI + 1], (image_kernel, image_kernel))

	image_for_center = i1 + i2 + i3

	components2 = ImageComponentAnalysis.label_components(image_for_center, comp_connect)
	components2 = Int.(components2 .> 0)
	components2 = ImageComponentAnalysis.label_components(components2, comp_connect)
	
	b1 = analyze_components(components2, BasicMeasurement(area=true, perimeter=true))
	b2 = analyze_components(components2, BoundingBox(box_area = true))
	df2 = leftjoin(b1, b2, on = :l)
	centroids2 = []
	for row in eachrow(df2)
		indices = row[:box_indices]
		x_point = ((indices[1][end] - indices[1][1]) ÷ 2) + indices[1][1]
		y_point = ((indices[2][end] - indices[2][1]) ÷ 2) + indices[2][1]
		push!(centroids2, (y_point, x_point))
	end
	
	output = length(unique(components2)), components2, df2, centroids2
	return output
end

# ╔═╡ 9fa7d638-d399-4887-81f9-609830969fe3
output = calc_output_sim(masked_array, header, slice_CCI1);

# ╔═╡ 8de5a772-7d6b-43ee-9c34-43195c172bfe
output

# ╔═╡ 508c758c-17ce-4a02-b464-32dd6f177f47
heatmap(transpose(output[2]))

# ╔═╡ baba875d-c791-4f12-a0d2-8160ba6a823e
heatmap(transpose(masked_array[:, :, 5]), colormap=:grays)

# ╔═╡ b3f82d69-6ccd-4760-af29-9122a0386c2a
md"""
## `center_points`
"""

# ╔═╡ cfaace51-8302-44cf-b619-ecba23e0b3f8
function center_points(dcm_array, output, header, tmp_center, CCI_slice)
	PixelSpacing = Phantoms.get_pixel_size(header)
	rows, cols = Int(header[(0x0028, 0x0010)]), Int(header[(0x0028, 0x0011)])
    sizes = []
    for row in eachrow(output[3])
		area = row[:area]
		append!(sizes, area)
	end

	centroids = output[4]
    largest = Dict()
    for index in 1:length(centroids)
		x = centroids[index][1]
		y = centroids[index][2]
		dist_loc = sqrt((tmp_center[2] - x)^2 + (tmp_center[1] - y)^2)
        dist_loc *= PixelSpacing[1]
        if dist_loc > 31
            largest[index] = [round(y), round(x)]
        else
            nothing
		end
	end

    max_dict = Dict()
	radius = round(2.5 / PixelSpacing[1], RoundUp)
    for key in largest
        tmp_arr = create_circular_mask(rows, cols, (key[2][2], key[2][1]), radius)
        tmp_arr = @. abs(tmp_arr * dcm_array[:,:,CCI_slice]) + abs(tmp_arr * dcm_array[:,:,CCI_slice - 1]) + abs(tmp_arr * dcm_array[:,:,CCI_slice + 1])
        tmp_arr = @. ifelse(tmp_arr == 0, missing, tmp_arr)
        max_dict[key[1]] = median(skipmissing(tmp_arr))
	end
    large1_index, large1_key = maximum(zip(values(max_dict), keys(max_dict)))
    pop!(max_dict, large1_key)
    large2_index, large2_key = maximum(zip(values(max_dict), keys(max_dict)))
    pop!(max_dict, large2_key)
    large3_index, large3_key = maximum(zip(values(max_dict), keys(max_dict)))

    center1 = largest[large1_key]
    center2 = largest[large2_key]  
    center3 = largest[large3_key]
	
    center = find_circle(center1, center2, center3)
	return center, center1, center2, center3
end

# ╔═╡ d5b7264b-0ca0-4b4d-abcd-e597dd4422fc
center_insert

# ╔═╡ fad7cd8d-0d3d-446c-a348-f61084d02cea
heatmap(transpose(masked_array[:, :, 5]), colormap=:grays)

# ╔═╡ 683207b8-c3d9-45be-a32b-a96fe8995400
center, center1, center2, center3 = center_points(dcm_array, output, header, center_insert, slice_CCI1)

# ╔═╡ 1ed8ee2e-3c5b-4515-bb45-8fb7b9c12d82
md"""
## `calc_centers`
"""

# ╔═╡ dc1a292c-ca0f-42b6-8f4c-30adbec04d76
function calc_centers(dcm_array, output, header, tmp_center, CCI_slice)
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

# ╔═╡ 66eb4749-fd37-4b1e-9405-2e4c80018b68
dict = calc_centers(dcm_array, output, header, center_insert, slice_CCI1)

# ╔═╡ 1416f1ff-9314-4ee4-9a27-074d9961d349
dict[:Large_MD]

# ╔═╡ dfd46611-7feb-40d2-9f85-0948ae61cddf
heatmap(transpose(masked_array[:, :, 5]), colormap=:grays)

# ╔═╡ 14c1613a-02eb-4e5e-9c78-677a4c3e8f76
md"""
## `mask_inserts`
"""

# ╔═╡ a94f2710-1953-4679-8b26-11c6ba990cfa
function mask_inserts(
	dcm_array, masked_array, header, CCI_slice, center_insert; 
	calcium_threshold=130, comp_connect=trues(3, 3)
	)
	
    output = calc_output_sim(masked_array, header, CCI_slice, calcium_threshold, comp_connect)
    insert_centers = calc_centers(dcm_array, output, header, center_insert, CCI_slice)

	PixelSpacing = Phantoms.get_pixel_size(header)
	rows, cols = Int(header[(0x0028, 0x0010)]), Int(header[(0x0028, 0x0011)])
	
    mask_L_HD = create_circular_mask(cols, rows, insert_centers[:Large_HD], (round(5 / PixelSpacing[1], RoundUp) / 2) + 1)
    mask_L_MD = create_circular_mask(cols, rows, insert_centers[:Large_MD], (round(5 / PixelSpacing[1], RoundUp) / 2) + 1)
    mask_L_LD = create_circular_mask(cols, rows, insert_centers[:Large_LD], (round(5 / PixelSpacing[1], RoundUp) / 2) + 1)   
    mask_M_HD = create_circular_mask(cols, rows, insert_centers[:Medium_HD], (round(3 / PixelSpacing[1], RoundUp) / 2) + 1)
    mask_M_MD = create_circular_mask(cols, rows, insert_centers[:Medium_MD], (round(3 / PixelSpacing[1], RoundUp) / 2) + 1)
    mask_M_LD = create_circular_mask(cols, rows, insert_centers[:Medium_LD], (round(3 / PixelSpacing[1], RoundUp) / 2) + 1) 
    mask_S_HD = create_circular_mask(cols, rows, insert_centers[:Small_HD], (round(1 / PixelSpacing[1], RoundUp) / 2) + 1)
    mask_S_MD = create_circular_mask(cols, rows, insert_centers[:Small_MD], (round(1 / PixelSpacing[1], RoundUp) / 2) + 1)
    mask_S_LD = create_circular_mask(cols, rows, insert_centers[:Small_LD], (round(1 / PixelSpacing[1], RoundUp) / 2) + 1) 

    return mask_L_HD, mask_M_HD, mask_S_HD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_LD, mask_M_LD, mask_S_LD
end

# ╔═╡ 22d60146-ca07-45cb-8eae-b8cc06b35504
mask_L_HD, mask_M_HD, mask_S_HD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_LD, mask_M_LD, mask_S_LD = mask_inserts(dcm_array, masked_array, header, slice_CCI1, center_insert);

# ╔═╡ 9dd44ca1-abd7-459e-87d4-09c576716133
@bind a3 PlutoUI.Slider(1:size(masked_array, 3), default=5, show_value=true)

# ╔═╡ 517d1578-d3c5-4748-a0ed-dc9fbf787dab
heatmap(transpose(masked_array[:, :, a3]), colormap=:grays)

# ╔═╡ a98fae0f-c107-49a7-9d60-ce9e7c767ae7
masks = transpose(mask_L_HD + mask_M_HD + mask_S_HD + mask_L_MD + mask_M_MD + mask_S_MD + mask_L_LD + mask_M_LD + mask_S_LD);

# ╔═╡ 700cfd05-93f6-4133-9e3c-457b078589c7
heatmap(transpose(masks), colormap=:grays)

# ╔═╡ 1957ae7a-e63a-420b-82cd-f7148115252c
md"""
## Visualize
"""

# ╔═╡ f3d3b979-3339-47a7-baa8-39d476bc2a28
function collect_tuple(tuple_array)
	row_num = size(tuple_array)
	col_num = length(tuple_array[1])
	container = zeros(Int64, row_num..., col_num)
	for i in 1:length(tuple_array)
		container[i,:] = collect(tuple_array[i])
	end
	return container
end

# ╔═╡ dbe32d0d-8e81-4aa0-bda7-52a132dee0be
function overlay_mask_bind(mask)
	indices = findall(x -> x == 1, mask)
	indices = Tuple.(indices)
	label_array = collect_tuple(indices)
	zs = unique(label_array[:,3])
	return PlutoUI.Slider(1:length(zs), default=25, show_value=true)
end

# ╔═╡ 2fe5f75a-f335-4b1d-83aa-e3f68af95b80
function overlay_mask_plot(array, mask, var, title::AbstractString)
	indices = findall(x -> x == 1, mask)
	indices = Tuple.(indices)
	label_array = collect_tuple(indices)
	zs = unique(label_array[:,3])
	indices_lbl = findall(x -> x == zs[var], label_array[:,3])
	
	fig = Figure()
	ax = Makie.Axis(fig[1, 1])
	ax.title = title
	heatmap!((array[:, :, zs[var]]), colormap=:grays)
	scatter!(label_array[:, 1][indices_lbl], label_array[:, 2][indices_lbl], markersize=1, color=:red)
	fig
end

# ╔═╡ b28baa54-34a0-4676-baef-c0ba4422762b
begin
	masks_3D = Array{Bool}(undef, size(dcm_array))
	for z in 1:size(dcm_array, 3)
		masks_3D[:, :, z] = masks
	end
end;

# ╔═╡ fa23d5ae-7d81-45cd-a30b-a226c9166822
@bind v1 overlay_mask_bind(masks_3D)

# ╔═╡ 0a9ba745-3850-4191-a41e-bfb49400bee2
overlay_mask_plot(dcm_array, masks_3D, v1, "masks overlayed")

# ╔═╡ Cell order:
# ╠═946debe7-078f-4db6-a383-d0b5a973bd4c
# ╠═601c4b6d-ca82-4fb0-b66b-5b50d507fa83
# ╠═3a43e3e5-adc3-4a0d-ad98-1a7e39dc6a69
# ╠═20e6e9ba-2d70-4d47-9187-95b903d95814
# ╟─65b807a2-8e49-4065-baf7-dca0534a7f20
# ╠═cfb8b5b9-cfaa-4b4f-81bd-3bc7b2eb30ba
# ╟─34221566-cdd5-480a-85c1-8262dda92a38
# ╠═defe2c1c-576f-451f-b4b9-ac88411bf45e
# ╟─cf853da4-c254-4a44-abc2-27b871bba172
# ╠═2f8f212d-d897-4034-a63c-e635de77a6fc
# ╠═f8c7a77e-b37f-4c6f-801c-5e5636cb46d5
# ╠═6bbd91e2-1832-4a01-87dc-0c253452e360
# ╠═d39cd4ac-8fa4-46d6-b0ae-0d99c912b8d7
# ╠═936c2be7-072f-4774-8261-3780e7e0897e
# ╠═fd9bba1e-27a7-48d9-b44d-14c213582e07
# ╟─d1f83e63-e42a-4a1c-a5eb-c4dfe93679ac
# ╠═5c18f738-2e97-4e05-9a4c-f688895f5d83
# ╠═21ad5a93-7121-4fff-95d7-c5a6b8bff336
# ╟─9ed39437-e066-46da-94bd-0cbacbbe9981
# ╠═5f389fbc-655f-4c32-97b8-1584f5990181
# ╠═26b57188-0d11-4f4d-bed3-37f7d8b37a96
# ╠═544734b7-a7de-41ac-b2ea-9587d1409bab
# ╟─f1f2f519-2ee3-4392-94f0-6114da0c98f4
# ╟─0b40f43c-78d1-4714-a840-3788ba8867b6
# ╟─6ce0efe4-6d4f-4ff4-983b-3932d75b5dc9
# ╟─3b0263c8-280d-4081-a1ef-b24bf1a841fe
# ╠═278f03d7-6f86-432f-9d48-6574e596fad4
# ╟─01e277ce-dfa7-4b8f-bd7e-2bd5c52c8def
# ╟─20c7b163-c2c2-4dd3-9cf2-1aeca69df7c0
# ╠═8363142b-2909-4133-9a50-a43232f3e82d
# ╠═a548b6b9-3777-437b-8575-1f91be77e4be
# ╠═6c87123f-769e-468c-b8c8-032b0a121e56
# ╟─16a9ac09-132a-43cb-a415-5c98493eba65
# ╠═3e25a510-1344-4858-b6b3-8787f0886cd2
# ╠═da0d4727-7726-4ae5-8e61-8ed7a48e5ee7
# ╟─afbe2b03-df48-45d6-92f0-f35a12614b02
# ╠═d77fd459-2207-4a61-825c-527469e5c623
# ╠═2aa8ab66-d297-4514-a8ae-56c19dfbb27c
# ╟─0802301f-864e-45a7-b39d-648b8fda1f03
# ╠═29b34c30-a280-48bd-893b-4a3fe788fd20
# ╠═df4c740d-cdd8-468a-8699-3e37669f7094
# ╠═8f91a763-3458-443a-a634-ee392ce8d130
# ╟─7dd55970-fca1-4cd6-87ed-aa5578be12b5
# ╠═57962d59-66ac-4e88-999d-d3ec825bb8ef
# ╟─85a3591e-60c2-4522-9e73-1be0d99845d6
# ╠═d46f59c3-a25c-4363-8563-853748a3ee84
# ╠═212ed8c0-6aee-448a-aef8-6db18ea94d18
# ╠═9e1eb1a4-5a86-4aca-9a5e-c7a75ebcc75b
# ╟─7b53db36-4bb5-4130-894b-b56be069ae6d
# ╠═fed039bf-2ace-4471-9d44-b55b2e6b5465
# ╠═e52eb19a-a8cb-4667-82b8-6fcf0db418ee
# ╟─b7e7d1d2-8213-4f0d-8fc0-4359c8268537
# ╟─6a4ee530-dc95-4732-b3b2-5171585bcab4
# ╠═005659e3-80e9-4308-8230-508bff92264f
# ╠═99299d76-9c48-4435-8d38-0abca676144a
# ╟─7ed5df88-34fd-4638-b25e-9f805e89e22e
# ╠═7c590aa2-cec9-4860-8cbb-ba44ddda81c1
# ╠═3cebbb7a-724e-4f2f-ac8e-3faeab27d76f
# ╠═f1b6675a-9755-48ff-9e49-c32b5e8c220b
# ╟─9255a372-52bd-4ea2-a9ce-171ab134e1e9
# ╠═72ce8586-1e3d-4da4-8e3f-d70f711e56ba
# ╠═9fa7d638-d399-4887-81f9-609830969fe3
# ╠═8de5a772-7d6b-43ee-9c34-43195c172bfe
# ╠═508c758c-17ce-4a02-b464-32dd6f177f47
# ╠═baba875d-c791-4f12-a0d2-8160ba6a823e
# ╟─b3f82d69-6ccd-4760-af29-9122a0386c2a
# ╠═cfaace51-8302-44cf-b619-ecba23e0b3f8
# ╠═d5b7264b-0ca0-4b4d-abcd-e597dd4422fc
# ╠═fad7cd8d-0d3d-446c-a348-f61084d02cea
# ╠═683207b8-c3d9-45be-a32b-a96fe8995400
# ╟─1ed8ee2e-3c5b-4515-bb45-8fb7b9c12d82
# ╠═dc1a292c-ca0f-42b6-8f4c-30adbec04d76
# ╠═66eb4749-fd37-4b1e-9405-2e4c80018b68
# ╠═1416f1ff-9314-4ee4-9a27-074d9961d349
# ╠═dfd46611-7feb-40d2-9f85-0948ae61cddf
# ╟─14c1613a-02eb-4e5e-9c78-677a4c3e8f76
# ╠═a94f2710-1953-4679-8b26-11c6ba990cfa
# ╠═22d60146-ca07-45cb-8eae-b8cc06b35504
# ╠═9dd44ca1-abd7-459e-87d4-09c576716133
# ╠═517d1578-d3c5-4748-a0ed-dc9fbf787dab
# ╠═a98fae0f-c107-49a7-9d60-ce9e7c767ae7
# ╠═700cfd05-93f6-4133-9e3c-457b078589c7
# ╟─1957ae7a-e63a-420b-82cd-f7148115252c
# ╟─f3d3b979-3339-47a7-baa8-39d476bc2a28
# ╟─dbe32d0d-8e81-4aa0-bda7-52a132dee0be
# ╠═2fe5f75a-f335-4b1d-83aa-e3f68af95b80
# ╠═b28baa54-34a0-4676-baef-c0ba4422762b
# ╠═fa23d5ae-7d81-45cd-a30b-a226c9166822
# ╠═0a9ba745-3850-4191-a41e-bfb49400bee2

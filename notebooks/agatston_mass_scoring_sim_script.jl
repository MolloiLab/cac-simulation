### A Pluto.jl notebook ###
# v0.19.8

using Markdown
using InteractiveUtils

# ╔═╡ fac1755c-3668-4892-a970-27d4341dae01
# ╠═╡ show_logs = false
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
		Pkg.add("Noise")
		Pkg.add(url="https://github.com/JuliaHealth/DICOM.jl")
		Pkg.add(url="https://github.com/Dale-Black/DICOMUtils.jl")
		Pkg.add(url="https://github.com/Dale-Black/PhantomSegmentation.jl")
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
	using Noise
	using DICOM
	using DICOMUtils
	using PhantomSegmentation
	using CalciumScoring
end

# ╔═╡ 8c8a2176-b524-44af-a5b2-1fd70e4bf3ad
TableOfContents()

# ╔═╡ 5433bb5e-062e-4bfe-9ff3-b19653b2574a
function create_mask(array, mask)
	@assert size(array) == size(mask)
	idxs = findall(x -> x == true, mask)
	overlayed_mask = zeros(size(array))
	for idx in idxs
		overlayed_mask[idx] = array[idx]
	end
	return overlayed_mask
end

# ╔═╡ 879a8e4c-1002-49f1-b4d2-e471e342c387
TYPE = "agatston"

# ╔═╡ 7189b87c-b0d2-4622-bd27-ddad13694c4c
VENDERS = ["80", "100", "120", "135"]

# ╔═╡ a044c5a1-3a40-4b00-9235-6e58351afb4b
SIZES = ["small", "medium", "large"]

# ╔═╡ ea981c62-fcec-4bfc-8068-af4397b7c499
DENSITIES = ["low", "normal"]

# ╔═╡ 0e0c5770-852d-4944-b387-5c0c7faa0062
blurs = [0, 0.5, 1, 1.5, 2]

# ╔═╡ 0f880597-e049-48b3-ad48-7dc17bd6e40e
begin
	dfs = []
	for blur in blurs
		for VENDER in VENDERS
			for SIZE in SIZES
				for DENSITY in DENSITIES
					SCAN_NUMBER = 1
					BASE_PATH = string("/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/images_new/", SIZE, "/", DENSITY, "/")
					root_path = string(BASE_PATH, VENDER)
					dcm_path_list = dcm_list_builder(root_path)
					pth = dcm_path_list[SCAN_NUMBER]
					scan = basename(pth)
					header, dcm_array, slice_thick_ori1 = dcm_reader(pth)
		
					# Segment Heart
					masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3)÷2)

					# Motion Blur
					if blur != 0
						for z in size(dcm_array, 3)
							dcm_array[:, :, z] = mult_gauss(dcm_array[:, :, z], blur)
						end
					end
				
					# Segment Calcium Rod
					local thresh
					if DENSITY == "low"
						thresh = 55
					elseif DENSITY ==  "normal"
						thresh = 130
					end
					
					# # Segment Calcium Rod
					# local thresh
					# if DENSITY == "low" && SIZE == "large"
					# 	thresh = 75
					# elseif DENSITY == "low" && SIZE == "medium"
					# 	thresh = 75
					# elseif DENSITY == "low"
					# 	thresh = 60
					# elseif DENSITY ==  "normal"
					# 	thresh = 130
					# end
	
					# # Segment Calcium Rod (reproducibility1)
					# local thresh
					# if DENSITY == "low" && SIZE == "large" && VENDER == "80"
					# 	thresh = 80
					# elseif DENSITY == "low" && SIZE == "large" && VENDER == "100"
					# 	thresh = 70
					# elseif DENSITY == "low" && SIZE == "large"
					# 	thresh = 75
					# elseif DENSITY == "low" && SIZE == "medium" && VENDER == "135"
					# 	thresh = 55
					# elseif DENSITY == "low" && SIZE == "medium"
					# 	thresh = 75
					# elseif DENSITY == "low"
					# 	thresh = 60
					# elseif DENSITY ==  "normal"
					# 	thresh = 130
					# end
	
					@info DENSITY, SIZE, VENDER
					calcium_image, slice_CCI, quality_slice, cal_rod_slice = mask_rod(masked_array, header; calcium_threshold=thresh)
			
					# Segment Calcium Inserts
					# mask_L_HD, mask_M_HD, mask_S_HD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_LD, mask_M_LD, mask_S_LD = mask_inserts_simulation(
					# 		dcm_array, masked_array, header, slice_CCI, center_insert
					# )
					root_new = string("/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/julia_arrays/", SIZE, "/") 
					mask_L_HD = Array(CSV.read(string(root_new, "mask_L_HD.csv"), DataFrame; header=false))
					mask_M_HD = Array(CSV.read(string(root_new, "mask_M_HD.csv"), DataFrame; header=false))
					mask_S_HD = Array(CSV.read(string(root_new, "mask_S_HD.csv"), DataFrame; header=false))
					mask_L_MD = Array(CSV.read(string(root_new, "mask_L_MD.csv"), DataFrame;
					header=false))
					mask_M_MD = Array(CSV.read(string(root_new, "mask_M_MD.csv"), DataFrame; header=false))
					mask_S_MD = Array(CSV.read(string(root_new, "mask_S_MD.csv"), DataFrame; header=false))
					mask_L_LD = Array(CSV.read(string(root_new, "mask_L_LD.csv"), DataFrame; header=false))
					mask_M_LD = Array(CSV.read(string(root_new, "mask_M_LD.csv"), DataFrame; header=false))
					mask_S_LD = Array(CSV.read(string(root_new, "mask_S_LD.csv"), DataFrame; header=false))
				
					# Mask Calibration Factor
					output = calc_output(masked_array, header, slice_CCI, thresh, trues(3, 3))
					insert_centers = calc_centers(dcm_array, output, header, center_insert, slice_CCI)
					rows, cols = Int(header[tag"Rows"]), Int(header[tag"Columns"])
					pixel_size = DICOMUtils.get_pixel_size(header)
					mass_cal_factor, angle_0_200HA, water_rod_metrics = mass_calibration(masked_array, insert_centers[:Large_LD], center_insert, 2, cols, rows, pixel_size)
	
					local agat_thresh
					agat_thresh = 130
					# if VENDER == "80"
					# 	agat_thresh = 110
					# elseif VENDER == "100"
					# 	agat_thresh = 87
					# elseif VENDER == "120"
					# 	agat_thresh = 77
					# else VENDER == "135"
					# 	agat_thresh = 72
					# end
					
					# Score Large Inserts
					arr = masked_array[:, :, 4:6]
					
					## High Density
					mask_L_HD_3D = Array{Bool}(undef, size(arr))
					for z in 1:size(arr, 3)
						mask_L_HD_3D[:, :, z] = mask_L_HD
					end
					dilated_mask_L_HD = dilate(dilate(mask_L_HD_3D))
					alg = Agatston()
					overlayed_mask_l_hd = create_mask(arr, dilated_mask_L_HD)
					agat_l_hd, mass_l_hd = score(overlayed_mask_l_hd, pixel_size, mass_cal_factor, alg; threshold=agat_thresh)
					
					## Medium Density
					mask_L_MD_3D = Array{Bool}(undef, size(arr))
					for z in 1:size(arr, 3)
						mask_L_MD_3D[:, :, z] = mask_L_MD
					end
					dilated_mask_L_MD = dilate(dilate(mask_L_MD_3D))
					overlayed_mask_l_md = create_mask(arr, dilated_mask_L_MD)
					agat_l_md, mass_l_md = score(overlayed_mask_l_md, pixel_size, mass_cal_factor, alg; threshold=agat_thresh)
					
					## Low Density
					mask_L_LD_3D = Array{Bool}(undef, size(arr))
					for z in 1:size(arr, 3)
						mask_L_LD_3D[:, :, z] = mask_L_LD
					end
					dilated_mask_L_LD = dilate(dilate(mask_L_LD_3D))
					overlayed_mask_l_ld = create_mask(arr, dilated_mask_L_LD)
					agat_l_ld, mass_l_ld = score(overlayed_mask_l_ld, pixel_size, mass_cal_factor, alg; threshold=agat_thresh)
				
				
					# Score Medium Inserts
					## High Density
					mask_M_HD_3D = Array{Bool}(undef, size(arr))
					for z in 1:size(arr, 3)
						mask_M_HD_3D[:, :, z] = mask_M_HD
					end
					dilated_mask_M_HD = dilate(dilate(dilate(dilate(mask_M_HD_3D))))
					overlayed_mask_m_hd = create_mask(arr, dilated_mask_M_HD)
					agat_m_hd, mass_m_hd = score(overlayed_mask_m_hd, pixel_size, mass_cal_factor, alg; threshold=agat_thresh)
					
					## Medium Density
					mask_M_MD_3D = Array{Bool}(undef, size(arr))
					for z in 1:size(arr, 3)
						mask_M_MD_3D[:, :, z] = mask_M_MD
					end
					dilated_mask_M_MD = dilate(dilate(dilate(dilate(mask_M_MD_3D))))
					overlayed_mask_m_md = create_mask(arr, dilated_mask_M_MD)
					agat_m_md, mass_m_md = score(overlayed_mask_m_md, pixel_size, mass_cal_factor, alg; threshold=agat_thresh)
					
					## Low Density
					mask_M_LD_3D = Array{Bool}(undef, size(arr))
					for z in 1:size(arr, 3)
						mask_M_LD_3D[:, :, z] = mask_M_LD
					end
					dilated_mask_M_LD = dilate(dilate(dilate(dilate(mask_M_LD_3D))))
					overlayed_mask_m_ld = create_mask(arr, dilated_mask_M_LD)
					agat_m_ld, mass_m_ld = score(overlayed_mask_m_ld, pixel_size, mass_cal_factor, alg; threshold=agat_thresh)
					
					# Score Small Inserts
					## High Density
					mask_S_HD_3D = Array{Bool}(undef, size(arr))
					for z in 1:size(arr, 3)
						mask_S_HD_3D[:, :, z] = mask_S_HD
					end
					dilated_mask_S_HD = dilate((dilate(dilate(dilate((mask_S_HD_3D))))))
					overlayed_mask_s_hd = create_mask(arr, dilated_mask_S_HD)
					agat_s_hd, mass_s_hd = score(overlayed_mask_s_hd, pixel_size, mass_cal_factor, alg; threshold=agat_thresh)
					
					## Medium Density
					mask_S_MD_3D = Array{Bool}(undef, size(arr))
					for z in 1:size(arr, 3)
						mask_S_MD_3D[:, :, z] = mask_S_MD
					end
					dilated_mask_S_MD = dilate((dilate(dilate(dilate(mask_S_MD_3D)))))
					overlayed_mask_s_md = create_mask(arr, dilated_mask_S_MD)
					agat_s_md, mass_s_md = score(overlayed_mask_s_md, pixel_size, mass_cal_factor, alg; threshold=agat_thresh)
					
					## Low Density
					mask_S_LD_3D = Array{Bool}(undef, size(arr))
					for z in 1:size(arr, 3)
						mask_S_LD_3D[:, :, z] = mask_S_LD
					end
					dilated_mask_S_LD = dilate((dilate(dilate(dilate(mask_S_LD_3D)))))
					overlayed_mask_s_ld = create_mask(arr, dilated_mask_S_LD)
					agat_s_ld, mass_s_ld = score(overlayed_mask_s_ld, pixel_size, mass_cal_factor, alg; threshold=agat_thresh)
					
					# Results
					local density_array
					if DENSITY == "low"
						density_array = [0, 25, 50, 100]
					elseif DENSITY == "normal"
						density_array = [0, 200, 400, 800]
					end
					inserts = [
						"Low Density",
						"Medium Density",
						"High Density"
					]
					
					## Agatston
					calculated_agat_large = [
						agat_l_ld,
						agat_l_md,
						agat_l_hd
					]
					calculated_agat_medium = [
						agat_m_ld,
						agat_m_md,
						agat_m_hd
					]
					calculated_agat_small = [
						agat_s_ld,
						agat_s_md,
						agat_s_hd
					]
					
					## Mass
					volume_gt = [
						7.065,
						63.585,
						176.625
					]
					ground_truth_mass_large = [
						volume_gt[3] * density_array[2] * 1e-3,
						volume_gt[3] * density_array[3] * 1e-3,
						volume_gt[3] * density_array[4] * 1e-3
					] # mg
					calculated_mass_large = [
						mass_l_ld,
						mass_l_md,
						mass_l_hd
					]
					ground_truth_mass_medium = [
						volume_gt[2] * density_array[2] * 1e-3,
						volume_gt[2] * density_array[3] * 1e-3,
						volume_gt[2] * density_array[4] * 1e-3
					]
					calculated_mass_medium = [
						mass_m_ld,
						mass_m_md,
						mass_m_hd
					]
					ground_truth_mass_small = [
						volume_gt[1] * density_array[2] * 1e-3,
						volume_gt[1] * density_array[3] * 1e-3,
						volume_gt[1] * density_array[4] * 1e-3
					]
					calculated_mass_small = [
						mass_s_ld,
						mass_s_md,
						mass_s_hd
					]
				
					df = DataFrame(
						VENDER = VENDER,
						SIZE = SIZE,
						DENSITY = DENSITY,
						scan = scan,
						inserts = inserts,
						ground_truth_mass_large = ground_truth_mass_large,
						calculated_mass_large = calculated_mass_large,
						ground_truth_mass_medium = ground_truth_mass_medium,
						calculated_mass_medium = calculated_mass_medium,
						ground_truth_mass_small = ground_truth_mass_small,
						calculated_mass_small = calculated_mass_small
					)
					push!(dfs, df)
				end
			end
		end
	end
end

# ╔═╡ 8a06a521-635c-4a08-935f-b07a4ea6e0d9
md"""
# Save results
"""

# ╔═╡ 9f8b8e1a-67ac-422b-88b6-a0db53198f88
dfs

# ╔═╡ e080cc84-7712-4c7e-ac69-3c08b7c51442
if ~isdir(string(cd(pwd, "..") , "/output_new/", TYPE))
	mkdir(string(cd(pwd, "..") , "/output_new/", TYPE))
end

# ╔═╡ 73de6943-c9b1-40c4-afb5-36aa9fc665a5
begin
	new_df = vcat(dfs[1:length(dfs)]...)
	output_path_new = string(cd(pwd, "..") , "/output_new/", TYPE, "/", "full2.csv")
	CSV.write(output_path_new, new_df)
end

# ╔═╡ Cell order:
# ╠═fac1755c-3668-4892-a970-27d4341dae01
# ╠═8c8a2176-b524-44af-a5b2-1fd70e4bf3ad
# ╠═5433bb5e-062e-4bfe-9ff3-b19653b2574a
# ╠═879a8e4c-1002-49f1-b4d2-e471e342c387
# ╠═7189b87c-b0d2-4622-bd27-ddad13694c4c
# ╠═a044c5a1-3a40-4b00-9235-6e58351afb4b
# ╠═ea981c62-fcec-4bfc-8068-af4397b7c499
# ╠═0e0c5770-852d-4944-b387-5c0c7faa0062
# ╠═0f880597-e049-48b3-ad48-7dc17bd6e40e
# ╟─8a06a521-635c-4a08-935f-b07a4ea6e0d9
# ╠═9f8b8e1a-67ac-422b-88b6-a0db53198f88
# ╠═e080cc84-7712-4c7e-ac69-3c08b7c51442
# ╠═73de6943-c9b1-40c4-afb5-36aa9fc665a5

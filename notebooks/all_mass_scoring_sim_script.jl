### A Pluto.jl notebook ###
# v0.19.2

using Markdown
using InteractiveUtils

# ╔═╡ 64565c8d-da53-40a9-a480-b47082ebe66c
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

# ╔═╡ 104f4083-33d6-4f79-8a3c-93a28cbc8fbb
TableOfContents()

# ╔═╡ 967a66af-9664-43d4-a06f-ccff9aee51c3
function create_mask(array, mask)
	@assert size(array) == size(mask)
	idxs = findall(x -> x == true, mask)
	overlayed_mask = zeros(size(array))
	for idx in idxs
		overlayed_mask[idx] = array[idx]
	end
	return overlayed_mask
end

# ╔═╡ 1bb68c5a-d6ee-472b-9d4d-05b7ad7215ea
begin
	TYPE = "all"
	BASE_PATH = "/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/images/"
end

# ╔═╡ cb9a71f5-4474-4b3a-ad2c-e40cfc3b2fb1
venders = ["80", "100", "120", "135"]

# ╔═╡ c1b12ea6-5495-477b-b3bb-5a72799b17af
sizes = ["small", "medium", "large"]

# ╔═╡ d4332b45-2f1d-4ee0-bb9d-fcdd55819bd1
densities = ["low", "normal"]

# ╔═╡ de5a56a4-536c-4aed-902c-0825c7cc4aae
# begin
# 	dfs = []
# 	for vender in venders
# 		for size_ in sizes
# 			for density in densities
# 				#### ---- AGATSTON ---- ####
# 				SCAN_NUMBER = 1
# 				VENDER = vender
# 				root_path = string(BASE_PATH, "/", size_, "/", density, "/", VENDER)
# 				dcm_path_list = dcm_list_builder(root_path)
# 				pth = dcm_path_list[SCAN_NUMBER]
# 				scan = basename(pth)
# 				header, dcm_array, slice_thick_ori1 = dcm_reader(pth)
# 				for z in size(dcm_array, 3)
# 					dcm_array[:, :, z] = mult_gauss(dcm_array[:, :, z], kern)
# 				end
	
# 				# Segment Heart
# 				masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3)÷2)
			
# 				# Segment Calcium Rod
# 				calcium_image, slice_CCI, quality_slice, cal_rod_slice = mask_rod(masked_array, header)
		
# 				# Segment Calcium Inserts
# 				mask_L_HD, mask_M_HD, mask_S_HD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_LD, mask_M_LD, mask_S_LD = mask_inserts_simulation(
# 						dcm_array, masked_array, header, slice_CCI, center_insert
# 				)
			
# 				# Mask Calibration Factor
# 				output = calc_output(masked_array, header, slice_CCI, 130, trues(3, 3))
# 				insert_centers = calc_centers(dcm_array, output, header, center_insert, slice_CCI)
# 				rows, cols = Int(header[tag"Rows"]), Int(header[tag"Columns"])
# 				pixel_size = DICOMUtils.get_pixel_size(header)
# 				mass_cal_factor, angle_0_200HA, water_rod_metrics = mass_calibration(masked_array, insert_centers[:Large_LD], center_insert, 2, cols, rows, pixel_size)
				
# 				# Score Large Inserts
# 				arr = masked_array[:, :, 4:6]
				
# 				## High Density
# 				mask_L_HD_3D = Array{Bool}(undef, size(arr))
# 				for z in 1:size(arr, 3)
# 					mask_L_HD_3D[:, :, z] = mask_L_HD
# 				end
# 				dilated_mask_L_HD = dilate(dilate(mask_L_HD_3D))
# 				alg = Agatston()
# 				overlayed_mask_l_hd = create_mask(arr, dilated_mask_L_HD)
# 				agat_l_hd, mass_l_hd = score(overlayed_mask_l_hd, pixel_size, mass_cal_factor, alg)
				
# 				## Medium Density
# 				mask_L_MD_3D = Array{Bool}(undef, size(arr))
# 				for z in 1:size(arr, 3)
# 					mask_L_MD_3D[:, :, z] = mask_L_MD
# 				end
# 				dilated_mask_L_MD = dilate(dilate(mask_L_MD_3D))
# 				overlayed_mask_l_md = create_mask(arr, dilated_mask_L_MD)
# 				agat_l_md, mass_l_md = score(overlayed_mask_l_md, pixel_size, mass_cal_factor, alg)
				
# 				## Low Density
# 				mask_L_LD_3D = Array{Bool}(undef, size(arr))
# 				for z in 1:size(arr, 3)
# 					mask_L_LD_3D[:, :, z] = mask_L_LD
# 				end
# 				dilated_mask_L_LD = dilate(dilate(mask_L_LD_3D))
# 				overlayed_mask_l_ld = create_mask(arr, dilated_mask_L_LD)
# 				agat_l_ld, mass_l_ld = score(overlayed_mask_l_ld, pixel_size, mass_cal_factor, alg)
			
			
# 				# Score Medium Inserts
# 				## High Density
# 				mask_M_HD_3D = Array{Bool}(undef, size(arr))
# 				for z in 1:size(arr, 3)
# 					mask_M_HD_3D[:, :, z] = mask_M_HD
# 				end
# 				dilated_mask_M_HD = dilate(dilate(dilate(dilate(mask_M_HD_3D))))
# 				overlayed_mask_m_hd = create_mask(arr, dilated_mask_M_HD)
# 				agat_m_hd, mass_m_hd = score(overlayed_mask_m_hd, pixel_size, mass_cal_factor, alg)
				
# 				## Medium Density
# 				mask_M_MD_3D = Array{Bool}(undef, size(arr))
# 				for z in 1:size(arr, 3)
# 					mask_M_MD_3D[:, :, z] = mask_M_MD
# 				end
# 				dilated_mask_M_MD = dilate(dilate(dilate(dilate(mask_M_MD_3D))))
# 				overlayed_mask_m_md = create_mask(arr, dilated_mask_M_MD)
# 				agat_m_md, mass_m_md = score(overlayed_mask_m_md, pixel_size, mass_cal_factor, alg)
				
# 				## Low Density
# 				mask_M_LD_3D = Array{Bool}(undef, size(arr))
# 				for z in 1:size(arr, 3)
# 					mask_M_LD_3D[:, :, z] = mask_M_LD
# 				end
# 				dilated_mask_M_LD = dilate(dilate(dilate(dilate(mask_M_LD_3D))))
# 				overlayed_mask_m_ld = create_mask(arr, dilated_mask_M_LD)
# 				agat_m_ld, mass_m_ld = score(overlayed_mask_m_ld, pixel_size, mass_cal_factor, alg)
				
# 				# Score Small Inserts
# 				## High Density
# 				mask_S_HD_3D = Array{Bool}(undef, size(arr))
# 				for z in 1:size(arr, 3)
# 					mask_S_HD_3D[:, :, z] = mask_S_HD
# 				end
# 				dilated_mask_S_HD = dilate((dilate(dilate(dilate((mask_S_HD_3D))))))
# 				overlayed_mask_s_hd = create_mask(arr, dilated_mask_S_HD)
# 				agat_s_hd, mass_s_hd = score(overlayed_mask_s_hd, pixel_size, mass_cal_factor, alg)
				
# 				## Medium Density
# 				mask_S_MD_3D = Array{Bool}(undef, size(arr))
# 				for z in 1:size(arr, 3)
# 					mask_S_MD_3D[:, :, z] = mask_S_MD
# 				end
# 				dilated_mask_S_MD = dilate((dilate(dilate(dilate(mask_S_MD_3D)))))
# 				overlayed_mask_s_md = create_mask(arr, dilated_mask_S_MD)
# 				agat_s_md, mass_s_md = score(overlayed_mask_s_md, pixel_size, mass_cal_factor, alg)
				
# 				## Low Density
# 				mask_S_LD_3D = Array{Bool}(undef, size(arr))
# 				for z in 1:size(arr, 3)
# 					mask_S_LD_3D[:, :, z] = mask_S_LD
# 				end
# 				dilated_mask_S_LD = dilate((dilate(dilate(dilate(mask_S_LD_3D)))))
# 				overlayed_mask_s_ld = create_mask(arr, dilated_mask_S_LD)
# 				agat_s_ld, mass_s_ld = score(overlayed_mask_s_ld, pixel_size, mass_cal_factor, alg)
				
# 				# Results
# 				density_array = [0, 200, 400, 800]
# 				inserts = [
# 					"Low Density",
# 					"Medium Density",
# 					"High Density"
# 				]
				
# 				## Agatston
# 				calculated_agat_large = [
# 					agat_l_ld,
# 					agat_l_md,
# 					agat_l_hd
# 				]
# 				calculated_agat_medium = [
# 					agat_m_ld,
# 					agat_m_md,
# 					agat_m_hd
# 				]
# 				calculated_agat_small = [
# 					agat_s_ld,
# 					agat_s_md,
# 					agat_s_hd
# 				]
				
# 				## Mass
# 				volume_gt = [
# 					7.065,
# 					63.585,
# 					176.625
# 				]
# 				ground_truth_mass_large = [
# 					volume_gt[3] * density_array[2] * 1e-3,
# 					volume_gt[3] * density_array[3] * 1e-3,
# 					volume_gt[3] * density_array[4] * 1e-3
# 				] # mg
# 				calculated_mass_large = [
# 					mass_l_ld,
# 					mass_l_md,
# 					mass_l_hd
# 				]
# 				ground_truth_mass_medium = [
# 					volume_gt[2] * density_array[2] * 1e-3,
# 					volume_gt[2] * density_array[3] * 1e-3,
# 					volume_gt[2] * density_array[4] * 1e-3
# 				]
# 				calculated_mass_medium = [
# 					mass_m_ld,
# 					mass_m_md,
# 					mass_m_hd
# 				]
# 				ground_truth_mass_small = [
# 					volume_gt[1] * density_array[2] * 1e-3,
# 					volume_gt[1] * density_array[3] * 1e-3,
# 					volume_gt[1] * density_array[4] * 1e-3
# 				]
# 				calculated_mass_small = [
# 					mass_s_ld,
# 					mass_s_md,
# 					mass_s_hd
# 				]
			
# 				df = DataFrame(
# 					kern = kern,
# 					size = size_,
# 					density = density,
# 					scan = scan,
# 					inserts = inserts,
# 					ground_truth_mass_large = ground_truth_mass_large,
# 					calculated_mass_large = calculated_mass_large,
# 					ground_truth_mass_medium = ground_truth_mass_medium,
# 					calculated_mass_medium = calculated_mass_medium,
# 					ground_truth_mass_small = ground_truth_mass_small,
# 					calculated_mass_small = calculated_mass_small
# 				)
# 				push!(dfs, df)

# 				#### ---- INTEGRATED ---- ####
# 				SCAN_NUMBER = 1
# 				VENDER = vender
# 				root_path = string(BASE_PATH, "/", size_, "/", density, "/", VENDER)
# 				dcm_path_list = dcm_list_builder(root_path)
# 				pth = dcm_path_list[SCAN_NUMBER]
# 				scan = basename(pth)
# 				header, dcm_array, slice_thick_ori1 = dcm_reader(pth)
# 				@info header
# 				for z in size(dcm_array, 3)
# 					dcm_array[:, :, z] = Noise.mult_gauss(dcm_array[:, :, z], kern)
# 				end
			
# 				# Segment Heart
# 				masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3)÷2)
			
# 				# Segment Calcium Rod
# 				calcium_image, slice_CCI, quality_slice, cal_rod_slice = mask_rod(masked_array, header)
			
# 				# Calibration Prep
# 				array_filtered = abs.(mapwindow(median, calcium_image[:, :, 2], (3, 3)))
# 				bool_arr = array_filtered .> 0
# 				bool_arr_erode = (((erode(erode(bool_arr)))))
# 				c_img = calcium_image[:, :, 1:3]
# 				mask_cal_3D = Array{Bool}(undef, size(c_img))
# 				for z in 1:size(c_img, 3)
# 					mask_cal_3D[:, :, z] = bool_arr_erode
# 				end
# 				cal_insert_mean = quantile!(c_img[mask_cal_3D], 0.7)
# 				# eroded_mask_L_HD = erode(erode(mask_L_HD_3D))
# 				# high_density_cal = mean(arr[eroded_mask_L_HD])
# 				# eroded_mask_L_MD = erode(erode(mask_L_MD_3D))
# 				# med_density_cal = mean(arr[eroded_mask_L_MD])
# 				# eroded_mask_L_LD = erode(erode(mask_L_LD_3D))
# 				# low_density_cal = mean(arr[eroded_mask_L_LD])
# 				density_array = [0, 200, 400, 800]
# 				density_array_cal = [0, 200]
# 				# intensity_array = [0, low_density_cal, med_density_cal, high_density_cal] # HU
# 				intensity_array_cal = [0, cal_insert_mean]
# 				df_cal = DataFrame(:density => density_array_cal, :intensity => intensity_array_cal)
# 				linearRegressor = lm(@formula(intensity ~ density), df_cal);
# 				linearFit = predict(linearRegressor)
# 				m = linearRegressor.model.pp.beta0[2]
# 				b = linearRegressor.model.rr.mu[1]
# 				density(intensity) = (intensity - b) / m
# 				intensity(ρ) = m*ρ + b
	
# 				mask_L_HD, mask_M_HD, mask_S_HD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_LD, mask_M_LD, mask_S_LD = mask_inserts_simulation(
# 					dcm_array, masked_array, header, slice_CCI, center_insert)
			
# 				arr = masked_array[:, :, 4:6]
# 				single_arr = masked_array[:, :, slice_CCI]
# 				pixel_size = DICOMUtils.get_pixel_size(header)
				
# 				# Score Large InsertS
# 				## High Density
# 				mask_L_HD_3D = Array{Bool}(undef, size(arr))
# 				for z in 1:size(arr, 3)
# 					mask_L_HD_3D[:, :, z] = mask_L_HD
# 				end
# 				dilated_mask_L_HD = dilate(dilate(mask_L_HD_3D))
# 				ring_mask_L_HD = dilate(dilate(dilate(dilate(mask_L_HD_3D)))) - dilate(dilate(dilate(mask_L_HD_3D)))
# 				single_ring_mask_L_HD = Bool.(ring_mask_L_HD[:, :, 3])
# 				s_bkg_L_HD = mean(single_arr[single_ring_mask_L_HD])
# 				S_Obj_HD = intensity(800)
# 				ρ_hd = 0.8 # mg/mm^3
# 				alg_L_HD = Integrated(arr[mask_L_HD_3D])
# 				mass_l_hd = score(s_bkg_L_HD, S_Obj_HD, pixel_size, ρ_hd, alg_L_HD)
			
# 				## Medium Density
# 				mask_L_MD_3D = Array{Bool}(undef, size(arr))
# 				for z in 1:size(arr, 3)
# 					mask_L_MD_3D[:, :, z] = mask_L_MD
# 				end
# 				dilated_mask_L_MD = dilate(dilate(mask_L_MD_3D))
# 				ring_mask_L_MD = dilate(dilate(dilate(dilate(mask_L_MD_3D)))) - dilate(dilate(dilate(mask_L_MD_3D)))
# 				single_ring_mask_L_MD = Bool.(ring_mask_L_MD[:, :, 3])
# 				s_bkg_L_MD = mean(single_arr[single_ring_mask_L_MD])
# 				S_Obj_MD = intensity(400)
# 				ρ_md = 0.4 # mg/mm^3
# 				alg_L_MD = Integrated(arr[mask_L_MD_3D])
# 				mass_l_md = score(s_bkg_L_MD, S_Obj_MD, pixel_size, ρ_md, alg_L_MD)
			
# 				## Low Density
# 				mask_L_LD_3D = Array{Bool}(undef, size(arr))
# 				for z in 1:size(arr, 3)
# 					mask_L_LD_3D[:, :, z] = mask_L_LD
# 				end
# 				dilated_mask_L_LD = dilate(dilate(mask_L_LD_3D))
# 				ring_mask_L_LD = dilate(dilate(dilate(dilate(mask_L_LD_3D)))) - dilate(dilate(dilate(mask_L_LD_3D)))
# 				single_ring_mask_L_LD = Bool.(ring_mask_L_LD[:, :, 3])
# 				s_bkg_L_LD = mean(single_arr[single_ring_mask_L_LD])
# 				S_Obj_LD = intensity(200)
# 				ρ_ld = 0.2 # mg/mm^3
# 				alg_L_LD = Integrated(arr[mask_L_LD_3D])
# 				mass_l_ld = score(s_bkg_L_LD, S_Obj_LD, pixel_size, ρ_ld, alg_L_LD)
				
# 				# Score Medium Inserts
# 				## High Density
# 				mask_M_HD_3D = Array{Bool}(undef, size(arr))
# 				for z in 1:size(arr, 3)
# 					mask_M_HD_3D[:, :, z] = mask_M_HD
# 				end
# 				dilated_mask_M_HD = dilate(dilate(dilate(dilate(mask_M_HD_3D))))
# 				ring_mask_M_HD = dilate(dilate(dilate(dilate(dilate(mask_M_HD_3D))))) - dilate(dilate(dilate(dilate(mask_M_HD_3D))))
# 				single_ring_mask_M_HD = Bool.(ring_mask_M_HD[:, :, 3])
# 				s_bkg_M_HD = mean(single_arr[single_ring_mask_M_HD])
# 				alg_M_HD = Integrated(arr[mask_M_HD_3D])
# 				mass_m_hd = score(s_bkg_M_HD, S_Obj_HD, pixel_size, ρ_hd, alg_M_HD)
				
# 				## Medium Density
# 				mask_M_MD_3D = Array{Bool}(undef, size(arr))
# 				for z in 1:size(arr, 3)
# 					mask_M_MD_3D[:, :, z] = mask_M_MD
# 				end
# 				dilated_mask_M_MD = dilate(dilate(dilate(dilate(mask_M_MD_3D))))
# 				ring_mask_M_MD = dilate(dilate(dilate(dilate(dilate(dilate(mask_M_MD_3D)))))) - dilate(dilate(dilate(dilate(dilate(mask_M_MD_3D)))))
# 				single_ring_mask_M_MD = Bool.(ring_mask_M_MD[:, :, 3])
# 				s_bkg_M_MD = mean(single_arr[single_ring_mask_M_MD])
# 				alg_M_MD = Integrated(arr[mask_M_MD_3D])
# 				mass_m_md = score(s_bkg_M_MD, S_Obj_MD, pixel_size, ρ_md, alg_M_MD)
			
# 				## Low Density
# 				mask_M_LD_3D = Array{Bool}(undef, size(arr))
# 				for z in 1:size(arr, 3)
# 					mask_M_LD_3D[:, :, z] = mask_M_LD
# 				end
# 				dilated_mask_M_LD = dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))))
# 				ring_mask_M_LD = dilate(dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))))) - dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))))
# 				single_ring_mask_M_LD = Bool.(ring_mask_M_LD[:, :, 3])
# 				s_bkg_M_LD = mean(single_arr[single_ring_mask_M_LD])
# 				alg_M_LD = Integrated(arr[mask_M_LD_3D])
# 				mass_m_ld = score(s_bkg_M_LD, S_Obj_LD, pixel_size, ρ_ld, alg_M_LD)
				
			
# 				# Score Small Inserts
# 				## High Density
# 				mask_S_HD_3D = Array{Bool}(undef, size(arr))
# 				for z in 1:size(arr, 3)
# 					mask_S_HD_3D[:, :, z] = mask_S_HD
# 				end
# 				dilated_mask_S_HD = dilate(dilate(dilate(dilate(dilate((mask_S_HD_3D))))))
# 				ring_mask_S_HD = dilate(dilate(dilate(dilate(dilate(mask_S_HD_3D))))) - dilate(dilate(dilate(dilate(mask_S_HD_3D))))
# 				single_ring_mask_S_HD = Bool.(ring_mask_S_HD[:, :, 3])
# 				s_bkg_S_HD = mean(single_arr[single_ring_mask_S_HD])
# 				alg_S_HD = Integrated(arr[mask_S_HD_3D])
# 				mass_s_hd = score(s_bkg_S_HD, S_Obj_HD, pixel_size, ρ_hd, alg_S_HD)
			
# 				## Medium Density
# 				mask_S_MD_3D = Array{Bool}(undef, size(arr))
# 				for z in 1:size(arr, 3)
# 					mask_S_MD_3D[:, :, z] = mask_S_MD
# 				end
# 				dilated_mask_S_MD = dilate(dilate(dilate(dilate(dilate(mask_S_MD_3D)))))
# 				ring_mask_S_MD = dilate(dilate(dilate(dilate(dilate(mask_S_MD_3D))))) - dilate(dilate(dilate(dilate(mask_S_MD_3D))))
# 				single_ring_mask_S_MD = Bool.(ring_mask_S_MD[:, :, 3])
# 				s_bkg_S_MD = mean(single_arr[single_ring_mask_S_MD])
# 				alg_S_MD = Integrated(arr[mask_S_MD_3D])
# 				mass_s_md = score(s_bkg_S_MD, S_Obj_MD, pixel_size, ρ_md, alg_S_MD)
				
# 				## Low Density
# 				mask_S_LD_3D = Array{Bool}(undef, size(arr))
# 				for z in 1:size(arr, 3)
# 					mask_S_LD_3D[:, :, z] = mask_S_LD
# 				end
# 				dilated_mask_S_LD = dilate(dilate(dilate(dilate(dilate(mask_S_LD_3D)))))
# 				ring_mask_S_LD = dilate(dilate(dilate(dilate(dilate(mask_S_LD_3D))))) - dilate(dilate(dilate(dilate(mask_S_LD_3D))))
# 				single_ring_mask_S_LD = Bool.(ring_mask_S_LD[:, :, 3])
# 				s_bkg_S_LD = mean(single_arr[single_ring_mask_S_LD])
# 				alg_S_LD = Integrated(arr[mask_S_LD_3D])
# 				mass_s_ld = score(s_bkg_S_LD, S_Obj_LD, pixel_size, ρ_ld, alg_S_LD)
				
# 				# Results
# 				density_array = [0, 200, 400, 800]
# 				inserts = [
# 					"Low Density",
# 					"Medium Density",
# 					"High Density"
# 				]
# 				volume_gt = [
# 					7.065,
# 					63.585,
# 					176.625
# 				]
# 				ground_truth_mass_large = [
# 					volume_gt[3] * density_array[2] * 1e-3,
# 					volume_gt[3] * density_array[3] * 1e-3,
# 					volume_gt[3] * density_array[4] * 1e-3
# 				] # mg
				
# 				calculated_mass_large = [
# 					mass_l_ld,
# 					mass_l_md,
# 					mass_l_hd
# 				]
# 				ground_truth_mass_medium = [
# 					volume_gt[2] * density_array[2] * 1e-3,
# 					volume_gt[2] * density_array[3] * 1e-3,
# 					volume_gt[2] * density_array[4] * 1e-3
# 				]
# 				calculated_mass_medium = [
# 					mass_m_ld,
# 					mass_m_md,
# 					mass_m_hd
# 				]
# 				ground_truth_mass_small = [
# 					volume_gt[1] * density_array[2] * 1e-3,
# 					volume_gt[1] * density_array[3] * 1e-3,
# 					volume_gt[1] * density_array[4] * 1e-3
# 				]
# 				calculated_mass_small = [
# 					mass_s_ld,
# 					mass_s_md,
# 					mass_s_hd
# 				]
				
# 				df = DataFrame(
# 					kern = kern,
# 					size = size_,
# 					density = density,
# 					scan = scan,
# 					inserts = inserts,
# 					ground_truth_mass_large = ground_truth_mass_large,
# 					calculated_mass_large = calculated_mass_large,
# 					ground_truth_mass_medium = ground_truth_mass_medium,
# 					calculated_mass_medium = calculated_mass_medium,
# 					ground_truth_mass_small = ground_truth_mass_small,
# 					calculated_mass_small = calculated_mass_small
# 				)
# 				push!(dfs, df)
# 			end
# 		end
# 	end
# end

# ╔═╡ b5e3929b-1bf8-46e4-9210-fc2250ea11bd
md"""
# Save results
"""

# ╔═╡ 15264a36-9e72-4266-945b-a0999f031675
# dfs

# ╔═╡ b2b13afb-8886-4875-977f-95508cd0fef0
# if ~isdir(string(cd(pwd, "..") , "/output/", TYPE))
# 	mkdir(string(cd(pwd, "..") , "/output/", TYPE))
# end

# ╔═╡ 45390d48-4a19-48a7-84e4-05b5cb7ccb34
# if length(dfs) == 12
# 	new_df = vcat(dfs[1:12]...)
# 	output_path_new = string(cd(pwd, "..") , "/output/", TYPE, "/", "full.csv")
# 	CSV.write(output_path_new, new_df)
# end

# ╔═╡ Cell order:
# ╠═64565c8d-da53-40a9-a480-b47082ebe66c
# ╠═104f4083-33d6-4f79-8a3c-93a28cbc8fbb
# ╟─967a66af-9664-43d4-a06f-ccff9aee51c3
# ╠═1bb68c5a-d6ee-472b-9d4d-05b7ad7215ea
# ╠═cb9a71f5-4474-4b3a-ad2c-e40cfc3b2fb1
# ╠═c1b12ea6-5495-477b-b3bb-5a72799b17af
# ╠═d4332b45-2f1d-4ee0-bb9d-fcdd55819bd1
# ╠═de5a56a4-536c-4aed-902c-0825c7cc4aae
# ╟─b5e3929b-1bf8-46e4-9210-fc2250ea11bd
# ╠═15264a36-9e72-4266-945b-a0999f031675
# ╠═b2b13afb-8886-4875-977f-95508cd0fef0
# ╠═45390d48-4a19-48a7-84e4-05b5cb7ccb34

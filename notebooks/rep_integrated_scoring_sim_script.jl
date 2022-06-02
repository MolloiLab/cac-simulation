### A Pluto.jl notebook ###
# v0.19.4

using Markdown
using InteractiveUtils

# ╔═╡ 70410feb-aa21-405d-8d26-9cb7c1ecdb9e
# ╠═╡ show_logs = false
begin
	let
		using Pkg
		Pkg.activate(mktempdir())
		Pkg.Registry.update()
		Pkg.add("PlutoUI")
		Pkg.add("Statistics")
		Pkg.add("StatsBase")
		Pkg.add("ImageMorphology")
		Pkg.add("ImageFiltering")
		Pkg.add("CSV")
		Pkg.add("DataFrames")
		Pkg.add("GLM")
		Pkg.add("Noise")
		Pkg.add(url="https://github.com/JuliaHealth/DICOM.jl")
		Pkg.add(url="https://github.com/Dale-Black/DICOMUtils.jl")
		Pkg.add(url="https://github.com/Dale-Black/PhantomSegmentation.jl")
		Pkg.add(url="https://github.com/Dale-Black/CalciumScoring.jl")
	end
	
	using PlutoUI
	using Statistics
	using StatsBase: quantile!
	using ImageMorphology
	using ImageFiltering
	using CSV
	using DataFrames
	using GLM
	using Noise
	using DICOM
	using DICOMUtils
	using PhantomSegmentation
	using CalciumScoring
end

# ╔═╡ f2845431-56a1-4fa4-a31f-a515c686d85a
TableOfContents()

# ╔═╡ 71b9dc8f-8879-43c3-a94b-c3efbeb672a6
function create_mask(array, mask)
	@assert size(array) == size(mask)
	idxs = findall(x -> x == true, mask)
	overlayed_mask = zeros(size(array))
	for idx in idxs
		overlayed_mask[idx] = array[idx]
	end
	return overlayed_mask
end

# ╔═╡ 26110c25-c288-4e5a-a3d1-260ca0c15a84
TYPE = "integrated_scoring"

# ╔═╡ 442943d0-9a1e-46ff-a924-91b5537a9b6a
VENDERS = ["80", "100", "120", "135"]

# ╔═╡ 57656cb2-146a-43ec-abbb-1b5a6a1f3afb
SIZES = ["small", "medium", "large"]

# ╔═╡ 6e98101e-f052-4b19-9915-429166605ced
DENSITIES = ["low", "normal"]

# ╔═╡ 40454c15-d8fa-4530-8b74-39a934afc257
begin
	dfs = []
	for VENDER in VENDERS
		for SIZE in SIZES
			for DENSITY in DENSITIES
				SCAN_NUMBER = 1
				BASE_PATH = string("/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/images_reproducibility1/", SIZE, "/", DENSITY, "/")
				root_path = string(BASE_PATH, VENDER)
				dcm_path_list = dcm_list_builder(root_path)
				pth = dcm_path_list[SCAN_NUMBER]
				scan = basename(pth)
				header, dcm_array, slice_thick_ori1 = dcm_reader(pth)
				
				# Segment Heart
				masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3)÷2)
			
				# Segment Calcium Rod
				local thresh
				if DENSITY == "low" && SIZE == "large" && VENDER == "80"
					thresh = 80
				elseif DENSITY == "low" && SIZE == "large" && VENDER == "100"
					thresh = 70
				elseif DENSITY == "low" && SIZE == "large"
					thresh = 75
				elseif DENSITY == "low" && SIZE == "medium" && VENDER == "135"
					thresh = 55
				elseif DENSITY == "low" && SIZE == "medium"
					thresh = 75
				elseif DENSITY == "low"
					thresh = 60
				elseif DENSITY ==  "normal"
					thresh = 130
				end

				@info DENSITY, SIZE, VENDER, thresh
				
				calcium_image, slice_CCI, quality_slice, cal_rod_slice = mask_rod(masked_array, header; calcium_threshold=thresh)
	
				# mask_L_HD, mask_M_HD, mask_S_HD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_LD, mask_M_LD, mask_S_LD = mask_inserts_simulation(
				# 	dcm_array, masked_array, header, slice_CCI, center_insert)

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

				# Calibration Prep
				local density_array
				if DENSITY == "low"
					density_array = [0, 25, 50, 100]
				elseif DENSITY == "normal"
					density_array = [0, 200, 400, 800]
				end
				array_filtered = abs.(mapwindow(median, calcium_image[:, :, 2], (3, 3)))
				bool_arr = array_filtered .> 0
				bool_arr_erode = (((erode(erode(bool_arr)))))
				c_img = calcium_image[:, :, 1:3]
				mask_cal_3D = Array{Bool}(undef, size(c_img))
				for z in 1:size(c_img, 3)
					mask_cal_3D[:, :, z] = bool_arr_erode
				end

				arr = masked_array[:, :, 4:6]
				single_arr = masked_array[:, :, 5]
				pixel_size = DICOMUtils.get_pixel_size(header)

				mask_L_HD_3D = Array{Bool}(undef, size(arr))
				for z in 1:size(arr, 3)
					mask_L_HD_3D[:, :, z] = mask_L_HD
				end

				mask_L_MD_3D = Array{Bool}(undef, size(arr))
				for z in 1:size(arr, 3)
					mask_L_MD_3D[:, :, z] = mask_L_MD
				end

				mask_L_LD_3D = Array{Bool}(undef, size(arr))
				for z in 1:size(arr, 3)
					mask_L_LD_3D[:, :, z] = mask_L_LD
				end
				
				cal_insert_mean = mean(c_img[mask_cal_3D])
				eroded_mask_L_HD = erode(erode(mask_L_HD_3D))
				high_density_cal = mean(arr[eroded_mask_L_HD])
				eroded_mask_L_MD = erode(erode(mask_L_MD_3D))
				med_density_cal = mean(arr[eroded_mask_L_MD])
				eroded_mask_L_LD = erode(erode(mask_L_LD_3D))
				low_density_cal = mean(arr[eroded_mask_L_LD])
				density_array_cal = [0, 200]
				intensity_array = [0, low_density_cal, med_density_cal, high_density_cal] # HU
				intensity_array_cal = [0, cal_insert_mean]
				local df_cal
				if DENSITY == "low"
					df_cal = DataFrame(:density => density_array[2:end], :intensity => intensity_array[2:end])
				elseif DENSITY == "normal"
					df_cal = DataFrame(:density => density_array, :intensity => intensity_array)
				end
				linearRegressor = lm(@formula(intensity ~ density), df_cal);
				linearFit = predict(linearRegressor)
				m = linearRegressor.model.pp.beta0[2]
				b = linearRegressor.model.rr.mu[1]
				density(intensity) = (intensity - b) / m
				intensity(ρ) = m*ρ + b


				# Score background
				background_mask = zeros(size(arr)...)
				background_mask[center_insert[1]-5:center_insert[1]+5, center_insert[2]-5:center_insert[2]+5, 2] .= 1
				background_mask = Bool.(background_mask)
				background_mask_dil = dilate(dilate(background_mask))

				single_bkg_center = Bool.(background_mask[:, :, 2])
				S_Obj_bkg = intensity(200)

				ring_background = background_mask_dil - background_mask
				single_ring_bkg = Bool.(ring_background[:, :, 2])
				s_bkg = mean(single_arr[single_ring_bkg])

				alg_bkg = Integrated(arr[background_mask])
				ρ_bkg = 0.2 # mg/mm^3
				mass_bkg = score(s_bkg, S_Obj_bkg, pixel_size, ρ_bkg, alg_bkg)
				
				# Score Large InsertS
				## High Density
				
				dilated_mask_L_HD = dilate(dilate(mask_L_HD_3D))
				ring_mask_L_HD = dilate(dilate(dilate(dilate(mask_L_HD_3D)))) - dilate(dilate(dilate(mask_L_HD_3D)))
				single_ring_mask_L_HD = Bool.(ring_mask_L_HD[:, :, 3])
				s_bkg_L_HD = mean(single_arr[single_ring_mask_L_HD])
				S_Obj_HD = intensity(800)
				ρ_hd = 0.8 # mg/mm^3
				alg_L_HD = Integrated(arr[mask_L_HD_3D])
				mass_l_hd = score(s_bkg_L_HD, S_Obj_HD, pixel_size, ρ_hd, alg_L_HD)
			
				## Medium Density
				
				dilated_mask_L_MD = dilate(dilate(mask_L_MD_3D))
				ring_mask_L_MD = dilate(dilate(dilate(dilate(mask_L_MD_3D)))) - dilate(dilate(dilate(mask_L_MD_3D)))
				single_ring_mask_L_MD = Bool.(ring_mask_L_MD[:, :, 3])
				s_bkg_L_MD = mean(single_arr[single_ring_mask_L_MD])
				S_Obj_MD = intensity(400)
				ρ_md = 0.4 # mg/mm^3
				alg_L_MD = Integrated(arr[mask_L_MD_3D])
				mass_l_md = score(s_bkg_L_MD, S_Obj_MD, pixel_size, ρ_md, alg_L_MD)
			
				## Low Density
				
				dilated_mask_L_LD = dilate(dilate(mask_L_LD_3D))
				ring_mask_L_LD = dilate(dilate(dilate(dilate(mask_L_LD_3D)))) - dilate(dilate(dilate(mask_L_LD_3D)))
				single_ring_mask_L_LD = Bool.(ring_mask_L_LD[:, :, 3])
				s_bkg_L_LD = mean(single_arr[single_ring_mask_L_LD])
				S_Obj_LD = intensity(200)
				ρ_ld = 0.2 # mg/mm^3
				alg_L_LD = Integrated(arr[mask_L_LD_3D])
				mass_l_ld = score(s_bkg_L_LD, S_Obj_LD, pixel_size, ρ_ld, alg_L_LD)
				
				# Score Medium Inserts
				## High Density
				mask_M_HD_3D = Array{Bool}(undef, size(arr))
				for z in 1:size(arr, 3)
					mask_M_HD_3D[:, :, z] = mask_M_HD
				end
				dilated_mask_M_HD = dilate(dilate(dilate(dilate(mask_M_HD_3D))))
				ring_mask_M_HD = dilate(dilate(dilate(dilate(dilate(mask_M_HD_3D))))) - dilate(dilate(dilate(dilate(mask_M_HD_3D))))
				single_ring_mask_M_HD = Bool.(ring_mask_M_HD[:, :, 3])
				s_bkg_M_HD = mean(single_arr[single_ring_mask_M_HD])
				alg_M_HD = Integrated(arr[mask_M_HD_3D])
				mass_m_hd = score(s_bkg_M_HD, S_Obj_HD, pixel_size, ρ_hd, alg_M_HD)
				
				## Medium Density
				mask_M_MD_3D = Array{Bool}(undef, size(arr))
				for z in 1:size(arr, 3)
					mask_M_MD_3D[:, :, z] = mask_M_MD
				end
				dilated_mask_M_MD = dilate(dilate(dilate(dilate(mask_M_MD_3D))))
				ring_mask_M_MD = dilate(dilate(dilate(dilate(dilate(dilate(mask_M_MD_3D)))))) - dilate(dilate(dilate(dilate(dilate(mask_M_MD_3D)))))
				single_ring_mask_M_MD = Bool.(ring_mask_M_MD[:, :, 3])
				s_bkg_M_MD = mean(single_arr[single_ring_mask_M_MD])
				alg_M_MD = Integrated(arr[mask_M_MD_3D])
				mass_m_md = score(s_bkg_M_MD, S_Obj_MD, pixel_size, ρ_md, alg_M_MD)
			
				## Low Density
				mask_M_LD_3D = Array{Bool}(undef, size(arr))
				for z in 1:size(arr, 3)
					mask_M_LD_3D[:, :, z] = mask_M_LD
				end
				dilated_mask_M_LD = dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))))
				ring_mask_M_LD = dilate(dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))))) - dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))))
				single_ring_mask_M_LD = Bool.(ring_mask_M_LD[:, :, 3])
				s_bkg_M_LD = mean(single_arr[single_ring_mask_M_LD])
				alg_M_LD = Integrated(arr[mask_M_LD_3D])
				mass_m_ld = score(s_bkg_M_LD, S_Obj_LD, pixel_size, ρ_ld, alg_M_LD)
				
			
				# Score Small Inserts
				## High Density
				mask_S_HD_3D = Array{Bool}(undef, size(arr))
				for z in 1:size(arr, 3)
					mask_S_HD_3D[:, :, z] = mask_S_HD
				end
				dilated_mask_S_HD = dilate(dilate(dilate(dilate(dilate((mask_S_HD_3D))))))
				ring_mask_S_HD = dilate(dilate(dilate(dilate(dilate(mask_S_HD_3D))))) - dilate(dilate(dilate(dilate(mask_S_HD_3D))))
				single_ring_mask_S_HD = Bool.(ring_mask_S_HD[:, :, 3])
				s_bkg_S_HD = mean(single_arr[single_ring_mask_S_HD])
				alg_S_HD = Integrated(arr[mask_S_HD_3D])
				mass_s_hd = score(s_bkg_S_HD, S_Obj_HD, pixel_size, ρ_hd, alg_S_HD)
			
				## Medium Density
				mask_S_MD_3D = Array{Bool}(undef, size(arr))
				for z in 1:size(arr, 3)
					mask_S_MD_3D[:, :, z] = mask_S_MD
				end
				dilated_mask_S_MD = dilate(dilate(dilate(dilate(dilate(mask_S_MD_3D)))))
				ring_mask_S_MD = dilate(dilate(dilate(dilate(dilate(mask_S_MD_3D))))) - dilate(dilate(dilate(dilate(mask_S_MD_3D))))
				single_ring_mask_S_MD = Bool.(ring_mask_S_MD[:, :, 3])
				s_bkg_S_MD = mean(single_arr[single_ring_mask_S_MD])
				alg_S_MD = Integrated(arr[mask_S_MD_3D])
				mass_s_md = score(s_bkg_S_MD, S_Obj_MD, pixel_size, ρ_md, alg_S_MD)
				
				## Low Density
				mask_S_LD_3D = Array{Bool}(undef, size(arr))
				for z in 1:size(arr, 3)
					mask_S_LD_3D[:, :, z] = mask_S_LD
				end
				dilated_mask_S_LD = dilate(dilate(dilate(dilate(dilate(mask_S_LD_3D)))))
				ring_mask_S_LD = dilate(dilate(dilate(dilate(dilate(mask_S_LD_3D))))) - dilate(dilate(dilate(dilate(mask_S_LD_3D))))
				single_ring_mask_S_LD = Bool.(ring_mask_S_LD[:, :, 3])
				s_bkg_S_LD = mean(single_arr[single_ring_mask_S_LD])
				alg_S_LD = Integrated(arr[mask_S_LD_3D])
				mass_s_ld = score(s_bkg_S_LD, S_Obj_LD, pixel_size, ρ_ld, alg_S_LD)
				
				# Results
				inserts = [
					"Low Density",
					"Medium Density",
					"High Density"
				]
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
					calculated_mass_small = calculated_mass_small,
					mass_bkg = mass_bkg
				)
				push!(dfs, df)
			end
		end
	end
end

# ╔═╡ 1b95b391-cad0-4aad-8885-fd0c79cd94c1
md"""
# Save Results
"""

# ╔═╡ cc20b566-f645-4ae5-8fe2-19ba5697aa6c
dfs

# ╔═╡ bb622f09-0da5-4c38-a7d3-de898af42490
if ~isdir(string(cd(pwd, "..") , "/output_repeated/", TYPE))
	mkdir(string(cd(pwd, "..") , "/output_repeated/", TYPE))
end

# ╔═╡ 30e1d9a6-868c-4d29-8752-48a3092d0759
if length(dfs) == 24
	new_df = vcat(dfs[1:24]...)
	output_path_new = string(cd(pwd, "..") , "/output_repeated/", TYPE, "/", "full.csv")
	CSV.write(output_path_new, new_df)
end

# ╔═╡ Cell order:
# ╠═70410feb-aa21-405d-8d26-9cb7c1ecdb9e
# ╠═f2845431-56a1-4fa4-a31f-a515c686d85a
# ╟─71b9dc8f-8879-43c3-a94b-c3efbeb672a6
# ╠═26110c25-c288-4e5a-a3d1-260ca0c15a84
# ╠═442943d0-9a1e-46ff-a924-91b5537a9b6a
# ╠═57656cb2-146a-43ec-abbb-1b5a6a1f3afb
# ╠═6e98101e-f052-4b19-9915-429166605ced
# ╠═40454c15-d8fa-4530-8b74-39a934afc257
# ╟─1b95b391-cad0-4aad-8885-fd0c79cd94c1
# ╠═cc20b566-f645-4ae5-8fe2-19ba5697aa6c
# ╠═bb622f09-0da5-4c38-a7d3-de898af42490
# ╠═30e1d9a6-868c-4d29-8752-48a3092d0759

### A Pluto.jl notebook ###
# v0.18.1

using Markdown
using InteractiveUtils

# ╔═╡ 70410feb-aa21-405d-8d26-9cb7c1ecdb9e
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
begin
	TYPE = "integrated_scoring"
	BASE_PATH = "/Users/daleblack/Google Drive/Datasets/Simulated/"
	SCAN_NUMBER = 1
end;

# ╔═╡ 442943d0-9a1e-46ff-a924-91b5537a9b6a
venders = ["80", "100", "120", "135"]

# ╔═╡ 83eba315-f05c-4e0b-82b9-63677be271ea
kernels = [0, 1, 2]

# ╔═╡ 40454c15-d8fa-4530-8b74-39a934afc257
begin
	dfs = []
	for kern in kernels
		for vender in venders
			SCAN_NUMBER = 1
			VENDER = vender
			root_path = string(BASE_PATH, VENDER)
			dcm_path_list = dcm_list_builder(root_path)
			pth = dcm_path_list[SCAN_NUMBER]
			scan = basename(pth)
			header, dcm_array, slice_thick_ori1 = dcm_reader(pth)
			for z in size(dcm_array, 3)
				dcm_array[:, :, z] = Noise.mult_gauss(dcm_array[:, :, z], kern)
			end
		
			# Segment Heart
			masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3)÷2)
		
			# Segment Calcium Rod
			calcium_image, slice_CCI, quality_slice, cal_rod_slice = mask_rod(masked_array, header)
		
			# Calibration Prep
			array_filtered = abs.(mapwindow(median, calcium_image[:, :, 2], (3, 3)))
			bool_arr = array_filtered .> 0
			bool_arr_erode = (((erode(erode(bool_arr)))))
			c_img = calcium_image[:, :, 1:3]
			mask_cal_3D = Array{Bool}(undef, size(c_img))
			for z in 1:size(c_img, 3)
				mask_cal_3D[:, :, z] = bool_arr_erode
			end
			cal_insert_mean = quantile!(c_img[mask_cal_3D], 0.7)
			# eroded_mask_L_HD = erode(erode(mask_L_HD_3D))
			# high_density_cal = mean(arr[eroded_mask_L_HD])
			# eroded_mask_L_MD = erode(erode(mask_L_MD_3D))
			# med_density_cal = mean(arr[eroded_mask_L_MD])
			# eroded_mask_L_LD = erode(erode(mask_L_LD_3D))
			# low_density_cal = mean(arr[eroded_mask_L_LD])
			density_array = [0, 200, 400, 800]
			density_array_cal = [0, 200]
			# intensity_array = [0, low_density_cal, med_density_cal, high_density_cal] # HU
			intensity_array_cal = [0, cal_insert_mean]
			df_cal = DataFrame(:density => density_array_cal, :intensity => intensity_array_cal)
			linearRegressor = lm(@formula(intensity ~ density), df_cal);
			linearFit = predict(linearRegressor)
			m = linearRegressor.model.pp.beta0[2]
			b = linearRegressor.model.rr.mu[1]
			density(intensity) = (intensity - b) / m
			intensity(ρ) = m*ρ + b

			mask_L_HD, mask_M_HD, mask_S_HD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_LD, mask_M_LD, mask_S_LD = mask_inserts_simulation(
				dcm_array, masked_array, header, slice_CCI, center_insert)
		
			arr = masked_array[:, :, 4:6]
			single_arr = masked_array[:, :, slice_CCI]
			pixel_size = DICOMUtils.get_pixel_size(header)
			
			# Score Large InsertS
			## High Density
			mask_L_HD_3D = Array{Bool}(undef, size(arr))
			for z in 1:size(arr, 3)
				mask_L_HD_3D[:, :, z] = mask_L_HD
			end
			dilated_mask_L_HD = dilate(dilate(mask_L_HD_3D))
			ring_mask_L_HD = dilate(dilate(dilate(dilate(mask_L_HD_3D)))) - dilate(dilate(dilate(mask_L_HD_3D)))
			single_ring_mask_L_HD = Bool.(ring_mask_L_HD[:, :, 3])
			s_bkg_L_HD = mean(single_arr[single_ring_mask_L_HD])
			S_Obj_HD = intensity(800)
			ρ_hd = 0.8 # mg/mm^3
			alg_L_HD = Integrated(arr[mask_L_HD_3D])
			mass_l_hd = score(s_bkg_L_HD, S_Obj_HD, pixel_size, ρ_hd, alg_L_HD)
		
			## Medium Density
			mask_L_MD_3D = Array{Bool}(undef, size(arr))
			for z in 1:size(arr, 3)
				mask_L_MD_3D[:, :, z] = mask_L_MD
			end
			dilated_mask_L_MD = dilate(dilate(mask_L_MD_3D))
			ring_mask_L_MD = dilate(dilate(dilate(dilate(mask_L_MD_3D)))) - dilate(dilate(dilate(mask_L_MD_3D)))
			single_ring_mask_L_MD = Bool.(ring_mask_L_MD[:, :, 3])
			s_bkg_L_MD = mean(single_arr[single_ring_mask_L_MD])
			S_Obj_MD = intensity(400)
			ρ_md = 0.4 # mg/mm^3
			alg_L_MD = Integrated(arr[mask_L_MD_3D])
			mass_l_md = score(s_bkg_L_MD, S_Obj_MD, pixel_size, ρ_md, alg_L_MD)
		
			## Low Density
			mask_L_LD_3D = Array{Bool}(undef, size(arr))
			for z in 1:size(arr, 3)
				mask_L_LD_3D[:, :, z] = mask_L_LD
			end
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
			density_array = [0, 200, 400, 800]
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
				kern = kern,
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

# ╔═╡ 1b95b391-cad0-4aad-8885-fd0c79cd94c1
md"""
# Save Results
"""

# ╔═╡ cc20b566-f645-4ae5-8fe2-19ba5697aa6c
dfs

# ╔═╡ bb622f09-0da5-4c38-a7d3-de898af42490
if ~isdir(string(cd(pwd, "..") , "/output/", TYPE))
	mkdir(string(cd(pwd, "..") , "/output/", TYPE))
end

# ╔═╡ 30e1d9a6-868c-4d29-8752-48a3092d0759
if length(dfs) == 12
	new_df = vcat(dfs[1:12]...)
	output_path_new = string(cd(pwd, "..") , "/output/", TYPE, "/", "full.csv")
	CSV.write(output_path_new, new_df)
end

# ╔═╡ Cell order:
# ╠═70410feb-aa21-405d-8d26-9cb7c1ecdb9e
# ╠═f2845431-56a1-4fa4-a31f-a515c686d85a
# ╟─71b9dc8f-8879-43c3-a94b-c3efbeb672a6
# ╠═26110c25-c288-4e5a-a3d1-260ca0c15a84
# ╠═442943d0-9a1e-46ff-a924-91b5537a9b6a
# ╠═83eba315-f05c-4e0b-82b9-63677be271ea
# ╠═40454c15-d8fa-4530-8b74-39a934afc257
# ╟─1b95b391-cad0-4aad-8885-fd0c79cd94c1
# ╠═cc20b566-f645-4ae5-8fe2-19ba5697aa6c
# ╠═bb622f09-0da5-4c38-a7d3-de898af42490
# ╠═30e1d9a6-868c-4d29-8752-48a3092d0759

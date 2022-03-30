### A Pluto.jl notebook ###
# v0.18.1

using Markdown
using InteractiveUtils

# ╔═╡ 6876b340-b075-11ec-385c-75b3de984c11
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
	using DICOM
	using DICOMUtils
	using PhantomSegmentation
	using CalciumScoring
end

# ╔═╡ bd839543-88ed-41c8-b6e1-0489a9d27714
TableOfContents()

# ╔═╡ 47510e0b-7af8-420a-8313-0581004c63e2
function create_mask(array, mask)
	@assert size(array) == size(mask)
	idxs = findall(x -> x == true, mask)
	overlayed_mask = zeros(size(array))
	for idx in idxs
		overlayed_mask[idx] = array[idx]
	end
	return overlayed_mask
end

# ╔═╡ 1bd15724-5c95-45ef-915f-27b91161f394
begin
	TYPE = "agatston"
	BASE_PATH = "/Users/daleblack/Google Drive/Datasets/Simulated/"
	SCAN_NUMBER = 1
end;

# ╔═╡ 6d44edf1-d97d-4333-9a6e-cda6b96d324f
venders = ["80", "100", "120", "135"]

# ╔═╡ 1d965bcc-1c10-486c-9218-7d3b1c089b20
kernels = [0, 3, 5]

# ╔═╡ 653f4b85-58f1-4bd7-a81e-86a8366c8395
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
			if kern != 0
				for z in size(dcm_array, 3)
					dcm_array[:, :, z] = imfilter(dcm_array[:, :, z], Kernel.gaussian(kern))
				end
			end

			# Segment Heart
			masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3)÷2)
		
			# Segment Calcium Rod
			calcium_image, slice_CCI, quality_slice, cal_rod_slice = mask_rod(masked_array, header)
	
			# Segment Calcium Inserts
			mask_L_HD, mask_M_HD, mask_S_HD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_LD, mask_M_LD, mask_S_LD = mask_inserts_simulation(
		            dcm_array, masked_array, header, slice_CCI, center_insert
			)
		
			# Mask Calibration Factor
			output = calc_output(masked_array, header, slice_CCI, 130, trues(3, 3))
		    insert_centers = calc_centers(dcm_array, output, header, center_insert, slice_CCI)
			rows, cols = Int(header[tag"Rows"]), Int(header[tag"Columns"])
			pixel_size = DICOMUtils.get_pixel_size(header)
			mass_cal_factor, angle_0_200HA, water_rod_metrics = mass_calibration(masked_array, insert_centers[:Large_LD], center_insert, 2, cols, rows, pixel_size)
			
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
			agat_l_hd, mass_l_hd = score(overlayed_mask_l_hd, pixel_size, mass_cal_factor, alg)
			
			## Medium Density
			mask_L_MD_3D = Array{Bool}(undef, size(arr))
			for z in 1:size(arr, 3)
				mask_L_MD_3D[:, :, z] = mask_L_MD
			end
			dilated_mask_L_MD = dilate(dilate(mask_L_MD_3D))
			overlayed_mask_l_md = create_mask(arr, dilated_mask_L_MD)
			agat_l_md, mass_l_md = score(overlayed_mask_l_md, pixel_size, mass_cal_factor, alg)
			
			## Low Density
			mask_L_LD_3D = Array{Bool}(undef, size(arr))
			for z in 1:size(arr, 3)
				mask_L_LD_3D[:, :, z] = mask_L_LD
			end
			dilated_mask_L_LD = dilate(dilate(mask_L_LD_3D))
			overlayed_mask_l_ld = create_mask(arr, dilated_mask_L_LD)
			agat_l_ld, mass_l_ld = score(overlayed_mask_l_ld, pixel_size, mass_cal_factor, alg)
		
		
			# Score Medium Inserts
			## High Density
			mask_M_HD_3D = Array{Bool}(undef, size(arr))
			for z in 1:size(arr, 3)
				mask_M_HD_3D[:, :, z] = mask_M_HD
			end
			dilated_mask_M_HD = dilate(dilate(dilate(dilate(mask_M_HD_3D))))
			overlayed_mask_m_hd = create_mask(arr, dilated_mask_M_HD)
			agat_m_hd, mass_m_hd = score(overlayed_mask_m_hd, pixel_size, mass_cal_factor, alg)
			
			## Medium Density
			mask_M_MD_3D = Array{Bool}(undef, size(arr))
			for z in 1:size(arr, 3)
				mask_M_MD_3D[:, :, z] = mask_M_MD
			end
			dilated_mask_M_MD = dilate(dilate(dilate(dilate(mask_M_MD_3D))))
			overlayed_mask_m_md = create_mask(arr, dilated_mask_M_MD)
			agat_m_md, mass_m_md = score(overlayed_mask_m_md, pixel_size, mass_cal_factor, alg)
			
			## Low Density
			mask_M_LD_3D = Array{Bool}(undef, size(arr))
			for z in 1:size(arr, 3)
				mask_M_LD_3D[:, :, z] = mask_M_LD
			end
			dilated_mask_M_LD = dilate(dilate(dilate(dilate(mask_M_LD_3D))))
			overlayed_mask_m_ld = create_mask(arr, dilated_mask_M_LD)
			agat_m_ld, mass_m_ld = score(overlayed_mask_m_ld, pixel_size, mass_cal_factor, alg)
			
			# Score Small Inserts
			## High Density
			mask_S_HD_3D = Array{Bool}(undef, size(arr))
			for z in 1:size(arr, 3)
				mask_S_HD_3D[:, :, z] = mask_S_HD
			end
			dilated_mask_S_HD = dilate((dilate(dilate(dilate((mask_S_HD_3D))))))
			overlayed_mask_s_hd = create_mask(arr, dilated_mask_S_HD)
			agat_s_hd, mass_s_hd = score(overlayed_mask_s_hd, pixel_size, mass_cal_factor, alg)
			
			## Medium Density
			mask_S_MD_3D = Array{Bool}(undef, size(arr))
			for z in 1:size(arr, 3)
				mask_S_MD_3D[:, :, z] = mask_S_MD
			end
			dilated_mask_S_MD = dilate((dilate(dilate(dilate(mask_S_MD_3D)))))
			overlayed_mask_s_md = create_mask(arr, dilated_mask_S_MD)
			agat_s_md, mass_s_md = score(overlayed_mask_s_md, pixel_size, mass_cal_factor, alg)
			
			## Low Density
			mask_S_LD_3D = Array{Bool}(undef, size(arr))
			for z in 1:size(arr, 3)
				mask_S_LD_3D[:, :, z] = mask_S_LD
			end
			dilated_mask_S_LD = dilate((dilate(dilate(dilate(mask_S_LD_3D)))))
			overlayed_mask_s_ld = create_mask(arr, dilated_mask_S_LD)
			agat_s_ld, mass_s_ld = score(overlayed_mask_s_ld, pixel_size, mass_cal_factor, alg)
			
			# Results
			density_array = [0, 200, 400, 800]
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
				kern = kern,
				scan = scan,
				inserts = inserts,
				calculated_agat_large = calculated_agat_large,
				calculated_agat_medium = calculated_agat_medium,
				calculated_agat_small = calculated_agat_small,
				ground_truth_mass_large = ground_truth_mass_large,
				calculated_mass_large = calculated_mass_large,
				ground_truth_mass_medium = ground_truth_mass_medium,
				calculated_mass_medium = calculated_mass_medium,
				ground_truth_mass_small = ground_truth_mass_small,
				calculated_mass_small = calculated_mass_small,
				mass_cal_factor = mass_cal_factor
			)
			push!(dfs, df)
		end
	end
end

# ╔═╡ ed2d91ef-73c7-48cd-bcd9-4bfb5d1f4ae5
md"""
# Save results
"""

# ╔═╡ a005fc54-44e3-43d1-8b4a-4767dc0e180e
dfs

# ╔═╡ ad45f9c8-b82f-4015-b912-f548ee988e05
if ~isdir(string(cd(pwd, "..") , "/output/", TYPE))
	mkdir(string(cd(pwd, "..") , "/output/", TYPE))
end

# ╔═╡ 99aee610-1b33-4b64-b752-ac8fe602bdcb
if length(dfs) == 12
	new_df = vcat(dfs[1:12]...)
	output_path_new = string(cd(pwd, "..") , "/output/", TYPE, "/", "full.csv")
	CSV.write(output_path_new, new_df)
end

# ╔═╡ Cell order:
# ╠═6876b340-b075-11ec-385c-75b3de984c11
# ╠═bd839543-88ed-41c8-b6e1-0489a9d27714
# ╟─47510e0b-7af8-420a-8313-0581004c63e2
# ╠═1bd15724-5c95-45ef-915f-27b91161f394
# ╠═6d44edf1-d97d-4333-9a6e-cda6b96d324f
# ╠═1d965bcc-1c10-486c-9218-7d3b1c089b20
# ╠═653f4b85-58f1-4bd7-a81e-86a8366c8395
# ╟─ed2d91ef-73c7-48cd-bcd9-4bfb5d1f4ae5
# ╠═a005fc54-44e3-43d1-8b4a-4767dc0e180e
# ╠═ad45f9c8-b82f-4015-b912-f548ee988e05
# ╠═99aee610-1b33-4b64-b752-ac8fe602bdcb

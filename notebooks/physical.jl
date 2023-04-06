### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ c7ef8f6c-7fb8-4bec-90e6-facbadcfca57
# ╠═╡ show_logs = false
begin
    using Pkg
    Pkg.activate(".")

    using PlutoUI, Statistics, CSV, DataFrames, GLM, CairoMakie, HypothesisTests, Colors, MLJBase, DICOM, DICOMUtils, PhantomSegmentation, CalciumScoring, ImageMorphology, ImageFiltering, Noise
    using StatsBase: quantile!, rmsd
end

# ╔═╡ 8cb5d250-92cc-4293-ae47-1ef4e115e37b
TableOfContents()

# ╔═╡ cf94ec93-180c-40d8-a2b7-7a9763b3479a
OUTPUT = "output_new"

# ╔═╡ 9f183350-beb0-4a45-a730-e1cd8140efea
SAVE_DF = "physical.csv"

# ╔═╡ 3ea0f3f3-b590-438e-96d1-c189b60f1196
VENDORS = ["Canon_Aquilion_One_Vision"]

# ╔═╡ 1298f72b-6c46-494c-91bf-34a3f7b4ca80
SCANS = [1, 3, 6, 8]
# SCANS = [1]

# ╔═╡ ad8d9fcb-ca8c-40e8-a331-469b5b68890c
function create_mask(array, mask)
    @assert size(array) == size(mask)
    idxs = findall(x -> x == true, mask)
    overlayed_mask = zeros(size(array))
    for idx in idxs
        overlayed_mask[idx] = array[idx]
    end
    return overlayed_mask
end

# ╔═╡ 3c2c7d5e-0a96-4008-8bdc-a5ea28114617
function dilate_mask_large(mask)
    return dilate(mask)
end

# ╔═╡ 528b8447-4f68-4128-978b-01d2df81045b
function ring_mask_large(dilated_mask)
    return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask)
end

# ╔═╡ 27794b80-3f8f-4eac-a7fc-9bfc25676f5d
function dilate_mask_medium(mask)
    return dilate(dilate(dilate(mask)))
end

# ╔═╡ 261851fd-aa68-475d-9552-91e94c242164
function ring_mask_medium(dilated_mask)
    return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask)
end

# ╔═╡ 0dc125bf-1447-4628-82dd-621afed62469
function dilate_mask_small(mask)
    return dilate(dilate(dilate(mask)))
end

# ╔═╡ bb619fa9-37a9-492a-9088-ba211833266d
function ring_mask_small(dilated_mask)
    return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask)
end

# ╔═╡ eba008b8-c288-4af5-b78d-48cf3d64e959
function dilate_mask_large_bkg(mask)
    return dilate(dilate(mask))
end

# ╔═╡ fde4f7f4-3abd-498f-92b9-867292c6d726
function dilate_mask_medium_bkg(mask)
    return dilate(mask)
end

# ╔═╡ 83e6b776-2a80-4e3b-9bba-045bd0ebab59
function dilate_mask_small_bkg(mask)
    return (mask)
end

# ╔═╡ 5d49cbd7-9d31-4e6a-98fd-4ef7cb87370e
joinpath(dirname(dirname(dirname(dirname(pwd())))), "Datasets", "CAC Data")

# ╔═╡ 28539bfd-9837-4ff8-90ff-b93633c5c282
"/Users/daleblack/Library/CloudStorage/GoogleDrive-djblack@uci.edu/My Drive/Datasets/CAC Data/Canon_Aquilion_One_Vision"

# ╔═╡ ca28ecc0-0324-43e7-8126-ccde7ef5f05b
begin
    dfs_vf = []
	dfs_s = []
	dfs_i = []
	dfs_a = []
    high_dens = []
    med_dens = []
    low_dens = []
    high_dens100 = []
    med_dens50 = []
    low_dens25 = []
    for VENDOR in VENDORS
		for SCAN in SCANS
			#---------------- Reusable Pieces ----------------#
			BASE_PATH = joinpath(dirname(dirname(dirname(dirname(pwd())))), "Datasets", "CAC Data")
			root_path = joinpath(BASE_PATH, VENDOR)
			dcm_path_list = dcm_list_builder(root_path)
			pth = dcm_path_list[SCAN]
			scan = basename(pth)
			header, dcm_array, slice_thick_ori1 = dcm_reader(pth)
			kV = header[tag"KVP"]
	
			# Segment Heart
			masked_array, center_insert, mask = mask_heart(
				header, dcm_array, size(dcm_array, 3) ÷ 2
			)
	
			# Segment Calcium Rod
			thresh = 130
	
			@info VENDOR, SCAN
	
			calcium_image, slice_CCI, quality_slice, cal_rod_slice = mask_rod(
				masked_array, header; calcium_threshold=thresh
			)

			# Segment Calcium Inserts
            mask_L_HD, mask_M_HD, mask_S_HD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_LD, mask_M_LD, mask_S_LD = mask_inserts(
                dcm_array, masked_array, header, slice_CCI, center_insert;
                calcium_threshold=thresh
            )
	
			# Prepare Mask Arrays
			arr = masked_array[:, :, slice_CCI-2:slice_CCI+2]
			pixel_size = DICOMUtils.get_pixel_size(header)
			voxel_size = pixel_size[1] * pixel_size[2] * pixel_size[3]

			## Background
			background_mask1 = zeros(size(arr)...)
			background_mask2 = zeros(size(arr)...)
			background_mask3 = zeros(size(arr)...)

			offs = 10
			rnge = 10
			background_mask1[
				(center_insert[1]-rnge+offs):(center_insert[1]+rnge+offs),
				(center_insert[2]-rnge):(center_insert[2]+rnge),
				2,
			] .= 1
			background_mask2[
				(center_insert[1]-rnge-offs):(center_insert[1]+rnge-offs),
				(center_insert[2]-rnge):(center_insert[2]+rnge),
				2,
			] .= 1
			background_mask3[
				(center_insert[1]-rnge):(center_insert[1]+rnge),
				(center_insert[2]-rnge+offs):(center_insert[2]+rnge+offs),
				2,
			] .= 1

			dilated_mask_L_bkg = dilate_mask_small_bkg(Bool.(background_mask1))
			ring_mask_L_bkg = ring_mask_small(dilated_mask_L_bkg)

			dilated_mask_M_bkg = dilate_mask_small_bkg(Bool.(background_mask2))
			ring_mask_M_bkg = ring_mask_small(dilated_mask_M_bkg)

			dilated_mask_S_bkg = dilate_mask_small_bkg(Bool.(background_mask3))
			ring_mask_S_bkg = ring_mask_small(dilated_mask_S_bkg)
			
			## Large
			### High Density
			mask_L_HD_3D = Array{Bool}(undef, size(arr))
			for z in axes(arr, 3)
				mask_L_HD_3D[:, :, z] = mask_L_HD
			end
			dilated_mask_L_HD = dilate_mask_large(mask_L_HD_3D)
			ring_mask_L_HD = ring_mask_large(dilated_mask_L_HD)
	
			### Medium Density
			mask_L_MD_3D = Array{Bool}(undef, size(arr))
			for z in axes(arr, 3)
				mask_L_MD_3D[:, :, z] = mask_L_MD
			end
			dilated_mask_L_MD = dilate_mask_large(mask_L_MD_3D)
			ring_mask_L_MD = ring_mask_large(dilated_mask_L_MD)
	
			### Low Density
			mask_L_LD_3D = Array{Bool}(undef, size(arr))
			for z in axes(arr, 3)
				mask_L_LD_3D[:, :, z] = mask_L_LD
			end
			dilated_mask_L_LD = dilate_mask_large(mask_L_LD_3D)
			ring_mask_L_LD = ring_mask_large(dilated_mask_L_LD)
	
	
			## Medium 
			### High Density
			mask_M_HD_3D = Array{Bool}(undef, size(arr))
			for z in axes(arr, 3)
				mask_M_HD_3D[:, :, z] = mask_M_HD
			end
			dilated_mask_M_HD = dilate_mask_medium(mask_M_HD_3D)
			ring_mask_M_HD = ring_mask_medium(dilated_mask_M_HD)
	
			### Medium Density
			mask_M_MD_3D = Array{Bool}(undef, size(arr))
			for z in axes(arr, 3)
				mask_M_MD_3D[:, :, z] = mask_M_MD
			end
			dilated_mask_M_MD = dilate_mask_medium(mask_M_MD_3D)
			ring_mask_M_MD = ring_mask_medium(dilated_mask_M_MD)
	
			### Low Density
			mask_M_LD_3D = Array{Bool}(undef, size(arr))
			for z in axes(arr, 3)
				mask_M_LD_3D[:, :, z] = mask_M_LD
			end
			dilated_mask_M_LD = dilate_mask_medium(mask_M_LD_3D)
			ring_mask_M_LD = ring_mask_medium(dilated_mask_M_LD)
	
			## Small
			### High Density
			mask_S_HD_3D = Array{Bool}(undef, size(arr))
			for z in axes(arr, 3)
				mask_S_HD_3D[:, :, z] = mask_S_HD
			end
			dilated_mask_S_HD = dilate_mask_small(mask_S_HD_3D)
			ring_mask_S_HD = ring_mask_small(dilated_mask_S_HD)
	
			### Medium Density
			mask_S_MD_3D = Array{Bool}(undef, size(arr))
			for z in axes(arr, 3)
				mask_S_MD_3D[:, :, z] = mask_S_MD
			end
			dilated_mask_S_MD = dilate_mask_small(mask_S_MD_3D)
			ring_mask_S_MD = ring_mask_small(dilated_mask_S_MD)
	
			### Low Density
			mask_S_LD_3D = Array{Bool}(undef, size(arr))
			for z in axes(arr, 3)
				mask_S_LD_3D[:, :, z] = mask_S_LD
			end
			dilated_mask_S_LD = dilate_mask_small(mask_S_LD_3D)
			ring_mask_S_LD = ring_mask_small(dilated_mask_S_LD)
	
			#---------------- Volume Fraction ----------------#
	
			## Calibration Prep
			density_array = [0.200, 0.400, 0.800]
			array_filtered = abs.(mapwindow(median, calcium_image[:, :, cal_rod_slice], (3, 3)))
			bool_arr = array_filtered .> 0
			bool_arr_erode = (erode(erode(erode(erode(bool_arr)))))
			c_img = calcium_image[:, :, cal_rod_slice-1:cal_rod_slice+2]
			mask_cal_3D = zeros(size(c_img))
			for z in axes(c_img, 3)
				mask_cal_3D[:, :, z] = Bool.(erode(bool_arr_erode))
			end
		
			hu_calcium = mean(c_img[Bool.(mask_cal_3D)])
			hu_heart_tissue = mean(arr[dilated_mask_L_bkg])
			ρ_calcium = 0.2

			# Background
			hu_heart_tissue_large_bkg = mean(arr[ring_mask_L_bkg])
			mass_large_bkg = score(arr[dilated_mask_L_bkg], hu_calcium, hu_heart_tissue_large_bkg, voxel_size, ρ_calcium, VolumeFraction())
		
			hu_heart_tissue_medium_bkg = mean(arr[ring_mask_M_bkg])
			mass_medium_bkg = score(arr[dilated_mask_M_bkg], hu_calcium, hu_heart_tissue_medium_bkg, voxel_size, ρ_calcium, VolumeFraction())
			
			hu_heart_tissue_small_bkg = mean(arr[ring_mask_S_bkg])
			mass_small_bkg = score(arr[dilated_mask_S_bkg], hu_calcium, hu_heart_tissue_small_bkg, voxel_size, ρ_calcium, VolumeFraction())
	
			mass_bkg = [mass_large_bkg, mass_medium_bkg, mass_small_bkg]
	
			# Score Large Inserts
			## High Density
			hu_heart_tissue_large_hd = mean(arr[ring_mask_L_HD])
			mass_large_hd = score(arr[dilated_mask_L_HD], hu_calcium, hu_heart_tissue_large_hd, voxel_size, ρ_calcium, VolumeFraction())
	
			## Medium Density
			hu_heart_tissue_large_md = mean(arr[ring_mask_L_MD])
			mass_large_md = score(arr[dilated_mask_L_MD], hu_calcium, hu_heart_tissue_large_md, voxel_size, ρ_calcium, VolumeFraction())
	
			## Low Density
			hu_heart_tissue_large_ld = mean(arr[ring_mask_L_LD])
			mass_large_ld = score(arr[dilated_mask_L_LD], hu_calcium, hu_heart_tissue_large_ld, voxel_size, ρ_calcium, VolumeFraction())
	
			# Score Medium Inserts
			## High Density
			hu_heart_tissue_medium_hd = mean(arr[ring_mask_M_HD])
			mass_medium_hd = score(arr[dilated_mask_M_HD], hu_calcium, hu_heart_tissue_medium_hd, voxel_size, ρ_calcium, VolumeFraction())
	
			## Medium Density
			hu_heart_tissue_medium_md = mean(arr[ring_mask_M_MD])
			mass_medium_md = score(arr[dilated_mask_M_MD], hu_calcium, hu_heart_tissue_medium_md, voxel_size, ρ_calcium, VolumeFraction())
	
			## Low Density
			hu_heart_tissue_medium_ld = mean(arr[ring_mask_M_LD])
			mass_medium_ld = score(arr[dilated_mask_M_LD], hu_calcium, hu_heart_tissue_medium_ld, voxel_size, ρ_calcium, VolumeFraction())
	
			# Score Small Inserts
			## High Density
			hu_heart_tissue_small_hd = mean(arr[ring_mask_S_HD])
			mass_small_hd = score(arr[dilated_mask_S_HD], hu_calcium, hu_heart_tissue_large_hd, voxel_size, ρ_calcium, VolumeFraction())
	
			## Medium Density
			hu_heart_tissue_small_md = mean(arr[ring_mask_S_MD])
			mass_small_md = score(arr[dilated_mask_S_MD], hu_calcium, hu_heart_tissue_large_md, voxel_size, ρ_calcium, VolumeFraction())
	
			## Low Density
			hu_heart_tissue_small_ld = mean(arr[ring_mask_S_LD])
			mass_small_ld = score(arr[dilated_mask_S_LD], hu_calcium, hu_heart_tissue_large_ld, voxel_size, ρ_calcium, VolumeFraction())
	
			# Results
	
			inserts = ["Low Density", "Medium Density", "High Density"]
			volume_gt = [7.065, 63.585, 176.625]
			ground_truth_mass_large = [
				volume_gt[3] * density_array[1],
				volume_gt[3] * density_array[2],
				volume_gt[3] * density_array[3],
			] # mg
	
			calculated_mass_large = [mass_large_ld, mass_large_md, mass_large_hd]
			ground_truth_mass_medium = [
				volume_gt[2] * density_array[1],
				volume_gt[2] * density_array[2],
				volume_gt[2] * density_array[3],
			] # mg
			calculated_mass_medium = [mass_medium_ld, mass_medium_md, mass_medium_hd]
			ground_truth_mass_small = [
				volume_gt[1] * density_array[1],
				volume_gt[1] * density_array[2],
				volume_gt[1] * density_array[3],
			] # mg
			calculated_mass_small = [mass_small_ld, mass_small_md, mass_small_hd]
	
			df = DataFrame(;
				VENDOR=VENDOR,
				scan=scan,
				inserts=inserts,
				ground_truth_mass_large=ground_truth_mass_large,
				calculated_mass_large=calculated_mass_large,
				ground_truth_mass_medium=ground_truth_mass_medium,
				calculated_mass_medium=calculated_mass_medium,
				ground_truth_mass_small=ground_truth_mass_small,
				calculated_mass_small=calculated_mass_small,
				mass_bkg = mass_bkg
			)
			push!(dfs_vf, df)
	
			#---------------- Integrated ----------------#
	
			intensity_array = [hu_heart_tissue, hu_calcium]
			df_cal = DataFrame(:density => [0, 0.200], :intensity => intensity_array)
	
			linearRegressor = lm(@formula(intensity ~ density), df_cal)
			linearFit = GLM.predict(linearRegressor)
			m = linearRegressor.model.pp.beta0[2]
			b = linearRegressor.model.rr.mu[1]
			density(intensity) = (intensity - b) / m
			intensity(ρ) = m * ρ + b
			
			ρ = 0.1 # mg/mm^3
			S_Obj = intensity(ρ)

			# Score background
			S_bkg_large = mean(arr[ring_mask_L_bkg])
			alg_bkg_large = Integrated(arr[dilated_mask_L_bkg])
			mass_large_bkg = score(S_bkg_large, S_Obj, pixel_size, ρ, alg_bkg_large)

			S_bkg_medium = mean(arr[ring_mask_M_bkg])
			alg_bkg_medium = Integrated(arr[dilated_mask_M_bkg])
			mass_medium_bkg = score(S_bkg_medium, S_Obj, pixel_size, ρ, alg_bkg_medium)

			S_bkg_small = mean(arr[ring_mask_S_bkg])
			alg_bkg_small = Integrated(arr[dilated_mask_S_bkg])
			mass_small_bkg = score(S_bkg_small, S_Obj, pixel_size, ρ, alg_bkg_small)

			mass_bkg = [mass_large_bkg, mass_medium_bkg, mass_small_bkg]
	
			# Score Large Insert
			## High Density
			s_bkg_L_HD = mean(arr[ring_mask_L_HD])
			alg_L_HD = Integrated(arr[dilated_mask_L_HD])
			mass_l_hd = score(s_bkg_L_HD, S_Obj, pixel_size, ρ, alg_L_HD)
	
			## Medium Density
			s_bkg_L_MD = mean(arr[ring_mask_L_MD])
			alg_L_MD = Integrated(arr[dilated_mask_L_MD])
			mass_l_md = score(s_bkg_L_MD, S_Obj, pixel_size, ρ, alg_L_MD)
	
			## Low Density
			s_bkg_L_LD = mean(arr[ring_mask_L_LD])
			alg_L_LD = Integrated(arr[dilated_mask_L_LD])
			mass_l_ld = score(s_bkg_L_LD, S_Obj, pixel_size, ρ, alg_L_LD)
	
			# Score Medium Inserts
			## High Density
			s_bkg_M_HD = mean(arr[ring_mask_M_HD])
			alg_M_HD = Integrated(arr[dilated_mask_M_HD])
			mass_m_hd = score(s_bkg_M_HD, S_Obj, pixel_size, ρ, alg_M_HD)
	
			## Medium Density
			s_bkg_M_MD = mean(arr[ring_mask_M_MD])
			alg_M_MD = Integrated(arr[dilated_mask_M_MD])
			mass_m_md = score(s_bkg_M_MD, S_Obj, pixel_size, ρ, alg_M_MD)
	
			## Low Density
			s_bkg_M_LD = mean(arr[ring_mask_M_LD])
			alg_M_LD = Integrated(arr[dilated_mask_M_LD])
			mass_m_ld = score(s_bkg_M_LD, S_Obj, pixel_size, ρ, alg_M_LD)
	
			# Score Small Inserts
			## High Density
			s_bkg_S_HD = mean(arr[ring_mask_S_HD])
			alg_S_HD = Integrated(arr[dilated_mask_S_HD])
			mass_s_hd = score(s_bkg_S_HD, S_Obj, pixel_size, ρ, alg_S_HD)
	
			## Medium Density
			s_bkg_S_MD = mean(arr[ring_mask_S_MD])
			alg_S_MD = Integrated(arr[dilated_mask_S_MD])
			mass_s_md = score(s_bkg_S_MD, S_Obj, pixel_size, ρ, alg_S_MD)
	
			## Low Density
			s_bkg_S_LD = mean(arr[ring_mask_S_LD])
			alg_S_LD = Integrated(arr[dilated_mask_S_LD])
			mass_s_ld = score(s_bkg_S_LD, S_Obj, pixel_size, ρ, alg_S_LD)
	
			# Results
	
			calculated_mass_large = [mass_l_ld, mass_l_md, mass_l_hd]
			calculated_mass_medium = [mass_m_ld, mass_m_md, mass_m_hd]
			calculated_mass_small = [mass_s_ld, mass_s_md, mass_s_hd]
	
			df = DataFrame(;
				VENDOR=VENDOR,
				scan=scan,
				inserts=inserts,
				ground_truth_mass_large=ground_truth_mass_large,
				calculated_mass_large=calculated_mass_large,
				ground_truth_mass_medium=ground_truth_mass_medium,
				calculated_mass_medium=calculated_mass_medium,
				ground_truth_mass_small=ground_truth_mass_small,
				calculated_mass_small=calculated_mass_small,
				mass_bkg=mass_bkg
			)
			push!(dfs_i, df)
	
			#---------------- SWCS ----------------#
			μ, σ = mean(c_img[Bool.(mask_cal_3D)]) / 2, std(c_img[Bool.(mask_cal_3D)])
	
			# Mask Calibration Factor
			output = calc_output(masked_array, header, slice_CCI, thresh, trues(3, 3))
            insert_centers = calc_centers(dcm_array, output, header, center_insert, slice_CCI)
            rows, cols = Int(header[tag"Rows"]), Int(header[tag"Columns"])
            pixel_size = DICOMUtils.get_pixel_size(header)
            mass_cal_factor, angle_0_200HA, water_rod_metrics = mass_calibration(masked_array, insert_centers[:Large_LD], center_insert, cal_rod_slice, rows, cols, pixel_size)
	
			# Background
			alg2 = SpatiallyWeighted()

			swcs_bkg_large = score(dilated_mask_L_bkg, μ, σ, alg2)
			swcs_bkg_medium = score(dilated_mask_M_bkg, μ, σ, alg2)
			swcs_bkg_small = score(dilated_mask_S_bkg, μ, σ, alg2)

			swcs_bkg = [swcs_bkg_large, swcs_bkg_medium, swcs_bkg_small]
	
			# Score Large Inserts
			## High Density
			overlayed_mask_l_hd = create_mask(arr, dilated_mask_L_HD)
			swcs_l_hd = score(overlayed_mask_l_hd, μ, σ, alg2)
	
			## Medium Density
			overlayed_mask_l_md = create_mask(arr, dilated_mask_L_MD)
			swcs_l_md = score(overlayed_mask_l_md, μ, σ, alg2)
	
			## Low Density
			overlayed_mask_l_ld = create_mask(arr, dilated_mask_L_LD)
			swcs_l_ld = score(overlayed_mask_l_ld, μ, σ, alg2)
	
			# Score Medium Inserts
			## High Density
			overlayed_mask_m_hd = create_mask(arr, dilated_mask_M_HD)
			swcs_m_hd = score(overlayed_mask_m_hd, μ, σ, alg2)
	
			## Medium Density
			overlayed_mask_m_md = create_mask(arr, dilated_mask_M_MD)
			swcs_m_md = score(overlayed_mask_m_md, μ, σ, alg2)
	
			## Low Density
			overlayed_mask_m_ld = create_mask(arr, dilated_mask_M_LD)
			swcs_m_ld = score(overlayed_mask_m_ld, μ, σ, alg2)
	
			# Score Small Inserts
			## High Density
			overlayed_mask_s_hd = create_mask(arr, dilated_mask_S_HD)
			swcs_s_hd = score(overlayed_mask_s_hd, μ, σ, alg2)
	
			## Medium Density
			overlayed_mask_s_md = create_mask(arr, dilated_mask_S_MD)
			swcs_s_md = score(overlayed_mask_s_md, μ, σ, alg2)
	
			## Low Density
			overlayed_mask_s_ld = create_mask(arr, dilated_mask_S_LD)
			swcs_s_ld = score(overlayed_mask_s_ld, μ, σ, alg2)
	
			# Results
			calculated_swcs_large = [swcs_l_ld, swcs_l_md, swcs_l_hd]
			calculated_swcs_medium = [swcs_m_ld, swcs_m_md, swcs_m_hd]
			calculated_swcs_small = [swcs_s_ld, swcs_s_md, swcs_s_hd]
	
			df = DataFrame(;
				VENDOR=VENDOR,
				scan=scan,
				inserts=inserts,
				ground_truth_mass_large=ground_truth_mass_large,
				ground_truth_mass_medium=ground_truth_mass_medium,
				ground_truth_mass_small=ground_truth_mass_small,
				calculated_swcs_large=calculated_swcs_large,
				calculated_swcs_medium=calculated_swcs_medium,
				calculated_swcs_small=calculated_swcs_small,
				swcs_bkg=swcs_bkg
			)
			push!(dfs_s, df)
	
			
	
			#---------------- Agatston ----------------#
			# Mask Calibration Factor
			output = calc_output(masked_array, header, slice_CCI, thresh, trues(3, 3))
            insert_centers = calc_centers(dcm_array, output, header, center_insert, slice_CCI)
            rows, cols = Int(header[tag"Rows"]), Int(header[tag"Columns"])
            pixel_size = DICOMUtils.get_pixel_size(header)
            mass_cal_factor, angle_0_200HA, water_rod_metrics = mass_calibration(masked_array, insert_centers[:Large_LD], center_insert, cal_rod_slice, rows, cols, pixel_size)
	
			# Background

			alg = Agatston()
			overlayed_bkg_mask_L = create_mask(arr, dilated_mask_L_bkg)
			overlayed_bkg_mask_M = create_mask(arr, dilated_mask_M_bkg)
			overlayed_bkg_mask_S = create_mask(arr, dilated_mask_S_bkg)

			agat_bkg, vol_bkg_large, mass_bkg_large = score(
				overlayed_bkg_mask_L,
				pixel_size,
				mass_cal_factor,
				alg;
				kV=kV
			)
			agat_bkg, vol_bkg_medium, mass_bkg_medium = score(
				overlayed_bkg_mask_M,
				pixel_size,
				mass_cal_factor,
				alg;
				kV=kV
			)
			agat_bkg, vol_bkg_small, mass_bkg_small = score(
				overlayed_bkg_mask_S,
				pixel_size,
				mass_cal_factor,
				alg;
				kV=kV
			)

			mass_bkg = [mass_bkg_large, mass_bkg_medium, mass_bkg_small]

			# Score Large Inserts
			## High Density
			alg = Agatston()
			overlayed_mask_l_hd = create_mask(arr, dilated_mask_L_HD)
			agat_l_hd, vol_l_hd, mass_l_hd = score(
				overlayed_mask_l_hd,
				pixel_size,
				mass_cal_factor,
				alg;
				kV=kV
			)

			## Medium Density
			overlayed_mask_l_md = create_mask(arr, dilated_mask_L_MD)
			agat_l_md, vol_l_md, mass_l_md = score(
				overlayed_mask_l_md,
				pixel_size,
				mass_cal_factor,
				alg;
				kV=kV
			)

			## Low Density
			overlayed_mask_l_ld = create_mask(arr, dilated_mask_L_LD)
			agat_l_ld, vol_l_ld, mass_l_ld = score(
				overlayed_mask_l_ld,
				pixel_size,
				mass_cal_factor,
				alg;
				kV=kV
			)

			# Score Medium Inserts
			## High Density
			overlayed_mask_m_hd = create_mask(arr, dilated_mask_M_HD)
			agat_m_hd, vol_m_hd, mass_m_hd = score(
				overlayed_mask_m_hd,
				pixel_size,
				mass_cal_factor,
				alg;
				kV=kV
			)

			## Medium Density
			overlayed_mask_m_md = create_mask(arr, dilated_mask_M_MD)
			agat_m_md, vol_m_md, mass_m_md = score(
				overlayed_mask_m_md,
				pixel_size,
				mass_cal_factor,
				alg;
				kV=kV
			)

			## Low Density
			overlayed_mask_m_ld = create_mask(arr, dilated_mask_M_LD)
			agat_m_ld, vol_m_ld, mass_m_ld = score(
				overlayed_mask_m_ld,
				pixel_size,
				mass_cal_factor,
				alg;
				kV=kV
			)

			# Score Small Inserts
			## High Density
			overlayed_mask_s_hd = create_mask(arr, dilated_mask_S_HD)
			agat_s_hd, vol_s_hd, mass_s_hd = score(
				overlayed_mask_s_hd,
				pixel_size,
				mass_cal_factor,
				alg;
				kV=kV
			)

			## Medium Density
			overlayed_mask_s_md = create_mask(arr, dilated_mask_S_MD)
			agat_s_md, vol_s_md, mass_s_md = score(
				overlayed_mask_s_md,
				pixel_size,
				mass_cal_factor,
				alg;
				kV=kV
			)

			## Low Density
			overlayed_mask_s_ld = create_mask(arr, dilated_mask_S_LD)
			agat_s_ld, vol_s_ld, mass_s_ld = score(
				overlayed_mask_s_ld,
				pixel_size,
				mass_cal_factor,
				alg;
				kV=kV
			)

			# Results

			## Agatston
			calculated_agat_large = [agat_l_ld, agat_l_md, agat_l_hd]
			calculated_agat_medium = [agat_m_ld, agat_m_md, agat_m_hd]
			calculated_agat_small = [agat_s_ld, agat_s_md, agat_s_hd]

			## Mass
			calculated_mass_large = [mass_l_ld, mass_l_md, mass_l_hd]
			calculated_mass_medium = [mass_m_ld, mass_m_md, mass_m_hd]
			calculated_mass_small = [mass_s_ld, mass_s_md, mass_s_hd]

			df = DataFrame(;
				VENDOR=VENDOR,
				scan=scan,
				inserts=inserts,
				ground_truth_mass_large=ground_truth_mass_large,
				calculated_mass_large=calculated_mass_large,
				ground_truth_mass_medium=ground_truth_mass_medium,
				calculated_agat_large=calculated_agat_large,
				calculated_agat_medium=calculated_agat_medium,
				calculated_agat_small=calculated_agat_small,
				calculated_mass_medium=calculated_mass_medium,
				ground_truth_mass_small=ground_truth_mass_small,
				calculated_mass_small=calculated_mass_small,
				mass_bkg=mass_bkg
			)
			push!(dfs_a, df)
		end
    end
end

# ╔═╡ d5f169dd-54a5-4829-a451-9be4c9ac2597
md"""
# Save Results
"""

# ╔═╡ 5aea6935-1a10-4040-8108-1dc8904462ea
dfs_vf

# ╔═╡ 91301347-5513-4759-8c13-510834c114b4
dfs_i

# ╔═╡ 72259444-e778-46f2-9a45-5b79d5d3c992
dfs_s

# ╔═╡ f501934b-afb3-48ac-8f37-1ab64d3e95bf
dfs_a

# ╔═╡ 29eff88e-0520-483c-b598-0d995934a460
begin
	vf_path = joinpath(cd(pwd, ".."), OUTPUT, "volume_fraction")
	if ~isdir(vf_path)
	    mkpath(vf_path)
	end

	s_path = joinpath(cd(pwd, ".."), OUTPUT, "swcs")
	if ~isdir(s_path)
	    mkpath(s_path)
	end

	i_path = joinpath(cd(pwd, ".."), OUTPUT, "integrated")
	if ~isdir(i_path)
	    mkpath(i_path)
	end

	a_path = joinpath(cd(pwd, ".."), OUTPUT, "agatston")
	if ~isdir(a_path)
	    mkpath(a_path)
	end
end

# ╔═╡ 0fed1ff3-2285-4d97-b4ab-f6a647383493
begin
    new_df = vcat(dfs_vf[1:length(dfs_vf)]...)
    output_path_new = joinpath(vf_path, SAVE_DF)
    CSV.write(output_path_new, new_df)

	new_df = vcat(dfs_s[1:length(dfs_s)]...)
    output_path_new = joinpath(s_path, SAVE_DF)
    CSV.write(output_path_new, new_df)

	new_df = vcat(dfs_i[1:length(dfs_i)]...)
    output_path_new = joinpath(i_path, SAVE_DF)
    CSV.write(output_path_new, new_df)

	new_df = vcat(dfs_a[1:length(dfs_a)]...)
    output_path_new = joinpath(a_path, SAVE_DF)
    CSV.write(output_path_new, new_df)
end

# ╔═╡ Cell order:
# ╠═c7ef8f6c-7fb8-4bec-90e6-facbadcfca57
# ╠═8cb5d250-92cc-4293-ae47-1ef4e115e37b
# ╠═cf94ec93-180c-40d8-a2b7-7a9763b3479a
# ╠═9f183350-beb0-4a45-a730-e1cd8140efea
# ╠═3ea0f3f3-b590-438e-96d1-c189b60f1196
# ╠═1298f72b-6c46-494c-91bf-34a3f7b4ca80
# ╟─ad8d9fcb-ca8c-40e8-a331-469b5b68890c
# ╟─3c2c7d5e-0a96-4008-8bdc-a5ea28114617
# ╟─528b8447-4f68-4128-978b-01d2df81045b
# ╟─27794b80-3f8f-4eac-a7fc-9bfc25676f5d
# ╟─261851fd-aa68-475d-9552-91e94c242164
# ╟─0dc125bf-1447-4628-82dd-621afed62469
# ╟─bb619fa9-37a9-492a-9088-ba211833266d
# ╟─eba008b8-c288-4af5-b78d-48cf3d64e959
# ╟─fde4f7f4-3abd-498f-92b9-867292c6d726
# ╟─83e6b776-2a80-4e3b-9bba-045bd0ebab59
# ╠═5d49cbd7-9d31-4e6a-98fd-4ef7cb87370e
# ╠═28539bfd-9837-4ff8-90ff-b93633c5c282
# ╠═ca28ecc0-0324-43e7-8126-ccde7ef5f05b
# ╟─d5f169dd-54a5-4829-a451-9be4c9ac2597
# ╠═5aea6935-1a10-4040-8108-1dc8904462ea
# ╠═91301347-5513-4759-8c13-510834c114b4
# ╠═72259444-e778-46f2-9a45-5b79d5d3c992
# ╠═f501934b-afb3-48ac-8f37-1ab64d3e95bf
# ╠═29eff88e-0520-483c-b598-0d995934a460
# ╠═0fed1ff3-2285-4d97-b4ab-f6a647383493

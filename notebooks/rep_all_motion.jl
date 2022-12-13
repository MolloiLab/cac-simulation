### A Pluto.jl notebook ###
# v0.19.16

using Markdown
using InteractiveUtils

# ╔═╡ c4e1c65a-f078-4633-b785-78f907bb9dcc
# ╠═╡ show_logs = false
begin
    using Pkg
    Pkg.activate(".")

    using PlutoUI, Statistics, CSV, DataFrames, GLM, CairoMakie, HypothesisTests, Colors, MLJBase, DICOM, DICOMUtils, PhantomSegmentation, CalciumScoring, ImageMorphology, ImageFiltering, Noise
    using StatsBase: quantile!, rmsd
end

# ╔═╡ 61474b4a-9139-4021-b41e-15e885f9466a
TableOfContents()

# ╔═╡ 0028f871-502e-43dd-a43a-624c032f5c8e
IMAGES = "images_reproducibility1"

# ╔═╡ 0f466539-11a6-4d6d-b579-53ffb9a25a39
OUTPUT = "output_repeated"

# ╔═╡ c2e4696a-e2e2-4f02-ab06-9ebe62732d0b
SAVE_DF = "motion.csv"

# ╔═╡ bcb96a7f-7553-426e-a1f5-317bf55d1856
VENDERS = ["80", "100", "120", "135"]

# ╔═╡ 567e50fb-3d49-4136-aa32-a0ff38a61c88
SIZES = ["small", "medium", "large"]

# ╔═╡ fdb032f0-3949-423f-9bd9-08b53d65f243
DENSITIES = ["low", "normal"]

# ╔═╡ e73c0eea-3ac4-42ad-a118-040729af29cc
function create_mask(array, mask)
    @assert size(array) == size(mask)
    idxs = findall(x -> x == true, mask)
    overlayed_mask = zeros(size(array))
    for idx in idxs
        overlayed_mask[idx] = array[idx]
    end
    return overlayed_mask
end

# ╔═╡ cb23368e-79f8-4db0-b272-430956751ff6
function dilate_mask_large(mask)
    return dilate(mask)
end

# ╔═╡ 10a9aa52-a8ed-4e76-ba90-43d8c77a33b7
function ring_mask_large(dilated_mask)
    return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask)
end

# ╔═╡ 715a8838-a7f0-437d-aa5f-68744e7ca668
function dilate_mask_medium(mask)
    return (mask)
end

# ╔═╡ facd5777-30ae-4869-9ccc-5e4128cc2de7
function ring_mask_medium(dilated_mask)
    return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask)
end

# ╔═╡ e7626765-d2ca-4129-8b47-f5925a2c0eb1
function dilate_mask_small(mask)
    return (mask)
end

# ╔═╡ e6387a3a-2377-4761-a6c4-77e375afb372
function ring_mask_small(dilated_mask)
    return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask)
end

# ╔═╡ b26eab05-14e1-40ad-980b-eb8f5206542b
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
    for VENDER in VENDERS
        for SIZE in SIZES
            for DENSITY in DENSITIES
				#---------------- Volume Fraction ----------------#
                SCAN_NUMBER = 1
                BASE_PATH = joinpath("/Users/daleblack/Google Drive/dev/MolloiLab/cac-simulation", IMAGES, SIZE, DENSITY)
                root_path = joinpath(BASE_PATH, string(VENDER * "-motion"))
                dcm_path_list = dcm_list_builder(root_path)
                pth = dcm_path_list[SCAN_NUMBER]
                scan = basename(pth)
                header, dcm_array, slice_thick_ori1 = dcm_reader(pth)
				kV = header[tag"KVP"]

                # Segment Heart
                masked_array, center_insert, mask = mask_heart(
                    header, dcm_array, size(dcm_array, 3) ÷ 2
                )

                # Segment Calcium Rod
                local thresh
                if DENSITY == "low"
                    thresh = 55
                elseif DENSITY == "normal"
                    thresh = 130
                end

                @info VENDER, SIZE, DENSITY

                calcium_image, slice_CCI, quality_slice, cal_rod_slice = mask_rod(
                    masked_array, header; calcium_threshold=thresh
                )

                root_new = string(
                    "/Users/daleblack/Google Drive/dev/MolloiLab/cac-simulation/julia_arrays/",
                    SIZE,
                    "/",
                )
                mask_L_HD = Array(
                    CSV.read(string(root_new, "mask_L_HD.csv"), DataFrame; header=false)
                )
                mask_M_HD = Array(
                    CSV.read(string(root_new, "mask_M_HD.csv"), DataFrame; header=false)
                )
                mask_S_HD = Array(
                    CSV.read(string(root_new, "mask_S_HD.csv"), DataFrame; header=false)
                )
                mask_L_MD = Array(
                    CSV.read(string(root_new, "mask_L_MD.csv"), DataFrame; header=false)
                )
                mask_M_MD = Array(
                    CSV.read(string(root_new, "mask_M_MD.csv"), DataFrame; header=false)
                )
                mask_S_MD = Array(
                    CSV.read(string(root_new, "mask_S_MD.csv"), DataFrame; header=false)
                )
                mask_L_LD = Array(
                    CSV.read(string(root_new, "mask_L_LD.csv"), DataFrame; header=false)
                )
                mask_M_LD = Array(
                    CSV.read(string(root_new, "mask_M_LD.csv"), DataFrame; header=false)
                )
                mask_S_LD = Array(
                    CSV.read(string(root_new, "mask_S_LD.csv"), DataFrame; header=false)
                )

                # Calibration Prep
                local density_array
                if DENSITY == "low"
                    density_array = [0.025, 0.050, 0.100]
                elseif DENSITY == "normal"
                    density_array = [0.200, 0.400, 0.800]
                end
                array_filtered = abs.(mapwindow(median, calcium_image[:, :, 2], (3, 3)))
                bool_arr = array_filtered .> 0
                bool_arr_erode = (((erode(erode(bool_arr)))))
                c_img = calcium_image[:, :, 1:3]
                mask_cal_3D = Array{Bool}(undef, size(c_img))
                for z in 1:size(c_img, 3)
                    mask_cal_3D[:, :, z] = bool_arr_erode
                end

                hu_calcium = mean(c_img[mask_cal_3D])
                ρ_calcium = 0.2

                arr = masked_array[:, :, 4:6]
                pixel_size = DICOMUtils.get_pixel_size(header)
                voxel_size = pixel_size[1] * pixel_size[2] * pixel_size[3]

                # Score Large Inserts
                ## High Density
                mask_L_HD_3D = Array{Bool}(undef, size(arr))
                for z in 1:size(arr, 3)
                    mask_L_HD_3D[:, :, z] = mask_L_HD
                end
                dilated_mask_L_HD = dilate_mask_large(mask_L_HD_3D)
                ring_mask_L_HD = ring_mask_large(dilated_mask_L_HD)
                hu_heart_tissue_large_hd = mean(arr[ring_mask_L_HD])
                mass_large_hd = score(arr[dilated_mask_L_HD], hu_calcium, hu_heart_tissue_large_hd, voxel_size, ρ_calcium, VolumeFraction())

                ## Medium Density
                mask_L_MD_3D = Array{Bool}(undef, size(arr))
                for z in 1:size(arr, 3)
                    mask_L_MD_3D[:, :, z] = mask_L_MD
                end
                dilated_mask_L_MD = dilate_mask_large(mask_L_MD_3D)
                ring_mask_L_MD = ring_mask_large(dilated_mask_L_MD)
                hu_heart_tissue_large_md = mean(arr[ring_mask_L_MD])
                mass_large_md = score(arr[dilated_mask_L_MD], hu_calcium, hu_heart_tissue_large_md, voxel_size, ρ_calcium, VolumeFraction())

                ## Low Density
                mask_L_LD_3D = Array{Bool}(undef, size(arr))
                for z in 1:size(arr, 3)
                    mask_L_LD_3D[:, :, z] = mask_L_LD
                end
                dilated_mask_L_LD = dilate_mask_large(mask_L_LD_3D)
                ring_mask_L_LD = ring_mask_large(dilated_mask_L_LD)
                hu_heart_tissue_large_ld = mean(arr[ring_mask_L_LD])
                mass_large_ld = score(arr[dilated_mask_L_LD], hu_calcium, hu_heart_tissue_large_ld, voxel_size, ρ_calcium, VolumeFraction())

                # Score Medium Inserts
                ## High Density
                mask_M_HD_3D = Array{Bool}(undef, size(arr))
                for z in 1:size(arr, 3)
                    mask_M_HD_3D[:, :, z] = mask_M_HD
                end
                dilated_mask_M_HD = dilate_mask_medium(mask_M_HD_3D)
                ring_mask_M_HD = ring_mask_medium(dilated_mask_M_HD)
                hu_heart_tissue_medium_hd = mean(arr[ring_mask_M_HD])
                mass_medium_hd = score(arr[dilated_mask_M_HD], hu_calcium, hu_heart_tissue_medium_hd, voxel_size, ρ_calcium, VolumeFraction())

                ## Medium Density
                mask_M_MD_3D = Array{Bool}(undef, size(arr))
                for z in 1:size(arr, 3)
                    mask_M_MD_3D[:, :, z] = mask_M_MD
                end
                dilated_mask_M_MD = dilate_mask_medium(mask_M_MD_3D)
                ring_mask_M_MD = ring_mask_medium(dilated_mask_M_MD)
                hu_heart_tissue_medium_md = mean(arr[ring_mask_M_MD])
                mass_medium_md = score(arr[dilated_mask_M_MD], hu_calcium, hu_heart_tissue_medium_md, voxel_size, ρ_calcium, VolumeFraction())

                ## Low Density
                mask_M_LD_3D = Array{Bool}(undef, size(arr))
                for z in 1:size(arr, 3)
                    mask_M_LD_3D[:, :, z] = mask_M_LD
                end
                dilated_mask_M_LD = dilate_mask_medium(mask_M_LD_3D)
                ring_mask_M_LD = ring_mask_medium(dilated_mask_M_LD)
                hu_heart_tissue_medium_ld = mean(arr[ring_mask_M_LD])
                mass_medium_ld = score(arr[dilated_mask_M_LD], hu_calcium, hu_heart_tissue_medium_ld, voxel_size, ρ_calcium, VolumeFraction())

                # Score Small Inserts
                ## High Density
                mask_S_HD_3D = Array{Bool}(undef, size(arr))
                for z in 1:size(arr, 3)
                    mask_S_HD_3D[:, :, z] = mask_S_HD
                end
                dilated_mask_S_HD = dilate_mask_small(mask_S_HD_3D)
                ring_mask_S_HD = ring_mask_small(dilated_mask_S_HD)
                hu_heart_tissue_small_hd = mean(arr[ring_mask_S_HD])
                mass_small_hd = score(arr[dilated_mask_S_HD], hu_calcium, hu_heart_tissue_large_hd, voxel_size, ρ_calcium, VolumeFraction())

                ## Medium Density
                mask_S_MD_3D = Array{Bool}(undef, size(arr))
                for z in 1:size(arr, 3)
                    mask_S_MD_3D[:, :, z] = mask_S_MD
                end
                dilated_mask_S_MD = dilate_mask_small(mask_S_MD_3D)
                ring_mask_S_MD = ring_mask_small(dilated_mask_S_MD)
                hu_heart_tissue_small_md = mean(arr[ring_mask_S_MD])
                mass_small_md = score(arr[dilated_mask_S_MD], hu_calcium, hu_heart_tissue_large_md, voxel_size, ρ_calcium, VolumeFraction())

                ## Low Density
                mask_S_LD_3D = Array{Bool}(undef, size(arr))
                for z in 1:size(arr, 3)
                    mask_S_LD_3D[:, :, z] = mask_S_LD
                end
                dilated_mask_S_LD = dilate_mask_small(mask_S_LD_3D)
                ring_mask_S_LD = ring_mask_small(dilated_mask_S_LD)
                hu_heart_tissue_small_ld = mean(arr[ring_mask_S_LD])
                mass_small_ld = score(arr[dilated_mask_S_LD], hu_calcium, hu_heart_tissue_large_ld, voxel_size, ρ_calcium, VolumeFraction())

                # Results
                local density_array_calculations
                if DENSITY == "low"
                    density_array_calculations = [0, 25, 50, 100]
                elseif DENSITY == "normal"
                    density_array_calculations = [0, 200, 400, 800]
                end

                inserts = ["Low Density", "Medium Density", "High Density"]
                volume_gt = [7.065, 63.585, 176.625]
                ground_truth_mass_large = [
                    volume_gt[3] * density_array_calculations[2] * 1e-3,
                    volume_gt[3] * density_array_calculations[3] * 1e-3,
                    volume_gt[3] * density_array_calculations[4] * 1e-3,
                ] # mg

                calculated_mass_large = [mass_large_ld, mass_large_md, mass_large_hd]
                ground_truth_mass_medium = [
                    volume_gt[2] * density_array_calculations[2] * 1e-3,
                    volume_gt[2] * density_array_calculations[3] * 1e-3,
                    volume_gt[2] * density_array_calculations[4] * 1e-3,
                ]
                calculated_mass_medium = [mass_medium_ld, mass_medium_md, mass_medium_hd]
                ground_truth_mass_small = [
                    volume_gt[1] * density_array_calculations[2] * 1e-3,
                    volume_gt[1] * density_array_calculations[3] * 1e-3,
                    volume_gt[1] * density_array_calculations[4] * 1e-3,
                ]
                calculated_mass_small = [mass_small_ld, mass_small_md, mass_small_hd]

                df = DataFrame(;
                    VENDER=VENDER,
                    SIZE=SIZE,
                    DENSITY=DENSITY,
                    scan=scan,
                    inserts=inserts,
                    ground_truth_mass_large=ground_truth_mass_large,
                    calculated_mass_large=calculated_mass_large,
                    ground_truth_mass_medium=ground_truth_mass_medium,
                    calculated_mass_medium=calculated_mass_medium,
                    ground_truth_mass_small=ground_truth_mass_small,
                    calculated_mass_small=calculated_mass_small
                )
                push!(dfs_vf, df)

                # if DENSITY == "normal"
                #     # Calculate mean intensities for calibration
                #     erode_mask_L_HD = erode(erode(mask_L_HD))
                #     mean_HD = mean(single_arr[erode_mask_L_HD])
                #     push!(high_dens, mean_HD)

                #     erode_mask_L_MD = erode(erode(mask_L_MD))
                #     mean_MD = mean(single_arr[erode_mask_L_MD])
                #     push!(med_dens, mean_MD)

                #     erode_mask_L_LD = erode(erode(mask_L_LD))
                #     mean_LD = mean(single_arr[erode_mask_L_LD])
                #     push!(low_dens, mean_LD)
                # else
                #     erode_mask_L_HD = erode(erode(mask_L_HD))
                #     mean_HD = mean(single_arr[erode_mask_L_HD])
                #     push!(high_dens100, mean_HD)

                #     erode_mask_L_MD = erode(erode(mask_L_MD))
                #     mean_MD = mean(single_arr[erode_mask_L_MD])
                #     push!(med_dens50, mean_MD)

                #     erode_mask_L_LD = erode(erode(mask_L_LD))
                #     mean_LD = mean(single_arr[erode_mask_L_LD])
                #     push!(low_dens25, mean_LD)
                # end

				#---------------- SWCS ----------------#
				local μ, σ
                if VENDER == "80"
                    μ, σ = 170, 40
                elseif VENDER == "100"
                    μ, σ = 165, 40
                elseif VENDER == "120"
                    μ, σ = 160, 40
                else
                    VENDER == "135"
                    μ, σ = 155, 40
                end

				# Mask Calibration Factor
                output = calc_output(
                    masked_array, header, slice_CCI, thresh, trues(3, 3)
                )
                insert_centers = calc_centers_simulation(
                    dcm_array, output, header, center_insert, slice_CCI
                )
                rows, cols = Int(header[tag"Rows"]), Int(header[tag"Columns"])
                pixel_size = DICOMUtils.get_pixel_size(header)
                mass_cal_factor, angle_0_200HA, water_rod_metrics = mass_calibration(
                    masked_array,
                    insert_centers[:Large_LD],
                    center_insert,
                    2,
                    cols,
                    rows,
                    pixel_size,
                )

                # Background
                arr = masked_array[:, :, 4:6]
                alg2 = SpatiallyWeighted()

                background_mask = zeros(size(arr)...)
                background_mask[
                    (center_insert[1]-5):(center_insert[1]+5),
                    (center_insert[2]-5):(center_insert[2]+5),
                    2,
                ] .= 1
                background_mask = Bool.(background_mask)
                swcs_bkg = score(background_mask, μ, σ, alg2)

                # Score Large Inserts
                ## High Density
                mask_L_HD_3D = Array{Bool}(undef, size(arr))
                for z in 1:size(arr, 3)
                    mask_L_HD_3D[:, :, z] = mask_L_HD
                end
                dilated_mask_L_HD = dilate(dilate(mask_L_HD_3D))
                alg = Agatston()
                overlayed_mask_l_hd = create_mask(arr, dilated_mask_L_HD)
                agat_l_hd, mass_l_hd = score(
                    overlayed_mask_l_hd, pixel_size, mass_cal_factor, alg
                )
                swcs_l_hd = score(overlayed_mask_l_hd, μ, σ, alg2)

                ## Medium Density
                mask_L_MD_3D = Array{Bool}(undef, size(arr))
                for z in 1:size(arr, 3)
                    mask_L_MD_3D[:, :, z] = mask_L_MD
                end
                dilated_mask_L_MD = dilate(dilate(mask_L_MD_3D))
                overlayed_mask_l_md = create_mask(arr, dilated_mask_L_MD)
                agat_l_md, mass_l_md = score(
                    overlayed_mask_l_md, pixel_size, mass_cal_factor, alg
                )
                swcs_l_md = score(overlayed_mask_l_md, μ, σ, alg2)

                ## Low Density
                mask_L_LD_3D = Array{Bool}(undef, size(arr))
                for z in 1:size(arr, 3)
                    mask_L_LD_3D[:, :, z] = mask_L_LD
                end
                dilated_mask_L_LD = dilate(dilate(mask_L_LD_3D))
                overlayed_mask_l_ld = create_mask(arr, dilated_mask_L_LD)
                agat_l_ld, mass_l_ld = score(
                    overlayed_mask_l_ld, pixel_size, mass_cal_factor, alg
                )
                swcs_l_ld = score(overlayed_mask_l_ld, μ, σ, alg2)

                # Score Medium Inserts
                ## High Density
                mask_M_HD_3D = Array{Bool}(undef, size(arr))
                for z in 1:size(arr, 3)
                    mask_M_HD_3D[:, :, z] = mask_M_HD
                end
                dilated_mask_M_HD = dilate(dilate(dilate(dilate(mask_M_HD_3D))))
                overlayed_mask_m_hd = create_mask(arr, dilated_mask_M_HD)
                agat_m_hd, mass_m_hd = score(
                    overlayed_mask_m_hd, pixel_size, mass_cal_factor, alg
                )
                swcs_m_hd = score(overlayed_mask_m_hd, μ, σ, alg2)

                ## Medium Density
                mask_M_MD_3D = Array{Bool}(undef, size(arr))
                for z in 1:size(arr, 3)
                    mask_M_MD_3D[:, :, z] = mask_M_MD
                end
                dilated_mask_M_MD = dilate(dilate(dilate(dilate(mask_M_MD_3D))))
                overlayed_mask_m_md = create_mask(arr, dilated_mask_M_MD)
                agat_m_md, mass_m_md = score(
                    overlayed_mask_m_md, pixel_size, mass_cal_factor, alg
                )
                swcs_m_md = score(overlayed_mask_m_md, μ, σ, alg2)

                ## Low Density
                mask_M_LD_3D = Array{Bool}(undef, size(arr))
                for z in 1:size(arr, 3)
                    mask_M_LD_3D[:, :, z] = mask_M_LD
                end
                dilated_mask_M_LD = dilate(dilate(dilate(dilate(mask_M_LD_3D))))
                overlayed_mask_m_ld = create_mask(arr, dilated_mask_M_LD)
                agat_m_ld, mass_m_ld = score(
                    overlayed_mask_m_ld, pixel_size, mass_cal_factor, alg
                )
                swcs_m_ld = score(overlayed_mask_m_ld, μ, σ, alg2)

                # Score Small Inserts
                ## High Density
                mask_S_HD_3D = Array{Bool}(undef, size(arr))
                for z in 1:size(arr, 3)
                    mask_S_HD_3D[:, :, z] = mask_S_HD
                end
                dilated_mask_S_HD = dilate((dilate(dilate(dilate((mask_S_HD_3D))))))
                overlayed_mask_s_hd = create_mask(arr, dilated_mask_S_HD)
                agat_s_hd, mass_s_hd = score(
                    overlayed_mask_s_hd, pixel_size, mass_cal_factor, alg
                )
                swcs_s_hd = score(overlayed_mask_s_hd, μ, σ, alg2)

                ## Medium Density
                mask_S_MD_3D = Array{Bool}(undef, size(arr))
                for z in 1:size(arr, 3)
                    mask_S_MD_3D[:, :, z] = mask_S_MD
                end
                dilated_mask_S_MD = dilate((dilate(dilate(dilate(mask_S_MD_3D)))))
                overlayed_mask_s_md = create_mask(arr, dilated_mask_S_MD)
                agat_s_md, mass_s_md = score(
                    overlayed_mask_s_md, pixel_size, mass_cal_factor, alg
                )
                swcs_s_md = score(overlayed_mask_s_md, μ, σ, alg2)

                ## Low Density
                mask_S_LD_3D = Array{Bool}(undef, size(arr))
                for z in 1:size(arr, 3)
                    mask_S_LD_3D[:, :, z] = mask_S_LD
                end
                dilated_mask_S_LD = dilate((dilate(dilate(dilate(mask_S_LD_3D)))))
                overlayed_mask_s_ld = create_mask(arr, dilated_mask_S_LD)
                agat_s_ld, mass_s_ld = score(
                    overlayed_mask_s_ld, pixel_size, mass_cal_factor, alg
                )
                swcs_s_ld = score(overlayed_mask_s_ld, μ, σ, alg2)

                # Results
                local density_array
                if DENSITY == "low"
                    density_array = [0, 25, 50, 100]
                elseif DENSITY == "normal"
                    density_array = [0, 200, 400, 800]
                end
                inserts = ["Low Density", "Medium Density", "High Density"]

                ## Agatston
                calculated_agat_large = [agat_l_ld, agat_l_md, agat_l_hd]
                calculated_agat_medium = [agat_m_ld, agat_m_md, agat_m_hd]
                calculated_agat_small = [agat_s_ld, agat_s_md, agat_s_hd]

                ## Mass
                volume_gt = [7.065, 63.585, 176.625]
                ground_truth_mass_large = [
                    volume_gt[3] * density_array[2] * 1e-3,
                    volume_gt[3] * density_array[3] * 1e-3,
                    volume_gt[3] * density_array[4] * 1e-3,
                ] # mg
                calculated_mass_large = [mass_l_ld, mass_l_md, mass_l_hd]
                ground_truth_mass_medium = [
                    volume_gt[2] * density_array[2] * 1e-3,
                    volume_gt[2] * density_array[3] * 1e-3,
                    volume_gt[2] * density_array[4] * 1e-3,
                ]
                calculated_mass_medium = [mass_m_ld, mass_m_md, mass_m_hd]
                ground_truth_mass_small = [
                    volume_gt[1] * density_array[2] * 1e-3,
                    volume_gt[1] * density_array[3] * 1e-3,
                    volume_gt[1] * density_array[4] * 1e-3,
                ]
                calculated_mass_small = [mass_s_ld, mass_s_md, mass_s_hd]
                calculated_swcs_large = [swcs_l_ld, swcs_l_md, swcs_l_hd]
                calculated_swcs_medium = [swcs_m_ld, swcs_m_md, swcs_m_hd]
                calculated_swcs_small = [swcs_s_ld, swcs_s_md, swcs_s_hd]

                df = DataFrame(;
                    VENDER=VENDER,
                    SIZE=SIZE,
                    DENSITY=DENSITY,
                    scan=scan,
                    inserts=inserts,
                    ground_truth_mass_large=ground_truth_mass_large,
                    calculated_mass_large=calculated_mass_large,
                    ground_truth_mass_medium=ground_truth_mass_medium,
                    calculated_mass_medium=calculated_mass_medium,
                    ground_truth_mass_small=ground_truth_mass_small,
                    calculated_mass_small=calculated_mass_small,
                    calculated_swcs_large=calculated_swcs_large,
                    calculated_swcs_medium=calculated_swcs_medium,
                    calculated_swcs_small=calculated_swcs_small,
                    swcs_bkg=swcs_bkg
                )
                push!(dfs_s, df)

				#---------------- Integrated ----------------#
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

                local intensity_array
                local density_array
                if VENDER == "80"
                    intensity_array = [0, 69.7, 377.3, 1375.3]
                    density_array = [0, 25, 200, 800]
                elseif VENDER == "100"
                    intensity_array = [0, 63.3, 326.0, 1179.4]
                    density_array = [0, 25, 200, 800]
                elseif VENDER == "120"
                    intensity_array = [0, 58.1, 296.9, 1072.1]
                    density_array = [0, 25, 200, 800]
                else
                    intensity_array = [0, 55.5, 282.7, 1019.7]
                    density_array = [0, 25, 200, 800]
                end

                df_cal = DataFrame(
                    :density => density_array, :intensity => intensity_array
                )

                linearRegressor = lm(@formula(intensity ~ density), df_cal)
                linearFit = GLM.predict(linearRegressor)
                m = linearRegressor.model.pp.beta0[2]
                b = linearRegressor.model.rr.mu[1]
                density(intensity) = (intensity - b) / m
                intensity(ρ) = m * ρ + b

                # Score background
                background_mask = zeros(size(arr)...)
                background_mask[
                    (center_insert[1]-5):(center_insert[1]+5),
                    (center_insert[2]-5):(center_insert[2]+5),
                    2,
                ] .= 1
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
                ring_mask_L_HD =
                    dilate(dilate(dilate(dilate(mask_L_HD_3D)))) -
                    dilate(dilate(dilate(mask_L_HD_3D)))
                single_ring_mask_L_HD = Bool.(ring_mask_L_HD[:, :, 3])
                s_bkg_L_HD = mean(single_arr[single_ring_mask_L_HD])
                S_Obj_HD = intensity(800)
                ρ_hd = 0.8 # mg/mm^3
                alg_L_HD = Integrated(arr[mask_L_HD_3D])
                mass_l_hd = score(s_bkg_L_HD, S_Obj_HD, pixel_size, ρ_hd, alg_L_HD)

                ## Medium Density

                dilated_mask_L_MD = dilate(dilate(mask_L_MD_3D))
                ring_mask_L_MD =
                    dilate(dilate(dilate(dilate(mask_L_MD_3D)))) -
                    dilate(dilate(dilate(mask_L_MD_3D)))
                single_ring_mask_L_MD = Bool.(ring_mask_L_MD[:, :, 3])
                s_bkg_L_MD = mean(single_arr[single_ring_mask_L_MD])
                S_Obj_MD = intensity(400)
                ρ_md = 0.4 # mg/mm^3
                alg_L_MD = Integrated(arr[mask_L_MD_3D])
                mass_l_md = score(s_bkg_L_MD, S_Obj_MD, pixel_size, ρ_md, alg_L_MD)

                ## Low Density

                dilated_mask_L_LD = dilate(dilate(mask_L_LD_3D))
                ring_mask_L_LD =
                    dilate(dilate(dilate(dilate(mask_L_LD_3D)))) -
                    dilate(dilate(dilate(mask_L_LD_3D)))
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
                ring_mask_M_HD =
                    dilate(dilate(dilate(dilate(dilate(mask_M_HD_3D))))) -
                    dilate(dilate(dilate(dilate(mask_M_HD_3D))))
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
                ring_mask_M_MD =
                    dilate(dilate(dilate(dilate(dilate(dilate(mask_M_MD_3D)))))) -
                    dilate(dilate(dilate(dilate(dilate(mask_M_MD_3D)))))
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
                ring_mask_M_LD =
                    dilate(dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))))) -
                    dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))))
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
                dilated_mask_S_HD = dilate(
                    dilate(dilate(dilate(dilate((mask_S_HD_3D)))))
                )
                ring_mask_S_HD =
                    dilate(dilate(dilate(dilate(dilate(mask_S_HD_3D))))) -
                    dilate(dilate(dilate(dilate(mask_S_HD_3D))))
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
                ring_mask_S_MD =
                    dilate(dilate(dilate(dilate(dilate(mask_S_MD_3D))))) -
                    dilate(dilate(dilate(dilate(mask_S_MD_3D))))
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
                ring_mask_S_LD =
                    dilate(dilate(dilate(dilate(dilate(mask_S_LD_3D))))) -
                    dilate(dilate(dilate(dilate(mask_S_LD_3D))))
                single_ring_mask_S_LD = Bool.(ring_mask_S_LD[:, :, 3])
                s_bkg_S_LD = mean(single_arr[single_ring_mask_S_LD])
                alg_S_LD = Integrated(arr[mask_S_LD_3D])
                mass_s_ld = score(s_bkg_S_LD, S_Obj_LD, pixel_size, ρ_ld, alg_S_LD)

                # Results
                local density_array_calculations
                if DENSITY == "low"
                    density_array_calculations = [0, 25, 50, 100]
                elseif DENSITY == "normal"
                    density_array_calculations = [0, 200, 400, 800]
                end

                inserts = ["Low Density", "Medium Density", "High Density"]
                volume_gt = [7.065, 63.585, 176.625]
                ground_truth_mass_large = [
                    volume_gt[3] * density_array_calculations[2] * 1e-3,
                    volume_gt[3] * density_array_calculations[3] * 1e-3,
                    volume_gt[3] * density_array_calculations[4] * 1e-3,
                ] # mg

                calculated_mass_large = [mass_l_ld, mass_l_md, mass_l_hd]
                ground_truth_mass_medium = [
                    volume_gt[2] * density_array_calculations[2] * 1e-3,
                    volume_gt[2] * density_array_calculations[3] * 1e-3,
                    volume_gt[2] * density_array_calculations[4] * 1e-3,
                ]
                calculated_mass_medium = [mass_m_ld, mass_m_md, mass_m_hd]
                ground_truth_mass_small = [
                    volume_gt[1] * density_array_calculations[2] * 1e-3,
                    volume_gt[1] * density_array_calculations[3] * 1e-3,
                    volume_gt[1] * density_array_calculations[4] * 1e-3,
                ]
                calculated_mass_small = [mass_s_ld, mass_s_md, mass_s_hd]

                df = DataFrame(;
                    VENDER=VENDER,
                    SIZE=SIZE,
                    DENSITY=DENSITY,
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

                if DENSITY == "normal"
                    # Calculate mean intensities for calibration
                    erode_mask_L_HD = erode(erode(mask_L_HD))
                    mean_HD = mean(single_arr[erode_mask_L_HD])
                    push!(high_dens, mean_HD)

                    erode_mask_L_MD = erode(erode(mask_L_MD))
                    mean_MD = mean(single_arr[erode_mask_L_MD])
                    push!(med_dens, mean_MD)

                    erode_mask_L_LD = erode(erode(mask_L_LD))
                    mean_LD = mean(single_arr[erode_mask_L_LD])
                    push!(low_dens, mean_LD)
                else
                    erode_mask_L_HD = erode(erode(mask_L_HD))
                    mean_HD = mean(single_arr[erode_mask_L_HD])
                    push!(high_dens100, mean_HD)

                    erode_mask_L_MD = erode(erode(mask_L_MD))
                    mean_MD = mean(single_arr[erode_mask_L_MD])
                    push!(med_dens50, mean_MD)

                    erode_mask_L_LD = erode(erode(mask_L_LD))
                    mean_LD = mean(single_arr[erode_mask_L_LD])
                    push!(low_dens25, mean_LD)
                end

				#---------------- Agatston ----------------#
				# Mask Calibration Factor
                output = calc_output(
                    masked_array, header, slice_CCI, thresh, trues(3, 3)
                )
                insert_centers = calc_centers_simulation(
                    dcm_array, output, header, center_insert, slice_CCI
                )
                rows, cols = Int(header[tag"Rows"]), Int(header[tag"Columns"])
                pixel_size = DICOMUtils.get_pixel_size(header)
                mass_cal_factor, angle_0_200HA, water_rod_metrics = mass_calibration(
                    masked_array,
                    insert_centers[:Large_LD],
                    center_insert,
                    2,
                    cols,
                    rows,
                    pixel_size,
                )

                local agat_thresh
                agat_thresh = 130

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
                agat_l_hd, mass_l_hd = score(
                    overlayed_mask_l_hd,
                    pixel_size,
                    mass_cal_factor,
                    alg;
                    kV=kV
                )

                ## Medium Density
                mask_L_MD_3D = Array{Bool}(undef, size(arr))
                for z in 1:size(arr, 3)
                    mask_L_MD_3D[:, :, z] = mask_L_MD
                end
                dilated_mask_L_MD = dilate(dilate(mask_L_MD_3D))
                overlayed_mask_l_md = create_mask(arr, dilated_mask_L_MD)
                agat_l_md, mass_l_md = score(
                    overlayed_mask_l_md,
                    pixel_size,
                    mass_cal_factor,
                    alg;
                    kV=kV
                )

                ## Low Density
                mask_L_LD_3D = Array{Bool}(undef, size(arr))
                for z in 1:size(arr, 3)
                    mask_L_LD_3D[:, :, z] = mask_L_LD
                end
                dilated_mask_L_LD = dilate(dilate(mask_L_LD_3D))
                overlayed_mask_l_ld = create_mask(arr, dilated_mask_L_LD)
                agat_l_ld, mass_l_ld = score(
                    overlayed_mask_l_ld,
                    pixel_size,
                    mass_cal_factor,
                    alg;
                    kV=kV
                )

                # Score Medium Inserts
                ## High Density
                mask_M_HD_3D = Array{Bool}(undef, size(arr))
                for z in 1:size(arr, 3)
                    mask_M_HD_3D[:, :, z] = mask_M_HD
                end
                dilated_mask_M_HD = dilate(dilate(dilate(dilate(mask_M_HD_3D))))
                overlayed_mask_m_hd = create_mask(arr, dilated_mask_M_HD)
                agat_m_hd, mass_m_hd = score(
                    overlayed_mask_m_hd,
                    pixel_size,
                    mass_cal_factor,
                    alg;
                    kV=kV
                )

                ## Medium Density
                mask_M_MD_3D = Array{Bool}(undef, size(arr))
                for z in 1:size(arr, 3)
                    mask_M_MD_3D[:, :, z] = mask_M_MD
                end
                dilated_mask_M_MD = dilate(dilate(dilate(dilate(mask_M_MD_3D))))
                overlayed_mask_m_md = create_mask(arr, dilated_mask_M_MD)
                agat_m_md, mass_m_md = score(
                    overlayed_mask_m_md,
                    pixel_size,
                    mass_cal_factor,
                    alg;
                    kV=kV
                )

                ## Low Density
                mask_M_LD_3D = Array{Bool}(undef, size(arr))
                for z in 1:size(arr, 3)
                    mask_M_LD_3D[:, :, z] = mask_M_LD
                end
                dilated_mask_M_LD = dilate(dilate(dilate(dilate(mask_M_LD_3D))))
                overlayed_mask_m_ld = create_mask(arr, dilated_mask_M_LD)
                agat_m_ld, mass_m_ld = score(
                    overlayed_mask_m_ld,
                    pixel_size,
                    mass_cal_factor,
                    alg;
                    kV=kV
                )

                # Score Small Inserts
                ## High Density
                mask_S_HD_3D = Array{Bool}(undef, size(arr))
                for z in 1:size(arr, 3)
                    mask_S_HD_3D[:, :, z] = mask_S_HD
                end
                dilated_mask_S_HD = dilate((dilate(dilate(dilate((mask_S_HD_3D))))))
                overlayed_mask_s_hd = create_mask(arr, dilated_mask_S_HD)
                agat_s_hd, mass_s_hd = score(
                    overlayed_mask_s_hd,
                    pixel_size,
                    mass_cal_factor,
                    alg;
                    kV=kV
                )

                ## Medium Density
                mask_S_MD_3D = Array{Bool}(undef, size(arr))
                for z in 1:size(arr, 3)
                    mask_S_MD_3D[:, :, z] = mask_S_MD
                end
                dilated_mask_S_MD = dilate((dilate(dilate(dilate(mask_S_MD_3D)))))
                overlayed_mask_s_md = create_mask(arr, dilated_mask_S_MD)
                agat_s_md, mass_s_md = score(
                    overlayed_mask_s_md,
                    pixel_size,
                    mass_cal_factor,
                    alg;
                    kV=kV
                )

                ## Low Density
                mask_S_LD_3D = Array{Bool}(undef, size(arr))
                for z in 1:size(arr, 3)
                    mask_S_LD_3D[:, :, z] = mask_S_LD
                end
                dilated_mask_S_LD = dilate((dilate(dilate(dilate(mask_S_LD_3D)))))
                overlayed_mask_s_ld = create_mask(arr, dilated_mask_S_LD)
                agat_s_ld, mass_s_ld = score(
                    overlayed_mask_s_ld,
                    pixel_size,
                    mass_cal_factor,
                    alg;
                    kV=kV
                )

                # Results
                local density_array
                if DENSITY == "low"
                    density_array = [0, 25, 50, 100]
                elseif DENSITY == "normal"
                    density_array = [0, 200, 400, 800]
                end
                inserts = ["Low Density", "Medium Density", "High Density"]

                ## Agatston
                calculated_agat_large = [agat_l_ld, agat_l_md, agat_l_hd]
                calculated_agat_medium = [agat_m_ld, agat_m_md, agat_m_hd]
                calculated_agat_small = [agat_s_ld, agat_s_md, agat_s_hd]

                ## Mass
                volume_gt = [7.065, 63.585, 176.625]
                ground_truth_mass_large = [
                    volume_gt[3] * density_array[2] * 1e-3,
                    volume_gt[3] * density_array[3] * 1e-3,
                    volume_gt[3] * density_array[4] * 1e-3,
                ] # mg
                calculated_mass_large = [mass_l_ld, mass_l_md, mass_l_hd]
                ground_truth_mass_medium = [
                    volume_gt[2] * density_array[2] * 1e-3,
                    volume_gt[2] * density_array[3] * 1e-3,
                    volume_gt[2] * density_array[4] * 1e-3,
                ]
                calculated_mass_medium = [mass_m_ld, mass_m_md, mass_m_hd]
                ground_truth_mass_small = [
                    volume_gt[1] * density_array[2] * 1e-3,
                    volume_gt[1] * density_array[3] * 1e-3,
                    volume_gt[1] * density_array[4] * 1e-3,
                ]
                calculated_mass_small = [mass_s_ld, mass_s_md, mass_s_hd]

                df = DataFrame(;
                    VENDER=VENDER,
                    SIZE=SIZE,
                    DENSITY=DENSITY,
                    scan=scan,
                    inserts=inserts,
                    ground_truth_mass_large=ground_truth_mass_large,
                    calculated_mass_large=calculated_mass_large,
                    ground_truth_mass_medium=ground_truth_mass_medium,
                    calculated_mass_medium=calculated_mass_medium,
                    ground_truth_mass_small=ground_truth_mass_small,
                    calculated_mass_small=calculated_mass_small
                )
                push!(dfs_a, df)
            end
        end
    end
end

# ╔═╡ 8e234939-154a-486d-a5ed-9fc9f2749965
md"""
# Save Results
"""

# ╔═╡ 1e6799f1-985c-46c0-adb3-089e40d19f0d
dfs_vf

# ╔═╡ dcc1ffd9-bd83-4ec8-8e0f-7d9aa659b2be
dfs_s

# ╔═╡ 6771b626-bb25-4d6b-945d-6b668c4043ee
dfs_i

# ╔═╡ 889d6ddd-97ba-49a3-940a-3f4a467ac721
dfs_a

# ╔═╡ 38e4f88d-9e5c-4d31-b0f6-4c79aa06baad
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

# ╔═╡ af147019-b7f3-4b05-ae10-32ed683ab97c
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
# ╠═c4e1c65a-f078-4633-b785-78f907bb9dcc
# ╠═61474b4a-9139-4021-b41e-15e885f9466a
# ╠═0028f871-502e-43dd-a43a-624c032f5c8e
# ╠═0f466539-11a6-4d6d-b579-53ffb9a25a39
# ╠═c2e4696a-e2e2-4f02-ab06-9ebe62732d0b
# ╠═bcb96a7f-7553-426e-a1f5-317bf55d1856
# ╠═567e50fb-3d49-4136-aa32-a0ff38a61c88
# ╠═fdb032f0-3949-423f-9bd9-08b53d65f243
# ╟─e73c0eea-3ac4-42ad-a118-040729af29cc
# ╟─cb23368e-79f8-4db0-b272-430956751ff6
# ╟─10a9aa52-a8ed-4e76-ba90-43d8c77a33b7
# ╟─715a8838-a7f0-437d-aa5f-68744e7ca668
# ╟─facd5777-30ae-4869-9ccc-5e4128cc2de7
# ╟─e7626765-d2ca-4129-8b47-f5925a2c0eb1
# ╟─e6387a3a-2377-4761-a6c4-77e375afb372
# ╠═b26eab05-14e1-40ad-980b-eb8f5206542b
# ╟─8e234939-154a-486d-a5ed-9fc9f2749965
# ╠═1e6799f1-985c-46c0-adb3-089e40d19f0d
# ╠═dcc1ffd9-bd83-4ec8-8e0f-7d9aa659b2be
# ╠═6771b626-bb25-4d6b-945d-6b668c4043ee
# ╠═889d6ddd-97ba-49a3-940a-3f4a467ac721
# ╠═38e4f88d-9e5c-4d31-b0f6-4c79aa06baad
# ╠═af147019-b7f3-4b05-ae10-32ed683ab97c

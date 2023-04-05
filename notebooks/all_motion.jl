### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 616a3e0e-82e8-4b08-ac82-1652403a693e
# ╠═╡ show_logs = false
begin
    using Pkg
    Pkg.activate(".")

    using PlutoUI, Statistics, CSV, DataFrames, GLM, CairoMakie, HypothesisTests, Colors, MLJBase, DICOM, DICOMUtils, PhantomSegmentation, CalciumScoring, ImageMorphology, ImageFiltering, Noise
    using StatsBase: quantile!, rmsd
end

# ╔═╡ 40c1b1e3-e103-48d5-8e9a-dd408ec162c6
TableOfContents()

# ╔═╡ e8ce5c83-98b6-4ec4-beac-677fd24ce62d
IMAGES = "images_new"

# ╔═╡ c7a8708a-14b3-418c-9a26-a45f8b84711b
OUTPUT = "output_new"

# ╔═╡ 992dca9e-3664-4798-a435-72ecd4c797d7
SAVE_DF = "motion.csv"

# ╔═╡ 1af7db2f-126d-4c5e-8091-7b937744b22c
VENDORS = ["80", "100", "120", "135"]

# ╔═╡ 3fa8dba1-e564-410c-9ce5-36870795f493
SIZES = ["small", "medium", "large"]

# ╔═╡ 0d2127b5-489a-401e-a7eb-9c2f5a0aa22a
DENSITIES = ["low", "normal"]

# ╔═╡ 1f88f342-7418-4d8a-af66-53a985f10209
function create_mask(array, mask)
    @assert size(array) == size(mask)
    idxs = findall(x -> x == true, mask)
    overlayed_mask = zeros(size(array))
    for idx in idxs
        overlayed_mask[idx] = array[idx]
    end
    return overlayed_mask
end

# ╔═╡ 2e03d03e-c048-48ce-8b5a-449bea68f258
function dilate_mask_large(mask)
    return dilate(mask)
end

# ╔═╡ e5e63f2b-38de-41bc-9785-576d8f3bb8d3
function ring_mask_large(dilated_mask)
    return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask)
end

# ╔═╡ 766dcf5a-b778-423a-b39b-610299fcf910
function dilate_mask_medium(mask)
    return (mask)
end

# ╔═╡ 39eebc5b-e16a-43ae-abdd-19e44813e9c8
function ring_mask_medium(dilated_mask)
    return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask)
end

# ╔═╡ 872ec5e3-6001-40a1-b5c6-075b3560d34d
function dilate_mask_small(mask)
    return (mask)
end

# ╔═╡ 26fee6aa-cbf1-475b-8553-8dab04e0f80d
function ring_mask_small(dilated_mask)
    return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask)
end

# ╔═╡ 066702c9-afa8-49e0-81bf-d2bfcdfd2c87
function dilate_mask_large_bkg(mask)
    return dilate(dilate(mask))
end

# ╔═╡ 485030f6-e3e8-477f-bda3-e1b0b321df54
function dilate_mask_medium_bkg(mask)
    return dilate(mask)
end

# ╔═╡ e5609112-e6ce-4fba-9237-c67941fc2fe0
function dilate_mask_small_bkg(mask)
    return (mask)
end

# ╔═╡ 73f69fdd-20e7-4aa6-9a00-942ab061433c
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
        for SIZE in SIZES
            for DENSITY in DENSITIES
                #---------------- Reusable Pieces ----------------#
				BASE_PATH = joinpath(dirname(pwd()), IMAGES, SIZE, DENSITY)
                root_path = joinpath(BASE_PATH, VENDOR * "-motion")
                dcm_path_list = dcm_list_builder(root_path)
                pth = dcm_path_list[1]
                scan = basename(pth)
                header, dcm_array, slice_thick_ori1 = dcm_reader(pth)
                kV = parse(Int64, VENDOR)

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

                @info VENDOR, SIZE, DENSITY

                calcium_image, slice_CCI, quality_slice, cal_rod_slice = mask_rod(
                    masked_array, header; calcium_threshold=thresh
                )

                # Load Masks
				root_new = joinpath(dirname(pwd()), "julia_arrays", SIZE)
                mask_L_HD = Array(CSV.read(joinpath(root_new, "mask_L_HD.csv"), DataFrame; header=false))
                mask_M_HD = Array(CSV.read(joinpath(root_new, "mask_M_HD.csv"), DataFrame; header=false))
                mask_S_HD = Array(CSV.read(joinpath(root_new, "mask_S_HD.csv"), DataFrame; header=false))
                mask_L_MD = Array(CSV.read(joinpath(root_new, "mask_L_MD.csv"), DataFrame; header=false))
                mask_M_MD = Array(CSV.read(joinpath(root_new, "mask_M_MD.csv"), DataFrame; header=false))
                mask_S_MD = Array(CSV.read(joinpath(root_new, "mask_S_MD.csv"), DataFrame; header=false))
                mask_L_LD = Array(CSV.read(joinpath(root_new, "mask_L_LD.csv"), DataFrame; header=false))
                mask_M_LD = Array(CSV.read(joinpath(root_new, "mask_M_LD.csv"), DataFrame; header=false))
                mask_S_LD = Array(CSV.read(joinpath(root_new, "mask_S_LD.csv"), DataFrame; header=false))

                # Prepare Mask Arrays
                arr = masked_array[:, :, 4:6]
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
                for z in 1:size(arr, 3)
                    mask_L_HD_3D[:, :, z] = mask_L_HD
                end
                dilated_mask_L_HD = dilate_mask_large(mask_L_HD_3D)
                ring_mask_L_HD = ring_mask_large(dilated_mask_L_HD)

                ### Medium Density
                mask_L_MD_3D = Array{Bool}(undef, size(arr))
                for z in 1:size(arr, 3)
                    mask_L_MD_3D[:, :, z] = mask_L_MD
                end
                dilated_mask_L_MD = dilate_mask_large(mask_L_MD_3D)
                ring_mask_L_MD = ring_mask_large(dilated_mask_L_MD)

                ### Low Density
                mask_L_LD_3D = Array{Bool}(undef, size(arr))
                for z in 1:size(arr, 3)
                    mask_L_LD_3D[:, :, z] = mask_L_LD
                end
                dilated_mask_L_LD = dilate_mask_large(mask_L_LD_3D)
                ring_mask_L_LD = ring_mask_large(dilated_mask_L_LD)


                ## Medium 
                ### High Density
                mask_M_HD_3D = Array{Bool}(undef, size(arr))
                for z in 1:size(arr, 3)
                    mask_M_HD_3D[:, :, z] = mask_M_HD
                end
                dilated_mask_M_HD = dilate_mask_medium(mask_M_HD_3D)
                ring_mask_M_HD = ring_mask_medium(dilated_mask_M_HD)

                ### Medium Density
                mask_M_MD_3D = Array{Bool}(undef, size(arr))
                for z in 1:size(arr, 3)
                    mask_M_MD_3D[:, :, z] = mask_M_MD
                end
                dilated_mask_M_MD = dilate_mask_medium(mask_M_MD_3D)
                ring_mask_M_MD = ring_mask_medium(dilated_mask_M_MD)

                ### Low Density
                mask_M_LD_3D = Array{Bool}(undef, size(arr))
                for z in 1:size(arr, 3)
                    mask_M_LD_3D[:, :, z] = mask_M_LD
                end
                dilated_mask_M_LD = dilate_mask_medium(mask_M_LD_3D)
                ring_mask_M_LD = ring_mask_medium(dilated_mask_M_LD)

                ## Small
                ### High Density
                mask_S_HD_3D = Array{Bool}(undef, size(arr))
                for z in 1:size(arr, 3)
                    mask_S_HD_3D[:, :, z] = mask_S_HD
                end
                dilated_mask_S_HD = dilate_mask_small(mask_S_HD_3D)
                ring_mask_S_HD = ring_mask_small(dilated_mask_S_HD)

                ### Medium Density
                mask_S_MD_3D = Array{Bool}(undef, size(arr))
                for z in 1:size(arr, 3)
                    mask_S_MD_3D[:, :, z] = mask_S_MD
                end
                dilated_mask_S_MD = dilate_mask_small(mask_S_MD_3D)
                ring_mask_S_MD = ring_mask_small(dilated_mask_S_MD)

                ### Low Density
                mask_S_LD_3D = Array{Bool}(undef, size(arr))
                for z in 1:size(arr, 3)
                    mask_S_LD_3D[:, :, z] = mask_S_LD
                end
                dilated_mask_S_LD = dilate_mask_small(mask_S_LD_3D)
                ring_mask_S_LD = ring_mask_small(dilated_mask_S_LD)

                #---------------- Volume Fraction ----------------#

                ## Calibration Prep
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
                mask_cal_3D = zeros(size(c_img))
                for z in 1:size(c_img, 3)
                    mask_cal_3D[:, :, z] = Bool.(erode(bool_arr_erode))
                end

                hu_calcium = mean(c_img[Bool.(mask_cal_3D)])
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
                push!(dfs_vf, df)

                #---------------- Integrated ----------------#

                local intensity_array
                if VENDOR == "80"
                    intensity_array = [23, 69.7, 377.3, 1375.3]
                elseif VENDOR == "100"
                    intensity_array = [23, 63.3, 326.0, 1179.4]
                elseif VENDOR == "120"
                    intensity_array = [20, 58.1, 296.9, 1072.1]
                else
                    intensity_array = [20, 55.5, 282.7, 1019.7]
                end

                # df_cal = DataFrame(:density => [0, 0.025, 0.200, 0.800], :intensity => intensity_array)
                df_cal = DataFrame(:density => [0, 0.2], :intensity => [hu_heart_tissue_large_bkg, hu_calcium])


                linearRegressor = lm(@formula(intensity ~ density), df_cal)
                linearFit = GLM.predict(linearRegressor)
                m = linearRegressor.model.pp.beta0[2]
                b = linearRegressor.model.rr.mu[1]
                density(intensity) = (intensity - b) / m
                intensity(ρ) = m * ρ + b

                # Score background
                S_Obj = intensity(ρ_calcium)
                ρ = ρ_calcium # mg/mm^3

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
                alg_L_HD = Integrated(arr[mask_L_HD_3D])
                mass_l_hd = score(s_bkg_L_HD, S_Obj, pixel_size, ρ, alg_L_HD)

                ## Medium Density
                s_bkg_L_MD = mean(arr[ring_mask_L_MD])
                alg_L_MD = Integrated(arr[mask_L_MD_3D])
                mass_l_md = score(s_bkg_L_MD, S_Obj, pixel_size, ρ, alg_L_MD)

                ## Low Density
                s_bkg_L_LD = mean(arr[ring_mask_L_LD])
                alg_L_LD = Integrated(arr[mask_L_LD_3D])
                mass_l_ld = score(s_bkg_L_LD, S_Obj, pixel_size, ρ, alg_L_LD)

                # Score Medium Inserts
                ## High Density
                s_bkg_M_HD = mean(arr[ring_mask_M_HD])
                alg_M_HD = Integrated(arr[mask_M_HD_3D])
                mass_m_hd = score(s_bkg_M_HD, S_Obj, pixel_size, ρ, alg_M_HD)

                ## Medium Density
                s_bkg_M_MD = mean(arr[ring_mask_M_MD])
                alg_M_MD = Integrated(arr[mask_M_MD_3D])
                mass_m_md = score(s_bkg_M_MD, S_Obj, pixel_size, ρ, alg_M_MD)

                ## Low Density
                s_bkg_M_LD = mean(arr[ring_mask_M_LD])
                alg_M_LD = Integrated(arr[mask_M_LD_3D])
                mass_m_ld = score(s_bkg_M_LD, S_Obj, pixel_size, ρ, alg_M_LD)

                # Score Small Inserts
                ## High Density
                s_bkg_S_HD = mean(arr[ring_mask_S_HD])
                alg_S_HD = Integrated(arr[mask_S_HD_3D])
                mass_s_hd = score(s_bkg_S_HD, S_Obj, pixel_size, ρ, alg_S_HD)

                ## Medium Density
                s_bkg_S_MD = mean(arr[ring_mask_S_MD])
                alg_S_MD = Integrated(arr[mask_S_MD_3D])
                mass_s_md = score(s_bkg_S_MD, S_Obj, pixel_size, ρ, alg_S_MD)

                ## Low Density
                s_bkg_S_LD = mean(arr[ring_mask_S_LD])
                alg_S_LD = Integrated(arr[mask_S_LD_3D])
                mass_s_ld = score(s_bkg_S_LD, S_Obj, pixel_size, ρ, alg_S_LD)

                # Results

                calculated_mass_large = [mass_l_ld, mass_l_md, mass_l_hd]
                calculated_mass_medium = [mass_m_ld, mass_m_md, mass_m_hd]
                calculated_mass_small = [mass_s_ld, mass_s_md, mass_s_hd]

                df = DataFrame(;
                    VENDOR=VENDOR,
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

                #---------------- SWCS ----------------#
                # local μ, σ
                # if VENDOR == "80"
                #     μ, σ = 170, 40
                # elseif VENDOR == "100"
                #     μ, σ = 165, 40
                # elseif VENDOR == "120"
                #     μ, σ = 160, 40
                # else
                #     VENDOR == "135"
                #     μ, σ = 155, 40
                # end
                # # μ, σ = mean(c_img[Bool.(mask_cal_3D)]), std(c_img[Bool.(mask_cal_3D)])
				μ, σ = mean(c_img[Bool.(mask_cal_3D)]) / 2, std(c_img[Bool.(mask_cal_3D)])


                # Mask Calibration Factor
                output = calc_output(
                    masked_array, header, slice_CCI, thresh, trues(3, 3)
                )
                insert_centers = calc_centers(
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
                    SIZE=SIZE,
                    DENSITY=DENSITY,
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
                output = calc_output(
                    masked_array, header, slice_CCI, thresh, trues(3, 3)
                )
                insert_centers = calc_centers(
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
                    SIZE=SIZE,
                    DENSITY=DENSITY,
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
end

# ╔═╡ f4e99c67-20e0-44d8-9bb0-e4ae7c162ad1
md"""
# Save Results
"""

# ╔═╡ 09d11f72-12aa-4c74-9f76-9eb2a29462b3
dfs_vf

# ╔═╡ 76a95026-f7ba-4804-93d2-c97fa4176669
dfs_s

# ╔═╡ 78c039dc-f066-4826-b2eb-50f0ff89b7df
dfs_i

# ╔═╡ aa89bea1-6a57-4267-9142-e552f9519b85
dfs_a

# ╔═╡ 4f1ce2d9-44dd-4567-b146-10795b1dabb9
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

# ╔═╡ de89f630-671d-4aa1-8d25-cd810397531e
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
# ╠═616a3e0e-82e8-4b08-ac82-1652403a693e
# ╠═40c1b1e3-e103-48d5-8e9a-dd408ec162c6
# ╠═e8ce5c83-98b6-4ec4-beac-677fd24ce62d
# ╠═c7a8708a-14b3-418c-9a26-a45f8b84711b
# ╠═992dca9e-3664-4798-a435-72ecd4c797d7
# ╠═1af7db2f-126d-4c5e-8091-7b937744b22c
# ╠═3fa8dba1-e564-410c-9ce5-36870795f493
# ╠═0d2127b5-489a-401e-a7eb-9c2f5a0aa22a
# ╟─1f88f342-7418-4d8a-af66-53a985f10209
# ╟─2e03d03e-c048-48ce-8b5a-449bea68f258
# ╟─e5e63f2b-38de-41bc-9785-576d8f3bb8d3
# ╟─766dcf5a-b778-423a-b39b-610299fcf910
# ╟─39eebc5b-e16a-43ae-abdd-19e44813e9c8
# ╟─872ec5e3-6001-40a1-b5c6-075b3560d34d
# ╟─26fee6aa-cbf1-475b-8553-8dab04e0f80d
# ╟─066702c9-afa8-49e0-81bf-d2bfcdfd2c87
# ╟─485030f6-e3e8-477f-bda3-e1b0b321df54
# ╟─e5609112-e6ce-4fba-9237-c67941fc2fe0
# ╠═73f69fdd-20e7-4aa6-9a00-942ab061433c
# ╟─f4e99c67-20e0-44d8-9bb0-e4ae7c162ad1
# ╠═09d11f72-12aa-4c74-9f76-9eb2a29462b3
# ╠═76a95026-f7ba-4804-93d2-c97fa4176669
# ╠═78c039dc-f066-4826-b2eb-50f0ff89b7df
# ╠═aa89bea1-6a57-4267-9142-e552f9519b85
# ╠═4f1ce2d9-44dd-4567-b146-10795b1dabb9
# ╠═de89f630-671d-4aa1-8d25-cd810397531e

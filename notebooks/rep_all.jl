### A Pluto.jl notebook ###
# v0.19.16

using Markdown
using InteractiveUtils

# ╔═╡ 5a293477-d953-403e-80d2-fb0bf5b9cb3a
# ╠═╡ show_logs = false
begin
    using Pkg
    Pkg.activate(".")

    using PlutoUI, Statistics, CSV, DataFrames, GLM, CairoMakie, HypothesisTests, Colors, MLJBase, DICOM, DICOMUtils, PhantomSegmentation, CalciumScoring, ImageMorphology, ImageFiltering, Noise
    using StatsBase: quantile!, rmsd
end

# ╔═╡ 2ede7b13-103e-439c-bc05-1bfb2292e55e
TableOfContents()

# ╔═╡ aeb1fa4d-affe-49c8-a555-770092794242
IMAGES = "images_reproducibility1"

# ╔═╡ a263b900-72a7-4914-9f71-b068b8ae910c
OUTPUT = "output_repeated"

# ╔═╡ 361b88e6-869c-4c79-985b-6cb40eb66cef
SAVE_DF = "full2.csv"

# ╔═╡ 50a74032-8b52-4803-b130-303a20e0c49c
VENDORS = ["80", "100", "120", "135"]

# ╔═╡ f4d7f4f6-e948-42de-b8d1-405a322dde86
SIZES = ["small", "medium", "large"]

# ╔═╡ df35ac7a-2a17-4fcf-a60c-869b654e53fd
DENSITIES = ["low", "normal"]

# ╔═╡ f314e335-df3a-4115-96f2-e7a8c63aeee6
function create_mask(array, mask)
    @assert size(array) == size(mask)
    idxs = findall(x -> x == true, mask)
    overlayed_mask = zeros(size(array))
    for idx in idxs
        overlayed_mask[idx] = array[idx]
    end
    return overlayed_mask
end

# ╔═╡ 7a5de9b2-25e3-4c70-981f-768b13113093
function dilate_mask_large(mask)
    return dilate(mask)
end

# ╔═╡ c77f25d2-6596-4b51-a586-6c75b8534d8a
function ring_mask_large(dilated_mask)
    return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask)
end

# ╔═╡ 223a3fa6-6244-4f4e-9366-37a0fad1474e
function dilate_mask_medium(mask)
    return (mask)
end

# ╔═╡ 7c0c6116-1675-4d43-befb-dd6a18e3454b
function ring_mask_medium(dilated_mask)
    return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask)
end

# ╔═╡ e1f68840-1eae-4939-9985-68cd2c8b680a
function dilate_mask_small(mask)
    return (mask)
end

# ╔═╡ 72677ba6-2004-4c56-8c66-563155256da3
function ring_mask_small(dilated_mask)
    return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask)
end

# ╔═╡ fc543228-f307-4aa5-b805-da855c62b609
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
                BASE_PATH = joinpath("/Users/daleblack/Google Drive/dev/MolloiLab/cac-simulation", IMAGES, SIZE, DENSITY)
                root_path = joinpath(BASE_PATH, VENDOR)
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
                root_new = joinpath("/Users/daleblack/Google Drive/dev/MolloiLab/cac-simulation/julia_arrays", SIZE)
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
                mask_cal_3D = Array{Bool}(undef, size(c_img))
                for z in 1:size(c_img, 3)
                    mask_cal_3D[:, :, z] = bool_arr_erode
                end

                if VENDOR == "80"
                    hu_calcium = 377.3
                elseif VENDOR == "100"
                    hu_calcium = 326.0
                elseif VENDOR == "120"
                    hu_calcium = 296.9
                else
                    hu_calcium = 282.7
                end
                ρ_calcium = 0.2


                # Background
                background_mask = zeros(size(arr)...)
                background_mask[
                    (center_insert[1]-5):(center_insert[1]+5),
                    (center_insert[2]-5):(center_insert[2]+5),
                    2,
                ] .= 1
                background_mask = Bool.(background_mask)
                background_ring = ring_mask_large(background_mask)
                hu_heart_tissue_bkg = mean(arr[background_ring])
                mass_bkg = score(arr[background_mask], hu_calcium, hu_heart_tissue_bkg, voxel_size, ρ_calcium, VolumeFraction())

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

                df_cal = DataFrame(:density => [0, 0.025, 0.200, 0.800], :intensity => intensity_array)

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
                S_Obj = intensity(ρ_calcium)
                ρ = ρ_calcium # mg/mm^3

                ring_background = Bool.(background_mask_dil - background_mask)
                s_bkg = mean(arr[ring_background])
                alg_bkg = Integrated(arr[background_mask])

                mass_bkg = score(s_bkg, S_Obj, pixel_size, ρ, alg_bkg)

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
                local μ, σ
                if VENDOR == "80"
                    μ, σ = 170, 40
                elseif VENDOR == "100"
                    μ, σ = 165, 40
                elseif VENDOR == "120"
                    μ, σ = 160, 40
                else
                    VENDOR == "135"
                    μ, σ = 155, 40
                end

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

                local agat_thresh
                agat_thresh = 130

                background_mask = zeros(size(arr)...)
                background_mask[
                    (center_insert[1]-5):(center_insert[1]+5),
                    (center_insert[2]-5):(center_insert[2]+5),
                    2,
                ] .= 1
                background_mask = Bool.(background_mask)
                alg = Agatston()
                overlayed_bkg_mask = create_mask(arr, background_mask)
                agat_bkg, mass_bkg = score(
                    overlayed_bkg_mask,
                    pixel_size,
                    mass_cal_factor,
                    alg;
                    kV=kV
                )

                # Score Large Inserts
                ## High Density
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
                overlayed_mask_l_md = create_mask(arr, dilated_mask_L_MD)
                agat_l_md, mass_l_md = score(
                    overlayed_mask_l_md,
                    pixel_size,
                    mass_cal_factor,
                    alg;
                    kV=kV
                )

                ## Low Density
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
                overlayed_mask_m_hd = create_mask(arr, dilated_mask_M_HD)
                agat_m_hd, mass_m_hd = score(
                    overlayed_mask_m_hd,
                    pixel_size,
                    mass_cal_factor,
                    alg;
                    kV=kV
                )

                ## Medium Density
                overlayed_mask_m_md = create_mask(arr, dilated_mask_M_MD)
                agat_m_md, mass_m_md = score(
                    overlayed_mask_m_md,
                    pixel_size,
                    mass_cal_factor,
                    alg;
                    kV=kV
                )

                ## Low Density
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
                overlayed_mask_s_hd = create_mask(arr, dilated_mask_S_HD)
                agat_s_hd, mass_s_hd = score(
                    overlayed_mask_s_hd,
                    pixel_size,
                    mass_cal_factor,
                    alg;
                    kV=kV
                )

                ## Medium Density
                overlayed_mask_s_md = create_mask(arr, dilated_mask_S_MD)
                agat_s_md, mass_s_md = score(
                    overlayed_mask_s_md,
                    pixel_size,
                    mass_cal_factor,
                    alg;
                    kV=kV
                )

                ## Low Density
                overlayed_mask_s_ld = create_mask(arr, dilated_mask_S_LD)
                agat_s_ld, mass_s_ld = score(
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

# ╔═╡ 0f7901f9-cab9-434a-ac18-fe57148f810b
md"""
# Save Results
"""

# ╔═╡ f5b3fbf7-2e82-4055-80ab-3a63c476e4bb
dfs_vf

# ╔═╡ 0d683ee1-9a1a-46df-a12f-d10f19ad1f44
dfs_s

# ╔═╡ a347fb07-a665-4671-ac16-f3e75e5043f1
dfs_i

# ╔═╡ f673d500-2031-4267-9cc5-285058a71b3d
dfs_a

# ╔═╡ a0b662f3-70ef-402c-b0a6-cd6fd2cf24ae
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

# ╔═╡ d6effd32-d81f-49bb-a666-334053e68014
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
# ╠═5a293477-d953-403e-80d2-fb0bf5b9cb3a
# ╠═2ede7b13-103e-439c-bc05-1bfb2292e55e
# ╠═aeb1fa4d-affe-49c8-a555-770092794242
# ╠═a263b900-72a7-4914-9f71-b068b8ae910c
# ╠═361b88e6-869c-4c79-985b-6cb40eb66cef
# ╠═50a74032-8b52-4803-b130-303a20e0c49c
# ╠═f4d7f4f6-e948-42de-b8d1-405a322dde86
# ╠═df35ac7a-2a17-4fcf-a60c-869b654e53fd
# ╟─f314e335-df3a-4115-96f2-e7a8c63aeee6
# ╟─7a5de9b2-25e3-4c70-981f-768b13113093
# ╟─c77f25d2-6596-4b51-a586-6c75b8534d8a
# ╟─223a3fa6-6244-4f4e-9366-37a0fad1474e
# ╟─7c0c6116-1675-4d43-befb-dd6a18e3454b
# ╟─e1f68840-1eae-4939-9985-68cd2c8b680a
# ╟─72677ba6-2004-4c56-8c66-563155256da3
# ╠═fc543228-f307-4aa5-b805-da855c62b609
# ╟─0f7901f9-cab9-434a-ac18-fe57148f810b
# ╠═f5b3fbf7-2e82-4055-80ab-3a63c476e4bb
# ╠═0d683ee1-9a1a-46df-a12f-d10f19ad1f44
# ╠═a347fb07-a665-4671-ac16-f3e75e5043f1
# ╠═f673d500-2031-4267-9cc5-285058a71b3d
# ╠═a0b662f3-70ef-402c-b0a6-cd6fd2cf24ae
# ╠═d6effd32-d81f-49bb-a666-334053e68014

### A Pluto.jl notebook ###
# v0.19.8

using Markdown
using InteractiveUtils

# ╔═╡ 6e3c5404-386a-4bbb-b911-2b314e589a36
# ╠═╡ show_logs = false
begin
    using Pkg
    Pkg.activate(".")

    using PlutoUI, Statistics, CSV, DataFrames, GLM, CairoMakie, HypothesisTests, Colors, MLJBase, DICOM, DICOMUtils, PhantomSegmentation, CalciumScoring, ImageMorphology, ImageFiltering, Noise
    using StatsBase: quantile!, rmsd
end

# ╔═╡ 6319ab98-eae2-45ce-b18c-29fda037fb77
TableOfContents()

# ╔═╡ 76d9f096-194f-4981-9b46-ea78b5a15ce6
function create_mask(array, mask)
    @assert size(array) == size(mask)
    idxs = findall(x -> x == true, mask)
    overlayed_mask = zeros(size(array))
    for idx in idxs
        overlayed_mask[idx] = array[idx]
    end
    return overlayed_mask
end

# ╔═╡ 1909e093-f52e-4d0d-ae99-bbd68f506883
function create_mask(array, mask)
    @assert size(array) == size(mask)
    idxs = findall(x -> x == true, mask)
    overlayed_mask = zeros(size(array))
    for idx in idxs
        overlayed_mask[idx] = array[idx]
    end
    return overlayed_mask
end

# ╔═╡ fb0d0ad3-1009-4f4f-9eb6-6835b7623aac
function dilate_mask_large(mask)
    return dilate(mask)
end

# ╔═╡ 08b8cce6-5cc0-4f5a-8a67-9cabbd01addb
function ring_mask_large(dilated_mask)
    return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask)
end

# ╔═╡ a1d6d9a3-5066-4eb9-a47a-a4b8266cfed3
function dilate_mask_medium(mask)
    return (mask)
end

# ╔═╡ 124d37bb-792c-469a-a6b1-b59b2a047dc8
function ring_mask_medium(dilated_mask)
    return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask)
end

# ╔═╡ ce37dd7d-a89b-4ad7-8385-194f653d7a4f
function dilate_mask_small(mask)
    return (mask)
end

# ╔═╡ 0d085820-dcbe-417d-87af-9a54a90340e9
function ring_mask_small(dilated_mask)
    return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask)
end

# ╔═╡ 5dda487e-0380-4c23-b340-5322ee3bd3ab
TYPE = "swcs"

# ╔═╡ 6f457845-bfdd-4bd9-8382-1279ab266b5a
VENDERS = ["80", "100", "120", "135"]

# ╔═╡ 1cd37003-cb90-44fc-8865-9f17d23e7fd4
SIZES = ["small", "medium", "large"]

# ╔═╡ ced8bb98-2590-4b52-b69b-3ea29fa60b99
DENSITIES = ["low", "normal"]


# ╔═╡ 9e8a0097-6680-4210-9cc7-582272de49c3
begin
    dfs = []
    for VENDER in VENDERS
        for SIZE in SIZES
            for DENSITY in DENSITIES
                SCAN_NUMBER = 1
                BASE_PATH = string(
                    "/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/images_new/",
                    SIZE,
                    "/",
                    DENSITY,
                    "/",
                )
                root_path = string(BASE_PATH, VENDER)
                dcm_path_list = dcm_list_builder(root_path)
                pth = dcm_path_list[SCAN_NUMBER]
                scan = basename(pth)
                header, dcm_array, slice_thick_ori1 = dcm_reader(pth)

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

                # Segment Calcium Inserts
                # mask_L_HD, mask_M_HD, mask_S_HD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_LD, mask_M_LD, mask_S_LD = mask_inserts_simulation(
                # 		dcm_array, masked_array, header, slice_CCI, center_insert
                # )
                root_new = string(
                    "/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/julia_arrays/",
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
                    blur=blur,
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
                push!(dfs, df)
            end
        end
    end
end

# ╔═╡ ee765a05-820e-441e-a000-64fed9346d25
md"""
# Save results
"""

# ╔═╡ d54279d4-494a-4430-9587-449a32e6cece
dfs

# ╔═╡ c7b29675-2392-40f9-b4c8-94afe85a9a9c
if ~isdir(string(cd(pwd, ".."), "/output_new/", TYPE))
    mkdir(string(cd(pwd, ".."), "/output_new/", TYPE))
end

# ╔═╡ b60b6d28-d1a3-46bd-b59a-37d1663452ca
begin
    new_df = vcat(dfs[1:length(dfs)]...)
    output_path_new = string(cd(pwd, ".."), "/output_new/", TYPE, "/", "full2.csv")
    CSV.write(output_path_new, new_df)
end

# ╔═╡ 99d772a2-e178-4f6b-94bc-1987c7a952ec

# ╔═╡ 7a3b633f-3045-42b9-9722-cd9a1b0b7217

# ╔═╡ dc4686c4-e59a-4b74-95ce-8ac66a4ed9cd

# ╔═╡ Cell order:
# ╠═6e3c5404-386a-4bbb-b911-2b314e589a36
# ╠═6319ab98-eae2-45ce-b18c-29fda037fb77
# ╟─76d9f096-194f-4981-9b46-ea78b5a15ce6
# ╠═0d10d0e3-5833-4330-9c6a-a6a620b56d9e
# ╠═af7aaaa6-4333-4157-b77c-dad9e07b3c46
# ╠═dd4d15fc-02af-42ab-a3b7-fdeaf1b1eb09
# ╠═02940e5d-22db-4f71-adba-b4d9807d5974
# ╠═c77f8ff1-ed5f-4d68-a9d1-b0ef304634f1
# ╠═9e8a0097-6680-4210-9cc7-582272de49c3
# ╟─ee765a05-820e-441e-a000-64fed9346d25
# ╠═d54279d4-494a-4430-9587-449a32e6cece
# ╠═c7b29675-2392-40f9-b4c8-94afe85a9a9c
# ╠═b60b6d28-d1a3-46bd-b59a-37d1663452ca
# ╠═99d772a2-e178-4f6b-94bc-1987c7a952ec
# ╠═7a3b633f-3045-42b9-9722-cd9a1b0b7217
# ╠═dc4686c4-e59a-4b74-95ce-8ac66a4ed9cd

### A Pluto.jl notebook ###
# v0.19.16

using Markdown
using InteractiveUtils

# ╔═╡ d8e1422a-f875-4bc8-8943-ec373a89f960
# ╠═╡ show_logs = false
begin
    using Pkg
    Pkg.activate(".")

    using PlutoUI, Statistics, CSV, DataFrames, GLM, CairoMakie, HypothesisTests, Colors, MLJBase, DICOM, DICOMUtils, PhantomSegmentation, CalciumScoring, ImageMorphology, ImageFiltering, Noise
    using StatsBase: quantile!, rmsd
end

# ╔═╡ be04238e-82c7-421e-b976-dc9420316921
TableOfContents()

# ╔═╡ ab1b36c0-e7aa-4b61-a2f0-1b7b9993b53e
function create_mask(array, mask)
    @assert size(array) == size(mask)
    idxs = findall(x -> x == true, mask)
    overlayed_mask = zeros(size(array))
    for idx in idxs
        overlayed_mask[idx] = array[idx]
    end
    return overlayed_mask
end

# ╔═╡ aded68a2-529e-47f3-8266-bfb5f2131f8f
function dilate_mask_medium(mask)
    return (mask)
end

# ╔═╡ 7f260aaa-8d45-4bd0-a654-19890878af28
function ring_mask_large(dilated_mask)
    return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask)
end

# ╔═╡ 577d6155-f097-4bff-8fb9-6a4371749a99
VENDERS = ["80", "100", "120", "135"]

# ╔═╡ c1be2978-a336-496b-b6bf-1474d50c4f59
function dilate_mask_large(mask)
    return dilate(mask)
end

# ╔═╡ 475ceb61-6133-4085-9421-52cd8cbc5671
function dilate_mask_small(mask)
    return (mask)
end

# ╔═╡ 8cbc5578-a465-4447-b366-504c6bf7e851
function ring_mask_medium(dilated_mask)
    return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask)
end

# ╔═╡ 8619ec6c-d6cc-40c5-9859-96c47881091d
function ring_mask_small(dilated_mask)
    return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask)
end

# ╔═╡ 6e9c7811-f168-4125-8661-e5be7c1022a5
SIZES = ["small", "medium", "large"]

# ╔═╡ c5c2c90b-6574-4b54-9e37-c32fee6f575f
DENSITIES = ["low", "normal"]

# ╔═╡ 48b0eced-4d4d-4c69-a188-ed6235ad42c7
TYPE = "swcs"

# ╔═╡ bb82eaa0-228e-42e5-8522-538d4cc96dfe
begin
    dfs = []
    for VENDER in VENDERS
        for SIZE in SIZES
            for DENSITY in DENSITIES
                SCAN_NUMBER = 1
                BASE_PATH = string(
                    "/Users/daleblack/Google Drive/dev/MolloiLab/cac-simulation/images_new/",
                    SIZE,
                    "/",
                    DENSITY,
                    "/",
                )
                root_path = string(BASE_PATH, string(VENDER * "-motion"))
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

                # Mask Calibration Factor
				global output
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

# ╔═╡ 1037f865-82ee-4227-bfc7-abdedf5a2518
output

# ╔═╡ 7a60d33c-c520-4c0c-b192-59435fab798f
md"""
# Save results
"""

# ╔═╡ e3c9c04c-ec92-4719-a9a5-adf8f87b504c
dfs

# ╔═╡ 73942ca6-429e-4a63-9a50-7534f72941bf
# if ~isdir(string(cd(pwd, ".."), "/output_new/", TYPE))
#     mkdir(string(cd(pwd, ".."), "/output_new/", TYPE))
# end

# ╔═╡ 1383f020-5f64-4239-a6a1-487ef6ee01f9
# begin
#     new_df = vcat(dfs[1:length(dfs)]...)
#     output_path_new = string(cd(pwd, ".."), "/output_new/", TYPE, "/", "full2.csv")
#     CSV.write(output_path_new, new_df)
# end

# ╔═╡ Cell order:
# ╠═d8e1422a-f875-4bc8-8943-ec373a89f960
# ╠═be04238e-82c7-421e-b976-dc9420316921
# ╠═ab1b36c0-e7aa-4b61-a2f0-1b7b9993b53e
# ╠═aded68a2-529e-47f3-8266-bfb5f2131f8f
# ╠═7f260aaa-8d45-4bd0-a654-19890878af28
# ╠═577d6155-f097-4bff-8fb9-6a4371749a99
# ╠═c1be2978-a336-496b-b6bf-1474d50c4f59
# ╠═475ceb61-6133-4085-9421-52cd8cbc5671
# ╠═8cbc5578-a465-4447-b366-504c6bf7e851
# ╠═8619ec6c-d6cc-40c5-9859-96c47881091d
# ╠═6e9c7811-f168-4125-8661-e5be7c1022a5
# ╠═c5c2c90b-6574-4b54-9e37-c32fee6f575f
# ╠═48b0eced-4d4d-4c69-a188-ed6235ad42c7
# ╠═bb82eaa0-228e-42e5-8522-538d4cc96dfe
# ╠═1037f865-82ee-4227-bfc7-abdedf5a2518
# ╟─7a60d33c-c520-4c0c-b192-59435fab798f
# ╠═e3c9c04c-ec92-4719-a9a5-adf8f87b504c
# ╠═73942ca6-429e-4a63-9a50-7534f72941bf
# ╠═1383f020-5f64-4239-a6a1-487ef6ee01f9

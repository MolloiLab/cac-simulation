### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ 17f7c8f3-b2bc-4c41-bec6-2e662b85dadf
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
        Pkg.add(; url="https://github.com/JuliaHealth/DICOM.jl")
        Pkg.add(; url="https://github.com/Dale-Black/DICOMUtils.jl")
        Pkg.add(; url="https://github.com/Dale-Black/PhantomSegmentation.jl")
        Pkg.add(; url="https://github.com/Dale-Black/CalciumScoring.jl")
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

# ╔═╡ 523393ee-0b74-4ece-84ef-ae6b13b70139
TableOfContents()

# ╔═╡ 67f2e345-171c-4ed4-8a1c-e1e0be3e33c5
function create_mask(array, mask)
    @assert size(array) == size(mask)
    idxs = findall(x -> x == true, mask)
    overlayed_mask = zeros(size(array))
    for idx in idxs
        overlayed_mask[idx] = array[idx]
    end
    return overlayed_mask
end

# ╔═╡ 2304c588-078f-48c3-9afc-ceceebb06f33
TYPE = "swcs"

# ╔═╡ 20a9edb7-5fa3-4175-861f-16bb15a1d030
VENDERS = ["80", "100", "120", "135"]

# ╔═╡ 901cf173-f94a-4fde-9014-f001e04138aa
SIZES = ["small", "medium", "large"]

# ╔═╡ 3bcf713b-61d7-4ebd-a282-d49c244e6723
DENSITIES = ["low", "normal"]

# ╔═╡ 1152699e-eae8-4d09-ba6c-d74eaaaf9ae1
blurs = [0, 0.5, 1, 1.5, 2]

# ╔═╡ 72bbc18a-3e45-479e-856b-4034de97ca12
begin
    dfs = []
    for blur in blurs
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

                    # Motion Blur
                    if blur != 0
                        for z in size(dcm_array, 3)
                            dcm_array[:, :, z] = mult_gauss(dcm_array[:, :, z], blur)
                        end
                    end

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

                    # # Segment Calcium Rod
                    # local thresh
                    # if DENSITY == "low" && SIZE == "large"
                    # 	thresh = 75
                    # elseif DENSITY == "low" && SIZE == "medium"
                    # 	thresh = 75
                    # elseif DENSITY == "low"
                    # 	thresh = 60
                    # elseif DENSITY ==  "normal"
                    # 	thresh = 130
                    # end

                    # # Segment Calcium Rod (reproducibility1)
                    # local thresh
                    # if DENSITY == "low" && SIZE == "large" && VENDER == "80"
                    # 	thresh = 80
                    # elseif DENSITY == "low" && SIZE == "large" && VENDER == "100"
                    # 	thresh = 70
                    # elseif DENSITY == "low" && SIZE == "large"
                    # 	thresh = 75
                    # elseif DENSITY == "low" && SIZE == "medium" && VENDER == "135"
                    # 	thresh = 55
                    # elseif DENSITY == "low" && SIZE == "medium"
                    # 	thresh = 75
                    # elseif DENSITY == "low"
                    # 	thresh = 60
                    # elseif DENSITY ==  "normal"
                    # 	thresh = 130
                    # end

                    @info blur, DENSITY, SIZE, VENDER
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
                        (center_insert[1] - 5):(center_insert[1] + 5),
                        (center_insert[2] - 5):(center_insert[2] + 5),
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
                        swcs_bkg=swcs_bkg,
                    )
                    push!(dfs, df)
                end
            end
        end
    end
end

# ╔═╡ a04b25ce-8322-45e8-8df3-ad1cac7d7f3e
md"""
# Save results
"""

# ╔═╡ 78ed7a73-d325-4dd4-bf1c-340300c1288a
dfs

# ╔═╡ 0d6776a5-9a21-4ce7-a61b-d47b97fdd6d8
if ~isdir(string(cd(pwd, ".."), "/output_repeated/", TYPE))
    mkdir(string(cd(pwd, ".."), "/output_repeated/", TYPE))
end

# ╔═╡ a3ab5f45-c061-4c15-8c5c-4d689dcdaa0a
begin
    new_df = vcat(dfs[1:length(dfs)]...)
    output_path_new = string(cd(pwd, ".."), "/output_repeated/", TYPE, "/", "full2.csv")
    CSV.write(output_path_new, new_df)
end

# ╔═╡ Cell order:
# ╠═17f7c8f3-b2bc-4c41-bec6-2e662b85dadf
# ╠═523393ee-0b74-4ece-84ef-ae6b13b70139
# ╟─67f2e345-171c-4ed4-8a1c-e1e0be3e33c5
# ╠═2304c588-078f-48c3-9afc-ceceebb06f33
# ╠═20a9edb7-5fa3-4175-861f-16bb15a1d030
# ╠═901cf173-f94a-4fde-9014-f001e04138aa
# ╠═3bcf713b-61d7-4ebd-a282-d49c244e6723
# ╠═1152699e-eae8-4d09-ba6c-d74eaaaf9ae1
# ╠═72bbc18a-3e45-479e-856b-4034de97ca12
# ╟─a04b25ce-8322-45e8-8df3-ad1cac7d7f3e
# ╠═78ed7a73-d325-4dd4-bf1c-340300c1288a
# ╠═0d6776a5-9a21-4ce7-a61b-d47b97fdd6d8
# ╠═a3ab5f45-c061-4c15-8c5c-4d689dcdaa0a

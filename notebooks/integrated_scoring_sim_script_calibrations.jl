### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ 0bba17a2-406d-4318-a417-30c090605b0d
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
        Pkg.add(; url="https://github.com/JuliaHealth/DICOM.jl")
        Pkg.add(; url="https://github.com/Dale-Black/DICOMUtils.jl")
        Pkg.add(; url="https://github.com/Dale-Black/PhantomSegmentation.jl")
        Pkg.add(; url="https://github.com/Dale-Black/CalciumScoring.jl")
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

# ╔═╡ 90e0891a-e9ed-4e7d-a410-4821bae07371
TableOfContents()

# ╔═╡ 3cebe9fa-3844-4567-89f3-807e956f2bc9
function create_mask(array, mask)
    @assert size(array) == size(mask)
    idxs = findall(x -> x == true, mask)
    overlayed_mask = zeros(size(array))
    for idx in idxs
        overlayed_mask[idx] = array[idx]
    end
    return overlayed_mask
end

# ╔═╡ fc679aac-09a1-48b6-aec8-29dc1642be9f
TYPE = "integrated_scoring"

# ╔═╡ 9ad8b7be-e72e-4ad4-98e8-213ed9734879
VENDERS = ["80", "100", "120", "135"]

# ╔═╡ 1952ad11-fdf7-44d4-9e79-de47021f7b77
SIZES = ["small", "medium", "large"]

# ╔═╡ a1af316d-8b6a-4092-9d67-1972494958f0
DENSITIES = ["low", "normal"]

# ╔═╡ feda751f-2c81-48ff-b414-643940dad296
cals = ["1pt", "3pt", "6pt", "specific"]

# ╔═╡ a0ad9779-9332-4a57-b325-cacf1af485d3
begin
    dfs = []
    high_dens = []
    med_dens = []
    low_dens = []
    high_dens100 = []
    med_dens50 = []
    low_dens25 = []
    for cal in cals
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

                    @info cal, DENSITY, SIZE, VENDER, thresh

                    calcium_image, slice_CCI, quality_slice, cal_rod_slice = mask_rod(
                        masked_array, header; calcium_threshold=thresh
                    )

                    # mask_L_HD, mask_M_HD, mask_S_HD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_LD, mask_M_LD, mask_S_LD = mask_inserts_simulation(
                    # 	dcm_array, masked_array, header, slice_CCI, center_insert)

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

                    # # Calibration Prep
                    # local density_array
                    # if DENSITY == "low"
                    # 	density_array = [0, 25, 50, 100]
                    # elseif DENSITY == "normal"
                    # 	density_array = [0, 200, 400, 800]
                    # end
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

                    local density_array
                    if cal == "1pt"
                        if VENDER == "80"
                            intensity_array = [0, 377.3]
                            density_array = [0, 200]
                        elseif VENDER == "100"
                            intensity_array = [0, 326.0]
                            density_array = [0, 200]
                        elseif VENDER == "120"
                            intensity_array = [0, 296.9]
                            density_array = [0, 200]
                        else
                            intensity_array = [0, 282.7]
                            density_array = [0, 200]
                        end
                    elseif cal == "3pt"
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
                    elseif cal == "6pt"
                        if VENDER == "80"
                            intensity_array = [0, 69.7, 114, 202, 377.3, 717, 1375.3]
                            density_array = [0, 25, 50, 100, 200, 400, 800]
                        elseif VENDER == "100"
                            intensity_array = [0, 63.3, 100, 176, 326.0, 617, 1179.4]
                            density_array = [0, 25, 50, 100, 200, 400, 800]
                        elseif VENDER == "120"
                            intensity_array = [0, 58.1, 92, 161, 296.9, 562, 1072.1]
                            density_array = [0, 25, 50, 100, 200, 400, 800]
                        else
                            intensity_array = [0, 55.5, 88, 153, 282.7, 534, 1019.7]
                            density_array = [0, 25, 50, 100, 200, 400, 800]
                        end
                    else
                        if DENSITY == "low"
                            intensity_array = [
                                0, low_density_cal, med_density_cal, high_density_cal
                            ]
                            density_array = [0, 25, 50, 100]
                        else
                            intensity_array = [
                                0, low_density_cal, med_density_cal, high_density_cal
                            ]
                            density_array = [0, 200, 400, 800]
                        end
                    end

                    # local df_cal
                    # if DENSITY == "low"
                    # 	df_cal = DataFrame(:density => density_array[2:end], :intensity => intensity_array[2:end])
                    # elseif DENSITY == "normal"
                    # 	df_cal = DataFrame(:density => density_array, :intensity => intensity_array)
                    # end

                    df_cal = DataFrame(
                        :density => density_array, :intensity => intensity_array
                    )

                    linearRegressor = lm(@formula(intensity ~ density), df_cal)
                    linearFit = predict(linearRegressor)
                    m = linearRegressor.model.pp.beta0[2]
                    b = linearRegressor.model.rr.mu[1]
                    density(intensity) = (intensity - b) / m
                    intensity(ρ) = m * ρ + b

                    # Score background
                    background_mask = zeros(size(arr)...)
                    background_mask[
                        (center_insert[1] - 5):(center_insert[1] + 5),
                        (center_insert[2] - 5):(center_insert[2] + 5),
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
                        cal=cal,
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
                        mass_bkg=mass_bkg,
                    )
                    push!(dfs, df)

                    # if DENSITY == "normal"
                    # 	# Calculate mean intensities for calibration
                    # 	erode_mask_L_HD = erode(erode(mask_L_HD))
                    # 	mean_HD = mean(single_arr[erode_mask_L_HD])
                    # 	push!(high_dens, mean_HD)

                    # 	erode_mask_L_MD = erode(erode(mask_L_MD))
                    # 	mean_MD = mean(single_arr[erode_mask_L_MD])
                    # 	push!(med_dens, mean_MD)

                    # 	erode_mask_L_LD = erode(erode(mask_L_LD))
                    # 	mean_LD = mean(single_arr[erode_mask_L_LD])
                    # 	push!(low_dens, mean_LD)
                    # else
                    # 	erode_mask_L_HD = erode(erode(mask_L_HD))
                    # 	mean_HD = mean(single_arr[erode_mask_L_HD])
                    # 	push!(high_dens100, mean_HD)

                    # 	erode_mask_L_MD = erode(erode(mask_L_MD))
                    # 	mean_MD = mean(single_arr[erode_mask_L_MD])
                    # 	push!(med_dens50, mean_MD)

                    # 	erode_mask_L_LD = erode(erode(mask_L_LD))
                    # 	mean_LD = mean(single_arr[erode_mask_L_LD])
                    # 	push!(low_dens25, mean_LD)
                    # end
                end
            end
        end
    end
end

# ╔═╡ cba1f065-cc60-4d2f-b770-afbec22ae251
# mean(high_dens), mean(med_dens), mean(low_dens), mean(high_dens100), mean(med_dens50), mean(low_dens25)

# ╔═╡ f2f8c290-29d1-4c26-aa23-a0ae4d236a12
md"""
# Save Results
"""

# ╔═╡ ec2383b0-596d-4d27-984b-78ab792cb27c
dfs

# ╔═╡ 57f81937-2af0-4fe8-aeda-da294c6f9acd
if ~isdir(string(cd(pwd, ".."), "/output_new/", TYPE))
    mkdir(string(cd(pwd, ".."), "/output_new/", TYPE))
end

# ╔═╡ d0db09cc-3f85-4b69-934d-d9ddda916d6d
begin
    new_df = vcat(dfs[1:length(dfs)]...)
    output_path_new = string(
        cd(pwd, ".."), "/output_new/", TYPE, "/", "full_calibrations.csv"
    )
    CSV.write(output_path_new, new_df)
end

# ╔═╡ Cell order:
# ╠═0bba17a2-406d-4318-a417-30c090605b0d
# ╠═90e0891a-e9ed-4e7d-a410-4821bae07371
# ╟─3cebe9fa-3844-4567-89f3-807e956f2bc9
# ╠═fc679aac-09a1-48b6-aec8-29dc1642be9f
# ╠═9ad8b7be-e72e-4ad4-98e8-213ed9734879
# ╠═1952ad11-fdf7-44d4-9e79-de47021f7b77
# ╠═a1af316d-8b6a-4092-9d67-1972494958f0
# ╠═feda751f-2c81-48ff-b414-643940dad296
# ╠═a0ad9779-9332-4a57-b325-cacf1af485d3
# ╠═cba1f065-cc60-4d2f-b770-afbec22ae251
# ╟─f2f8c290-29d1-4c26-aa23-a0ae4d236a12
# ╠═ec2383b0-596d-4d27-984b-78ab792cb27c
# ╠═57f81937-2af0-4fe8-aeda-da294c6f9acd
# ╠═d0db09cc-3f85-4b69-934d-d9ddda916d6d

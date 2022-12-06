### A Pluto.jl notebook ###
# v0.19.16

using Markdown
using InteractiveUtils

# ╔═╡ 10334dff-74c1-47cc-9eb8-389438e6eccb
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.activate(".")

    using PlutoUI, Statistics, CSV, DataFrames, GLM, CairoMakie, HypothesisTests, Colors, MLJBase, DICOM, DICOMUtils, PhantomSegmentation, CalciumScoring, ImageMorphology, ImageFiltering, Noise
    using StatsBase: quantile!, rmsd
end

# ╔═╡ 2c51cad4-845b-4da4-8ee2-fd69a186fc32
TableOfContents()

# ╔═╡ f6c04411-8086-420a-98a8-31dad838d7e8
function create_mask(array, mask)
    @assert size(array) == size(mask)
    idxs = findall(x -> x == true, mask)
    overlayed_mask = zeros(size(array))
    for idx in idxs
        overlayed_mask[idx] = array[idx]
    end
    return overlayed_mask
end

# ╔═╡ 55c4156c-ceef-48e0-ab9e-c744fdcb5a2d
function dilate_mask_large(mask)
	return dilate(mask)
end

# ╔═╡ 9e4bfdcd-2899-4eb9-9824-5929a986c0d2
function ring_mask_large(dilated_mask)
	return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask);
end

# ╔═╡ 7019b86e-7bd2-4fa5-95e9-bf883b802940
function dilate_mask_medium(mask)
	return (mask)
end

# ╔═╡ 2484db7c-4a69-4ede-bee0-aa8d798e2717
function ring_mask_medium(dilated_mask)
	return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask);
end

# ╔═╡ 871a62c7-9ee1-42bc-9cd5-5302b8cdd932
function dilate_mask_small(mask)
	return (mask)
end

# ╔═╡ d61bcd0c-ba3e-4386-a3eb-e84e1b64cc4d
function ring_mask_small(dilated_mask)
	return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask);
end

# ╔═╡ 06ca1788-fd1c-43c3-8c5b-0dde8ff6c152
TYPE = "volume_fraction"

# ╔═╡ 37c85960-4254-463e-b5cf-30395307936a
VENDERS = ["80", "100", "120", "135"]

# ╔═╡ 7b04f6ab-595f-4df8-9ca3-1539d24b2df2
SIZES = ["small", "medium", "large"]

# ╔═╡ 4ac607b3-d0c2-4186-b572-95aa255d52df
DENSITIES = ["low", "normal"]

# ╔═╡ e8f6f8c3-a066-4122-867e-0ef55a53baf4
blurs = [0, 0.5, 1, 1.5, 2]

# ╔═╡ e963e6e8-bf77-4836-ba4a-0c6b2bcf2e3e
begin
    dfs = []
    high_dens = []
    med_dens = []
    low_dens = []
    high_dens100 = []
    med_dens50 = []
    low_dens25 = []
    for blur in blurs
        for VENDER in VENDERS
            for SIZE in SIZES
                for DENSITY in DENSITIES
                    SCAN_NUMBER = 1
                    BASE_PATH = string(
                        "/Users/daleblack/Google Drive/dev/MolloiLab/cac-simulation/images_reproducibility1/",
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

                    @info blur, DENSITY, SIZE, VENDER, thresh

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
                    )
                    push!(dfs, df)

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
                end
            end
        end
    end
end

# ╔═╡ 2ed22766-a87e-48b1-8a4e-63e5a954d4f0
md"""
# Save Results
"""

# ╔═╡ 434db292-e5a2-4ce4-b1dd-43e2d6b1f513
dfs

# ╔═╡ 896d5af4-ebbb-469a-baff-bc78b49fb0cc
if ~isdir(string(cd(pwd, ".."), "/output_repeated/", TYPE))
    mkpath(string(cd(pwd, ".."), "/output_repeated/", TYPE))
end

# ╔═╡ 24246499-8340-4dac-8473-dca542992260
begin
    new_df = vcat(dfs[1:length(dfs)]...)
    output_path_new = string(cd(pwd, ".."), "/output_repeated/", TYPE, "/", "full2.csv")
    CSV.write(output_path_new, new_df)
end

# ╔═╡ Cell order:
# ╠═10334dff-74c1-47cc-9eb8-389438e6eccb
# ╠═2c51cad4-845b-4da4-8ee2-fd69a186fc32
# ╟─f6c04411-8086-420a-98a8-31dad838d7e8
# ╟─55c4156c-ceef-48e0-ab9e-c744fdcb5a2d
# ╟─9e4bfdcd-2899-4eb9-9824-5929a986c0d2
# ╟─7019b86e-7bd2-4fa5-95e9-bf883b802940
# ╟─2484db7c-4a69-4ede-bee0-aa8d798e2717
# ╟─871a62c7-9ee1-42bc-9cd5-5302b8cdd932
# ╟─d61bcd0c-ba3e-4386-a3eb-e84e1b64cc4d
# ╠═06ca1788-fd1c-43c3-8c5b-0dde8ff6c152
# ╠═37c85960-4254-463e-b5cf-30395307936a
# ╠═7b04f6ab-595f-4df8-9ca3-1539d24b2df2
# ╠═4ac607b3-d0c2-4186-b572-95aa255d52df
# ╠═e8f6f8c3-a066-4122-867e-0ef55a53baf4
# ╠═e963e6e8-bf77-4836-ba4a-0c6b2bcf2e3e
# ╠═2ed22766-a87e-48b1-8a4e-63e5a954d4f0
# ╠═434db292-e5a2-4ce4-b1dd-43e2d6b1f513
# ╠═896d5af4-ebbb-469a-baff-bc78b49fb0cc
# ╠═24246499-8340-4dac-8473-dca542992260

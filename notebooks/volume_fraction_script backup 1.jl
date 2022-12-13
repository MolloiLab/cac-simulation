### A Pluto.jl notebook ###
# v0.19.16

using Markdown
using InteractiveUtils

# ╔═╡ ad9ad67a-bce4-4170-aae6-e28189d56080
# ╠═╡ show_logs = false
begin
    using Pkg
    Pkg.activate(".")

    using PlutoUI, Statistics, CSV, DataFrames, GLM, CairoMakie, HypothesisTests, Colors, MLJBase, DICOM, DICOMUtils, PhantomSegmentation, CalciumScoring, ImageMorphology, ImageFiltering, Noise
    using StatsBase: quantile!, rmsd
end

# ╔═╡ 87d5f049-e9b2-4e73-b9d4-132bdaec8bbf
TableOfContents()

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
TYPE = "volume_fraction"

# ╔═╡ 6f457845-bfdd-4bd9-8382-1279ab266b5a
VENDERS = ["80", "100", "120", "135"]

# ╔═╡ 1cd37003-cb90-44fc-8865-9f17d23e7fd4
SIZES = ["small", "medium", "large"]

# ╔═╡ ced8bb98-2590-4b52-b69b-3ea29fa60b99
DENSITIES = ["low", "normal"]

# ╔═╡ a6050593-6070-415f-8a6b-46c3776be999
begin
    dfs = []
    high_dens = []
    med_dens = []
    low_dens = []
    high_dens100 = []
    med_dens50 = []
    low_dens25 = []
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

# ╔═╡ 312f1b9d-c2c2-48a4-b6d3-19cbb8e800af
md"""
# Save Results
"""

# ╔═╡ cef74cb4-596e-4222-8a1e-ee4faabbfc60
dfs

# ╔═╡ 677102fb-39e8-4ba4-9211-17849277e52b
if ~isdir(string(cd(pwd, ".."), "/output_new/", TYPE))
    mkpath(string(cd(pwd, ".."), "/output_new/", TYPE))
end

# ╔═╡ c504b84a-2c28-4a52-b74d-67314d0f6c9c
begin
    new_df = vcat(dfs[1:length(dfs)]...)
    output_path_new = string(cd(pwd, ".."), "/output_new/", TYPE, "/", "full2.csv")
    CSV.write(output_path_new, new_df)
end

# ╔═╡ Cell order:
# ╠═ad9ad67a-bce4-4170-aae6-e28189d56080
# ╠═87d5f049-e9b2-4e73-b9d4-132bdaec8bbf
# ╟─1909e093-f52e-4d0d-ae99-bbd68f506883
# ╠═fb0d0ad3-1009-4f4f-9eb6-6835b7623aac
# ╠═08b8cce6-5cc0-4f5a-8a67-9cabbd01addb
# ╠═a1d6d9a3-5066-4eb9-a47a-a4b8266cfed3
# ╠═124d37bb-792c-469a-a6b1-b59b2a047dc8
# ╠═ce37dd7d-a89b-4ad7-8385-194f653d7a4f
# ╠═0d085820-dcbe-417d-87af-9a54a90340e9
# ╠═5dda487e-0380-4c23-b340-5322ee3bd3ab
# ╠═6f457845-bfdd-4bd9-8382-1279ab266b5a
# ╠═1cd37003-cb90-44fc-8865-9f17d23e7fd4
# ╠═ced8bb98-2590-4b52-b69b-3ea29fa60b99
# ╠═61f0444a-e1a6-46f4-8caf-0b3759f64077
# ╠═a6050593-6070-415f-8a6b-46c3776be999
# ╟─312f1b9d-c2c2-48a4-b6d3-19cbb8e800af
# ╠═cef74cb4-596e-4222-8a1e-ee4faabbfc60
# ╠═677102fb-39e8-4ba4-9211-17849277e52b
# ╠═c504b84a-2c28-4a52-b74d-67314d0f6c9c

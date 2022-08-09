### A Pluto.jl notebook ###
# v0.19.4

using Markdown
using InteractiveUtils

# ╔═╡ 4bf9f49b-f23f-45ea-83ce-2f71a8a47e5d
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

# ╔═╡ 32c6e7e9-21a4-4b15-9b17-57851d599ded
TableOfContents()

# ╔═╡ 356e5a8c-9bb5-487e-b4fb-ab9a54e44e30
function create_mask(array, mask)
    @assert size(array) == size(mask)
    idxs = findall(x -> x == true, mask)
    overlayed_mask = zeros(size(array))
    for idx in idxs
        overlayed_mask[idx] = array[idx]
    end
    return overlayed_mask
end

# ╔═╡ cbe7d757-d1e0-407d-9ed6-c9b929db4e8d
TYPE = "integrated_scoring"

# ╔═╡ 79d4812d-f3be-44e8-914f-54de012f9fe7
VENDERS = ["80", "100", "120", "135"]

# ╔═╡ 88c3c336-7533-4319-ad3f-e7ae2616c112
SIZES = ["small", "medium", "large"]

# ╔═╡ 26704a7d-f695-4ce2-892b-20f606462b68
DENSITIES = ["low", "normal"]

# ╔═╡ d24d719b-2805-41d8-ba50-b0fa237d55b2
begin
    dfs = []
    for VENDER in VENDERS
        for SIZE in SIZES
            for DENSITY in DENSITIES
                SCAN_NUMBER = 1
                BASE_PATH = string(
                    "/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/images/",
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
                if DENSITY == "low" &&
                    SIZE == "large" &&
                    (VENDER == "120" || VENDER == "135")
                    thresh = 75
                elseif DENSITY == "low" && SIZE == "medium"
                    thresh = 75
                elseif DENSITY == "low"
                    thresh = 60
                elseif DENSITY == "normal"
                    thresh = 130
                end

                @info DENSITY, SIZE, VENDER

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
                single_arr = masked_array[:, :, slice_CCI]
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

                df = DataFrame(;
                    VENDER=VENDER,
                    SIZE=SIZE,
                    DENSITY=DENSITY,
                    high_density_cal=high_density_cal,
                    med_density_cal=med_density_cal,
                    low_density_cal=low_density_cal,
                )
                push!(dfs, df)
            end
        end
    end
end

# ╔═╡ f8d2544d-5b4d-4ece-954a-d4f883da3619
dfs

# ╔═╡ de01bd04-0abb-4d90-bbec-334df299b94a
if length(dfs) == 24
    global new_df = vcat(dfs[1:24]...)
end

# ╔═╡ f49a9b42-822c-4e28-a32f-9f765936cf30
df_low, df_normal = groupby(new_df, :DENSITY);

# ╔═╡ 8db8f72c-1ced-4cca-a482-7ece3cf154bf
df_normal_small, df_normal_medium, df_normal_large = groupby(df_normal, :SIZE);

# ╔═╡ 66b0ec94-4da4-4c77-8fea-071d1e164b8d
df_low_small, df_low_medium, df_low_large = groupby(df_low, :SIZE);

# ╔═╡ ab1fae31-cc88-4633-8452-fbed2b4d792c
md"""
## Low Density, Small, 80kV
"""

# ╔═╡ c1b081db-5fe7-4dda-bf92-4c71d0269f6d
df_low_small80, df_low_small100, df_low_small120, df_low_small130 = groupby(
    df_low_small, :VENDER
)

# ╔═╡ 94348b42-9c5d-4a44-aa4f-cd285a0d560e
df_low_small80

# ╔═╡ Cell order:
# ╠═4bf9f49b-f23f-45ea-83ce-2f71a8a47e5d
# ╠═32c6e7e9-21a4-4b15-9b17-57851d599ded
# ╠═356e5a8c-9bb5-487e-b4fb-ab9a54e44e30
# ╠═cbe7d757-d1e0-407d-9ed6-c9b929db4e8d
# ╠═79d4812d-f3be-44e8-914f-54de012f9fe7
# ╠═88c3c336-7533-4319-ad3f-e7ae2616c112
# ╠═26704a7d-f695-4ce2-892b-20f606462b68
# ╠═d24d719b-2805-41d8-ba50-b0fa237d55b2
# ╠═f8d2544d-5b4d-4ece-954a-d4f883da3619
# ╠═de01bd04-0abb-4d90-bbec-334df299b94a
# ╠═f49a9b42-822c-4e28-a32f-9f765936cf30
# ╠═8db8f72c-1ced-4cca-a482-7ece3cf154bf
# ╠═66b0ec94-4da4-4c77-8fea-071d1e164b8d
# ╟─ab1fae31-cc88-4633-8452-fbed2b4d792c
# ╠═c1b081db-5fe7-4dda-bf92-4c71d0269f6d
# ╠═94348b42-9c5d-4a44-aa4f-cd285a0d560e

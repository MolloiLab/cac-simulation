### A Pluto.jl notebook ###
# v0.19.4

using Markdown
using InteractiveUtils

# ╔═╡ eb0d0540-d2f6-11ec-290a-d728e7dffe5c
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
        Pkg.add("Distributions")
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
    using Distributions
    using DICOM
    using DICOMUtils
    using PhantomSegmentation
    using CalciumScoring
end

# ╔═╡ dea06aa6-20e3-4099-8d06-060e7970a535
TableOfContents()

# ╔═╡ b76f95ef-cd93-4e50-994d-ad3fa8433441
function create_mask(array, mask)
    @assert size(array) == size(mask)
    idxs = findall(x -> x == true, mask)
    overlayed_mask = zeros(size(array))
    for idx in idxs
        overlayed_mask[idx] = array[idx]
    end
    return overlayed_mask
end

# ╔═╡ 36b5fe08-1ac1-4734-ba7e-64157c2be717
TYPE = "swcs"

# ╔═╡ a859a9dc-8cb0-4ade-b340-83edfc471927
VENDERS = ["80", "100", "120", "135"]

# ╔═╡ 57e199ba-86a7-4ae2-8ee0-284fdce05dea
SIZES = ["small", "medium", "large"]

# ╔═╡ 75800753-24ca-4fb2-8016-6ae43efad58d
DENSITIES = ["low", "normal"]

# ╔═╡ f0b6724e-9c9f-4998-bb09-5c66345640a4
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
                if DENSITY == "low" && SIZE == "small"
                    thresh = 60
                elseif DENSITY == "low" &&
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

                array_filtered = abs.(mapwindow(median, calcium_image[:, :, 2], (3, 3)))
                bool_arr = array_filtered .> 0
                bool_arr_erode = ((erode(erode(erode(bool_arr)))))

                c_img = calcium_image[:, :, 1:3]
                mask_cal_3D = Array{Bool}(undef, size(c_img))
                for z in 1:size(c_img, 3)
                    mask_cal_3D[:, :, z] = bool_arr_erode
                end
                calibration = c_img[mask_cal_3D]

                df = DataFrame(;
                    VENDER=VENDER,
                    SIZE=SIZE,
                    DENSITY=DENSITY,
                    scan=scan,
                    mean=mean(calibration),
                    std=std(calibration),
                )
                push!(dfs, df)
            end
        end
    end
end

# ╔═╡ 8d20427c-d5a5-46c2-ba78-97138b966776
md"""
## Load DataFrame
"""

# ╔═╡ 45e74fa4-f0c5-457a-8a4b-b433a381b930
if length(dfs) == 24
    new_df = vcat(dfs[1:24]...)
end

# ╔═╡ c190bc16-2a17-4466-9335-0371e7213326
df_80, df_100, df_120, df_135 = groupby(new_df, :scan)

# ╔═╡ 96c5bc65-53db-4aa5-b1f1-8050f95a0a4e
md"""
## Calculate Means and Stds
"""

# ╔═╡ 60480856-d759-47af-8d45-4112b1c194f6
mean_tot, std_tot = mean(new_df[!, :mean]), mean(new_df[!, :std])

# ╔═╡ 3ad16684-e053-48a7-afe9-e45a4231c15c
mean_80, std_80 = mean(df_80[!, :mean]), mean(df_80[!, :std])

# ╔═╡ b61f1e11-5391-4746-926d-82a67087826d
mean_100, std_100 = mean(df_100[!, :mean]), mean(df_100[!, :std])

# ╔═╡ d19413f5-be8a-40fb-9dd2-4571e4bc74cd
mean_120, std_120 = mean(df_120[!, :mean]), mean(df_120[!, :std])

# ╔═╡ 0f06986b-aee4-4c30-b89f-988ef9745355
mean_135, std_135 = mean(df_135[!, :mean]), mean(df_135[!, :std])

# ╔═╡ d0746fef-10f3-4e77-ad62-99e2c6da5b2f
md"""
## Save Means and Stds
"""

# ╔═╡ 183a40b4-04e6-4215-92d9-c6f820c90482
df_save = DataFrame(;
    mean_80=mean_80,
    std_80=std_80,
    mean_100=mean_100,
    std_100=std_100,
    mean_120=mean_120,
    std_120=std_120,
    mean_135=mean_135,
    std_135=std_135,
)

# ╔═╡ Cell order:
# ╠═eb0d0540-d2f6-11ec-290a-d728e7dffe5c
# ╠═dea06aa6-20e3-4099-8d06-060e7970a535
# ╟─b76f95ef-cd93-4e50-994d-ad3fa8433441
# ╠═36b5fe08-1ac1-4734-ba7e-64157c2be717
# ╠═a859a9dc-8cb0-4ade-b340-83edfc471927
# ╠═57e199ba-86a7-4ae2-8ee0-284fdce05dea
# ╠═75800753-24ca-4fb2-8016-6ae43efad58d
# ╠═f0b6724e-9c9f-4998-bb09-5c66345640a4
# ╟─8d20427c-d5a5-46c2-ba78-97138b966776
# ╠═45e74fa4-f0c5-457a-8a4b-b433a381b930
# ╠═c190bc16-2a17-4466-9335-0371e7213326
# ╟─96c5bc65-53db-4aa5-b1f1-8050f95a0a4e
# ╠═60480856-d759-47af-8d45-4112b1c194f6
# ╠═3ad16684-e053-48a7-afe9-e45a4231c15c
# ╠═b61f1e11-5391-4746-926d-82a67087826d
# ╠═d19413f5-be8a-40fb-9dd2-4571e4bc74cd
# ╠═0f06986b-aee4-4c30-b89f-988ef9745355
# ╟─d0746fef-10f3-4e77-ad62-99e2c6da5b2f
# ╠═183a40b4-04e6-4215-92d9-c6f820c90482

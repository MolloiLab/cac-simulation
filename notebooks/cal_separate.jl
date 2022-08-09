### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ 5b6528ec-f640-11ec-1a28-7b6a8e869196
# ╠═╡ show_logs = false
begin
    let
        using Pkg
        Pkg.activate(mktempdir())
        Pkg.Registry.update()
        Pkg.add("PlutoUI")
        Pkg.add("CairoMakie")
        Pkg.add("MAT")
        Pkg.add("Statistics")
        Pkg.add("DataFrames")
        Pkg.add("CSV")
        Pkg.add(; url="https://github.com/JuliaHealth/DICOM.jl")
        Pkg.add(; url="https://github.com/Dale-Black/DICOMUtils.jl")
    end

    using PlutoUI
    using CairoMakie
    using MAT
    using Statistics
    using DataFrames
    using CSV
    using DICOM
    using DICOMUtils
end

# ╔═╡ ec8aedd4-ab5c-495b-916a-2784cacad883
TableOfContents()

# ╔═╡ e9f3b804-eb8d-4a05-b902-66b9c0a28d78
md"""
## Save Calibrations
"""

# ╔═╡ 28af2f7b-e201-4004-ad20-7fd8573568c6
energies = [80, 100, 120, 135]

# ╔═╡ 733205c8-3fc7-4ecf-9454-90800788271b
densities = [25, 50, 100, 200, 400, 800]

# ╔═╡ a6d8c5ba-80c2-4366-a285-3d782b8c4c91
for e in energies
    for d in densities
        path = string(
            "/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/mat_files/calibration/",
            d,
            "rod",
            e,
            ".mat",
        )
        vars1 = matread(path)
        array1 = vars1[string("I")]
        array1 = Int16.(round.(array1))

        dcm_path = "/Users/daleblack/Google Drive/Datasets/Canon_Aquilion_One_Vision/Large_rep1/96E1EB4F"

        dcm = dcm_parse(dcm_path)
        dcm[tag"Pixel Data"] = array1
        dcm[tag"Instance Number"] = d
        dcm[tag"Rows"] = size(array1, 1)
        dcm[tag"Columns"] = size(array1, 2)

        output_root1 = string(
            "/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/cal_images"
        )
        if !isdir(output_root1)
            mkdir(output_root1)
        end

        output_root2 = string(output_root1, "/", e)
        if !isdir(output_root2)
            mkdir(output_root2)
        end

        output_dcm = string(output_root2, "/", d, ".dcm")
        dcm_write(output_dcm, dcm)
    end
end

# ╔═╡ 87bcb701-d356-4ccb-9660-4a6440c22207
md"""
## Check DICOMs
"""

# ╔═╡ 9189c9f6-924e-4cf8-9225-3361fdea10c9
root = "/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/cal_images/"

# ╔═╡ 6af05761-1931-4251-9396-df3194b55b58
pt = [185, 315]

# ╔═╡ 809ea1d6-c149-4080-a62a-29208125cbe4
pad = 10

# ╔═╡ 955a54d3-1b98-4baf-9d0a-42f97c28051d
begin
    dfs = []
    for e in energies
        pth = string(root, e, "/")
        header, dcm_array, slice_thick_ori1 = dcm_reader(pth)
        slice1 = dcm_array[:, :, 1]
        slice2 = dcm_array[:, :, 2]
        slice3 = dcm_array[:, :, 3]
        slice4 = dcm_array[:, :, 4]
        slice5 = dcm_array[:, :, 5]
        slice6 = dcm_array[:, :, 6]

        _25 = mean(slice1[(pt[1] - pad):(pt[1] + pad), (pt[2] - pad):(pt[2] + pad)])
        _50 = mean(slice2[(pt[1] - pad):(pt[1] + pad), (pt[2] - pad):(pt[2] + pad)])
        _100 = mean(slice3[(pt[1] - pad):(pt[1] + pad), (pt[2] - pad):(pt[2] + pad)])
        _200 = mean(slice4[(pt[1] - pad):(pt[1] + pad), (pt[2] - pad):(pt[2] + pad)])
        _400 = mean(slice5[(pt[1] - pad):(pt[1] + pad), (pt[2] - pad):(pt[2] + pad)])
        _800 = mean(slice6[(pt[1] - pad):(pt[1] + pad), (pt[2] - pad):(pt[2] + pad)])

        df = DataFrame(; kV=e, _25=_25, _50=_50, _100=_100, _200=_200, _400=_400, _800=_800)
        push!(dfs, df)
    end
end

# ╔═╡ e1c47aa3-e36b-4793-af2a-dd69f5c4e58b
new_df = vcat(dfs[1:length(dfs)]...)

# ╔═╡ 9b4d575c-2842-4256-8fe4-e3c653bbc1c9
if ~isdir(string(cd(pwd, ".."), "/output_new"))
    mkdir(string(cd(pwd, ".."), "/output_new"))
end

# ╔═╡ bb95d8a6-caa4-42e9-99bc-8971b05f3748
output_path = string(cd(pwd, ".."), "/output_new", "/calibrations.csv")

# ╔═╡ 42dae7f5-1c94-4635-8b1f-bdfc8b00e603
CSV.write(output_path, new_df)

# ╔═╡ Cell order:
# ╠═5b6528ec-f640-11ec-1a28-7b6a8e869196
# ╠═ec8aedd4-ab5c-495b-916a-2784cacad883
# ╟─e9f3b804-eb8d-4a05-b902-66b9c0a28d78
# ╠═28af2f7b-e201-4004-ad20-7fd8573568c6
# ╠═733205c8-3fc7-4ecf-9454-90800788271b
# ╠═a6d8c5ba-80c2-4366-a285-3d782b8c4c91
# ╟─87bcb701-d356-4ccb-9660-4a6440c22207
# ╠═9189c9f6-924e-4cf8-9225-3361fdea10c9
# ╠═6af05761-1931-4251-9396-df3194b55b58
# ╠═809ea1d6-c149-4080-a62a-29208125cbe4
# ╠═955a54d3-1b98-4baf-9d0a-42f97c28051d
# ╠═e1c47aa3-e36b-4793-af2a-dd69f5c4e58b
# ╠═9b4d575c-2842-4256-8fe4-e3c653bbc1c9
# ╠═bb95d8a6-caa4-42e9-99bc-8971b05f3748
# ╠═42dae7f5-1c94-4635-8b1f-bdfc8b00e603

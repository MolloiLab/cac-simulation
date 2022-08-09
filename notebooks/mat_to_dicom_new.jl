### A Pluto.jl notebook ###
# v0.19.8

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try
            Base.loaded_modules[Base.PkgId(
                Base.UUID("6e696c72-6542-2067-7265-42206c756150"),
                "AbstractPlutoDingetjes",
            )].Bonds.initial_value
        catch
            b -> missing
        end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 4427b686-e049-11ec-3271-ab17a8f19bfc
# ╠═╡ show_logs = false
begin
    let
        using Pkg
        Pkg.activate(mktempdir())
        Pkg.Registry.update()
        Pkg.add("PlutoUI")
        Pkg.add("CairoMakie")
        Pkg.add("MAT")
        Pkg.add(; url="https://github.com/JuliaHealth/DICOM.jl")
        Pkg.add(; url="https://github.com/Dale-Black/DICOMUtils.jl")
    end

    using PlutoUI
    using CairoMakie
    using MAT
    using DICOM
    using DICOMUtils
end

# ╔═╡ 20bf9bda-4a51-4f64-a299-86d6cb6ea22f
TableOfContents()

# ╔═╡ 97d46d78-70ff-42b0-8e99-51721390e171
densities = ["low", "normal"]

# ╔═╡ 726188c6-0938-4647-8aa1-fca66c983ee3
sizes = ["small", "medium", "large"]

# ╔═╡ ac21a7ea-a2c5-401e-affd-43152cdc85aa
energies = [80, 100, 120, 135]

# ╔═╡ f02d5b27-9184-45e8-b23f-df6a86867fe9
files = [1, 2, 3]

# ╔═╡ 8d9777f7-22c6-4c61-b29d-85e3300b5c45
for density in densities
    for _size in sizes
        for energy in energies
            for file in files
                file_inserts = file + 3
                ENERGY = energy
                ROD = string("QRM", ENERGY, "rod.mat")
                VESSEL = string("QRM", ENERGY, "vessels.mat")
                BASE_PATH = string(
                    "/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/mat_files/new_exposure/",
                    _size,
                    "/",
                    density,
                    "/",
                )

                path = string(BASE_PATH, ROD)
                vars1 = matread(path)
                array1 = vars1[string("I")]
                array1 = Int16.(round.(array1))

                path2 = string(BASE_PATH, VESSEL)
                vars2 = matread(path2)
                array2 = vars2[string("I")]
                array2 = Int16.(round.(array2))

                dcm_path = "/Users/daleblack/Google Drive/Datasets/Canon_Aquilion_One_Vision/Large_rep1/96E1EB4F"

                dcm = dcm_parse(dcm_path)
                dcm[tag"Pixel Data"] = array1
                dcm[tag"Instance Number"] = file
                dcm[tag"Rows"] = size(array1, 1)
                dcm[tag"Columns"] = size(array1, 2)

                output_root1 = string(
                    "/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/images_new/"
                )
                if !isdir(output_root1)
                    mkdir(output_root1)
                end
                output_root2 = string(
                    "/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/images_new/",
                    _size,
                )
                if !isdir(output_root2)
                    mkdir(output_root2)
                end
                output_root3 = string(
                    "/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/images_new/",
                    _size,
                    "/",
                    density,
                )
                if !isdir(output_root3)
                    mkdir(output_root3)
                end
                output_root = string(
                    "/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/images_new/",
                    _size,
                    "/",
                    density,
                    "/",
                    ENERGY,
                )
                if !isdir(output_root)
                    mkdir(output_root)
                end

                output_path_rod = string(output_root, "/", file, ".dcm")
                dcm_write(output_path_rod, dcm)

                dcm2 = dcm_parse(dcm_path)
                dcm2[tag"Pixel Data"] = array2
                dcm2[tag"Instance Number"] = file_inserts
                dcm2[tag"Rows"] = size(array2, 1)
                dcm2[tag"Columns"] = size(array2, 2)

                output_path_inserts = string(output_root, "/", file_inserts, ".dcm")
                dcm_write(output_path_inserts, dcm2)
            end
        end
    end
end

# ╔═╡ 258103d8-16ff-46ec-9a93-0fb526f17d3f
md"""
## Check DICOM image(s)
"""

# ╔═╡ b4a9e9b5-fcd2-4d62-a320-b8651674064a
path = "/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/images_new/large/normal/80"

# ╔═╡ 7eff2120-6f38-43c0-bf57-6ccdc49a9f7a
dcmdir_combined = dcmdir_parse(path);

# ╔═╡ 9b7bb48e-6e79-4c06-824e-88f45bde68f0
vol_combined = load_dcm_array(dcmdir_combined);

# ╔═╡ 0b0590d0-8e45-48d4-bc70-9c7c33b3fed1
@bind c PlutoUI.Slider(1:size(vol_combined, 3); default=1, show_value=true)

# ╔═╡ a07b20bc-3891-4d43-b63b-340c50d3183f
heatmap(transpose(vol_combined[:, :, c]); colormap=:grays)

# ╔═╡ Cell order:
# ╠═4427b686-e049-11ec-3271-ab17a8f19bfc
# ╠═20bf9bda-4a51-4f64-a299-86d6cb6ea22f
# ╠═97d46d78-70ff-42b0-8e99-51721390e171
# ╠═726188c6-0938-4647-8aa1-fca66c983ee3
# ╠═ac21a7ea-a2c5-401e-affd-43152cdc85aa
# ╠═f02d5b27-9184-45e8-b23f-df6a86867fe9
# ╠═8d9777f7-22c6-4c61-b29d-85e3300b5c45
# ╟─258103d8-16ff-46ec-9a93-0fb526f17d3f
# ╠═b4a9e9b5-fcd2-4d62-a320-b8651674064a
# ╠═7eff2120-6f38-43c0-bf57-6ccdc49a9f7a
# ╠═9b7bb48e-6e79-4c06-824e-88f45bde68f0
# ╟─0b0590d0-8e45-48d4-bc70-9c7c33b3fed1
# ╠═a07b20bc-3891-4d43-b63b-340c50d3183f

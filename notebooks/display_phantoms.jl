### A Pluto.jl notebook ###
# v0.19.16

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 953b1bfb-56db-4230-ac06-b8148999058e
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.activate(".")

    using PlutoUI, Statistics, CSV, DataFrames, GLM, CairoMakie, HypothesisTests, Colors, MLJBase, DICOM, DICOMUtils, PhantomSegmentation, CalciumScoring, ImageMorphology, ImageFiltering, Noise
    using StatsBase: quantile!, rmsd
end

# ╔═╡ e32075ee-7eea-4500-8f24-c45d306cfcf4
TableOfContents()

# ╔═╡ 3eb54d12-3616-47d6-ab2c-5de69bc23b18
md"""
## Load DICOMS

All you need to do is set `base_path` once and leave it. After that, the only thing that should change is the `VENDOR`, once for every set, and the `SCAN_NUMBER`, once for each scan.
"""

# ╔═╡ afd61c12-9156-4a0e-8db6-1b95f482a5c7
begin
	VENDORS = ["80", "100", "120", "135"]
	SIZES = ["small", "medium", "large"]
	DENSITIES = ["low", "normal"]

	IMAGES = "images_new"
end

# ╔═╡ 2fe41f74-1500-4c3a-a20d-4bec392c3db0
@bind v PlutoUI.Slider(1:4; default=1, show_value=true)

# ╔═╡ 948ade07-a248-4218-859b-69d8573a2994
@bind d PlutoUI.Slider(1:2; default=1, show_value=true)

# ╔═╡ a6d483ae-c3d8-41cc-a855-b67ff2906a69
@bind s PlutoUI.Slider(1:3; default=1, show_value=true)

# ╔═╡ 577cfad8-b66e-442c-b74f-efdb68b6f647
begin
	VENDOR = VENDORS[v]
	SIZE = SIZES[s]
	DENSITY = DENSITIES[d]
end

# ╔═╡ b7d8fb07-960f-4147-8d76-32fe519dafd7
begin
	    BASE_PATH = joinpath("/Users/daleblack/Google Drive/dev/MolloiLab/cac-simulation", IMAGES, SIZE, DENSITY)
		root_path = joinpath(BASE_PATH, VENDOR)
		dcm_path_list = dcm_list_builder(root_path)
		pth = dcm_path_list[1]
		header, dcm_array, slice_thick_ori1 = dcm_reader(pth)
end;

# ╔═╡ 59615deb-e6fe-4d66-9564-3c8a0178db2f
begin
		root_path_motion = joinpath(BASE_PATH, VENDOR * "-motion")
		dcm_path_list_motion = dcm_list_builder(root_path_motion)
		pth_motion = dcm_path_list_motion[1]
		_, dcm_array_motion, _ = dcm_reader(pth_motion)
end;

# ╔═╡ 7782ea73-65b4-4745-908a-38a714398376
let
    f = Figure()
	
    ax1 = Axis(f[1, 1])
	heatmap!(transpose(dcm_array[:, :, 5]); colormap=:grays)
	hidexdecorations!(ax1)
	hideydecorations!(ax1)
	ax1.title = "Stationary"

	ax2 = Axis(f[2, 1])
	heatmap!(transpose(dcm_array_motion[:, :, 5]); colormap=:grays)
	hidexdecorations!(ax2)
	hideydecorations!(ax2)
	ax2.title = "Motion"

	save(joinpath(dirname(pwd()),"figures", "phantoms_test.png"), f, resolution = (600, 900))

	f
end

# ╔═╡ 1ea1517c-87b0-4511-b49f-822cf4be64c5
function load_dcm_inputs(directions::AbstractVector)
	
	return PlutoUI.combine() do Child
		inputs = [
			md"""
				$(directions[1]): $(
				Child(directions[1], PlutoUI.Slider(1:4))

				$(directions[2]): $(
				Child(directions[2], PlutoUI.Slider(1:3))

				$(directions[3]): $(
				Child(directions[3], PlutoUI.Slider(1:2))
			"""
		]
	end
end

# ╔═╡ 94f8d4e8-b114-418e-9798-235fb526320e
@bind speeds load_dcm_inputs(["North", "East", "South"])

# ╔═╡ Cell order:
# ╠═953b1bfb-56db-4230-ac06-b8148999058e
# ╠═e32075ee-7eea-4500-8f24-c45d306cfcf4
# ╟─3eb54d12-3616-47d6-ab2c-5de69bc23b18
# ╠═afd61c12-9156-4a0e-8db6-1b95f482a5c7
# ╠═577cfad8-b66e-442c-b74f-efdb68b6f647
# ╠═b7d8fb07-960f-4147-8d76-32fe519dafd7
# ╠═59615deb-e6fe-4d66-9564-3c8a0178db2f
# ╠═2fe41f74-1500-4c3a-a20d-4bec392c3db0
# ╠═948ade07-a248-4218-859b-69d8573a2994
# ╠═a6d483ae-c3d8-41cc-a855-b67ff2906a69
# ╠═7782ea73-65b4-4745-908a-38a714398376
# ╠═94f8d4e8-b114-418e-9798-235fb526320e
# ╠═1ea1517c-87b0-4511-b49f-822cf4be64c5

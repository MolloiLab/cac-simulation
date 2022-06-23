### A Pluto.jl notebook ###
# v0.19.8

using Markdown
using InteractiveUtils

# ╔═╡ e4efd54a-f30b-11ec-0893-1d541db42058
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
		Pkg.add(url="https://github.com/JuliaHealth/DICOM.jl")
		Pkg.add(url="https://github.com/Dale-Black/DICOMUtils.jl")
		Pkg.add(url="https://github.com/Dale-Black/PhantomSegmentation.jl")
		Pkg.add(url="https://github.com/Dale-Black/CalciumScoring.jl")
		Pkg.add("OrderedCollections")
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
	using OrderedCollections
end

# ╔═╡ a19a47e4-9601-4d1a-98b7-f632bf879ca8
TableOfContents()

# ╔═╡ 06cfea73-89c1-4908-b40c-573a3799ee46
TYPE = "calibration"

# ╔═╡ d1d6a762-5e45-4a51-a9bf-053d671ef634
VENDERS = ["80", "100", "120", "135"]

# ╔═╡ 06bf1275-374b-4520-afbf-46231e57345a
SIZES = ["small"]

# ╔═╡ ecd0330c-472d-447b-991e-827dfbc2177f
DENSITIES = ["low", "normal"]

# ╔═╡ 01e699f9-4d9b-407b-89c1-8f6e0da1ef1f
blurs = [0]

# ╔═╡ f063b6dc-b3a2-42fc-b605-e278636d413a
begin
	dfs = []
	cals = OrderedDict()
	for VENDER in VENDERS
		high_dens = []
		med_dens = []
		low_dens = []
		high_dens100 = []
		med_dens50 = []
		low_dens25 = []
		for blur in blurs
			for SIZE in SIZES
				for DENSITY in DENSITIES
					SCAN_NUMBER = 1
					BASE_PATH = string("/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/images_new/", SIZE, "/", DENSITY, "/")
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
					masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3)÷2)
				
					# Segment Calcium Rod
					local thresh
					if DENSITY == "low"
						thresh = 55
					elseif DENSITY ==  "normal"
						thresh = 130
					end
	
					@info blur, DENSITY, SIZE, VENDER, thresh
					
					calcium_image, slice_CCI, quality_slice, cal_rod_slice = mask_rod(masked_array, header; calcium_threshold=thresh)
	
					root_new = string("/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/julia_arrays/", SIZE, "/") 
					mask_L_HD = Array(CSV.read(string(root_new, "mask_L_HD.csv"), DataFrame; header=false))
					mask_M_HD = Array(CSV.read(string(root_new, "mask_M_HD.csv"), DataFrame; header=false))
					mask_S_HD = Array(CSV.read(string(root_new, "mask_S_HD.csv"), DataFrame; header=false))
					mask_L_MD = Array(CSV.read(string(root_new, "mask_L_MD.csv"), DataFrame;
					header=false))
					mask_M_MD = Array(CSV.read(string(root_new, "mask_M_MD.csv"), DataFrame; header=false))
					mask_S_MD = Array(CSV.read(string(root_new, "mask_S_MD.csv"), DataFrame; header=false))
					mask_L_LD = Array(CSV.read(string(root_new, "mask_L_LD.csv"), DataFrame; header=false))
					mask_M_LD = Array(CSV.read(string(root_new, "mask_M_LD.csv"), DataFrame; header=false))
					mask_S_LD = Array(CSV.read(string(root_new, "mask_S_LD.csv"), DataFrame; header=false))
	
					arr = masked_array[:, :, 4:6]
					single_arr = masked_array[:, :, 5]
					pixel_size = DICOMUtils.get_pixel_size(header)
	
					if DENSITY == "normal"
						# Calculate mean intensities for calibration
						erode_mask_L_HD = erode(erode(mask_L_HD))
						mean_HD = mean(single_arr[erode_mask_L_HD])
						push!(high_dens, mean_HD)
		
						erode_mask_L_MD = erode(erode(mask_L_MD))
						mean_MD = mean(single_arr[erode_mask_L_MD])
						push!(med_dens, mean_MD)
		
						erode_mask_L_LD = erode(erode(mask_L_LD))
						mean_LD = mean(single_arr[erode_mask_L_LD])
						push!(low_dens, mean_LD)
					else
						erode_mask_L_HD = erode(erode(mask_L_HD))
						mean_HD = mean(single_arr[erode_mask_L_HD])
						push!(high_dens100, mean_HD)
		
						erode_mask_L_MD = erode(erode(mask_L_MD))
						mean_MD = mean(single_arr[erode_mask_L_MD])
						push!(med_dens50, mean_MD)
		
						erode_mask_L_LD = erode(erode(mask_L_LD))
						mean_LD = mean(single_arr[erode_mask_L_LD])
						push!(low_dens25, mean_LD)
					end
				end
			end
		end
		cal = (mean(low_dens25), mean(med_dens50), mean(high_dens100), mean(low_dens), mean(med_dens), mean(high_dens))
		cals["$(VENDER)"] = cal
	end
end

# ╔═╡ b4cfa24f-0f28-4c75-b26f-bfcd1c56e773
cals

# ╔═╡ Cell order:
# ╠═e4efd54a-f30b-11ec-0893-1d541db42058
# ╠═a19a47e4-9601-4d1a-98b7-f632bf879ca8
# ╠═06cfea73-89c1-4908-b40c-573a3799ee46
# ╠═d1d6a762-5e45-4a51-a9bf-053d671ef634
# ╠═06bf1275-374b-4520-afbf-46231e57345a
# ╠═ecd0330c-472d-447b-991e-827dfbc2177f
# ╠═01e699f9-4d9b-407b-89c1-8f6e0da1ef1f
# ╠═f063b6dc-b3a2-42fc-b605-e278636d413a
# ╠═b4cfa24f-0f28-4c75-b26f-bfcd1c56e773

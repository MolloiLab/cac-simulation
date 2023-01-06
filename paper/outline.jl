### A Pluto.jl notebook ###
# v0.19.18

using Markdown
using InteractiveUtils

# ╔═╡ 48f50ffe-6a1b-4f4d-a515-d6a63d78eef8
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using Images, FileIO, PlutoUI, CSV, DataFrames, PrettyTables, Latexify
end

# ╔═╡ d3043d93-7792-4005-ad8b-9b981a6f5aa0
TableOfContents()

# ╔═╡ 8ffa0667-835f-4503-afb9-6dc5772d5775
md"""
# Abstract
"""

# ╔═╡ 8ace0d7f-1863-49a0-897c-5fabc85b1a13
md"""
# 1. Introduction
"""

# ╔═╡ 2453932a-7525-419d-8ba2-f8d1289a0c0f
md"""
# 2. Method
"""

# ╔═╡ 946fdf4b-519a-4b7c-a025-3f8c0e3d594d
md"""
## 2.1 - Simulation
"""

# ╔═╡ 208eaf04-df19-4c39-afb2-3a78a07486f7
load(joinpath(dirname(pwd()), "figures", "phantom materials.png"))

# ╔═╡ ac21720a-ea27-490a-9bcd-0ef4972bd5d3
md"""
Fig. 1 Shows a sketch of the simulated phantom with the colors highlighting the different materials in the simulated phantoms.
"""

# ╔═╡ bc049335-7083-4626-9625-b82b16abc1ad
md"""
## 2.2 - Physical Phantom
"""

# ╔═╡ b599acd8-4808-4e02-8604-281f11ebdf89
load(joinpath(dirname(pwd()), "figures", "phantom slices.png"))

# ╔═╡ dac84773-ef02-42aa-885e-cc6415174840
md"""
Fig. 2 Shows axial slice views of a of various phantoms included in this study. (A) shows a simulated small phantom with normal-density inserts (200, 400, 800 mg/cm^3) and a tube voltage of 120 kV. (B) shows a simulated large phantom with low-density inserts (25, 50, 100 mg/cm^3) and a tube voltage of 120 kV with simulated motion.(C) shows an axial slice of the calibration rod (200 mg/cm^3) for a simulated medium phantom with a tube voltage of 100 kV. (D) shows an axial slice view of a physical QRM Thorax Phantom with a Cardio Calcification Insert. The red arrow shows simulated beam hardening artifact, and the blue arrows show simulated streaking artifact caused by motion.
"""

# ╔═╡ 6cf1b2c9-325e-4237-9f78-4c1266b9e521
md"""
## 2.3 - Agatston Scoring
"""

# ╔═╡ f0be8fb5-f867-4348-9a3c-af3578b26285
md"""
## 2.4 - Integrated Calcium Mass
"""

# ╔═╡ b2f7d165-4850-47ee-b596-618d26bd333f
load(joinpath(dirname(pwd()), "figures", "integrated-volume fraction methods.png"))

# ╔═╡ 7fd911bf-7e9e-4705-8217-b90af882aeae
md"""
Fig 3. Shows two identical simulated vessel lumens with ROIs for calcium measurement.(A) Shows the ROIs needed for the integrated calcium mass technique. The central ROI that yields ``S_{Obj}`` and the background ROI that yields ``S_{Bkg}`` are unaffected by the partial volume effect, while the object ROI used to calculate ``I`` is affected by the partial volume effect. (B) Shows the ROIs needed for the volume fraction technique and a zoomed-in portion of the simulated vessel lumen, where ``i`` is the individual voxel intensity. ``S_{Bkg}`` is a ring-like background ROI unaffected by the partial volume effect. ``S_{Obj}`` is the mean intensity of known calcium density obtained from a calibration rod unaffected by the partial volume effect.
"""

# ╔═╡ c9377e1e-769c-430d-8914-4168b8750362
md"""
## 2.5 - Volume Fraction Mass
"""

# ╔═╡ 461c7e3f-2945-487c-b3d0-6b173b962129
md"""
## 2.6 - Spatially Weighted Calcium Scoring
"""

# ╔═╡ 65fa25fa-852a-4ff9-ac51-ec072112ef90
md"""
## 2.7 - Statistical Analysis
"""

# ╔═╡ 345bc398-5cd7-48ca-9e96-192893213266
md"""
# 3. Results
"""

# ╔═╡ f740d3da-21a5-4cc6-b235-2bd96bd90594
md"""
## 3.1 - Simulated Phantom
"""

# ╔═╡ c368e80b-a9ef-4141-94a2-0cbdbd794ae6
md"""
### 3.1.1 - Accuracy
"""

# ╔═╡ b60e850e-22fb-43b9-839a-bebdfa90e018
# load(joinpath(dirname(pwd()), "figures", "stationary", "accuracy_normal.png"))

# ╔═╡ 4c043244-e54c-4f3d-9e79-ee445faf6525
# md"""
# Fig. 5 Shows the linear regression analysis comparing measured calcium to the known calcium for the normal density (200, 400, 800 mg/ cm3) stationary phantoms. Fig 5A shows the results of integrated calcium mass. Fig 5B shows the results of the volume fraction method. Fig 5C shows the results of Agatston mass scoring. The best fit line along with the root mean squared error and root mean squared deviation values are shown in each plot. 
# """

# ╔═╡ 95730773-705e-46fb-a8fe-ab2d4333dbd2
load(joinpath(dirname(pwd()), "figures", "stationary", "accuracy_low.png"))

# ╔═╡ 5bb0f56b-7558-40fc-bed0-4d25211aa9ff
md"""
Fig. 4 Shows the linear regression analysis comparing measured calcium to the known calcium for the low-density (25, 50, 100 mg/ cm3) stationary phantoms. Every tube voltage (80, 100, 120, 135 kV) and size (small, medium, large) is included in the analysis. (A) shows the results of integrated calcium mass. (B) shows the results of the volume fraction method. (C) shows the results of Agatston mass scoring. The best fit line, along with the root mean squared error and root mean squared deviation values are shown in each plot. 
"""

# ╔═╡ 70ff5fd1-da44-4cb8-95b6-6680ef4002f0
md"""
### 3.1.2 - Reproducibility
"""

# ╔═╡ 6bd02498-6cc2-4311-9370-78bfd9dca0da
load(joinpath(dirname(pwd()), "figures", "stationary", "reproducibility.png"))

# ╔═╡ baa3757e-f1bb-416c-83bf-8977d0a93cdc
md"""
Fig. 5 Shows reproducibility measurements for all four scoring techniques on the stationary phantoms. Every tube voltage (80, 100, 120, 135 kV), size (small, medium, large), and density (low, normal) is included in the analysis. The measurements from the first set of images were plotted against the second set of images for each technique. Integrated mass (A), Volume fraction (B), Agatston mass scoring (C), and spatially weighted calcium scoring (D) are shown along with the root mean squared error and root mean squared deviation values.
"""

# ╔═╡ 9b2b3a55-905f-4aa8-8d6b-136ac8f140de
md"""
### 3.1.3 - Sensitivity and Specificity
"""

# ╔═╡ ad82f85c-fe65-44a8-a1e8-b9e7665fcdb7
load(joinpath(dirname(pwd()), "figures", "stationary", "sensitivity_specificity.png"))

# ╔═╡ 36f219cb-560e-49fc-b8fc-07f1c19bf500
md"""
Fig 6. Shows the percentage of false-negative (CAC=0) and false-positive (CAC>0) scores on the stationary phantoms. Every tube voltage (80, 100, 120, 135 kV), size (small, medium, large), and density (low, normal) is included in the analysis.
"""

# ╔═╡ f6c2118f-29c5-4dfa-bf8a-4ee30f8be196
# md"""
# ## 3.2 - Motion
# """

# ╔═╡ 03b59f85-f72b-4100-8aa7-3767fca347c4
# md"""
# ### 3.2.1 - Accuracy
# """

# ╔═╡ 6f614da3-da89-41d2-bf1b-9cc7761ccc99
# load(joinpath(dirname(pwd()), "figures", "motion", "accuracy_normal.png"))

# ╔═╡ fd92053f-c8fd-46f7-a6f3-4a5f41531a98
# md"""
# Fig. 9 Shows the linear regression analysis comparing measured calcium to the known calcium for the normal density (200, 400, 800 mg/ cm3) motion-affected phantoms. Fig 9A shows the results of integrated calcium mass. Fig 9B shows the results of the volume fraction method. Fig 9C shows the results of Agatston mass scoring. The best fit line, along with the root mean squared error and root mean squared deviation values, are shown in each plot. 
# """

# ╔═╡ 9520d761-f743-4fcd-860c-fdf80b996d23
# load(joinpath(dirname(pwd()), "figures", "motion", "accuracy_low.png"))

# ╔═╡ 9e42ea7f-3c72-4cdc-b6be-177df1e46f17
# md"""
# Fig. 10 Shows the linear regression analysis comparing measured calcium to the known calcium for the low-density (25, 50, 100 mg/ cm3) motion-affected phantoms. Fig 10A shows the results of integrated calcium mass. Fig 10B shows the results of the volume fraction method. Fig 10C shows the results of Agatston mass scoring. The best fit line, along with the root mean squared error and root mean squared deviation values, are shown in each plot. 
# """

# ╔═╡ caf409b6-eac7-422e-8480-04a55ad19d3c
# md"""
# ### 3.1.2 - Reproducibility
# """

# ╔═╡ 831bb50b-e268-476f-a87b-f532c847aa60
# load(joinpath(dirname(pwd()), "figures", "motion", "reproducibility.png"))

# ╔═╡ 70c9020a-79e7-488c-ad97-f853e659b322
# md"""
# Fig. 11 Shows reproducibility measurements for all four scoring techniques on the motion-affected phantoms. The measurements from the first set of images were plotted against the second set of images for each technique. Integrated mass (A), Volume fraction (B), Agatston mass scoring (C), and spatially weighted calcium scoring (D) are shown along with the root mean squared error and root mean squared deviation values.
# """

# ╔═╡ 62532e8f-5e71-489a-89a1-70a56f1b1d7a
# md"""
# ### 3.1.3 - Sensitivity and Specificity
# """

# ╔═╡ b55ce920-715b-4549-a6d0-2ecce40d21bd
# load(joinpath(dirname(pwd()), "figures", "motion", "sensitivity_specificity.png"))

# ╔═╡ b58573a5-7308-4ce1-b091-cbc447e134b8
# md"""
# Fig 12. Shows the percentage of false-negative (CAC=0) and false-positive (CAC>0) scores for a region of pure background on the motion-affected phantoms.
# """

# ╔═╡ 36f6c9c0-8963-449a-ad4d-33bc6fa6d44a
md"""
## 3.2 - Physical Phantom
"""

# ╔═╡ 812cb8f8-0ae6-4cb8-a15c-b34a6fc0e196
md"""
### 3.2.1 - Accuracy
"""

# ╔═╡ 81744560-96ef-468d-a17b-91598f1129fe
load(joinpath(dirname(pwd()), "figures", "physical", "accuracy.png"))

# ╔═╡ 4908af35-e8a0-4a78-89ac-a8280a3512ff
md"""
Fig. 7 Shows the linear regression analysis comparing measured calcium to the known calcium for the normal density (200, 400, 800 mg/ cm3) physical phantom scans at a tube voltage of 120 kV. (A) shows the results of integrated calcium mass. (B)shows the results of the volume fraction method. (C) shows the results of Agatston mass scoring. The best fit line, along with the root mean squared error and root mean squared deviation values, are shown in each plot. 
"""

# ╔═╡ c66b50d3-8a0b-4f78-8e1d-324ece091d98
md"""
### 3.2.2 - Reproducibility

*Note: The reproducibility results show a lower RMSE and RMSD for Agatston compared to Integrated and Volume Fraction. This is misleading because Agatston underestimates the mass by almost 50%, so the reproducibility comparison between smaller numbers is of course going to result in a smaller RMSE/D than Integrated/Volume Fraction techniques.*
"""

# ╔═╡ 21c2ea30-2b7c-4e10-b8f6-b616be8e8ce2
load(joinpath(dirname(pwd()), "figures", "physical", "reproducibility.png"))

# ╔═╡ aae2d4d7-248c-4ec6-a258-f88ab26dd134
md"""
Fig. 8 Shows reproducibility measurements on the physical phantom scans at a tube voltage of 120 kV. Measurements from two scans were plotted against a second set of measurements from two different scans at different angles. Integrated mass (A), Volume fraction (B), Agatston mass scoring (C) are shown along with the root mean squared error and root mean squared deviation values.
"""

# ╔═╡ c2da7cda-f636-4c97-aca3-09a649e29731
md"""
### 3.2.3 - Sensitivity and Specificity
"""

# ╔═╡ b009fdb8-f059-4d8b-8069-158daf50820c
load(joinpath(dirname(pwd()), "figures", "physical", "sensitivity_specificity.png"))

# ╔═╡ 583f9090-a2b6-4c9c-bcd3-cd588018a30f
md"""
Fig 9. Shows the percentage of false-negative (CAC=0) and false-positive (CAC>0) scores for the physical phantom scans at a tube voltage of 120 kV.
"""

# ╔═╡ a5db1e0e-19f8-48c0-8fd9-006e115c056e
begin
	df_stationary = CSV.read(joinpath(dirname(pwd()), "figures", "stationary", "summary.csv"), DataFrame)
	insertcols!(df_stationary, 1, "Phantom"=>repeat(["Stationary"], size(df_stationary)[1]))
	
	df_motion = CSV.read(joinpath(dirname(pwd()), "figures", "motion", "summary.csv"), DataFrame)
	insertcols!(df_motion, 1, "Phantom"=>repeat(["Motion"], size(df_motion)[1]))

	df_physical = CSV.read(joinpath(dirname(pwd()), "figures", "physical", "summary.csv"), DataFrame)
	insertcols!(df_physical, 1, "Phantom"=>repeat(["Physical"], size(df_physical)[1]))

	_dfs = append!(df_stationary, df_motion, cols=:union)
	dfs = append!(_dfs, df_physical, cols=:union)
end

# ╔═╡ e2a5e3e7-6088-47cb-a475-ae4cf75a2597
CSV.write(joinpath(pwd(), "results-total-summary.csv"), dfs)

# ╔═╡ 75aaf5be-9013-4b92-99a4-a862d2e5343a
md"""
# 4. Discussion
"""

# ╔═╡ 269f772d-bbd2-4224-aeea-711a308b1eb8
md"""
# 5. Conclusion
"""

# ╔═╡ fa9db86b-0894-4ce5-b431-8ae4b6dacf24
md"""
# Acknowledgements
"""

# ╔═╡ 1214a60f-d5d8-4a26-a284-464605dfccbc
md"""
# References
"""

# ╔═╡ Cell order:
# ╠═48f50ffe-6a1b-4f4d-a515-d6a63d78eef8
# ╠═d3043d93-7792-4005-ad8b-9b981a6f5aa0
# ╟─8ffa0667-835f-4503-afb9-6dc5772d5775
# ╟─8ace0d7f-1863-49a0-897c-5fabc85b1a13
# ╟─2453932a-7525-419d-8ba2-f8d1289a0c0f
# ╟─946fdf4b-519a-4b7c-a025-3f8c0e3d594d
# ╟─208eaf04-df19-4c39-afb2-3a78a07486f7
# ╟─ac21720a-ea27-490a-9bcd-0ef4972bd5d3
# ╟─bc049335-7083-4626-9625-b82b16abc1ad
# ╟─b599acd8-4808-4e02-8604-281f11ebdf89
# ╟─dac84773-ef02-42aa-885e-cc6415174840
# ╟─6cf1b2c9-325e-4237-9f78-4c1266b9e521
# ╟─f0be8fb5-f867-4348-9a3c-af3578b26285
# ╟─b2f7d165-4850-47ee-b596-618d26bd333f
# ╟─7fd911bf-7e9e-4705-8217-b90af882aeae
# ╟─c9377e1e-769c-430d-8914-4168b8750362
# ╟─461c7e3f-2945-487c-b3d0-6b173b962129
# ╟─65fa25fa-852a-4ff9-ac51-ec072112ef90
# ╟─345bc398-5cd7-48ca-9e96-192893213266
# ╟─f740d3da-21a5-4cc6-b235-2bd96bd90594
# ╟─c368e80b-a9ef-4141-94a2-0cbdbd794ae6
# ╟─b60e850e-22fb-43b9-839a-bebdfa90e018
# ╟─4c043244-e54c-4f3d-9e79-ee445faf6525
# ╟─95730773-705e-46fb-a8fe-ab2d4333dbd2
# ╟─5bb0f56b-7558-40fc-bed0-4d25211aa9ff
# ╟─70ff5fd1-da44-4cb8-95b6-6680ef4002f0
# ╟─6bd02498-6cc2-4311-9370-78bfd9dca0da
# ╟─baa3757e-f1bb-416c-83bf-8977d0a93cdc
# ╟─9b2b3a55-905f-4aa8-8d6b-136ac8f140de
# ╟─ad82f85c-fe65-44a8-a1e8-b9e7665fcdb7
# ╟─36f219cb-560e-49fc-b8fc-07f1c19bf500
# ╟─f6c2118f-29c5-4dfa-bf8a-4ee30f8be196
# ╟─03b59f85-f72b-4100-8aa7-3767fca347c4
# ╟─6f614da3-da89-41d2-bf1b-9cc7761ccc99
# ╟─fd92053f-c8fd-46f7-a6f3-4a5f41531a98
# ╟─9520d761-f743-4fcd-860c-fdf80b996d23
# ╟─9e42ea7f-3c72-4cdc-b6be-177df1e46f17
# ╟─caf409b6-eac7-422e-8480-04a55ad19d3c
# ╟─831bb50b-e268-476f-a87b-f532c847aa60
# ╟─70c9020a-79e7-488c-ad97-f853e659b322
# ╟─62532e8f-5e71-489a-89a1-70a56f1b1d7a
# ╟─b55ce920-715b-4549-a6d0-2ecce40d21bd
# ╟─b58573a5-7308-4ce1-b091-cbc447e134b8
# ╟─36f6c9c0-8963-449a-ad4d-33bc6fa6d44a
# ╟─812cb8f8-0ae6-4cb8-a15c-b34a6fc0e196
# ╟─81744560-96ef-468d-a17b-91598f1129fe
# ╟─4908af35-e8a0-4a78-89ac-a8280a3512ff
# ╟─c66b50d3-8a0b-4f78-8e1d-324ece091d98
# ╟─21c2ea30-2b7c-4e10-b8f6-b616be8e8ce2
# ╟─aae2d4d7-248c-4ec6-a258-f88ab26dd134
# ╟─c2da7cda-f636-4c97-aca3-09a649e29731
# ╠═b009fdb8-f059-4d8b-8069-158daf50820c
# ╟─583f9090-a2b6-4c9c-bcd3-cd588018a30f
# ╟─a5db1e0e-19f8-48c0-8fd9-006e115c056e
# ╟─e2a5e3e7-6088-47cb-a475-ae4cf75a2597
# ╟─75aaf5be-9013-4b92-99a4-a862d2e5343a
# ╟─269f772d-bbd2-4224-aeea-711a308b1eb8
# ╟─fa9db86b-0894-4ce5-b431-8ae4b6dacf24
# ╟─1214a60f-d5d8-4a26-a284-464605dfccbc

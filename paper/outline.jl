### A Pluto.jl notebook ###
# v0.19.16

using Markdown
using InteractiveUtils

# ╔═╡ 48f50ffe-6a1b-4f4d-a515-d6a63d78eef8
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using Images, FileIO, PlutoUI
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

# ╔═╡ 052c0b93-c3b7-40c1-8aec-062c3e9a224d
load(joinpath(dirname(pwd()), "figures", "phantom geometries.png"))

# ╔═╡ f722fa66-c6c7-4d5b-8ee9-c552f703a2ba
md"""
### 2.1.1 - Motion Blurring
"""

# ╔═╡ bc049335-7083-4626-9625-b82b16abc1ad
md"""
## 2.2 - Physical Phantom
"""

# ╔═╡ 6cf1b2c9-325e-4237-9f78-4c1266b9e521
md"""
## 2.3 - Agatston Scoring
"""

# ╔═╡ f0be8fb5-f867-4348-9a3c-af3578b26285
md"""
## 2.4 Integrated Calcium Mass
"""

# ╔═╡ 63b474d3-619c-44a5-81ea-e1a0baa69a6d
md"""
### 2.4.1 - Calibration
"""

# ╔═╡ c9377e1e-769c-430d-8914-4168b8750362
md"""
## 2.5 - Volume Fraction Mass
"""

# ╔═╡ 65fa25fa-852a-4ff9-ac51-ec072112ef90
md"""
## 2.6 - Statistical Analysis
"""

# ╔═╡ 345bc398-5cd7-48ca-9e96-192893213266
md"""
# 3. Results
"""

# ╔═╡ f740d3da-21a5-4cc6-b235-2bd96bd90594
md"""
## 3.1 - Stationary
"""

# ╔═╡ c368e80b-a9ef-4141-94a2-0cbdbd794ae6
md"""
### 3.1.1 - Accuracy
"""

# ╔═╡ b60e850e-22fb-43b9-839a-bebdfa90e018
load(joinpath(dirname(pwd()), "figures", "stationary", "accuracy_normal.png"))

# ╔═╡ 95730773-705e-46fb-a8fe-ab2d4333dbd2
load(joinpath(dirname(pwd()), "figures", "stationary", "accuracy_low.png"))

# ╔═╡ 70ff5fd1-da44-4cb8-95b6-6680ef4002f0
md"""
### 3.1.2 - Reproducibility
"""

# ╔═╡ 6bd02498-6cc2-4311-9370-78bfd9dca0da
load(joinpath(dirname(pwd()), "figures", "stationary", "reproducibility.png"))

# ╔═╡ 9b2b3a55-905f-4aa8-8d6b-136ac8f140de
md"""
### 3.1.3 - Sensitivity and Specificity
"""

# ╔═╡ ad82f85c-fe65-44a8-a1e8-b9e7665fcdb7
load(joinpath(dirname(pwd()), "figures", "stationary", "sensitivity_specificity.png"))

# ╔═╡ f6c2118f-29c5-4dfa-bf8a-4ee30f8be196
md"""
## 3.2 - Motion
"""

# ╔═╡ 03b59f85-f72b-4100-8aa7-3767fca347c4
md"""
### 3.2.1 - Accuracy
"""

# ╔═╡ 6f614da3-da89-41d2-bf1b-9cc7761ccc99
load(joinpath(dirname(pwd()), "figures", "motion", "accuracy_normal.png"))

# ╔═╡ 9520d761-f743-4fcd-860c-fdf80b996d23
load(joinpath(dirname(pwd()), "figures", "motion", "accuracy_low.png"))

# ╔═╡ caf409b6-eac7-422e-8480-04a55ad19d3c
md"""
### 3.1.2 - Reproducibility
"""

# ╔═╡ 831bb50b-e268-476f-a87b-f532c847aa60
load(joinpath(dirname(pwd()), "figures", "motion", "reproducibility.png"))

# ╔═╡ 62532e8f-5e71-489a-89a1-70a56f1b1d7a
md"""
### 3.1.3 - Sensitivity and Specificity
"""

# ╔═╡ b55ce920-715b-4549-a6d0-2ecce40d21bd
load(joinpath(dirname(pwd()), "figures", "motion", "sensitivity_specificity.png"))

# ╔═╡ 36f6c9c0-8963-449a-ad4d-33bc6fa6d44a
md"""
## 3.2 - Physical Phantom
"""

# ╔═╡ 812cb8f8-0ae6-4cb8-a15c-b34a6fc0e196
md"""
### 3.1.1 - Accuracy
"""

# ╔═╡ 81744560-96ef-468d-a17b-91598f1129fe
load(joinpath(dirname(pwd()), "figures", "physical", "accuracy.png"))

# ╔═╡ c66b50d3-8a0b-4f78-8e1d-324ece091d98
md"""
### 3.1.2 - Reproducibility
"""

# ╔═╡ 21c2ea30-2b7c-4e10-b8f6-b616be8e8ce2
load(joinpath(dirname(pwd()), "figures", "physical", "reproducibility.png"))

# ╔═╡ c2da7cda-f636-4c97-aca3-09a649e29731
md"""
### 3.1.3 - Sensitivity and Specificity
"""

# ╔═╡ b009fdb8-f059-4d8b-8069-158daf50820c
load(joinpath(dirname(pwd()), "figures", "physical", "sensitivity_specificity.png"))

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
# ╟─052c0b93-c3b7-40c1-8aec-062c3e9a224d
# ╟─f722fa66-c6c7-4d5b-8ee9-c552f703a2ba
# ╟─bc049335-7083-4626-9625-b82b16abc1ad
# ╟─6cf1b2c9-325e-4237-9f78-4c1266b9e521
# ╟─f0be8fb5-f867-4348-9a3c-af3578b26285
# ╟─63b474d3-619c-44a5-81ea-e1a0baa69a6d
# ╟─c9377e1e-769c-430d-8914-4168b8750362
# ╟─65fa25fa-852a-4ff9-ac51-ec072112ef90
# ╟─345bc398-5cd7-48ca-9e96-192893213266
# ╟─f740d3da-21a5-4cc6-b235-2bd96bd90594
# ╟─c368e80b-a9ef-4141-94a2-0cbdbd794ae6
# ╟─b60e850e-22fb-43b9-839a-bebdfa90e018
# ╟─95730773-705e-46fb-a8fe-ab2d4333dbd2
# ╟─70ff5fd1-da44-4cb8-95b6-6680ef4002f0
# ╟─6bd02498-6cc2-4311-9370-78bfd9dca0da
# ╟─9b2b3a55-905f-4aa8-8d6b-136ac8f140de
# ╟─ad82f85c-fe65-44a8-a1e8-b9e7665fcdb7
# ╟─f6c2118f-29c5-4dfa-bf8a-4ee30f8be196
# ╟─03b59f85-f72b-4100-8aa7-3767fca347c4
# ╟─6f614da3-da89-41d2-bf1b-9cc7761ccc99
# ╟─9520d761-f743-4fcd-860c-fdf80b996d23
# ╟─caf409b6-eac7-422e-8480-04a55ad19d3c
# ╟─831bb50b-e268-476f-a87b-f532c847aa60
# ╟─62532e8f-5e71-489a-89a1-70a56f1b1d7a
# ╟─b55ce920-715b-4549-a6d0-2ecce40d21bd
# ╟─36f6c9c0-8963-449a-ad4d-33bc6fa6d44a
# ╟─812cb8f8-0ae6-4cb8-a15c-b34a6fc0e196
# ╟─81744560-96ef-468d-a17b-91598f1129fe
# ╟─c66b50d3-8a0b-4f78-8e1d-324ece091d98
# ╟─21c2ea30-2b7c-4e10-b8f6-b616be8e8ce2
# ╟─c2da7cda-f636-4c97-aca3-09a649e29731
# ╟─b009fdb8-f059-4d8b-8069-158daf50820c
# ╟─75aaf5be-9013-4b92-99a4-a862d2e5343a
# ╟─269f772d-bbd2-4224-aeea-711a308b1eb8
# ╟─fa9db86b-0894-4ce5-b431-8ae4b6dacf24
# ╟─1214a60f-d5d8-4a26-a284-464605dfccbc

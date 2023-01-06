### A Pluto.jl notebook ###
# v0.19.18

using Markdown
using InteractiveUtils

# ╔═╡ 34651281-01a4-46cb-b094-5d0582519c30
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using Images, FileIO, PlutoUI, CSV, DataFrames, PrettyTables, Latexify
end

# ╔═╡ ed03ef4a-e103-40c2-8a83-cbd565f9a35c
TableOfContents()

# ╔═╡ 29fdd303-5953-47c0-a789-6a8eaaa2dfc7
md"""
# Abstract
"""

# ╔═╡ 4df8f7e4-b5dc-4c97-aa62-032137e83ec4
md"""
Agatston scoring is limited by factors such as arbitrary thresholding and does not include all the calcium information in computed tomography scans. Agatston scoring’s limitations lead to problems with reproducibility and accuracy and, therefore, may not correlate well with known calcium mass. A scoring technique that removes the need for thresholding and quantifies calcium mass more accurately and reproducibly is needed.


This study adapted two techniques to create two novel calcium scoring approaches that include all the calcium information within an image. The integrated Hounsfield technique was adjusted for calcium mass quantification and the volume fraction technique. Integrated calcium mass and volume fraction calcium mass were calculated and compared to known calcium mass, along with Agatston scoring and spatially weighted calcium scoring on simulated CT scans across different patient sizes, energies, and motion levels. 

The simulation was created, based on previously validated software, to match the scanning parameters of a 320-slice CT scanner. To simulate different patient sizes, additional fat rings were added, which resulted in a small-sized phantom of 30x20cm, a medium-sized phantom of 35x25cm, and a large-sized phantom of 40x30cm. Three calcification inserts of different diameters (1, 3, and 5 mm) and different hydroxyapatite (HA) densities were placed within the phantom. All calcium scoring measurements were repeated across different kVs (80-135 kV), patient sizes (small, medium, large), calcium insert sizes (1-5 mm), and calcium insert densities (25-800 mg/cm3). Physical phantom scans, acquired by a previous group, were then used to validate the results of the simulation study.

Both integrated and volume fraction calcium mass yielded lower root mean squared error (RMSE) and deviation (RMSD) values than Agatston scoring in all accuracy comparisons. Integrated and volume fraction calcium mass also yielded lower RMSE and RMSD values than Agatston or spatially weighted calcium scoring in every reproducibility comparison. Specifically, the RMSE and RMSD of the low-density calcifications for the integrated calcium mass technique vs. known calcium mass were 0.495 mg and 0.494 mg, respectively. The RMSE and RMSD values for the volume fraction calcium mass technique were 0.585 mg and 0.575 mg, respectively. Compared to Agatston scoring, which produced RMSE and RMSD values of 3.509 mg and 2.240 mg, respectively. Reproducibility across all patient sizes, calcium densities, and tube voltages, integrated calcium mass produced RMSE and RMSD values of 0.708 mg and 0.692 mg, respectively. Volume fraction calcium mass produced RMSE and RMSD values of 0.783 mg and 0.776 mg, respectively. Spatially weighted calcium scoring produced RMSE and RMSD values of 4.722 mg and 4.714 mg, respectively. Agatston scoring produced RMSE and RMSD values of 0.982 mg and 0.980 mg, respectively. The percentage of false-negative (CAC=0) scores produced by integrated calcium mass, volume fraction calcium mass, spatially weighted calcium scoring, and Agatston scoring were 11.111%, 11.111%, 42.593%, and 38.889%, respectively. The percentage of false-positive (CAC>0) scores produced by integrated calcium mass, volume fraction calcium mass, spatially weighted calcium scoring, and Agatston scoring were 11.111%, 11.111%, 16.667, and 0.0%, respectively. 
"""

# ╔═╡ 11cff734-d9b3-48c5-b0c8-3d5ac163962c
md"""
# 1. Introduction
"""

# ╔═╡ c73b73b7-973b-487d-80a1-8dce586cc3c8
md"""
The coronary artery calcification (CAC) score is a standard atherosclerotic marker [CITE]. CAC scoring is a test performed using computed tomography (CT) that measures the amount of calcium buildup within the walls of the heart's arteries and is an essential predictor of coronary heart disease [CITE]. The leading cause of death in the United States is coronary heart disease, killing 659,000 people annually, and it is the third leading cause of mortality worldwide [CITE].

Agatston scoring is the most common CAC scoring technique [CITE] and is a good predictor of major adverse cardiac events (MACE) [CITE]. Although, studies have shown that a significant number of patients have been characterized as having no calcium (CAC=0) while still developing MACE [CITE]. This is possibly due in part to the intensity thresholding requirements associated with the Agatston scoring technique, so other approaches like spatially weighted calcium scoring have been put forth as an alternative to Agatston scoring. Spatially weighted calcium scoring improves upon traditional calcium scoring by avoiding thresholding and has been shown to predict MACE more accurately [CITE], [CITE]. Spatially weighted calcium scoring is still limited in distinguishing non-zero CAC from noise. It also lacks quantitative insight as it is an arbitrary score without any direct physical association.

The Agatston technique can be used to estimate calcium volume and calcium mass. The calcium mass score acquired via the Agaston technique is fundamentally similar to traditional Agatston calcium scoring and suffers from many of the same limitations inherent to the Agatston scoring approach [CITE], [CITE]. A previous study also showed that calcium mass quantification via the Agatston approach produces up to 50% underestimation of calcium mass for large patients [CITE]. 

A CAC measurement technique that improves the accuracy, reproducibility, sensitivity, and specificity of the measured calcium mass would help to improve coronary heart disease diagnosis and stratify risk more accurately, especially for undetectable and nearly undetectable calcium levels. A more quantitative calcium mass calculation would likewise help improve risk stratification and diagnosis of coronary heart disease. The integrated intensity or Hounsfield unit (HU) technique recently demonstrated the ability to accurately assess coronary artery cross-sectional area, even past the visible threshold [CITE]. Similarly, volume fraction mass quantification has demonstrated the ability to accurately measure the mass of lung tissue in pulmonary computed tomography (CT) scans. These techniques have not been adapted for calcium mass quantification and CAC scoring. 


This study adapted the integrated HU and volume fraction mass quantification techniques to quantify coronary artery calcium in CT scans by accounting for the partial volume effect. These new calcium mass quantification techniques, integrated calcium mass, and volume fraction calcium mass were evaluated on simulated CT scans where the total calcium mass was known. These techniques were compared to the known calcium mass to determine which technique was the most robust and accurate compared to the gold standard Agatston scoring technique and the recently introduced spatially weighted calcium scoring technique. 

Throughout this study, Agatston scoring refers to the calcium mass score derived from the Agatston approach. Spatially weighted calcium scoring refers only to the calcium score without associated physical units. Integrated calcium mass and volume fraction calcium mass refer only to the calcium mass calculated via the integrated HU approach and volume fraction mass quantification technique, respectively. This study's calcium scoring algorithms are publicly available at https://github.com/Dale-Black/CalciumScoring.jl [13].

"""

# ╔═╡ 980ddc9c-1b2e-4aa4-8441-e49d5eb5ba49
md"""
# 2. Method
"""

# ╔═╡ 03fab3b3-98b8-408e-b857-387ea9764ff0
md"""
## 2.1 - Simulation
"""

# ╔═╡ fc429bd1-4318-4975-b004-b0b031927f3e
md"""
The simulation study was set to match the scanning parameters of the 320-slice CT scanner (Canon Aquilion One, Canon America Medical Systems, Tustin, CA) as previously reported [CITE]. The X-ray spectrum was created with an interpolating polynomial model [CITE]. The linear attenuation coefficients were created from their chemical composition [CITE]. Poisson noise was added to simulate quantum noise. The simulation did not include Compton scatter, but beam hardening was included. A 3200x2200 pixel digital phantom was designed based on an anthropomorphic thorax phantom with a size of 30x20cm (QRM-Thorax, QRM, Mӧhrendorf, Germany). To simulate different patient sizes, additional fat rings emulated by a mixture of 20% water and 80% lipid were added, which resulted in a medium-sized phantom of 35x25cm and a large-sized phantom of 40x30cm. There were nine calcification inserts within the thorax with different densities and sizes. Three calcification inserts of different diameters (1, 3, and 5 mm) and different hydroxyapatite (HA) densities were placed within each phantom. For the normal density study, the HA densities in the inserts were 200, 400, and 800 mg/cm3. For the low-density study, the densities were changed to 25, 50, and 100 mg/cm3. Each phantom also contained a 10 mm diameter calibration rod. All phantom sizes and density levels were scanned using 80, 100, 120, and 135 kV tube voltages. For small, medium, and large patient sizes, the exposure value was adjusted to 0.9, 2.0, and 5.4 mR, respectively, resulting in similar noise levels for different-sized phantoms.

Simulation materials and geometries are shown in Figure 1. Acquisition and reconstruction parameters for the simulated and physical phantoms are shown in Table 1. The calibration rods were all 10 mm in diameter. All calcium scoring measurements were repeated across each kV, patient size, calcium insert size, and calcium insert density. 

Segmenting regions of interest (ROIs) is an important step in calcium measurement. For this study, segmentations were done automatically based on previous work by van Praagh et al. [CITE] and adapted for simulated phantoms. The automatic segmentation approach effectively segments calcium inserts based on the known geometry of the simulated phantom without requiring manual intervention. This is important because much of this study is focused on the accurate measurement of calcium mass for inserts below the threshold of visual detectability. In patients, coronary artery centerline estimates can be extracted in non-contrast CT scans by automatic methods, such as atlas-based approaches [CITE].  Recently, deep learning has shown promise in automatically segmenting cardiac anatomy and has the potential to accurately segment coronary artery centerlines in non-contrast CT scans using a supervised learning approach in patient images like the OrCaScore dataset [CITE]. The coronary artery centerlines can then automatically generate ROIs for calcium mass measurement of coronary artery calcification.
"""

# ╔═╡ a05b6404-db94-4074-9ae7-e39e9da5c10d
load(joinpath(dirname(pwd()), "figures", "phantom materials.png"))

# ╔═╡ 20866c33-94d0-451a-8fb0-0b758bf2c1dc
md"""
Fig. 1 Shows a sketch of the simulated phantom with the colors highlighting the different materials in the simulated phantoms.
"""

# ╔═╡ 810d38ae-d826-4f6d-b095-e80f8ebb766e
df_sim = DataFrame(
	"Parameter" => [
		"Manufacturer",
		"CT System",
		"Reconstruction",
		"Tube Voltage (kV)",
		"Exposure Small (mR)",
		"Exposure Medium (mR)",
		"Exposure Large (mR)",
		"Slice Thickness (mm)",
		"Matrix Size Small (pixels)",
		"Matrix Size Medium (pixels)",
		"Matrix Size Large (pixels)",
		"Detector Element Width (mm)",
		"Detector Thickness (mm)",
	],
	"Simulation" => [
		"Canon",
		"Aquilion One Vision",
		"FBP",
		"80, 100, 120, 135",
		"0.9",
		"2",
		"5.4",
		"0.5",
		"640 x 440",
		"740 x 540",
		"840 x 640",
		"0.5",
		"0.5",
	],
	"Physical" => [
		"Canon"
		"Aquilion One Vision"
		"FBP"
		"120"
		"15"
		"N/A"
		"84"
		"3.0"
		"512 x 512"
		"N/A"
		"512 x 512"
		"UNSURE"
		"UNSURE"
	]
)

# ╔═╡ 606bf394-2371-4550-ae21-10b4f7f1f751
md"""
Table 1. Acquisition and reconstruction parameters for the simulated and physical phantoms.
"""

# ╔═╡ 339b7757-334d-427f-964b-594dd1e9c5b6
md"""
To study the effect that motion has on Agatston scoring, integrated calcium mass, volume fraction calcium mass, and spatially weighted calcium scoring, we simulated motion using a previously reported random motion filter [CITE: http://proceedings.mlr.press/v102/shaw19a.html] and [CITE: https://www.sciencedirect.com/science/article/pii/S0169260721003102]. We then performed repeated analysis on this simulated motion data to better understand the effect that motion has on all four calcium scoring techniques.
"""

# ╔═╡ 87739b9f-b6f3-4274-aed8-910ca66d4259
md"""
## 2.2 - Physical Phantom
"""

# ╔═╡ 41ef0933-6d79-4d29-aa3d-27335c265aad
md"""
This study also utilized an anthropomorphic thorax phantom (QRM-Thorax, QRM, Möhrendorf, Germany) with an insert containing calcium (Cardiac Calcification Insert (CCI), QRM, Möhrendorf, Germany). All images were acquired by G van Praagh et al. [CITE], and a limited subset of scans was included in this analysis. The cardiac calcification insert phantom consisted of nine calcification inserts made up of hydroxyapatite (HA). Within the cardiac calcification insert phantom, two calibration rods consisting of water-equivalent material and 200 mgHAcm-3 were also present. The calcifications had diameters and lengths of 1.0, 3.0, and 5.0 mm. Three different densities were present in the phantom for each calcification size: 200, 400, and 800 mgHAcm-3. Two different patient sizes, small and large, were included in this analysis by the addition of a fat ring. This fat ring increased the phantom size from 300 × 200 mm to 400 × 300 mm. Four scans, two small and two large patient scans were included in this analysis. Figure 2 Shows axial slice views of various phantoms included in this study. Fig. 2A shows a simulated small phantom with normal-density inserts (200, 400, 800 mg/cm^3) and a tube voltage of 120 kV. Fig. 2B shows a simulated large phantom with low-density inserts (25, 50, 100 mg/cm^3) and a tube voltage of 120 kV with simulated motion. Fig. 2C shows an axial slice of the calibration rod (200 mg/cm^3) for a simulated medium phantom with a tube voltage of 100 kV. Fig. 2D shows an axial slice view of a physical QRM Thorax Phantom with a Cardio Calcification Insert. The red arrow shows simulated beam hardening artifact, and the blue arrows show simulated streaking artifact caused by motion.

Segmenting regions of interest is an important step in calcium measurement. For this study, segmentations were done automatically based on previous work by van Praagh et al. [CITE] and adapted for the Julia programming language [CITE]. 
"""

# ╔═╡ 400b0f06-2906-4ec7-bbdb-d4c17b2edf55
load(joinpath(dirname(pwd()), "figures", "phantom slices.png"))

# ╔═╡ 3c409f13-8f26-48cd-9bb0-408f8e9b2785
md"""
Fig. 2 Shows axial slice views of various phantoms included in this study. (A) shows a simulated small phantom with normal-density inserts (200, 400, 800 mg/cm^3) and a tube voltage of 120 kV. (B) shows a simulated large phantom with low-density inserts (25, 50, 100 mg/cm^3) and a tube voltage of 120 kV with simulated motion. (C) shows an axial slice of the calibration rod (200 mg/cm^3) for a simulated medium phantom with a tube voltage of 100 kV. (D) shows an axial slice view of a physical QRM Thorax Phantom with a Cardio Calcification Insert. The red arrow shows simulated beam hardening artifact, and the blue arrows show simulated streaking artifact caused by motion.
"""

# ╔═╡ 174c4db1-352b-438d-a79b-ee51e026f704
md"""
## 2.3 - Agatston Scoring
"""

# ╔═╡ d045cdb6-7f03-4da1-8a50-8db1ccbee43c
md"""
Agatston scoring is defined at a tube voltage of 120 kV [CITE], but recent papers have shown how Agatston scoring can be adjusted for use in low-dose scans (70, 80, and 100 kV) [CITE]. For this study, we assumed an exponentially decreasing trendline and extrapolated beyond to account for a higher tube voltage of 135 kV. This results in an extrapolation formula shown in Equation 1, where ``y_{thresh}`` corresponds to the extrapolated threshold (HU) and ``TV`` corresponds to the tube voltage (kV). All kV-specific thresholds used in this study are shown in Table 2.
"""

# ╔═╡ aa0e180a-ff22-48f5-8799-e895c1208b2d
md"""

```math
\begin{equation}
y_{thresh}=(378) e^{-0.009(TV)}
\end{equation}
\tag{1}
```

"""

# ╔═╡ 5eb54efb-a2dd-447d-a4db-5cdaae8d11d7
df_agat_thresh = DataFrame(
	"Tube Voltage (kV)" => [80, 100, 120, 135],
	"Threshold (HU)" => [177, 145, 130, 112]
)

# ╔═╡ c9ead9dd-2747-483b-a98c-805fddfeb5b6
md"""
Table 2. Tube voltage adapted thresholds for Agatston scoring.
"""

# ╔═╡ 3957e508-bc7e-4415-991a-09e74e6824d4
md"""
## 2.4 - Integrated Calcium Mass
"""

# ╔═╡ f33a8610-e75a-4e82-9fa7-057f6eca70fb
md"""
As previously reported [CITE], the integrated Hounsfield technique assumes that although the partial volume effect influences the HU of a particular voxel, the total integrated HU within an ROI is conserved. This study addresses the issue of non-detectable CAC by adjusting the integrated HU technique for use in CAC scoring, calling it integrated calcium mass. The cross-sectional area equation (Eq. 2) can be adjusted for use in three-dimensional regions (Eq. 3) and applied to calcium mass quantification. Fig. 3A shows the cross-section of a simulated coronary plaque with the measurements that need to be computed for each image. ``S_{Obj}`` is the intensity (HU) of a region within the plaque containing pure calcium with no partial volume affected voxels. ``S_{Bkg}`` is the measured intensity (HU) of a ring-like section of the coronary lesion with no calcium. ``I`` is the sum of every voxel intensity (HU), including those voxels affected by the partial volume effect. ``V`` is the total number of voxels in the cross-section multiplied by the voxel size. Eq. 4 shows how to convert the volume of the object (``V_{Obj}``) to mass (``M_{Obj}``), where ``\rho_{S_{Obj}}`` is the density of the calcification, specific to the object intensity (``S_{Obj}``).

It is impractical to accurately measure the intensity of a calcification with no partial volume effect (``S_{Obj}``) for small calcifications. Therefore, this study utilized a calibration rod and computed a volume based on the measured intensity of this calibration rod ``S_{Obj}``, then converted the volume to mass by multiplying by the known density of that rod. Integrating over the entire ROI should remove any effect from noise, assuming noise only affects individual voxels.

"""

# ╔═╡ e2d7f8ed-ff9d-489e-ba8e-397ef1b70e06
md"""
```math
\begin{aligned}
I &= [(A-CSA) \times S_{Bkg} ] + [CSA \times S_{Obj}]
\\
CSA &= \frac{I - (A \times S_{Bkg})}{S_{Obj}-S_{Bkg}}
\end{aligned}
\tag{2}
```
\

```math
\begin{aligned}
I &= [(V - V_{Obj}) \times S_{Bkg} ] + [ V_{Obj} \times S_{Obj}]
\\
 V_{Obj} &= \frac{I - (V \times S_{Bkg})}{S_{Obj}-S_{Bkg}}
\end{aligned}
\tag{3}
```

\

```math
\begin{aligned}
M_{Obj} = V_{Obj} \times \rho_{S_{Obj}}
\end{aligned}
\tag{4}
```
"""

# ╔═╡ 7a489daa-64b0-441a-8091-d7493e2c4cf4
md"""
## 2.5 - Volume Fraction Mass
"""

# ╔═╡ 7ad72b43-f604-4953-9309-e7e226687f82
md"""
The volume fraction calcium mass technique is similar to integrated calcium mass, but instead of calculating the calcium within an entire ROI, the percent of calcium contained within one voxel is calculated. The percent calcium contained within each voxel is then summed up within a ROI to obtain the total percentage of calcium for that ROI (Eq. 5). Given the known size of the ROI and density of the calibration rod, the volume and mass of calcium can be calculated (Eq. 6). 

Equation 5 shows how to calculate the percentage of calcium contained within each voxel. ``i`` is the intensity of one voxel of interest (HU), ``S_{Obj}`` is the intensity of pure background, which can be obtained from a ring-like object as seen in Fig. 3B, and ``S_{Obj}`` is the intensity of pure calcium which can be obtained from a calibration rod of any density. The result, ``k_{i}``, is then the percentage of calcium within one voxel ``i``. The entire region of voxels is then summed to give ``K`` which is the total percentage of calcium contained within the ROI.

Equation 6 shows how to calculate the volume of calcium ``V_{Obj}``, given the total percentage of calcium within an ROI ``K`` and the known volume of that ROI ``V_{ROI}``. This can then be converted into a calcium mass ``M_{Obj}`` given the density of the calibration rod. Fig. 3 shows the parameters needed to calculate mass via the integrated calcium mass technique (Fig. 3A) and the volume fraction calcium mass technique (Fig. 3B).
"""

# ╔═╡ e7bfebba-6597-409a-a89c-a81f9bd5706c
md"""
```math
\begin{aligned}
k_{i} &= \frac{i - S_{Bkg}}{S_{Obj} - S_{Bkg}} \\
K &= \sum_{i} k_{i}
\end{aligned}
\tag{5}
```

\

```math
\begin{aligned}
V_{Obj} &= K \times V_{ROI} \\
M_{Obj} &= V_{Obj} \times \rho_{S_{Obj}}
\end{aligned}
\tag{6}
```
"""

# ╔═╡ dc703cec-e97e-4600-994e-48cea5a67392
load(joinpath(dirname(pwd()), "figures", "integrated-volume fraction methods.png"))

# ╔═╡ d822ff70-0c79-48dd-a15f-9b29078f6f2c
md"""
Fig 3. Shows two identical simulated vessel lumens with ROIs for calcium measurement.(A) Shows the ROIs needed for the integrated calcium mass technique. The central ROI that yields ``S_{Obj}`` and the background ROI that yields ``S_{Bkg}`` are unaffected by the partial volume effect, while the object ROI used to calculate ``I`` is affected by the partial volume effect. (B) Shows the ROIs needed for the volume fraction technique and a zoomed-in portion of the simulated vessel lumen, where ``i`` is the individual voxel intensity. ``S_{Bkg}`` is a ring-like background ROI unaffected by the partial volume effect. ``S_{Obj}`` is the mean intensity of known calcium density obtained from a calibration rod unaffected by the partial volume effect.
"""

# ╔═╡ 64078475-d2b7-43b5-8396-fee783d17389
md"""
## 2.6 - Spatially Weighted Calcium Scoring
"""

# ╔═╡ 44bccf94-4896-407a-a7ce-422d019ec14e
md"""
This study reimplemented the previously described spatially weighted calcium scoring technique [CITE] to compare integrated calcium mass and volume fraction calcium mass against a more recently proposed calcium scoring approach. The same set of voxels included in the integrated calcium mass and volume fraction calcium mass regions of interest were included in the spatially weighted calcium scoring calculations. Each voxel was assigned a weight using a normal distribution function with means and standard deviations obtained from 100 mg/cm3 calcium rod measurements across multiple simulated scans. These measurements are shown in Table 3. Only small patient sizes were used for measurements, and the signal-to-noise ratio was kept constant between different energies.
"""

# ╔═╡ f7a9f808-00e0-4c01-86e6-213012c5a14c
df_swcs_cal = DataFrame(
	"kV" => [80, 100, 120, 135],
	"Mean" => [201.490, 174.366, 158.264, 152.881],
	"Std" => [34.304, 28.171, 22.966, 21.078]
)

# ╔═╡ 2dceaac8-39c6-471c-bc40-ecf1082199d5
md"""
Table 3. Means and standard deviations of the 100 mg/cm^3 calcium calibration rods.
"""

# ╔═╡ 377b892d-6b63-4b5a-bc2d-bf3490285f24
md"""
## 2.7 - Statistical Analysis
"""

# ╔═╡ 75c3bc65-717e-4164-8247-3757088ce28a
md"""
All calcium scoring calculations and analyses were performed using the programming language Julia [CITE]. Root-mean-square error (RMSE) and root-mean-square deviation (RMSD) were calculated for all linear regression measurements to test for accuracy (RMSE) and precision (RMSD). Equation 7 shows how to calculate RMSE and RMSD. ``N`` is the total number of data points, ``\hat{y}`` is the calculated calcium masses, and ``y`` is either the ground truth calcium masses (RMSE) or the linear regression-based calcium masses (RMSD), which is computed based on the calculated calcium masses.

"""

# ╔═╡ 81e50a3f-6996-4d84-a320-dbb05251514f
md"""
```math
\begin{equation}
RMS = \sqrt{\frac{\sum{|y-{\hat{y}|}^2}}{N}}
\end{equation}
\tag{7}
```
"""

# ╔═╡ 3eb8006c-e5cb-426a-8083-aa5510761050
md"""
# 3. Results
"""

# ╔═╡ 51fc6aac-6ba0-417b-8825-8ce173dcbb74
md"""
## 3.1 - Simulated Phantom
"""

# ╔═╡ ae521a11-98ed-44f0-9feb-921a8a2f742e
md"""
### 3.1.1 - Accuracy
"""

# ╔═╡ 9bd0f11a-908b-4063-ad2c-3dc27835af65
md"""
Agatson scoring, integrated calcium mass, and volume fraction calcium mass can calculate mass directly, making the comparison between known mass and those techniques straightforward. The spatially weighted approach produces an arbitrary score correlated with the Agatston score [CITE]. A mass calculation derived from spatially weighted calcium scoring is not practical due to the lack of thresholding in the spatially weighted technique. 

Linear regression was performed for integrated calcium mass, volume fraction calcium mass, and Agatston scoring against known calcium mass on all stationary and motion-affected simulated phantoms. The low-density and normal-density phantoms were examined separately. 

The correlation coefficient (r), RMSE, and RMSD values were calculated for Agatston scoring, integrated calcium mass, and volume fraction calcium mass. Linear regression used the known calcium mass as the reference. Integrated calcium mass on normal-density phantoms produced an r-correlation coefficient, RMSE, and RMSD value of 1.000, 0.677 mg, and 0.602 mg, respectively. Volume fraction calcium mass on normal-density phantoms produced an r-correlation coefficient, RMSE, and RMSD value of 1.000, 0.609 mg, and 0.599 mg, respectively. Agatston scoring on normal-density phantoms produced an r-correlation coefficient, RMSE, and RMSD value of 1.000, 1.650 mg, and 0.906 mg, respectively. All normal-density accuracy measurements are shown in Table 4.

Integrated calcium mass on low-density phantoms (Fig. 4A) produced an r-correlation coefficient, RMSE, and RMSD value of 0.991, 0.495 mg, and 0.494 mg, respectively. Volume fraction calcium mass on low-density phantoms (Fig. 4B) produced an r-correlation coefficient, RMSE, and RMSD value of 0.989, 0.585 mg, and 0.575 mg, respectively. Agatston scoring on low-density phantoms (Fig. 4C) produced an r-correlation coefficient, RMSE, and RMSD value of 0.789, 3.509 mg, and 2.240 mg, respectively. 

A similar analysis was performed for all motion-affected phantoms, and the results can be seen in Table 4 The trend continued to show that integrated calcium mass and volume fraction calcium mass outperformed Agatston scoring for both low-density and normal-density phantoms.
"""

# ╔═╡ 76e59ea9-445e-4db4-9d3c-232a7aa57a85
load(joinpath(dirname(pwd()), "figures", "stationary", "accuracy_low.png"))

# ╔═╡ fbb70982-2828-4383-889f-fb195f9f9b94
md"""
Fig. 4 Shows the linear regression analysis comparing measured calcium to the known calcium for the low-density (25, 50, 100 mg/ cm3) stationary phantoms. Every tube voltage (80, 100, 120, 135 kV) and size (small, medium, large) is included in the analysis. (A) shows the results of integrated calcium mass. (B) shows the results of the volume fraction method. (C) shows the results of Agatston mass scoring. The best fit line, along with the root mean squared error and root mean squared deviation values are shown in each plot. 
"""

# ╔═╡ 233a15d1-2e71-4a93-ab5f-c9599a43baf2
md"""
### 3.1.2 - Reproducibility
"""

# ╔═╡ a906e0da-6f16-494c-8325-6a1bf38c2c2d
md"""
This study included another set of simulated images with identical geometry to understand how repeatable the four different calcium scoring techniques are. The only variation between the two groups of images comes from the random quantum noise associated with the simulation itself. First, all false-negative values were removed from the analysis and then the reproducibility of the results was calculated. Then the results for the stationary phantoms were plotted, comparing measurements from the first set of images against the second set of images for all four scoring techniques (Fig. 5). The RMSE and RMSD values for integrated calcium mass were 0.708 mg and 0.692 mg, respectively (Fig. 5A). The RMSE and RMSD values for volume fraction calcium mass were 0.783 mg and 0.776 mg, respectively (Fig. 5B). The RMSE and RMSD values for Agatston scoring were 0.982 mg and 0.980 mg, respectively (Fig. 5C). The RMSE and RMSD values for spatially weighted calcium scoring were 4.722 mg and 4.714 mg, respectively (Fig. 5D). 

The correlation coefficient and best-fit line were also calulated for each calcium scoring technique and are shown in Table 4. The results were then repeated for the motion-affected simulated phantom data. Again, similar trends were seen in reproducibility analysis for both the stationary scans and motion-affected scans. These results are shown in Table 4.
"""

# ╔═╡ dbd35cc5-6528-4056-b153-fc34b77ba76c
load(joinpath(dirname(pwd()), "figures", "stationary", "reproducibility.png"))

# ╔═╡ 20cc068b-e06b-4679-977c-3e0e1faf492f
md"""
Fig. 5 Shows reproducibility measurements for all four scoring techniques on the stationary phantoms. Every tube voltage (80, 100, 120, 135 kV), size (small, medium, large), and density (low, normal) is included in the analysis. The measurements from the first set of images were plotted against the second set of images for each technique. Integrated mass (A), Volume fraction (B), Agatston mass scoring (C), and spatially weighted calcium scoring (D) are shown along with the root mean squared error and root mean squared deviation values.
"""

# ╔═╡ 13573036-576c-4214-b1ba-b0f3ff764b48
md"""
### 3.1.3 - Sensitivity and Specificity
"""

# ╔═╡ d58aa466-3044-4608-a48a-d64fd46549d3
md"""
The percentage of false-negative (CAC=0) and false-positive (CAC>0) scores was also calculated to understand the sensitivity and specificity of all four calcium scoring techniques. Any region containing known calcium that resulted in a CAC score of zero was determined to be a false-negative (CAC=0) score, and any region of pure background that resulted in a positive calcium score was determined to be a false-positive (CAC>0) score. This was simple to obtain for Agatston scoring, as Agatston returns zero as a possible score. 

For spatially weighted calcium scoring, the score is always greater than zero, making a score of zero impossible. This is unhelpful since noise will eventually become a dominating factor for low-density and small calcifications making it difficult to distinguish between a small amount of calcium and noise, so spatially weighted calcium scores for regions of pure background were calculated across all simulated phantoms, and the mean and standard deviation of these scores were calculated. Any spatially weighted calcium score containing regions of known calcium that resulted in a score less than the mean background score plus one standard deviation was then considered a false-negative score. Likewise, any spatially weighted calcium score greater than the mean score of pure background plus one standard deviation was considered a false-positive (CAC>0) score for regions containing only background.

A similar approach for determining false-negative (CAC=0) and false-positive (CAC>0) scores was applied to integrated calcium mass and volume fraction calcium mass. Both calcium mass techniques result in real masses from ``(-\infty,  \infty)``. For low-density and small calcifications, the noise will dominate the calcium mass results making it difficult to distinguish between a small amount of calcium and noise, so a mean and standard deviation mass was obtained in the same way as spatially weighted calcium scoring. Any integrated calcium mass or volume fraction calcium mass result containing regions of known calcium that resulted in a mass less than the mean background mass plus one standard deviation was then considered a false-negative score. Likewise, integrated calcium mass or volume fraction calcium mass result greater than the mean mass of pure background plus one standard deviation was considered a false-positive (CAC>0) mass for regions containing only background.

The percentage of false-negative (CAC=0) and false-positive (CAC>0) scores were then calculated for all four techniques. For the stationary simulated phantoms, integrated calcium mass and volume fraction calcium mass techniques produced the fewest false-negative (CAC=0) results (21 out of 216 scores), and spatially weighted calcium scoring produced the most (92 out of 216 scores). While Agatson scoring produced a total of 84 false-negative (CAC=0) out of 216 scores. The integrated calcium mass and volume fraction calcium mass techniques produced 8 false-positive (CAC>0) results out of 72 scores. Spatially weighted calcium scoring produced 12 false-positive (CAC>0) out of 72 total scores, and Agatston scoring produced zero false-positive scores (CAC>0). These results are shown in Fig. 6. Results for the sensitivity and specificity of each technique, evaluated on the motion-affected phantoms, are shown in Table 4.

"""

# ╔═╡ 7639ce2b-48e7-45aa-80c9-de20cd486714
load(joinpath(dirname(pwd()), "figures", "stationary", "sensitivity_specificity.png"))

# ╔═╡ bae2f0a9-3819-421a-9660-775ea2e43d42
md"""
Fig 6. Shows the percentage of false-negative (CAC=0) and false-positive (CAC>0) scores on the stationary phantoms. Every tube voltage (80, 100, 120, 135 kV), size (small, medium, large), and density (low, normal) is included in the analysis.
"""

# ╔═╡ e17c61ba-91a0-4dbb-9f5a-2182280aec0b
md"""
## 3.2 - Physical Phantom
"""

# ╔═╡ 0ad389de-3aed-41a0-9680-430ddfac54ad
md"""
The physical phantom scans acquired by Praagh et al. [CITE] were also used for analysis. Two small and two large phantom scans acquired by the CANON scanner (Table 1) were included in this analysis to validate the simulation results. Spatially weighted calcium scoring requires a collection of 100 mg/cm^3 calibration rod measurements for the weighting function. The physical QRM phantom scans acquired by Praagh et al. [CITE] only contained a 200 mg/cm^3 calibration rod. Therefore, spatially weighted calcium scoring was excluded from the analysis of the physical phantom.
"""

# ╔═╡ 3a80288d-cf55-4f2e-87e3-c456bd483954
md"""
### 3.2.1 - Accuracy
"""

# ╔═╡ 3471dcb8-f2c7-4215-bb9f-66febdb4d28d
md"""
Linear regression was performed for integrated calcium mass, volume fraction calcium mass, and Agatston scoring against known calcium mass on all four physical phantom scans. Only normal-density (200, 400, 800 mg/cm^3) is included in the physical QRM phantom.

The correlation coefficient (r), RMSE, and RMSD values were calculated for Agatston scoring, integrated calcium mass, and volume fraction calcium mass. Linear regression used the known calcium mass as the reference. Integrated calcium mass on the normal-density physical phantom (Fig. 7A) produced an r-correlation coefficient, RMSE, and RMSD value of 0.981, 8.811 mg, and 5.752 mg, respectively. Volume fraction calcium mass on the normal-density physical phantom (Fig. 7B) produced an r-correlation coefficient, RMSE, and RMSD value of 0.975, 8.559 mg, and 6.647 mg, respectively. Agatston scoring on the normal-density physical phantom (Fig. 7C) produced an r-correlation coefficient, RMSE, and RMSD value of 0.973, 22.485 mg, and 4.653 mg, respectively.

All normal-density physical phantom accuracy measurements, including the best-fit line and the r-correlation coefficient, are shown in Table 4.
"""

# ╔═╡ 20ca113e-c67f-46a5-8dae-dc00434be100
load(joinpath(dirname(pwd()), "figures", "physical", "accuracy.png"))

# ╔═╡ 2393e3f9-5431-4d27-8337-f412e411d737
md"""
Fig. 7 Shows the linear regression analysis comparing measured calcium to the known calcium for the normal density (200, 400, 800 mg/ cm3) physical phantom scans at a tube voltage of 120 kV. (A) shows the results of integrated calcium mass. (B)shows the results of the volume fraction method. (C) shows the results of Agatston mass scoring. The best fit line, along with the root mean squared error and root mean squared deviation values, are shown in each plot. 
"""

# ╔═╡ 6ff085e7-daa8-45a5-b96c-7cd1863509c3
md"""
### 3.2.2 - Reproducibility
"""

# ╔═╡ 30f8b807-26ed-48e2-9009-64db01467194
md"""
To analyze the reproducibility of integrated calcium mass, volume fraction calcium mass, and Agatston scoring on the physical phantom, the first small phantom, and large phantom scans were compared to the second small phantom and large phantom scans. These scans were acquired under identical settings, with the only variation coming from a rotation of the phantom. All false-negative values were removed from the analysis, and then the reproducibility of each technique was calculated. 

The RMSE and RMSD values for integrated calcium mass were 2.838 mg and 1.906 mg, respectively (Fig. 8A). The RMSE and RMSD values for volume fraction calcium mass were 3.203 mg and 1.777 mg, respectively (Fig. 8B). The RMSE and RMSD values for Agatston scoring were 2.144 mg and 1.792 mg, respectively (Fig. 8C).

The up to 50% underestimation for the Agatston scoring masses (Fig. 7C) means the calculated masses are, on average, smaller than the calculated masses from the integrated and volume fraction calcium mass techniques. Since RMS values are sensitive to the magnitude of the inputs, the smaller RMSE and RMSD values for Agatston scoring, relative to integrated calcium mass and volume fraction calcium mass techniques, are explained by the underestimated mass calculations of the Agatston technique.

The correlation coefficient and best-fit line were also calculated for each calcium scoring technique and are shown in Table 4.
"""

# ╔═╡ a184dfb4-1a2d-4b38-bb9f-937fcfea9e62
load(joinpath(dirname(pwd()), "figures", "physical", "reproducibility.png"))

# ╔═╡ 60657f1d-8c4a-41a5-961e-5425a5bfb9b8
md"""
Fig. 8 Shows reproducibility measurements on the physical phantom scans at a tube voltage of 120 kV. Measurements from two scans were plotted against a second set of measurements from two different scans at different angles. Integrated mass (A), Volume fraction (B), Agatston mass scoring (C) are shown along with the root mean squared error and root mean squared deviation values.
"""

# ╔═╡ daf56f50-e577-48ec-80a5-33e7fbeb2a38
md"""
### 3.2.3 - Sensitivity and Specificity
"""

# ╔═╡ dcdf781f-18c5-406f-8b71-dc21aa462ba3
md"""
The approach used on the simulated phantoms that determined the cutoff value for false-negative (CAC=0) and false-positive (CAC>0) scores, based on the mean and standard deviation of pure background, was also utilized for physical phantom sensitivity and specificity analysis.

For the physical phantom scans, integrated calcium mass and volume fraction calcium mass techniques produced six false-negative (CAC=0) scores out of 36 total scores. While Agatson scoring produced a total of 7 false-negative (CAC=0) out of 216 scores. 

The integrated calcium mass and volume fraction calcium mass techniques produced one false-positive (CAC>0) result out of 12 scores. Agatston scoring produced zero false-positive scores (CAC>0). These results are shown in Fig. 9. and are summarized in Table 4.
"""

# ╔═╡ a7c9ce18-1bf8-45af-8298-2099044c8a4f
load(joinpath(dirname(pwd()), "figures", "physical", "sensitivity_specificity.png"))

# ╔═╡ 72de39a4-66e2-4180-a5f3-85fae564a17c
md"""
Fig 9. Shows the percentage of false-negative (CAC=0) and false-positive (CAC>0) scores for the physical phantom scans at a tube voltage of 120 kV.
"""

# ╔═╡ 0612338d-ca4c-4504-894f-7b47cfcb0949
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

# ╔═╡ 01e48aa6-5606-44ad-b624-dbbb84c39295
md"""
Table 4. Summary of all the accuracy, reproducibily, sensitivity, and specificity measurements across every type of phantom (stationary, motion-affected, and physical phantom)
"""

# ╔═╡ 18d5e6e7-3efb-4afc-87a6-e6721679b951
md"""
# 4. Discussion
"""

# ╔═╡ 4833c0c6-d751-4ab8-949d-86c2d17547ec
md"""
Calcium mass was measured for different patient sizes, kVps, calcium sizes, and calcium densities, with and without motion on simulated phantoms, and validated in a subset of physical phantom scans acquired by Praagh et al. and shared with us for analysis. The calculated calcium mass was compared against the known mass of the calcium inserts.

The results indicate that integrated calcium mass and volume fraction calcium mass are more accurate, reproducible, sensitive, and specific than Agatston scoring (Table 4).

Agatston scoring has commonly been used in the past for predicting patient outcomes. However, a limitation of Agatston scoring is that it's only defined at 120 kVp, and a threshold of 130 HU is commonly used for calcium detection. However, the calcium attenuation coefficient is energy dependent, which makes scoring challenging when images are acquired at lower kVps, to reduce patient radiation dose. Recent reports have introduced correction factors for calcium measurements at lower kVps [CITE]. Another limitation of a thresholding approach for calcium measurement is that it is affected by partial volume effect and motion. We have introduced two new methods for calcium mass quantifications based on integrated intensity (HU) and volume fraction that address the above limitations.

Previous studies have shown that up to 5% of patients with suspected zero CAC (CAC=0) will experience MACE, despite no detectible calcium by Agatston scoring [CITE]. One question arises as to whether these patients had calcium that is not detectible by traditional Agatston scoring or simply no calcium. Integrated calcium mass, volume fraction calcium mass, and spatially weighted calcium scoring attempt to address this concern by removing the intensity thresholding requirements of standard Agatston scoring. This study shows that integrated calcium mass and volume fraction calcium mass are more sensitive to low-density calcifications than spatially weighted calcium scoring and Agatston scoring. The results showed that the percentage of false-negative (CAC=0) scores on the stationary simulated phantom were 11.111, 11.111, 42.593, and 38.889 for integrated calcium mass, volume fraction calcium mass, spatially weighted calcium scoring, and Agatston scoring, respectively (Figure 6). The substantial reduction in false-negative zero calcium scores for the integrated calcium mass and volume fraction calcium mass techniques compared to the existing techniques will help address the current limitation for patients with false-negative (CAC=0) scores.

The integrated calcium mass and volume fraction calcium mass techniques can detect calcifications that are currently indetectable by the Agatston scoring approach due to its thresholding requirement. Furthermore, a previous study has shown that calcium volume was positively and independently associated with major adverse cardiac event risk, and calcium density was inversely associated with major adverse cardiac event risk [CITE]. Another study has shown that calcium density score was the strongest positive independent predictor of major adverse cardiac events, compared to Agatston score, mass score, and volume score [CITE]. Disagreements between these studies are possibly related to the thresholding approach of Agatston scoring and poor reproducibility of Agatston scoring, which is also a limitation of all the traditional calcium (mass, volume, and density) scoring approaches based on the Agatston technique. Integrated calcium mass and volume fraction calcium mass provide more accurate, reproducible, and quantitative approaches to calcium measurement. Future studies on patient data comparing Agatston scoring with integrated calcium mass and volume fraction calcium mass might help explain these seemingly contradictory results better

When acquiring the mass of calcium in the QRM phantom using Agatston scoring, it has been shown the mass scores can be significantly different from the known physical mass [CITE]. Our results (Figures 4 and 7) show that integrated calcium mass and volume fraction calcium mass measure calcium mass more accurately than Agatston scoring, on both simulated and physical phantoms. We also see an improvement in reproducibility when comparing integrated calcium mass and volume fraction calcium mass against Agatston scoring and spatially weighted calcium scoring (Fig. 5). 

Although many parameters were considered in the simulation, only one type of reconstruction was utilized, which prevented us from addressing the effect that reconstruction type has on these calcium scoring techniques. Previous studies have shown good agreement between Agatston scores based on filtered back-projection, hybrid iterative reconstruction, and deep learning-based reconstruction [CITE]. 

Slice thickness plays an important role in calcium scoring, and traditional Agatston scoring is only defined at a slice thickness of 3 mm. Recent studies have shown that the accuracy and sensitivity of Agatston scoring are improved when slice thickness is decreased [CITE]. Our simulation was limited to 0.5 mm slice thickness which is expected to provide more accurate and sensitive comparisons for Agatston scoring. Nonetheless, future studies might provide insights by varying the slice thickness. This study was also limited by the lack of realistic cardiac hardware, such as stents, that are known to cause blooming or motion artifacts [CITE]. Future studies should account for this by including physical phantoms with realistic cardiac hardware or a robust patient data set. 

Fewer total scans were included in the physical phantom analysis compared to the simulated analysis, resulting in fewer calcium scores. This is a limitation and becomes more pronounced in the sensitivity and specificity analysis, where false-negative (CAC=0) and false-positive (CAC>0) results are rare, to begin with. In addition, motion was not included in the physical phantom analysis. Likewise, the simulated motion is based on motion in magnetic resonance imaging. More physical phantom scans, including scans affected by motion, need to be included in future studies. This would provide more robust accuracy, reproducibility, sensitivity, and specificity measurements and would result in a more reliable expected percentage of false-negative (CAC=0) and false-positive (CAC>0) scores for each calcium quantification technique.

Previous studies show that Agatston scoring consistently underestimates calcium density and volume, with even further underestimation for low-density and motion-affected plaques [CITE]. Werf et el. indicates that low-density calcifications might fall below the 130 HU threshold because of blurring from motion which artificially reduces the Agatston score [CITE]. Our study is consistent with these results (Fig. 4C and 7C); we showed that Agatston scoring produced the most false-negative (CAC=0) classifications for the simulated data and the physical phantom scans. Future studies are warranted in physical phantoms with lower-density calcification inserts (< 200 mg/cm3) to understand how integrated calcium mass and volume fraction calcium mass compares to Agatston mass within the low-density regime on physical data. Very high coronary artery calcium density (> 1000 mg/cm3) is quite rare, and Agatston scoring has already been shown to be a good predictor of cardiovascular disease within this subset of patients [CITE]. 

Tzolos et al. showed that Agatston scoring struggles to detect small calcifications in the coronary arteries of patients due to the threshold requirement of 130 HU and the minimum connected component requirement of 1 mm [33]. Our results are consistent with this study and indicate that the size of the calcium insert is a critical variable in accounting for false-negative (CAC=0) Agatston scores. The small inserts resulted in 40 false-negative (CAC=0) scores out of 216 total scores, whereas the medium and large inserts accounted for only 24 and 20 false-negative (CAC=0) scores, respectively. Density was also an important factor. Low-density calcifications resulted in 40 false-negative (CAC=0) scores out of 216 total scores, whereas the medium-density and high-density calcifications resulted in 32 and 12 false-negative (CAC=0) scores, respectively. Based on our results, integrated calcium mass and volume fraction calcium mass improved upon many of the issues associated with Agatston scoring and resulted in fewer false-negative (CAC=0) scores.
"""

# ╔═╡ a78dc4b4-5feb-445b-9738-86ea383d8d47
md"""
# 5. Conclusion
"""

# ╔═╡ 0d3b8a40-ecf1-4e89-9bd4-57840bef7625
md"""
This study shows that both the integrated calcium mass and volume fraction calcium mass techniques improve sensitivity, accuracy, and reproducibility of calcium mass measurement as compared to traditional Agatston scoring. This improvement in calcium scoring can potentially improve risk stratification for patients undergoing calcium scoring and further improve outcomes compared to Agatston scoring.
"""

# ╔═╡ 7285058e-0e48-43d8-8fde-45670bab2c43
md"""
# Acknowledgements
"""

# ╔═╡ f20e768c-5a93-4055-99ee-40a045a6a188
md"""
# References
"""

# ╔═╡ Cell order:
# ╠═34651281-01a4-46cb-b094-5d0582519c30
# ╟─ed03ef4a-e103-40c2-8a83-cbd565f9a35c
# ╟─29fdd303-5953-47c0-a789-6a8eaaa2dfc7
# ╟─4df8f7e4-b5dc-4c97-aa62-032137e83ec4
# ╟─11cff734-d9b3-48c5-b0c8-3d5ac163962c
# ╟─c73b73b7-973b-487d-80a1-8dce586cc3c8
# ╟─980ddc9c-1b2e-4aa4-8441-e49d5eb5ba49
# ╟─03fab3b3-98b8-408e-b857-387ea9764ff0
# ╟─fc429bd1-4318-4975-b004-b0b031927f3e
# ╟─a05b6404-db94-4074-9ae7-e39e9da5c10d
# ╟─20866c33-94d0-451a-8fb0-0b758bf2c1dc
# ╟─810d38ae-d826-4f6d-b095-e80f8ebb766e
# ╟─606bf394-2371-4550-ae21-10b4f7f1f751
# ╟─339b7757-334d-427f-964b-594dd1e9c5b6
# ╟─87739b9f-b6f3-4274-aed8-910ca66d4259
# ╟─41ef0933-6d79-4d29-aa3d-27335c265aad
# ╟─400b0f06-2906-4ec7-bbdb-d4c17b2edf55
# ╟─3c409f13-8f26-48cd-9bb0-408f8e9b2785
# ╟─174c4db1-352b-438d-a79b-ee51e026f704
# ╟─d045cdb6-7f03-4da1-8a50-8db1ccbee43c
# ╟─aa0e180a-ff22-48f5-8799-e895c1208b2d
# ╟─5eb54efb-a2dd-447d-a4db-5cdaae8d11d7
# ╟─c9ead9dd-2747-483b-a98c-805fddfeb5b6
# ╟─3957e508-bc7e-4415-991a-09e74e6824d4
# ╟─f33a8610-e75a-4e82-9fa7-057f6eca70fb
# ╟─e2d7f8ed-ff9d-489e-ba8e-397ef1b70e06
# ╟─7a489daa-64b0-441a-8091-d7493e2c4cf4
# ╟─7ad72b43-f604-4953-9309-e7e226687f82
# ╟─e7bfebba-6597-409a-a89c-a81f9bd5706c
# ╟─dc703cec-e97e-4600-994e-48cea5a67392
# ╟─d822ff70-0c79-48dd-a15f-9b29078f6f2c
# ╟─64078475-d2b7-43b5-8396-fee783d17389
# ╟─44bccf94-4896-407a-a7ce-422d019ec14e
# ╟─f7a9f808-00e0-4c01-86e6-213012c5a14c
# ╟─2dceaac8-39c6-471c-bc40-ecf1082199d5
# ╟─377b892d-6b63-4b5a-bc2d-bf3490285f24
# ╟─75c3bc65-717e-4164-8247-3757088ce28a
# ╟─81e50a3f-6996-4d84-a320-dbb05251514f
# ╟─3eb8006c-e5cb-426a-8083-aa5510761050
# ╟─51fc6aac-6ba0-417b-8825-8ce173dcbb74
# ╟─ae521a11-98ed-44f0-9feb-921a8a2f742e
# ╟─9bd0f11a-908b-4063-ad2c-3dc27835af65
# ╟─76e59ea9-445e-4db4-9d3c-232a7aa57a85
# ╟─fbb70982-2828-4383-889f-fb195f9f9b94
# ╟─233a15d1-2e71-4a93-ab5f-c9599a43baf2
# ╟─a906e0da-6f16-494c-8325-6a1bf38c2c2d
# ╟─dbd35cc5-6528-4056-b153-fc34b77ba76c
# ╟─20cc068b-e06b-4679-977c-3e0e1faf492f
# ╟─13573036-576c-4214-b1ba-b0f3ff764b48
# ╟─d58aa466-3044-4608-a48a-d64fd46549d3
# ╟─7639ce2b-48e7-45aa-80c9-de20cd486714
# ╟─bae2f0a9-3819-421a-9660-775ea2e43d42
# ╟─e17c61ba-91a0-4dbb-9f5a-2182280aec0b
# ╟─0ad389de-3aed-41a0-9680-430ddfac54ad
# ╟─3a80288d-cf55-4f2e-87e3-c456bd483954
# ╟─3471dcb8-f2c7-4215-bb9f-66febdb4d28d
# ╟─20ca113e-c67f-46a5-8dae-dc00434be100
# ╟─2393e3f9-5431-4d27-8337-f412e411d737
# ╟─6ff085e7-daa8-45a5-b96c-7cd1863509c3
# ╟─30f8b807-26ed-48e2-9009-64db01467194
# ╟─a184dfb4-1a2d-4b38-bb9f-937fcfea9e62
# ╟─60657f1d-8c4a-41a5-961e-5425a5bfb9b8
# ╟─daf56f50-e577-48ec-80a5-33e7fbeb2a38
# ╟─dcdf781f-18c5-406f-8b71-dc21aa462ba3
# ╟─a7c9ce18-1bf8-45af-8298-2099044c8a4f
# ╟─72de39a4-66e2-4180-a5f3-85fae564a17c
# ╟─0612338d-ca4c-4504-894f-7b47cfcb0949
# ╟─01e48aa6-5606-44ad-b624-dbbb84c39295
# ╟─18d5e6e7-3efb-4afc-87a6-e6721679b951
# ╟─4833c0c6-d751-4ab8-949d-86c2d17547ec
# ╟─a78dc4b4-5feb-445b-9738-86ea383d8d47
# ╟─0d3b8a40-ecf1-4e89-9bd4-57840bef7625
# ╟─7285058e-0e48-43d8-8fde-45670bab2c43
# ╟─f20e768c-5a93-4055-99ee-40a045a6a188

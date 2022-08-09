### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ e2b53960-0e99-11ed-27bd-5f46a57ba02f
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
        Pkg.add(; url="https://github.com/JuliaHealth/DICOM.jl")
        Pkg.add(; url="https://github.com/Dale-Black/DICOMUtils.jl")
        Pkg.add(; url="https://github.com/Dale-Black/PhantomSegmentation.jl")
        Pkg.add(; url="https://github.com/Dale-Black/CalciumScoring.jl")
        Pkg.add("CairoMakie")
        Pkg.add("HypothesisTests")
        Pkg.add("Colors")
        Pkg.add("MLJBase")
        Pkg.add("AlgebraOfGraphics")
    end

    using PlutoUI
    using Statistics
    using StatsBase: quantile!, rmsd
    using ImageMorphology
    using ImageFiltering
    using CSV
    using DataFrames
    using GLM
    using DICOM
    using DICOMUtils
    using PhantomSegmentation
    using CalciumScoring
    using CairoMakie
    using HypothesisTests
    using Colors
    using MLJBase
    using AlgebraOfGraphics
end

# ╔═╡ b39114c0-182e-44dc-b00b-906c81921ebd
TableOfContents()

# ╔═╡ 3389da26-b3b3-455c-a382-6bb44437cd05
md"""
## TODO
- Assess calibrations (1pt, 3pt, 6pt) for simulated data
- Add motion reproducibility (SWCS, ICS, AS) to simulated data
"""

# ╔═╡ 88d08c61-2f3b-475f-b392-031cffbea40b
md"""
## Integrated
"""

# ╔═╡ a8eec50e-2993-4bfe-94da-b943a919e6e4
path_integrated = "/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/output_new/integrated_scoring";

# ╔═╡ 76df194b-31bc-4db6-ba99-f6a35b3c0a65
df_i = CSV.read(string(path_integrated, "/full2.csv"), DataFrame)

# ╔═╡ 450c6c6c-96ae-42cf-bc6d-5c80e99a2fa7
df_i_low, df_i_normal = groupby(df_i, :DENSITY);

# ╔═╡ dde42a35-ef12-45d7-a0ae-ede776fb555a
df_i_low_small, df_i_low_medium, df_i_low_large = groupby(df_i_low, :SIZE);

# ╔═╡ 97e09092-c961-47de-95c0-91e791a938a0
df_i_normal_small, df_i_normal_medium, df_i_normal_large = groupby(df_i_normal, :SIZE);

# ╔═╡ b14db409-696e-44e2-a4bc-afe263c3647d
df_i

# ╔═╡ 8b77c804-5f9a-4919-a55f-5ab7762d291e
md"""
## Agatston
"""

# ╔═╡ 69468127-7c21-4d34-9f97-f59cd301ea3d
path_agat = "/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/output_new/agatston";

# ╔═╡ 9b5ed20b-36b6-4bfb-9cac-97a5f1a0b947
df_a = CSV.read(string(path_agat, "/full2.csv"), DataFrame);

# ╔═╡ 1ad42be8-6538-45c0-bed3-cf6f5da8b082
df_a_low, df_a_normal = groupby(df_a, :DENSITY);

# ╔═╡ d692b7b0-a12b-4b3a-a040-15634bcd252c
df_a_low_small, df_a_low_medium, df_a_low_large = groupby(df_a_low, :SIZE);

# ╔═╡ c3583be4-9add-4cf5-a738-6b63bd3a3712
df_a_normal_small, df_a_normal_medium, df_a_normal_large = groupby(df_a_normal, :SIZE);

# ╔═╡ 1d65a547-c69e-48d5-8ce2-d0827ee87adc
md"""
## Spatially Weighted
"""

# ╔═╡ d934cfb8-78ae-4fcd-a568-e5b74409d1c8
path_swcs = "/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/output_new/swcs";

# ╔═╡ e7be547b-88a0-47b2-a84b-eaa48ec7e4eb
df_s = CSV.read(string(path_swcs, "/full2.csv"), DataFrame);

# ╔═╡ 0687e798-a30c-44cf-a03c-99308a561113
df_s_low, df_s_normal = groupby(df_s, :DENSITY);

# ╔═╡ 97afb112-d5fe-4014-92de-16efb1d3cb1d
df_s_low_small, df_s_low_medium, df_s_low_large = groupby(df_s_low, :SIZE);

# ╔═╡ 535f8216-9aeb-43c9-8e51-a4a93eeda076
df_s_normal_small, df_s_normal_medium, df_s_normal_large = groupby(df_s_normal, :SIZE);

# ╔═╡ ae26e32a-dff5-4e02-b03e-bd13cf79b14b
md"""
# Figures
"""

# ╔═╡ 20b127df-a98e-4657-b5e2-dbd297a10c81
md"""
## Zero CAC Scores (FIG)
"""

# ╔═╡ cb05cac7-ecd9-44b2-8fe7-f80d272b482a
medphys_theme = Theme(;
    Axis=(
        backgroundcolor=:white,
        xgridcolor=:gray,
        xgridwidth=0.5,
        xlabelfont=:Helvetica,
        xticklabelfont=:Helvetica,
        xlabelsize=20,
        xticklabelsize=20,
        # xminorticksvisible = true,
        ygridcolor=:gray,
        ygridwidth=0.5,
        ylabelfont=:Helvetica,
        yticklabelfont=:Helvetica,
        ylabelsize=20,
        yticklabelsize=20,
        # yminortickvisible = true,
        bottomsplinecolor=:black,
        leftspinecolor=:black,
        titlefont=:Helvetica,
        titlesize=30,
    ),
)

# ╔═╡ 8ae299c5-f2a2-4ed4-bd57-45cc065ade7a
md"""
#### SWCS
"""

# ╔═╡ 2865797a-c029-4c7a-80de-e717b33f07fe
df_s_80, df_s_100, df_s_120, df_s_135, = groupby(df_s, :scan);

# ╔═╡ 76ce94f1-3fbf-461a-8dbc-7fffd7e3fabb
begin
    mean_swcs_80, std_swcs_80 = mean(df_s_80[!, :swcs_bkg]), std(df_s_80[!, :swcs_bkg])
    mean_swcs_100, std_swcs_100 = mean(df_s_100[!, :swcs_bkg]), std(df_s_100[!, :swcs_bkg])
    mean_swcs_120, std_swcs_120 = mean(df_s_120[!, :swcs_bkg]), std(df_s_120[!, :swcs_bkg])
    mean_swcs_135, std_swcs_135 = mean(df_s_135[!, :swcs_bkg]), std(df_s_135[!, :swcs_bkg])
end;

# ╔═╡ db624048-33f3-4f4c-aa67-cfacd72777d5
df_s_80

# ╔═╡ 82fa5c00-6801-400b-b2eb-b3b0f45672ea
begin
    array_s_80 = hcat(
        df_s_80[!, :calculated_swcs_large],
        df_s_80[!, :calculated_swcs_medium],
        df_s_80[!, :calculated_swcs_small],
    )
    array_s_100 = hcat(
        df_s_100[!, :calculated_swcs_large],
        df_s_100[!, :calculated_swcs_medium],
        df_s_100[!, :calculated_swcs_small],
    )
    array_s_120 = hcat(
        df_s_120[!, :calculated_swcs_large],
        df_s_120[!, :calculated_swcs_medium],
        df_s_120[!, :calculated_swcs_small],
    )
    array_s_135 = hcat(
        df_s_135[!, :calculated_swcs_large],
        df_s_135[!, :calculated_swcs_medium],
        df_s_135[!, :calculated_swcs_small],
    )
end;

# ╔═╡ f0240c36-8b76-4116-ac55-43e099e6cd62
begin
    num_zeroCAC_80 = length(findall(x -> x < mean_swcs_80, array_s_80))
    num_zeroCAC_100 = length(findall(x -> x < mean_swcs_100, array_s_100))
    num_zeroCAC_120 = length(findall(x -> x < mean_swcs_120, array_s_120))
    num_zeroCAC_135 = length(findall(x -> x < mean_swcs_135, array_s_135))

    total_zero_s = num_zeroCAC_80 + num_zeroCAC_100 + num_zeroCAC_120 + num_zeroCAC_135
    total_cac = length(array_s_80) * 4
end;

# ╔═╡ 71626005-b2d5-4afb-ab65-6606e6f88a59
md"""
#### Agatston
"""

# ╔═╡ 632d337c-5af7-45c2-a7ba-d3a757768361
array_a = hcat(
    df_a[!, :calculated_mass_large],
    df_a[!, :calculated_mass_medium],
    df_a[!, :calculated_mass_small],
);

# ╔═╡ 9bbfd561-a51d-44af-ab06-f9afe555fc98
num_zero_a = length(findall(x -> x == 0, array_a));

# ╔═╡ 509ae5c2-4423-4a61-87b7-f119615cd19c
md"""
#### Integrated
"""

# ╔═╡ 4d00a89d-2467-4505-9d87-e4066e107a6e
df_i_80, df_i_100, df_i_120, df_i_135, = groupby(df_i, :scan);

# ╔═╡ c69f965d-196e-481a-8622-b7e988d17a2b
begin
    mean_i_80, std_i_80 = mean(df_i_80[!, :mass_bkg]), std(df_i_80[!, :mass_bkg])
    mean_i_100, std_i_100 = mean(df_i_100[!, :mass_bkg]), std(df_i_100[!, :mass_bkg])
    mean_i_120, std_i_120 = mean(df_i_120[!, :mass_bkg]), std(df_i_120[!, :mass_bkg])
    mean_i_135, std_i_135 = mean(df_i_135[!, :mass_bkg]), std(df_i_135[!, :mass_bkg])
end

# ╔═╡ 2d3d0da0-587d-4a85-9360-82c7749a1945
begin
    array_i_80 = hcat(
        df_i_80[!, :calculated_mass_large],
        df_i_80[!, :calculated_mass_medium],
        df_i_80[!, :calculated_mass_small],
    )
    array_i_100 = hcat(
        df_i_100[!, :calculated_mass_large],
        df_i_100[!, :calculated_mass_medium],
        df_i_100[!, :calculated_mass_small],
    )
    array_i_120 = hcat(
        df_i_120[!, :calculated_mass_large],
        df_i_120[!, :calculated_mass_medium],
        df_i_120[!, :calculated_mass_small],
    )
    array_i_135 = hcat(
        df_i_135[!, :calculated_mass_large],
        df_i_135[!, :calculated_mass_medium],
        df_i_135[!, :calculated_mass_small],
    )
end;

# ╔═╡ 50852830-cddf-4208-a823-2cfb86c48ec6
begin
    num_zeroCAC_80_i = length(findall(x -> x < mean_i_80, array_i_80))
    num_zeroCAC_100_i = length(findall(x -> x < mean_i_100, array_i_100))
    num_zeroCAC_120_i = length(findall(x -> x < mean_i_120, array_i_120))
    num_zeroCAC_135_i = length(findall(x -> x < mean_i_135, array_i_135))

    total_zero_i =
        num_zeroCAC_80_i + num_zeroCAC_100_i + num_zeroCAC_120_i + num_zeroCAC_135_i
end;

# ╔═╡ ba9b6361-115d-47df-808f-e67147041427
function zero_cac_plot()
    f = Figure()
    colors = Makie.wong_colors()

    ##-- TOP --##
    axtop = Axis(f[1, 1]; xticks=(1:3, ["Integrated", "Spatially Weighted", "Agatston"]))

    table = [1, 2, 3]
    heights1 = [
        (total_zero_i / total_cac) * 100,
        (total_zero_s / total_cac) * 100,
        (num_zero_a / total_cac) * 100,
    ]
    barplot!(axtop, table, heights1; color=colors[1:3], bar_labels=:y)

    axtop.title = "Zero CAC Scores"
    axtop.ylabel = "% False-negative zero CAC scores"
    ylims!(axtop; low=0, high=100)
    axtop.yticks = [0, 25, 50, 75, 100]

    save(
        "/Users/daleblack/Google Drive/Research/Papers/My Papers/cac-simulation/figures/zero_cac_motion.png",
        f,
    )
    return f
end

# ╔═╡ 10fe5219-e93b-446a-be5d-9da6201b7fbd
with_theme(medphys_theme) do
    zero_cac_plot()
end

# ╔═╡ 33fc9890-72d0-4103-b393-89aae9e216ea
total_zero_i, total_zero_s, num_zero_a

# ╔═╡ 2fe1b4a3-a2ac-40ee-b90e-1906a33a5438
total_cac

# ╔═╡ dcdc9952-0fff-420b-8f30-b4da6ab6a5f9
md"""
## Linear Regression Mass Score (FIG)
"""

# ╔═╡ 14ddfbb8-1181-4843-b446-9352415d05e9
md"""
#### Normal Density
"""

# ╔═╡ c26af14f-0cf7-4a1e-8ed8-abd3eb247d73
let
    df = df_i_normal
    gt_array = vec(
        hcat(
            df[!, :ground_truth_mass_large],
            df[!, :ground_truth_mass_medium],
            df[!, :ground_truth_mass_small],
        ),
    )
    calc_array = vec(
        hcat(
            df[!, :calculated_mass_large],
            df[!, :calculated_mass_medium],
            df[!, :calculated_mass_small],
        ),
    )
    data = DataFrame(; X=gt_array, Y=calc_array)
    global model_i_normal
    model_i_normal = lm(@formula(Y ~ X), data)
    global r2_1
    r2_1 = GLM.r2(model_i_normal)
    global rms_values1
    rms_values1 = [
        rms(data[!, :X], data[!, :Y]), rmsd(data[!, :Y], GLM.predict(model_i_normal))
    ]
end

# ╔═╡ aadb907e-dead-482f-8768-f6851c08a18c
begin
    newX1 = DataFrame(; X=collect(1:1000))
    pred_i_norm = GLM.predict(model_i_normal, newX1)
end

# ╔═╡ a1d4e6f8-b43c-40bd-a3f4-fba304871cb9
co1 = coef(model_i_normal)

# ╔═╡ 33fbcb79-10a4-454b-8114-1f0ffd0904a2
let
    df = df_a_normal
    gt_array = vec(
        hcat(
            df[!, :ground_truth_mass_large],
            df[!, :ground_truth_mass_medium],
            df[!, :ground_truth_mass_small],
        ),
    )
    calc_array = vec(
        hcat(
            df[!, :calculated_mass_large],
            df[!, :calculated_mass_medium],
            df[!, :calculated_mass_small],
        ),
    )
    data = DataFrame(; X=gt_array, Y=calc_array)
    global model_a_normal
    model_a_normal = lm(@formula(Y ~ X), data)
    global r2_3
    r2_3 = GLM.r2(model_a_normal)
    global rms_values3
    rms_values3 = [
        rms(data[!, :X], data[!, :Y]), rmsd(data[!, :Y], GLM.predict(model_a_normal))
    ]
end

# ╔═╡ 4b1ebeb0-1797-457d-b52f-0ea7d8e57318
begin
    newX3 = DataFrame(; X=collect(1:1000))
    pred_a_norm = GLM.predict(model_a_normal, newX3)
end

# ╔═╡ 1f196a66-9c0a-43d2-b41f-b302954753c2
co3 = coef(model_a_normal)

# ╔═╡ aeb864db-fba4-4d9c-8120-672339aaa90a
function lin_reg_norm()
    f = Figure()

    ##-- A --##
    ax1 = Axis(f[1, 1])

    df = df_i_normal
    scatter!(ax1, df[!, :ground_truth_mass_large], df[!, :calculated_mass_large])
    errorbars!(
        ax1,
        df[!, :ground_truth_mass_large],
        df[!, :calculated_mass_large],
        rms(df[!, :ground_truth_mass_large], df[!, :calculated_mass_large]),
    )
    scatter!(ax1, df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium])
    errorbars!(
        ax1,
        df[!, :ground_truth_mass_medium],
        df[!, :calculated_mass_medium],
        rms(df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium]),
    )
    scatter!(
        ax1, df[!, :ground_truth_mass_small], df[!, :calculated_mass_small]; color=:red
    )
    errorbars!(
        ax1,
        df[!, :ground_truth_mass_small],
        df[!, :calculated_mass_small],
        rms(df[!, :ground_truth_mass_small], df[!, :calculated_mass_small]),
    )
    lines!(ax1, [-1000, 1000], [-1000, 1000])
    lines!(ax1, collect(1:1000), pred_i_norm; linestyle=:dashdot)
    Textbox(
        f[1, 1];
        placeholder="y = $(trunc(co1[2]; digits=3))x + $(trunc(co1[1]; digits=3)) \nr = $(trunc(r2_1; digits=3)) \nRMSE: $(trunc(rms_values1[1]; digits=3)) \nRMSD: $(trunc(rms_values1[2]; digits=3))",
        tellheight=false,
        tellwidth=false,
        boxcolor=:white,
        halign=:left,
        valign=:top,
        textsize=12,
    )

    xlims!(ax1; low=0, high=200)
    ylims!(ax1; low=0, high=200)
    ax1.xticks = [0, 50, 100, 150, 200]
    ax1.yticks = [0, 50, 100, 150, 200]
    ax1.xlabel = "Known Mass (mg)"
    ax1.ylabel = "Calculated Mass (mg)"
    ax1.title = "Integrated (Normal-Density)"

    ##-- B --##
    ax2 = Axis(f[2, 1])

    df3 = df_a_normal
    sc1 = scatter!(ax2, df3[!, :ground_truth_mass_large], df3[!, :calculated_mass_large])
    errorbars!(
        ax2,
        df3[!, :ground_truth_mass_large],
        df3[!, :calculated_mass_large],
        rms(df3[!, :ground_truth_mass_large], df3[!, :calculated_mass_large]),
    )
    sc2 = scatter!(ax2, df3[!, :ground_truth_mass_medium], df3[!, :calculated_mass_medium])
    errorbars!(
        ax2,
        df3[!, :ground_truth_mass_medium],
        df3[!, :calculated_mass_medium],
        rms(df3[!, :ground_truth_mass_medium], df3[!, :calculated_mass_medium]),
    )
    sc3 = scatter!(
        ax2, df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small]; color=:red
    )
    errorbars!(
        ax2,
        df3[!, :ground_truth_mass_small],
        df3[!, :calculated_mass_small],
        rms(df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small]),
    )
    ln1 = lines!(ax2, [-1000, 1000], [-1000, 1000])
    ln2 = lines!(ax2, collect(1:1000), pred_a_norm; linestyle=:dashdot)
    Textbox(
        f[2, 1];
        placeholder="y = $(trunc(co3[2]; digits=3))x + $(trunc(co3[1]; digits=3)) \nr = $(trunc(r2_3; digits=3)) \nRMSE: $(trunc(rms_values3[1]; digits=3)) \nRMSD: $(trunc(rms_values3[2]; digits=3))",
        tellheight=false,
        tellwidth=false,
        boxcolor=:white,
        halign=:left,
        valign=:top,
        textsize=12,
    )

    xlims!(ax2; low=0, high=200)
    ylims!(ax2; low=0, high=200)
    ax2.xticks = [0, 50, 100, 150, 200]
    ax2.yticks = [0, 50, 100, 150, 200]
    ax2.xlabel = "Known Mass (mg)"
    ax2.ylabel = "Calculated Mass (mg)"
    ax2.title = "Agatston (Normal-Density)"

    ##-- LABELS --##
    f[1:2, 2] = Legend(
        f,
        [sc1, sc2, sc3, ln1, ln2],
        ["Large Inserts", "Medium Inserts", "Small Inserts", "Unity", "Fitted Line"];
        framevisible=false,
    )

    for (label, layout) in zip(["A", "B"], [f[1, 1], f[2, 1]])
        Label(
            layout[1, 1, TopLeft()],
            label;
            textsize=25,
            padding=(0, 60, 25, 0),
            halign=:right,
        )
    end

    save(
        "/Users/daleblack/Google Drive/Research/Papers/My Papers/cac-simulation/figures/linear_reg_norm_motion.png",
        f,
    )
    return f
end

# ╔═╡ ad497883-8eb1-463b-8aa5-096b3a7b0a56
with_theme(medphys_theme) do
    lin_reg_norm()
end

# ╔═╡ b25416e5-9312-4ed9-a814-d643c35bcf63
md"""
#### Low Density
"""

# ╔═╡ 76845429-d904-4fae-aa79-bb6341eeac61
let
    df = df_i_low
    gt_array = vec(
        hcat(
            df[!, :ground_truth_mass_large],
            df[!, :ground_truth_mass_medium],
            df[!, :ground_truth_mass_small],
        ),
    )
    calc_array = vec(
        hcat(
            df[!, :calculated_mass_large],
            df[!, :calculated_mass_medium],
            df[!, :calculated_mass_small],
        ),
    )
    data = DataFrame(; X=gt_array, Y=calc_array)
    global model_i_low
    model_i_low = lm(@formula(Y ~ X), data)
    global r2_2
    r2_2 = GLM.r2(model_i_low)
    global rms_values2
    rms_values2 = [
        rms(data[!, :X], data[!, :Y]), rmsd(data[!, :Y], GLM.predict(model_i_low))
    ]
end

# ╔═╡ 22d70538-4e53-4564-aa5c-66ece0d4cd43
begin
    newX2 = DataFrame(; X=collect(1:1000))
    pred_i_low = GLM.predict(model_i_low, newX2)
end

# ╔═╡ db8b8c00-3d75-4a7e-8903-00041b4112fd
co2 = coef(model_i_low)

# ╔═╡ 05febcf3-2ef8-461c-a6da-26705000e3cb
let
    df = df_a_low
    gt_array = vec(
        hcat(
            df[!, :ground_truth_mass_large],
            df[!, :ground_truth_mass_medium],
            df[!, :ground_truth_mass_small],
        ),
    )
    calc_array = vec(
        hcat(
            df[!, :calculated_mass_large],
            df[!, :calculated_mass_medium],
            df[!, :calculated_mass_small],
        ),
    )
    data = DataFrame(; X=gt_array, Y=calc_array)
    global model_a_low
    model_a_low = lm(@formula(Y ~ X), data)
    global r2_4
    r2_4 = GLM.r2(model_a_low)
    global rms_values4
    rms_values4 = [
        rms(data[!, :X], data[!, :Y]), rmsd(data[!, :Y], GLM.predict(model_a_low))
    ]
end

# ╔═╡ 4ccd5566-6555-479c-a36f-ad2b8400fa6d
begin
    newX4 = DataFrame(; X=collect(1:1000))
    pred_a_low = GLM.predict(model_a_low, newX4)
end

# ╔═╡ cd966539-e89a-4f9a-b748-4e906c81829e
co4 = coef(model_a_low)

# ╔═╡ 808845d8-8830-49b5-b511-2bc1b45593e4
function lin_reg_low()
    f = Figure()
    ##-- A --##
    ax1 = Axis(f[1, 1])

    df2 = df_i_low
    sc1 = scatter!(ax1, df2[!, :ground_truth_mass_large], df2[!, :calculated_mass_large])
    errorbars!(
        ax1,
        df2[!, :ground_truth_mass_large],
        df2[!, :calculated_mass_large],
        rms(df2[!, :ground_truth_mass_large], df2[!, :calculated_mass_large]),
    )
    sc2 = scatter!(ax1, df2[!, :ground_truth_mass_medium], df2[!, :calculated_mass_medium])
    errorbars!(
        ax1,
        df2[!, :ground_truth_mass_medium],
        df2[!, :calculated_mass_medium],
        rms(df2[!, :ground_truth_mass_medium], df2[!, :calculated_mass_medium]),
    )
    sc3 = scatter!(
        ax1, df2[!, :ground_truth_mass_small], df2[!, :calculated_mass_small]; color=:red
    )
    errorbars!(
        ax1,
        df2[!, :ground_truth_mass_small],
        df2[!, :calculated_mass_small],
        rms(df2[!, :ground_truth_mass_small], df2[!, :calculated_mass_small]),
    )
    ln1 = lines!(ax1, [-1000, 1000], [-1000, 1000])
    ln2 = lines!(ax1, collect(1:1000), pred_i_low; linestyle=:dashdot)
    Textbox(
        f[1, 1];
        placeholder="y = $(trunc(co2[2]; digits=3))x + $(trunc(co2[1]; digits=3)) \nr = $(trunc(r2_2; digits=3)) \nRMSE: $(trunc(rms_values2[1]; digits=3)) \nRMSD: $(trunc(rms_values2[2]; digits=3))",
        tellheight=false,
        tellwidth=false,
        boxcolor=:white,
        halign=:left,
        valign=:top,
        textsize=12,
    )

    xlims!(ax1; low=0, high=25)
    ylims!(ax1; low=0, high=25)
    ax1.xticks = [0, 5, 10, 15, 20, 25]
    ax1.yticks = [0, 5, 10, 15, 20, 25]
    ax1.xlabel = "Known Mass (mg)"
    ax1.ylabel = "Calculated Mass (mg)"
    ax1.title = "Integrated (Low-Density)"
    # hidedecorations!(ax1, ticklabels=false, ticks=false, label=false)

    ##-- B --##
    ax2 = Axis(f[2, 1])

    df4 = df_a_low
    sc1 = scatter!(ax2, df4[!, :ground_truth_mass_large], df4[!, :calculated_mass_large])
    errorbars!(
        ax2,
        df4[!, :ground_truth_mass_large],
        df4[!, :calculated_mass_large],
        rms(df4[!, :ground_truth_mass_large], df4[!, :calculated_mass_large]),
    )
    sc2 = scatter!(ax2, df4[!, :ground_truth_mass_medium], df4[!, :calculated_mass_medium])
    errorbars!(
        ax2,
        df4[!, :ground_truth_mass_medium],
        df4[!, :calculated_mass_medium],
        rms(df4[!, :ground_truth_mass_medium], df4[!, :calculated_mass_medium]),
    )
    sc3 = scatter!(
        ax2, df4[!, :ground_truth_mass_small], df4[!, :calculated_mass_small]; color=:red
    )
    errorbars!(
        ax2,
        df4[!, :ground_truth_mass_small],
        df4[!, :calculated_mass_small],
        rms(df4[!, :ground_truth_mass_small], df4[!, :calculated_mass_small]),
    )
    ln1 = lines!(ax2, [-1000, 1000], [-1000, 1000])
    ln2 = lines!(ax2, collect(1:1000), pred_a_low; linestyle=:dashdot)
    Textbox(
        f[2, 1];
        placeholder="y = $(trunc(co4[2]; digits=3))x + $(trunc(co4[1]; digits=3)) \nr = $(trunc(r2_4; digits=3)) \nRMSE: $(trunc(rms_values4[1]; digits=3)) \nRMSD: $(trunc(rms_values4[2]; digits=3))",
        tellheight=false,
        tellwidth=false,
        boxcolor=:white,
        halign=:left,
        valign=:top,
        textsize=12,
    )

    xlims!(ax2; low=0, high=25)
    ylims!(ax2; low=0, high=25)
    ax2.xticks = [0, 5, 10, 15, 20, 25]
    ax2.yticks = [0, 5, 10, 15, 20, 25]
    ax2.xlabel = "Known Mass (mg)"
    ax2.ylabel = "Calculated Mass (mg)"
    ax2.title = "Agatston (Low-Density)"
    # hidedecorations!(ax2, ticklabels=false, ticks=false, label=false)

    ##-- LABELS --##

    f[1:2, 2] = Legend(
        f,
        [sc1, sc2, sc3, ln1, ln2],
        ["Large Inserts", "Medium Inserts", "Small Inserts", "Unity", "Fitted Line"];
        framevisible=false,
    )

    for (label, layout) in zip(["A", "B"], [f[1, 1], f[2, 1]])
        Label(
            layout[1, 1, TopLeft()],
            label;
            textsize=25,
            padding=(0, 60, 25, 0),
            halign=:right,
        )
    end

    save(
        "/Users/daleblack/Google Drive/Research/Papers/My Papers/cac-simulation/figures/linear_reg_low_motion.png",
        f,
    )
    return f
end

# ╔═╡ c4d460d1-d2a0-436e-afec-82a7a7c34249
with_theme(medphys_theme) do
    lin_reg_low()
end

# ╔═╡ eb7b6dde-cfda-429d-951a-61eb26703d0f
md"""
## Reproducibility (FIG)
"""

# ╔═╡ 632786a6-a7f5-4a2e-9919-c37d9463d74c
md"""
#### SWCS
"""

# ╔═╡ 1d1c1b47-8806-4d1e-acd9-7707827b512e
array_s = hcat(
    df_s[!, :calculated_swcs_large],
    df_s[!, :calculated_swcs_medium],
    df_s[!, :calculated_swcs_small],
)

# ╔═╡ 7757de4b-e93c-4576-8d66-1e49566abc22
vec_s = vec(array_s);

# ╔═╡ a413ea64-33f1-4de2-84ec-64b4471274c5
md"""
#### Agatston
"""

# ╔═╡ e8bbc77f-d269-4ef3-889d-a0f0ed8d230c
vec_a = vec(array_a);

# ╔═╡ 64c83e13-68fd-40b8-95bc-c3a2b053d710
md"""
#### Integrated
"""

# ╔═╡ cbafab0a-31d4-44df-83ed-8d25113c91d1
array_i = hcat(
    df_i[!, :calculated_mass_large],
    df_i[!, :calculated_mass_medium],
    df_i[!, :calculated_mass_small],
)

# ╔═╡ 6c8bcd01-7e14-4e5b-b9a4-b0317f5ffe0f
vec_i = vec(array_i);

# ╔═╡ e17998a3-b47d-4f96-be7f-1c4f741bdb42
md"""
#### Known mass
"""

# ╔═╡ b75523b9-cdb7-4e89-bfd1-4231bd9764a5
array_k = hcat(
    df_a[!, :ground_truth_mass_large],
    df_a[!, :ground_truth_mass_medium],
    df_a[!, :ground_truth_mass_small],
);

# ╔═╡ 5e51727d-1cff-43f4-bb3f-b2d90f6630ce
vec_k = vec(array_k);

# ╔═╡ 847841f5-e161-43ea-a02d-dab852f4ca97
md"""
#### Integrated
"""

# ╔═╡ 1679bc7b-5d86-4eb0-b3dc-004165d3638c
path_integrated_r = "/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/output_repeated/integrated_scoring";

# ╔═╡ 06325405-93c1-4aeb-b8b5-c3a667e595ba
df_i_r = CSV.read(string(path_integrated_r, "/full2.csv"), DataFrame);

# ╔═╡ 70d951c6-4ad8-4645-b187-ca7fe4f9ac0d
df_i_low_r, df_i_normal_r = groupby(df_i_r, :DENSITY);

# ╔═╡ 3be56210-24cf-4153-bd3e-05ecb50ab783
df_i_low_small_r, df_i_low_medium_r, df_i_low_large_r = groupby(df_i_low_r, :SIZE);

# ╔═╡ 479f5aea-eb39-4eba-994c-5949c4cca64b
df_i_normal_small_r, df_i_normal_medium_r, df_i_normal_large_r = groupby(
    df_i_normal_r, :SIZE
);

# ╔═╡ 04dd49f0-e64c-431d-915f-b2e59f8ba3a0
array_i_r = hcat(
    df_i_r[!, :calculated_mass_large],
    df_i_r[!, :calculated_mass_medium],
    df_i_r[!, :calculated_mass_small],
)

# ╔═╡ 917ecb3c-704c-42f7-baca-d90b50d3b5b6
begin
    mean_bkg_i = mean(df_i[!, :mass_bkg])
    mean_bkg_i_r = mean(df_i_r[!, :mass_bkg])
end

# ╔═╡ 86e05664-74d6-4181-bca5-d82d4719ebac
begin
    idxs_i_r_large = Tuple(findall(x -> x < mean_bkg_i_r, array_i_r[:, 1]))
    idxs_i_large = Tuple(findall(x -> x < mean_bkg_i, array_i[:, 1]))
    indxs_i_large_tot = (idxs_i_r_large..., idxs_i_large...)
    indx_i_large_tot = Tuple(unique(indxs_i_large_tot))
end

# ╔═╡ 2a659abb-ea56-437d-a66e-dd509f107946
begin
    idxs_i_r_med = Tuple(findall(x -> x < mean_bkg_i_r, array_i_r[:, 2]))
    idxs_i_med = Tuple(findall(x -> x < mean_bkg_i, array_i[:, 2]))
    indxs_i_med_tot = (idxs_i_r_med..., idxs_i_med...)
    indx_i_med_tot = Tuple(unique(indxs_i_med_tot))
end

# ╔═╡ e1ddab59-ce43-4d26-afdf-ffbd764f34c0
begin
    idxs_i_r_small = Tuple(findall(x -> x < mean_bkg_i_r, array_i_r[:, 3]))
    idxs_i_small = Tuple(findall(x -> x < mean_bkg_i, array_i[:, 3]))
    indxs_i_small_tot = (idxs_i_r_small..., idxs_i_small...)
    indx_i_small_tot = Tuple(unique(indxs_i_small_tot))
end

# ╔═╡ dab35b27-5585-4a67-a5b5-536b4699fb8c
begin
    if length(indx_i_large_tot) > 0
        df_i_large = df_i[Not(indx_i_large_tot...), :]
        df_i_r_large = df_i_r[Not(indx_i_large_tot...), :]
    end
    if length(indx_i_med_tot) > 0
        df_i_med = df_i[Not(indx_i_med_tot...), :]
        df_i_r_med = df_i_r[Not(indx_i_med_tot...), :]
    end
    if length(indx_i_small_tot) > 0
        df_i_small = df_i[Not(indx_i_small_tot...), :]
        df_i_r_small = df_i_r[Not(indx_i_small_tot...), :]
    end
end;

# ╔═╡ 11d27f35-5250-4a5a-9486-fc490a3f673e
let
    r_array = [
        Tuple(df_i_r[!, :calculated_mass_large])...,
        Tuple(df_i_r[!, :calculated_mass_medium])...,
        Tuple(df_i_r_small[!, :calculated_mass_small])...,
    ]
    calc_array = [
        Tuple(df_i[!, :calculated_mass_large])...,
        Tuple(df_i[!, :calculated_mass_medium])...,
        Tuple(df_i_small[!, :calculated_mass_small])...,
    ]
    data = DataFrame(; X=r_array, Y=calc_array)
    global model_ir
    model_ir = lm(@formula(Y ~ X), data)
    global r2_1r
    r2_1r = GLM.r2(model_ir)
    global rms_values1r
    rms_values1r = [rms(data[!, :X], data[!, :Y]), rmsd(data[!, :Y], GLM.predict(model_ir))]
end

# ╔═╡ dd4bdc8e-8a19-4d40-9e83-2de28e3e527d
begin
    newX1r = DataFrame(; X=collect(1:1000))
    pred_ir = GLM.predict(model_ir, newX1r)
end

# ╔═╡ 5798dc6c-c43c-492c-b8ae-2ab9a7bb5d02
co1r = coef(model_ir)

# ╔═╡ 0d5255e0-bb2a-44c5-a68e-7da1d1620126
md"""
#### Agatston
"""

# ╔═╡ 2f6e13cb-d23f-4e42-9dcc-952ecabbccc9
path_agat_r = "/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/output_repeated/agatston";

# ╔═╡ a58871c6-7633-4e61-beb3-840c40ddfbd2
df_a_r = CSV.read(string(path_agat_r, "/full2.csv"), DataFrame);

# ╔═╡ 6d8e3c08-6c20-4d28-8294-01fd62d592ce
df_a_low_r, df_a_normal_r = groupby(df_a_r, :DENSITY);

# ╔═╡ 3fa6abba-b1d8-4e64-a4cc-d01dba6ec581
df_a_low_small_r, df_a_low_medium_r, df_a_low_large_r = groupby(df_a_low_r, :SIZE);

# ╔═╡ 992d5958-28fb-4b8b-b333-979544048d0e
df_a_normal_small_r, df_a_normal_medium_r, df_a_normal_large_r = groupby(
    df_a_normal_r, :SIZE
);

# ╔═╡ dfec184d-d5c4-481a-88a6-99327b257a3a
array_a_r = hcat(
    df_a_r[!, :calculated_mass_large],
    df_a_r[!, :calculated_mass_medium],
    df_a_r[!, :calculated_mass_large],
);

# ╔═╡ 2ee363a5-9847-4af0-a793-68e4f2d5c875
begin
    idxs_a_r_large = Tuple(findall(x -> x == 0, array_a_r[:, 1]))
    idxs_a_large = Tuple(findall(x -> x == 0, array_a[:, 1]))
    indxs_a_large_tot = (idxs_a_r_large..., idxs_a_large...)
    indx_a_large_tot = Tuple(unique(indxs_a_large_tot))
end

# ╔═╡ d1f6aa00-cb54-435f-a793-22312a5d6c79
begin
    idxs_a_r_med = Tuple(findall(x -> x == 0, array_a_r[:, 2]))
    idxs_a_med = Tuple(findall(x -> x == 0, array_a[:, 2]))
    indxs_a_med_tot = (idxs_a_r_med..., idxs_a_med...)
    indx_a_med_tot = Tuple(unique(indxs_a_med_tot))
end

# ╔═╡ f07cf34d-225e-4416-b19c-db8d1a39c6c9
begin
    idxs_a_r_small = Tuple(findall(x -> x == 0, array_a_r[:, 3]))
    idxs_a_small = Tuple(findall(x -> x == 0, array_a[:, 3]))
    indxs_a_small_tot = (idxs_a_r_small..., idxs_a_small...)
    indx_a_small_tot = Tuple(unique(indxs_a_small_tot))
end

# ╔═╡ 24b81917-f6a3-4b64-b5e0-274193d05201
begin
    if length(indx_a_large_tot) > 0
        df_a_large = df_a[Not(indx_a_large_tot...), :]
        df_a_r_large = df_a_r[Not(indx_a_large_tot...), :]
    end
    if length(indx_a_med_tot) > 0
        df_a_med = df_a[Not(indx_a_med_tot...), :]
        df_a_r_med = df_a_r[Not(indx_a_med_tot...), :]
    end
    if length(indx_a_small_tot) > 0
        df_a_small = df_a[Not(indx_a_small_tot...), :]
        df_a_r_small = df_a_r[Not(indx_a_small_tot...), :]
    end
end;

# ╔═╡ fb448e70-45ac-470b-90a5-9f45ca2a2d28
let
    r_array = [
        Tuple(df_a_r_large[!, :calculated_mass_large])...,
        Tuple(df_a_r_med[!, :calculated_mass_medium])...,
        Tuple(df_a_r_small[!, :calculated_mass_small])...,
    ]
    calc_array = [
        Tuple(df_a_large[!, :calculated_mass_large])...,
        Tuple(df_a_med[!, :calculated_mass_medium])...,
        Tuple(df_a_small[!, :calculated_mass_small])...,
    ]
    data = DataFrame(; X=r_array, Y=calc_array)
    global model_ar
    model_ar = lm(@formula(Y ~ X), data)
    global r2_2r
    r2_2r = GLM.r2(model_ar)
    global rms_values2r
    rms_values2r = [rms(data[!, :X], data[!, :Y]), rmsd(data[!, :Y], GLM.predict(model_ar))]
end

# ╔═╡ 1e0a0439-f21b-464e-a137-b147c448792a
begin
    newX2r = DataFrame(; X=collect(1:1000))
    pred_ar = GLM.predict(model_ar, newX2r)
end

# ╔═╡ ce90f9ee-97d9-4bb8-8a07-dccbdf98675e
co2r = coef(model_ar)

# ╔═╡ 2fd240f7-8f4e-4eca-b83c-89611b0b5dac
md"""
#### SWCS
"""

# ╔═╡ ed434ee4-6288-4d62-b3a0-bec7fd3ffad8
path_swcs_r = "/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/output_repeated/swcs";

# ╔═╡ 8107c0a6-24c8-4e8d-9092-7b86af925cc2
df_s_r = CSV.read(string(path_swcs_r, "/full2.csv"), DataFrame);

# ╔═╡ 5230514e-4ee0-4bee-bab5-5846b320c466
df_s_low_r, df_s_normal_r = groupby(df_s_r, :DENSITY);

# ╔═╡ cea55edf-d2ea-4107-89d0-7a2aa760f6df
df_s_low_small_r, df_s_low_medium_r, df_s_low_large_r = groupby(df_s_low_r, :SIZE);

# ╔═╡ 94d401c5-f160-4113-bc2d-a8617d0cc0da
df_s_normal_small_r, df_s_normal_medium_r, df_s_normal_large_r = groupby(
    df_s_normal_r, :SIZE
);

# ╔═╡ ba2f0bdb-97c4-4db7-b614-e7f1a392cdb7
array_s_r = hcat(
    df_s_r[!, :calculated_swcs_large],
    df_s_r[!, :calculated_swcs_medium],
    df_s_r[!, :calculated_swcs_small],
);

# ╔═╡ d5220791-a3e8-44b6-aae3-cfe2a52838df
begin
    mean_bkg_s = mean(df_s[!, :swcs_bkg])
    mean_bkg_s_r = mean(df_s_r[!, :swcs_bkg])
end

# ╔═╡ 9996d074-7e32-467e-9b15-e806f9237606
begin
    idxs_s_r_large = Tuple(findall(x -> x < mean_bkg_s_r, array_s_r[:, 1]))
    idxs_s_large = Tuple(findall(x -> x < mean_bkg_s, array_s[:, 1]))
    indxs_s_large_tot = (idxs_s_r_large..., idxs_s_large...)
    indx_s_large_tot = Tuple(unique(indxs_s_large_tot))
end

# ╔═╡ a1aa6205-b436-459a-a60e-a566e784457b
begin
    idxs_s_r_med = Tuple(findall(x -> x < mean_bkg_s_r, array_s_r[:, 2]))
    idxs_s_med = Tuple(findall(x -> x < mean_bkg_s, array_s[:, 2]))
    indxs_s_med_tot = (idxs_s_r_med..., idxs_s_med...)
    indx_s_med_tot = Tuple(unique(indxs_s_med_tot))
end

# ╔═╡ 4ed80ff1-1e52-4f2d-a758-73cd40c75bf0
begin
    idxs_s_r_small = Tuple(findall(x -> x < mean_bkg_s_r, array_s_r[:, 3]))
    idxs_s_small = Tuple(findall(x -> x < mean_bkg_s, array_s[:, 3]))
    indxs_s_small_tot = (idxs_s_r_small..., idxs_s_small...)
    indx_s_small_tot = Tuple(unique(indxs_s_small_tot))
end

# ╔═╡ 94ad5025-2c74-4970-9ce7-28ec698e254a
begin
    if length(indx_s_large_tot) > 0
        df_s_large = df_s[Not(indx_s_large_tot...), :]
        df_s_r_large = df_s_r[Not(indx_s_large_tot...), :]
    end
    if length(indx_s_med_tot) > 0
        df_s_med = df_s[Not(indx_s_med_tot...), :]
        df_s_r_med = df_s_r[Not(indx_s_med_tot...), :]
    end
    if length(indx_s_small_tot) > 0
        df_s_small = df_s[Not(indx_s_small_tot...), :]
        df_s_r_small = df_s_r[Not(indx_s_small_tot...), :]
    end
end;

# ╔═╡ fa33e295-cb79-4107-86ae-6524938bcf96
let
    r_array = [
        Tuple(df_s_r_large[!, :calculated_mass_large])...,
        Tuple(df_s_r_med[!, :calculated_mass_medium])...,
        Tuple(df_s_r_small[!, :calculated_mass_small])...,
    ]
    calc_array = [
        Tuple(df_s_large[!, :calculated_mass_large])...,
        Tuple(df_s_med[!, :calculated_mass_medium])...,
        Tuple(df_s_small[!, :calculated_mass_small])...,
    ]
    data = DataFrame(; X=r_array, Y=calc_array)
    global model_sr
    model_sr = lm(@formula(Y ~ X), data)
    global r2_3r
    r2_3r = GLM.r2(model_sr)
    global rms_values3r
    rms_values3r = [rms(data[!, :X], data[!, :Y]), rmsd(data[!, :Y], GLM.predict(model_sr))]
end

# ╔═╡ 9d4c20c5-61ad-49ba-ac18-64db25038d1e
begin
    newX3r = DataFrame(; X=collect(1:1000))
    pred_sr = GLM.predict(model_sr, newX3r)
end

# ╔═╡ 813b929c-ccd8-4ee7-a496-87647444e9ac
co3r = coef(model_sr)

# ╔═╡ dbbfd45b-03a5-42eb-b155-545f5de072ba
function reprod()
    f = Figure()

    ##-- A --##
    ax1 = Axis(f[1, 1])
    scatter!(ax1, df_i_r[!, :calculated_mass_large], df_i[!, :calculated_mass_large])
    scatter!(ax1, df_i_r[!, :calculated_mass_medium], df_i[!, :calculated_mass_medium])
    scatter!(
        ax1,
        df_i_r_small[!, :calculated_mass_small],
        df_i_small[!, :calculated_mass_small];
        color=:red,
    )
    lines!(ax1, [-1000, 1000], [-1000, 1000]; label="Unity")
    lines!(ax1, collect(1:1000), pred_ir; linestyle=:dashdot)
    if co1r[1] > 0
        Textbox(
            f[1, 1];
            placeholder="y = $(trunc(co1r[2]; digits=3))x+$(trunc(co1r[1]; digits=3)) \nr = $(trunc(r2_1r; digits=3)) \nRMSE: $(trunc(rms_values1r[1]; digits=3)) \nRMSD: $(trunc(rms_values1r[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            textsize=12,
        )
    else
        Textbox(
            f[1, 1];
            placeholder="y = $(trunc(co1r[2]; digits=3))x$(trunc(co1r[1]; digits=3)) \nr = $(trunc(r2_1r; digits=3)) \nRMSE: $(trunc(rms_values1r[1]; digits=3)) \nRMSD: $(trunc(rms_values1r[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            textsize=12,
        )
    end

    xlims!(ax1; low=0, high=200)
    ylims!(ax1; low=0, high=200)
    ax1.xticks = [0, 50, 100, 150, 200]
    ax1.yticks = [0, 50, 100, 150, 200]
    ax1.xlabel = "Mass 1 (mg)"
    ax1.ylabel = "Mass 2 (mg)"
    ax1.title = "Integrated"

    ##-- B --##
    ax2 = Axis(f[2, 1])
    scatter!(
        ax2,
        df_a_r_large[!, :calculated_mass_large],
        df_a_large[!, :calculated_mass_large];
        label="Large Inserts",
    )
    scatter!(
        ax2,
        df_a_r_med[!, :calculated_mass_medium],
        df_a_med[!, :calculated_mass_medium];
        label="Medium Inserts",
    )
    scatter!(
        ax2,
        df_a_r_small[!, :calculated_mass_small],
        df_a_small[!, :calculated_mass_small];
        label="Small Inserts",
        color=:red,
    )
    lines!(ax2, [-1000, 1000], [-1000, 1000]; label="Unity")
    lines!(ax2, collect(1:1000), pred_ar; linestyle=:dashdot)
    if co2r[1] > 0
        Textbox(
            f[2, 1];
            placeholder="y = $(trunc(co2r[2]; digits=3))x+$(trunc(co2r[1]; digits=3)) \nr = $(trunc(r2_2r; digits=3)) \nRMSE: $(trunc(rms_values2r[1]; digits=3)) \nRMSD: $(trunc(rms_values2r[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            textsize=12,
        )
    else
        Textbox(
            f[2, 1];
            placeholder="y = $(trunc(co2r[2]; digits=3))x$(trunc(co2r[1]; digits=3)) \nr = $(trunc(r2_2r; digits=3)) \nRMSE: $(trunc(rms_values2r[1]; digits=3)) \nRMSD: $(trunc(rms_values2r[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            textsize=12,
        )
    end

    xlims!(ax2; low=0, high=200)
    ylims!(ax2; low=0, high=200)
    ax2.xticks = [0, 50, 100, 150, 200]
    ax2.yticks = [0, 50, 100, 150, 200]
    ax2.xlabel = "Mass 1 (mg)"
    ax2.ylabel = "Mass 2 (mg)"
    ax2.title = "Agatston"

    # ##-- C --##
    ax3 = Axis(f[3, 1])
    scatter!(
        ax3,
        df_s_r_large[!, :calculated_swcs_large],
        df_s_large[!, :calculated_swcs_large];
        label="Large Inserts",
    )
    scatter!(
        ax3,
        df_s_r_med[!, :calculated_swcs_medium],
        df_s_med[!, :calculated_swcs_medium];
        label="Medium Inserts",
    )
    scatter!(
        ax3,
        df_s_r_small[!, :calculated_swcs_small],
        df_s_small[!, :calculated_swcs_small];
        label="Small Inserts",
        color=:red,
    )
    lines!(ax3, [-1000, 1000], [-1000, 1000]; label="Unity")
    lines!(ax3, collect(1:1000), pred_sr; linestyle=:dashdot, label="Fitted Line")
    if co3r[1] > 0
        Textbox(
            f[3, 1];
            placeholder="y = $(trunc(co3r[2]; digits=3))x+$(trunc(co3r[1]; digits=3)) \nr = $(trunc(r2_3r; digits=3)) \nRMSE: $(trunc(rms_values3r[1]; digits=3)) \nRMSD: $(trunc(rms_values3r[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            textsize=12,
        )
    else
        Textbox(
            f[3, 1];
            placeholder="y = $(trunc(co3r[2]; digits=3))x$(trunc(co3r[1]; digits=3)) \nr = $(trunc(r2_3r; digits=3)) \nRMSE: $(trunc(rms_values3r[1]; digits=3)) \nRMSD: $(trunc(rms_values3r[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            textsize=12,
        )
    end

    xlims!(ax3; low=0, high=500)
    ylims!(ax3; low=0, high=500)
    ax3.xticks = [0, 125, 250, 375, 500]
    ax3.yticks = [0, 125, 250, 375, 500]
    ax3.xlabel = "SWCS 1"
    ax3.ylabel = "SWCS 2"
    ax3.title = "Spatially Weighted"

    ##-- LABELS --##
    f[2, 2] = Legend(f, ax3; framevisible=false)
    for (label, layout) in zip(["A", "B", "C"], [f[1, 1], f[2, 1], f[3, 1]])
        Label(
            layout[1, 1, TopLeft()],
            label;
            textsize=25,
            padding=(0, -10, 5, 0),
            halign=:right,
        )
    end

    save(
        "/Users/daleblack/Google Drive/Research/Papers/My Papers/cac-simulation/figures/reprod_motion.png",
        f,
    )
    return f
end

# ╔═╡ 85f53aca-3136-4144-8992-966b435d6c8d
with_theme(medphys_theme) do
    reprod()
end

# ╔═╡ fd50c12c-3b8b-4168-853d-e7236478cba5
# md"""
# ## Motion Blurring (FIG)
# """

# ╔═╡ b4b95e73-a1a5-4491-bbe0-18c0cf70ccbf
# function motion_blur()
# 	f = Figure()

# 	##-- A --##
# 	ax1 = Axis(f[1, 1])

# 	df = df_i
# 	scatter!(ax1, df[!, :ground_truth_mass_large], df[!, :calculated_mass_large])
# 	errorbars!(ax1, df[!, :ground_truth_mass_large], df[!, :calculated_mass_large], rms(df[!, :ground_truth_mass_large], df[!, :calculated_mass_large]))
# 	scatter!(ax1, df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium])
# 	errorbars!(ax1, df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium], rms(df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium]))
# 	scatter!(ax1, df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], color=:red)
# 	errorbars!(ax1, df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], rms(df[!, :ground_truth_mass_small], df[!, :calculated_mass_small]))
# 	lines!(ax1, [-1000, 1000], [-1000, 1000],)
# 	lines!(ax1, collect(1:1000), pred_ii, linestyle=:dashdot)
# 	Textbox(
# 		f[1, 1], 
# 		placeholder = "y = $(trunc(coii[2]; digits=3))x + $(trunc(coii[1]; digits=3)) \nr = $(trunc(r2_ii; digits=3)) \nRMSE: $(trunc(rms_valuesii[1]; digits=3)) \nRMSD: $(trunc(rms_valuesii[2]; digits=3))", 
# 		tellheight = false,
#         tellwidth = false,
# 		boxcolor=:white,
# 		halign=:left,
# 		valign=:top,
# 		textsize=12
# 	)

# 	xlims!(ax1, low=0, high=200)
# 	ylims!(ax1, low=0, high=200)
# 	ax1.xticks = [0, 50, 100, 150, 200]
# 	ax1.yticks = [0, 50, 100, 150, 200]
# 	ax1.xlabel = "Known Mass (mg)"
# 	ax1.ylabel = "Calculated Mass (mg)"
# 	ax1.title = "Integrated"

# 	##-- B --##
# 	ax2 = Axis(f[2, 1])

# 	df3 = df_a
# 	sc1=scatter!(ax2, df3[!, :ground_truth_mass_large], df3[!, :calculated_mass_large])
# 	errorbars!(ax2, df3[!, :ground_truth_mass_large], df3[!, :calculated_mass_large], rms(df3[!, :ground_truth_mass_large], df3[!, :calculated_mass_large]))
# 	sc2=scatter!(ax2, df3[!, :ground_truth_mass_medium], df3[!, :calculated_mass_medium])
# 	errorbars!(ax2, df3[!, :ground_truth_mass_medium], df3[!, :calculated_mass_medium], rms(df3[!, :ground_truth_mass_medium], df3[!, :calculated_mass_medium]))
# 	sc3=scatter!(ax2, df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small], color=:red)
# 	errorbars!(ax2, df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small], rms(df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small]))
# 	ln1=lines!(ax2, [-1000, 1000], [-1000, 1000])
# 	ln2=lines!(ax2, collect(1:1000), pred_aa, linestyle=:dashdot)
# 	Textbox(
# 		f[2, 1], 
# 		placeholder = "y = $(trunc(coaa[2]; digits=3))x + $(trunc(coaa[1]; digits=3)) \nr = $(trunc(r2_aa; digits=3)) \nRMSE: $(trunc(rms_valuesaa[1]; digits=3)) \nRMSD: $(trunc(rms_valuesaa[2]; digits=3))", 
# 		tellheight = false,
#         tellwidth = false,
# 		boxcolor=:white,
# 		halign=:left,
# 		valign=:top,
# 		textsize=12
# 	)

# 	xlims!(ax2, low=0, high=200)
# 	ylims!(ax2, low=0, high=200)
# 	ax2.xticks = [0, 50, 100, 150, 200]
# 	ax2.yticks = [0, 50, 100, 150, 200]
# 	ax2.xlabel = "Known Mass (mg)"
# 	ax2.ylabel = "Calculated Mass (mg)"
# 	ax2.title = "Agatston"

# 	##-- LABELS --##
# 	f[1:2, 2] = Legend(f, [sc1, sc2, sc3, ln1, ln2], ["Large Inserts", "Medium Inserts", "Small Inserts", "Unity", "Fitted Line"], framevisible = false)
# 	for (label, layout) in zip(["A", "B"], [f[1,1], f[2,1]])
# 	    Label(layout[1, 1, TopLeft()], label,
# 	        textsize = 25,
# 	        padding = (0, 60, 25, 0),
# 	        halign = :right)
# 	end

# 	save("/Users/daleblack/Google Drive/Research/Papers/My Papers/cac-simulation/figures/motionblur.png", f)
# 	f
# end

# ╔═╡ 57b206a6-3e9a-44af-b063-de358b597cab
# with_theme(medphys_theme) do
# 	motion_blur()
# end

# ╔═╡ e092697c-cf33-4e7a-ba31-32d0f74576e8
# md"""
# #### Integrated
# """

# ╔═╡ 3e85daec-fbe8-4b4c-9b01-37d345e6b762
# df_ii0, df_ii05, df_ii1, df_ii15, df_ii2 = groupby(df_i, :blur);

# ╔═╡ c4590c18-7b89-47b4-9788-5959099d8096
# let
# 	df = df_i
# 	gt_array = vec(hcat(df[!, :ground_truth_mass_large], df[!, :ground_truth_mass_medium], df[!, :ground_truth_mass_small]))
# 	calc_array = vec(hcat(df[!, :calculated_mass_large], df[!, :calculated_mass_medium], df[!, :calculated_mass_small]))
# 	data = DataFrame(
# 		X = gt_array,
# 		Y= calc_array
# 	)
# 	global model_ii
# 	model_ii = lm(@formula(Y ~ X), data)
# 	global r2_ii
# 	r2_ii = GLM.r2(model_ii)
# 	global rms_valuesii
# 	rms_valuesii = [
# 		rms(data[!, :X], data[!, :Y]),
# 		rmsd(data[!, :Y], GLM.predict(model_ii))
# 	]
# end

# ╔═╡ 6e6c5d71-4a1a-4c31-9a23-96d99ac5f1f0
# begin
# 	newXii = DataFrame(X=collect(1:1000));
# 	pred_ii = GLM.predict(model_ii, newXii)
# end

# ╔═╡ b03ab3c2-598b-440e-a00d-a93bc483828f
# coii = coef(model_ii)

# ╔═╡ 9bf77bda-92dc-43c5-aaa0-711e5f127868
# md"""
# #### Agatston
# """

# ╔═╡ 86f392b6-9332-43b6-b09d-aecbdf72933a
# let
# 	df = df_a
# 	gt_array = vec(hcat(df[!, :ground_truth_mass_large], df[!, :ground_truth_mass_medium], df[!, :ground_truth_mass_small]))
# 	calc_array = vec(hcat(df[!, :calculated_mass_large], df[!, :calculated_mass_medium], df[!, :calculated_mass_small]))
# 	data = DataFrame(
# 		X = gt_array,
# 		Y= calc_array
# 	)
# 	global model_aa
# 	model_aa = lm(@formula(Y ~ X), data)
# 	global r2_aa
# 	r2_aa = GLM.r2(model_aa)
# 	global rms_valuesaa
# 	rms_valuesaa = [
# 		rms(data[!, :X], data[!, :Y]),
# 		rmsd(data[!, :Y], GLM.predict(model_aa))
# 	]
# end

# ╔═╡ 45fc61c9-1173-4f59-a883-974e1029edc2
# begin
# 	newXaa = DataFrame(X=collect(1:1000));
# 	pred_aa = GLM.predict(model_aa, newXaa)
# end

# ╔═╡ 7c31bf02-dcd3-4880-85eb-0247e646177d
# coaa = coef(model_aa)

# ╔═╡ 7ec6d916-685b-4044-bca6-cca08b51def8
# df_aa0, df_aa05, df_aa1, df_aa15, df_aa2 = groupby(df_a, :blur);

# ╔═╡ Cell order:
# ╠═e2b53960-0e99-11ed-27bd-5f46a57ba02f
# ╠═b39114c0-182e-44dc-b00b-906c81921ebd
# ╟─3389da26-b3b3-455c-a382-6bb44437cd05
# ╟─88d08c61-2f3b-475f-b392-031cffbea40b
# ╠═a8eec50e-2993-4bfe-94da-b943a919e6e4
# ╠═76df194b-31bc-4db6-ba99-f6a35b3c0a65
# ╠═450c6c6c-96ae-42cf-bc6d-5c80e99a2fa7
# ╠═dde42a35-ef12-45d7-a0ae-ede776fb555a
# ╠═97e09092-c961-47de-95c0-91e791a938a0
# ╠═b14db409-696e-44e2-a4bc-afe263c3647d
# ╟─8b77c804-5f9a-4919-a55f-5ab7762d291e
# ╠═69468127-7c21-4d34-9f97-f59cd301ea3d
# ╠═9b5ed20b-36b6-4bfb-9cac-97a5f1a0b947
# ╠═1ad42be8-6538-45c0-bed3-cf6f5da8b082
# ╠═d692b7b0-a12b-4b3a-a040-15634bcd252c
# ╠═c3583be4-9add-4cf5-a738-6b63bd3a3712
# ╟─1d65a547-c69e-48d5-8ce2-d0827ee87adc
# ╠═d934cfb8-78ae-4fcd-a568-e5b74409d1c8
# ╠═e7be547b-88a0-47b2-a84b-eaa48ec7e4eb
# ╠═0687e798-a30c-44cf-a03c-99308a561113
# ╠═97afb112-d5fe-4014-92de-16efb1d3cb1d
# ╠═535f8216-9aeb-43c9-8e51-a4a93eeda076
# ╟─ae26e32a-dff5-4e02-b03e-bd13cf79b14b
# ╟─20b127df-a98e-4657-b5e2-dbd297a10c81
# ╟─cb05cac7-ecd9-44b2-8fe7-f80d272b482a
# ╟─ba9b6361-115d-47df-808f-e67147041427
# ╟─10fe5219-e93b-446a-be5d-9da6201b7fbd
# ╟─8ae299c5-f2a2-4ed4-bd57-45cc065ade7a
# ╠═2865797a-c029-4c7a-80de-e717b33f07fe
# ╠═76ce94f1-3fbf-461a-8dbc-7fffd7e3fabb
# ╠═db624048-33f3-4f4c-aa67-cfacd72777d5
# ╠═82fa5c00-6801-400b-b2eb-b3b0f45672ea
# ╠═f0240c36-8b76-4116-ac55-43e099e6cd62
# ╟─71626005-b2d5-4afb-ab65-6606e6f88a59
# ╠═632d337c-5af7-45c2-a7ba-d3a757768361
# ╠═9bbfd561-a51d-44af-ab06-f9afe555fc98
# ╟─509ae5c2-4423-4a61-87b7-f119615cd19c
# ╠═4d00a89d-2467-4505-9d87-e4066e107a6e
# ╠═c69f965d-196e-481a-8622-b7e988d17a2b
# ╠═2d3d0da0-587d-4a85-9360-82c7749a1945
# ╠═50852830-cddf-4208-a823-2cfb86c48ec6
# ╠═33fc9890-72d0-4103-b393-89aae9e216ea
# ╠═2fe1b4a3-a2ac-40ee-b90e-1906a33a5438
# ╟─dcdc9952-0fff-420b-8f30-b4da6ab6a5f9
# ╟─aeb864db-fba4-4d9c-8120-672339aaa90a
# ╟─808845d8-8830-49b5-b511-2bc1b45593e4
# ╟─ad497883-8eb1-463b-8aa5-096b3a7b0a56
# ╟─c4d460d1-d2a0-436e-afec-82a7a7c34249
# ╟─14ddfbb8-1181-4843-b446-9352415d05e9
# ╠═c26af14f-0cf7-4a1e-8ed8-abd3eb247d73
# ╠═aadb907e-dead-482f-8768-f6851c08a18c
# ╠═a1d4e6f8-b43c-40bd-a3f4-fba304871cb9
# ╠═33fbcb79-10a4-454b-8114-1f0ffd0904a2
# ╠═4b1ebeb0-1797-457d-b52f-0ea7d8e57318
# ╠═1f196a66-9c0a-43d2-b41f-b302954753c2
# ╟─b25416e5-9312-4ed9-a814-d643c35bcf63
# ╠═76845429-d904-4fae-aa79-bb6341eeac61
# ╠═22d70538-4e53-4564-aa5c-66ece0d4cd43
# ╠═db8b8c00-3d75-4a7e-8903-00041b4112fd
# ╠═05febcf3-2ef8-461c-a6da-26705000e3cb
# ╠═4ccd5566-6555-479c-a36f-ad2b8400fa6d
# ╠═cd966539-e89a-4f9a-b748-4e906c81829e
# ╟─eb7b6dde-cfda-429d-951a-61eb26703d0f
# ╟─dbbfd45b-03a5-42eb-b155-545f5de072ba
# ╟─85f53aca-3136-4144-8992-966b435d6c8d
# ╟─632786a6-a7f5-4a2e-9919-c37d9463d74c
# ╠═1d1c1b47-8806-4d1e-acd9-7707827b512e
# ╠═7757de4b-e93c-4576-8d66-1e49566abc22
# ╟─a413ea64-33f1-4de2-84ec-64b4471274c5
# ╠═e8bbc77f-d269-4ef3-889d-a0f0ed8d230c
# ╟─64c83e13-68fd-40b8-95bc-c3a2b053d710
# ╠═cbafab0a-31d4-44df-83ed-8d25113c91d1
# ╠═6c8bcd01-7e14-4e5b-b9a4-b0317f5ffe0f
# ╟─e17998a3-b47d-4f96-be7f-1c4f741bdb42
# ╠═b75523b9-cdb7-4e89-bfd1-4231bd9764a5
# ╠═5e51727d-1cff-43f4-bb3f-b2d90f6630ce
# ╟─847841f5-e161-43ea-a02d-dab852f4ca97
# ╠═1679bc7b-5d86-4eb0-b3dc-004165d3638c
# ╠═06325405-93c1-4aeb-b8b5-c3a667e595ba
# ╠═70d951c6-4ad8-4645-b187-ca7fe4f9ac0d
# ╠═3be56210-24cf-4153-bd3e-05ecb50ab783
# ╠═479f5aea-eb39-4eba-994c-5949c4cca64b
# ╠═04dd49f0-e64c-431d-915f-b2e59f8ba3a0
# ╠═917ecb3c-704c-42f7-baca-d90b50d3b5b6
# ╠═86e05664-74d6-4181-bca5-d82d4719ebac
# ╠═2a659abb-ea56-437d-a66e-dd509f107946
# ╠═e1ddab59-ce43-4d26-afdf-ffbd764f34c0
# ╠═dab35b27-5585-4a67-a5b5-536b4699fb8c
# ╠═11d27f35-5250-4a5a-9486-fc490a3f673e
# ╠═dd4bdc8e-8a19-4d40-9e83-2de28e3e527d
# ╠═5798dc6c-c43c-492c-b8ae-2ab9a7bb5d02
# ╟─0d5255e0-bb2a-44c5-a68e-7da1d1620126
# ╠═2f6e13cb-d23f-4e42-9dcc-952ecabbccc9
# ╠═a58871c6-7633-4e61-beb3-840c40ddfbd2
# ╠═6d8e3c08-6c20-4d28-8294-01fd62d592ce
# ╠═3fa6abba-b1d8-4e64-a4cc-d01dba6ec581
# ╠═992d5958-28fb-4b8b-b333-979544048d0e
# ╠═dfec184d-d5c4-481a-88a6-99327b257a3a
# ╠═2ee363a5-9847-4af0-a793-68e4f2d5c875
# ╠═d1f6aa00-cb54-435f-a793-22312a5d6c79
# ╠═f07cf34d-225e-4416-b19c-db8d1a39c6c9
# ╠═24b81917-f6a3-4b64-b5e0-274193d05201
# ╠═fb448e70-45ac-470b-90a5-9f45ca2a2d28
# ╠═1e0a0439-f21b-464e-a137-b147c448792a
# ╠═ce90f9ee-97d9-4bb8-8a07-dccbdf98675e
# ╟─2fd240f7-8f4e-4eca-b83c-89611b0b5dac
# ╠═ed434ee4-6288-4d62-b3a0-bec7fd3ffad8
# ╠═8107c0a6-24c8-4e8d-9092-7b86af925cc2
# ╠═5230514e-4ee0-4bee-bab5-5846b320c466
# ╠═cea55edf-d2ea-4107-89d0-7a2aa760f6df
# ╠═94d401c5-f160-4113-bc2d-a8617d0cc0da
# ╠═ba2f0bdb-97c4-4db7-b614-e7f1a392cdb7
# ╠═d5220791-a3e8-44b6-aae3-cfe2a52838df
# ╠═9996d074-7e32-467e-9b15-e806f9237606
# ╠═a1aa6205-b436-459a-a60e-a566e784457b
# ╠═4ed80ff1-1e52-4f2d-a758-73cd40c75bf0
# ╠═94ad5025-2c74-4970-9ce7-28ec698e254a
# ╠═fa33e295-cb79-4107-86ae-6524938bcf96
# ╠═9d4c20c5-61ad-49ba-ac18-64db25038d1e
# ╠═813b929c-ccd8-4ee7-a496-87647444e9ac
# ╠═fd50c12c-3b8b-4168-853d-e7236478cba5
# ╠═b4b95e73-a1a5-4491-bbe0-18c0cf70ccbf
# ╠═57b206a6-3e9a-44af-b063-de358b597cab
# ╠═e092697c-cf33-4e7a-ba31-32d0f74576e8
# ╠═3e85daec-fbe8-4b4c-9b01-37d345e6b762
# ╠═c4590c18-7b89-47b4-9788-5959099d8096
# ╠═6e6c5d71-4a1a-4c31-9a23-96d99ac5f1f0
# ╠═b03ab3c2-598b-440e-a00d-a93bc483828f
# ╠═9bf77bda-92dc-43c5-aaa0-711e5f127868
# ╠═86f392b6-9332-43b6-b09d-aecbdf72933a
# ╠═45fc61c9-1173-4f59-a883-974e1029edc2
# ╠═7c31bf02-dcd3-4880-85eb-0247e646177d
# ╠═7ec6d916-685b-4044-bca6-cca08b51def8

### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ c82bb3e6-e052-11ec-3392-eba4a4b464ba
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

# ╔═╡ 2df207ed-cf9c-4d7d-8354-4da14f93276c
TableOfContents()

# ╔═╡ 672bf97d-e5f8-4b06-bff4-9ae468501c78
md"""
## TODO
- Assess calibrations (1pt, 3pt, 6pt) for simulated data
- Add motion reproducibility (SWCS, ICS, AS) to simulated data
"""

# ╔═╡ 4cf61487-3913-4e11-970e-4ca31a5ffc8d
md"""
## Integrated
"""

# ╔═╡ 57d46998-368a-4916-8380-ee49d5473a49
path_integrated = "/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/output_new/integrated_scoring";

# ╔═╡ 2c960bd8-ae64-453f-b29f-275bf5263774
begin
    df_i = CSV.read(string(path_integrated, "/full2.csv"), DataFrame)
    df_i, df_irest = groupby(df_i, :blur)
    df_i = df_i[!, 2:end]
end;

# ╔═╡ 7f4fae09-7916-4a98-a102-7f861900c457
df_i_low, df_i_normal = groupby(df_i, :DENSITY);

# ╔═╡ c987e55d-9ffd-40a5-9d64-1653d4000837
df_i_low_small, df_i_low_medium, df_i_low_large = groupby(df_i_low, :SIZE);

# ╔═╡ 0c367cb1-d959-4c10-b31c-18d3cb396255
df_i_normal_small, df_i_normal_medium, df_i_normal_large = groupby(df_i_normal, :SIZE);

# ╔═╡ a58beea4-b186-44c6-bf75-98d368f7ebe4
md"""
## Agatston
"""

# ╔═╡ 4fd5b649-cea1-4d1f-95ac-35aa12f8fbff
path_agat = "/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/output_new/agatston";

# ╔═╡ 3fa3748a-22d6-49d8-9888-a749dca99ce2
begin
    df_a = CSV.read(string(path_agat, "/full2.csv"), DataFrame)
    df_a, df_arest = groupby(df_a, :blur)
    df_a = df_a[!, 2:end]
end;

# ╔═╡ c38b9994-7698-4c81-b65b-534b994d8ff0
df_a_low, df_a_normal = groupby(df_a, :DENSITY);

# ╔═╡ 79952ee5-3ccd-4d85-95b4-f7610c90d63d
df_a_low_small, df_a_low_medium, df_a_low_large = groupby(df_a_low, :SIZE);

# ╔═╡ 7223ae32-7241-463f-be13-2d2ce6bc7de7
df_a_normal_small, df_a_normal_medium, df_a_normal_large = groupby(df_a_normal, :SIZE);

# ╔═╡ f39a9ed5-c897-4b90-8494-bf231312970e
md"""
## Spatially Weighted
"""

# ╔═╡ 32623778-9a19-4318-81c6-1fc7aa0157fe
path_swcs = "/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/output_new/swcs";

# ╔═╡ 317f8c26-7b39-45c3-b9d3-02925d4b0514
begin
    df_s = CSV.read(string(path_swcs, "/full2.csv"), DataFrame)
    df_s, df_srest = groupby(df_s, :blur)
    df_s = df_s[!, 2:end]
end;

# ╔═╡ 06deefa7-9170-4e2f-ac49-d6dc65c0af76
df_s_low, df_s_normal = groupby(df_s, :DENSITY);

# ╔═╡ 5a4cc121-fd02-49bf-bdf5-07b71b48ce19
df_s_low_small, df_s_low_medium, df_s_low_large = groupby(df_s_low, :SIZE);

# ╔═╡ d5d8e9da-eba8-4788-b787-95b3a2b41329
df_s_normal_small, df_s_normal_medium, df_s_normal_large = groupby(df_s_normal, :SIZE);

# ╔═╡ 310c1999-77a6-432e-850c-4111411cceb0
md"""
# Figures
"""

# ╔═╡ 3a8c6f3c-0fa6-477c-85e2-7739cf2537ae
md"""
## Zero CAC Scores (FIG)
"""

# ╔═╡ 3c40ca1d-4111-4bf1-941e-30b59cd1a264
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

# ╔═╡ 932e1a11-ec52-48bc-af27-bb3e09e6b27f
md"""
#### SWCS
"""

# ╔═╡ dc8805a5-2b4c-47fe-88c8-0f4434ddd5cf
df_s_80, df_s_100, df_s_120, df_s_135, = groupby(df_s, :scan);

# ╔═╡ 2c78dcc4-03c7-49b5-ac45-30c198671761
begin
    mean_swcs_80, std_swcs_80 = mean(df_s_80[!, :swcs_bkg]), std(df_s_80[!, :swcs_bkg])
    mean_swcs_100, std_swcs_100 = mean(df_s_100[!, :swcs_bkg]), std(df_s_100[!, :swcs_bkg])
    mean_swcs_120, std_swcs_120 = mean(df_s_120[!, :swcs_bkg]), std(df_s_120[!, :swcs_bkg])
    mean_swcs_135, std_swcs_135 = mean(df_s_135[!, :swcs_bkg]), std(df_s_135[!, :swcs_bkg])
end;

# ╔═╡ d3e13dcd-201a-49ab-86f7-a7d3320088c2
begin
    array_s_80 = Array(df_s_80[!, 12:14])
    array_s_100 = Array(df_s_100[!, 12:14])
    array_s_120 = Array(df_s_120[!, 12:14])
    array_s_135 = Array(df_s_135[!, 12:14])
end;

# ╔═╡ 4b97d892-b0db-4f75-9239-9da215f905f0
begin
    num_zeroCAC_80 = length(findall(x -> x < mean_swcs_80, array_s_80))
    num_zeroCAC_100 = length(findall(x -> x < mean_swcs_100, array_s_100))
    num_zeroCAC_120 = length(findall(x -> x < mean_swcs_120, array_s_120))
    num_zeroCAC_135 = length(findall(x -> x < mean_swcs_135, array_s_135))

    total_zero_s = num_zeroCAC_80 + num_zeroCAC_100 + num_zeroCAC_120 + num_zeroCAC_135
    total_cac = length(array_s_80) * 4
end;

# ╔═╡ da2b55dc-2fe9-468f-a53e-baf232ca79cd
total_cac

# ╔═╡ 01d903d8-fca2-4317-a2d7-93429f06d47c
md"""
#### Agatston
"""

# ╔═╡ 97a612d5-44e7-43f1-b61e-8518c7e6aa28
array_a = hcat(df_a[!, 7], df_a[!, 9], df_a[!, 11]);

# ╔═╡ a31dfad2-e420-4c41-b805-227778a39fc9
num_zero_a = length(findall(x -> x == 0, array_a));

# ╔═╡ 86d5f7fd-6a1b-4fdc-9e06-6cd66f8a169c
md"""
#### Integrated
"""

# ╔═╡ aabe8e7b-faed-499d-aea2-fd44836a3cfb
df_i_80, df_i_100, df_i_120, df_i_135, = groupby(df_i, :scan);

# ╔═╡ dd644657-51e9-4b27-a186-6efba0bb83f8
begin
    mean_i_80, std_i_80 = mean(df_i_80[!, :mass_bkg]), std(df_i_80[!, :mass_bkg])
    mean_i_100, std_i_100 = mean(df_i_100[!, :mass_bkg]), std(df_i_100[!, :mass_bkg])
    mean_i_120, std_i_120 = mean(df_i_120[!, :mass_bkg]), std(df_i_120[!, :mass_bkg])
    mean_i_135, std_i_135 = mean(df_i_135[!, :mass_bkg]), std(df_i_135[!, :mass_bkg])
end

# ╔═╡ d79d6478-0837-40f0-a7a5-270bf44d6166
begin
    array_i_80 = hcat(df_i_80[!, 7], df_i_80[!, 9], df_i_80[!, 11])
    array_i_100 = hcat(df_i_100[!, 7], df_i_100[!, 9], df_i_100[!, 11])
    array_i_120 = hcat(df_i_120[!, 7], df_i_120[!, 9], df_i_120[!, 11])
    array_i_135 = hcat(df_i_135[!, 7], df_i_135[!, 9], df_i_135[!, 11])
end;

# ╔═╡ 2b5a8dca-faae-4163-ab14-ae4f98b2a7cc
begin
    num_zeroCAC_80_i = length(findall(x -> x < mean_i_80, array_i_80))
    num_zeroCAC_100_i = length(findall(x -> x < mean_i_100, array_i_100))
    num_zeroCAC_120_i = length(findall(x -> x < mean_i_120, array_i_120))
    num_zeroCAC_135_i = length(findall(x -> x < mean_i_135, array_i_135))

    total_zero_i =
        num_zeroCAC_80_i + num_zeroCAC_100_i + num_zeroCAC_120_i + num_zeroCAC_135_i
end;

# ╔═╡ 01ecd9f2-4aa6-4f9e-b1ec-e51e3aa99571
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
        "/Users/daleblack/Google Drive/Research/Papers/My Papers/cac-simulation/figures/zero_cac.png",
        f,
    )
    return f
end

# ╔═╡ cc88bab2-8b5e-474a-889d-f0fbfec790a2
with_theme(medphys_theme) do
    zero_cac_plot()
end

# ╔═╡ f4b51729-d53a-40e2-919b-5d4562cffbc0
total_zero_i, total_zero_s, num_zero_a

# ╔═╡ 7dfc24a4-e006-45f4-b5b9-977a7c3c0b7c
md"""
## Linear Regression Mass Score (FIG)
"""

# ╔═╡ 33e21324-b4f8-4b5a-b759-dc38057a612d
md"""
#### Normal Density
"""

# ╔═╡ 1ef8e352-7e89-461f-9ac8-8a39cf2c14d7
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

# ╔═╡ ece4a68f-1f3b-4140-bc34-064cd3953763
begin
    newX1 = DataFrame(; X=collect(1:1000))
    pred_i_norm = GLM.predict(model_i_normal, newX1)
end

# ╔═╡ 731edfcf-6c6b-4bbc-b23f-74be29d80513
co1 = coef(model_i_normal)

# ╔═╡ c86c3d16-ace0-4de2-a7fb-267ea2925302
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

# ╔═╡ 36718970-38e9-4497-b9dd-d9947b1717e4
begin
    newX3 = DataFrame(; X=collect(1:1000))
    pred_a_norm = GLM.predict(model_a_normal, newX3)
end

# ╔═╡ e7bb990c-3414-480a-b9f4-1b532bc006c9
co3 = coef(model_a_normal)

# ╔═╡ c93ba92a-ccc7-4ae0-9207-300715821fc5
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
        "/Users/daleblack/Google Drive/Research/Papers/My Papers/cac-simulation/figures/linear_reg_norm.png",
        f,
    )
    return f
end

# ╔═╡ a4bc0a8a-4904-4217-97a9-44158f99ae70
with_theme(medphys_theme) do
    lin_reg_norm()
end

# ╔═╡ cf59cd3f-026c-4321-a1ac-de217177b52e
md"""
#### Low Density
"""

# ╔═╡ 0ce89635-b285-4ad7-a936-134eaac981ca
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

# ╔═╡ 6a163966-b24c-4e6a-b2a1-795bbf54bc36
begin
    newX2 = DataFrame(; X=collect(1:1000))
    pred_i_low = GLM.predict(model_i_low, newX2)
end

# ╔═╡ 27e191ad-b04c-4e48-9763-128637fdb4b6
co2 = coef(model_i_low)

# ╔═╡ d4c89efd-fb4b-4553-8920-847c393cb4bc
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
        rms(data[!, :X], data[!, :Y]), rms(GLM.predict(model_a_low), data[!, :Y])
    ]
end

# ╔═╡ 5ccbe974-afb9-454f-bd1a-88879a507565
begin
    newX4 = DataFrame(; X=collect(1:1000))
    pred_a_low = GLM.predict(model_a_low, newX4)
end

# ╔═╡ d43c412e-22b5-4406-a883-3410641b5b16
co4 = coef(model_a_low)

# ╔═╡ 2e04217d-8dfb-48c9-85dc-9fb42cfd3039
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
    ylims!(ax2; low=-10, high=30)
    ax2.xticks = [0, 5, 10, 15, 20, 25]
    ax2.yticks = [-10, 0, 10, 20, 30]
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
        "/Users/daleblack/Google Drive/Research/Papers/My Papers/cac-simulation/figures/linear_reg_low.png",
        f,
    )
    return f
end

# ╔═╡ 98c3eaef-502e-4a20-b5b6-46a3d6b394d3
with_theme(medphys_theme) do
    lin_reg_low()
end

# ╔═╡ 69a5a71a-b05f-4c71-bfb2-4d085e48c645
md"""
#### SWCS
"""

# ╔═╡ d98ac8e2-0d07-4194-8bd6-2644a27f9b48
array_s = Array(hcat(df_s[!, (end - 3):(end - 1)]));

# ╔═╡ 2a790b87-d76e-4a66-8add-4c776b0f1a93
vec_s = vec(array_s);

# ╔═╡ ea7e4439-7c29-4ae9-8038-0dd760cf9545
md"""
#### Agatston
"""

# ╔═╡ 6e0c83a5-45b2-4966-a7e4-721b2b403ec7
vec_a = vec(array_a);

# ╔═╡ 1abbc7c8-46cd-42ad-bf2a-2f2f4f883f79
md"""
#### Integrated
"""

# ╔═╡ 6372b044-a4c1-4a67-888e-1d473f3d6a87
array_i = hcat(df_i[!, 7], df_i[!, 9], df_i[!, 11]);

# ╔═╡ 1bca72fc-adb9-4d54-8e30-8b23ba4da329
vec_i = vec(array_i);

# ╔═╡ 4fde5c38-784d-4461-ab16-525ffb678a42
md"""
#### Known mass
"""

# ╔═╡ 3e9a9103-5d84-44de-aeab-9adce3adb876
array_k = hcat(df_a[!, 6], df_a[!, 8], df_a[!, 10]);

# ╔═╡ 49a9cc2f-bedd-4992-99d0-08f8001c5ba0
vec_k = vec(array_k);

# ╔═╡ ee93ecc4-a964-4b0d-8673-da2b2f1430b7
md"""
#### Correlations
"""

# ╔═╡ b3d14f1b-d538-415a-bb8d-df9abb36be70
pearson_a = cor(hcat(vec_a, vec_k));

# ╔═╡ e7c525f4-e9da-424e-9483-e2d97832e19d
pearson_i = cor(hcat(vec_i, vec_k));

# ╔═╡ 7ddde041-6e01-4f4c-b930-c2b5143fd51f
pearson_s = cor(hcat(vec_s, vec_k));

# ╔═╡ c05610e1-ecc6-4391-b6cf-964bb4696362
p_a = EqualVarianceTTest(vec_a, vec_k);

# ╔═╡ 8fd6011f-24b4-43a8-809e-371ef33ceedc
p_i = EqualVarianceTTest(vec_i, vec_k);

# ╔═╡ e89f1021-7a7a-4d32-82b3-dd93d2823388
p_s = EqualVarianceTTest(vec_s, vec_k);

# ╔═╡ 870ba06e-5f39-4200-a29d-9a771396c9f8
p_corr_coef = [
    (pearson_a[2], "p = ", pvalue(p_a))
    (pearson_i[2], "p = ", pvalue(p_i))
    (pearson_s[2], "p = ", pvalue(p_s))
];

# ╔═╡ 70dcfdd9-b6c2-4b1c-8b6f-78724050acaf
let
    df = df_a
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
    model = lm(@formula(Y ~ X), data)
    global r_squared_a
    r_squared_a = r2(model)
end;

# ╔═╡ ad15a980-a7a3-4118-9152-6903bf24b605
let
    df = df_i
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
    model = lm(@formula(Y ~ X), data)
    global r_squared_i
    r_squared_i = r2(model)
end;

# ╔═╡ b11aadd1-6124-4155-9ce1-c2f162053bf9
r_squared = [r_squared_a, r_squared_i, "NA"]

# ╔═╡ 636cfdaa-4b53-426b-9465-aa348f4bee7f
RMS_values = [rms(vec_a, vec_k), rms(vec_i, vec_k), "NA"]

# ╔═╡ 4dcbfbdb-d620-438b-b4c5-c288563d63a8
t_test_rmse1 = OneSampleTTest(Float64.(RMS_values[1:2]));

# ╔═╡ 0f59b081-4fcd-43ba-9d79-83af3ae2b3ac
pvalue(t_test_rmse1)

# ╔═╡ 20dc54a6-6351-440f-bc32-9b95e962f970
DataFrame(; Pearson_Corr_Coeff=p_corr_coef, R_Squared=r_squared, RMS_values=RMS_values)

# ╔═╡ 7ec78f71-6b20-4c8c-a9da-a216404bee72
md"""
## Reproducibility (FIG)
"""

# ╔═╡ 88df8b9d-ff41-41d3-99fe-8ab9a050a803
md"""
#### Integrated
"""

# ╔═╡ 53c1b176-e2f7-4cb9-baa9-26d61ab8c18f
path_integrated_r = "/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/output_repeated/integrated_scoring";

# ╔═╡ 89b28ac5-dd69-4812-84fd-64b54606d146
begin
    df_i_r = CSV.read(string(path_integrated_r, "/full2.csv"), DataFrame)
    df_i_r, df_i_r_rest = groupby(df_i_r, :blur)
    df_i_r = df_i_r[!, 2:end]
end;

# ╔═╡ 4f5328e7-66dd-433c-a7f6-5dc461c83a4c
df_i_low_r, df_i_normal_r = groupby(df_i_r, :DENSITY);

# ╔═╡ 6e515ee8-2836-4285-a0a9-131e1e7b4b7f
df_i_low_small_r, df_i_low_medium_r, df_i_low_large_r = groupby(df_i_low_r, :SIZE);

# ╔═╡ 9ae3b4cb-fac5-4168-868d-724de9b0c9d2
df_i_normal_small_r, df_i_normal_medium_r, df_i_normal_large_r = groupby(
    df_i_normal_r, :SIZE
);

# ╔═╡ d9e26f65-6d40-48e5-b567-4cce5a2c530b
array_i_r = hcat(
    df_i_r[!, :calculated_mass_large],
    df_i_r[!, :calculated_mass_medium],
    df_i_r[!, :calculated_mass_small],
)

# ╔═╡ 5273506c-7380-45cc-ac3d-bd48998fe6aa
begin
    mean_bkg_i = mean(df_i[!, :mass_bkg])
    mean_bkg_i_r = mean(df_i_r[!, :mass_bkg])
end

# ╔═╡ 27f1aaa0-7a3e-4cb1-abc8-7a52105c994b
begin
    idxs_i_r_large = Tuple(findall(x -> x < mean_bkg_i_r, array_i_r[:, 1]))
    idxs_i_large = Tuple(findall(x -> x < mean_bkg_i, array_i[:, 1]))
    indxs_i_large_tot = (idxs_i_r_large..., idxs_i_large...)
    indx_i_large_tot = Tuple(unique(indxs_i_large_tot))
end

# ╔═╡ 8e2a5e37-a32a-418a-a392-465ba617bc03
begin
    idxs_i_r_med = Tuple(findall(x -> x < mean_bkg_i_r, array_i_r[:, 2]))
    idxs_i_med = Tuple(findall(x -> x < mean_bkg_i, array_i[:, 2]))
    indxs_i_med_tot = (idxs_i_r_med..., idxs_i_med...)
    indx_i_med_tot = Tuple(unique(indxs_i_med_tot))
end

# ╔═╡ 1cf7d56a-a350-458a-8466-4f6b821ea44f
begin
    idxs_i_r_small = Tuple(findall(x -> x < mean_bkg_i_r, array_i_r[:, 3]))
    idxs_i_small = Tuple(findall(x -> x < mean_bkg_i, array_i[:, 3]))
    indxs_i_small_tot = (idxs_i_r_small..., idxs_i_small...)
    indx_i_small_tot = Tuple(unique(indxs_i_small_tot))
end

# ╔═╡ 8d758456-4fd9-4fd8-be4e-ea8d1dfdf020
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

# ╔═╡ 4b82a62c-9f2c-4e61-b4de-f61fa8292dad
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

# ╔═╡ 5fb6004b-e9b4-4da5-961b-51170ab8e57e
begin
    newX1r = DataFrame(; X=collect(1:1000))
    pred_ir = GLM.predict(model_ir, newX1r)
end

# ╔═╡ 437ae6df-aeb3-457e-ae3a-3c3701e53f34
co1r = coef(model_ir)

# ╔═╡ 8eaec985-641a-498f-9cd9-c2c85d877142
md"""
#### Agatston
"""

# ╔═╡ fcd83805-06e2-4cda-aa56-0c13c69424d8
path_agat_r = "/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/output_repeated/agatston";

# ╔═╡ 5f2eb2d2-84c4-4160-a92a-f181b4126450
begin
    df_a_r = CSV.read(string(path_agat_r, "/full2.csv"), DataFrame)
    df_a_r, df_a_r_rest = groupby(df_a_r, :blur)
    df_a_r = df_a_r[!, 2:end]
end;

# ╔═╡ 869de9cc-6ab0-4e0b-ad1f-59db1fd2f69f
df_a_low_r, df_a_normal_r = groupby(df_a_r, :DENSITY);

# ╔═╡ 43b3ba9b-247c-41cd-968e-3a27736426b0
df_a_low_small_r, df_a_low_medium_r, df_a_low_large_r = groupby(df_a_low_r, :SIZE);

# ╔═╡ b16090b5-a232-4e71-9148-b71296efa999
df_a_normal_small_r, df_a_normal_medium_r, df_a_normal_large_r = groupby(
    df_a_normal_r, :SIZE
);

# ╔═╡ 24034fa9-f9d4-4bd0-a736-488b06403adc
array_a_r = hcat(
    df_a_r[!, :calculated_mass_large],
    df_a_r[!, :calculated_mass_medium],
    df_a_r[!, :calculated_mass_large],
);

# ╔═╡ 0c77926d-f5f6-4d5d-9a68-711dbd963431
begin
    idxs_a_r_large = Tuple(findall(x -> x == 0, array_a_r[:, 1]))
    idxs_a_large = Tuple(findall(x -> x == 0, array_a[:, 1]))
    indxs_a_large_tot = (idxs_a_r_large..., idxs_a_large...)
    indx_a_large_tot = Tuple(unique(indxs_a_large_tot))
end

# ╔═╡ 07270ee9-5d00-4600-875e-7be2af3cd889
begin
    idxs_a_r_med = Tuple(findall(x -> x == 0, array_a_r[:, 2]))
    idxs_a_med = Tuple(findall(x -> x == 0, array_a[:, 2]))
    indxs_a_med_tot = (idxs_a_r_med..., idxs_a_med...)
    indx_a_med_tot = Tuple(unique(indxs_a_med_tot))
end

# ╔═╡ aea5b2cb-dc21-41a5-804d-54136c79b7c2
begin
    idxs_a_r_small = Tuple(findall(x -> x == 0, array_a_r[:, 3]))
    idxs_a_small = Tuple(findall(x -> x == 0, array_a[:, 3]))
    indxs_a_small_tot = (idxs_a_r_small..., idxs_a_small...)
    indx_a_small_tot = Tuple(unique(indxs_a_small_tot))
end

# ╔═╡ be382251-024f-49fb-a6bb-901a441fb5a2
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

# ╔═╡ a3996a98-d95f-4a51-8776-d6286f5c8ed5
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

# ╔═╡ 664e7901-31bc-4c49-8d6a-e770e0553c31
begin
    newX2r = DataFrame(; X=collect(1:1000))
    pred_ar = GLM.predict(model_ar, newX2r)
end

# ╔═╡ f278a992-1527-466c-909b-8ff83fed0a61
co2r = coef(model_ar)

# ╔═╡ 5bc912d4-be86-4cca-9e8b-8af31d269edf
md"""
#### SWCS
"""

# ╔═╡ 594b8d0e-f65c-4fe3-b571-e2ee68a9ab5f
path_swcs_r = "/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/output_repeated/swcs";

# ╔═╡ 1d208fd8-b847-4afe-84b5-3a553f50a858
begin
    df_s_r = CSV.read(string(path_swcs_r, "/full2.csv"), DataFrame)
    df_s_r, df_s_r_rest = groupby(df_s_r, :blur)
    df_s_r = df_s_r[!, 2:end]
end;

# ╔═╡ a9f68d7e-7818-47d5-ba1f-965634225b30
df_s_low_r, df_s_normal_r = groupby(df_s_r, :DENSITY);

# ╔═╡ f2177e1f-2a55-4f81-9a69-ab985ae3b7c2
df_s_low_small_r, df_s_low_medium_r, df_s_low_large_r = groupby(df_s_low_r, :SIZE);

# ╔═╡ 35bc90e2-9d70-4daf-a825-b462831c5bf6
df_s_normal_small_r, df_s_normal_medium_r, df_s_normal_large_r = groupby(
    df_s_normal_r, :SIZE
);

# ╔═╡ 624074ec-572a-45ab-932a-aed7edb1f846
array_s_r = hcat(
    df_s_r[!, :calculated_swcs_large],
    df_s_r[!, :calculated_swcs_medium],
    df_s_r[!, :calculated_swcs_small],
);

# ╔═╡ e8b56e73-f22b-4587-9141-558f5f545a5e
begin
    mean_bkg_s = mean(df_s[!, :swcs_bkg])
    mean_bkg_s_r = mean(df_s_r[!, :swcs_bkg])
end

# ╔═╡ 22d79cce-54c3-4df3-a269-78584a3afc7a
begin
    idxs_s_r_large = Tuple(findall(x -> x < mean_bkg_s_r, array_s_r[:, 1]))
    idxs_s_large = Tuple(findall(x -> x < mean_bkg_s, array_s[:, 1]))
    indxs_s_large_tot = (idxs_s_r_large..., idxs_s_large...)
    indx_s_large_tot = Tuple(unique(indxs_s_large_tot))
end

# ╔═╡ 00a151d1-f576-4f80-bcc2-634db90c5841
begin
    idxs_s_r_med = Tuple(findall(x -> x < mean_bkg_s_r, array_s_r[:, 2]))
    idxs_s_med = Tuple(findall(x -> x < mean_bkg_s, array_s[:, 2]))
    indxs_s_med_tot = (idxs_s_r_med..., idxs_s_med...)
    indx_s_med_tot = Tuple(unique(indxs_s_med_tot))
end

# ╔═╡ 8221a487-cb47-48ab-a91a-96e6cb911a52
begin
    idxs_s_r_small = Tuple(findall(x -> x < mean_bkg_s_r, array_s_r[:, 3]))
    idxs_s_small = Tuple(findall(x -> x < mean_bkg_s, array_s[:, 3]))
    indxs_s_small_tot = (idxs_s_r_small..., idxs_s_small...)
    indx_s_small_tot = Tuple(unique(indxs_s_small_tot))
end

# ╔═╡ 95652ae8-c26e-4c4f-b2c6-5134fdfa9b02
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

# ╔═╡ a385671b-af3a-480b-bb0c-18d1c13cf943
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

# ╔═╡ fd610b17-74b8-44e1-bcac-5ac910db05e5
begin
    newX3r = DataFrame(; X=collect(1:1000))
    pred_sr = GLM.predict(model_sr, newX3r)
end

# ╔═╡ 8f6d3806-cc69-4eb2-97bb-cf0cd0e0b565
co3r = coef(model_sr)

# ╔═╡ a85f777a-76f2-4c64-9973-ea9dec245600
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
        "/Users/daleblack/Google Drive/Research/Papers/My Papers/cac-simulation/figures/reprod.png",
        f,
    )
    return f
end

# ╔═╡ d2b90d91-4e63-45b0-8273-5231dbf2778e
with_theme(medphys_theme) do
    reprod()
end

# ╔═╡ 3e361d54-19f3-48d0-8e85-7351ec2f3335
# md"""
# ## Motion Blurring (FIG)
# """

# ╔═╡ f00cfcc7-9ecc-4a84-a23d-2a30845be8f8
# function motion_blur()
# 	f = Figure()

# 	##-- A --##
# 	ax1 = Axis(f[1, 1])

# 	df = df_ii
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

# 	df3 = df_aa
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

# ╔═╡ 043dfeac-a7c3-4937-8ba4-49e5749c0ca1
# with_theme(medphys_theme) do
# 	motion_blur()
# end

# ╔═╡ f8c2e340-8aa4-4af6-b4c8-43cb33bef721
# md"""
# #### Integrated
# """

# ╔═╡ 620ad304-cb77-4a7c-afc7-76a49d220cfd
# df_ii = CSV.read(string(path_integrated, "/full2.csv"), DataFrame);

# ╔═╡ 663a9cab-3235-486b-97c2-23a4a3c86c10
# df_ii0, df_ii05, df_ii1, df_ii15, df_ii2 = groupby(df_ii, :blur);

# ╔═╡ 909f4f32-cfe6-4c82-bfa0-cea38a8de5ca
# let
# 	df = df_ii
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

# ╔═╡ 9ffe97f5-624d-459e-a064-be4e45123b85
# begin
# 	newXii = DataFrame(X=collect(1:1000));
# 	pred_ii = GLM.predict(model_ii, newXii)
# end

# ╔═╡ 67f78b69-7202-4128-8918-ee65fac779bb
# coii = coef(model_ii)

# ╔═╡ 56efa0dd-c7dc-4364-a840-26ceb052826c
# md"""
# #### Agatston
# """

# ╔═╡ bd7e30d2-916e-49b5-905e-7caa2d9080b4
# df_aa = CSV.read(string(path_agat, "/full2.csv"), DataFrame);

# ╔═╡ 8e381b33-f2f5-476f-9a9f-d0b7ec1bb774
# let
# 	df = df_aa
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

# ╔═╡ 41946bb8-fe10-4086-831d-3a68878fe944
# begin
# 	newXaa = DataFrame(X=collect(1:1000));
# 	pred_aa = GLM.predict(model_aa, newXaa)
# end

# ╔═╡ ca8f3808-2173-4a6e-b248-bafc0473bbdb
# coaa = coef(model_aa)

# ╔═╡ bd115bb7-065c-41af-b332-4b701bd6120d
# df_aa0, df_aa05, df_aa1, df_aa15, df_aa2 = groupby(df_aa, :blur);

# ╔═╡ a6f0d36e-823d-4198-b4ea-d95ce16ac65c
md"""
# Tables
"""

# ╔═╡ b7064bb1-f21a-4bb6-a395-9f1c73102058
md"""
## Simulation Parameters
"""

# ╔═╡ 21b03bf8-47ba-4915-936a-92c488671bb1
parameters = [
	"Manufacturer",
	"CT System",
	"Reconstruction",
	"Tube Voltage (kVp)",
	"Exposure Small (mR)",
	"Exposure Medium (mR)",
	"Exposure Large (mR)",
	"Slice Thickness (mm)",
	"Matrix Size Small",
	"Mtrix Size Medium",
	"Matrix Size Large",
	"Detector Element Width (mm)",
	"Detector Thickness (mm)",
]

# ╔═╡ e0a66c8a-e1b6-4c70-9996-ec3fe62786f3
simulations = [
	"Canon"
	"Aquilion One Vision"
	"FBP"
	"120"
	"0.9"
	"2"
	"5.4"
	"0.5"
	"640x440"
	"740x540"
	"840x640"
	"0.5"
	"0.5"
]

# ╔═╡ bf3b0e00-5f2b-435b-9717-96052f7e9071
simulation_parameters = DataFrame(
	"Parameter" => parameters,
	"Simulation" => simulations
)

# ╔═╡ 09e62a9b-4f42-4285-88de-3b90a8a4003a
md"""
## Spatially Weighted Calcium Scoring Means & Stds
"""

# ╔═╡ d74aa830-5a48-4535-8579-39ce2dc8142d
kvp = [
	80
	100
	120
	135
]

# ╔═╡ 6eed3e33-f8f1-41f8-96b1-6e38fa129733
means = [
	201.4898
	174.3658
	158.2645
	152.8815
]

# ╔═╡ 36962426-647f-4194-92cc-aed51c628484
stds = [
	34.3038
	28.1708
	22.9656
	21.0778
]

# ╔═╡ 512f3683-3231-4140-96d8-dd2378f5094d
means_stds = DataFrame(
	"kVp" => kvp,
	"Mean" => means,
	"Std" => stds
)

# ╔═╡ 5d7dfbc7-7a69-4acb-93e4-a8ebf0f8b787


# ╔═╡ 235bc94e-166d-4df1-9ed6-f331e6e1cf0c
md"""
## Summaries
"""

# ╔═╡ e9dc2696-d2a8-4c4f-8665-da526e242570
r_vals_rep = [
	r2_1r
	r2_2r
	r2_3r
]

# ╔═╡ 574b38a3-964a-4d7a-adf4-e36198c003a3
rmse_vals_rep = [
	rms_values1r[1]
	rms_values2r[1]
	rms_values3r[1]
]

# ╔═╡ bc0ede74-4122-4d56-8d99-0043013752df
rmsd_vals_rep = [
	rms_values1r[2]
	rms_values2r[2]
	rms_values3r[2]
]

# ╔═╡ 165f2e6f-afd4-4f65-aae7-0023a531dbd0
summ_reprod = DataFrame(
	"R Correlation Coefficient" => r_vals_rep,
	"RMSE" => rmse_vals_rep,
	"RMSD" => rmsd_vals_rep
)

# ╔═╡ 7c314dcf-aa06-491d-bb45-89c5c8006688
r_vals_norm = [
	r2_1,
	r2_3,
]

# ╔═╡ 5daffb8c-9387-4bdf-9fd1-46dd8bd1b873
rmse_vals_norm = [
	rms_values1[1],
	rms_values3[1]
]

# ╔═╡ 9f1b2d37-1d54-4339-869a-a927c66618c4
rmsd_vals_norm = [
	rms_values1[2],
	rms_values3[2]
]

# ╔═╡ a6c8ebe4-e598-4cb7-a2c2-568f9ddd1b3a
summ_regression_norm = DataFrame(
	"R Correlation Coefficient" => r_vals_norm,
	"RMSE" => rmse_vals_norm,
	"RMSD" => rmsd_vals_norm
)

# ╔═╡ 1fc7c060-f348-4be7-8b61-7ccf0412a38b
r_vals_low = [
	r2_2,
	r2_4,
]

# ╔═╡ 4ba8c595-cdcc-4892-ba2a-591f6baf5c2e
rmse_vals_low = [
	rms_values2[1],
	rms_values4[1]
]

# ╔═╡ 1752eeef-ce73-471e-b5b5-f02f498239e1
rmsd_vals_low = [
	rms_values2[2],
	rms_values4[2]
]

# ╔═╡ 1069b0be-5a04-4ef2-bdb0-1016edbca301
summ_regression_low = DataFrame(
	"R Correlation Coefficient" => r_vals_low,
	"RMSE" => rmse_vals_low,
	"RMSD" => rmsd_vals_low
)

# ╔═╡ 8ef3a331-fd8f-4b61-b355-b133f0480432
zero_cac_num = [
	string(total_zero_i, "/", total_cac)
	string(total_zero_s, "/", total_cac)
	string(num_zero_a, "/", total_cac)
]

# ╔═╡ 6cf69eaf-aaf2-4bc1-aa10-e22efae2b175
zero_cac_perc = [
	round(total_zero_i / total_cac * 100, digits=2)
	round(total_zero_s / total_cac * 100, digits=2)
	round(num_zero_a / total_cac * 100, digits=2)
]

# ╔═╡ 38f87d7b-6d22-41de-97da-927a478d8e8b
summ_zero_cac = DataFrame(
	"False Negatives (CAC=0)" => zero_cac_num,
	"Percentage False Negatives (%)" => zero_cac_perc
)

# ╔═╡ Cell order:
# ╠═c82bb3e6-e052-11ec-3392-eba4a4b464ba
# ╠═2df207ed-cf9c-4d7d-8354-4da14f93276c
# ╟─672bf97d-e5f8-4b06-bff4-9ae468501c78
# ╟─4cf61487-3913-4e11-970e-4ca31a5ffc8d
# ╠═57d46998-368a-4916-8380-ee49d5473a49
# ╠═2c960bd8-ae64-453f-b29f-275bf5263774
# ╠═7f4fae09-7916-4a98-a102-7f861900c457
# ╠═c987e55d-9ffd-40a5-9d64-1653d4000837
# ╠═0c367cb1-d959-4c10-b31c-18d3cb396255
# ╟─a58beea4-b186-44c6-bf75-98d368f7ebe4
# ╠═4fd5b649-cea1-4d1f-95ac-35aa12f8fbff
# ╠═3fa3748a-22d6-49d8-9888-a749dca99ce2
# ╠═c38b9994-7698-4c81-b65b-534b994d8ff0
# ╠═79952ee5-3ccd-4d85-95b4-f7610c90d63d
# ╠═7223ae32-7241-463f-be13-2d2ce6bc7de7
# ╟─f39a9ed5-c897-4b90-8494-bf231312970e
# ╠═32623778-9a19-4318-81c6-1fc7aa0157fe
# ╠═317f8c26-7b39-45c3-b9d3-02925d4b0514
# ╠═06deefa7-9170-4e2f-ac49-d6dc65c0af76
# ╠═5a4cc121-fd02-49bf-bdf5-07b71b48ce19
# ╠═d5d8e9da-eba8-4788-b787-95b3a2b41329
# ╟─310c1999-77a6-432e-850c-4111411cceb0
# ╟─3a8c6f3c-0fa6-477c-85e2-7739cf2537ae
# ╟─3c40ca1d-4111-4bf1-941e-30b59cd1a264
# ╟─01ecd9f2-4aa6-4f9e-b1ec-e51e3aa99571
# ╟─cc88bab2-8b5e-474a-889d-f0fbfec790a2
# ╠═f4b51729-d53a-40e2-919b-5d4562cffbc0
# ╠═da2b55dc-2fe9-468f-a53e-baf232ca79cd
# ╟─932e1a11-ec52-48bc-af27-bb3e09e6b27f
# ╠═dc8805a5-2b4c-47fe-88c8-0f4434ddd5cf
# ╠═2c78dcc4-03c7-49b5-ac45-30c198671761
# ╠═d3e13dcd-201a-49ab-86f7-a7d3320088c2
# ╠═4b97d892-b0db-4f75-9239-9da215f905f0
# ╟─01d903d8-fca2-4317-a2d7-93429f06d47c
# ╠═97a612d5-44e7-43f1-b61e-8518c7e6aa28
# ╠═a31dfad2-e420-4c41-b805-227778a39fc9
# ╟─86d5f7fd-6a1b-4fdc-9e06-6cd66f8a169c
# ╠═aabe8e7b-faed-499d-aea2-fd44836a3cfb
# ╠═dd644657-51e9-4b27-a186-6efba0bb83f8
# ╠═d79d6478-0837-40f0-a7a5-270bf44d6166
# ╠═2b5a8dca-faae-4163-ab14-ae4f98b2a7cc
# ╟─7dfc24a4-e006-45f4-b5b9-977a7c3c0b7c
# ╠═c93ba92a-ccc7-4ae0-9207-300715821fc5
# ╠═2e04217d-8dfb-48c9-85dc-9fb42cfd3039
# ╟─a4bc0a8a-4904-4217-97a9-44158f99ae70
# ╟─98c3eaef-502e-4a20-b5b6-46a3d6b394d3
# ╟─33e21324-b4f8-4b5a-b759-dc38057a612d
# ╠═1ef8e352-7e89-461f-9ac8-8a39cf2c14d7
# ╠═ece4a68f-1f3b-4140-bc34-064cd3953763
# ╠═731edfcf-6c6b-4bbc-b23f-74be29d80513
# ╠═c86c3d16-ace0-4de2-a7fb-267ea2925302
# ╠═36718970-38e9-4497-b9dd-d9947b1717e4
# ╠═e7bb990c-3414-480a-b9f4-1b532bc006c9
# ╟─cf59cd3f-026c-4321-a1ac-de217177b52e
# ╠═0ce89635-b285-4ad7-a936-134eaac981ca
# ╠═6a163966-b24c-4e6a-b2a1-795bbf54bc36
# ╠═27e191ad-b04c-4e48-9763-128637fdb4b6
# ╠═d4c89efd-fb4b-4553-8920-847c393cb4bc
# ╠═5ccbe974-afb9-454f-bd1a-88879a507565
# ╠═d43c412e-22b5-4406-a883-3410641b5b16
# ╟─69a5a71a-b05f-4c71-bfb2-4d085e48c645
# ╠═d98ac8e2-0d07-4194-8bd6-2644a27f9b48
# ╠═2a790b87-d76e-4a66-8add-4c776b0f1a93
# ╟─ea7e4439-7c29-4ae9-8038-0dd760cf9545
# ╠═6e0c83a5-45b2-4966-a7e4-721b2b403ec7
# ╟─1abbc7c8-46cd-42ad-bf2a-2f2f4f883f79
# ╠═6372b044-a4c1-4a67-888e-1d473f3d6a87
# ╠═1bca72fc-adb9-4d54-8e30-8b23ba4da329
# ╟─4fde5c38-784d-4461-ab16-525ffb678a42
# ╠═3e9a9103-5d84-44de-aeab-9adce3adb876
# ╠═49a9cc2f-bedd-4992-99d0-08f8001c5ba0
# ╟─ee93ecc4-a964-4b0d-8673-da2b2f1430b7
# ╠═b3d14f1b-d538-415a-bb8d-df9abb36be70
# ╠═e7c525f4-e9da-424e-9483-e2d97832e19d
# ╠═7ddde041-6e01-4f4c-b930-c2b5143fd51f
# ╠═c05610e1-ecc6-4391-b6cf-964bb4696362
# ╠═8fd6011f-24b4-43a8-809e-371ef33ceedc
# ╠═e89f1021-7a7a-4d32-82b3-dd93d2823388
# ╠═870ba06e-5f39-4200-a29d-9a771396c9f8
# ╠═70dcfdd9-b6c2-4b1c-8b6f-78724050acaf
# ╠═ad15a980-a7a3-4118-9152-6903bf24b605
# ╠═b11aadd1-6124-4155-9ce1-c2f162053bf9
# ╠═636cfdaa-4b53-426b-9465-aa348f4bee7f
# ╠═4dcbfbdb-d620-438b-b4c5-c288563d63a8
# ╠═0f59b081-4fcd-43ba-9d79-83af3ae2b3ac
# ╠═20dc54a6-6351-440f-bc32-9b95e962f970
# ╟─7ec78f71-6b20-4c8c-a9da-a216404bee72
# ╟─a85f777a-76f2-4c64-9973-ea9dec245600
# ╟─d2b90d91-4e63-45b0-8273-5231dbf2778e
# ╟─88df8b9d-ff41-41d3-99fe-8ab9a050a803
# ╠═53c1b176-e2f7-4cb9-baa9-26d61ab8c18f
# ╠═89b28ac5-dd69-4812-84fd-64b54606d146
# ╠═4f5328e7-66dd-433c-a7f6-5dc461c83a4c
# ╠═6e515ee8-2836-4285-a0a9-131e1e7b4b7f
# ╠═9ae3b4cb-fac5-4168-868d-724de9b0c9d2
# ╠═d9e26f65-6d40-48e5-b567-4cce5a2c530b
# ╠═5273506c-7380-45cc-ac3d-bd48998fe6aa
# ╠═27f1aaa0-7a3e-4cb1-abc8-7a52105c994b
# ╠═8e2a5e37-a32a-418a-a392-465ba617bc03
# ╠═1cf7d56a-a350-458a-8466-4f6b821ea44f
# ╠═8d758456-4fd9-4fd8-be4e-ea8d1dfdf020
# ╠═4b82a62c-9f2c-4e61-b4de-f61fa8292dad
# ╠═5fb6004b-e9b4-4da5-961b-51170ab8e57e
# ╠═437ae6df-aeb3-457e-ae3a-3c3701e53f34
# ╟─8eaec985-641a-498f-9cd9-c2c85d877142
# ╠═fcd83805-06e2-4cda-aa56-0c13c69424d8
# ╠═5f2eb2d2-84c4-4160-a92a-f181b4126450
# ╠═869de9cc-6ab0-4e0b-ad1f-59db1fd2f69f
# ╠═43b3ba9b-247c-41cd-968e-3a27736426b0
# ╠═b16090b5-a232-4e71-9148-b71296efa999
# ╠═24034fa9-f9d4-4bd0-a736-488b06403adc
# ╠═0c77926d-f5f6-4d5d-9a68-711dbd963431
# ╠═07270ee9-5d00-4600-875e-7be2af3cd889
# ╠═aea5b2cb-dc21-41a5-804d-54136c79b7c2
# ╠═be382251-024f-49fb-a6bb-901a441fb5a2
# ╠═a3996a98-d95f-4a51-8776-d6286f5c8ed5
# ╠═664e7901-31bc-4c49-8d6a-e770e0553c31
# ╠═f278a992-1527-466c-909b-8ff83fed0a61
# ╟─5bc912d4-be86-4cca-9e8b-8af31d269edf
# ╠═594b8d0e-f65c-4fe3-b571-e2ee68a9ab5f
# ╠═1d208fd8-b847-4afe-84b5-3a553f50a858
# ╠═a9f68d7e-7818-47d5-ba1f-965634225b30
# ╠═f2177e1f-2a55-4f81-9a69-ab985ae3b7c2
# ╠═35bc90e2-9d70-4daf-a825-b462831c5bf6
# ╠═624074ec-572a-45ab-932a-aed7edb1f846
# ╠═e8b56e73-f22b-4587-9141-558f5f545a5e
# ╠═22d79cce-54c3-4df3-a269-78584a3afc7a
# ╠═00a151d1-f576-4f80-bcc2-634db90c5841
# ╠═8221a487-cb47-48ab-a91a-96e6cb911a52
# ╠═95652ae8-c26e-4c4f-b2c6-5134fdfa9b02
# ╠═a385671b-af3a-480b-bb0c-18d1c13cf943
# ╠═fd610b17-74b8-44e1-bcac-5ac910db05e5
# ╠═8f6d3806-cc69-4eb2-97bb-cf0cd0e0b565
# ╠═3e361d54-19f3-48d0-8e85-7351ec2f3335
# ╟─f00cfcc7-9ecc-4a84-a23d-2a30845be8f8
# ╠═043dfeac-a7c3-4937-8ba4-49e5749c0ca1
# ╠═f8c2e340-8aa4-4af6-b4c8-43cb33bef721
# ╠═620ad304-cb77-4a7c-afc7-76a49d220cfd
# ╠═663a9cab-3235-486b-97c2-23a4a3c86c10
# ╠═909f4f32-cfe6-4c82-bfa0-cea38a8de5ca
# ╠═9ffe97f5-624d-459e-a064-be4e45123b85
# ╠═67f78b69-7202-4128-8918-ee65fac779bb
# ╠═56efa0dd-c7dc-4364-a840-26ceb052826c
# ╠═bd7e30d2-916e-49b5-905e-7caa2d9080b4
# ╠═8e381b33-f2f5-476f-9a9f-d0b7ec1bb774
# ╠═41946bb8-fe10-4086-831d-3a68878fe944
# ╠═ca8f3808-2173-4a6e-b248-bafc0473bbdb
# ╠═bd115bb7-065c-41af-b332-4b701bd6120d
# ╟─a6f0d36e-823d-4198-b4ea-d95ce16ac65c
# ╟─b7064bb1-f21a-4bb6-a395-9f1c73102058
# ╟─21b03bf8-47ba-4915-936a-92c488671bb1
# ╟─e0a66c8a-e1b6-4c70-9996-ec3fe62786f3
# ╟─bf3b0e00-5f2b-435b-9717-96052f7e9071
# ╟─09e62a9b-4f42-4285-88de-3b90a8a4003a
# ╟─d74aa830-5a48-4535-8579-39ce2dc8142d
# ╟─6eed3e33-f8f1-41f8-96b1-6e38fa129733
# ╟─36962426-647f-4194-92cc-aed51c628484
# ╟─512f3683-3231-4140-96d8-dd2378f5094d
# ╠═5d7dfbc7-7a69-4acb-93e4-a8ebf0f8b787
# ╟─235bc94e-166d-4df1-9ed6-f331e6e1cf0c
# ╟─e9dc2696-d2a8-4c4f-8665-da526e242570
# ╟─574b38a3-964a-4d7a-adf4-e36198c003a3
# ╟─bc0ede74-4122-4d56-8d99-0043013752df
# ╟─165f2e6f-afd4-4f65-aae7-0023a531dbd0
# ╟─7c314dcf-aa06-491d-bb45-89c5c8006688
# ╟─5daffb8c-9387-4bdf-9fd1-46dd8bd1b873
# ╟─9f1b2d37-1d54-4339-869a-a927c66618c4
# ╠═a6c8ebe4-e598-4cb7-a2c2-568f9ddd1b3a
# ╟─1fc7c060-f348-4be7-8b61-7ccf0412a38b
# ╟─4ba8c595-cdcc-4892-ba2a-591f6baf5c2e
# ╟─1752eeef-ce73-471e-b5b5-f02f498239e1
# ╠═1069b0be-5a04-4ef2-bdb0-1016edbca301
# ╟─8ef3a331-fd8f-4b61-b355-b133f0480432
# ╟─6cf69eaf-aaf2-4bc1-aa10-e22efae2b175
# ╠═38f87d7b-6d22-41de-97da-927a478d8e8b

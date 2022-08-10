### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ 88588807-9e94-4845-a92a-c45043c8863f
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

# ╔═╡ 90c88c16-b5a2-44ae-b131-f84b77d0ad33
TableOfContents()

# ╔═╡ cafdaf2f-565e-416c-8ec8-5a5aff875f72
md"""
## Integrated
"""

# ╔═╡ ce61c241-be65-4884-938b-64c176a4d67e
path_integrated = "/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/output_new/integrated_scoring";

# ╔═╡ 29737c66-02db-4190-92e9-27b143460fac
begin
    df_i = CSV.read(string(path_integrated, "/full_calibrations.csv"), DataFrame)
    df_i1, df_i3, df_i6, df_ispec = groupby(df_i, :cal)
end;

# ╔═╡ 569443ef-eaf9-412a-b6e1-df242bdca4db
md"""
## Linear Regression
"""

# ╔═╡ eb186bb7-79f7-4a61-9cec-d0c94a9f3488
md"""
## Linear Regression Mass Score (FIG)
"""

# ╔═╡ b1093c89-449a-47c9-b793-b0e3c0dd08b7
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

# ╔═╡ 5e19343c-56f7-4ead-a2fa-ac25cc41bd6a
md"""
#### 1 Point
"""

# ╔═╡ 0bb8427e-53d3-4765-b361-1e0b747a94c0
let
    df = df_i1
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

# ╔═╡ 6f35362b-f441-4fab-9778-dd5d97fe2e84
begin
    newX1 = DataFrame(; X=collect(1:1000))
    pred_i_norm = GLM.predict(model_i_normal, newX1)
end

# ╔═╡ a09a334a-f175-408a-9282-184661053f82
co1 = coef(model_i_normal)

# ╔═╡ 03c80fba-7061-48e3-9f08-920f9a6d2bec
md"""
#### 3 Point
"""

# ╔═╡ 7c919af2-3043-498d-b96e-872fe0e18b54
let
    df = df_i3
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

# ╔═╡ 714429cc-8c2a-4c2d-9250-af4a05e8a2c8
begin
    newX3 = DataFrame(; X=collect(1:1000))
    pred_a_norm = GLM.predict(model_a_normal, newX3)
end

# ╔═╡ 1678d74c-0362-4535-8629-9a560c57cf21
co3 = coef(model_a_normal)

# ╔═╡ e50909ca-fae2-4914-b5de-71553b7e04b0
function lin_reg_norm()
    f = Figure()

    ##-- A --##
    ax1 = Axis(f[1, 1])

    df = df_i1
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
    ax1.title = "Integrated (1 Point)"

    ##-- B --##
    ax2 = Axis(f[2, 1])

    df3 = df_i3
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
    ax2.title = "Integrated (3 Point)"

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
        "/Users/daleblack/Google Drive/Research/Papers/My Papers/cac-simulation/figures/calibration_compare_1_3.png",
        f,
    )
    return f
end

# ╔═╡ 80133408-248f-4ea5-9f1f-d93a6dd6ec23
with_theme(medphys_theme) do
    lin_reg_norm()
end

# ╔═╡ fddd9210-9b8a-463f-99a5-6148ee2ead8b
md"""
#### 6 Point
"""

# ╔═╡ a7b409b8-0d0e-4e4d-aa69-5d72ab34e2c4
let
    df = df_i6
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

# ╔═╡ 7ea608e1-36c5-48b4-87f7-71ddcafe118c
begin
    newX2 = DataFrame(; X=collect(1:1000))
    pred_i_low = GLM.predict(model_i_low, newX2)
end

# ╔═╡ d2a98d9a-ad56-4d94-aa0f-079987f94854
co2 = coef(model_i_low)

# ╔═╡ df760f9e-8f84-453a-92fe-b717307c40ab
md"""
#### Specific
"""

# ╔═╡ 7645093a-e3c8-417b-91f3-e42ffca00e76
let
    df = df_ispec
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

# ╔═╡ 0c9d689b-5096-43c0-a1f0-df62e0ea261d
begin
    newX4 = DataFrame(; X=collect(1:1000))
    pred_a_low = GLM.predict(model_a_low, newX4)
end

# ╔═╡ 3f34747c-9c34-4408-b4fb-2bff5b9ebc0f
co4 = coef(model_a_low)

# ╔═╡ 1ad0639b-7481-4221-89a1-4c77070cb81e
function lin_reg()
    f = Figure()

    ##-- A --##
    ax1 = Axis(f[1, 1])

    df = df_i1
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
    ax1.title = "1 Point"

    ##-- B --##
    ax2 = Axis(f[2, 1])

    df3 = df_i3
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
    ax2.title = "3 Points"

	##-- C --##
	ax3 = Axis(f[1, 2])
	
	df2 = df_i6
	sc1 = scatter!(ax3, df2[!, :ground_truth_mass_large], df2[!, :calculated_mass_large])
	errorbars!(
	    ax3,
	    df2[!, :ground_truth_mass_large],
	    df2[!, :calculated_mass_large],
	    rms(df2[!, :ground_truth_mass_large], df2[!, :calculated_mass_large]),
	)
	sc2 = scatter!(ax3, df2[!, :ground_truth_mass_medium], df2[!, :calculated_mass_medium])
	errorbars!(
	    ax3,
	    df2[!, :ground_truth_mass_medium],
	    df2[!, :calculated_mass_medium],
	    rms(df2[!, :ground_truth_mass_medium], df2[!, :calculated_mass_medium]),
	)
	sc3 = scatter!(
	    ax3, df2[!, :ground_truth_mass_small], df2[!, :calculated_mass_small]; color=:red
	)
	errorbars!(
	    ax3,
	    df2[!, :ground_truth_mass_small],
	    df2[!, :calculated_mass_small],
	    rms(df2[!, :ground_truth_mass_small], df2[!, :calculated_mass_small]),
	)
	ln1 = lines!(ax3, [-1000, 1000], [-1000, 1000])
	ln2 = lines!(ax3, collect(1:1000), pred_i_low; linestyle=:dashdot)
	Textbox(
	    f[1, 2];
	    placeholder="y = $(trunc(co2[2]; digits=3))x + $(trunc(co2[1]; digits=3)) \nr = $(trunc(r2_2; digits=3)) \nRMSE: $(trunc(rms_values2[1]; digits=3)) \nRMSD: $(trunc(rms_values2[2]; digits=3))",
	    tellheight=false,
	    tellwidth=false,
	    boxcolor=:white,
	    halign=:left,
	    valign=:top,
	    textsize=12,
	)
	
	xlims!(ax3; low=0, high=25)
	ylims!(ax3; low=0, high=25)
	ax3.xticks = [0, 5, 10, 15, 20, 25]
	ax3.yticks = [0, 5, 10, 15, 20, 25]
	ax3.xlabel = "Known Mass (mg)"
	ax3.ylabel = "Calculated Mass (mg)"
	ax3.title = "6 Points"
	# hidedecorations!(ax3, ticklabels=false, ticks=false, label=false)

	##-- D --##
	ax4 = Axis(f[2, 2])
	
	df4 = df_ispec
	sc1 = scatter!(ax4, df4[!, :ground_truth_mass_large], df4[!, :calculated_mass_large])
	errorbars!(
	    ax4,
	    df4[!, :ground_truth_mass_large],
	    df4[!, :calculated_mass_large],
	    rms(df4[!, :ground_truth_mass_large], df4[!, :calculated_mass_large]),
	)
	sc2 = scatter!(ax4, df4[!, :ground_truth_mass_medium], df4[!, :calculated_mass_medium])
	errorbars!(
	    ax4,
	    df4[!, :ground_truth_mass_medium],
	    df4[!, :calculated_mass_medium],
	    rms(df4[!, :ground_truth_mass_medium], df4[!, :calculated_mass_medium]),
	)
	sc3 = scatter!(
	    ax4, df4[!, :ground_truth_mass_small], df4[!, :calculated_mass_small]; color=:red
	)
	errorbars!(
	    ax4,
	    df4[!, :ground_truth_mass_small],
	    df4[!, :calculated_mass_small],
	    rms(df4[!, :ground_truth_mass_small], df4[!, :calculated_mass_small]),
	)
	ln1 = lines!(ax4, [-1000, 1000], [-1000, 1000])
	ln2 = lines!(ax4, collect(1:1000), pred_a_low; linestyle=:dashdot)
	Textbox(
	    f[2, 2];
	    placeholder="y = $(trunc(co4[2]; digits=3))x + $(trunc(co4[1]; digits=3)) \nr = $(trunc(r2_4; digits=3)) \nRMSE: $(trunc(rms_values4[1]; digits=3)) \nRMSD: $(trunc(rms_values4[2]; digits=3))",
	    tellheight=false,
	    tellwidth=false,
	    boxcolor=:white,
	    halign=:left,
	    valign=:top,
	    textsize=12,
	)
	
	xlims!(ax4; low=0, high=25)
	ylims!(ax4; low=0, high=25)
	ax4.xticks = [0, 5, 10, 15, 20, 25]
	ax4.yticks = [0, 5, 10, 15, 20, 25]
	ax4.xlabel = "Known Mass (mg)"
	ax4.ylabel = "Calculated Mass (mg)"
	ax4.title = "Scan Specific"
	# hidedecorations!(ax4, ticklabels=false, ticks=false, label=false)

    ##-- LABELS --##
    f[1:2, 3] = Legend(
        f,
        [sc1, sc2, sc3, ln1, ln2],
        ["Large Inserts", "Medium Inserts", "Small Inserts", "Unity", "Fitted Line"];
        framevisible=false,
    )

    for (label, layout) in zip(["A", "B", "C", "D"], [f[1, 1], f[2, 1], f[1, 2], f[2, 2]])
        Label(
            layout[1, 1, TopLeft()],
            label;
            textsize=25,
            padding=(0, 60, 25, 0),
            halign=:right,
        )
    end

    save(
        "/Users/daleblack/Google Drive/Research/Papers/My Papers/cac-simulation/figures/calibration_comparisons.png",
        f,
    )
    return f
end

# ╔═╡ d1b68d52-cd19-4b18-bc84-f9e672c38797
with_theme(medphys_theme) do
    lin_reg()
end

# ╔═╡ 7132d524-6901-4f7b-aec2-7d18e5e90ca2
function lin_reg_low()
    f = Figure()
    ##-- A --##
    ax1 = Axis(f[1, 1])

    df2 = df_i6
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
    ax1.title = "Integrated (6 Point)"
    # hidedecorations!(ax1, ticklabels=false, ticks=false, label=false)

    ##-- B --##
    ax2 = Axis(f[2, 1])

    df4 = df_ispec
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
    ax2.title = "Integrated (Scan Specific)"
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
        "/Users/daleblack/Google Drive/Research/Papers/My Papers/cac-simulation/figures/calibration_compare_6_spec.png",
        f,
    )
    return f
end

# ╔═╡ 48340c73-054a-449a-95a7-f197e09f2836
with_theme(medphys_theme) do
    lin_reg_low()
end

# ╔═╡ Cell order:
# ╠═88588807-9e94-4845-a92a-c45043c8863f
# ╠═90c88c16-b5a2-44ae-b131-f84b77d0ad33
# ╟─cafdaf2f-565e-416c-8ec8-5a5aff875f72
# ╠═ce61c241-be65-4884-938b-64c176a4d67e
# ╠═29737c66-02db-4190-92e9-27b143460fac
# ╟─569443ef-eaf9-412a-b6e1-df242bdca4db
# ╟─1ad0639b-7481-4221-89a1-4c77070cb81e
# ╟─d1b68d52-cd19-4b18-bc84-f9e672c38797
# ╟─eb186bb7-79f7-4a61-9cec-d0c94a9f3488
# ╟─b1093c89-449a-47c9-b793-b0e3c0dd08b7
# ╟─e50909ca-fae2-4914-b5de-71553b7e04b0
# ╠═7132d524-6901-4f7b-aec2-7d18e5e90ca2
# ╟─80133408-248f-4ea5-9f1f-d93a6dd6ec23
# ╠═48340c73-054a-449a-95a7-f197e09f2836
# ╟─5e19343c-56f7-4ead-a2fa-ac25cc41bd6a
# ╠═0bb8427e-53d3-4765-b361-1e0b747a94c0
# ╠═6f35362b-f441-4fab-9778-dd5d97fe2e84
# ╠═a09a334a-f175-408a-9282-184661053f82
# ╟─03c80fba-7061-48e3-9f08-920f9a6d2bec
# ╠═7c919af2-3043-498d-b96e-872fe0e18b54
# ╠═714429cc-8c2a-4c2d-9250-af4a05e8a2c8
# ╠═1678d74c-0362-4535-8629-9a560c57cf21
# ╟─fddd9210-9b8a-463f-99a5-6148ee2ead8b
# ╠═a7b409b8-0d0e-4e4d-aa69-5d72ab34e2c4
# ╠═7ea608e1-36c5-48b4-87f7-71ddcafe118c
# ╠═d2a98d9a-ad56-4d94-aa0f-079987f94854
# ╟─df760f9e-8f84-453a-92fe-b717307c40ab
# ╠═7645093a-e3c8-417b-91f3-e42ffca00e76
# ╠═0c9d689b-5096-43c0-a1f0-df62e0ea261d
# ╠═3f34747c-9c34-4408-b4fb-2bff5b9ebc0f

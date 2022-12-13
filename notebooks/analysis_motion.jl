### A Pluto.jl notebook ###
# v0.19.16

using Markdown
using InteractiveUtils

# ╔═╡ c82bb3e6-e052-11ec-3392-eba4a4b464ba
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.activate(".")

    using PlutoUI, Statistics, CSV, DataFrames, GLM, CairoMakie, HypothesisTests, Colors, MLJBase
	using StatsBase: quantile!, rmsd
end

# ╔═╡ 2df207ed-cf9c-4d7d-8354-4da14f93276c
TableOfContents()

# ╔═╡ d72c94c2-a208-430c-9d76-85d8cfb33a22
md"""
#### Helper Functions
"""

# ╔═╡ cabe7e9a-e932-406b-95dd-2c9128decdc7
function prepare_linear_regression(df)
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
    r_squared = GLM.r2(model)
    rms_values = [rms(data[!, :X], data[!, :Y]), rmsd(data[!, :Y], GLM.predict(model))]

	x = DataFrame(; X=collect(1:1000))
    fitted_line = GLM.predict(model, x)
	coefficient = coef(model)

	return r_squared, rms_values, fitted_line, coefficient
end

# ╔═╡ 4a255a58-b9a0-4175-a750-b2562361631d
function prepare_linear_regression(df_large, df_medium, df_small, df_reprod_large, df_reprod_medium, df_reprod_small)
    r_array = [
        Tuple(df_reprod_large[!, :calculated_mass_large])...,
        Tuple(df_reprod_medium[!, :calculated_mass_medium])...,
        Tuple(df_reprod_small[!, :calculated_mass_small])...,
    ]
    calc_array = [
        Tuple(df_large[!, :calculated_mass_large])...,
        Tuple(df_medium[!, :calculated_mass_medium])...,
        Tuple(df_small[!, :calculated_mass_small])...,
    ]
    data = DataFrame(; X=r_array, Y=calc_array)
    model = lm(@formula(Y ~ X), data)
    r_squared = GLM.r2(model)
    rms_values = [rms(data[!, :X], data[!, :Y]), rmsd(data[!, :Y], GLM.predict(model))]

	x = DataFrame(; X=collect(1:1000))
    fitted_line = GLM.predict(model, x)
	coefficient = coef(model)

	return r_squared, rms_values, fitted_line, coefficient
end

# ╔═╡ 7fc65de7-7514-4532-803e-1120150fd2fb
function remove_false_negatives(df, array, df_reprod, array_reprod)
	
	idxs_r_large = Tuple(findall(x -> x <= 0, array_reprod[:, 1]))
    idxs_large = Tuple(findall(x -> x <= 0, array[:, 1]))
    indxs_large_tot = Tuple(unique([idxs_r_large..., idxs_large...]))

	idxs_r_med = Tuple(findall(x -> x <= 0, array_reprod[:, 2]))
    idxs_med = Tuple(findall(x -> x <= 0, array[:, 2]))
    indxs_med_tot = Tuple(unique([idxs_r_med..., idxs_med...]))

	idxs_r_small = Tuple(findall(x -> x <= 0, array_reprod[:, 3]))
    idxs_small = Tuple(findall(x -> x <= 0, array[:, 3]))
    indxs_small_tot = Tuple(unique([idxs_r_small..., idxs_small...]))

	if length(indxs_large_tot) > 0
        df_large = df[Not(indxs_large_tot...), :]
        df_r_large = df_reprod[Not(indxs_large_tot...), :]
	else
		df_large = copy(df)
		df_r_large = copy(df_reprod)
    end
    if length(indxs_med_tot) > 0
        df_med = df[Not(indxs_med_tot...), :]
        df_r_med = df_reprod[Not(indxs_med_tot...), :]
    else
		df_med = copy(df)
		df_r_med = copy(df_reprod)
    end
    if length(indxs_small_tot) > 0
        df_small = df[Not(indxs_small_tot...), :]
        df_r_small = df_reprod[Not(indxs_small_tot...), :]
	else
		df_small = copy(df)
		df_r_small= copy(df_reprod)
    end

	return df_large, df_r_large, df_med, df_r_med, df_small, df_r_small
end

# ╔═╡ 043dec5f-07a8-472d-8082-c7a771f93270
FIGURE_PATH = "motion"

# ╔═╡ e79c34d0-bf06-4a7a-93b7-e49f09ce3f4e
md"""
# Load CSVs
"""

# ╔═╡ 4cf61487-3913-4e11-970e-4ca31a5ffc8d
md"""
## Integrated
"""

# ╔═╡ 57d46998-368a-4916-8380-ee49d5473a49
path_integrated = "/Users/daleblack/Google Drive/dev/MolloiLab/cac-simulation/output_new/integrated";

# ╔═╡ 2c960bd8-ae64-453f-b29f-275bf5263774
df_i = CSV.read(string(path_integrated, "/motion.csv"), DataFrame);

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
path_agat = "/Users/daleblack/Google Drive/dev/MolloiLab/cac-simulation/output_new/agatston";

# ╔═╡ 3fa3748a-22d6-49d8-9888-a749dca99ce2
df_a = CSV.read(string(path_agat, "/motion.csv"), DataFrame);

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
path_swcs = "/Users/daleblack/Google Drive/dev/MolloiLab/cac-simulation/output_new/swcs";

# ╔═╡ 317f8c26-7b39-45c3-b9d3-02925d4b0514
df_s = CSV.read(string(path_swcs, "/motion.csv"), DataFrame);

# ╔═╡ 06deefa7-9170-4e2f-ac49-d6dc65c0af76
df_s_low, df_s_normal = groupby(df_s, :DENSITY);

# ╔═╡ 5a4cc121-fd02-49bf-bdf5-07b71b48ce19
df_s_low_small, df_s_low_medium, df_s_low_large = groupby(df_s_low, :SIZE);

# ╔═╡ d5d8e9da-eba8-4788-b787-95b3a2b41329
df_s_normal_small, df_s_normal_medium, df_s_normal_large = groupby(df_s_normal, :SIZE);

# ╔═╡ 090bddd5-081a-43c8-b8d1-23fce50331e8
md"""
## Volume Fraction
"""

# ╔═╡ 993c899f-d741-40d5-81e0-ea7fd2871146
path_vf = "/Users/daleblack/Google Drive/dev/MolloiLab/cac-simulation/output_new/volume_fraction";

# ╔═╡ 5e803aad-9d82-4801-aff6-3c1c2f748b02
df_vf = CSV.read(string(path_vf, "/motion.csv"), DataFrame);

# ╔═╡ 5c76c3c7-bb56-4d90-8f0f-7f7d0e271d3c
df_vf_low, df_vf_normal = groupby(df_vf, :DENSITY);

# ╔═╡ 0085693a-157d-42a8-9d35-ec0feec24c9c
df_vf_low_small, df_vf_low_medium, df_vf_low_large = groupby(df_vf_low, :SIZE);

# ╔═╡ a03343e0-a9f4-4a49-8490-1d542bbe42fd
df_vf_normal_small, df_vf_normal_medium, df_vf_normal_large = groupby(df_vf_normal, :SIZE);

# ╔═╡ 310c1999-77a6-432e-850c-4111411cceb0
md"""
# Figures
"""

# ╔═╡ 3c40ca1d-4111-4bf1-941e-30b59cd1a264
medphys_theme = Theme(;
    Axis=(
        backgroundcolor=:white,
        xgridcolor=:gray,
        xgridwidth=0.5,
        xlabelfont=:Helvetica,
        xticklabelfont=:Helvetica,
        xlabelsize=16,
        xticklabelsize=20,
        # xminorticksvisible = true,
        ygridcolor=:gray,
        ygridwidth=0.5,
        ylabelfont=:Helvetica,
        yticklabelfont=:Helvetica,
        ylabelsize=16,
        yticklabelsize=20,
        # yminortickvisible = true,
        bottomsplinecolor=:black,
        leftspinecolor=:black,
        titlefont=:Helvetica,
        titlesize=30,
    ),
)

# ╔═╡ 7dfc24a4-e006-45f4-b5b9-977a7c3c0b7c
md"""
## Accuracy
"""

# ╔═╡ 79a09453-d375-4f1d-9318-887d806dbb14
md"""
#### Normal Density
"""

# ╔═╡ 166718a5-1403-4f1b-b0f7-539fa006997a
r_squared_normal_i, rms_values_normal_i, fitted_line_normal_i, coefficient_normal_i = prepare_linear_regression(df_i_normal)

# ╔═╡ bcedd090-80d5-480a-aad7-97a5fb8131f5
r_squared_normal_vf, rms_values_normal_vf, fitted_line_normal_vf, coefficient_normal_vf = prepare_linear_regression(df_vf_normal)

# ╔═╡ e0175e63-0ba9-4763-8ce6-e1df05966497
r_squared_normal_a, rms_values_normal_a, fitted_line_normal_a, coefficient_normal_a = prepare_linear_regression(df_a_normal)

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
    lines!(ax1, collect(1:1000), fitted_line_normal_i; linestyle=:dashdot)
    # Textbox(
    #     f[1, 1];
    #     placeholder="y = $(trunc(co1[2]; digits=3))x + $(trunc(co1[1]; digits=3)) \nr = $(trunc(r2_1; digits=3)) \nRMSE: $(trunc(rms_values1[1]; digits=3)) \nRMSD: $(trunc(rms_values1[2]; digits=3))",
    #     tellheight=false,
    #     tellwidth=false,
    #     boxcolor=:white,
    #     halign=:left,
    #     valign=:top,
    #     textsize=12,
    # )
	Textbox(
        f[1, 1];
        placeholder="RMSE: $(trunc(rms_values_normal_i[1]; digits=3)) \nRMSD: $(trunc(rms_values_normal_i[2]; digits=3))",
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

    df3 = df_vf_normal
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
    ln2 = lines!(ax2, collect(1:1000), fitted_line_normal_vf; linestyle=:dashdot)
    # Textbox(
    #     f[2, 1];
    #     placeholder="y = $(trunc(co1_vf[2]; digits=3))x + $(trunc(co1_vf[1]; digits=3)) \nr = $(trunc(r2_1_vf; digits=3)) \nRMSE: $(trunc(rms_values1vf[1]; digits=3)) \nRMSD: $(trunc(rms_values1vf[2]; digits=3))",
    #     tellheight=false,
    #     tellwidth=false,
    #     boxcolor=:white,
    #     halign=:left,
    #     valign=:top,
    #     textsize=12,
    # )
	Textbox(
        f[2, 1];
        placeholder="RMSE: $(trunc(rms_values_normal_vf[1]; digits=3)) \nRMSD: $(trunc(rms_values_normal_vf[2]; digits=3))",
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
    ax2.title = "Volume Fraction (Normal-Density)"


    ##-- C --##
    ax2 = Axis(f[3, 1])

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
    ln2 = lines!(ax2, collect(1:1000), fitted_line_normal_a; linestyle=:dashdot)
    # Textbox(
    #     f[3, 1];
    #     placeholder="y = $(trunc(co3[2]; digits=3))x + $(trunc(co3[1]; digits=3)) \nr = $(trunc(r2_3; digits=3)) \nRMSE: $(trunc(rms_values3[1]; digits=3)) \nRMSD: $(trunc(rms_values3[2]; digits=3))",
    #     tellheight=false,
    #     tellwidth=false,
    #     boxcolor=:white,
    #     halign=:left,
    #     valign=:top,
    #     textsize=12,
    # )
	Textbox(
        f[3, 1];
        placeholder="RMSE: $(trunc(rms_values_normal_a[1]; digits=3)) \nRMSD: $(trunc(rms_values_normal_a[2]; digits=3))",
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
    f[2, 2] = Legend(
        f,
        [sc1, sc2, sc3, ln1, ln2],
        ["Large Inserts", "Medium Inserts", "Small Inserts", "Unity", "Fitted Line"];
        framevisible=false,
    )

    for (label, layout) in zip(["A", "B", "C"], [f[1, 1], f[2, 1], f[3, 1]])
        Label(
            layout[1, 1, TopLeft()],
            label;
            textsize=25,
            padding=(0, 90, 25, 0),
            halign=:right,
        )
    end

    save(joinpath(dirname(pwd()),"figures", FIGURE_PATH, "accuracy_normal.png"), f)
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

# ╔═╡ 5d7ea4ec-088c-4f94-bb64-531673836115
r_squared_low_i, rms_values_low_i, fitted_line_low_i, coefficient_low_i = prepare_linear_regression(df_i_low)

# ╔═╡ 1886135c-4ecc-4f27-9a77-7ead70d7735d
r_squared_low_vf, rms_values_low_vf, fitted_line_low_vf, coefficient_low_vf = prepare_linear_regression(df_vf_low)

# ╔═╡ 8b5c75a7-8973-46bf-88bf-0fe65f4d20ce
r_squared_low_a, rms_values_low_a, fitted_line_low_a, coefficient_low_a = prepare_linear_regression(df_a_low)

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
    ln2 = lines!(ax1, collect(1:1000), fitted_line_low_i; linestyle=:dashdot)
    # Textbox(
    #     f[1, 1];
    #     placeholder="y = $(trunc(co2[2]; digits=3))x + $(trunc(co2[1]; digits=3)) \nr = $(trunc(r2_2; digits=3)) \nRMSE: $(trunc(rms_values2[1]; digits=3)) \nRMSD: $(trunc(rms_values2[2]; digits=3))",
    #     tellheight=false,
    #     tellwidth=false,
    #     boxcolor=:white,
    #     halign=:left,
    #     valign=:top,
    #     textsize=12,
    # )
	Textbox(
        f[1, 1];
        placeholder="RMSE: $(trunc(rms_values_low_i[1]; digits=3)) \nRMSD: $(trunc(rms_values_low_i[2]; digits=3))",
        tellheight=false,
        tellwidth=false,
        boxcolor=:white,
        halign=:left,
        valign=:top,
        textsize=12,
    )

    xlims!(ax1; low=0, high=25)
    ylims!(ax1; low=-10, high=30)
    ax1.xticks = [0, 5, 10, 15, 20, 25]
    ax1.yticks = [-10, 0, 10, 20, 30]
    ax1.xlabel = "Known Mass (mg)"
    ax1.ylabel = "Calculated Mass (mg)"
    ax1.title = "Integrated (Low-Density)"
    # hidedecorations!(ax1, ticklabels=false, ticks=false, label=false)

	##-- B --##
    ax2 = Axis(f[2, 1])

    df4 = df_vf_low
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
    ln2 = lines!(ax2, collect(1:1000), fitted_line_low_vf; linestyle=:dashdot)
    # Textbox(
    #     f[2, 1];
    #     placeholder="y = $(trunc(co2_vf[2]; digits=3))x + $(trunc(co2_vf[1]; digits=3)) \nr = $(trunc(r2_2_vf; digits=3)) \nRMSE: $(trunc(rms_values2_vf[1]; digits=3)) \nRMSD: $(trunc(rms_values2_vf[2]; digits=3))",
    #     tellheight=false,
    #     tellwidth=false,
    #     boxcolor=:white,
    #     halign=:left,
    #     valign=:top,
    #     textsize=12,
    # )
	Textbox(
        f[2, 1];
        placeholder="RMSE: $(trunc(rms_values_low_vf[1]; digits=3)) \nRMSD: $(trunc(rms_values_low_vf[2]; digits=3))",
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
    ax2.title = "Volume Fraction (Low-Density)"
    # hidedecorations!(ax2, ticklabels=false, ticks=false, label=false)

    ##-- C --##
    ax2 = Axis(f[3, 1])

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
    ln2 = lines!(ax2, collect(1:1000), fitted_line_low_a; linestyle=:dashdot)
    # Textbox(
    #     f[3, 1];
    #     placeholder="y = $(trunc(co4[2]; digits=3))x + $(trunc(co4[1]; digits=3)) \nr = $(trunc(r2_4; digits=3)) \nRMSE: $(trunc(rms_values4[1]; digits=3)) \nRMSD: $(trunc(rms_values4[2]; digits=3))",
    #     tellheight=false,
    #     tellwidth=false,
    #     boxcolor=:white,
    #     halign=:left,
    #     valign=:top,
    #     textsize=12,
    # )
	Textbox(
        f[3, 1];
        placeholder="RMSE: $(trunc(rms_values_low_a[1]; digits=3)) \nRMSD: $(trunc(rms_values_low_a[2]; digits=3))",
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

    f[2, 2] = Legend(
        f,
        [sc1, sc2, sc3, ln1, ln2],
        ["Large Inserts", "Medium Inserts", "Small Inserts", "Unity", "Fitted Line"];
        framevisible=false,
    )

    for (label, layout) in zip(["A", "B", "C"], [f[1, 1], f[2, 1], f[3, 1]])
        Label(
            layout[1, 1, TopLeft()],
            label;
            textsize=25,
            padding=(0, 90, 25, 0),
            halign=:right,
        )
    end

    save(joinpath(dirname(pwd()),"figures", FIGURE_PATH, "accuracy_low.png"), f)
    return f
end

# ╔═╡ 98c3eaef-502e-4a20-b5b6-46a3d6b394d3
with_theme(medphys_theme) do
    lin_reg_low()
end

# ╔═╡ 7ec78f71-6b20-4c8c-a9da-a216404bee72
md"""
## Reproducibility
"""

# ╔═╡ 88df8b9d-ff41-41d3-99fe-8ab9a050a803
md"""
#### Integrated
"""

# ╔═╡ 53c1b176-e2f7-4cb9-baa9-26d61ab8c18f
path_integrated_r = "/Users/daleblack/Google Drive/dev/MolloiLab/cac-simulation/output_repeated/integrated";

# ╔═╡ 89b28ac5-dd69-4812-84fd-64b54606d146
df_i_r = CSV.read(string(path_integrated_r, "/motion.csv"), DataFrame);

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

# ╔═╡ 4c9feda4-9fc1-42c7-8cac-a17d51627c6d
md"""
#### Volume Fraction
"""

# ╔═╡ e261c635-98d7-4357-a429-a34b73b47bb5
path_volume_fraction_r = "/Users/daleblack/Google Drive/dev/MolloiLab/cac-simulation/output_repeated/volume_fraction";

# ╔═╡ d2e4fa84-1f12-4924-823a-6ea7ab0d3d8e
df_vf_r = CSV.read(string(path_volume_fraction_r, "/motion.csv"), DataFrame);

# ╔═╡ 101b02a5-9e67-460d-9db0-140684bd1e11
df_vf_low_r, df_vf_normal_r = groupby(df_vf_r, :DENSITY);

# ╔═╡ fe275fdc-c833-43ac-a4c7-7a49b2b728f8
df_vf_low_small_r, df_vf_low_medium_r, df_vf_low_large_r = groupby(df_vf_low_r, :SIZE);

# ╔═╡ 8babe305-6bf0-43e3-9fe7-ded195d2addf
df_vf_normal_small_r, df_vf_normal_medium_r, df_vf_normal_large_r = groupby(
    df_vf_normal_r, :SIZE
);

# ╔═╡ ea17e8f9-c06e-439e-82f9-dbfdef759998
array_vf_r = hcat(
    df_vf_r[!, :calculated_mass_large],
    df_vf_r[!, :calculated_mass_medium],
    df_vf_r[!, :calculated_mass_small],
)

# ╔═╡ 8eaec985-641a-498f-9cd9-c2c85d877142
md"""
#### Agatston
"""

# ╔═╡ fcd83805-06e2-4cda-aa56-0c13c69424d8
path_agat_r = "/Users/daleblack/Google Drive/dev/MolloiLab/cac-simulation/output_repeated/agatston";

# ╔═╡ 5f2eb2d2-84c4-4160-a92a-f181b4126450
df_a_r = CSV.read(string(path_agat_r, "/motion.csv"), DataFrame);

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

# ╔═╡ 5bc912d4-be86-4cca-9e8b-8af31d269edf
md"""
#### SWCS
"""

# ╔═╡ 594b8d0e-f65c-4fe3-b571-e2ee68a9ab5f
path_swcs_r = "/Users/daleblack/Google Drive/dev/MolloiLab/cac-simulation/output_repeated/swcs";

# ╔═╡ 1d208fd8-b847-4afe-84b5-3a553f50a858
df_s_r = CSV.read(string(path_swcs_r, "/motion.csv"), DataFrame);

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

# ╔═╡ 1ecb8fd3-a36b-4b06-8daf-e598e297ecf8
function remove_false_negatives(df, array, df_reprod, array_reprod, swcs::Bool)
	if swcs == true
		idxs_r_large = Tuple(findall(x -> x <= mean_bkg_s_r, array_reprod[:, 1]))
	    idxs_large = Tuple(findall(x -> x <= mean_bkg_s, array[:, 1]))
	    global indxs_large_tot = Tuple(unique([idxs_r_large..., idxs_large...]))
	
		idxs_r_med = Tuple(findall(x -> x <= mean_bkg_s_r, array_reprod[:, 2]))
	    idxs_med = Tuple(findall(x -> x <= mean_bkg_s, array[:, 2]))
	    indxs_med_tot = Tuple(unique([idxs_r_med..., idxs_med...]))
	
		idxs_r_small = Tuple(findall(x -> x <= mean_bkg_s_r, array_reprod[:, 3]))
	    idxs_small = Tuple(findall(x -> x <= mean_bkg_s, array[:, 3]))
	    indxs_small_tot = Tuple(unique([idxs_r_small..., idxs_small...]))
	
		if length(indxs_large_tot) > 0
	        df_large = df[Not(indxs_large_tot...), :]
	        df_r_large = df_reprod[Not(indxs_large_tot...), :]
		else
			df_large = copy(df)
			df_r_large = copy(df_reprod)
	    end
	    if length(indxs_med_tot) > 0
	        df_med = df[Not(indxs_med_tot...), :]
	        df_r_med = df_reprod[Not(indxs_med_tot...), :]
	    else
			df_med = copy(df)
			df_r_med = copy(df_reprod)
	    end
	    if length(indxs_small_tot) > 0
	        df_small = df[Not(indxs_small_tot...), :]
	        df_r_small = df_reprod[Not(indxs_small_tot...), :]
		else
			df_small = copy(df)
			df_r_small= copy(df_reprod)
	    end
	
		return df_large, df_r_large, df_med, df_r_med, df_small, df_r_small
	else
		error("Not using SWCS)")
	end
end

# ╔═╡ 3a8c6f3c-0fa6-477c-85e2-7739cf2537ae
md"""
## Sensitivity
"""

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
	array_s = Array(df_s[!, 12:14])
    array_s_80 = Array(df_s_80[!, 12:14])
    array_s_100 = Array(df_s_100[!, 12:14])
    array_s_120 = Array(df_s_120[!, 12:14])
    array_s_135 = Array(df_s_135[!, 12:14])
end;

# ╔═╡ 97f07af4-a4af-4a26-9f81-5acb836c9f2c
df_s_large, df_s_r_large, df_s_medium, df_s_r_medium, df_s_small, df_s_r_small = remove_false_negatives(df_s, array_s, df_s_r, array_s_r, true);

# ╔═╡ a5dbb0ed-f3d9-4a7d-b830-6b6acf2d685b
r_squared_reprod_s, rms_values_reprod_s, fitted_line_reprod_s, coefficient_reprod_s = prepare_linear_regression(df_s_r_large, df_s_r_medium, df_s_r_small, df_s_large, df_s_medium, df_s_small)

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
array_a = hcat(df_a[!, :calculated_mass_large], df_a[!, :calculated_mass_medium], df_a[!, :calculated_mass_small]);

# ╔═╡ ef060aa5-860a-471d-a7d2-9a6d72e0c3fe
df_a_large, df_a_r_large, df_a_medium, df_a_r_medium, df_a_small, df_a_r_small = remove_false_negatives(df_a, array_a, df_a_r, array_a_r);

# ╔═╡ 3d987e71-9742-4b70-be10-2ed82a6df5f4
r_squared_reprod_a, rms_values_reprod_a, fitted_line_reprod_a, coefficient_reprod_a = prepare_linear_regression(df_a_r_large, df_a_r_medium, df_a_r_small, df_a_large, df_a_medium, df_a_small)

# ╔═╡ a31dfad2-e420-4c41-b805-227778a39fc9
num_zero_a = length(findall(x -> x <= 0, array_a))

# ╔═╡ 86d5f7fd-6a1b-4fdc-9e06-6cd66f8a169c
md"""
#### Integrated
"""

# ╔═╡ 79ad1bda-835e-4da2-8eae-6047a25ed99e
array_i = hcat(df_i[!, :calculated_mass_large], df_i[!, :calculated_mass_medium], df_i[!, :calculated_mass_small]);

# ╔═╡ 4b79a10e-d861-44ee-8b07-877e51145d5d
df_i_large, df_i_r_large, df_i_medium, df_i_r_medium, df_i_small, df_i_r_small = remove_false_negatives(df_i, array_i, df_i_r, array_i_r);

# ╔═╡ 14cdef5f-c9ba-4810-93fe-0f0915b803a2
r_squared_reprod_i, rms_values_reprod_i, fitted_line_reprod_i, coefficient_reprod_i = prepare_linear_regression(df_i_r_large, df_i_r_medium, df_i_r_small, df_i_large, df_i_medium, df_i_small);

# ╔═╡ fb50a01e-819f-4e4b-aef4-39d8255443a3
total_zero_i = length(findall(x -> x <= 0, array_i))

# ╔═╡ 99162b79-82a7-4293-b616-5e7669d04235
md"""
#### Volume Fraction
"""

# ╔═╡ 9abd5004-6d24-49e0-95a1-e62d1ce853cf
array_vf = hcat(df_vf[!, :calculated_mass_large], df_vf[!, :calculated_mass_medium], df_vf[!, :calculated_mass_small]);

# ╔═╡ e095f5ff-0687-4cb8-8075-5f98d1f4dd0e
df_vf_large, df_vf_r_large, df_vf_medium, df_vf_r_medium, df_vf_small, df_vf_r_small = remove_false_negatives(df_vf, array_vf, df_vf_r, array_vf_r);

# ╔═╡ 54e384c2-ac35-4679-8be5-ac96fa778c81
r_squared_reprod_vf, rms_values_reprod_vf, fitted_line_reprod_vf, coefficient_reprod_vf = prepare_linear_regression(df_vf_r_large, df_vf_r_medium, df_vf_r_small, df_vf_large, df_vf_medium, df_vf_small);

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
    lines!(ax1, collect(1:1000), fitted_line_reprod_i; linestyle=:dashdot)
    # if co1r[1] > 0
    #     Textbox(
    #         f[1, 1];
    #         placeholder="y = $(trunc(co1r[2]; digits=3))x+$(trunc(co1r[1]; digits=3)) \nr = $(trunc(r2_1r; digits=3)) \nRMSE: $(trunc(rms_values1r[1]; digits=3)) \nRMSD: $(trunc(rms_values1r[2]; digits=3))",
    #         tellheight=false,
    #         tellwidth=false,
    #         boxcolor=:white,
    #         halign=:left,
    #         valign=:top,
    #         textsize=12,
    #     )
    # else
    #     Textbox(
    #         f[1, 1];
    #         placeholder="y = $(trunc(co1r[2]; digits=3))x$(trunc(co1r[1]; digits=3)) \nr = $(trunc(r2_1r; digits=3)) \nRMSE: $(trunc(rms_values1r[1]; digits=3)) \nRMSD: $(trunc(rms_values1r[2]; digits=3))",
    #         tellheight=false,
    #         tellwidth=false,
    #         boxcolor=:white,
    #         halign=:left,
    #         valign=:top,
    #         textsize=12,
    #     )
    # end

	if coefficient_reprod_i[1] > 0
        Textbox(
            f[1, 1];
            placeholder="RMSE: $(trunc(rms_values_reprod_i[1]; digits=3)) \nRMSD: $(trunc(rms_values_reprod_i[2]; digits=3))",
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
            placeholder="RMSE: $(trunc(rms_values_reprod_i[1]; digits=3)) \nRMSD: $(trunc(rms_values_reprod_i[2]; digits=3))",
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
        df_vf_r_large[!, :calculated_mass_large],
        df_vf_large[!, :calculated_mass_large];
        label="Large Inserts",
    )
    scatter!(
        ax2,
        df_vf_r_medium[!, :calculated_mass_medium],
        df_vf_medium[!, :calculated_mass_medium];
        label="Medium Inserts",
    )
    scatter!(
        ax2,
        df_vf_r_small[!, :calculated_mass_small],
        df_vf_small[!, :calculated_mass_small];
        label="Small Inserts",
        color=:red,
    )
    lines!(ax2, [-1000, 1000], [-1000, 1000]; label="Unity")
    lines!(ax2, collect(1:1000), fitted_line_reprod_vf; linestyle=:dashdot)
    # if co2r[1] > 0
    #     Textbox(
    #         f[2, 1];
    #         placeholder="y = $(trunc(co2r[2]; digits=3))x+$(trunc(co2r[1]; digits=3)) \nr = $(trunc(r2_2r; digits=3)) \nRMSE: $(trunc(rms_values2r[1]; digits=3)) \nRMSD: $(trunc(rms_values2r[2]; digits=3))",
    #         tellheight=false,
    #         tellwidth=false,
    #         boxcolor=:white,
    #         halign=:left,
    #         valign=:top,
    #         textsize=12,
    #     )
    # else
    #     Textbox(
    #         f[2, 1];
    #         placeholder="y = $(trunc(co2r[2]; digits=3))x$(trunc(co2r[1]; digits=3)) \nr = $(trunc(r2_2r; digits=3)) \nRMSE: $(trunc(rms_values2r[1]; digits=3)) \nRMSD: $(trunc(rms_values2r[2]; digits=3))",
    #         tellheight=false,
    #         tellwidth=false,
    #         boxcolor=:white,
    #         halign=:left,
    #         valign=:top,
    #         textsize=12,
    #     )
    # end

	if coefficient_reprod_vf[1] > 0
        Textbox(
            f[2, 1];
            placeholder="RMSE: $(trunc(rms_values_reprod_vf[1]; digits=3)) \nRMSD: $(trunc(rms_values_reprod_vf[2]; digits=3))",
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
            placeholder="nRMSE: $(trunc(rms_values_reprod_vf[1]; digits=3)) \nRMSD: $(trunc(rms_values_reprod_vf[2]; digits=3))",
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
    ax2.title = "Volume Fraction"

    ##-- C --##
    ax3 = Axis(f[1, 2])
    scatter!(
        ax3,
        df_a_r_large[!, :calculated_mass_large],
        df_a_large[!, :calculated_mass_large];
        label="Large Inserts",
    )
    scatter!(
        ax3,
        df_a_r_medium[!, :calculated_mass_medium],
        df_a_medium[!, :calculated_mass_medium];
        label="Medium Inserts",
    )
    scatter!(
        ax3,
        df_a_r_small[!, :calculated_mass_small],
        df_a_small[!, :calculated_mass_small];
        label="Small Inserts",
        color=:red,
    )
    lines!(ax3, [-1000, 1000], [-1000, 1000]; label="Unity")
    lines!(ax3, collect(1:1000), fitted_line_reprod_a; linestyle=:dashdot)
    # if co2r[1] > 0
    #     Textbox(
    #         f[2, 1];
    #         placeholder="y = $(trunc(co2r[2]; digits=3))x+$(trunc(co2r[1]; digits=3)) \nr = $(trunc(r2_2r; digits=3)) \nRMSE: $(trunc(rms_values2r[1]; digits=3)) \nRMSD: $(trunc(rms_values2r[2]; digits=3))",
    #         tellheight=false,
    #         tellwidth=false,
    #         boxcolor=:white,
    #         halign=:left,
    #         valign=:top,
    #         textsize=12,
    #     )
    # else
    #     Textbox(
    #         f[2, 1];
    #         placeholder="y = $(trunc(co2r[2]; digits=3))x$(trunc(co2r[1]; digits=3)) \nr = $(trunc(r2_2r; digits=3)) \nRMSE: $(trunc(rms_values2r[1]; digits=3)) \nRMSD: $(trunc(rms_values2r[2]; digits=3))",
    #         tellheight=false,
    #         tellwidth=false,
    #         boxcolor=:white,
    #         halign=:left,
    #         valign=:top,
    #         textsize=12,
    #     )
    # end

	if coefficient_reprod_a[1] > 0
        Textbox(
            f[1, 2];
            placeholder="RMSE: $(trunc(rms_values_reprod_a[1]; digits=3)) \nRMSD: $(trunc(rms_values_reprod_a[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            textsize=12,
        )
    else
        Textbox(
            f[1, 2];
            placeholder="nRMSE: $(trunc(rms_values_reprod_a[1]; digits=3)) \nRMSD: $(trunc(rms_values_reprod_a[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            textsize=12,
        )
    end

    xlims!(ax3; low=0, high=200)
    ylims!(ax3; low=0, high=200)
    ax3.xticks = [0, 50, 100, 150, 200]
    ax3.yticks = [0, 50, 100, 150, 200]
    ax3.xlabel = "Mass 1 (mg)"
    ax3.ylabel = "Mass 2 (mg)"
    ax3.title = "Agatston"

    # ##-- C --##
    ax4 = Axis(f[2, 2])
    scatter!(
        ax4,
        df_s_r_large[!, :calculated_swcs_large],
        df_s_large[!, :calculated_swcs_large];
        label="Large Inserts",
    )
    scatter!(
        ax4,
        df_s_r_medium[!, :calculated_swcs_medium],
        df_s_medium[!, :calculated_swcs_medium];
        label="Medium Inserts",
    )
    scatter!(
        ax4,
        df_s_r_small[!, :calculated_swcs_small],
        df_s_small[!, :calculated_swcs_small];
        label="Small Inserts",
        color=:red,
    )
    lines!(ax4, [-1000, 1000], [-1000, 1000]; label="Unity")
    lines!(ax4, collect(1:1000), fitted_line_reprod_s; linestyle=:dashdot, label="Fitted Line")
    # if co3r[1] > 0
    #     Textbox(
    #         f[3, 1];
    #         placeholder="y = $(trunc(co3r[2]; digits=3))x+$(trunc(co3r[1]; digits=3)) \nr = $(trunc(r2_3r; digits=3)) \nRMSE: $(trunc(rms_values3r[1]; digits=3)) \nRMSD: $(trunc(rms_values3r[2]; digits=3))",
    #         tellheight=false,
    #         tellwidth=false,
    #         boxcolor=:white,
    #         halign=:left,
    #         valign=:top,
    #         textsize=12,
    #     )
    # else
    #     Textbox(
    #         f[3, 1];
    #         placeholder="y = $(trunc(co3r[2]; digits=3))x$(trunc(co3r[1]; digits=3)) \nr = $(trunc(r2_3r; digits=3)) \nRMSE: $(trunc(rms_values3r[1]; digits=3)) \nRMSD: $(trunc(rms_values3r[2]; digits=3))",
    #         tellheight=false,
    #         tellwidth=false,
    #         boxcolor=:white,
    #         halign=:left,
    #         valign=:top,
    #         textsize=12,
    #     )
    # end

	if coefficient_reprod_s[1] > 0
        Textbox(
            f[2, 2];
            placeholder="RMSE: $(trunc(rms_values_reprod_s[1]; digits=3)) \nRMSD: $(trunc(rms_values_reprod_s[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            textsize=12,
        )
    else
        Textbox(
            f[2, 2];
            placeholder="RMSE: $(trunc(rms_values_reprod_s[1]; digits=3)) \nRMSD: $(trunc(rms_values_reprod_s[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            textsize=12,
        )
    end

    xlims!(ax4; low=0, high=500)
    ylims!(ax4; low=0, high=500)
    ax4.xticks = [0, 125, 250, 375, 500]
    ax4.yticks = [0, 125, 250, 375, 500]
    ax4.xlabel = "SWCS 1"
    ax4.ylabel = "SWCS 2"
    ax4.title = "Spatially Weighted"

    ##-- LABELS --##
    f[1:2, 3] = Legend(f, ax3; framevisible=false)
    for (label, layout) in zip(["A", "B", "C", "D"], [f[1, 1], f[2, 1], f[1, 2], f[2, 2]])
        Label(
            layout[1, 1, TopLeft()],
            label;
            textsize=25,
            padding=(0, 0, 40, 0),
            halign=:right,
        )
    end

    save(joinpath(dirname(pwd()),"figures", FIGURE_PATH, "reproducibility.png"), f)
    return f
end

# ╔═╡ d2b90d91-4e63-45b0-8273-5231dbf2778e
with_theme(medphys_theme) do
    reprod()
end

# ╔═╡ 965264df-f9d8-4730-b7bb-d7f5334bb5c5
total_zero_vf = length(findall(x -> x <= 0, array_vf))

# ╔═╡ 01ecd9f2-4aa6-4f9e-b1ec-e51e3aa99571
function false_negative()
    f = Figure()
    colors = Makie.wong_colors()

    ##-- TOP --##
    axtop = Axis(f[1, 1]; xticks=(1:4, ["Integrated", "Volume Fraction", "Spatially Weighted", "Agatston"]))

    table = [1, 2, 3, 4]
    heights1 = [
		(total_zero_i / total_cac) * 100,
		(total_zero_vf / total_cac) * 100,
        (total_zero_s / total_cac) * 100,
        (num_zero_a / total_cac) * 100,
    ]
    barplot!(axtop, table, heights1; color=colors[1:4], bar_labels=:y)

    axtop.title = "False-Negative Scores (CAC=0)"
    axtop.ylabel = "% False-Negative Zero CAC Scores"
    ylims!(axtop; low=0, high=100)
    axtop.yticks = [0, 25, 50, 75, 100]

    save(joinpath(dirname(pwd()),"figures", FIGURE_PATH, "false_negative.png"), f)
    return f
end

# ╔═╡ cc88bab2-8b5e-474a-889d-f0fbfec790a2
with_theme(medphys_theme) do
    false_negative()
end

# ╔═╡ f4b51729-d53a-40e2-919b-5d4562cffbc0
total_zero_i, total_zero_vf, total_zero_s, num_zero_a

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

# ╔═╡ 235bc94e-166d-4df1-9ed6-f331e6e1cf0c
md"""
## Summaries
"""

# ╔═╡ 509ef00f-0a2f-4fc2-84b4-e20d43f3401f
md"""
### Accuracy
"""

# ╔═╡ 3d970695-2ec7-4223-9ebd-5060cade30c5
md"""
#### Normal
"""

# ╔═╡ 7c314dcf-aa06-491d-bb45-89c5c8006688
r_squared_values_normal = [
	r_squared_normal_i,
	r_squared_normal_vf,
	r_squared_normal_a,
]

# ╔═╡ 5daffb8c-9387-4bdf-9fd1-46dd8bd1b873
rmse_values_normal = [
	rms_values_normal_i[1],
	rms_values_normal_vf[1],
	rms_values_normal_a[1]
]

# ╔═╡ 9f1b2d37-1d54-4339-869a-a927c66618c4
rmsd_values_normal = [
	rms_values_normal_i[2],
	rms_values_normal_vf[2],
	rms_values_normal_a[2]
]

# ╔═╡ a6c8ebe4-e598-4cb7-a2c2-568f9ddd1b3a
summ_regression_normal = DataFrame(
	"Techniques" => ["Integrated", "Volume Fraction", "Agatston"],
	"R Correlation Coefficient" => r_squared_values_normal,
	"RMSE" => rmse_values_normal,
	"RMSD" => rmsd_values_normal
)

# ╔═╡ 39f1f521-990a-4646-94bd-93def74bdafb
md"""
#### Low
"""

# ╔═╡ 387c9e36-b67a-4aaa-8138-aef34951702e
r_squared_values_low = [
	r_squared_low_i,
	r_squared_low_vf,
	r_squared_low_a,
]

# ╔═╡ 7f5848a9-8d1b-4b81-b90f-d1c88a7f1fee
rmse_values_low = [
	rms_values_low_i[1],
	rms_values_low_vf[1],
	rms_values_low_a[1]
]

# ╔═╡ 142a8cb0-e296-47d1-b584-24a35b8d3339
rmsd_values_low = [
	rms_values_low_i[2],
	rms_values_low_vf[2],
	rms_values_low_a[2]
]

# ╔═╡ d3ebd2b6-c3c6-415a-b656-9025c2a66836
summ_regression_low = DataFrame(
	"Techniques" => ["Integrated", "Volume Fraction", "Agatston"],
	"R Correlation Coefficient" => r_squared_values_low,
	"RMSE" => rmse_values_low,
	"RMSD" => rmsd_values_low
)

# ╔═╡ 4ab40ed8-eaa8-42a8-a9ef-5cd06161a9c2
md"""
### Reproducibility
"""

# ╔═╡ e9dc2696-d2a8-4c4f-8665-da526e242570
r_squared_values_reprod = [
	r_squared_reprod_i
	r_squared_reprod_a
	r_squared_reprod_s
]

# ╔═╡ 574b38a3-964a-4d7a-adf4-e36198c003a3
rmse_values_reprod = [
	rms_values_reprod_i[1]
	rms_values_reprod_a[1]
	rms_values_reprod_s[1]
]

# ╔═╡ bc0ede74-4122-4d56-8d99-0043013752df
rmsd_values_reprod = [
	rms_values_reprod_i[2]
	rms_values_reprod_a[2]
	rms_values_reprod_s[2]
]

# ╔═╡ 165f2e6f-afd4-4f65-aae7-0023a531dbd0
summ_reprod = DataFrame(
	"Techniques" => ["Integrated", "SWCS", "Agatston"],
	"R Correlation Coefficient" => r_squared_values_reprod,
	"RMSE" => rmse_values_reprod,
	"RMSD" => rmsd_values_reprod
)

# ╔═╡ 70043dad-f0db-4ebc-a9a6-57ba93bb9ae2
md"""
### Sensitivity
"""

# ╔═╡ 8ef3a331-fd8f-4b61-b355-b133f0480432
zero_cac_num = [
	string(total_zero_i, "/", total_cac)
	string(total_zero_vf, "/", total_cac)
	string(total_zero_s, "/", total_cac)
	string(num_zero_a, "/", total_cac)
]

# ╔═╡ 6cf69eaf-aaf2-4bc1-aa10-e22efae2b175
zero_cac_perc = [
	round(total_zero_i / total_cac * 100, digits=2)
	round(total_zero_vf / total_cac * 100, digits=2)
	round(total_zero_s / total_cac * 100, digits=2)
	round(num_zero_a / total_cac * 100, digits=2)
]

# ╔═╡ 38f87d7b-6d22-41de-97da-927a478d8e8b
summ_zero_cac = DataFrame(
	"Techniques" => ["Integrated", "Volume Fraction", "Spatially Weighted", "Agatston"],
	"False Negatives (CAC=0)" => zero_cac_num,
	"Percentage False Negatives (%)" => zero_cac_perc
)

# ╔═╡ Cell order:
# ╠═c82bb3e6-e052-11ec-3392-eba4a4b464ba
# ╠═2df207ed-cf9c-4d7d-8354-4da14f93276c
# ╟─d72c94c2-a208-430c-9d76-85d8cfb33a22
# ╟─cabe7e9a-e932-406b-95dd-2c9128decdc7
# ╟─4a255a58-b9a0-4175-a750-b2562361631d
# ╟─7fc65de7-7514-4532-803e-1120150fd2fb
# ╟─1ecb8fd3-a36b-4b06-8daf-e598e297ecf8
# ╠═043dec5f-07a8-472d-8082-c7a771f93270
# ╟─e79c34d0-bf06-4a7a-93b7-e49f09ce3f4e
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
# ╟─090bddd5-081a-43c8-b8d1-23fce50331e8
# ╠═993c899f-d741-40d5-81e0-ea7fd2871146
# ╠═5e803aad-9d82-4801-aff6-3c1c2f748b02
# ╠═5c76c3c7-bb56-4d90-8f0f-7f7d0e271d3c
# ╠═0085693a-157d-42a8-9d35-ec0feec24c9c
# ╠═a03343e0-a9f4-4a49-8490-1d542bbe42fd
# ╟─310c1999-77a6-432e-850c-4111411cceb0
# ╠═3c40ca1d-4111-4bf1-941e-30b59cd1a264
# ╟─7dfc24a4-e006-45f4-b5b9-977a7c3c0b7c
# ╟─c93ba92a-ccc7-4ae0-9207-300715821fc5
# ╟─a4bc0a8a-4904-4217-97a9-44158f99ae70
# ╟─2e04217d-8dfb-48c9-85dc-9fb42cfd3039
# ╟─98c3eaef-502e-4a20-b5b6-46a3d6b394d3
# ╟─79a09453-d375-4f1d-9318-887d806dbb14
# ╠═166718a5-1403-4f1b-b0f7-539fa006997a
# ╠═bcedd090-80d5-480a-aad7-97a5fb8131f5
# ╠═e0175e63-0ba9-4763-8ce6-e1df05966497
# ╟─cf59cd3f-026c-4321-a1ac-de217177b52e
# ╠═5d7ea4ec-088c-4f94-bb64-531673836115
# ╠═1886135c-4ecc-4f27-9a77-7ead70d7735d
# ╠═8b5c75a7-8973-46bf-88bf-0fe65f4d20ce
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
# ╠═4b79a10e-d861-44ee-8b07-877e51145d5d
# ╠═14cdef5f-c9ba-4810-93fe-0f0915b803a2
# ╟─4c9feda4-9fc1-42c7-8cac-a17d51627c6d
# ╠═e261c635-98d7-4357-a429-a34b73b47bb5
# ╠═d2e4fa84-1f12-4924-823a-6ea7ab0d3d8e
# ╠═101b02a5-9e67-460d-9db0-140684bd1e11
# ╠═fe275fdc-c833-43ac-a4c7-7a49b2b728f8
# ╠═8babe305-6bf0-43e3-9fe7-ded195d2addf
# ╠═ea17e8f9-c06e-439e-82f9-dbfdef759998
# ╠═e095f5ff-0687-4cb8-8075-5f98d1f4dd0e
# ╠═54e384c2-ac35-4679-8be5-ac96fa778c81
# ╟─8eaec985-641a-498f-9cd9-c2c85d877142
# ╠═fcd83805-06e2-4cda-aa56-0c13c69424d8
# ╠═5f2eb2d2-84c4-4160-a92a-f181b4126450
# ╠═869de9cc-6ab0-4e0b-ad1f-59db1fd2f69f
# ╠═43b3ba9b-247c-41cd-968e-3a27736426b0
# ╠═b16090b5-a232-4e71-9148-b71296efa999
# ╠═24034fa9-f9d4-4bd0-a736-488b06403adc
# ╠═ef060aa5-860a-471d-a7d2-9a6d72e0c3fe
# ╠═3d987e71-9742-4b70-be10-2ed82a6df5f4
# ╟─5bc912d4-be86-4cca-9e8b-8af31d269edf
# ╠═594b8d0e-f65c-4fe3-b571-e2ee68a9ab5f
# ╠═1d208fd8-b847-4afe-84b5-3a553f50a858
# ╠═a9f68d7e-7818-47d5-ba1f-965634225b30
# ╠═f2177e1f-2a55-4f81-9a69-ab985ae3b7c2
# ╠═35bc90e2-9d70-4daf-a825-b462831c5bf6
# ╠═624074ec-572a-45ab-932a-aed7edb1f846
# ╠═e8b56e73-f22b-4587-9141-558f5f545a5e
# ╠═97f07af4-a4af-4a26-9f81-5acb836c9f2c
# ╠═a5dbb0ed-f3d9-4a7d-b830-6b6acf2d685b
# ╟─3a8c6f3c-0fa6-477c-85e2-7739cf2537ae
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
# ╠═79ad1bda-835e-4da2-8eae-6047a25ed99e
# ╠═fb50a01e-819f-4e4b-aef4-39d8255443a3
# ╟─99162b79-82a7-4293-b616-5e7669d04235
# ╠═9abd5004-6d24-49e0-95a1-e62d1ce853cf
# ╠═965264df-f9d8-4730-b7bb-d7f5334bb5c5
# ╟─a6f0d36e-823d-4198-b4ea-d95ce16ac65c
# ╟─b7064bb1-f21a-4bb6-a395-9f1c73102058
# ╟─21b03bf8-47ba-4915-936a-92c488671bb1
# ╟─e0a66c8a-e1b6-4c70-9996-ec3fe62786f3
# ╠═bf3b0e00-5f2b-435b-9717-96052f7e9071
# ╟─09e62a9b-4f42-4285-88de-3b90a8a4003a
# ╟─d74aa830-5a48-4535-8579-39ce2dc8142d
# ╟─6eed3e33-f8f1-41f8-96b1-6e38fa129733
# ╟─36962426-647f-4194-92cc-aed51c628484
# ╠═512f3683-3231-4140-96d8-dd2378f5094d
# ╟─235bc94e-166d-4df1-9ed6-f331e6e1cf0c
# ╟─509ef00f-0a2f-4fc2-84b4-e20d43f3401f
# ╟─3d970695-2ec7-4223-9ebd-5060cade30c5
# ╟─7c314dcf-aa06-491d-bb45-89c5c8006688
# ╟─5daffb8c-9387-4bdf-9fd1-46dd8bd1b873
# ╟─9f1b2d37-1d54-4339-869a-a927c66618c4
# ╠═a6c8ebe4-e598-4cb7-a2c2-568f9ddd1b3a
# ╟─39f1f521-990a-4646-94bd-93def74bdafb
# ╟─387c9e36-b67a-4aaa-8138-aef34951702e
# ╟─7f5848a9-8d1b-4b81-b90f-d1c88a7f1fee
# ╟─142a8cb0-e296-47d1-b584-24a35b8d3339
# ╠═d3ebd2b6-c3c6-415a-b656-9025c2a66836
# ╟─4ab40ed8-eaa8-42a8-a9ef-5cd06161a9c2
# ╟─e9dc2696-d2a8-4c4f-8665-da526e242570
# ╟─574b38a3-964a-4d7a-adf4-e36198c003a3
# ╟─bc0ede74-4122-4d56-8d99-0043013752df
# ╠═165f2e6f-afd4-4f65-aae7-0023a531dbd0
# ╟─70043dad-f0db-4ebc-a9a6-57ba93bb9ae2
# ╟─8ef3a331-fd8f-4b61-b355-b133f0480432
# ╟─6cf69eaf-aaf2-4bc1-aa10-e22efae2b175
# ╠═38f87d7b-6d22-41de-97da-927a478d8e8b

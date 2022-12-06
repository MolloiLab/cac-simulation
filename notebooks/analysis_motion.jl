### A Pluto.jl notebook ###
# v0.19.16

using Markdown
using InteractiveUtils

# ╔═╡ c7043dc2-8d49-45db-8f89-2490094b60a6
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.activate(".")

    using PlutoUI, Statistics, CSV, DataFrames, GLM, CairoMakie, HypothesisTests, Colors, MLJBase
	using StatsBase: quantile!, rmsd
end

# ╔═╡ fda1f661-91bb-4abe-83ba-200095e7efc6
TableOfContents()

# ╔═╡ 6658974b-e177-45dd-ba28-03e54440296e
md"""
#### Helper Functions
"""

# ╔═╡ 44804164-d4fd-42a6-8018-3032f2929a1a
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

# ╔═╡ 7f10366d-cc73-4ca7-9ebe-14c360cb8d34
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

# ╔═╡ cd1d88db-fc7b-4790-b57e-d5c6d34a9657
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

# ╔═╡ 62e3c411-70e5-4453-9226-2d11cb3c91fd
md"""
# Load CSVs
"""

# ╔═╡ f4927d17-9916-432d-9880-47ac4f83e859
md"""
## Integrated
"""

# ╔═╡ fcd9c21d-2f6c-4421-bad1-bb5c239329fb
path_integrated = "/Users/daleblack/Google Drive/dev/MolloiLab/cac-simulation/output_new/integrated_scoring";

# ╔═╡ f8f7f2bf-7afd-474a-ba0d-9c45f16c5ddd
df_i = CSV.read(string(path_integrated, "/full2.csv"), DataFrame)

# ╔═╡ fe082455-1797-4d28-bd05-621fc9dcc540
df_i_low, df_i_normal = groupby(df_i, :DENSITY);

# ╔═╡ 460410a9-2dc4-42ef-b1a8-997e04727e0f
df_i_low_small, df_i_low_medium, df_i_low_large = groupby(df_i_low, :SIZE);

# ╔═╡ c7abbd34-d2df-403e-96a4-764ebcb57702
df_i_normal_small, df_i_normal_medium, df_i_normal_large = groupby(df_i_normal, :SIZE);

# ╔═╡ d9fb7c5c-16df-4c97-92c7-baabe1905c0f
md"""
## Agatston
"""

# ╔═╡ cbde2904-fbcc-4863-9bce-fe448fdc5ea7
path_agat = "/Users/daleblack/Google Drive/dev/MolloiLab/cac-simulation/output_new/agatston";

# ╔═╡ b3bf2fb7-fa89-4efb-b6e2-65c5baff8b4d
df_a = CSV.read(string(path_agat, "/full2.csv"), DataFrame)

# ╔═╡ 3bc73f74-e009-4f2d-ba1b-32539dfa8505
df_a_low, df_a_normal = groupby(df_a, :DENSITY);

# ╔═╡ 9cf1580b-8fd6-4881-88fd-25fd60f48807
df_a_low_small, df_a_low_medium, df_a_low_large = groupby(df_a_low, :SIZE);

# ╔═╡ f85b9c18-8c65-47d7-9949-e424e6cb2e29
df_a_normal_small, df_a_normal_medium, df_a_normal_large = groupby(df_a_normal, :SIZE);

# ╔═╡ cf8b2459-1204-4035-8f8c-79e92a4bab3e
md"""
## Spatially Weighted
"""

# ╔═╡ 50413fca-f3ec-415d-a83e-d118beaf5009
path_swcs = "/Users/daleblack/Google Drive/dev/MolloiLab/cac-simulation/output_new/swcs";

# ╔═╡ 0fdf77f9-7906-49d3-8556-0903ba76301f
df_s = CSV.read(string(path_swcs, "/full2.csv"), DataFrame)

# ╔═╡ f9267f12-3480-4b9d-ab22-68daf4b8077e
df_s_low, df_s_normal = groupby(df_s, :DENSITY);

# ╔═╡ 5aae2475-846f-49ea-bbe8-87f3af101999
df_s_low_small, df_s_low_medium, df_s_low_large = groupby(df_s_low, :SIZE);

# ╔═╡ 9a21414e-d55f-4094-b28e-bcd72b14c3e3
df_s_normal_small, df_s_normal_medium, df_s_normal_large = groupby(df_s_normal, :SIZE);

# ╔═╡ ede89d4a-54b8-46c8-b6a4-aa9c4c7137df
md"""
## Volume Fraction
"""

# ╔═╡ b95045b3-e108-4464-8fb3-1bf3ca4eea45
path_vf = "/Users/daleblack/Google Drive/dev/MolloiLab/cac-simulation/output_new/volume_fraction";

# ╔═╡ fa637952-26e2-4bf3-8184-79200cc0621f
df_vf = CSV.read(string(path_vf, "/full2.csv"), DataFrame)

# ╔═╡ e7c5bc2b-d6fd-493b-8aec-32167c9701ca
df_vf_low, df_vf_normal = groupby(df_vf, :DENSITY);

# ╔═╡ a8d566eb-e4a8-4ab5-959e-9ae9915267c1
df_vf_low_small, df_vf_low_medium, df_vf_low_large = groupby(df_vf_low, :SIZE);

# ╔═╡ b28a4846-7ca3-4827-b1be-c0f6951050fc
df_vf_normal_small, df_vf_normal_medium, df_vf_normal_large = groupby(df_vf_normal, :SIZE);

# ╔═╡ 9db9fad7-8982-4f8c-8f4f-3196fc3f1af0
md"""
# Figures
"""

# ╔═╡ f7bb8133-456f-49ad-9690-a9b42e39144f
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

# ╔═╡ 7e55db13-2b25-4283-bc42-411af0af0a00
md"""
## Accuracy
"""

# ╔═╡ f2777887-fe96-4795-b0d8-8430749846fe
md"""
#### Normal Density
"""

# ╔═╡ 7d551a01-ecdf-473b-84d5-db310e65d4a1
r_squared_normal_i, rms_values_normal_i, fitted_line_normal_i, coefficient_normal_i = prepare_linear_regression(df_i_normal)

# ╔═╡ 5986e921-5e2a-4084-92a2-96aa2297edfd
r_squared_normal_vf, rms_values_normal_vf, fitted_line_normal_vf, coefficient_normal_vf = prepare_linear_regression(df_vf_normal)

# ╔═╡ 1e551f4e-8381-4ae1-980a-a660db4dd333
r_squared_normal_a, rms_values_normal_a, fitted_line_normal_a, coefficient_normal_a = prepare_linear_regression(df_a_normal)

# ╔═╡ e5366f83-eefc-41e1-9cda-ec9303808c20
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

    # save(
    #     "/Users/daleblack/Google Drive/Research/Papers/My Papers/cac-simulation/figures/linear_reg_norm.png",
    #     f,
    # )
    return f
end

# ╔═╡ 7af4578d-f999-4ee1-9d83-994539399b92
with_theme(medphys_theme) do
    lin_reg_norm()
end

# ╔═╡ 9e1ac65a-927a-45f9-a061-1f8ce137f1bb
md"""
#### Low Density
"""

# ╔═╡ 704bf7f2-e2ac-4ccf-9efe-dd5bc939263e
r_squared_low_i, rms_values_low_i, fitted_line_low_i, coefficient_low_i = prepare_linear_regression(df_i_low)

# ╔═╡ 9d12a490-4225-45cf-8315-5275d01fbc7d
r_squared_low_vf, rms_values_low_vf, fitted_line_low_vf, coefficient_low_vf = prepare_linear_regression(df_vf_low)

# ╔═╡ 55a71f6d-e3af-4fce-bb49-f9a68831d8f0
r_squared_low_a, rms_values_low_a, fitted_line_low_a, coefficient_low_a = prepare_linear_regression(df_a_low)

# ╔═╡ 0f086c11-fb86-48c3-a6bd-844207054dc9
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

    save(
        "/Users/daleblack/Google Drive/Research/Papers/My Papers/cac-simulation/figures/linear_reg_low.png",
        f,
    )
    return f
end

# ╔═╡ 0220f0a6-426a-4021-82c0-218793708ede
with_theme(medphys_theme) do
    lin_reg_low()
end

# ╔═╡ 0bd48352-a9a6-4b9c-a85d-31d2f99e96c6
md"""
## Reproducibility
"""

# ╔═╡ 9d631928-295f-449d-91f5-7b52e2b059e2
md"""
#### Integrated
"""

# ╔═╡ 5f5017c0-26f6-408a-98a6-f17fd436bd4e
path_integrated_r = "/Users/daleblack/Google Drive/dev/MolloiLab/cac-simulation/output_repeated/integrated_scoring";

# ╔═╡ bb042948-668f-45b1-a954-4b694e64ceaf
df_i_r = CSV.read(string(path_integrated_r, "/full2.csv"), DataFrame)

# ╔═╡ 168756a1-1965-40a4-ae5a-c11af0624c0a
df_i_low_r, df_i_normal_r = groupby(df_i_r, :DENSITY);

# ╔═╡ 67c093c7-4352-48d1-9fd1-634cefbea552
df_i_low_small_r, df_i_low_medium_r, df_i_low_large_r = groupby(df_i_low_r, :SIZE);

# ╔═╡ 6ea59087-3a3f-41c9-8b74-5f892d83b46d
df_i_normal_small_r, df_i_normal_medium_r, df_i_normal_large_r = groupby(
    df_i_normal_r, :SIZE
);

# ╔═╡ fd331a11-7ec5-45d1-aad1-521824c886d1
array_i_r = hcat(
    df_i_r[!, :calculated_mass_large],
    df_i_r[!, :calculated_mass_medium],
    df_i_r[!, :calculated_mass_small],
)

# ╔═╡ a42a3104-1858-4ab0-9de0-79c0d85bdaa8
md"""
#### Volume Fraction
"""

# ╔═╡ 471b8219-1268-4efd-8a94-0292336938f4
path_volume_fraction_r = "/Users/daleblack/Google Drive/dev/MolloiLab/cac-simulation/output_repeated/volume_fraction";

# ╔═╡ f95a03f3-be7d-4ff0-af37-165a5405c552
df_vf_r = CSV.read(string(path_volume_fraction_r, "/full2.csv"), DataFrame)

# ╔═╡ 5ce9729f-82fc-43a9-b809-699702255a15
df_vf_low_r, df_vf_normal_r = groupby(df_vf_r, :DENSITY);

# ╔═╡ 54f85688-77d1-4e91-acec-55c4c03c721a
df_vf_low_small_r, df_vf_low_medium_r, df_vf_low_large_r = groupby(df_vf_low_r, :SIZE);

# ╔═╡ 28488ced-f472-4e36-956e-988ef3b06361
df_vf_normal_small_r, df_vf_normal_medium_r, df_vf_normal_large_r = groupby(
    df_vf_normal_r, :SIZE
);

# ╔═╡ d48c2606-4b79-4150-9669-2a6efc036280
array_vf_r = hcat(
    df_vf_r[!, :calculated_mass_large],
    df_vf_r[!, :calculated_mass_medium],
    df_vf_r[!, :calculated_mass_small],
)

# ╔═╡ c86b93dc-18da-47b7-8099-b27396bdd48a
md"""
#### Agatston
"""

# ╔═╡ be39ddd8-16ee-4049-b61c-e16f85e6f289
path_agat_r = "/Users/daleblack/Google Drive/dev/MolloiLab/cac-simulation/output_repeated/agatston";

# ╔═╡ 8b208615-24dc-4c53-8cfa-f3ffd419dea3
df_a_r = CSV.read(string(path_agat_r, "/full2.csv"), DataFrame)

# ╔═╡ 316ce953-f138-4d32-a7bb-385ef932c8b2
df_a_low_r, df_a_normal_r = groupby(df_a_r, :DENSITY);

# ╔═╡ 03c97478-5f88-42fb-bf29-f01d1effd3eb
df_a_low_small_r, df_a_low_medium_r, df_a_low_large_r = groupby(df_a_low_r, :SIZE);

# ╔═╡ bea7c372-8aa8-4feb-ba6b-0bb3e624aa2c
df_a_normal_small_r, df_a_normal_medium_r, df_a_normal_large_r = groupby(
    df_a_normal_r, :SIZE
);

# ╔═╡ d49df961-f841-435c-9c2a-c4cf293cff32
array_a_r = hcat(
    df_a_r[!, :calculated_mass_large],
    df_a_r[!, :calculated_mass_medium],
    df_a_r[!, :calculated_mass_large],
);

# ╔═╡ d5a746f9-1351-40fd-b0c6-c7f5c9b6fdd4
md"""
#### SWCS
"""

# ╔═╡ a6ab095b-f34b-4e54-b752-d606ccd20b4a
path_swcs_r = "/Users/daleblack/Google Drive/dev/MolloiLab/cac-simulation/output_repeated/swcs";

# ╔═╡ cd62f88b-2054-4eb6-8c76-042af27d3b68
df_s_r = CSV.read(string(path_swcs_r, "/full2.csv"), DataFrame)

# ╔═╡ e9f69c55-3a91-4aa6-9113-143ccc495f4d
df_s_low_r, df_s_normal_r = groupby(df_s_r, :DENSITY);

# ╔═╡ 0bd80638-ed4c-4e13-bdcb-77e880453429
df_s_low_small_r, df_s_low_medium_r, df_s_low_large_r = groupby(df_s_low_r, :SIZE);

# ╔═╡ 26ce1629-f38c-4dd9-8ea4-5d49dff5eb7c
df_s_normal_small_r, df_s_normal_medium_r, df_s_normal_large_r = groupby(
    df_s_normal_r, :SIZE
);

# ╔═╡ 45dc885c-f00d-4873-939b-2f5f70c3fa18
array_s_r = hcat(
    df_s_r[!, :calculated_swcs_large],
    df_s_r[!, :calculated_swcs_medium],
    df_s_r[!, :calculated_swcs_small],
);

# ╔═╡ 3544bf0e-d154-4ccc-8696-4d9a7c954405
begin
    mean_bkg_s = mean(df_s[!, :swcs_bkg])
    mean_bkg_s_r = mean(df_s_r[!, :swcs_bkg])
end

# ╔═╡ e201ebfc-17f2-4be8-83a6-06578ae0ee77
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

# ╔═╡ 3114afe2-6d7f-466b-bbac-b2f8b9ed79a5
md"""
## Sensitivity
"""

# ╔═╡ 0c282db8-3e45-49bc-9183-f27412fdd273
md"""
#### SWCS
"""

# ╔═╡ b0ec93ef-5833-4d20-a8fe-959e97588816
df_s_80, df_s_100, df_s_120, df_s_135, = groupby(df_s, :scan);

# ╔═╡ 3c27fc44-e4bb-440b-b88b-3da8284dc329
begin
    mean_swcs_80, std_swcs_80 = mean(df_s_80[!, :swcs_bkg]), std(df_s_80[!, :swcs_bkg])
    mean_swcs_100, std_swcs_100 = mean(df_s_100[!, :swcs_bkg]), std(df_s_100[!, :swcs_bkg])
    mean_swcs_120, std_swcs_120 = mean(df_s_120[!, :swcs_bkg]), std(df_s_120[!, :swcs_bkg])
    mean_swcs_135, std_swcs_135 = mean(df_s_135[!, :swcs_bkg]), std(df_s_135[!, :swcs_bkg])
end;

# ╔═╡ bece0a05-e843-4be6-94d7-dd46c127cc92
mean_swcs_80, mean_swcs_100, mean_swcs_120, mean_swcs_135

# ╔═╡ 4f934e78-117e-4054-a101-ea04731b769f
begin
	array_s = Array(df_s[!, end-3:end-1])
    array_s_80 = Array(df_s_80[!, end-3:end-1])
    array_s_100 = Array(df_s_100[!, end-3:end-1])
    array_s_120 = Array(df_s_120[!, end-3:end-1])
    array_s_135 = Array(df_s_135[!, end-3:end-1])
end;

# ╔═╡ 03c62590-3615-4c06-b9f0-77267308c1a5
df_s_large, df_s_r_large, df_s_medium, df_s_r_medium, df_s_small, df_s_r_small = remove_false_negatives(df_s, array_s, df_s_r, array_s_r, true);

# ╔═╡ 99426eaa-b5a0-48aa-9731-843679fc74ab
df_s_r_large

# ╔═╡ 26904b95-eec0-4c10-b846-f18b66322a24
r_squared_reprod_s, rms_values_reprod_s, fitted_line_reprod_s, coefficient_reprod_s = prepare_linear_regression(df_s_r_large, df_s_r_medium, df_s_r_small, df_s_large, df_s_medium, df_s_small)

# ╔═╡ 26be4298-b7c3-4426-9ee9-aaaaa89b8d87
begin
    num_zeroCAC_80 = length(findall(x -> x < mean_swcs_80, array_s_80))
    num_zeroCAC_100 = length(findall(x -> x < mean_swcs_100, array_s_100))
    num_zeroCAC_120 = length(findall(x -> x < mean_swcs_120, array_s_120))
    num_zeroCAC_135 = length(findall(x -> x < mean_swcs_135, array_s_135))

    total_zero_s = num_zeroCAC_80 + num_zeroCAC_100 + num_zeroCAC_120 + num_zeroCAC_135
    total_cac = length(array_s_80) * 4
end;

# ╔═╡ e78bbfc6-8699-4488-bd3a-488486b3c921
total_cac

# ╔═╡ f2a8bd87-67c2-4e7e-8bce-6bf7e65cbe48
total_zero_s

# ╔═╡ 58447cc6-384a-4dc5-b733-1f8a6bfdbea1
md"""
#### Agatston
"""

# ╔═╡ 6328e32a-5fa3-4de8-b39e-58f800730fd5
array_a = hcat(df_a[!, :calculated_mass_large], df_a[!, :calculated_mass_medium], df_a[!, :calculated_mass_small]);

# ╔═╡ 4a88479f-f53a-4923-9a1b-141878a2b331
df_a_large, df_a_r_large, df_a_medium, df_a_r_medium, df_a_small, df_a_r_small = remove_false_negatives(df_a, array_a, df_a_r, array_a_r);

# ╔═╡ ca125df2-4b7a-4056-8f61-bba9b3dc79eb
r_squared_reprod_a, rms_values_reprod_a, fitted_line_reprod_a, coefficient_reprod_a = prepare_linear_regression(df_a_r_large, df_a_r_medium, df_a_r_small, df_a_large, df_a_medium, df_a_small)

# ╔═╡ 2e170bad-e903-483f-9cc9-5f7b2b8b2b17
num_zero_a = length(findall(x -> x <= 0, array_a))

# ╔═╡ 19aa4b7b-08aa-48fb-8c0c-5d85fdfd9f1d
md"""
#### Integrated
"""

# ╔═╡ 15f718ae-39d5-4aea-bf35-513d67d66710
array_i = hcat(df_i[!, :calculated_mass_large], df_i[!, :calculated_mass_medium], df_i[!, :calculated_mass_small]);

# ╔═╡ 9ee37984-894a-43cf-bc2e-e2a04db22bb2
df_i_large, df_i_r_large, df_i_medium, df_i_r_medium, df_i_small, df_i_r_small = remove_false_negatives(df_i, array_i, df_i_r, array_i_r);

# ╔═╡ ed01e6bb-3a7d-4bb7-8995-920a22fc5f48
r_squared_reprod_i, rms_values_reprod_i, fitted_line_reprod_i, coefficient_reprod_i = prepare_linear_regression(df_i_r_large, df_i_r_medium, df_i_r_small, df_i_large, df_i_medium, df_i_small);

# ╔═╡ ad844c27-cdd1-49d3-b71b-08bc38b628cf
total_zero_i = length(findall(x -> x <= 0, array_i))

# ╔═╡ 9af1bd1a-1ce4-4795-bf64-232fdabd3727
md"""
#### Volume Fraction
"""

# ╔═╡ d0a14e60-7f95-4fad-b132-8a3e90b1bd2a
array_vf = hcat(df_vf[!, :calculated_mass_large], df_vf[!, :calculated_mass_medium], df_vf[!, :calculated_mass_small]);

# ╔═╡ 42b5e2a8-71fe-4813-be24-255b9ed64330
df_vf_large, df_vf_r_large, df_vf_medium, df_vf_r_medium, df_vf_small, df_vf_r_small = remove_false_negatives(df_vf, array_vf, df_vf_r, array_vf_r);

# ╔═╡ 2b839a05-7007-4f9f-ac3f-6bbf51c5071c
r_squared_reprod_vf, rms_values_reprod_vf, fitted_line_reprod_vf, coefficient_reprod_vf = prepare_linear_regression(df_vf_r_large, df_vf_r_medium, df_vf_r_small, df_vf_large, df_vf_medium, df_vf_small);

# ╔═╡ 047a921c-2c65-4a1a-8f9c-a7d1ade380ce
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
            placeholder="RMSE: $(trunc(rms_values_reprod_vf[1]; digits=3)) \nRMSD: $(trunc(rms_values_reprod_vf[2]; digits=3))",
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
            placeholder="RMSE: $(trunc(rms_values_reprod_a[1]; digits=3)) \nRMSD: $(trunc(rms_values_reprod_a[2]; digits=3))",
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

    save(
        "/Users/daleblack/Google Drive/Research/Papers/My Papers/cac-simulation/figures-review/reprod.png",
        f,
    )
    return f
end

# ╔═╡ bc23bb1f-cb6d-4a80-82bd-9c587d29d39c
with_theme(medphys_theme) do
    reprod()
end

# ╔═╡ 8389c664-7da1-4df1-abd2-3f9dc4ab2775
total_zero_vf = length(findall(x -> x <= 0, array_vf))

# ╔═╡ 668b632c-bb5c-4b79-bd4f-69f6642a0bc8
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

    save(
        "/Users/daleblack/Google Drive/Research/Papers/My Papers/cac-simulation/figures-review/zero_cac.png",
        f,
    )
    return f
end

# ╔═╡ eb5dd9cf-70fb-415e-930a-e95393e15985
with_theme(medphys_theme) do
    false_negative()
end

# ╔═╡ fb91c84e-06f9-49aa-bf53-d6f8df85cdb4
total_zero_i, total_zero_vf, total_zero_s, num_zero_a

# ╔═╡ 47ae09c1-ddf3-487b-9463-626c15938320
md"""
# Tables
"""

# ╔═╡ d2829dd3-8e08-48ea-8636-f78bc203541d
md"""
## Simulation Parameters
"""

# ╔═╡ 0823bd5e-c106-4ef9-839c-8452a4ee8568
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

# ╔═╡ de822771-48f9-4620-a998-3d4235a5501d
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

# ╔═╡ 0c95bec7-f6cb-4c2b-91a2-275e05d9e572
simulation_parameters = DataFrame(
	"Parameter" => parameters,
	"Simulation" => simulations
)

# ╔═╡ 4eabcaae-55ea-4c3c-b276-15ee6f7ff9ec
md"""
## Spatially Weighted Calcium Scoring Means & Stds
"""

# ╔═╡ b05be082-6174-469b-be11-0b0c7d95f1f5
kvp = [
	80
	100
	120
	135
]

# ╔═╡ 494435e9-e008-4284-9132-bdc78cf999e3
means = [
	201.4898
	174.3658
	158.2645
	152.8815
]

# ╔═╡ d0b3d0fd-4db6-4d09-8005-733191b70076
stds = [
	34.3038
	28.1708
	22.9656
	21.0778
]

# ╔═╡ 97798f97-3912-48a8-8ffe-c05cdb736e09
means_stds = DataFrame(
	"kVp" => kvp,
	"Mean" => means,
	"Std" => stds
)

# ╔═╡ dbcfe49b-03c6-41f9-8003-9d054e782127
md"""
## Summaries
"""

# ╔═╡ 88b83264-0e98-43cc-be25-b0ad0311aebc
md"""
### Accuracy
"""

# ╔═╡ 1a12b990-db62-475f-819d-c76c9deadba7
md"""
#### Normal
"""

# ╔═╡ 7cabd9f5-e450-46e4-ad32-97e9c2b143d0
r_squared_values_normal = [
	r_squared_normal_i,
	r_squared_normal_vf,
	r_squared_normal_a,
]

# ╔═╡ e7461e35-bfe2-4b71-b785-02532e5ecf48
rmse_values_normal = [
	rms_values_normal_i[1],
	rms_values_normal_vf[1],
	rms_values_normal_a[1]
]

# ╔═╡ 596f29f3-fe65-4349-95e1-e2e42dd22c91
rmsd_values_normal = [
	rms_values_normal_i[2],
	rms_values_normal_vf[2],
	rms_values_normal_a[2]
]

# ╔═╡ ab6f7193-d8f2-4288-9c65-41e392552e84
summ_regression_normal = DataFrame(
	"Techniques" => ["Integrated", "Volume Fraction", "Agatston"],
	"R Correlation Coefficient" => r_squared_values_normal,
	"RMSE" => rmse_values_normal,
	"RMSD" => rmsd_values_normal
)

# ╔═╡ 90593bbe-4143-450e-b780-1bdd5ee88512
md"""
#### Low
"""

# ╔═╡ 244f3279-7a1e-45e3-b046-665331f28274
r_squared_values_low = [
	r_squared_low_i,
	r_squared_low_vf,
	r_squared_low_a,
]

# ╔═╡ c314868f-5f42-4211-9e35-20895f55c494
rmse_values_low = [
	rms_values_low_i[1],
	rms_values_low_vf[1],
	rms_values_low_a[1]
]

# ╔═╡ e14ed023-83de-43ac-9db2-38985ccfb4cc
rmsd_values_low = [
	rms_values_low_i[2],
	rms_values_low_vf[2],
	rms_values_low_a[2]
]

# ╔═╡ d245eb36-d4ba-4f43-a119-a49defdc85d8
summ_regression_low = DataFrame(
	"Techniques" => ["Integrated", "Volume Fraction", "Agatston"],
	"R Correlation Coefficient" => r_squared_values_low,
	"RMSE" => rmse_values_low,
	"RMSD" => rmsd_values_low
)

# ╔═╡ a17cd0a7-5164-446c-bfc6-4d097bd77063
md"""
### Reproducibility
"""

# ╔═╡ e83e986c-3840-4ece-a578-c78ae06a1df9
r_squared_values_reprod = [
	r_squared_reprod_i
	r_squared_reprod_a
	r_squared_reprod_s
]

# ╔═╡ 438d943e-0f45-4bd5-843d-58d92aa1e33e
rmse_values_reprod = [
	rms_values_reprod_i[1]
	rms_values_reprod_a[1]
	rms_values_reprod_s[1]
]

# ╔═╡ 0c73850e-c673-43b4-b623-288ef3084de3
rmsd_values_reprod = [
	rms_values_reprod_i[2]
	rms_values_reprod_a[2]
	rms_values_reprod_s[2]
]

# ╔═╡ eac67a26-c28d-4607-a4f2-e586d75837a1
summ_reprod = DataFrame(
	"Techniques" => ["Integrated", "SWCS", "Agatston"],
	"R Correlation Coefficient" => r_squared_values_reprod,
	"RMSE" => rmse_values_reprod,
	"RMSD" => rmsd_values_reprod
)

# ╔═╡ fb3306a7-b593-408d-94f9-0827bcfa5275
md"""
### Sensitivity
"""

# ╔═╡ 001f2158-e8a7-47fa-9e30-74aec8f89f5c
zero_cac_num = [
	string(total_zero_i, "/", total_cac)
	string(total_zero_vf, "/", total_cac)
	string(total_zero_s, "/", total_cac)
	string(num_zero_a, "/", total_cac)
]

# ╔═╡ a9dbe651-2067-4926-8fd1-4b07dc26e0df
zero_cac_perc = [
	round(total_zero_i / total_cac * 100, digits=2)
	round(total_zero_vf / total_cac * 100, digits=2)
	round(total_zero_s / total_cac * 100, digits=2)
	round(num_zero_a / total_cac * 100, digits=2)
]

# ╔═╡ 2f6d2c8f-b3aa-4303-8e74-3b8daa550058
summ_zero_cac = DataFrame(
	"Techniques" => ["Integrated", "Volume Fraction", "Spatially Weighted", "Agatston"],
	"False Negatives (CAC=0)" => zero_cac_num,
	"Percentage False Negatives (%)" => zero_cac_perc
)

# ╔═╡ Cell order:
# ╠═c7043dc2-8d49-45db-8f89-2490094b60a6
# ╠═fda1f661-91bb-4abe-83ba-200095e7efc6
# ╟─6658974b-e177-45dd-ba28-03e54440296e
# ╟─44804164-d4fd-42a6-8018-3032f2929a1a
# ╟─7f10366d-cc73-4ca7-9ebe-14c360cb8d34
# ╟─cd1d88db-fc7b-4790-b57e-d5c6d34a9657
# ╟─e201ebfc-17f2-4be8-83a6-06578ae0ee77
# ╟─62e3c411-70e5-4453-9226-2d11cb3c91fd
# ╟─f4927d17-9916-432d-9880-47ac4f83e859
# ╠═fcd9c21d-2f6c-4421-bad1-bb5c239329fb
# ╠═f8f7f2bf-7afd-474a-ba0d-9c45f16c5ddd
# ╠═fe082455-1797-4d28-bd05-621fc9dcc540
# ╠═460410a9-2dc4-42ef-b1a8-997e04727e0f
# ╠═c7abbd34-d2df-403e-96a4-764ebcb57702
# ╟─d9fb7c5c-16df-4c97-92c7-baabe1905c0f
# ╠═cbde2904-fbcc-4863-9bce-fe448fdc5ea7
# ╠═b3bf2fb7-fa89-4efb-b6e2-65c5baff8b4d
# ╠═3bc73f74-e009-4f2d-ba1b-32539dfa8505
# ╠═9cf1580b-8fd6-4881-88fd-25fd60f48807
# ╠═f85b9c18-8c65-47d7-9949-e424e6cb2e29
# ╟─cf8b2459-1204-4035-8f8c-79e92a4bab3e
# ╠═50413fca-f3ec-415d-a83e-d118beaf5009
# ╠═0fdf77f9-7906-49d3-8556-0903ba76301f
# ╠═f9267f12-3480-4b9d-ab22-68daf4b8077e
# ╠═5aae2475-846f-49ea-bbe8-87f3af101999
# ╠═9a21414e-d55f-4094-b28e-bcd72b14c3e3
# ╟─ede89d4a-54b8-46c8-b6a4-aa9c4c7137df
# ╠═b95045b3-e108-4464-8fb3-1bf3ca4eea45
# ╠═fa637952-26e2-4bf3-8184-79200cc0621f
# ╠═e7c5bc2b-d6fd-493b-8aec-32167c9701ca
# ╠═a8d566eb-e4a8-4ab5-959e-9ae9915267c1
# ╠═b28a4846-7ca3-4827-b1be-c0f6951050fc
# ╟─9db9fad7-8982-4f8c-8f4f-3196fc3f1af0
# ╠═f7bb8133-456f-49ad-9690-a9b42e39144f
# ╟─7e55db13-2b25-4283-bc42-411af0af0a00
# ╟─e5366f83-eefc-41e1-9cda-ec9303808c20
# ╟─7af4578d-f999-4ee1-9d83-994539399b92
# ╟─0f086c11-fb86-48c3-a6bd-844207054dc9
# ╟─0220f0a6-426a-4021-82c0-218793708ede
# ╟─f2777887-fe96-4795-b0d8-8430749846fe
# ╠═7d551a01-ecdf-473b-84d5-db310e65d4a1
# ╠═5986e921-5e2a-4084-92a2-96aa2297edfd
# ╠═1e551f4e-8381-4ae1-980a-a660db4dd333
# ╟─9e1ac65a-927a-45f9-a061-1f8ce137f1bb
# ╠═704bf7f2-e2ac-4ccf-9efe-dd5bc939263e
# ╠═9d12a490-4225-45cf-8315-5275d01fbc7d
# ╠═55a71f6d-e3af-4fce-bb49-f9a68831d8f0
# ╟─0bd48352-a9a6-4b9c-a85d-31d2f99e96c6
# ╟─047a921c-2c65-4a1a-8f9c-a7d1ade380ce
# ╠═99426eaa-b5a0-48aa-9731-843679fc74ab
# ╠═bc23bb1f-cb6d-4a80-82bd-9c587d29d39c
# ╟─9d631928-295f-449d-91f5-7b52e2b059e2
# ╠═5f5017c0-26f6-408a-98a6-f17fd436bd4e
# ╠═bb042948-668f-45b1-a954-4b694e64ceaf
# ╠═168756a1-1965-40a4-ae5a-c11af0624c0a
# ╠═67c093c7-4352-48d1-9fd1-634cefbea552
# ╠═6ea59087-3a3f-41c9-8b74-5f892d83b46d
# ╠═fd331a11-7ec5-45d1-aad1-521824c886d1
# ╠═9ee37984-894a-43cf-bc2e-e2a04db22bb2
# ╠═ed01e6bb-3a7d-4bb7-8995-920a22fc5f48
# ╟─a42a3104-1858-4ab0-9de0-79c0d85bdaa8
# ╠═471b8219-1268-4efd-8a94-0292336938f4
# ╠═f95a03f3-be7d-4ff0-af37-165a5405c552
# ╠═5ce9729f-82fc-43a9-b809-699702255a15
# ╠═54f85688-77d1-4e91-acec-55c4c03c721a
# ╠═28488ced-f472-4e36-956e-988ef3b06361
# ╠═d48c2606-4b79-4150-9669-2a6efc036280
# ╠═42b5e2a8-71fe-4813-be24-255b9ed64330
# ╠═2b839a05-7007-4f9f-ac3f-6bbf51c5071c
# ╟─c86b93dc-18da-47b7-8099-b27396bdd48a
# ╠═be39ddd8-16ee-4049-b61c-e16f85e6f289
# ╠═8b208615-24dc-4c53-8cfa-f3ffd419dea3
# ╠═316ce953-f138-4d32-a7bb-385ef932c8b2
# ╠═03c97478-5f88-42fb-bf29-f01d1effd3eb
# ╠═bea7c372-8aa8-4feb-ba6b-0bb3e624aa2c
# ╠═d49df961-f841-435c-9c2a-c4cf293cff32
# ╠═4a88479f-f53a-4923-9a1b-141878a2b331
# ╠═ca125df2-4b7a-4056-8f61-bba9b3dc79eb
# ╟─d5a746f9-1351-40fd-b0c6-c7f5c9b6fdd4
# ╠═a6ab095b-f34b-4e54-b752-d606ccd20b4a
# ╠═cd62f88b-2054-4eb6-8c76-042af27d3b68
# ╠═e9f69c55-3a91-4aa6-9113-143ccc495f4d
# ╠═0bd80638-ed4c-4e13-bdcb-77e880453429
# ╠═26ce1629-f38c-4dd9-8ea4-5d49dff5eb7c
# ╠═45dc885c-f00d-4873-939b-2f5f70c3fa18
# ╠═3544bf0e-d154-4ccc-8696-4d9a7c954405
# ╠═03c62590-3615-4c06-b9f0-77267308c1a5
# ╠═26904b95-eec0-4c10-b846-f18b66322a24
# ╟─3114afe2-6d7f-466b-bbac-b2f8b9ed79a5
# ╟─668b632c-bb5c-4b79-bd4f-69f6642a0bc8
# ╠═eb5dd9cf-70fb-415e-930a-e95393e15985
# ╠═fb91c84e-06f9-49aa-bf53-d6f8df85cdb4
# ╠═e78bbfc6-8699-4488-bd3a-488486b3c921
# ╟─0c282db8-3e45-49bc-9183-f27412fdd273
# ╠═b0ec93ef-5833-4d20-a8fe-959e97588816
# ╠═3c27fc44-e4bb-440b-b88b-3da8284dc329
# ╠═bece0a05-e843-4be6-94d7-dd46c127cc92
# ╠═4f934e78-117e-4054-a101-ea04731b769f
# ╠═26be4298-b7c3-4426-9ee9-aaaaa89b8d87
# ╠═f2a8bd87-67c2-4e7e-8bce-6bf7e65cbe48
# ╟─58447cc6-384a-4dc5-b733-1f8a6bfdbea1
# ╠═6328e32a-5fa3-4de8-b39e-58f800730fd5
# ╠═2e170bad-e903-483f-9cc9-5f7b2b8b2b17
# ╟─19aa4b7b-08aa-48fb-8c0c-5d85fdfd9f1d
# ╠═15f718ae-39d5-4aea-bf35-513d67d66710
# ╠═ad844c27-cdd1-49d3-b71b-08bc38b628cf
# ╟─9af1bd1a-1ce4-4795-bf64-232fdabd3727
# ╠═d0a14e60-7f95-4fad-b132-8a3e90b1bd2a
# ╠═8389c664-7da1-4df1-abd2-3f9dc4ab2775
# ╟─47ae09c1-ddf3-487b-9463-626c15938320
# ╟─d2829dd3-8e08-48ea-8636-f78bc203541d
# ╟─0823bd5e-c106-4ef9-839c-8452a4ee8568
# ╟─de822771-48f9-4620-a998-3d4235a5501d
# ╠═0c95bec7-f6cb-4c2b-91a2-275e05d9e572
# ╟─4eabcaae-55ea-4c3c-b276-15ee6f7ff9ec
# ╟─b05be082-6174-469b-be11-0b0c7d95f1f5
# ╟─494435e9-e008-4284-9132-bdc78cf999e3
# ╟─d0b3d0fd-4db6-4d09-8005-733191b70076
# ╠═97798f97-3912-48a8-8ffe-c05cdb736e09
# ╟─dbcfe49b-03c6-41f9-8003-9d054e782127
# ╟─88b83264-0e98-43cc-be25-b0ad0311aebc
# ╟─1a12b990-db62-475f-819d-c76c9deadba7
# ╟─7cabd9f5-e450-46e4-ad32-97e9c2b143d0
# ╟─e7461e35-bfe2-4b71-b785-02532e5ecf48
# ╟─596f29f3-fe65-4349-95e1-e2e42dd22c91
# ╠═ab6f7193-d8f2-4288-9c65-41e392552e84
# ╟─90593bbe-4143-450e-b780-1bdd5ee88512
# ╟─244f3279-7a1e-45e3-b046-665331f28274
# ╟─c314868f-5f42-4211-9e35-20895f55c494
# ╟─e14ed023-83de-43ac-9db2-38985ccfb4cc
# ╠═d245eb36-d4ba-4f43-a119-a49defdc85d8
# ╟─a17cd0a7-5164-446c-bfc6-4d097bd77063
# ╟─e83e986c-3840-4ece-a578-c78ae06a1df9
# ╟─438d943e-0f45-4bd5-843d-58d92aa1e33e
# ╟─0c73850e-c673-43b4-b623-288ef3084de3
# ╠═eac67a26-c28d-4607-a4f2-e586d75837a1
# ╟─fb3306a7-b593-408d-94f9-0827bcfa5275
# ╟─001f2158-e8a7-47fa-9e30-74aec8f89f5c
# ╟─a9dbe651-2067-4926-8fd1-4b07dc26e0df
# ╠═2f6d2c8f-b3aa-4303-8e74-3b8daa550058

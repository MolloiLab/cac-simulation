### A Pluto.jl notebook ###
# v0.19.16

using Markdown
using InteractiveUtils

# ╔═╡ 82a02e07-06d4-4e15-924d-c515b76d17f8
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.activate(".")

    using PlutoUI, Statistics, CSV, DataFrames, GLM, CairoMakie, HypothesisTests, Colors, MLJBase
	using StatsBase: quantile!, rmsd
end

# ╔═╡ 90bd3c98-d3c1-4209-8692-256b91679a62
TableOfContents()

# ╔═╡ fc91cc88-b137-49aa-a101-1049b8721a94
md"""
#### Helper Functions
"""

# ╔═╡ d68438a7-5b1c-4c49-8fcf-9e60b68f0aee
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

# ╔═╡ 522a9c57-d9b4-4c03-889e-d2a78a21ca5d
function prepare_linear_regression(df_large, df_medium, df_small, df_reprod_large, df_reprod_medium, df_reprod_small)
	local r_array
	local calc_array
	try
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
	catch
		r_array = [
	        Tuple(df_reprod_large[!, :calculated_swcs_large])...,
	        Tuple(df_reprod_medium[!, :calculated_swcs_medium])...,
	        Tuple(df_reprod_small[!, :calculated_swcs_small])...,
	    ]
	    calc_array = [
	        Tuple(df_large[!, :calculated_swcs_large])...,
	        Tuple(df_medium[!, :calculated_swcs_medium])...,
	        Tuple(df_small[!, :calculated_swcs_small])...,
	    ]
	end
    data = DataFrame(; X=r_array, Y=calc_array)
    model = lm(@formula(Y ~ X), data)
    r_squared = GLM.r2(model)
    rms_values = [rms(data[!, :X], data[!, :Y]), rmsd(data[!, :Y], GLM.predict(model))]

	x = DataFrame(; X=collect(1:1000))
    fitted_line = GLM.predict(model, x)
	coefficient = coef(model)

	return r_squared, rms_values, fitted_line, coefficient
end

# ╔═╡ b302881f-99cd-4ac6-b22f-6d6d5c098e27
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

# ╔═╡ 5d6c0929-6608-4a01-877b-b70a3a8dbdfe
FIGURE_PATH = "physical"

# ╔═╡ 2c8d827a-10bc-48a9-8657-c7c8d83fd5df
md"""
# Load CSVs
"""

# ╔═╡ 312b7073-caa1-44eb-9c53-7d2508a24603
md"""
## Integrated
"""

# ╔═╡ 2e67c175-5dfa-4024-a990-5dfcc85c127b
path_integrated = "/Users/daleblack/Google Drive/dev/MolloiLab/cac-simulation/output_new/integrated";

# ╔═╡ 96d711cf-d365-405c-98ec-b43f5b329d24
df_i = CSV.read(string(path_integrated, "/physical.csv"), DataFrame);

# ╔═╡ 523f9106-0431-45c1-9e6e-f4b8a84c069e
(df_i_small_1, df_i_small_2), (df_i_large_1, df_i_large_2) = groupby(df_i[1:6, :], :scan), groupby(df_i[7:end, :], :scan);

# ╔═╡ b98030b9-047e-4f60-9d20-c6cadf2aac99
df_i_1, df_i_2 = vcat(df_i_small_1, df_i_large_1, cols=:union), vcat(df_i_small_2, df_i_large_2, cols=:union);

# ╔═╡ 205dd89c-f180-4733-8b79-b38b66047689
md"""
## Agatston
"""

# ╔═╡ f69fc8e8-0f63-4a45-83cd-538e92c8f838
path_agat = "/Users/daleblack/Google Drive/dev/MolloiLab/cac-simulation/output_new/agatston";

# ╔═╡ 93fbe90d-7c8a-42c7-a184-82b039437e14
df_a = CSV.read(string(path_agat, "/physical.csv"), DataFrame);

# ╔═╡ ba0f04d1-0f01-46a0-9d6d-affddc9601a1
(df_a_small_1, df_a_small_2), (df_a_large_1, df_a_large_2) = groupby(df_a[1:6, :], :scan), groupby(df_a[7:end, :], :scan);

# ╔═╡ 6fbd654f-1bc2-49ee-86a6-e3ad4ecf38c8
df_a_1, df_a_2 = vcat(df_a_small_1, df_a_large_1, cols=:union), vcat(df_a_small_2, df_a_large_2, cols=:union);

# ╔═╡ 5eab9d9b-d3f6-4701-bfbd-b27221183f40
md"""
## Spatially Weighted
"""

# ╔═╡ 4e9a8c62-a80f-4e6a-8da9-b007ccce6fef
path_swcs = "/Users/daleblack/Google Drive/dev/MolloiLab/cac-simulation/output_new/swcs";

# ╔═╡ 4acec5a7-2ab3-4fb4-9c5e-85e6b36aa31e
df_s = CSV.read(string(path_swcs, "/physical.csv"), DataFrame);

# ╔═╡ 1b166b53-0ed1-4889-a6b6-af57cf73e6ac
(df_s_small_1, df_s_small_2), (df_s_large_1, df_s_large_2) = groupby(df_s[1:6, :], :scan), groupby(df_s[7:end, :], :scan);

# ╔═╡ 61830308-aad1-425b-842e-9ecff697da8f
df_s_1, df_s_2 = vcat(df_s_small_1, df_s_large_1, cols=:union), vcat(df_s_small_2, df_s_large_2, cols=:union);

# ╔═╡ 1c88082d-186c-43df-b4df-a8f574a94dbe
md"""
## Volume Fraction
"""

# ╔═╡ c38d5c24-3a96-4310-a682-d969345b647a
path_vf = "/Users/daleblack/Google Drive/dev/MolloiLab/cac-simulation/output_new/volume_fraction";

# ╔═╡ 73e6d4e4-c5cd-4d3f-8429-caccc72990af
df_vf = CSV.read(string(path_vf, "/physical.csv"), DataFrame);

# ╔═╡ 68fdecfe-4a7a-4d5e-bd84-e676689ab40d
(df_vf_small_1, df_vf_small_2), (df_vf_large_1, df_vf_large_2) = groupby(df_vf[1:6, :], :scan), groupby(df_vf[7:end, :], :scan);

# ╔═╡ 0d9f55fd-8248-4bc1-8423-888476760961
df_vf_1, df_vf_2 = vcat(df_vf_small_1, df_vf_large_1, cols=:union), vcat(df_vf_small_2, df_vf_large_2, cols=:union);

# ╔═╡ a4fb145a-80eb-4118-804b-c72bba493b3b
md"""
# Figures
"""

# ╔═╡ b8c20e8b-e6a8-45df-9867-b345b330b831
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

# ╔═╡ b0c5bee4-a399-4415-a8a6-5313cc31e1d4
md"""
## Accuracy
"""

# ╔═╡ 9f072412-aab7-4c53-b5aa-b62aaf3ebda1
r_squared_normal_i, rms_values_normal_i, fitted_line_normal_i, coefficient_normal_i = prepare_linear_regression(df_i)

# ╔═╡ 878c24a3-42f5-473f-8a71-3200ad1f21b3
r_squared_normal_vf, rms_values_normal_vf, fitted_line_normal_vf, coefficient_normal_vf = prepare_linear_regression(df_vf)

# ╔═╡ b6632788-a235-4831-ae74-5820b97ec8dc
r_squared_normal_a, rms_values_normal_a, fitted_line_normal_a, coefficient_normal_a = prepare_linear_regression(df_a)

# ╔═╡ 9248087e-2e21-4bd5-aed7-c5eedbcdf0c1
function lin_reg_norm()
    f = Figure()

    ##-- A --##
    ax1 = Axis(f[1, 1])

    df = df_i
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
    ax1.title = "Integrated"
	
    ##-- B --##
    ax2 = Axis(f[2, 1])

    df3 = df_vf
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
    ax2.title = "Volume Fraction"


    ##-- C --##
    ax2 = Axis(f[3, 1])

    df3 = df_a
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
    ax2.title = "Agatston"

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

    save(joinpath(dirname(pwd()),"figures", FIGURE_PATH, "accuracy.png"), f)
    return f
end

# ╔═╡ e7e38a45-aed8-4e25-8fdc-d08cb1b3f197
with_theme(medphys_theme) do
    lin_reg_norm()
end

# ╔═╡ f011edbc-2547-4be8-895b-5ebe611be122
md"""
## Reproducibility
"""

# ╔═╡ 3e5c3d00-fcd3-4fc6-9a52-4f728a5ba69c
md"""
#### Integrated
"""

# ╔═╡ 3ce29e9d-3be6-41a2-9cff-2399aff77b0a
array_i_1 = hcat(
    df_i_1[!, :calculated_mass_large],
    df_i_1[!, :calculated_mass_medium],
    df_i_1[!, :calculated_mass_small],
)

# ╔═╡ 0f58b7de-28e1-4704-9a31-f56d0c7c1b18
array_i_2 = hcat(
    df_i_2[!, :calculated_mass_large],
    df_i_2[!, :calculated_mass_medium],
    df_i_2[!, :calculated_mass_small],
)

# ╔═╡ dfcbccc0-3a6f-491c-8e6f-f68a96bf8a7c
md"""
#### Volume Fraction
"""

# ╔═╡ c49cb0e1-700e-4c77-bcbb-1b815e78c902
array_vf_1 = hcat(
    df_vf_1[!, :calculated_mass_large],
    df_vf_1[!, :calculated_mass_medium],
    df_vf_1[!, :calculated_mass_small],
)

# ╔═╡ c7719179-5c77-4e6b-ad1d-86ea90ff78f9
array_vf_2 = hcat(
    df_vf_2[!, :calculated_mass_large],
    df_vf_2[!, :calculated_mass_medium],
    df_vf_2[!, :calculated_mass_small],
)

# ╔═╡ 3daa60e3-b1ba-4106-9256-3e7919c926d0
md"""
#### Agatston
"""

# ╔═╡ 37581950-10e8-4be0-95b0-a6cd53e2f556
array_a_1 = hcat(
    df_a_1[!, :calculated_mass_large],
    df_a_1[!, :calculated_mass_medium],
    df_a_1[!, :calculated_mass_small],
)

# ╔═╡ a144cec4-0c90-4713-abbc-23795b3d8584
array_a_2 = hcat(
    df_a_2[!, :calculated_mass_large],
    df_a_2[!, :calculated_mass_medium],
    df_a_2[!, :calculated_mass_small],
)

# ╔═╡ 666b7001-c521-4bb2-ac40-20966fb9ada5
md"""
#### SWCS
"""

# ╔═╡ 7e250de1-bb94-4a8a-8d1c-bf0479cc40a5
mean_background_s = mean(df_s[:, :swcs_bkg])

# ╔═╡ e912ad90-9477-4253-97db-f679d2c9abb3
function remove_false_negatives(df, array, df_reprod, array_reprod, swcs::Bool)
	if swcs == true
		idxs_r_large = Tuple(findall(x -> x <= mean_background_s, array_reprod[:, 1]))
	    idxs_large = Tuple(findall(x -> x <= mean_background_s, array[:, 1]))
	    global indxs_large_tot = Tuple(unique([idxs_r_large..., idxs_large...]))
	
		idxs_r_med = Tuple(findall(x -> x <= mean_background_s, array_reprod[:, 2]))
	    idxs_med = Tuple(findall(x -> x <= mean_background_s, array[:, 2]))
	    indxs_med_tot = Tuple(unique([idxs_r_med..., idxs_med...]))
	
		idxs_r_small = Tuple(findall(x -> x <= mean_background_s, array_reprod[:, 3]))
	    idxs_small = Tuple(findall(x -> x <= mean_background_s, array[:, 3]))
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

# ╔═╡ 7e65b2f7-6e24-4f2c-8c2d-3c6eebc41c75
df_i_large, df_i_r_large, df_i_medium, df_i_r_medium, df_i_small, df_i_r_small = remove_false_negatives(df_i_1, array_i_1, df_i_2, array_i_2);

# ╔═╡ b7749170-0017-4e65-8bb6-2d95da3c783d
r_squared_reprod_i, rms_values_reprod_i, fitted_line_reprod_i, coefficient_reprod_i = prepare_linear_regression(df_i_r_large, df_i_r_medium, df_i_r_small, df_i_large, df_i_medium, df_i_small);

# ╔═╡ eaf446e3-503c-41b9-b746-1073f11786ed
df_vf_large, df_vf_r_large, df_vf_medium, df_vf_r_medium, df_vf_small, df_vf_r_small = remove_false_negatives(df_vf_1, array_vf_1, df_vf_2, array_vf_2);

# ╔═╡ 60844187-29ab-43a9-9970-4a41f9af3f31
r_squared_reprod_vf, rms_values_reprod_vf, fitted_line_reprod_vf, coefficient_reprod_vf = prepare_linear_regression(df_vf_r_large, df_vf_r_medium, df_vf_r_small, df_vf_large, df_vf_medium, df_vf_small);

# ╔═╡ b37915b4-aed4-4d9a-bffd-7d95387c71f2
df_a_large, df_a_r_large, df_a_medium, df_a_r_medium, df_a_small, df_a_r_small = remove_false_negatives(df_a_1, array_a_1, df_a_2, array_a_2);

# ╔═╡ 5e34bd7f-7e56-4a7c-b899-6ad334c093ed
r_squared_reprod_a, rms_values_reprod_a, fitted_line_reprod_a, coefficient_reprod_a = prepare_linear_regression(df_a_r_large, df_a_r_medium, df_a_r_small, df_a_large, df_a_medium, df_a_small)

# ╔═╡ aac1e177-5f76-496d-b6a1-254c4338f25c
array_s_1 = hcat(
    df_s_1[!, :calculated_swcs_large],
    df_s_1[!, :calculated_swcs_medium],
    df_s_1[!, :calculated_swcs_small],
)

# ╔═╡ 7c8dfdba-7b75-4812-b5f7-b0ed7ed08b31
array_s_2 = hcat(
    df_s_2[!, :calculated_swcs_large],
    df_s_2[!, :calculated_swcs_medium],
    df_s_2[!, :calculated_swcs_small],
)

# ╔═╡ 6103e58a-f454-4894-aea9-afa706ad11a6
df_s_large, df_s_r_large, df_s_medium, df_s_r_medium, df_s_small, df_s_r_small = remove_false_negatives(df_s_1, array_s_1, df_s_2, array_s_2, true);

# ╔═╡ 8c28ca4e-eff1-4921-94e1-2a04f3a06a1a
r_squared_reprod_s, rms_values_reprod_s, fitted_line_reprod_s, coefficient_reprod_s = prepare_linear_regression(df_s_r_large, df_s_r_medium, df_s_r_small, df_s_large, df_s_medium, df_s_small)

# ╔═╡ f534bb1e-3c46-476a-ae05-1384fcab1260
function reprod()
    f = Figure()

    ##-- A --##
    ax1 = Axis(f[1, 1])
    scatter!(ax1, df_i_r_large[!, :calculated_mass_large], df_i_large[!, :calculated_mass_large])
    scatter!(ax1, df_i_r_medium[!, :calculated_mass_medium], df_i_medium[!, :calculated_mass_medium])
    scatter!(
        ax1,
        df_i_r_small[!, :calculated_mass_small],
        df_i_small[!, :calculated_mass_small];
        color=:red,
    )
    lines!(ax1, [-1000, 1000], [-1000, 1000]; label="Unity")
    lines!(ax1, collect(1:1000), fitted_line_reprod_i; linestyle=:dashdot)

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

    # ##-- D --##
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

# ╔═╡ 86b5f204-d004-47dd-9709-28965692bc4f
with_theme(medphys_theme) do
    reprod()
end

# ╔═╡ b6d6696c-9d9c-4944-b936-51e877a36663
md"""
## Sensitivity and Specificity
"""

# ╔═╡ 2051d1c0-239e-4899-ba4d-9d9918ecc57f
md"""
### False Negative
"""

# ╔═╡ 3edf9057-a2c4-4696-a5f6-95b3e079ac39
# function false_negative()
#     f = Figure()
#     colors = Makie.wong_colors()

#     ##-- TOP --##
#     axtop = Axis(f[1, 1]; xticks=(1:4, ["Integrated", "Volume Fraction", "Spatially Weighted", "Agatston"]))

#     table = [1, 2, 3, 4]
#     heights1 = [
# 		(total_zero_i / total_cac) * 100,
# 		(total_zero_vf / total_cac) * 100,
#         (total_zero_s / total_cac) * 100,
#         (num_zero_a / total_cac) * 100,
#     ]
#     barplot!(axtop, table, heights1; color=colors[1:4], bar_labels=:y)

#     axtop.title = "False-Negative Scores (CAC=0)"
#     axtop.ylabel = "% False-Negative Zero CAC Scores"
#     ylims!(axtop; low=0, high=100)
#     axtop.yticks = [0, 25, 50, 75, 100]

#     save(joinpath(dirname(pwd()),"figures", FIGURE_PATH, "false_negative.png"), f)
#     return f
# end

# ╔═╡ fb286880-9f84-4bcf-acd4-69ea0afd7eb9
# with_theme(medphys_theme) do
#     false_negative()
# end

# ╔═╡ 28e60f4d-a0d0-44d9-81f2-68f945b0b103
md"""
#### SWCS
"""

# ╔═╡ fc4eb83f-9882-4873-b2de-90db0ab77b68
array_s = Array(hcat(df_s[!, :calculated_swcs_large], df_s[!, :calculated_swcs_medium], df_s[!, :calculated_swcs_small]));

# ╔═╡ 029f2477-dc9a-477d-97e1-ea2f20e92082
total_cac = length(array_s)

# ╔═╡ d332aabf-c47a-4287-be7f-1c13d67a5c54
mean_s, std_s = mean(df_s[!, :swcs_bkg]), std(df_s[!, :swcs_bkg])

# ╔═╡ 579f650b-fb90-4c92-96ed-15fbf9f70626
total_zero_s = length(findall(x -> x <= mean_s + std_s, array_s))

# ╔═╡ 64f9e589-6a08-4f13-8d96-4e025fa44eb8
md"""
#### Agatston
"""

# ╔═╡ b58db5d8-1cc5-4ba7-9712-22a436a43fdf
array_a = hcat(df_a[!, :calculated_mass_large], df_a[!, :calculated_mass_medium], df_a[!, :calculated_mass_small]);

# ╔═╡ 04f4bbf5-1c40-4b6d-9c63-936d7e13d68e
num_zero_a = length(findall(x -> x <= 0, array_a))

# ╔═╡ 716b8ee3-bad8-4e04-bc66-8f332f8e7d67
md"""
#### Integrated
"""

# ╔═╡ dbcc68d3-89cf-464b-98f0-1babcc608cc2
array_i = hcat(df_i[!, :calculated_mass_large], df_i[!, :calculated_mass_medium], df_i[!, :calculated_mass_small]);

# ╔═╡ 6053ad1a-90a4-4d7f-a96a-ed31487906da
mean_i, std_i = mean(df_i[!, :mass_bkg]), std(df_i[!, :mass_bkg])

# ╔═╡ 6cee3262-94f6-424d-9d2c-dd4482c72bd1
total_zero_i = length(findall(x -> x <= mean_i + std_i, array_i))

# ╔═╡ d15a6d6c-cae7-47ac-920e-e483dc8ef824
md"""
#### Volume Fraction
"""

# ╔═╡ 7ef0bb08-5a6a-412c-b8eb-9830aa1cee3c
array_vf = hcat(df_vf[!, :calculated_mass_large], df_vf[!, :calculated_mass_medium], df_vf[!, :calculated_mass_small]);

# ╔═╡ b4510b14-68bd-40b9-84d8-37d0645d27b1
mean_vf, std_vf = mean(df_vf[!, :mass_bkg]), std(df_vf[!, :mass_bkg])

# ╔═╡ 66c2f000-669d-4ec4-bdd1-1c585646157b
total_zero_vf = length(findall(x -> x <= mean_i + std_i, array_i))

# ╔═╡ 3536bec9-a79b-4473-b98d-d30b714cf844
total_zero_i, total_zero_vf, total_zero_s, num_zero_a

# ╔═╡ 7349fb1a-5686-40b0-bb11-c48ab4c54630
md"""
### False Positive
"""

# ╔═╡ 5e2164bc-945e-48c2-b58f-2009102eeff5
# function false_positive()
#     f = Figure()
#     colors = Makie.wong_colors()

#     ##-- TOP --##
#     axtop = Axis(f[1, 1]; xticks=(1:4, ["Integrated", "Volume Fraction", "Spatially Weighted", "Agatston"]))

#     table = [1, 2, 3, 4]
#     heights1 = [
# 		(total_zero_i_pos / total_cac_pos) * 100,
# 		(total_zero_vf_pos / total_cac_pos) * 100,
#         (total_zero_s_pos / total_cac_pos) * 100,
#         (total_zero_a_pos / total_cac_pos) * 100,
#     ]
#     barplot!(axtop, table, heights1; color=colors[1:4], bar_labels=:y)

#     axtop.title = "False-Positive Scores (CAC>0)"
#     axtop.ylabel = "% False-Positive CAC Scores"
#     ylims!(axtop; low=0, high=100)
#     axtop.yticks = [0, 25, 50, 75, 100]

#     save(joinpath(dirname(pwd()),"figures", FIGURE_PATH, "false_positive.png"), f)
#     return f
# end

# ╔═╡ 95f8d3bb-3a0b-409e-a6f0-1795ff15e3d5
# with_theme(medphys_theme) do
#     false_positive()
# end

# ╔═╡ 0a54b027-4617-47c8-a49d-a6ff5908f24b
md"""
#### SWCS
"""

# ╔═╡ 2fe30c9a-0df9-43c6-b0fb-41b7a14a8cc6
array_s_pos = df_s[!, :swcs_bkg]

# ╔═╡ ec8e9942-5c19-4be0-ada8-d534785d2e13
total_cac_pos = length(array_s_pos)

# ╔═╡ 37d70a7a-ad14-4cbb-9dc6-00d4e5f42064
total_zero_s_pos = length(findall(x -> x >= mean_s + std_s, array_s_pos))

# ╔═╡ 755e06f9-e3a3-4d9e-a264-420ef669aa01
md"""
#### Agatston
"""

# ╔═╡ 87e15e87-9fc6-41c3-94e4-8e1dc1135bb4
array_a_pos = df_a[!, :mass_bkg]

# ╔═╡ fc881c54-ff32-422a-8d24-c0e963526bad
total_zero_a_pos = length(findall(x -> x > 0, array_a_pos))

# ╔═╡ 83b4cfe0-7a04-446f-8fb8-05138b5788a4
md"""
#### Integrated
"""

# ╔═╡ ee252446-3a17-4f24-aab8-44b94cd11916
array_i_pos = df_i[!, :mass_bkg]

# ╔═╡ 001b1139-ad4b-47ce-aa00-41f4896e272f
total_zero_i_pos = length(findall(x -> x >= mean_i + std_i, array_i_pos))

# ╔═╡ 82c9c6b6-97fd-4407-b6ec-049a30d319dd
md"""
#### Volume Fraction
"""

# ╔═╡ 85ba8eb8-3ee8-49e1-926d-fe4c59f51cb7
array_vf_pos = df_vf[!, :mass_bkg]

# ╔═╡ 916425e6-b276-47f3-a5a8-c51854ef5a44
total_zero_vf_pos = length(findall(x -> x >= mean_vf + std_vf, array_vf_pos))

# ╔═╡ 6bdc877f-b3ad-476e-82c6-5164fd744586
function sensitivity_specificity()
    f = Figure()
    colors = Makie.wong_colors()

    ##-- TOP --##
    ax1 = Axis(f[1, 1]; xticks=(1:4, ["Integrated", "Volume Fraction", "Spatially Weighted", "Agatston"]))

    table = [1, 2, 3, 4]
    heights1 = [
		(total_zero_i / total_cac) * 100,
		(total_zero_vf / total_cac) * 100,
        (total_zero_s / total_cac) * 100,
        (num_zero_a / total_cac) * 100,
    ]
    barplot!(ax1, table, heights1; color=colors[1:4], bar_labels=:y)

    ax1.title = "False-Negative Scores (CAC=0)"
    ax1.ylabel = "% False-Negative Zero CAC Scores"
    ylims!(ax1; low=0, high=100)
    ax1.yticks = [0, 25, 50, 75, 100]

	##-- TOP --##
    ax2 = Axis(f[2, 1]; xticks=(1:4, ["Integrated", "Volume Fraction", "Spatially Weighted", "Agatston"]))

    table = [1, 2, 3, 4]
    heights1 = [
		(total_zero_i_pos / total_cac_pos) * 100,
		(total_zero_vf_pos / total_cac_pos) * 100,
        (total_zero_s_pos / total_cac_pos) * 100,
        (total_zero_a_pos / total_cac_pos) * 100,
    ]
    barplot!(ax2, table, heights1; color=colors[1:4], bar_labels=:y)

    ax2.title = "False-Positive Scores (CAC>0)"
    ax2.ylabel = "% False-Positive CAC Scores"
    ylims!(ax2; low=0, high=100)
    ax2.yticks = [0, 25, 50, 75, 100]

    save(joinpath(dirname(pwd()),"figures", FIGURE_PATH, "sensitivity_specificity.png"), f, resolution=(800, 800))

	
    return f
end

# ╔═╡ 92ba508c-c7cb-40be-b558-5f71f24fdea2
with_theme(medphys_theme) do
    sensitivity_specificity()
end

# ╔═╡ 0362dbe6-6f9c-4b7c-9ff4-d900efc659fb
md"""
# Tables
"""

# ╔═╡ 0a054afc-0787-4956-95d9-135bc50d8dfc
md"""
## Simulation Parameters
"""

# ╔═╡ ef3b994c-25e0-48b4-8757-61df5ab1871d
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

# ╔═╡ 342fb4d9-af12-49ca-8180-74f907d46d7e
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

# ╔═╡ cbb22572-c63e-4cbf-9e70-c01665a86e32
simulation_parameters = DataFrame(
	"Parameter" => parameters,
	"Simulation" => simulations
)

# ╔═╡ 717f0e18-40bd-4c1c-810f-3303fae5c511
md"""
## Spatially Weighted Calcium Scoring Means & Stds
"""

# ╔═╡ e683b065-9558-4c39-901c-5a3679e3ce54
kvp = [
	80
	100
	120
	135
]

# ╔═╡ c08c96d9-10c0-4f58-9875-73049a177c7e
means = [
	201.4898
	174.3658
	158.2645
	152.8815
]

# ╔═╡ afaa324a-89bb-4536-9a7d-260b7e250645
stds = [
	34.3038
	28.1708
	22.9656
	21.0778
]

# ╔═╡ cf68ca50-8bd9-47cf-a940-96a099ad5e7a
means_stds = DataFrame(
	"kVp" => kvp,
	"Mean" => means,
	"Std" => stds
)

# ╔═╡ 1929ca2c-c587-4c9c-9ed2-2beaf5824978
md"""
## Summaries
"""

# ╔═╡ 01e3fccb-31c6-454e-874d-f61c5dd1af94
md"""
### Accuracy
"""

# ╔═╡ 38abb4c4-5616-4c17-9d2e-be2a93faf5e5
r_squared_values_normal = [
	r_squared_normal_i,
	r_squared_normal_vf,
	r_squared_normal_a,
]

# ╔═╡ c391965f-b842-4988-a93e-b3dde2aa9796
rmse_values_normal = [
	rms_values_normal_i[1],
	rms_values_normal_vf[1],
	rms_values_normal_a[1]
]

# ╔═╡ c4c4707e-4a51-489b-b0eb-f2e3639205c3
rmsd_values_normal = [
	rms_values_normal_i[2],
	rms_values_normal_vf[2],
	rms_values_normal_a[2]
]

# ╔═╡ b110bc4c-db2a-4011-8d8b-b4ada36385bb
summ_regression_normal = DataFrame(
	"Techniques" => ["Integrated", "Volume Fraction", "Agatston"],
	"R Correlation Coefficient" => r_squared_values_normal,
	"RMSE" => rmse_values_normal,
	"RMSD" => rmsd_values_normal
)

# ╔═╡ df1f26e8-4e2f-4bec-87cf-7f25fb4f0071
md"""
### Reproducibility
"""

# ╔═╡ ae9a858d-2866-4109-82fc-e764d6381833
r_squared_values_reprod = [
	r_squared_reprod_i
	r_squared_reprod_a
	r_squared_reprod_s
]

# ╔═╡ ed56b7c2-a4ab-4b06-9025-bb54dc150743
rmse_values_reprod = [
	rms_values_reprod_i[1]
	rms_values_reprod_a[1]
	rms_values_reprod_s[1]
]

# ╔═╡ 507a2819-ab1b-4db2-a54f-44e75aa91f02
rmsd_values_reprod = [
	rms_values_reprod_i[2]
	rms_values_reprod_a[2]
	rms_values_reprod_s[2]
]

# ╔═╡ 56906cf8-2473-41f5-914a-c93d61e4642f
summ_reprod = DataFrame(
	"Techniques" => ["Integrated", "SWCS", "Agatston"],
	"R Correlation Coefficient" => r_squared_values_reprod,
	"RMSE" => rmse_values_reprod,
	"RMSD" => rmsd_values_reprod
)

# ╔═╡ 108a5151-e07a-4dc3-91f7-c8e355617dfd
md"""
### Sensitivity
"""

# ╔═╡ e1cb16f1-c4ea-412c-931b-2761f39765dc
zero_cac_num = [
	string(total_zero_i, "/", total_cac)
	string(total_zero_vf, "/", total_cac)
	string(total_zero_s, "/", total_cac)
	string(num_zero_a, "/", total_cac)
]

# ╔═╡ 3ed1fc03-2f3a-4b22-9e9a-293670b43e63
zero_cac_perc = [
	round(total_zero_i / total_cac * 100, digits=2)
	round(total_zero_vf / total_cac * 100, digits=2)
	round(total_zero_s / total_cac * 100, digits=2)
	round(num_zero_a / total_cac * 100, digits=2)
]

# ╔═╡ a8e1fa2a-6e51-491a-b805-c7fdc17f2dda
summ_zero_cac = DataFrame(
	"Techniques" => ["Integrated", "Volume Fraction", "Spatially Weighted", "Agatston"],
	"False Negatives (CAC=0)" => zero_cac_num,
	"Percentage False Negatives (%)" => zero_cac_perc
)


# ╔═╡ Cell order:
# ╠═82a02e07-06d4-4e15-924d-c515b76d17f8
# ╠═90bd3c98-d3c1-4209-8692-256b91679a62
# ╟─fc91cc88-b137-49aa-a101-1049b8721a94
# ╟─d68438a7-5b1c-4c49-8fcf-9e60b68f0aee
# ╟─522a9c57-d9b4-4c03-889e-d2a78a21ca5d
# ╟─b302881f-99cd-4ac6-b22f-6d6d5c098e27
# ╟─e912ad90-9477-4253-97db-f679d2c9abb3
# ╠═5d6c0929-6608-4a01-877b-b70a3a8dbdfe
# ╟─2c8d827a-10bc-48a9-8657-c7c8d83fd5df
# ╟─312b7073-caa1-44eb-9c53-7d2508a24603
# ╠═2e67c175-5dfa-4024-a990-5dfcc85c127b
# ╠═96d711cf-d365-405c-98ec-b43f5b329d24
# ╠═523f9106-0431-45c1-9e6e-f4b8a84c069e
# ╠═b98030b9-047e-4f60-9d20-c6cadf2aac99
# ╟─205dd89c-f180-4733-8b79-b38b66047689
# ╠═f69fc8e8-0f63-4a45-83cd-538e92c8f838
# ╠═93fbe90d-7c8a-42c7-a184-82b039437e14
# ╠═ba0f04d1-0f01-46a0-9d6d-affddc9601a1
# ╠═6fbd654f-1bc2-49ee-86a6-e3ad4ecf38c8
# ╟─5eab9d9b-d3f6-4701-bfbd-b27221183f40
# ╠═4e9a8c62-a80f-4e6a-8da9-b007ccce6fef
# ╠═4acec5a7-2ab3-4fb4-9c5e-85e6b36aa31e
# ╠═1b166b53-0ed1-4889-a6b6-af57cf73e6ac
# ╠═61830308-aad1-425b-842e-9ecff697da8f
# ╟─1c88082d-186c-43df-b4df-a8f574a94dbe
# ╠═c38d5c24-3a96-4310-a682-d969345b647a
# ╠═73e6d4e4-c5cd-4d3f-8429-caccc72990af
# ╠═68fdecfe-4a7a-4d5e-bd84-e676689ab40d
# ╠═0d9f55fd-8248-4bc1-8423-888476760961
# ╟─a4fb145a-80eb-4118-804b-c72bba493b3b
# ╠═b8c20e8b-e6a8-45df-9867-b345b330b831
# ╟─b0c5bee4-a399-4415-a8a6-5313cc31e1d4
# ╠═9f072412-aab7-4c53-b5aa-b62aaf3ebda1
# ╠═878c24a3-42f5-473f-8a71-3200ad1f21b3
# ╠═b6632788-a235-4831-ae74-5820b97ec8dc
# ╟─9248087e-2e21-4bd5-aed7-c5eedbcdf0c1
# ╟─e7e38a45-aed8-4e25-8fdc-d08cb1b3f197
# ╟─f011edbc-2547-4be8-895b-5ebe611be122
# ╟─f534bb1e-3c46-476a-ae05-1384fcab1260
# ╟─86b5f204-d004-47dd-9709-28965692bc4f
# ╟─3e5c3d00-fcd3-4fc6-9a52-4f728a5ba69c
# ╠═3ce29e9d-3be6-41a2-9cff-2399aff77b0a
# ╠═0f58b7de-28e1-4704-9a31-f56d0c7c1b18
# ╠═7e65b2f7-6e24-4f2c-8c2d-3c6eebc41c75
# ╠═b7749170-0017-4e65-8bb6-2d95da3c783d
# ╟─dfcbccc0-3a6f-491c-8e6f-f68a96bf8a7c
# ╠═c49cb0e1-700e-4c77-bcbb-1b815e78c902
# ╠═c7719179-5c77-4e6b-ad1d-86ea90ff78f9
# ╠═eaf446e3-503c-41b9-b746-1073f11786ed
# ╠═60844187-29ab-43a9-9970-4a41f9af3f31
# ╟─3daa60e3-b1ba-4106-9256-3e7919c926d0
# ╠═37581950-10e8-4be0-95b0-a6cd53e2f556
# ╠═a144cec4-0c90-4713-abbc-23795b3d8584
# ╠═b37915b4-aed4-4d9a-bffd-7d95387c71f2
# ╠═5e34bd7f-7e56-4a7c-b899-6ad334c093ed
# ╟─666b7001-c521-4bb2-ac40-20966fb9ada5
# ╠═7e250de1-bb94-4a8a-8d1c-bf0479cc40a5
# ╠═aac1e177-5f76-496d-b6a1-254c4338f25c
# ╠═7c8dfdba-7b75-4812-b5f7-b0ed7ed08b31
# ╠═6103e58a-f454-4894-aea9-afa706ad11a6
# ╠═8c28ca4e-eff1-4921-94e1-2a04f3a06a1a
# ╟─b6d6696c-9d9c-4944-b936-51e877a36663
# ╟─6bdc877f-b3ad-476e-82c6-5164fd744586
# ╟─92ba508c-c7cb-40be-b558-5f71f24fdea2
# ╟─2051d1c0-239e-4899-ba4d-9d9918ecc57f
# ╟─3edf9057-a2c4-4696-a5f6-95b3e079ac39
# ╟─fb286880-9f84-4bcf-acd4-69ea0afd7eb9
# ╠═3536bec9-a79b-4473-b98d-d30b714cf844
# ╟─28e60f4d-a0d0-44d9-81f2-68f945b0b103
# ╠═fc4eb83f-9882-4873-b2de-90db0ab77b68
# ╠═029f2477-dc9a-477d-97e1-ea2f20e92082
# ╠═d332aabf-c47a-4287-be7f-1c13d67a5c54
# ╠═579f650b-fb90-4c92-96ed-15fbf9f70626
# ╟─64f9e589-6a08-4f13-8d96-4e025fa44eb8
# ╠═b58db5d8-1cc5-4ba7-9712-22a436a43fdf
# ╠═04f4bbf5-1c40-4b6d-9c63-936d7e13d68e
# ╟─716b8ee3-bad8-4e04-bc66-8f332f8e7d67
# ╠═dbcc68d3-89cf-464b-98f0-1babcc608cc2
# ╠═6053ad1a-90a4-4d7f-a96a-ed31487906da
# ╠═6cee3262-94f6-424d-9d2c-dd4482c72bd1
# ╟─d15a6d6c-cae7-47ac-920e-e483dc8ef824
# ╠═7ef0bb08-5a6a-412c-b8eb-9830aa1cee3c
# ╠═b4510b14-68bd-40b9-84d8-37d0645d27b1
# ╠═66c2f000-669d-4ec4-bdd1-1c585646157b
# ╟─7349fb1a-5686-40b0-bb11-c48ab4c54630
# ╟─5e2164bc-945e-48c2-b58f-2009102eeff5
# ╟─95f8d3bb-3a0b-409e-a6f0-1795ff15e3d5
# ╟─0a54b027-4617-47c8-a49d-a6ff5908f24b
# ╠═2fe30c9a-0df9-43c6-b0fb-41b7a14a8cc6
# ╠═ec8e9942-5c19-4be0-ada8-d534785d2e13
# ╠═37d70a7a-ad14-4cbb-9dc6-00d4e5f42064
# ╟─755e06f9-e3a3-4d9e-a264-420ef669aa01
# ╠═87e15e87-9fc6-41c3-94e4-8e1dc1135bb4
# ╠═fc881c54-ff32-422a-8d24-c0e963526bad
# ╟─83b4cfe0-7a04-446f-8fb8-05138b5788a4
# ╠═ee252446-3a17-4f24-aab8-44b94cd11916
# ╠═001b1139-ad4b-47ce-aa00-41f4896e272f
# ╟─82c9c6b6-97fd-4407-b6ec-049a30d319dd
# ╠═85ba8eb8-3ee8-49e1-926d-fe4c59f51cb7
# ╠═916425e6-b276-47f3-a5a8-c51854ef5a44
# ╟─0362dbe6-6f9c-4b7c-9ff4-d900efc659fb
# ╟─0a054afc-0787-4956-95d9-135bc50d8dfc
# ╠═ef3b994c-25e0-48b4-8757-61df5ab1871d
# ╠═342fb4d9-af12-49ca-8180-74f907d46d7e
# ╠═cbb22572-c63e-4cbf-9e70-c01665a86e32
# ╟─717f0e18-40bd-4c1c-810f-3303fae5c511
# ╠═e683b065-9558-4c39-901c-5a3679e3ce54
# ╠═c08c96d9-10c0-4f58-9875-73049a177c7e
# ╠═afaa324a-89bb-4536-9a7d-260b7e250645
# ╠═cf68ca50-8bd9-47cf-a940-96a099ad5e7a
# ╟─1929ca2c-c587-4c9c-9ed2-2beaf5824978
# ╟─01e3fccb-31c6-454e-874d-f61c5dd1af94
# ╠═38abb4c4-5616-4c17-9d2e-be2a93faf5e5
# ╠═c391965f-b842-4988-a93e-b3dde2aa9796
# ╠═c4c4707e-4a51-489b-b0eb-f2e3639205c3
# ╠═b110bc4c-db2a-4011-8d8b-b4ada36385bb
# ╟─df1f26e8-4e2f-4bec-87cf-7f25fb4f0071
# ╠═ae9a858d-2866-4109-82fc-e764d6381833
# ╠═ed56b7c2-a4ab-4b06-9025-bb54dc150743
# ╠═507a2819-ab1b-4db2-a54f-44e75aa91f02
# ╠═56906cf8-2473-41f5-914a-c93d61e4642f
# ╟─108a5151-e07a-4dc3-91f7-c8e355617dfd
# ╠═e1cb16f1-c4ea-412c-931b-2761f39765dc
# ╠═3ed1fc03-2f3a-4b22-9e9a-293670b43e63
# ╠═a8e1fa2a-6e51-491a-b805-c7fdc17f2dda

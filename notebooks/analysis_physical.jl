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

# ╔═╡ 957ee3a3-52a4-47c3-9aaa-708c8391bdfd
md"""
## Sensitivity and Specificity
"""

# ╔═╡ 48569e15-133e-4cc1-8a48-a6d7cc6d013a
md"""
### False Negative
"""

# ╔═╡ bfbd8d42-773b-4d43-83d2-91f6f14f65a8
# ## ---- FIND OPTIMAL THRESHOLD FOR SENS/SPEC ---- ##
# begin
# 	cutoffs = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45]
# 	for cutoff in cutoffs
# 		global total_zero_i, total_zero_vf, total_zero_i_pos, total_zero_vf_pos

# 		total_zero_i = length(findall(x -> x <= cutoff, array_i))
# 		total_zero_vf = length(findall(x -> x <= cutoff, array_vf))
# 		total_zero_i_pos = length(findall(x -> x > cutoff, array_i_pos))
# 		total_zero_vf_pos = length(findall(x -> x > cutoff, array_vf_pos))
		
		
# 		rms_i_neg = rmsd([0], [total_zero_i])
# 		rms_vf_neg = rmsd([0], [total_zero_vf])
# 		rms_i_pos = rmsd([0], [total_zero_i_pos])
# 		rms_vf_pos = rmsd([0], [total_zero_vf_pos])

# 		tot_rms = mean([rms_i_neg, rms_vf_neg, rms_i_pos, rms_vf_pos])
# 		@show tot_rms, cutoff
# 	end
# end

# ╔═╡ b393c657-95fb-44ec-aa9a-603969a880ce
cutoff = 0.25
# cutoff = 0

# ╔═╡ 5687d65a-37fb-43f7-b8b9-eb2620fc95bf
md"""
#### SWCS
"""

# ╔═╡ c6976b6f-72f8-4e86-97d3-21de16d8dd26
array_s = hcat(df_s[!, :calculated_swcs_large], df_s[!, :calculated_swcs_medium], df_s[!, :calculated_swcs_small]);

# ╔═╡ 3503ebb9-6998-441a-8e9e-02fa75f77822
total_cac = length(array_s)

# ╔═╡ 4f517bf7-2aaf-4a92-987b-bf4c6469a777
mean_swcs = mean(df_s[:, :swcs_bkg])

# ╔═╡ dd27edbc-b64c-4b78-a715-32f8b5716d6f
begin
    total_zero_s = length(findall(x -> x < mean_swcs, array_s))
end;

# ╔═╡ 66c1ae2a-0ec4-4115-ad52-f79de9e3a1f9
md"""
#### Agatston
"""

# ╔═╡ 745360d6-9208-4609-902a-4c7c7d811331
array_a = hcat(df_a[!, :calculated_mass_large], df_a[!, :calculated_mass_medium], df_a[!, :calculated_mass_small]);

# ╔═╡ 62c5a925-6600-4dc4-9f17-99106b54648d
num_zero_a = length(findall(x -> x <= 0, array_a))

# ╔═╡ 339e6ba8-6214-4996-abc2-ea71d3ac35c5
md"""
#### Integrated
"""

# ╔═╡ 2482e89e-e8b3-48a3-b60a-04251c9f6bfe
array_i = hcat(df_i[!, :calculated_mass_large], df_i[!, :calculated_mass_medium], df_i[!, :calculated_mass_small]);

# ╔═╡ 5e763173-51af-4604-995f-57f508da4be4
total_zero_i = length(findall(x -> x <= cutoff, array_i))

# ╔═╡ 960f1234-aff7-48de-9213-6cfa850ad55f
md"""
#### Volume Fraction
"""

# ╔═╡ 3f60e3a3-9c9b-4d2c-8aeb-271a8a1ede23
array_vf = hcat(df_vf[!, :calculated_mass_large], df_vf[!, :calculated_mass_medium], df_vf[!, :calculated_mass_small]);

# ╔═╡ bf373eef-6ce0-443a-9407-b95668585ae7
total_zero_vf = length(findall(x -> x <= cutoff, array_vf))

# ╔═╡ a04874c9-f1da-42d0-b2b1-52344a60a9a2
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

# ╔═╡ 5b92674d-6628-463e-b4f5-d25aec8a23d1
with_theme(medphys_theme) do
    false_negative()
end

# ╔═╡ fd73f26b-04f3-41d1-b77a-ba9124d84a1b
total_zero_i, total_zero_vf, total_zero_s, num_zero_a

# ╔═╡ 6f9380b0-383d-4f78-9214-eb89763d6082
md"""
### False Positive
"""

# ╔═╡ 06a41250-a9c8-45b5-9198-ee01c1fde29a
md"""
#### SWCS
"""

# ╔═╡ 773eff0d-47dd-45e6-965e-d8a5e730b3a1
df_s

# ╔═╡ 4805a652-b27e-49ec-8fd5-309c845d70d3
mean_swcs

# ╔═╡ 762ef6fa-a5f6-48a7-a870-ebc906339b6c
total_zero_s_pos = length(findall(x -> x > mean_swcs, df_s[!, :swcs_bkg]))

# ╔═╡ da14ac16-ef28-4c84-8139-1d3215d1a795
total_cac_pos = length(df_s[!, :swcs_bkg])

# ╔═╡ 7c8c5693-a976-4f0f-a36f-71f54f5baf86
total_zero_s_pos, total_cac_pos

# ╔═╡ 9ddd8afd-a8d3-4e1f-a609-b84a45baab06
total_zero_s_pos / total_cac_pos

# ╔═╡ 16759fcd-65b6-4d82-8a05-781bd300637d
md"""
#### Agatston
"""

# ╔═╡ 8d4efd07-3625-402e-8dff-ae47018ba0f5
array_a_pos = df_a[!, :mass_bkg]

# ╔═╡ acf9a966-0f82-4b99-985b-3839be57aaa5
total_zero_a_pos = length(findall(x -> x > 0, array_a_pos))

# ╔═╡ 73e9a2f3-9471-4423-aabb-7d791078b116
md"""
#### Integrated
"""

# ╔═╡ 3ceedf0c-5723-4547-ba2e-58fb0b5ec26c
array_i_pos = df_i[!, :mass_bkg]

# ╔═╡ f83c7d00-81c7-41df-a8fb-75ab96700412
total_zero_i_pos = length(findall(x -> x > cutoff, array_i_pos))

# ╔═╡ 8a23605d-c130-4a05-bd1a-beac82d0351a
md"""
#### Volume Fraction
"""

# ╔═╡ 15f2614f-bdad-446a-995f-5b0759249c01
array_vf_pos = df_vf[!, :mass_bkg]

# ╔═╡ 8481f39b-e278-4d2e-93e9-a2e251739f02
total_zero_vf_pos = length(findall(x -> x > cutoff, array_vf_pos))

# ╔═╡ e6f32ffd-be2f-4281-aea2-b006a695b1ff
function false_positive()
    f = Figure()
    colors = Makie.wong_colors()

    ##-- TOP --##
    axtop = Axis(f[1, 1]; xticks=(1:4, ["Integrated", "Volume Fraction", "Spatially Weighted", "Agatston"]))

    table = [1, 2, 3, 4]
    heights1 = [
		(total_zero_i_pos / total_cac_pos) * 100,
		(total_zero_vf_pos / total_cac_pos) * 100,
        (total_zero_s_pos / total_cac_pos) * 100,
        (total_zero_a_pos / total_cac_pos) * 100,
    ]
    barplot!(axtop, table, heights1; color=colors[1:4], bar_labels=:y)

    axtop.title = "False-Positive Scores (CAC>0)"
    axtop.ylabel = "% False-Positive CAC Scores"
    ylims!(axtop; low=0, high=100)
    axtop.yticks = [0, 25, 50, 75, 100]

    save(joinpath(dirname(pwd()),"figures", FIGURE_PATH, "false_positive.png"), f)
    return f
end

# ╔═╡ 169a7728-3af1-4d02-b879-84d5f02cd37c
with_theme(medphys_theme) do
    false_positive()
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
# ╟─957ee3a3-52a4-47c3-9aaa-708c8391bdfd
# ╟─5b92674d-6628-463e-b4f5-d25aec8a23d1
# ╟─169a7728-3af1-4d02-b879-84d5f02cd37c
# ╟─48569e15-133e-4cc1-8a48-a6d7cc6d013a
# ╟─a04874c9-f1da-42d0-b2b1-52344a60a9a2
# ╠═fd73f26b-04f3-41d1-b77a-ba9124d84a1b
# ╠═3503ebb9-6998-441a-8e9e-02fa75f77822
# ╠═bfbd8d42-773b-4d43-83d2-91f6f14f65a8
# ╠═b393c657-95fb-44ec-aa9a-603969a880ce
# ╟─5687d65a-37fb-43f7-b8b9-eb2620fc95bf
# ╠═c6976b6f-72f8-4e86-97d3-21de16d8dd26
# ╠═4f517bf7-2aaf-4a92-987b-bf4c6469a777
# ╠═dd27edbc-b64c-4b78-a715-32f8b5716d6f
# ╟─66c1ae2a-0ec4-4115-ad52-f79de9e3a1f9
# ╠═745360d6-9208-4609-902a-4c7c7d811331
# ╠═62c5a925-6600-4dc4-9f17-99106b54648d
# ╟─339e6ba8-6214-4996-abc2-ea71d3ac35c5
# ╠═2482e89e-e8b3-48a3-b60a-04251c9f6bfe
# ╠═5e763173-51af-4604-995f-57f508da4be4
# ╟─960f1234-aff7-48de-9213-6cfa850ad55f
# ╠═3f60e3a3-9c9b-4d2c-8aeb-271a8a1ede23
# ╠═bf373eef-6ce0-443a-9407-b95668585ae7
# ╟─6f9380b0-383d-4f78-9214-eb89763d6082
# ╟─e6f32ffd-be2f-4281-aea2-b006a695b1ff
# ╟─06a41250-a9c8-45b5-9198-ee01c1fde29a
# ╠═773eff0d-47dd-45e6-965e-d8a5e730b3a1
# ╠═4805a652-b27e-49ec-8fd5-309c845d70d3
# ╠═762ef6fa-a5f6-48a7-a870-ebc906339b6c
# ╠═da14ac16-ef28-4c84-8139-1d3215d1a795
# ╠═7c8c5693-a976-4f0f-a36f-71f54f5baf86
# ╠═9ddd8afd-a8d3-4e1f-a609-b84a45baab06
# ╟─16759fcd-65b6-4d82-8a05-781bd300637d
# ╠═8d4efd07-3625-402e-8dff-ae47018ba0f5
# ╠═acf9a966-0f82-4b99-985b-3839be57aaa5
# ╟─73e9a2f3-9471-4423-aabb-7d791078b116
# ╠═3ceedf0c-5723-4547-ba2e-58fb0b5ec26c
# ╠═f83c7d00-81c7-41df-a8fb-75ab96700412
# ╟─8a23605d-c130-4a05-bd1a-beac82d0351a
# ╠═15f2614f-bdad-446a-995f-5b0759249c01
# ╠═8481f39b-e278-4d2e-93e9-a2e251739f02
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

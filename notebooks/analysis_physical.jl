### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 82a02e07-06d4-4e15-924d-c515b76d17f8
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.activate(".")

    using PlutoUI, Statistics, CSV, DataFrames, GLM, CairoMakie, HypothesisTests, Colors, MLJBase, PrettyTables, LinearAlgebra, Printf
	using StatsBase: quantile!, rmsd, zscore
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

# ╔═╡ 2e862b70-3c59-4572-8e71-f3233b4363fa
function create_textbox(loc, coef, r, rms; p=5)
    local tb
    if coef[1] > 0
        tb = Label(
            loc;
            # text="y = $(@sprintf "%.2f" coef[2])x+$(@sprintf "%.2f" coef[1]) \nr = $(@sprintf "%.2f" r) \nRMSE: $(@sprintf "%.2f" rms[1]) \nRMSD: $(@sprintf "%.2f" rms[2])",
			text="RMSE: $(@sprintf "%.2f" rms[1]) \nRMSD: $(@sprintf "%.2f" rms[2])",
            padding=(p, p, p, p),
            tellheight=false,
            tellwidth=false,
            halign=:left,
            justification=:left,
            valign=:top,
            fontsize=12
        )
    else
        tb = Label(
            loc,
            # text="y = $(@sprintf "%.2f" coef[2])x$(@sprintf "%.2f" coef[1]) \nr = $(@sprintf "%.2f" r) \nRMSE: $(@sprintf "%.2f" rms[1]) \nRMSD: $(@sprintf "%.2f" rms[2])",
			text="RMSE: $(@sprintf "%.2f" rms[1]) \nRMSD: $(@sprintf "%.2f" rms[2])",
            padding=(p, p, p, p),
            tellheight=false,
            tellwidth=false,
            halign=:left,
            justification=:left,
            valign=:top,
            fontsize=12
        )
    end
    return tb
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

# ╔═╡ bdecc190-cecc-4f4a-b865-742ef553ca6a
root_dir = joinpath(dirname(pwd()), "output_new")

# ╔═╡ 2e67c175-5dfa-4024-a990-5dfcc85c127b
path_integrated = joinpath(root_dir, "integrated")

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
path_agat = joinpath(root_dir, "agatston")

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
path_swcs = joinpath(root_dir, "swcs")

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
path_vf = joinpath(root_dir, "volume_fraction")

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
medphys_theme = Theme(
    Axis = (
        backgroundcolor = :white,
		xgridcolor = :gray,
		xgridwidth = 0.1,
		xlabelsize = 15,
		xticklabelsize = 20,
		ygridcolor = :gray,
		ygridwidth = 0.1,
		ylabelsize = 15,
		yticklabelsize = 20,
		bottomsplinecolor = :black,
		leftspinecolor = :black,
		titlesize = 20
	)
);

# ╔═╡ b0c5bee4-a399-4415-a8a6-5313cc31e1d4
md"""
## Accuracy
"""

# ╔═╡ d4ddfb59-5f1b-414c-87aa-aed494e008d9
begin
	r_squared_normal_i, rms_values_normal_i, fitted_line_normal_i, coefficient_normal_i = prepare_linear_regression(df_i)
	r_squared_normal_i = round.(r_squared_normal_i; digits=3)
	rms_values_normal_i = round.(rms_values_normal_i; digits=3)
	coefficient_normal_i = round.(coefficient_normal_i; digits=3)
end

# ╔═╡ fa3a6b92-d83a-443e-81f1-c6bd48f55d60
begin
	r_squared_normal_vf, rms_values_normal_vf, fitted_line_normal_vf, coefficient_normal_vf = prepare_linear_regression(df_vf)
	r_squared_normal_vf = round.(r_squared_normal_vf; digits=3)
	rms_values_normal_vf = round.(rms_values_normal_vf; digits=3)
	coefficient_normal_vf = round.(coefficient_normal_vf; digits=3)
end

# ╔═╡ 1e50121c-ac1f-4107-863b-d2da79c77700
begin
	r_squared_normal_a, rms_values_normal_a, fitted_line_normal_a, coefficient_normal_a = prepare_linear_regression(df_a)
	r_squared_normal_a = round.(r_squared_normal_a; digits=3)
	rms_values_normal_a = round.(rms_values_normal_a; digits=3)
	coefficient_normal_a = round.(coefficient_normal_a; digits=3)
end

# ╔═╡ 9248087e-2e21-4bd5-aed7-c5eedbcdf0c1
function lin_reg()
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
	create_textbox(f[1, 1], coefficient_normal_i, r_squared_normal_i, rms_values_normal_i)

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
	create_textbox(f[2, 1], coefficient_normal_vf, r_squared_normal_vf, rms_values_normal_vf)
	

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
	create_textbox(f[3, 1], coefficient_normal_a, r_squared_normal_a, rms_values_normal_a)
	

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
            fontsize=25,
            padding=(0, 90, 25, 0),
            halign=:right,
        )
    end

    save(joinpath(dirname(pwd()),"figures", FIGURE_PATH, "accuracy.png"), f)
    return f
end

# ╔═╡ e7e38a45-aed8-4e25-8fdc-d08cb1b3f197
with_theme(lin_reg, medphys_theme)

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
begin
	df_i_large, df_i_r_large, df_i_medium, df_i_r_medium, df_i_small, df_i_r_small = remove_false_negatives(df_i_1, array_i_1, df_i_2, array_i_2)

	df_i_large_norm, df_i_r_large_norm, df_i_medium_norm, df_i_r_medium_norm, df_i_small_norm, df_i_r_small_norm = mapcols(zscore, df_i_large[:, 4:end]), mapcols(zscore, df_i_r_large[:, 4:end]), mapcols(zscore, df_i_medium[:, 4:end]), mapcols(zscore, df_i_r_medium[:, 4:end]), mapcols(zscore, df_i_small[:, 4:end]), mapcols(zscore, df_i_r_small[:, 4:end])
end;

# ╔═╡ b7749170-0017-4e65-8bb6-2d95da3c783d
begin
	r_squared_reprod_i, rms_values_reprod_i, fitted_line_reprod_i, coefficient_reprod_i = prepare_linear_regression(df_i_r_large, df_i_r_medium, df_i_r_small, df_i_large, df_i_medium, df_i_small)
	
	r_squared_reprod_i = round.(r_squared_reprod_i; digits=3)
	rms_values_reprod_i = round.(rms_values_reprod_i; digits=3)
	coefficient_reprod_i = round.(coefficient_reprod_i; digits=3)
end

# ╔═╡ 0afaa9ac-0d8f-4474-ad90-aae299a873f0
begin
	r_squared_reprod_i_norm, rms_values_reprod_i_norm, fitted_line_reprod_i_norm, coefficient_reprod_i_norm = prepare_linear_regression(df_i_r_large_norm, df_i_r_medium_norm, df_i_r_small_norm, df_i_large_norm, df_i_medium_norm, df_i_small_norm)
	
	r_squared_reprod_i_norm = round.(r_squared_reprod_i_norm; digits=3)
	rms_values_reprod_i_norm = round.(rms_values_reprod_i_norm; digits=3)
	coefficient_reprod_i_norm = round.(coefficient_reprod_i_norm; digits=3)
end

# ╔═╡ a1b4c6e2-662a-4b0f-bda1-b154d4630aa8
rms_values_reprod_i_norm

# ╔═╡ 54cde88d-515b-4304-a382-165542a396c3
begin
	df_vf_large, df_vf_r_large, df_vf_medium, df_vf_r_medium, df_vf_small, df_vf_r_small = remove_false_negatives(df_vf_1, array_vf_1, df_vf_2, array_vf_2)

	df_vf_large_norm, df_vf_r_large_norm, df_vf_medium_norm, df_vf_r_medium_norm, df_vf_small_norm, df_vf_r_small_norm = mapcols(zscore, df_vf_large[:, 4:end]), mapcols(zscore, df_vf_r_large[:, 4:end]), mapcols(zscore, df_vf_medium[:, 4:end]), mapcols(zscore, df_vf_r_medium[:, 4:end]), mapcols(zscore, df_vf_small[:, 4:end]), mapcols(zscore, df_vf_r_small[:, 4:end])
end;

# ╔═╡ 2974a158-a1c7-4509-b609-60711195df12
begin
	r_squared_reprod_vf, rms_values_reprod_vf, fitted_line_reprod_vf, coefficient_reprod_vf = prepare_linear_regression(df_vf_r_large, df_vf_r_medium, df_vf_r_small, df_vf_large, df_vf_medium, df_vf_small)
	
	r_squared_reprod_vf = round.(r_squared_reprod_vf; digits=3)
	rms_values_reprod_vf = round.(rms_values_reprod_vf; digits=3)
	coefficient_reprod_vf = round.(coefficient_reprod_vf; digits=3)
end

# ╔═╡ 99e53995-61a3-45f3-a699-8bdef0dfe152
begin
	r_squared_reprod_vf_norm, rms_values_reprod_vf_norm, fitted_line_reprod_vf_norm, coefficient_reprod_vf_norm = prepare_linear_regression(df_vf_r_large_norm, df_vf_r_medium_norm, df_vf_r_small_norm, df_vf_large_norm, df_vf_medium_norm, df_vf_small_norm)
	
	r_squared_reprod_vf_norm = round.(r_squared_reprod_vf_norm; digits=3)
	rms_values_reprod_vf_norm = round.(rms_values_reprod_vf_norm; digits=3)
	coefficient_reprod_vf_norm = round.(coefficient_reprod_vf_norm; digits=3)
end

# ╔═╡ eaae39f7-21d8-488d-841f-e8a4e57175bd
rms_values_reprod_vf_norm

# ╔═╡ a087b301-3b51-44b0-8afd-74f983c26b06
begin
	df_a_large, df_a_r_large, df_a_medium, df_a_r_medium, df_a_small, df_a_r_small = remove_false_negatives(df_a_1, array_a_1, df_a_2, array_a_2)

	df_a_large_norm, df_a_r_large_norm, df_a_medium_norm, df_a_r_medium_norm, df_a_small_norm, df_a_r_small_norm = mapcols(zscore, df_a_large[:, 4:end]), mapcols(zscore, df_a_r_large[:, 4:end]), mapcols(zscore, df_a_medium[:, 4:end]), mapcols(zscore, df_a_r_medium[:, 4:end]), mapcols(zscore, df_a_small[:, 4:end]), mapcols(zscore, df_a_r_small[:, 4:end])
end;

# ╔═╡ c482a5e0-aeb4-4300-8dd1-79e99ac688ee
begin
	r_squared_reprod_a, rms_values_reprod_a, fitted_line_reprod_a, coefficient_reprod_a = prepare_linear_regression(df_a_r_large, df_a_r_medium, df_a_r_small, df_a_large, df_a_medium, df_a_small)
	
	r_squared_reprod_a = round.(r_squared_reprod_a; digits=3)
	rms_values_reprod_a = round.(rms_values_reprod_a; digits=3)
	coefficient_reprod_a = round.(coefficient_reprod_a; digits=3)
end

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
	create_textbox(f[1, 1], coefficient_reprod_i, r_squared_reprod_i, rms_values_reprod_i)

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
	create_textbox(f[2, 1], coefficient_reprod_vf, r_squared_reprod_vf, rms_values_reprod_vf)

    xlims!(ax2; low=0, high=200)
    ylims!(ax2; low=0, high=200)
    ax2.xticks = [0, 50, 100, 150, 200]
    ax2.yticks = [0, 50, 100, 150, 200]
    ax2.xlabel = "Mass 1 (mg)"
    ax2.ylabel = "Mass 2 (mg)"
    ax2.title = "Volume Fraction"

    ##-- C --##
    ax3 = Axis(f[3, 1])
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
	create_textbox(f[3, 1], coefficient_reprod_a, r_squared_reprod_a, rms_values_reprod_a)

    xlims!(ax3; low=0, high=200)
    ylims!(ax3; low=0, high=200)
    ax3.xticks = [0, 50, 100, 150, 200]
    ax3.yticks = [0, 50, 100, 150, 200]
    ax3.xlabel = "Mass 1 (mg)"
    ax3.ylabel = "Mass 2 (mg)"
    ax3.title = "Agatston"

    ##-- LABELS --##
    f[2, 2] = Legend(f, ax3; framevisible=false)
    for (label, layout) in zip(["A", "B", "C"], [f[1, 1], f[2, 1], f[3, 1]])
        Label(
            layout[1, 1, TopLeft()],
            label;
            fontsize=25,
            padding=(0, 90, 25, 0),
            halign=:right,
        )
    end

    save(joinpath(dirname(pwd()),"figures", FIGURE_PATH, "reproducibility.png"), f)
    return f
end

# ╔═╡ 86b5f204-d004-47dd-9709-28965692bc4f
with_theme(reprod, medphys_theme)

# ╔═╡ 1bf8321e-b661-45b2-b191-cc83d302422a
begin
	r_squared_reprod_a_norm, rms_values_reprod_a_norm, fitted_line_reprod_a_norm, coefficient_reprod_a_norm = prepare_linear_regression(df_a_r_large_norm, df_a_r_medium_norm, df_a_r_small_norm, df_a_large_norm, df_a_medium_norm, df_a_small_norm)
	
	r_squared_reprod_a_norm = round.(r_squared_reprod_a_norm; digits=3)
	rms_values_reprod_a_norm = round.(rms_values_reprod_a_norm; digits=3)
	coefficient_reprod_a_norm = round.(coefficient_reprod_a_norm; digits=3)
end

# ╔═╡ ab0eb1de-75ea-42d5-8cbc-175286da9777
rms_values_reprod_a_norm

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
begin
	r_squared_reprod_s, rms_values_reprod_s, fitted_line_reprod_s, coefficient_reprod_s = prepare_linear_regression(df_s_r_large, df_s_r_medium, df_s_r_small, df_s_large, df_s_medium, df_s_small)
	
	r_squared_reprod_s = round.(r_squared_reprod_s; digits=3)
	rms_values_reprod_s = round.(rms_values_reprod_s; digits=3)
	coefficient_reprod_s = round.(coefficient_reprod_s; digits=3)
end

# ╔═╡ de5afff2-5383-4b32-8c14-116f8ca548fb
md"""
## Sensitivity and Specificity
"""

# ╔═╡ c3106ba6-5cbd-4f24-b2fa-49b0882396c0
std_level = 1

# ╔═╡ 97985b6c-da75-4c0d-b4a6-2b3402fa7d11
md"""
#### False Negative
"""

# ╔═╡ f4be0691-2f1d-4130-857b-bd8f0948d4fe
md"""
#### SWCS
"""

# ╔═╡ 1b086ed5-fdde-4ea4-8994-ab92cbe871d2
array_s = Array(hcat(df_s[!, :calculated_swcs_large], df_s[!, :calculated_swcs_medium], df_s[!, :calculated_swcs_small]));

# ╔═╡ 81d34483-1bf0-47e9-a929-a8e414147613
total_cac = length(array_s)

# ╔═╡ 99fa2181-393e-4f74-8f76-2dc5df7cc328
mean_s, std_s = mean(df_s[!, :swcs_bkg]), std(df_s[!, :swcs_bkg]) * std_level

# ╔═╡ 88d035df-ffd9-411d-87a0-cbb512dc3f34
total_zero_s = length(findall(x -> x <= mean_s + std_s, array_s))

# ╔═╡ 8a665678-d2a8-403f-9823-cc90c32f3760
md"""
#### Agatston
"""

# ╔═╡ 187ea5d8-a86d-44ed-ae3b-5d284c71f04d
array_a = hcat(df_a[!, :calculated_mass_large], df_a[!, :calculated_mass_medium], df_a[!, :calculated_mass_small]);

# ╔═╡ 2b32141d-6ff4-423a-85e7-d9219311cbb9
num_zero_a = length(findall(x -> x <= 0, array_a))

# ╔═╡ f5dbe1d6-ae15-4921-bce7-d2fa53765558
md"""
#### Integrated
"""

# ╔═╡ 28733b96-971e-4f77-a019-139a68f9ced8
array_i = hcat(df_i[!, :calculated_mass_large], df_i[!, :calculated_mass_medium], df_i[!, :calculated_mass_small]);

# ╔═╡ 0e82f73d-a46f-4fa0-815f-0c221bfec789
mean_i, std_i = mean(df_i[!, :mass_bkg]), std(df_i[!, :mass_bkg]) * std_level

# ╔═╡ 6f684d6b-91ce-40d8-9b9c-5da4d18d9925
total_zero_i = length(findall(x -> x <= mean_i + std_i, array_i))

# ╔═╡ 32c62fd9-01fa-4b60-a0bb-addbf37f5bb4
md"""
#### Volume Fraction
"""

# ╔═╡ 0e5f0094-c24a-4bb1-ba22-8d3c1da77903
array_vf = hcat(df_vf[!, :calculated_mass_large], df_vf[!, :calculated_mass_medium], df_vf[!, :calculated_mass_small]);

# ╔═╡ c6c9e7d0-9846-4303-b8c5-aa262a428623
mean_vf, std_vf = mean(df_vf[!, :mass_bkg]), std(df_vf[!, :mass_bkg]) * std_level

# ╔═╡ 6d842833-aa03-4e72-9e57-8bfa6e0a0b6e
total_zero_vf = length(findall(x -> x <= mean_i + std_i, array_i))

# ╔═╡ f4ff80f2-9146-48c6-90eb-fd5cbbcbb1b7
md"""
#### False Positive
"""

# ╔═╡ 8f45145a-dc1c-4291-947d-9e6eb8349b6d
md"""
#### SWCS
"""

# ╔═╡ 278a4320-673a-45ba-9de6-7f525aab406f
array_s_pos = df_s[!, :swcs_bkg]

# ╔═╡ c3c30ec3-052d-46bb-b6d7-a73d162c9f46
total_cac_pos = length(array_s_pos)

# ╔═╡ 3c4aff49-c42b-4f5d-a6a7-af4f7fa03197
total_zero_s_pos = length(findall(x -> x >= mean_s + std_s, array_s_pos))

# ╔═╡ 7b17e354-11b0-473e-b0b5-0999af7dc022
md"""
#### Agatston
"""

# ╔═╡ 98e91087-bbb1-4215-86db-9c6f6ea7c5b2
array_a_pos = df_a[!, :mass_bkg]

# ╔═╡ 0e02b18f-21a8-4dd2-b8c9-5d9a553d219d
total_zero_a_pos = length(findall(x -> x > 0, array_a_pos))

# ╔═╡ da61c2d8-7ffe-4c13-a955-d9fe9001997f
md"""
#### Integrated
"""

# ╔═╡ a5f6027f-16ca-4788-ba52-a885d8dce519
array_i_pos = df_i[!, :mass_bkg]

# ╔═╡ 46cef33e-5ef3-4689-8461-3f227b3cc584
total_zero_i_pos = length(findall(x -> x >= mean_i + std_i, array_i_pos))

# ╔═╡ 773f86e3-7d02-4f20-a9d6-b4bad563b2c4
md"""
#### Volume Fraction
"""

# ╔═╡ cb3c6875-9886-4134-aa5b-a96303bb75c2
array_vf_pos = df_vf[!, :mass_bkg]

# ╔═╡ 62d9b719-6e29-4def-a404-26645e132180
total_zero_vf_pos = length(findall(x -> x >= mean_vf + std_vf, array_vf_pos))

# ╔═╡ 6590c12f-2345-488f-a3c8-ffa289b9438c
function sensitivity_specificity()
    f = Figure()
    colors = Makie.wong_colors()

    ##-- TOP --##
    ax1 = Axis(f[1, 1]; xticks=(1:3, ["Integrated", "Volume Fraction", "Agatston"]))

    table = [1, 2, 3]
    heights1 = [
		(total_zero_i / total_cac) * 100,
		(total_zero_vf / total_cac) * 100,
        (num_zero_a / total_cac) * 100,
    ]
    barplot!(ax1, table, heights1; color=colors[1:3], bar_labels=:y)

    ax1.title = "False-Negative Scores (CAC=0)"
    ax1.ylabel = "% False-Negative Zero CAC Scores"
    ylims!(ax1; low=0, high=100)
    ax1.yticks = [0, 25, 50, 75, 100]

	##-- TOP --##
    ax2 = Axis(f[2, 1]; xticks=(1:3, ["Integrated", "Volume Fraction", "Agatston"]))

    table = [1, 2, 3]
    heights1 = [
		(total_zero_i_pos / total_cac_pos) * 100,
		(total_zero_vf_pos / total_cac_pos) * 100,
        (total_zero_a_pos / total_cac_pos) * 100,
    ]
    barplot!(ax2, table, heights1; color=colors[1:3], bar_labels=:y)

    ax2.title = "False-Positive Scores (CAC>0)"
    ax2.ylabel = "% False-Positive CAC Scores"
    ylims!(ax2; low=0, high=100)
    ax2.yticks = [0, 25, 50, 75, 100]

    save(joinpath(dirname(pwd()),"figures", FIGURE_PATH, "sensitivity_specificity.png"), f, resolution=(800, 800))

	
    return f
end

# ╔═╡ 0cbe1116-c343-4f93-abac-4e44a8f9c3f9
with_theme(medphys_theme) do
    sensitivity_specificity()
end

# ╔═╡ 4d8d3832-3b61-4f88-b4fa-0c0861d7f5a1
md"""
# Tables
"""

# ╔═╡ 58729782-e6df-4e98-8a4e-9cfdcba81e72
md"""
## Accuracy
"""

# ╔═╡ e4e8f32a-e6d8-43a2-bdfc-8fe51c6bd82d
integrated_normal_accuracy = [r_squared_normal_i, rms_values_normal_i..., "y = $(coefficient_normal_i[2])x + $(coefficient_normal_i[1])"]

# ╔═╡ 3dcc90ca-bbf4-4228-ac2f-eeb2bf31c12d
volume_fraction_normal_accuracy = [r_squared_normal_vf, rms_values_normal_vf..., "y = $(coefficient_normal_vf[2])x + $(coefficient_normal_vf[1])"]

# ╔═╡ 0e0c5392-9971-4eea-aac0-fc2d6f0e897b
agatston_normal_accuracy = [r_squared_normal_a, rms_values_normal_a..., "y = $(coefficient_normal_a[2])x + $(coefficient_normal_a[1])"]

# ╔═╡ b08e7070-28fe-43c9-bfb5-6812feda20d7
md"""
## Reproducibility
"""

# ╔═╡ 4c57e3d6-daf3-4dcb-9330-2dff0074b19d
integrated_reprod_accuracy = [r_squared_reprod_i, rms_values_reprod_i..., "y = $(coefficient_reprod_i[2])x + $(coefficient_reprod_i[1])"]

# ╔═╡ 127ef411-2ef8-41ff-8a05-93c9240bd7b1
volume_fraction_reprod_accuracy = [r_squared_reprod_vf, rms_values_reprod_vf..., "y = $(coefficient_reprod_vf[2])x + $(coefficient_reprod_vf[1])"]

# ╔═╡ 98d6d8ed-e3f0-47b6-8183-3325053e6ca2
agatston_reprod_accuracy = [r_squared_reprod_a, rms_values_reprod_a..., "y = $(coefficient_reprod_a[2])x + $(coefficient_reprod_a[1])"]

# ╔═╡ af6b8453-c2f7-4dec-8dca-911e08bd0f2b
swcs_reprod_accuracy = [r_squared_reprod_s, rms_values_reprod_s..., "y = $(coefficient_reprod_s[2])x + $(coefficient_reprod_s[1])"]

# ╔═╡ 422aa020-585f-45be-ab30-4870689a5a68
md"""
## Sensitivity and Specificity
"""

# ╔═╡ fd91903c-ede5-43d3-8772-dca309503cdf
begin
	integrated_sens_spec = ["$total_zero_i / $total_cac", "$total_zero_i_pos / $total_cac_pos"]
	volume_fraction_sens_spec = ["$total_zero_vf / $total_cac", "$total_zero_vf_pos / $total_cac_pos"]
	agatston_sens_spec = ["$num_zero_a / $total_cac", "$total_zero_a_pos / $total_cac_pos"]
end

# ╔═╡ 20a25474-5bf5-4267-83ff-32dde8b54371
md"""
# Summaries
"""

# ╔═╡ da4a557d-550b-4450-8022-6a70d3c9105a
begin
	df_accuracy = DataFrame(
		"Metrics" => ["r^2", "RMSE", "RMSD", "Best Fit Line"],
		"Integrated" => integrated_normal_accuracy,
		"Volume Fraction" => volume_fraction_normal_accuracy,
		"Agatston" => agatston_normal_accuracy,
	)
	df_reprod = DataFrame(
		"Metrics" => ["r^2", "RMSE", "RMSD", "Best Fit Line"],
		"Integrated" => integrated_reprod_accuracy,
		"Volume Fraction" => volume_fraction_reprod_accuracy,
		"Agatston" => agatston_reprod_accuracy,
	)
	df_sens_spec = DataFrame(
		"Metrics" => ["False Negatives", "False Positives"],
		"Integrated" => integrated_sens_spec,
		"Volume Fraction" => volume_fraction_sens_spec,
		"Agatston" => agatston_sens_spec
	)
	_dfs = append!(df_accuracy, df_reprod, cols=:union)
	dfs = append!(_dfs, df_sens_spec, cols=:union)
	
	results = [repeat(["Accuracy (Normal Density)"], 4)..., repeat(["Reproducibility"], 4)..., "Specificity", "Sensitivity"]
	insertcols!(dfs, 1, :Results=>results)
end;

# ╔═╡ 6d466a96-601e-40bb-aeff-41a867048118
dfs

# ╔═╡ e7e38213-3f4d-4563-8875-bf8e52bcd129
CSV.write(joinpath(dirname(pwd()),"figures", FIGURE_PATH, "summary.csv"), dfs)

# ╔═╡ 374a8805-e0c9-47d0-a14b-b7d585013827
md"""
## Simulation Parameters
"""

# ╔═╡ 5e8c5d10-6b80-4093-94f1-101d80b73615
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

# ╔═╡ 40783122-036b-4fa8-8ccf-b7ea741152ae
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

# ╔═╡ 28e88ef6-dc45-4297-8348-79ac44daaf17
simulation_parameters = DataFrame(
	"Parameter" => parameters,
	"Simulation" => simulations
)

# ╔═╡ 6182ddf5-97bc-47b8-9735-40d1c5e07144
md"""
## Spatially Weighted Calcium Scoring Means & Stds
"""

# ╔═╡ b6e8232e-a37d-44ad-9712-6a7e6ec20808
kvp = [
	80
	100
	120
	135
]

# ╔═╡ 42ef44d2-f2e6-4030-971b-bb6187eec577
means = [
	201.4898
	174.3658
	158.2645
	152.8815
]

# ╔═╡ 266447ef-e127-475b-9bcb-549df9d974ed
stds = [
	34.3038
	28.1708
	22.9656
	21.0778
]

# ╔═╡ 4fdc8ce0-2ae2-40bc-9eeb-419277787991
means_stds = DataFrame(
	"kVp" => kvp,
	"Mean" => means,
	"Std" => stds
)

# ╔═╡ Cell order:
# ╠═82a02e07-06d4-4e15-924d-c515b76d17f8
# ╠═90bd3c98-d3c1-4209-8692-256b91679a62
# ╟─fc91cc88-b137-49aa-a101-1049b8721a94
# ╟─d68438a7-5b1c-4c49-8fcf-9e60b68f0aee
# ╟─522a9c57-d9b4-4c03-889e-d2a78a21ca5d
# ╟─b302881f-99cd-4ac6-b22f-6d6d5c098e27
# ╟─e912ad90-9477-4253-97db-f679d2c9abb3
# ╟─2e862b70-3c59-4572-8e71-f3233b4363fa
# ╠═5d6c0929-6608-4a01-877b-b70a3a8dbdfe
# ╟─2c8d827a-10bc-48a9-8657-c7c8d83fd5df
# ╟─312b7073-caa1-44eb-9c53-7d2508a24603
# ╠═bdecc190-cecc-4f4a-b865-742ef553ca6a
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
# ╟─9248087e-2e21-4bd5-aed7-c5eedbcdf0c1
# ╠═e7e38a45-aed8-4e25-8fdc-d08cb1b3f197
# ╠═d4ddfb59-5f1b-414c-87aa-aed494e008d9
# ╠═fa3a6b92-d83a-443e-81f1-c6bd48f55d60
# ╠═1e50121c-ac1f-4107-863b-d2da79c77700
# ╟─f011edbc-2547-4be8-895b-5ebe611be122
# ╟─f534bb1e-3c46-476a-ae05-1384fcab1260
# ╟─86b5f204-d004-47dd-9709-28965692bc4f
# ╟─3e5c3d00-fcd3-4fc6-9a52-4f728a5ba69c
# ╠═3ce29e9d-3be6-41a2-9cff-2399aff77b0a
# ╠═0f58b7de-28e1-4704-9a31-f56d0c7c1b18
# ╠═7e65b2f7-6e24-4f2c-8c2d-3c6eebc41c75
# ╠═b7749170-0017-4e65-8bb6-2d95da3c783d
# ╠═0afaa9ac-0d8f-4474-ad90-aae299a873f0
# ╠═a1b4c6e2-662a-4b0f-bda1-b154d4630aa8
# ╟─dfcbccc0-3a6f-491c-8e6f-f68a96bf8a7c
# ╠═c49cb0e1-700e-4c77-bcbb-1b815e78c902
# ╠═c7719179-5c77-4e6b-ad1d-86ea90ff78f9
# ╠═54cde88d-515b-4304-a382-165542a396c3
# ╠═2974a158-a1c7-4509-b609-60711195df12
# ╠═99e53995-61a3-45f3-a699-8bdef0dfe152
# ╠═eaae39f7-21d8-488d-841f-e8a4e57175bd
# ╟─3daa60e3-b1ba-4106-9256-3e7919c926d0
# ╠═37581950-10e8-4be0-95b0-a6cd53e2f556
# ╠═a144cec4-0c90-4713-abbc-23795b3d8584
# ╠═a087b301-3b51-44b0-8afd-74f983c26b06
# ╠═c482a5e0-aeb4-4300-8dd1-79e99ac688ee
# ╠═1bf8321e-b661-45b2-b191-cc83d302422a
# ╠═ab0eb1de-75ea-42d5-8cbc-175286da9777
# ╟─666b7001-c521-4bb2-ac40-20966fb9ada5
# ╠═7e250de1-bb94-4a8a-8d1c-bf0479cc40a5
# ╠═aac1e177-5f76-496d-b6a1-254c4338f25c
# ╠═7c8dfdba-7b75-4812-b5f7-b0ed7ed08b31
# ╠═6103e58a-f454-4894-aea9-afa706ad11a6
# ╠═8c28ca4e-eff1-4921-94e1-2a04f3a06a1a
# ╟─de5afff2-5383-4b32-8c14-116f8ca548fb
# ╟─6590c12f-2345-488f-a3c8-ffa289b9438c
# ╠═0cbe1116-c343-4f93-abac-4e44a8f9c3f9
# ╠═c3106ba6-5cbd-4f24-b2fa-49b0882396c0
# ╟─97985b6c-da75-4c0d-b4a6-2b3402fa7d11
# ╟─f4be0691-2f1d-4130-857b-bd8f0948d4fe
# ╠═1b086ed5-fdde-4ea4-8994-ab92cbe871d2
# ╠═81d34483-1bf0-47e9-a929-a8e414147613
# ╠═99fa2181-393e-4f74-8f76-2dc5df7cc328
# ╠═88d035df-ffd9-411d-87a0-cbb512dc3f34
# ╟─8a665678-d2a8-403f-9823-cc90c32f3760
# ╠═187ea5d8-a86d-44ed-ae3b-5d284c71f04d
# ╠═2b32141d-6ff4-423a-85e7-d9219311cbb9
# ╟─f5dbe1d6-ae15-4921-bce7-d2fa53765558
# ╠═28733b96-971e-4f77-a019-139a68f9ced8
# ╠═0e82f73d-a46f-4fa0-815f-0c221bfec789
# ╠═6f684d6b-91ce-40d8-9b9c-5da4d18d9925
# ╟─32c62fd9-01fa-4b60-a0bb-addbf37f5bb4
# ╠═0e5f0094-c24a-4bb1-ba22-8d3c1da77903
# ╠═c6c9e7d0-9846-4303-b8c5-aa262a428623
# ╠═6d842833-aa03-4e72-9e57-8bfa6e0a0b6e
# ╟─f4ff80f2-9146-48c6-90eb-fd5cbbcbb1b7
# ╟─8f45145a-dc1c-4291-947d-9e6eb8349b6d
# ╠═278a4320-673a-45ba-9de6-7f525aab406f
# ╠═c3c30ec3-052d-46bb-b6d7-a73d162c9f46
# ╠═3c4aff49-c42b-4f5d-a6a7-af4f7fa03197
# ╟─7b17e354-11b0-473e-b0b5-0999af7dc022
# ╠═98e91087-bbb1-4215-86db-9c6f6ea7c5b2
# ╠═0e02b18f-21a8-4dd2-b8c9-5d9a553d219d
# ╟─da61c2d8-7ffe-4c13-a955-d9fe9001997f
# ╠═a5f6027f-16ca-4788-ba52-a885d8dce519
# ╠═46cef33e-5ef3-4689-8461-3f227b3cc584
# ╟─773f86e3-7d02-4f20-a9d6-b4bad563b2c4
# ╠═cb3c6875-9886-4134-aa5b-a96303bb75c2
# ╠═62d9b719-6e29-4def-a404-26645e132180
# ╟─4d8d3832-3b61-4f88-b4fa-0c0861d7f5a1
# ╟─58729782-e6df-4e98-8a4e-9cfdcba81e72
# ╠═e4e8f32a-e6d8-43a2-bdfc-8fe51c6bd82d
# ╠═3dcc90ca-bbf4-4228-ac2f-eeb2bf31c12d
# ╠═0e0c5392-9971-4eea-aac0-fc2d6f0e897b
# ╟─b08e7070-28fe-43c9-bfb5-6812feda20d7
# ╠═4c57e3d6-daf3-4dcb-9330-2dff0074b19d
# ╠═127ef411-2ef8-41ff-8a05-93c9240bd7b1
# ╠═98d6d8ed-e3f0-47b6-8183-3325053e6ca2
# ╠═af6b8453-c2f7-4dec-8dca-911e08bd0f2b
# ╟─422aa020-585f-45be-ab30-4870689a5a68
# ╠═fd91903c-ede5-43d3-8772-dca309503cdf
# ╟─20a25474-5bf5-4267-83ff-32dde8b54371
# ╠═da4a557d-550b-4450-8022-6a70d3c9105a
# ╠═6d466a96-601e-40bb-aeff-41a867048118
# ╠═e7e38213-3f4d-4563-8875-bf8e52bcd129
# ╟─374a8805-e0c9-47d0-a14b-b7d585013827
# ╟─5e8c5d10-6b80-4093-94f1-101d80b73615
# ╟─40783122-036b-4fa8-8ccf-b7ea741152ae
# ╠═28e88ef6-dc45-4297-8348-79ac44daaf17
# ╟─6182ddf5-97bc-47b8-9735-40d1c5e07144
# ╠═b6e8232e-a37d-44ad-9712-6a7e6ec20808
# ╠═42ef44d2-f2e6-4030-971b-bb6187eec577
# ╠═266447ef-e127-475b-9bcb-549df9d974ed
# ╠═4fdc8ce0-2ae2-40bc-9eeb-419277787991

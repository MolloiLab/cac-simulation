### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ c82bb3e6-e052-11ec-3392-eba4a4b464ba
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.activate(".")

    using PlutoUI, Statistics, CSV, DataFrames, GLM, CairoMakie, HypothesisTests, Colors, MLJBase, Printf
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

# ╔═╡ c5567a17-2345-470f-97c8-9383617298ac
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

# ╔═╡ 043dec5f-07a8-472d-8082-c7a771f93270
FIGURE_PATH = "stationary"

# ╔═╡ e79c34d0-bf06-4a7a-93b7-e49f09ce3f4e
md"""
# Load CSVs
"""

# ╔═╡ 4cf61487-3913-4e11-970e-4ca31a5ffc8d
md"""
## Integrated
"""

# ╔═╡ 7e1b5b89-ed32-42df-aad9-f38524d36988
root_dir = joinpath(dirname(pwd()), "output_new")

# ╔═╡ 57d46998-368a-4916-8380-ee49d5473a49
path_integrated = joinpath(root_dir, "integrated")

# ╔═╡ 2c960bd8-ae64-453f-b29f-275bf5263774
df_i = CSV.read(string(path_integrated, "/full2.csv"), DataFrame);

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
path_agat = joinpath(root_dir, "agatston")

# ╔═╡ 3fa3748a-22d6-49d8-9888-a749dca99ce2
df_a = CSV.read(string(path_agat, "/full2.csv"), DataFrame);

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
path_swcs = joinpath(root_dir, "swcs")

# ╔═╡ 317f8c26-7b39-45c3-b9d3-02925d4b0514
df_s = CSV.read(string(path_swcs, "/full2.csv"), DataFrame);

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
path_vf = joinpath(root_dir, "volume_fraction")

# ╔═╡ 5e803aad-9d82-4801-aff6-3c1c2f748b02
df_vf = CSV.read(string(path_vf, "/full2.csv"), DataFrame);

# ╔═╡ 5c76c3c7-bb56-4d90-8f0f-7f7d0e271d3c
df_vf_low, df_vf_normal = groupby(df_vf, :DENSITY);

# ╔═╡ 0085693a-157d-42a8-9d35-ec0feec24c9c
df_vf_low_small, df_vf_low_medium, df_vf_low_large = groupby(df_vf_low, :SIZE);

# ╔═╡ a03343e0-a9f4-4a49-8490-1d542bbe42fd
df_vf_normal_small, df_vf_normal_medium, df_vf_normal_large = groupby(df_vf_normal, :SIZE);

# ╔═╡ 593d7b00-94c4-49c2-87a7-723dbef64667
md"""
# Figures
"""

# ╔═╡ 63f97b37-dc0e-456d-811f-409f8360b2da
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

# ╔═╡ ac937ee5-e820-4960-836f-50df8c7ae3ee
md"""
## Accuracy
"""

# ╔═╡ 9b3c6bda-5c21-4c07-be6f-939a1ac0496a
md"""
#### Normal Density
"""

# ╔═╡ b4e25c03-ac59-4104-acf5-2f049497d460
begin
	r_squared_normal_i, rms_values_normal_i, fitted_line_normal_i, coefficient_normal_i = prepare_linear_regression(df_i_normal)
	r_squared_normal_i = round.(r_squared_normal_i; digits=3)
	rms_values_normal_i = round.(rms_values_normal_i; digits=3)
	coefficient_normal_i = round.(coefficient_normal_i; digits=3)
end

# ╔═╡ 074a6437-9745-4680-b9b5-baf9d767afc9
begin
	r_squared_normal_vf, rms_values_normal_vf, fitted_line_normal_vf, coefficient_normal_vf = prepare_linear_regression(df_vf_normal)
	r_squared_normal_vf = round.(r_squared_normal_vf; digits=3)
	rms_values_normal_vf = round.(rms_values_normal_vf; digits=3)
	coefficient_normal_vf = round.(coefficient_normal_vf; digits=3)
end

# ╔═╡ f51b5a1e-a688-470d-8a83-b65f03fff875
begin
	r_squared_normal_a, rms_values_normal_a, fitted_line_normal_a, coefficient_normal_a = prepare_linear_regression(df_a_normal)
	r_squared_normal_a = round.(r_squared_normal_a; digits=3)
	rms_values_normal_a = round.(rms_values_normal_a; digits=3)
	coefficient_normal_a = round.(coefficient_normal_a; digits=3)
end

# ╔═╡ decd4320-112a-4c4c-9ed5-680a721b81c1
function lin_reg_norm()
    f = Figure()

    ##-- A --##
    ax1 = Axis(f[1, 1])
    df = df_i_normal
    scatter!(
		df[!, :ground_truth_mass_large], df[!, :calculated_mass_large]
	)
    errorbars!(
		df[!, :ground_truth_mass_large], df[!, :calculated_mass_large], rms(df[!, :ground_truth_mass_large], df[!, :calculated_mass_large]),
    )
    scatter!(
		df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium]
	)
    errorbars!(
		df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium], rms(df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium]),
    )
    scatter!(
		df[!, :ground_truth_mass_small], df[!, :calculated_mass_small]; color=:red
	)
    errorbars!(
		df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], rms(df[!, :ground_truth_mass_small], df[!, :calculated_mass_small]),
    )
    lines!(ax1, [-1000, 1000], [-1000, 1000])
    lines!(ax1, collect(1:1000), fitted_line_normal_i; linestyle=:dashdot)
	create_textbox(f[1, 1], coefficient_normal_i, r_squared_normal_i, rms_values_normal_i)

    xlims!(ax1; low=0, high=200)
    ylims!(ax1; low=0, high=200)
    ax1.xticks = [0, 50, 100, 150, 200]
    ax1.yticks = [0, 100, 200]
    ax1.xlabel = "Known Mass (mg)"
    ax1.ylabel = "Calculated Mass (mg)"
    ax1.title = "Integrated (Normal-Density)"
	
    ##-- B --##
    ax2 = Axis(f[2, 1])
    df3 = df_vf_normal
    sc1 = scatter!(
		df3[!, :ground_truth_mass_large], df3[!, :calculated_mass_large]
	)
    errorbars!(
        df3[!, :ground_truth_mass_large],
        df3[!, :calculated_mass_large],
        rms(df3[!, :ground_truth_mass_large], df3[!, :calculated_mass_large]),
    )
    sc2 = scatter!(df3[!, :ground_truth_mass_medium], df3[!, :calculated_mass_medium]
	)
    errorbars!(
        df3[!, :ground_truth_mass_medium],
        df3[!, :calculated_mass_medium],
        rms(df3[!, :ground_truth_mass_medium], df3[!, :calculated_mass_medium]),
    )
    sc3 = scatter!(
		df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small]; color=:red
    )
    errorbars!(
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
    ax2.yticks = [0, 100, 200]
    ax2.xlabel = "Known Mass (mg)"
    ax2.ylabel = "Calculated Mass (mg)"
    ax2.title = "Volume Fraction (Normal-Density)"


    ##-- C --##
    ax2 = Axis(f[3, 1])
    df3 = df_a_normal
    sc1 = scatter!(
		df3[!, :ground_truth_mass_large], df3[!, :calculated_mass_large]
	)
    errorbars!(
        df3[!, :ground_truth_mass_large],
        df3[!, :calculated_mass_large],
        rms(df3[!, :ground_truth_mass_large], df3[!, :calculated_mass_large]),
    )
    sc2 = scatter!(
		df3[!, :ground_truth_mass_medium], df3[!, :calculated_mass_medium]
	)
    errorbars!(
        df3[!, :ground_truth_mass_medium],
        df3[!, :calculated_mass_medium],
        rms(df3[!, :ground_truth_mass_medium], df3[!, :calculated_mass_medium]),
    )
    sc3 = scatter!(
		df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small]; color=:red
    )
    errorbars!(
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
    ax2.yticks = [0, 100, 200]
    ax2.xlabel = "Known Mass (mg)"
    ax2.ylabel = "Calculated Mass (mg)"
    ax2.title = "Agatston (Normal-Density)"

    ##-- LABELS --##
    f[2, 2] = Legend(
        f,
        [sc1, sc2, sc3, ln1, ln2],
        ["Large Inserts", "Medium Inserts", "Small Inserts", "Unity", "Fitted Line"];
        framevisible=false,
		padding=(10, 10, 8, 8)
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

	# Label(
	# 	f[1:end, 0],
	# 	"Calculated Mass (mg)";
	# 	fontsize=20,
	# 	rotation = pi/2
	# )

    save(joinpath(dirname(pwd()),"figures", FIGURE_PATH, "accuracy_normal.png"), f)
    return f
end

# ╔═╡ 0fa91a0b-c205-4b10-9876-9a9e163f3d7e
with_theme(lin_reg_norm, medphys_theme)

# ╔═╡ 91a157fb-a515-4ab0-be1a-8d9148db6d67
md"""
#### Low Density
"""

# ╔═╡ 1478eef5-bcde-4513-8794-9b9cbe388ee2
begin
	r_squared_low_i, rms_values_low_i, fitted_line_low_i, coefficient_low_i = prepare_linear_regression(df_i_low)
	r_squared_low_i = round.(r_squared_low_i; digits=3)
	rms_values_low_i = round.(rms_values_low_i; digits=3)
	coefficient_low_i = round.(coefficient_low_i; digits=3)
end

# ╔═╡ ffeb8498-d7a3-4122-b71f-56326490042b
begin
	r_squared_low_vf, rms_values_low_vf, fitted_line_low_vf, coefficient_low_vf = prepare_linear_regression(df_vf_low)
	r_squared_low_vf = round.(r_squared_low_vf; digits=3)
	rms_values_low_vf = round.(rms_values_low_vf; digits=3)
	coefficient_low_vf = round.(coefficient_low_vf; digits=3)
end

# ╔═╡ ad784be1-5252-4914-bdd2-085edf5a95e0
begin
	r_squared_low_a, rms_values_low_a, fitted_line_low_a, coefficient_low_a = prepare_linear_regression(df_a_low)
	r_squared_low_a = round.(r_squared_low_a; digits=3)
	rms_values_low_a = round.(rms_values_low_a; digits=3)
	coefficient_low_a = round.(coefficient_low_a; digits=3)
end

# ╔═╡ bc9939d7-f3a8-48e5-ada4-8fcf0fc1137a
function lin_reg_low()
    f = Figure()
    ##-- A --##
    ax = Axis(
		f[1, 1],
		xticks = [0, 5, 10, 15, 20, 25],
		yticks = [0, 15, 30],
		xlabel = "Known Mass (mg)",
    	ylabel = "Calculated Mass (mg)",
    	title = "Integrated (Low-Density)"
	)

    df2 = df_i_low
    sc1 = scatter!(
		df2[!, :ground_truth_mass_large], df2[!, :calculated_mass_large]
	)
    errorbars!(
        df2[!, :ground_truth_mass_large],
        df2[!, :calculated_mass_large],
        rms(df2[!, :ground_truth_mass_large], df2[!, :calculated_mass_large]),
    )
    sc2 = scatter!(
		df2[!, :ground_truth_mass_medium], df2[!, :calculated_mass_medium]
	)
    errorbars!(
        df2[!, :ground_truth_mass_medium],
        df2[!, :calculated_mass_medium],
        rms(df2[!, :ground_truth_mass_medium], df2[!, :calculated_mass_medium]),
    )
    sc3 = scatter!(
		df2[!, :ground_truth_mass_small], df2[!, :calculated_mass_small]; color=:red
    )
    errorbars!(
        df2[!, :ground_truth_mass_small],
        df2[!, :calculated_mass_small],
        rms(df2[!, :ground_truth_mass_small], df2[!, :calculated_mass_small]),
    )
    ln1 = lines!([-1000, 1000], [-1000, 1000])
    ln2 = lines!(collect(1:1000), fitted_line_low_i; linestyle=:dashdot)
	create_textbox(f[1, 1], coefficient_low_i, r_squared_low_i, rms_values_low_i)

    xlims!(ax; low=0, high=25)
    ylims!(ax; low=-5, high=30)

	##-- B --##
    ax = Axis(
		f[2, 1],
		xticks = [0, 5, 10, 15, 20, 25],
		yticks = [0, 15, 30],
		xlabel = "Known Mass (mg)",
    	ylabel = "Calculated Mass (mg)",
    	title = "Volume Fraction (Low-Density)"
	)


    df4 = df_vf_low
    sc1 = scatter!(
		df4[!, :ground_truth_mass_large], df4[!, :calculated_mass_large]
	)
    errorbars!(
        df4[!, :ground_truth_mass_large],
        df4[!, :calculated_mass_large],
        rms(df4[!, :ground_truth_mass_large], df4[!, :calculated_mass_large]),
    )
    sc2 = scatter!(df4[!, :ground_truth_mass_medium], df4[!, :calculated_mass_medium]
	)
    errorbars!(
        df4[!, :ground_truth_mass_medium],
        df4[!, :calculated_mass_medium],
        rms(df4[!, :ground_truth_mass_medium], df4[!, :calculated_mass_medium]),
    )
    sc3 = scatter!(
		df4[!, :ground_truth_mass_small], df4[!, :calculated_mass_small]; color=:red
    )
    errorbars!(
        df4[!, :ground_truth_mass_small],
        df4[!, :calculated_mass_small],
        rms(df4[!, :ground_truth_mass_small], df4[!, :calculated_mass_small]),
    )
    ln1 = lines!([-1000, 1000], [-1000, 1000])
    ln2 = lines!(collect(1:1000), fitted_line_low_vf; linestyle=:dashdot)
	create_textbox(f[2, 1], coefficient_low_vf, r_squared_low_vf, rms_values_low_vf)
	
    xlims!(ax; low=0, high=25)
    ylims!(ax; low=-5, high=30)

    ##-- C --##
    ax = Axis(
		f[3, 1],
		xticks = [0, 5, 10, 15, 20, 25],
		yticks = [0, 15, 30],
		xlabel = "Known Mass (mg)",
    	ylabel = "Calculated Mass (mg)",
    	title = "Agatston (Low-Density)"
	)

    df4 = df_a_low
    sc1 = scatter!(
		df4[!, :ground_truth_mass_large], df4[!, :calculated_mass_large]
	)
    errorbars!(
        df4[!, :ground_truth_mass_large],
        df4[!, :calculated_mass_large],
        rms(df4[!, :ground_truth_mass_large], df4[!, :calculated_mass_large]),
    )
    sc2 = scatter!(
		df4[!, :ground_truth_mass_medium], df4[!, :calculated_mass_medium]
	)
    errorbars!(
        df4[!, :ground_truth_mass_medium],
        df4[!, :calculated_mass_medium],
        rms(df4[!, :ground_truth_mass_medium], df4[!, :calculated_mass_medium]),
    )
    sc3 = scatter!(
		df4[!, :ground_truth_mass_small], df4[!, :calculated_mass_small]; color=:red
    )
    errorbars!(
        df4[!, :ground_truth_mass_small],
        df4[!, :calculated_mass_small],
        rms(df4[!, :ground_truth_mass_small], df4[!, :calculated_mass_small]),
    )
    ln1 = lines!([-1000, 1000], [-1000, 1000])
    ln2 = lines!(collect(1:1000), fitted_line_low_a; linestyle=:dashdot)
	create_textbox(f[3, 1], coefficient_low_a, r_squared_low_a, rms_values_low_a)

    xlims!(ax; low=0, high=25)
    ylims!(ax; low=-5, high=30)

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

    save(joinpath(dirname(pwd()),"figures", FIGURE_PATH, "accuracy_low.png"), f)
    return f
end

# ╔═╡ 3f02194d-9c69-45d1-80fc-4cc50b1f3d8f
with_theme(lin_reg_low, medphys_theme)

# ╔═╡ d4753644-29fb-41a9-b8d1-5f46bb87a613
md"""
## Reproducibility
"""

# ╔═╡ aef525b3-be1a-4a80-ba18-7ade91baea2d
root_dir_r = joinpath(dirname(pwd()), "output_repeated")

# ╔═╡ f8a69b2e-1cab-4281-9eb5-22bed0046c49
md"""
#### Integrated
"""

# ╔═╡ af6d4416-6750-4674-8d5d-88456d58f841
path_integrated_r = joinpath(root_dir_r, "integrated")

# ╔═╡ 2e531f2f-6d07-4818-b9f8-ecfa8a5d419a
df_i_r = CSV.read(string(path_integrated_r, "/full2.csv"), DataFrame);

# ╔═╡ b7ad23f1-11d4-4146-8eae-37c0856b6e01
df_i_low_r, df_i_normal_r = groupby(df_i_r, :DENSITY);

# ╔═╡ a1a417b0-2f19-4911-bec7-4bac0d2358b7
df_i_low_small_r, df_i_low_medium_r, df_i_low_large_r = groupby(df_i_low_r, :SIZE);

# ╔═╡ c31d03f6-3846-4f78-aa76-97a397af500a
df_i_normal_small_r, df_i_normal_medium_r, df_i_normal_large_r = groupby(
    df_i_normal_r, :SIZE
);

# ╔═╡ 74c1846e-4e69-4434-a695-7b197ed0df0a
array_i_r = hcat(
    df_i_r[!, :calculated_mass_large],
    df_i_r[!, :calculated_mass_medium],
    df_i_r[!, :calculated_mass_small],
)

# ╔═╡ 52261c32-415d-425f-898b-d1bbf66052cf
array_i = hcat(df_i[!, :calculated_mass_large], df_i[!, :calculated_mass_medium], df_i[!, :calculated_mass_small]);

# ╔═╡ f7b5997f-8aba-4311-b120-836e808fee2f
md"""
#### Volume Fraction
"""

# ╔═╡ 47ff39ea-1f55-4e7d-94de-5d4ded126648
path_volume_fraction_r = joinpath(root_dir_r, "volume_fraction")

# ╔═╡ 24e3e4d1-91b4-4ae0-abf8-174118926b6b
df_vf_r = CSV.read(string(path_volume_fraction_r, "/full2.csv"), DataFrame);

# ╔═╡ fc1e3216-6929-4f8f-b1c9-b700f3e35181
df_vf_low_r, df_vf_normal_r = groupby(df_vf_r, :DENSITY);

# ╔═╡ 6ea76780-2726-4d61-bff1-fb590f7ba2f5
df_vf_low_small_r, df_vf_low_medium_r, df_vf_low_large_r = groupby(df_vf_low_r, :SIZE);

# ╔═╡ 288f8d4c-323a-4962-8aff-026ef1bfa77d
df_vf_normal_small_r, df_vf_normal_medium_r, df_vf_normal_large_r = groupby(
    df_vf_normal_r, :SIZE
);

# ╔═╡ fd9e7b94-03c5-4ea9-8405-2218e07df3b5
array_vf_r = hcat(
    df_vf_r[!, :calculated_mass_large],
    df_vf_r[!, :calculated_mass_medium],
    df_vf_r[!, :calculated_mass_small],
)

# ╔═╡ 3be21a82-9524-4edb-93d3-4d515b8b8758
md"""
#### Agatston
"""

# ╔═╡ 04afe0fa-c290-48b3-aa68-c8f74346a321
path_agat_r = joinpath(root_dir_r, "agatston")

# ╔═╡ d181bc5e-468d-42e0-8284-986756a1368a
df_a_r = CSV.read(string(path_agat_r, "/full2.csv"), DataFrame);

# ╔═╡ 814ff1b6-103c-44b9-b9fd-f334507bf998
df_a_low_r, df_a_normal_r = groupby(df_a_r, :DENSITY);

# ╔═╡ 92daece0-acb9-415d-9dbc-d773bffdc18c
df_a_low_small_r, df_a_low_medium_r, df_a_low_large_r = groupby(df_a_low_r, :SIZE);

# ╔═╡ 305a1661-e34d-41ff-912a-386be27aa76d
df_a_normal_small_r, df_a_normal_medium_r, df_a_normal_large_r = groupby(
    df_a_normal_r, :SIZE
);

# ╔═╡ ce630749-7c8e-41c7-853e-422092d64958
array_a_r = hcat(
    df_a_r[!, :calculated_mass_large],
    df_a_r[!, :calculated_mass_medium],
    df_a_r[!, :calculated_mass_large],
);

# ╔═╡ 1e03184d-ec78-4cd2-8b96-6b37f7406049
md"""
#### SWCS
"""

# ╔═╡ b0c4b9bb-ef7f-4764-84bb-854cf6fd86db
path_swcs_r = joinpath(root_dir_r, "swcs")

# ╔═╡ 0a6e6032-049d-4f1a-af5b-92fc9b463016
df_s_r = CSV.read(string(path_swcs_r, "/full2.csv"), DataFrame);

# ╔═╡ 77eec308-31c7-4acb-b1ba-63cd0d52eef4
df_s_low_r, df_s_normal_r = groupby(df_s_r, :DENSITY);

# ╔═╡ 274372b0-5a2d-4976-863a-2b184dd1b9f0
df_s_low_small_r, df_s_low_medium_r, df_s_low_large_r = groupby(df_s_low_r, :SIZE);

# ╔═╡ c9824d60-dd8d-4aa2-9ad4-83b9487a182d
df_s_normal_small_r, df_s_normal_medium_r, df_s_normal_large_r = groupby(
    df_s_normal_r, :SIZE
);

# ╔═╡ 4e1d4a2c-a6d9-4133-9950-db034a12211a
array_s_r = hcat(
    df_s_r[!, :calculated_swcs_large],
    df_s_r[!, :calculated_swcs_medium],
    df_s_r[!, :calculated_swcs_small],
);

# ╔═╡ d543c8b6-a53f-4c0d-9374-0b0e162bd448
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

# ╔═╡ eead6f88-82d7-4b0a-9028-d34a4af13ef8
df_i_large, df_i_r_large, df_i_medium, df_i_r_medium, df_i_small, df_i_r_small = remove_false_negatives(df_i, array_i, df_i_r, array_i_r);

# ╔═╡ 9050cef9-d977-404c-b301-c44b4a535bcc
begin
	r_squared_reprod_i, rms_values_reprod_i, fitted_line_reprod_i, coefficient_reprod_i = prepare_linear_regression(df_i_r_large, df_i_r_medium, df_i_r_small, df_i_large, df_i_medium, df_i_small)
	
	r_squared_reprod_i = round.(r_squared_reprod_i; digits=3)
	rms_values_reprod_i = round.(rms_values_reprod_i; digits=3)
	coefficient_reprod_i = round.(coefficient_reprod_i; digits=3)
end

# ╔═╡ 80bdd86d-930f-4d6c-8845-dd3680390af0
array_s = hcat(df_s[!, :calculated_swcs_large], df_s[!, :calculated_swcs_medium], df_s[!, :calculated_swcs_small]);

# ╔═╡ 7dcd6c10-9206-4f3d-8c85-11faacf63084
df_s_large, df_s_r_large, df_s_medium, df_s_r_medium, df_s_small, df_s_r_small = remove_false_negatives(df_s, array_s, df_s_r, array_s_r, true);

# ╔═╡ 51649ae7-9035-493f-a7a8-64ce35c0fafd
begin
	r_squared_reprod_s, rms_values_reprod_s, fitted_line_reprod_s, coefficient_reprod_s = prepare_linear_regression(df_s_r_large, df_s_r_medium, df_s_r_small, df_s_large, df_s_medium, df_s_small)
	
	r_squared_reprod_s = round.(r_squared_reprod_s; digits=3)
	rms_values_reprod_s = round.(rms_values_reprod_s; digits=3)
	coefficient_reprod_s = round.(coefficient_reprod_s; digits=3)
end

# ╔═╡ c50117bb-bbb7-4347-b4bb-d446677f4509
md"""
## Sensitivity and Specificity
"""

# ╔═╡ 75e0f9d7-ef66-42bd-bf06-cebecb09a526
std_level = 1.5

# ╔═╡ 9ed516b4-232b-45c5-af90-2964ed2e6a36
md"""
### False Negative
"""

# ╔═╡ 10725f84-8ce4-454a-8c6e-41cdaa25dd1a
md"""
#### SWCS
"""

# ╔═╡ eb4c455f-120b-4c2b-95d6-3022aef521d8
rand() * 10

# ╔═╡ 02dd53b4-3f0f-437f-885b-1c7f4c78d8e5
begin
	false_negative_s = []
	for i in 1:3:nrow(df_s)-2
		mean_s, std_s = mean(df_s[i:i+2, :swcs_bkg]), std(df_s[i:i+2, :swcs_bkg])*std_level 
		# @info mean_s, 10
		array_s = hcat(df_s[i:i+2, :calculated_swcs_large], df_s[i:i+2, :calculated_swcs_medium], df_s[i:i+2, :calculated_swcs_small]);
		neg = length(findall(x -> x <= (mean_s + std_s), array_s))
		push!(false_negative_s, neg)
	end
end

# ╔═╡ 608e3046-5df2-4491-992c-46233e93f4c8
total_zero_s = sum(false_negative_s)

# ╔═╡ 65dc9755-c46f-49ea-8ecf-d68bcb42a628
md"""
#### Agatston
"""

# ╔═╡ 8508a5ed-2646-4b2b-9a3b-42300aa1f4ed
array_a = hcat(df_a[!, :calculated_mass_large], df_a[!, :calculated_mass_medium], df_a[!, :calculated_mass_small]);

# ╔═╡ cdcadf95-cf49-41b4-9936-db56052fea52
df_a_large, df_a_r_large, df_a_medium, df_a_r_medium, df_a_small, df_a_r_small = remove_false_negatives(df_a, array_a, df_a_r, array_a_r);

# ╔═╡ 0769da37-e229-4dc6-b877-eb7e23872176
begin
	r_squared_reprod_a, rms_values_reprod_a, fitted_line_reprod_a, coefficient_reprod_a = prepare_linear_regression(df_a_r_large, df_a_r_medium, df_a_r_small, df_a_large, df_a_medium, df_a_small)
	
	r_squared_reprod_a = round.(r_squared_reprod_a; digits=3)
	rms_values_reprod_a = round.(rms_values_reprod_a; digits=3)
	coefficient_reprod_a = round.(coefficient_reprod_a; digits=3)
end

# ╔═╡ f59e8199-6291-4825-8741-a5c83bca77d5
num_zero_a = length(findall(x -> x <= 0, array_a))

# ╔═╡ 61c58af3-c6af-4eaf-a34c-a4541aa28e02
total_cac = length(array_a)

# ╔═╡ a80dd9c0-fa52-4133-bda7-74a304c02842
length(findall(x -> x <= 0, df_a[!, :calculated_mass_large])), length(findall(x -> x <= 0, df_a[!, :calculated_mass_medium])), length(findall(x -> x <= 0, df_a[!, :calculated_mass_small]))

# ╔═╡ 1561e1a0-2f39-44a0-9c35-e5939ae4ab8d
df_a_ld, df_a_md, df_a_hd = groupby(df_a, :inserts);

# ╔═╡ e64f797b-f282-496d-89f8-fb4ed826323a
length(findall(x -> x <= 0, hcat(df_a_ld[!, :calculated_mass_large], df_a_ld[!, :calculated_mass_medium], df_a_ld[!, :calculated_mass_small]))), length(findall(x -> x <= 0, hcat(df_a_md[!, :calculated_mass_large], df_a_md[!, :calculated_mass_medium], df_a_md[!, :calculated_mass_small]))), length(findall(x -> x <= 0, hcat(df_a_hd[!, :calculated_mass_large], df_a_hd[!, :calculated_mass_medium], df_a_hd[!, :calculated_mass_small])))

# ╔═╡ 84f7225d-f8c6-4754-8635-faacd768abf9
md"""
#### Integrated
"""

# ╔═╡ 958cc5b3-1434-4eed-b7e7-7b3e36fe848c
begin
	false_negative_i = []
	for i in 1:3:nrow(df_i)-2
		mean_i, std_i = mean(df_i[i:i+2, :mass_bkg]), std(df_i[i:i+2, :mass_bkg])*std_level 
		@info mean_i
		array_i = hcat(df_i[i:i+2, :calculated_mass_large], df_i[i:i+2, :calculated_mass_medium], df_i[i:i+2, :calculated_mass_small]);
		neg = length(findall(x -> x <= mean_i + std_i, array_i))
		push!(false_negative_i, neg)
	end
end

# ╔═╡ f029a8c7-7cc1-49c1-8337-c76ac0d585cd
total_zero_i = sum(false_negative_i)

# ╔═╡ 73852846-3fd1-4246-8bd3-afabd7b8e480
md"""
#### Volume Fraction
"""

# ╔═╡ e1665711-ce73-48f3-9321-d9bda350f43b
array_vf = hcat(df_vf[!, :calculated_mass_large], df_vf[!, :calculated_mass_medium], df_vf[!, :calculated_mass_small]);

# ╔═╡ 4ecaa1c8-ca20-455f-8856-24982ce1f687
df_vf_large, df_vf_r_large, df_vf_medium, df_vf_r_medium, df_vf_small, df_vf_r_small = remove_false_negatives(df_vf, array_vf, df_vf_r, array_vf_r);

# ╔═╡ 44c71fcc-f4a1-4ca7-9552-305dfc4e6f01
begin
	r_squared_reprod_vf, rms_values_reprod_vf, fitted_line_reprod_vf, coefficient_reprod_vf = prepare_linear_regression(df_vf_r_large, df_vf_r_medium, df_vf_r_small, df_vf_large, df_vf_medium, df_vf_small)
	
	r_squared_reprod_vf = round.(r_squared_reprod_vf; digits=3)
	rms_values_reprod_vf = round.(rms_values_reprod_vf; digits=3)
	coefficient_reprod_vf = round.(coefficient_reprod_vf; digits=3)
end

# ╔═╡ b6fe3edc-205b-46a2-8cf3-e309bcaec383
function reprod()
    f = Figure()

    ##-- A --##
    ax = Axis(
		f[1, 1],
	    xticks = [0, 50, 100, 150, 200],
	    yticks = [0, 50, 100, 150, 200],
	    xlabel = "Mass 1 (mg)",
	    ylabel = "Mass 2 (mg)",
	    title = "Integrated",
		)
    scatter!(df_i_r[!, :calculated_mass_large], df_i[!, :calculated_mass_large])
    scatter!(df_i_r[!, :calculated_mass_medium], df_i[!, :calculated_mass_medium])
    scatter!(
        df_i_r_small[!, :calculated_mass_small],
        df_i_small[!, :calculated_mass_small];
        color=:red,
    )
    lines!([-1000, 1000], [-1000, 1000]; label="Unity")
    lines!(collect(1:1000), fitted_line_reprod_i; linestyle=:dashdot)
	create_textbox(f[1, 1], coefficient_reprod_i, r_squared_reprod_i, rms_values_reprod_i)

    xlims!(ax; low=0, high=200)
    ylims!(ax; low=0, high=200)

	##-- B --##
    ax = Axis(
		f[2, 1],
		xticks = [0, 50, 100, 150, 200],
    	yticks = [0, 50, 100, 150, 200],
    	xlabel = "Mass 1 (mg)",
    	ylabel = "Mass 2 (mg)",
    	title = "Volume Fraction",
	)
    scatter!(
        df_vf_r_large[!, :calculated_mass_large],
        df_vf_large[!, :calculated_mass_large];
        label="Large Inserts",
    )
    scatter!(
        df_vf_r_medium[!, :calculated_mass_medium],
        df_vf_medium[!, :calculated_mass_medium];
        label="Medium Inserts",
    )
    scatter!(
        df_vf_r_small[!, :calculated_mass_small],
        df_vf_small[!, :calculated_mass_small];
        label="Small Inserts",
        color=:red,
    )
    lines!([-1000, 1000], [-1000, 1000]; label="Unity")
    lines!(collect(1:1000), fitted_line_reprod_vf; linestyle=:dashdot)
	create_textbox(f[2, 1], coefficient_reprod_vf, r_squared_reprod_vf, rms_values_reprod_vf)

    xlims!(ax; low=0, high=200)
    ylims!(ax; low=0, high=200)


    ##-- C --##
    ax = Axis(
		f[1, 2],
		xticks = [0, 50, 100, 150, 200],
    	yticks = [0, 50, 100, 150, 200],
    	xlabel = "Mass 1 (mg)",
    	ylabel = "Mass 2 (mg)",
    	title = "Agatston",
	)
    scatter!(
        df_a_r_large[!, :calculated_mass_large],
        df_a_large[!, :calculated_mass_large];
        label="Large Inserts",
    )
    scatter!(
        df_a_r_medium[!, :calculated_mass_medium],
        df_a_medium[!, :calculated_mass_medium];
        label="Medium Inserts",
    )
    scatter!(
        df_a_r_small[!, :calculated_mass_small],
        df_a_small[!, :calculated_mass_small];
        label="Small Inserts",
        color=:red,
    )
    lines!([-1000, 1000], [-1000, 1000]; label="Unity")
    lines!(collect(1:1000), fitted_line_reprod_a; linestyle=:dashdot)
	create_textbox(f[1, 2], coefficient_reprod_a, r_squared_reprod_a, rms_values_reprod_a)
	
    xlims!(ax; low=0, high=200)
    ylims!(ax; low=0, high=200)


    # ##-- C --##
    ax = Axis(
		f[2, 2],
		xticks = [0, 125, 250, 375, 500],
    	yticks = [0, 125, 250, 375, 500],
    	xlabel = "SWCS 1",
    	ylabel = "SWCS 2",
    	title = "Spatially Weighted",
	)
    scatter!(
        df_s_r_large[!, :calculated_swcs_large],
        df_s_large[!, :calculated_swcs_large];
        label="Large Inserts",
    )
    scatter!(
        df_s_r_medium[!, :calculated_swcs_medium],
        df_s_medium[!, :calculated_swcs_medium];
        label="Medium Inserts",
    )
    scatter!(
        df_s_r_small[!, :calculated_swcs_small],
        df_s_small[!, :calculated_swcs_small];
        label="Small Inserts",
        color=:red,
    )
    lines!([-1000, 1000], [-1000, 1000]; label="Unity")
    lines!(collect(1:1000), fitted_line_reprod_s; linestyle=:dashdot, label="Fitted Line")
	create_textbox(f[2, 2], coefficient_reprod_s, r_squared_reprod_s, rms_values_reprod_s)
	

    xlims!(ax; low=0, high=500)
    ylims!(ax; low=0, high=500)


    ##-- LABELS --##
    f[1:2, 3] = Legend(f, ax; framevisible=false)
    for (label, layout) in zip(["A", "B", "C", "D"], [f[1, 1], f[2, 1], f[1, 2], f[2, 2]])
        Label(
            layout[1, 1, TopLeft()],
            label;
            fontsize=25,
            padding=(0, 0, 40, 0),
            halign=:right,
        )
    end

    save(joinpath(dirname(pwd()),"figures", FIGURE_PATH, "reproducibility.png"), f)
    return f
end

# ╔═╡ e401dc35-3e36-4e69-803d-457b3212db21
reprod()

# ╔═╡ 899c93e4-cd46-4618-89a3-e25d5c80bc5c
mean_vf, std_vf = mean(df_vf[!, :mass_bkg]), std(df_vf[!, :mass_bkg]) * std_level

# ╔═╡ beaa94a0-0de3-40a3-86a1-912cbc018775
begin
	false_negative_vf = []
	for i in 1:3:nrow(df_i)-2
		mean_vf, std_vf = mean(df_vf[i:i+2, :mass_bkg]), std(df_vf[i:i+2, :mass_bkg])*std_level 
		array_vf = hcat(df_vf[i:i+2, :calculated_mass_large], df_vf[i:i+2, :calculated_mass_medium], df_vf[i:i+2, :calculated_mass_small]);
		neg = length(findall(x -> x <= mean_vf + std_vf, array_vf))
		push!(false_negative_vf, neg)
	end
end

# ╔═╡ b708d2c9-f868-4814-b1b7-6908317aa5ba
total_zero_vf = sum(false_negative_vf)

# ╔═╡ 5e32f34c-2f32-457e-a9cf-59184170b6ed
total_zero_i, total_zero_vf, total_zero_s, num_zero_a, total_cac

# ╔═╡ 1d04d29f-401d-45c4-99ff-3bbdb0fb4f6b
md"""
### False Positive
"""

# ╔═╡ 559f77d3-69c0-4155-8cc0-1f669f929bc3
md"""
#### SWCS
"""

# ╔═╡ b540e5f0-d271-4d42-a167-bf79766a53ac
begin
	false_positive_s = []
	for i in 1:3:nrow(df_s)-2
		mean_s, std_s = mean(df_s[i:i+2, :swcs_bkg]), std(df_s[i:i+2, :swcs_bkg])*std_level
		array_s_pos = df_s[i:i+2, :swcs_bkg]
		pos = length(findall(x -> x > (mean_s + std_s), array_s_pos))
		push!(false_positive_s, pos)
	end
end

# ╔═╡ 045cb9a5-12ac-4016-a39c-10b2b15e90f9
total_zero_s_pos = sum(false_positive_s)

# ╔═╡ 1f0df2eb-a434-41fa-b690-4fcdae434a30
md"""
#### Agatston
"""

# ╔═╡ 4a0c4183-3ccc-468a-aa04-75719c21750d
array_a_pos = df_a[!, :mass_bkg]

# ╔═╡ 9d8911ab-05e3-4a07-8f2e-d9e1540447a5
total_cac_pos = length(array_a_pos)

# ╔═╡ f542ed82-d2f3-4033-8083-b2e7a3a54650
total_zero_a_pos = length(findall(x -> x > 0, array_a_pos))

# ╔═╡ 029cffc3-e313-46e6-9313-b50900b4dd1a
md"""
#### Integrated
"""

# ╔═╡ 3512ff80-5109-473b-9c25-7cebe090068a
begin
	false_positive_i = []
	for i in 1:3:nrow(df_i)-2
		mean_i, std_i = mean(df_i[i:i+2, :mass_bkg]), std(df_i[i:i+2, :mass_bkg])*std_level
		array_i_pos = df_i[i:i+2, :mass_bkg]
		pos = length(findall(x -> x > (mean_i + std_level), array_i_pos))
		push!(false_positive_i, pos)
	end
end

# ╔═╡ 3dd31e91-908c-47f1-ae0f-d1031aca6a93
total_zero_i_pos = sum(false_positive_i)

# ╔═╡ c4cff6d0-0d3f-490a-9c35-cb34167e6db3
md"""
#### Volume Fraction
"""

# ╔═╡ 393dbe06-232e-4240-b22d-90016a428ec5
begin
	false_positive_vf = []
	for i in 1:3:nrow(df_vf)-2
		mean_vf, std_vf = mean(df_vf[i:i+2, :mass_bkg]), std(df_vf[i:i+2, :mass_bkg])*std_level
		array_vf_pos = df_vf[i:i+2, :mass_bkg]
		pos = length(findall(x -> x > (mean_vf + std_vf), array_vf_pos))
		push!(false_positive_vf, pos)
	end
end

# ╔═╡ 42c24484-9a06-4a98-9825-43d78ac4e970
total_zero_vf_pos = sum(false_positive_vf)

# ╔═╡ 34991466-181f-4036-8f47-12e20f2144bc
function sensitivity_specificity()
    f = Figure(resolution=(800, 800))
    colors = Makie.wong_colors()

    ##-- A --##
    ax = Axis(
		f[1, 1]; 
		xticks = (1:4, ["Integrated", "Volume Fraction", "Spatially Weighted", "Agatston"]),
		title = "False-Negative (CAC=0)",
		ylabel = "False-Negative (%)",
		yticks = [0, 25, 50, 75, 100]
	)

    table = [1, 2, 3, 4]
	h1 = (total_zero_i / total_cac) * 100
	h2 = (total_zero_vf / total_cac) * 100
    h3 = (total_zero_s / total_cac) * 100
    h4 = (num_zero_a / total_cac) * 100
	heights1 = [h1, h2, h3, h4]
	l1 = @sprintf "%.2f" h1
	l2 = @sprintf "%.2f" h2
	l3 = @sprintf "%.2f" h3
	l4 = @sprintf "%.2f" h4
    barplot!(table, heights1; color=colors[1:4], bar_labels=[l1, l2, l3, l4])

    ylims!(ax; low=0, high=100)

	##-- B --##
	ax = Axis(
		f[2, 1]; 
		xticks = (1:4, ["Integrated", "Volume Fraction", "Spatially Weighted", "Agatston"]),
		title = "False-Positive (CAC>0)",
		ylabel = "False-Positive (%)",
		yticks = [0, 25, 50, 75, 100]
	)

    table = [1, 2, 3, 4]
	h1 = (total_zero_i_pos / total_cac_pos) * 100
	h2 = (total_zero_vf_pos / total_cac_pos) * 100
    h3 = (total_zero_s_pos / total_cac_pos) * 100
    h4 = (total_zero_a_pos / total_cac_pos) * 100
    heights1 = [h1, h2, h3, h4]
	l1 = @sprintf "%.2f" h1
	l2 = @sprintf "%.2f" h2
	l3 = @sprintf "%.2f" h3
	l4 = @sprintf "%.2f" h4
    barplot!(table, heights1; color=colors[1:4], bar_labels=[l1, l2, l3, l4])

    ylims!(ax; low=0, high=100)

	
	for (label, layout) in zip(["A", "B"], [f[1,1], f[2,1]])
	    Label(layout[1, 1, TopLeft()], label,
	        fontsize = 25,
	        padding = (0, 60, 25, 0),
	        halign = :right)
	end

    save(joinpath(dirname(pwd()),"figures", FIGURE_PATH, "sensitivity_specificity.png"), f)


    return f
end


# ╔═╡ bf35fb67-20b9-4b74-bd15-3315d5e37855
with_theme(sensitivity_specificity, medphys_theme)

# ╔═╡ 819c99a9-5a22-447b-9136-0fb3372af1e4
total_zero_i_pos, total_zero_vf_pos, total_zero_s_pos, total_zero_a_pos, total_cac_pos

# ╔═╡ 63a8d206-6b10-4241-b061-50151df24b3b
md"""
# Tables
"""

# ╔═╡ 0998feb7-0b76-4b77-8477-df571a3c4802
md"""
## Accuracy
"""

# ╔═╡ 57c6e0d6-989e-411d-83ad-1c3aaecf9da1
integrated_normal_accuracy = [r_squared_normal_i, rms_values_normal_i..., "y = $(coefficient_normal_i[2])x + $(coefficient_normal_i[1])"]

# ╔═╡ aea852a4-edc6-4291-ac77-945aeea8179b
volume_fraction_normal_accuracy = [r_squared_normal_vf, rms_values_normal_vf..., "y = $(coefficient_normal_vf[2])x + $(coefficient_normal_vf[1])"]

# ╔═╡ 45fb414e-2980-44ba-ac33-7721c2d6ce3d
agatston_normal_accuracy = [r_squared_normal_a, rms_values_normal_a..., "y = $(coefficient_normal_a[2])x + $(coefficient_normal_a[1])"]

# ╔═╡ 15725860-c3f5-4f23-ac6a-67d37690d5cc
integrated_low_accuracy = [r_squared_low_i, rms_values_low_i..., "y = $(coefficient_low_i[2])x + $(coefficient_low_i[1])"]

# ╔═╡ 4e917952-9da3-417e-9fb9-6c05136e1988
volume_fraction_low_accuracy = [r_squared_low_vf, rms_values_low_vf..., "y = $(coefficient_low_vf[2])x + $(coefficient_low_vf[1])"]

# ╔═╡ 33a42461-8497-4782-b091-44eedd1e8539
agatston_low_accuracy = [r_squared_low_a, rms_values_low_a..., "y = $(coefficient_low_a[2])x + $(coefficient_low_a[1])"]

# ╔═╡ 1842bf61-be48-44d4-9711-b37eadb60a47
md"""
## Reproducibility
"""

# ╔═╡ ab2267d7-1b70-4351-bef8-9c1ea1ef08c3
integrated_reprod_accuracy = [r_squared_reprod_i, rms_values_reprod_i..., "y = $(coefficient_reprod_i[2])x + $(coefficient_reprod_i[1])"]

# ╔═╡ aceea802-4159-4488-99e0-293404745464
volume_fraction_reprod_accuracy = [r_squared_reprod_vf, rms_values_reprod_vf..., "y = $(coefficient_reprod_vf[2])x + $(coefficient_reprod_vf[1])"]

# ╔═╡ 39a9de4b-cf74-4bb7-b7f3-b941923d4487
agatston_reprod_accuracy = [r_squared_reprod_a, rms_values_reprod_a..., "y = $(coefficient_reprod_a[2])x + $(coefficient_reprod_a[1])"]

# ╔═╡ 9e6c3f7b-6cb3-42b9-a8b6-8eef2f58da90
swcs_reprod_accuracy = [r_squared_reprod_s, rms_values_reprod_s..., "y = $(coefficient_reprod_s[2])x + $(coefficient_reprod_s[1])"]

# ╔═╡ 85729c9f-90f4-4103-988d-bd04d9ca47a6
md"""
## Sensitivity and Specificity
"""

# ╔═╡ 71d323aa-cba2-46f2-b3a6-8453077fc474
begin
	integrated_sens_spec = ["$total_zero_i / $total_cac", "$total_zero_i_pos / $total_cac_pos"]
	volume_fraction_sens_spec = ["$total_zero_vf / $total_cac", "$total_zero_vf_pos / $total_cac_pos"]
	swcs_sens_spec = ["$total_zero_s / $total_cac", "$total_zero_s_pos / $total_cac_pos"]
	agatston_sens_spec = ["$num_zero_a / $total_cac", "$total_zero_a_pos / $total_cac_pos"]
end

# ╔═╡ 52eb314c-23aa-4288-aa6f-3e1790832d70
swcs_sens_spec

# ╔═╡ 1f3405c3-c76d-4d64-85ba-97057cc9c582
md"""
# Summaries
"""

# ╔═╡ 4bf820cc-5c76-4a95-b041-650909abedcf
begin
	df_accuracy_normal = DataFrame(
		"Metrics" => ["r^2", "RMSE", "RMSD", "Best Fit Line"],
		"Integrated" => integrated_normal_accuracy,
		"Volume Fraction" => volume_fraction_normal_accuracy,
		"Agatston" => agatston_normal_accuracy,
	)
	df_accuracy_low = DataFrame(
		"Metrics" => ["r^2", "RMSE", "RMSD", "Best Fit Line"],
		"Integrated" =>  integrated_low_accuracy,
		"Volume Fraction" => volume_fraction_low_accuracy,
		"Agatston" => agatston_low_accuracy,
	)
	df_reprod = DataFrame(
		"Metrics" => ["r^2", "RMSE", "RMSD", "Best Fit Line"],
		"Integrated" => integrated_reprod_accuracy,
		"Volume Fraction" => volume_fraction_reprod_accuracy,
		"SWCS" => swcs_reprod_accuracy,
		"Agatston" => agatston_reprod_accuracy,
	)
	df_sens_spec = DataFrame(
		"Metrics" => ["False Negatives", "False Positives"],
		"Integrated" => integrated_sens_spec,
		"Volume Fraction" => volume_fraction_sens_spec,
		"SWCS" => swcs_sens_spec,
		"Agatston" => agatston_sens_spec
	)
	_dfs = append!(df_accuracy_normal, df_accuracy_low, cols=:union)
	_dfs = append!(_dfs, df_reprod, cols=:union)
	dfs = append!(_dfs, df_sens_spec, cols=:union)
	
	results = [repeat(["Accuracy (Normal Density)"], 4)..., repeat(["Accuracy (Low Density)"], 4)..., repeat(["Reproducibility"], 4)..., "Specificity", "Sensitivity"]
	insertcols!(dfs, 1, :Results=>results)
end;

# ╔═╡ e1c698de-8acb-44e9-acf7-514c7e5d806e
dfs

# ╔═╡ f832169d-fa15-4507-8074-748a77c6d53d
CSV.write(joinpath(dirname(pwd()),"figures", FIGURE_PATH, "summary.csv"), dfs)

# ╔═╡ 2eb8e060-0086-4b88-b6d6-f851aa04f057
md"""
## Simulation Parameters
"""

# ╔═╡ 949b37e4-a965-4543-ab29-1050bdccd605
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

# ╔═╡ bf93f3c0-24a8-4caa-b721-2632ce52dfbe
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

# ╔═╡ 72c24c61-9b9e-44b4-894b-f08e51f84b01
simulation_parameters = DataFrame(
	"Parameter" => parameters,
	"Simulation" => simulations
)

# ╔═╡ 1f14d642-6bcd-4564-a995-985892360e59
md"""
## Spatially Weighted Calcium Scoring Means & Stds
"""

# ╔═╡ 28852770-07b7-435c-a9c6-bc8162b4d437
kvp = [
	80
	100
	120
	135
]

# ╔═╡ 07a850f2-7476-47d9-a778-3f50a1fbd862
means = [
	201.4898
	174.3658
	158.2645
	152.8815
]

# ╔═╡ 23432875-ffb4-4dfb-8a8c-dfe5c4476a1d
stds = [
	34.3038
	28.1708
	22.9656
	21.0778
]

# ╔═╡ da160a8a-6519-4e4e-90cf-962781f244de
means_stds = DataFrame(
	"kVp" => kvp,
	"Mean" => means,
	"Std" => stds
)

# ╔═╡ Cell order:
# ╠═c82bb3e6-e052-11ec-3392-eba4a4b464ba
# ╠═2df207ed-cf9c-4d7d-8354-4da14f93276c
# ╟─d72c94c2-a208-430c-9d76-85d8cfb33a22
# ╟─cabe7e9a-e932-406b-95dd-2c9128decdc7
# ╟─4a255a58-b9a0-4175-a750-b2562361631d
# ╟─7fc65de7-7514-4532-803e-1120150fd2fb
# ╟─1ecb8fd3-a36b-4b06-8daf-e598e297ecf8
# ╟─c5567a17-2345-470f-97c8-9383617298ac
# ╠═043dec5f-07a8-472d-8082-c7a771f93270
# ╟─e79c34d0-bf06-4a7a-93b7-e49f09ce3f4e
# ╟─4cf61487-3913-4e11-970e-4ca31a5ffc8d
# ╠═7e1b5b89-ed32-42df-aad9-f38524d36988
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
# ╟─593d7b00-94c4-49c2-87a7-723dbef64667
# ╠═63f97b37-dc0e-456d-811f-409f8360b2da
# ╟─ac937ee5-e820-4960-836f-50df8c7ae3ee
# ╟─decd4320-112a-4c4c-9ed5-680a721b81c1
# ╟─0fa91a0b-c205-4b10-9876-9a9e163f3d7e
# ╟─bc9939d7-f3a8-48e5-ada4-8fcf0fc1137a
# ╟─3f02194d-9c69-45d1-80fc-4cc50b1f3d8f
# ╟─9b3c6bda-5c21-4c07-be6f-939a1ac0496a
# ╠═b4e25c03-ac59-4104-acf5-2f049497d460
# ╠═074a6437-9745-4680-b9b5-baf9d767afc9
# ╠═f51b5a1e-a688-470d-8a83-b65f03fff875
# ╟─91a157fb-a515-4ab0-be1a-8d9148db6d67
# ╠═1478eef5-bcde-4513-8794-9b9cbe388ee2
# ╠═ffeb8498-d7a3-4122-b71f-56326490042b
# ╠═ad784be1-5252-4914-bdd2-085edf5a95e0
# ╟─d4753644-29fb-41a9-b8d1-5f46bb87a613
# ╠═aef525b3-be1a-4a80-ba18-7ade91baea2d
# ╟─b6fe3edc-205b-46a2-8cf3-e309bcaec383
# ╟─e401dc35-3e36-4e69-803d-457b3212db21
# ╟─f8a69b2e-1cab-4281-9eb5-22bed0046c49
# ╠═af6d4416-6750-4674-8d5d-88456d58f841
# ╠═2e531f2f-6d07-4818-b9f8-ecfa8a5d419a
# ╠═b7ad23f1-11d4-4146-8eae-37c0856b6e01
# ╠═a1a417b0-2f19-4911-bec7-4bac0d2358b7
# ╠═c31d03f6-3846-4f78-aa76-97a397af500a
# ╠═74c1846e-4e69-4434-a695-7b197ed0df0a
# ╠═52261c32-415d-425f-898b-d1bbf66052cf
# ╠═eead6f88-82d7-4b0a-9028-d34a4af13ef8
# ╠═9050cef9-d977-404c-b301-c44b4a535bcc
# ╟─f7b5997f-8aba-4311-b120-836e808fee2f
# ╠═47ff39ea-1f55-4e7d-94de-5d4ded126648
# ╠═24e3e4d1-91b4-4ae0-abf8-174118926b6b
# ╠═fc1e3216-6929-4f8f-b1c9-b700f3e35181
# ╠═6ea76780-2726-4d61-bff1-fb590f7ba2f5
# ╠═288f8d4c-323a-4962-8aff-026ef1bfa77d
# ╠═fd9e7b94-03c5-4ea9-8405-2218e07df3b5
# ╠═4ecaa1c8-ca20-455f-8856-24982ce1f687
# ╠═44c71fcc-f4a1-4ca7-9552-305dfc4e6f01
# ╟─3be21a82-9524-4edb-93d3-4d515b8b8758
# ╠═04afe0fa-c290-48b3-aa68-c8f74346a321
# ╠═d181bc5e-468d-42e0-8284-986756a1368a
# ╠═814ff1b6-103c-44b9-b9fd-f334507bf998
# ╠═92daece0-acb9-415d-9dbc-d773bffdc18c
# ╠═305a1661-e34d-41ff-912a-386be27aa76d
# ╠═ce630749-7c8e-41c7-853e-422092d64958
# ╠═cdcadf95-cf49-41b4-9936-db56052fea52
# ╠═0769da37-e229-4dc6-b877-eb7e23872176
# ╟─1e03184d-ec78-4cd2-8b96-6b37f7406049
# ╠═b0c4b9bb-ef7f-4764-84bb-854cf6fd86db
# ╠═0a6e6032-049d-4f1a-af5b-92fc9b463016
# ╠═77eec308-31c7-4acb-b1ba-63cd0d52eef4
# ╠═274372b0-5a2d-4976-863a-2b184dd1b9f0
# ╠═c9824d60-dd8d-4aa2-9ad4-83b9487a182d
# ╠═4e1d4a2c-a6d9-4133-9950-db034a12211a
# ╠═d543c8b6-a53f-4c0d-9374-0b0e162bd448
# ╠═80bdd86d-930f-4d6c-8845-dd3680390af0
# ╠═7dcd6c10-9206-4f3d-8c85-11faacf63084
# ╠═51649ae7-9035-493f-a7a8-64ce35c0fafd
# ╟─c50117bb-bbb7-4347-b4bb-d446677f4509
# ╠═34991466-181f-4036-8f47-12e20f2144bc
# ╠═bf35fb67-20b9-4b74-bd15-3315d5e37855
# ╠═75e0f9d7-ef66-42bd-bf06-cebecb09a526
# ╟─9ed516b4-232b-45c5-af90-2964ed2e6a36
# ╠═5e32f34c-2f32-457e-a9cf-59184170b6ed
# ╠═819c99a9-5a22-447b-9136-0fb3372af1e4
# ╟─10725f84-8ce4-454a-8c6e-41cdaa25dd1a
# ╠═eb4c455f-120b-4c2b-95d6-3022aef521d8
# ╠═02dd53b4-3f0f-437f-885b-1c7f4c78d8e5
# ╠═608e3046-5df2-4491-992c-46233e93f4c8
# ╟─65dc9755-c46f-49ea-8ecf-d68bcb42a628
# ╠═8508a5ed-2646-4b2b-9a3b-42300aa1f4ed
# ╠═f59e8199-6291-4825-8741-a5c83bca77d5
# ╠═61c58af3-c6af-4eaf-a34c-a4541aa28e02
# ╠═a80dd9c0-fa52-4133-bda7-74a304c02842
# ╠═1561e1a0-2f39-44a0-9c35-e5939ae4ab8d
# ╠═e64f797b-f282-496d-89f8-fb4ed826323a
# ╟─84f7225d-f8c6-4754-8635-faacd768abf9
# ╠═958cc5b3-1434-4eed-b7e7-7b3e36fe848c
# ╠═f029a8c7-7cc1-49c1-8337-c76ac0d585cd
# ╟─73852846-3fd1-4246-8bd3-afabd7b8e480
# ╠═e1665711-ce73-48f3-9321-d9bda350f43b
# ╠═899c93e4-cd46-4618-89a3-e25d5c80bc5c
# ╠═beaa94a0-0de3-40a3-86a1-912cbc018775
# ╠═b708d2c9-f868-4814-b1b7-6908317aa5ba
# ╟─1d04d29f-401d-45c4-99ff-3bbdb0fb4f6b
# ╟─559f77d3-69c0-4155-8cc0-1f669f929bc3
# ╠═b540e5f0-d271-4d42-a167-bf79766a53ac
# ╠═045cb9a5-12ac-4016-a39c-10b2b15e90f9
# ╟─1f0df2eb-a434-41fa-b690-4fcdae434a30
# ╠═4a0c4183-3ccc-468a-aa04-75719c21750d
# ╠═9d8911ab-05e3-4a07-8f2e-d9e1540447a5
# ╠═f542ed82-d2f3-4033-8083-b2e7a3a54650
# ╟─029cffc3-e313-46e6-9313-b50900b4dd1a
# ╠═3512ff80-5109-473b-9c25-7cebe090068a
# ╠═3dd31e91-908c-47f1-ae0f-d1031aca6a93
# ╟─c4cff6d0-0d3f-490a-9c35-cb34167e6db3
# ╠═393dbe06-232e-4240-b22d-90016a428ec5
# ╠═42c24484-9a06-4a98-9825-43d78ac4e970
# ╟─63a8d206-6b10-4241-b061-50151df24b3b
# ╟─0998feb7-0b76-4b77-8477-df571a3c4802
# ╠═57c6e0d6-989e-411d-83ad-1c3aaecf9da1
# ╠═aea852a4-edc6-4291-ac77-945aeea8179b
# ╠═45fb414e-2980-44ba-ac33-7721c2d6ce3d
# ╠═15725860-c3f5-4f23-ac6a-67d37690d5cc
# ╠═4e917952-9da3-417e-9fb9-6c05136e1988
# ╠═33a42461-8497-4782-b091-44eedd1e8539
# ╟─1842bf61-be48-44d4-9711-b37eadb60a47
# ╠═ab2267d7-1b70-4351-bef8-9c1ea1ef08c3
# ╠═aceea802-4159-4488-99e0-293404745464
# ╠═39a9de4b-cf74-4bb7-b7f3-b941923d4487
# ╠═9e6c3f7b-6cb3-42b9-a8b6-8eef2f58da90
# ╟─85729c9f-90f4-4103-988d-bd04d9ca47a6
# ╠═71d323aa-cba2-46f2-b3a6-8453077fc474
# ╠═52eb314c-23aa-4288-aa6f-3e1790832d70
# ╟─1f3405c3-c76d-4d64-85ba-97057cc9c582
# ╠═4bf820cc-5c76-4a95-b041-650909abedcf
# ╠═e1c698de-8acb-44e9-acf7-514c7e5d806e
# ╠═f832169d-fa15-4507-8074-748a77c6d53d
# ╟─2eb8e060-0086-4b88-b6d6-f851aa04f057
# ╟─949b37e4-a965-4543-ab29-1050bdccd605
# ╟─bf93f3c0-24a8-4caa-b721-2632ce52dfbe
# ╠═72c24c61-9b9e-44b4-894b-f08e51f84b01
# ╟─1f14d642-6bcd-4564-a995-985892360e59
# ╠═28852770-07b7-435c-a9c6-bc8162b4d437
# ╠═07a850f2-7476-47d9-a778-3f50a1fbd862
# ╠═23432875-ffb4-4dfb-8a8c-dfe5c4476a1d
# ╠═da160a8a-6519-4e4e-90cf-962781f244de

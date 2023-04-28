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

# ╔═╡ aaa022d9-2d23-4bc3-96fc-dc349848820e
md"""
### Density
"""

# ╔═╡ e3c5f720-adec-4b59-b8a7-a26385c002ff
function lin_reg_agat_swcs()
    f = Figure()

    ##-- A --##
    ax = Axis(
		f[1, 1],
    	xlabel = "SWCS Score",
    	ylabel = "Agatston Score",
    	title = "Agatston vs SWCS (Normal Density)"
	)
    sc1 = scatter!(
		df_s_normal[!, :calculated_swcs_large], df_a_normal[!, :calculated_agat_large]
	)
    sc1 = scatter!(
		df_s_normal[!, :calculated_swcs_medium], df_a_normal[!, :calculated_agat_medium]
	)
    sc1 = scatter!(
		df_s_normal[!, :calculated_swcs_small], df_a_normal[!, :calculated_agat_small]; color=:red
	)

	
    ##-- B --##
    ax = Axis(
		f[2, 1],
    	xlabel = "SWCS Score",
    	ylabel = "Agatston Score",
    	title = "Agatston vs SWCS (Low Density)"
	)
    sc1 = scatter!(
		df_s_low[!, :calculated_swcs_large], df_a_low[!, :calculated_agat_large]
	)
    sc2 = scatter!(
		df_s_low[!, :calculated_swcs_medium], df_a_low[!, :calculated_agat_medium]
	)
    sc3 = scatter!(
		df_s_low[!, :calculated_swcs_small], df_a_low[!, :calculated_agat_small]; color=:red
	)

    ##-- LABELS --##
    f[2, 2] = Legend(
        f,
        [sc1, sc2, sc3],
        ["Large Inserts", "Medium Inserts", "Small Inserts"];
        framevisible=false,
		padding=(10, 10, 8, 8)
    )

    for (label, layout) in zip(["A", "B"], [f[1, 1], f[2, 1]])
        Label(
            layout[1, 1, TopLeft()],
            label;
            fontsize=25,
            padding=(0, 90, 25, 0),
            halign=:right,
        )
    end

	f
end

# ╔═╡ 0822b8d6-6285-4b69-80e1-23538044e6cc
with_theme(lin_reg_agat_swcs, medphys_theme)

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

# ╔═╡ f8e0735d-7854-48d0-85a9-ccfb3f8fc898
r_squared_normal_i, r_squared_normal_vf, r_squared_normal_a 

# ╔═╡ 44341e09-1593-4d24-8c77-1e3fceb219ee
coefficient_normal_i, coefficient_normal_vf, coefficient_normal_a

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

# ╔═╡ 766474b7-8158-48b5-a4fb-6ed4ab055338
round.((r_squared_low_i, r_squared_low_vf, r_squared_low_a), digits=2)

# ╔═╡ 3cf4a3a1-0997-425c-9f4f-0785a45eea1c
coefficient_low_i, coefficient_low_vf, coefficient_low_a

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

# ╔═╡ 77f50463-59bc-477d-8a26-97aa08b0e21f
md"""
### Patient Size
"""

# ╔═╡ 086bdda7-cdf5-4510-95e9-c405a2acb9d2
begin 
	df_i_small_patient, df_i_medium_patient, df_i_large_patient = groupby(df_i, :SIZE)
	df_vf_small_patient, df_vf_medium_patient, df_vf_large_patient = groupby(df_vf, :SIZE)
	df_a_small_patient, df_a_medium_patient, df_a_large_patient = groupby(df_a, :SIZE)
	df_s_small_patient, df_s_medium_patient, df_s_large_patient = groupby(df_s, :SIZE)
end;

# ╔═╡ 5b2b1e40-d75e-4281-9ddd-4445775b79ba
md"""
#### Small Patient
"""

# ╔═╡ 8eac91f6-7872-4065-bc78-6ed611b2b192
begin
	r_squared_normal_i_sp, rms_values_normal_i_sp, fitted_line_normal_i_sp, coefficient_normal_i_sp = prepare_linear_regression(df_i_small_patient)
	r_squared_normal_i_sp = round.(r_squared_normal_i_sp; digits=3)
	rms_values_normal_i_sp = round.(rms_values_normal_i_sp; digits=3)
	coefficient_normal_i_sp = round.(coefficient_normal_i_sp; digits=3)
end

# ╔═╡ 46a5b2bd-3f00-4e3a-b799-f1a9f1041928
begin
	r_squared_normal_vf_sp, rms_values_normal_vf_sp, fitted_line_normal_vf_sp, coefficient_normal_vf_sp = prepare_linear_regression(df_vf_small_patient)
	r_squared_normal_vf_sp = round.(r_squared_normal_vf_sp; digits=3)
	rms_values_normal_vf_sp = round.(rms_values_normal_vf_sp; digits=3)
	coefficient_normal_vf_sp = round.(coefficient_normal_vf_sp; digits=3)
end

# ╔═╡ 7dd8bbb5-80dc-490d-b518-d7b14452e305
begin
	r_squared_normal_a_sp, rms_values_normal_a_sp, fitted_line_normal_a_sp, coefficient_normal_a_sp = prepare_linear_regression(df_a_small_patient)
	r_squared_normal_a_sp = round.(r_squared_normal_a_sp; digits=3)
	rms_values_normal_a_sp = round.(rms_values_normal_a_sp; digits=3)
	coefficient_normal_a_sp = round.(coefficient_normal_a_sp; digits=3)
end

# ╔═╡ c41d53f9-a74b-4096-8709-5a92c2bc7eaa
function lin_reg_small_patient()
    f = Figure()

    ##-- A --##
    ax1 = Axis(f[1, 1])
    df = df_i_small_patient
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
    lines!(ax1, collect(1:1000), fitted_line_normal_i_sp; linestyle=:dashdot)
	create_textbox(f[1, 1], coefficient_normal_i_sp, r_squared_normal_i_sp, rms_values_normal_i_sp)

    xlims!(ax1; low=0, high=200)
    ylims!(ax1; low=0, high=200)
    ax1.xticks = [0, 50, 100, 150, 200]
    ax1.yticks = [0, 100, 200]
    ax1.xlabel = "Known Mass (mg)"
    ax1.ylabel = "Calculated Mass (mg)"
    ax1.title = "Integrated (Small Patient)"
	
    ##-- B --##
    ax2 = Axis(f[2, 1])
    df3 = df_vf_small_patient
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
    ln2 = lines!(ax2, collect(1:1000), fitted_line_normal_vf_sp; linestyle=:dashdot)
	create_textbox(f[2, 1], coefficient_normal_vf_sp, r_squared_normal_vf_sp, rms_values_normal_vf_sp)
	
    xlims!(ax2; low=0, high=200)
    ylims!(ax2; low=0, high=200)
    ax2.xticks = [0, 50, 100, 150, 200]
    ax2.yticks = [0, 100, 200]
    ax2.xlabel = "Known Mass (mg)"
    ax2.ylabel = "Calculated Mass (mg)"
    ax2.title = "Volume Fraction (Small Patient)"


    ##-- C --##
    ax2 = Axis(f[3, 1])
    df3 = df_vf_small_patient
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
    ln2 = lines!(ax2, collect(1:1000), fitted_line_normal_a_sp; linestyle=:dashdot)
	create_textbox(f[3, 1], coefficient_normal_a_sp, r_squared_normal_a_sp, rms_values_normal_a_sp)

    xlims!(ax2; low=0, high=200)
    ylims!(ax2; low=0, high=200)
    ax2.xticks = [0, 50, 100, 150, 200]
    ax2.yticks = [0, 100, 200]
    ax2.xlabel = "Known Mass (mg)"
    ax2.ylabel = "Calculated Mass (mg)"
    ax2.title = "Agatston (Small Patient)"

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

    # save(joinpath(dirname(pwd()),"figures", FIGURE_PATH, "accuracy_normal.png"), f)
    return f
end

# ╔═╡ 971ba7ff-c14a-4b24-8e3c-5df3f68a05fd
with_theme(lin_reg_small_patient, medphys_theme)

# ╔═╡ 4cbeb9d2-4ed1-4190-8fd0-aeff30f3c81d
md"""
#### Medium Patient
"""

# ╔═╡ 398e4d28-d01b-4bd7-8145-7d78a3de672f
begin
	r_squared_normal_i_mp, rms_values_normal_i_mp, fitted_line_normal_i_mp, coefficient_normal_i_mp = prepare_linear_regression(df_i_medium_patient)
	r_squared_normal_i_mp = round.(r_squared_normal_i_mp; digits=3)
	rms_values_normal_i_mp = round.(rms_values_normal_i_mp; digits=3)
	coefficient_normal_i_mp = round.(coefficient_normal_i_mp; digits=3)
end

# ╔═╡ a6f72994-a3e1-4bb0-963f-12f8c2788ad8
begin
	r_squared_normal_vf_mp, rms_values_normal_vf_mp, fitted_line_normal_vf_mp, coefficient_normal_vf_mp = prepare_linear_regression(df_vf_medium_patient)
	r_squared_normal_vf_mp = round.(r_squared_normal_vf_mp; digits=3)
	rms_values_normal_vf_mp = round.(rms_values_normal_vf_mp; digits=3)
	coefficient_normal_vf_mp = round.(coefficient_normal_vf_mp; digits=3)
end

# ╔═╡ 134fc008-1eb8-496c-a079-0b6cf177a3d1
begin
	r_squared_normal_a_mp, rms_values_normal_a_mp, fitted_line_normal_a_mp, coefficient_normal_a_mp = prepare_linear_regression(df_a_medium_patient)
	r_squared_normal_a_mp = round.(r_squared_normal_a_mp; digits=3)
	rms_values_normal_a_mp = round.(rms_values_normal_a_mp; digits=3)
	coefficient_normal_a_mp = round.(coefficient_normal_a_mp; digits=3)
end

# ╔═╡ d5a323d0-5ca7-46f5-b169-1e17c4ec0413
function lin_reg_medium_patient()
    f = Figure()

    ##-- A --##
    ax1 = Axis(f[1, 1])
    df = df_i_medium_patient
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
    lines!(ax1, collect(1:1000), fitted_line_normal_i_mp; linestyle=:dashdot)
	create_textbox(f[1, 1], coefficient_normal_i_mp, r_squared_normal_i_mp, rms_values_normal_i_mp)

    xlims!(ax1; low=0, high=200)
    ylims!(ax1; low=0, high=200)
    ax1.xticks = [0, 50, 100, 150, 200]
    ax1.yticks = [0, 100, 200]
    ax1.xlabel = "Known Mass (mg)"
    ax1.ylabel = "Calculated Mass (mg)"
    ax1.title = "Integrated (Medium Patient)"
	
    ##-- B --##
    ax2 = Axis(f[2, 1])
    df3 = df_vf_medium_patient
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
    ln2 = lines!(ax2, collect(1:1000), fitted_line_normal_vf_mp; linestyle=:dashdot)
	create_textbox(f[2, 1], coefficient_normal_vf_mp, r_squared_normal_vf_mp, rms_values_normal_vf_mp)
	
    xlims!(ax2; low=0, high=200)
    ylims!(ax2; low=0, high=200)
    ax2.xticks = [0, 50, 100, 150, 200]
    ax2.yticks = [0, 100, 200]
    ax2.xlabel = "Known Mass (mg)"
    ax2.ylabel = "Calculated Mass (mg)"
    ax2.title = "Volume Fraction (Medium Patient)"


    ##-- C --##
    ax2 = Axis(f[3, 1])
    df3 = df_vf_medium_patient
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
    ln2 = lines!(ax2, collect(1:1000), fitted_line_normal_a_mp; linestyle=:dashdot)
	create_textbox(f[3, 1], coefficient_normal_a_mp, r_squared_normal_a_mp, rms_values_normal_a_mp)

    xlims!(ax2; low=0, high=200)
    ylims!(ax2; low=0, high=200)
    ax2.xticks = [0, 50, 100, 150, 200]
    ax2.yticks = [0, 100, 200]
    ax2.xlabel = "Known Mass (mg)"
    ax2.ylabel = "Calculated Mass (mg)"
    ax2.title = "Agatston (Medium Patient)"

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

    # save(joinpath(dirname(pwd()),"figures", FIGURE_PATH, "accuracy_normal.png"), f)
    return f
end

# ╔═╡ 89ace836-d183-4872-8037-c48a5db02b89
with_theme(lin_reg_medium_patient, medphys_theme)

# ╔═╡ 7ae5193f-2d29-4152-b062-791629a5f1fa
md"""
#### Large Patient
"""

# ╔═╡ 6e2f2819-6421-416f-a097-df1faefd9748
begin
	r_squared_normal_i_lp, rms_values_normal_i_lp, fitted_line_normal_i_lp, coefficient_normal_i_lp = prepare_linear_regression(df_i_large_patient)
	r_squared_normal_i_lp = round.(r_squared_normal_i_lp; digits=3)
	rms_values_normal_i_lp = round.(rms_values_normal_i_lp; digits=3)
	coefficient_normal_i_lp = round.(coefficient_normal_i_lp; digits=3)
end

# ╔═╡ e6f08aec-06e2-4dd4-86c6-0e87db2d6f88
begin
	r_squared_normal_vf_lp, rms_values_normal_vf_lp, fitted_line_normal_vf_lp, coefficient_normal_vf_lp = prepare_linear_regression(df_vf_large_patient)
	r_squared_normal_vf_lp = round.(r_squared_normal_vf_lp; digits=3)
	rms_values_normal_vf_lp = round.(rms_values_normal_vf_lp; digits=3)
	coefficient_normal_vf_lp = round.(coefficient_normal_vf_lp; digits=3)
end

# ╔═╡ 384f8f09-6df0-4fd8-91a5-f901001afea2
begin
	r_squared_normal_a_lp, rms_values_normal_a_lp, fitted_line_normal_a_lp, coefficient_normal_a_lp = prepare_linear_regression(df_a_large_patient)
	r_squared_normal_a_lp = round.(r_squared_normal_a_lp; digits=3)
	rms_values_normal_a_lp = round.(rms_values_normal_a_lp; digits=3)
	coefficient_normal_a_lp = round.(coefficient_normal_a_lp; digits=3)
end

# ╔═╡ a28a1259-9161-4c18-82f2-9bfb26a50690
function lin_reg_large_patient()
    f = Figure()

    ##-- A --##
    ax1 = Axis(f[1, 1])
    df = df_i_large_patient
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
    lines!(ax1, collect(1:1000), fitted_line_normal_i_lp; linestyle=:dashdot)
	create_textbox(f[1, 1], coefficient_normal_i_lp, r_squared_normal_i_lp, rms_values_normal_i_lp)

    xlims!(ax1; low=0, high=200)
    ylims!(ax1; low=0, high=200)
    ax1.xticks = [0, 50, 100, 150, 200]
    ax1.yticks = [0, 100, 200]
    ax1.xlabel = "Known Mass (mg)"
    ax1.ylabel = "Calculated Mass (mg)"
    ax1.title = "Integrated (Large Patient)"
	
    ##-- B --##
    ax2 = Axis(f[2, 1])
    df3 = df_vf_large_patient
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
    ln2 = lines!(ax2, collect(1:1000), fitted_line_normal_vf_lp; linestyle=:dashdot)
	create_textbox(f[2, 1], coefficient_normal_vf_lp, r_squared_normal_vf_lp, rms_values_normal_vf_lp)
	
    xlims!(ax2; low=0, high=200)
    ylims!(ax2; low=0, high=200)
    ax2.xticks = [0, 50, 100, 150, 200]
    ax2.yticks = [0, 100, 200]
    ax2.xlabel = "Known Mass (mg)"
    ax2.ylabel = "Calculated Mass (mg)"
    ax2.title = "Volume Fraction (Large Patient)"


    ##-- C --##
    ax2 = Axis(f[3, 1])
    df3 = df_vf_large_patient
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
    ln2 = lines!(ax2, collect(1:1000), fitted_line_normal_a_lp; linestyle=:dashdot)
	create_textbox(f[3, 1], coefficient_normal_a_lp, r_squared_normal_a_lp, rms_values_normal_a_lp)

    xlims!(ax2; low=0, high=200)
    ylims!(ax2; low=0, high=200)
    ax2.xticks = [0, 50, 100, 150, 200]
    ax2.yticks = [0, 100, 200]
    ax2.xlabel = "Known Mass (mg)"
    ax2.ylabel = "Calculated Mass (mg)"
    ax2.title = "Agatston (Large Patient)"

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

    # save(joinpath(dirname(pwd()),"figures", FIGURE_PATH, "accuracy_normal.png"), f)
    return f
end

# ╔═╡ 2555e6e0-2712-4ba0-a1f8-38a8e89b025b
with_theme(lin_reg_large_patient, medphys_theme)

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

# ╔═╡ 598e1797-837f-4bef-8fac-716f0690142a
array_vf = hcat(df_vf[:, :calculated_mass_large], df_vf[:, :calculated_mass_medium], df_vf[:, :calculated_mass_small]);

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

# ╔═╡ 4ecaa1c8-ca20-455f-8856-24982ce1f687
df_vf_large, df_vf_r_large, df_vf_medium, df_vf_r_medium, df_vf_small, df_vf_r_small = remove_false_negatives(df_vf, array_vf, df_vf_r, array_vf_r);

# ╔═╡ 44c71fcc-f4a1-4ca7-9552-305dfc4e6f01
begin
	r_squared_reprod_vf, rms_values_reprod_vf, fitted_line_reprod_vf, coefficient_reprod_vf = prepare_linear_regression(df_vf_r_large, df_vf_r_medium, df_vf_r_small, df_vf_large, df_vf_medium, df_vf_small)
	
	r_squared_reprod_vf = round.(r_squared_reprod_vf; digits=3)
	rms_values_reprod_vf = round.(rms_values_reprod_vf; digits=3)
	coefficient_reprod_vf = round.(coefficient_reprod_vf; digits=3)
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

# ╔═╡ 473918d5-df19-439d-a125-e0241495d29d
md"""
### Density
"""

# ╔═╡ fc3c711f-dee1-4e88-9337-9d930a40e3ea
md"""
#### Normal Density
"""

# ╔═╡ a50a0ed7-5733-44ec-857b-710156ce40e3
md"""
##### Integrated
"""

# ╔═╡ 4877368a-2491-4f45-8d52-5050e335ad4d
array_i_nd = hcat(
	df_i_normal[!, :calculated_mass_large], 
	df_i_normal[!, :calculated_mass_medium], 
	df_i_normal[!, :calculated_mass_small]
);

# ╔═╡ 6720e3d7-6e6a-4a7e-b749-f4777cba587c
array_i_r_nd = hcat(
    df_i_normal_r[!, :calculated_mass_large],
    df_i_normal_r[!, :calculated_mass_medium],
    df_i_normal_r[!, :calculated_mass_small],
);

# ╔═╡ a5409851-8bc1-434a-be9f-d6a5bd8e392a
df_i_large_nd, df_i_r_large_nd, df_i_medium_nd, df_i_r_medium_nd, df_i_small_nd, df_i_r_small_nd = remove_false_negatives(df_i_normal, array_i_nd, df_i_normal_r, array_i_r_nd);

# ╔═╡ 1309668f-d000-4062-a582-a90e23ef43f8
begin
	r_squared_reprod_i_nd, rms_values_reprod_i_nd, fitted_line_reprod_i_nd, coefficient_reprod_i_nd = prepare_linear_regression(df_i_r_large_nd, df_i_r_medium_nd, df_i_r_small_nd, df_i_large_nd, df_i_medium_nd, df_i_small_nd)
	
	r_squared_reprod_i_nd = round.(r_squared_reprod_i_nd; digits=3)
	rms_values_reprod_i_nd = round.(rms_values_reprod_i_nd; digits=3)
	coefficient_reprod_i_nd = round.(coefficient_reprod_i_nd; digits=3)
end

# ╔═╡ 768e6d0a-9c75-4b2f-b52a-bb68b834f68e
md"""
##### Volume Fraction
"""

# ╔═╡ ac366a6b-8384-4735-9aa5-51c7f7604877
array_vf_nd = hcat(
	df_vf_normal[!, :calculated_mass_large], 
	df_vf_normal[!, :calculated_mass_medium], 
	df_vf_normal[!, :calculated_mass_small]
);

# ╔═╡ 3802ae31-ae12-45ba-b4b5-2854729b1ff0
array_vf_r_nd = hcat(
    df_vf_normal_r[!, :calculated_mass_large],
    df_vf_normal_r[!, :calculated_mass_medium],
    df_vf_normal_r[!, :calculated_mass_small],
);

# ╔═╡ ae19ce13-bba7-41cb-8be2-2e56ce7fa879
df_vf_large_nd, df_vf_r_large_nd, df_vf_medium_nd, df_vf_r_medium_nd, df_vf_small_nd, df_vf_r_small_nd = remove_false_negatives(df_vf_normal, array_vf_nd, df_vf_normal_r, array_vf_r_nd);

# ╔═╡ 2011184c-54ec-4243-a42b-af2a6eee5552
begin
	r_squared_reprod_vf_nd, rms_values_reprod_vf_nd, fitted_line_reprod_vf_nd, coefficient_reprod_vf_nd = prepare_linear_regression(df_vf_r_large_nd, df_vf_r_medium_nd, df_vf_r_small_nd, df_vf_large_nd, df_vf_medium_nd, df_vf_small_nd)
	
	r_squared_reprod_vf_nd = round.(r_squared_reprod_vf_nd; digits=3)
	rms_values_reprod_vf_nd = round.(rms_values_reprod_vf_nd; digits=3)
	coefficient_reprod_vf_nd = round.(coefficient_reprod_vf_nd; digits=3)
end

# ╔═╡ 8ac18600-6073-4f2f-8709-43101f264a22
md"""
##### Agatston
"""

# ╔═╡ 9c4fc2f5-3b58-4f8a-a681-e246ac757dc2
array_a_nd = hcat(
	df_a_normal[!, :calculated_mass_large], 
	df_a_normal[!, :calculated_mass_medium], 
	df_a_normal[!, :calculated_mass_small]
);

# ╔═╡ cb402344-d304-4ec8-9fe4-2bda0a0d7bfd
array_a_r_nd = hcat(
    df_a_normal_r[!, :calculated_mass_large],
    df_a_normal_r[!, :calculated_mass_medium],
    df_a_normal_r[!, :calculated_mass_small],
);

# ╔═╡ f40d6a0e-384d-4d61-af93-e7398046d364
df_a_large_nd, df_a_r_large_nd, df_a_medium_nd, df_a_r_medium_nd, df_a_small_nd, df_a_r_small_nd = remove_false_negatives(df_a_normal, array_a_nd, df_a_normal_r, array_a_r_nd);

# ╔═╡ 022674ae-64de-4e3a-b70d-89ef428cc780
begin
	r_squared_reprod_a_nd, rms_values_reprod_a_nd, fitted_line_reprod_a_nd, coefficient_reprod_a_nd = prepare_linear_regression(df_a_r_large_nd, df_a_r_medium_nd, df_a_r_small_nd, df_a_large_nd, df_a_medium_nd, df_a_small_nd)
	
	r_squared_reprod_a_nd = round.(r_squared_reprod_a_nd; digits=3)
	rms_values_reprod_a_nd = round.(rms_values_reprod_a_nd; digits=3)
	coefficient_reprod_a_nd = round.(coefficient_reprod_a_nd; digits=3)
end

# ╔═╡ 16812354-9c81-48c7-83c0-97c0c7c6963f
md"""
##### Spatially Weighted
"""

# ╔═╡ dfe97497-f0b5-4cce-a6b8-b68382a1142a
array_s_nd = hcat(
	df_s_normal[!, :calculated_swcs_large], 
	df_s_normal[!, :calculated_swcs_medium], 
	df_s_normal[!, :calculated_swcs_small]
);

# ╔═╡ 832051a2-8290-49eb-be5a-cfba36f19238
array_s_r_nd = hcat(
    df_s_normal_r[!, :calculated_swcs_large],
    df_s_normal_r[!, :calculated_swcs_medium],
    df_s_normal_r[!, :calculated_swcs_small],
);

# ╔═╡ 9492ef06-4986-415b-9bec-cf3a84d686c6
df_s_large_nd, df_s_r_large_nd, df_s_medium_nd, df_s_r_medium_nd, df_s_small_nd, df_s_r_small_nd = remove_false_negatives(df_s_normal, array_s_nd, df_s_normal_r, array_s_r_nd);

# ╔═╡ d3f7af8b-4a3a-4d1d-81ff-e529bd8994c5
begin
	r_squared_reprod_s_nd, rms_values_reprod_s_nd, fitted_line_reprod_s_nd, coefficient_reprod_s_nd = prepare_linear_regression(df_s_r_large_nd, df_s_r_medium_nd, df_s_r_small_nd, df_s_large_nd, df_s_medium_nd, df_s_small_nd)
	
	r_squared_reprod_s_nd = round.(r_squared_reprod_s_nd; digits=3)
	rms_values_reprod_s_nd = round.(rms_values_reprod_s_nd; digits=3)
	coefficient_reprod_s_nd = round.(coefficient_reprod_s_nd; digits=3)
end

# ╔═╡ b7c5e3f3-d0a5-4501-b260-73960c00f160
function reprod_nd()
    f = Figure()

    ##-- A --##
    ax = Axis(
		f[1, 1],
	    xticks = [0, 50, 100, 150, 200],
	    yticks = [0, 50, 100, 150, 200],
	    xlabel = "Mass 1 (mg)",
	    ylabel = "Mass 2 (mg)",
	    title = "Integrated (Normal Density)",
		)
    scatter!(
		df_i_normal_r[!, :calculated_mass_large], 
		df_i_normal[!, :calculated_mass_large]
	)
    scatter!(
		df_i_normal_r[!, :calculated_mass_medium], 
		df_i_normal[!, :calculated_mass_medium]
	)
    scatter!(
        df_i_normal_r[!, :calculated_mass_small],
        df_i_normal[!, :calculated_mass_small];
        color=:red,
    )
    lines!([-1000, 1000], [-1000, 1000]; label="Unity")
    lines!(collect(1:1000), fitted_line_reprod_i_nd; linestyle=:dashdot)
	create_textbox(f[1, 1], coefficient_reprod_i_nd, r_squared_reprod_i_nd, rms_values_reprod_i_nd)

    xlims!(ax; low=0, high=200)
    ylims!(ax; low=0, high=200)

	##-- B --##
    ax = Axis(
		f[2, 1],
		xticks = [0, 50, 100, 150, 200],
    	yticks = [0, 50, 100, 150, 200],
    	xlabel = "Mass 1 (mg)",
    	ylabel = "Mass 2 (mg)",
    	title = "Volume Fraction (Normal Density)",
	)
    scatter!(
        df_vf_normal_r[!, :calculated_mass_large],
        df_vf_normal[!, :calculated_mass_large];
        label="Large Inserts",
    )
    scatter!(
        df_vf_normal_r[!, :calculated_mass_medium],
        df_vf_normal[!, :calculated_mass_medium];
        label="Medium Inserts",
    )
    scatter!(
        df_vf_normal_r[!, :calculated_mass_small],
        df_vf_normal[!, :calculated_mass_small];
        label="Small Inserts",
        color=:red,
    )
    lines!([-1000, 1000], [-1000, 1000]; label="Unity")
    lines!(collect(1:1000), fitted_line_reprod_vf_nd; linestyle=:dashdot)
	create_textbox(f[2, 1], coefficient_reprod_vf_nd, r_squared_reprod_vf_nd, rms_values_reprod_vf_nd)

    xlims!(ax; low=0, high=200)
    ylims!(ax; low=0, high=200)


    ##-- C --##
    ax = Axis(
		f[1, 2],
		xticks = [0, 50, 100, 150, 200],
    	yticks = [0, 50, 100, 150, 200],
    	xlabel = "Mass 1 (mg)",
    	ylabel = "Mass 2 (mg)",
    	title = "Agatston (Normal Density)",
	)
    scatter!(
        df_a_normal_r[!, :calculated_mass_large],
        df_a_normal[!, :calculated_mass_large];
        label="Large Inserts",
    )
    scatter!(
        df_a_normal_r[!, :calculated_mass_medium],
        df_a_normal[!, :calculated_mass_medium];
        label="Medium Inserts",
    )
    scatter!(
        df_a_normal_r[!, :calculated_mass_small],
        df_a_normal[!, :calculated_mass_small];
        label="Small Inserts",
        color=:red,
    )
    lines!([-1000, 1000], [-1000, 1000]; label="Unity")
    lines!(collect(1:1000), fitted_line_reprod_a_nd; linestyle=:dashdot)
	create_textbox(f[1, 2], coefficient_reprod_a_nd, r_squared_reprod_a_nd, rms_values_reprod_a_nd)
	
    xlims!(ax; low=0, high=200)
    ylims!(ax; low=0, high=200)


    # ##-- D --##
    ax = Axis(
		f[2, 2],
		xticks = [0, 125, 250, 375, 500],
    	yticks = [0, 125, 250, 375, 500],
    	xlabel = "SWCS 1",
    	ylabel = "SWCS 2",
    	title = "Spatially Weighted (Normal Density)",
	)
    scatter!(
        df_s_normal_r[!, :calculated_swcs_large],
        df_s_normal[!, :calculated_swcs_large];
        label="Large Inserts",
    )
    scatter!(
        df_s_normal_r[!, :calculated_swcs_medium],
        df_s_normal[!, :calculated_swcs_medium];
        label="Medium Inserts",
    )
    scatter!(
        df_s_normal_r[!, :calculated_swcs_small],
        df_s_normal[!, :calculated_swcs_small];
        label="Small Inserts",
        color=:red,
    )
    lines!([-1000, 1000], [-1000, 1000]; label="Unity")
    lines!(collect(1:1000), fitted_line_reprod_s_nd; linestyle=:dashdot, label="Fitted Line")
	create_textbox(f[2, 2], coefficient_reprod_s_nd, r_squared_reprod_s_nd, rms_values_reprod_s_nd)
	

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

    # save(joinpath(dirname(pwd()),"figures", FIGURE_PATH, "reproducibility.png"), f)
    return f
end

# ╔═╡ 8af3d98c-f33c-4cb2-b388-387bf2385530
with_theme(reprod_nd, medphys_theme)

# ╔═╡ 0c5c5a7b-3aa4-4193-b54a-7c753b819a60
md"""
#### Low Density
"""

# ╔═╡ 1eea3eb6-8006-4a9f-adc2-00e19aec2af9
md"""
##### Integrated
"""

# ╔═╡ a72aed45-966b-48fd-a1f5-ae39c1f1ca05
array_i_ld = hcat(
	df_i_low[!, :calculated_mass_large], 
	df_i_low[!, :calculated_mass_medium], 
	df_i_low[!, :calculated_mass_small]
);

# ╔═╡ 6e5c06d3-7b5e-4aae-985d-e7458a3cb263
array_i_r_ld = hcat(
    df_i_low_r[!, :calculated_mass_large],
    df_i_low_r[!, :calculated_mass_medium],
    df_i_low_r[!, :calculated_mass_small],
);

# ╔═╡ 8c68de3d-8126-47e4-bdbd-982ed0244b65
df_i_large_ld, df_i_r_large_ld, df_i_medium_ld, df_i_r_medium_ld, df_i_small_ld, df_i_r_small_ld = remove_false_negatives(df_i_low, array_i_ld, df_i_low_r, array_i_r_ld);

# ╔═╡ 3863e320-368c-4816-82d1-eeccc8ed004d
begin
	r_squared_reprod_i_ld, rms_values_reprod_i_ld, fitted_line_reprod_i_ld, coefficient_reprod_i_ld = prepare_linear_regression(df_i_r_large_ld, df_i_r_medium_ld, df_i_r_small_ld, df_i_large_ld, df_i_medium_ld, df_i_small_ld)
	
	r_squared_reprod_i_ld = round.(r_squared_reprod_i_ld; digits=3)
	rms_values_reprod_i_ld = round.(rms_values_reprod_i_ld; digits=3)
	coefficient_reprod_i_ld = round.(coefficient_reprod_i_ld; digits=3)
end

# ╔═╡ d3278fe2-e57b-4983-a6eb-8dead9a70c94
md"""
##### Volume Fraction
"""

# ╔═╡ d205dd2f-b1a2-4285-ae27-31ac9cb5bfe2
array_vf_ld = hcat(
	df_vf_low[!, :calculated_mass_large], 
	df_vf_low[!, :calculated_mass_medium], 
	df_vf_low[!, :calculated_mass_small]
);

# ╔═╡ 967baf16-5b6b-4611-9703-30990bf05852
array_vf_r_ld = hcat(
    df_vf_low_r[!, :calculated_mass_large],
    df_vf_low_r[!, :calculated_mass_medium],
    df_vf_low_r[!, :calculated_mass_small],
);

# ╔═╡ b83a626e-d45d-40ba-89ef-4787d44332ad
df_vf_large_ld, df_vf_r_large_ld, df_vf_medium_ld, df_vf_r_medium_ld, df_vf_small_ld, df_vf_r_small_ld = remove_false_negatives(df_vf_low, array_vf_ld, df_vf_low_r, array_vf_r_ld);

# ╔═╡ 7969e05f-5148-4d1b-ac77-f8b9dde09e64
begin
	r_squared_reprod_vf_ld, rms_values_reprod_vf_ld, fitted_line_reprod_vf_ld, coefficient_reprod_vf_ld = prepare_linear_regression(df_vf_r_large_ld, df_vf_r_medium_ld, df_vf_r_small_ld, df_vf_large_ld, df_vf_medium_ld, df_vf_small_ld)
	
	r_squared_reprod_vf_ld = round.(r_squared_reprod_vf_ld; digits=3)
	rms_values_reprod_vf_ld = round.(rms_values_reprod_vf_ld; digits=3)
	coefficient_reprod_vf_ld = round.(coefficient_reprod_vf_ld; digits=3)
end

# ╔═╡ ef9255e1-3300-4354-8c9f-6226e3f32145
md"""
##### Agatston
"""

# ╔═╡ c77c6d5c-a989-42fd-99ef-4adbf6aa6db5
array_a_ld = hcat(
	df_a_low[!, :calculated_mass_large], 
	df_a_low[!, :calculated_mass_medium], 
	df_a_low[!, :calculated_mass_small]
);

# ╔═╡ adcc69bc-c067-49cc-a017-5dc11ad92063
array_a_r_ld = hcat(
    df_a_low_r[!, :calculated_mass_large],
    df_a_low_r[!, :calculated_mass_medium],
    df_a_low_r[!, :calculated_mass_small],
);

# ╔═╡ ed7cad03-0e74-441f-b5d2-2426c0cf5767
df_a_large_ld, df_a_r_large_ld, df_a_medium_ld, df_a_r_medium_ld, df_a_small_ld, df_a_r_small_ld = remove_false_negatives(df_a_low, array_a_ld, df_a_low_r, array_a_r_ld);

# ╔═╡ 27abe3ed-eaa6-4acf-b83f-56a0b7838501
begin
	r_squared_reprod_a_ld, rms_values_reprod_a_ld, fitted_line_reprod_a_ld, coefficient_reprod_a_ld = prepare_linear_regression(df_a_r_large_ld, df_a_r_medium_ld, df_a_r_small_ld, df_a_large_ld, df_a_medium_ld, df_a_small_ld)
	
	r_squared_reprod_a_ld = round.(r_squared_reprod_a_ld; digits=3)
	rms_values_reprod_a_ld = round.(rms_values_reprod_a_ld; digits=3)
	coefficient_reprod_a_ld = round.(coefficient_reprod_a_ld; digits=3)
end

# ╔═╡ 8b39efb6-b0dc-418c-9e1f-7bf7231db159
md"""
##### Spatially Weighted
"""

# ╔═╡ a2d4a55c-da21-43b8-aec3-d80b5e29bc04
array_s_ld = hcat(
	df_s_low[!, :calculated_swcs_large], 
	df_s_low[!, :calculated_swcs_medium], 
	df_s_low[!, :calculated_swcs_small]
);

# ╔═╡ 919efe63-7c77-45b8-b8d6-ed7f6e0bf694
array_s_r_ld = hcat(
    df_s_low_r[!, :calculated_swcs_large],
    df_s_low_r[!, :calculated_swcs_medium],
    df_s_low_r[!, :calculated_swcs_small],
);

# ╔═╡ 51bc03d9-6821-47cd-9435-1d79aa37f779
df_s_large_ld, df_s_r_large_ld, df_s_medium_ld, df_s_r_medium_ld, df_s_small_ld, df_s_r_small_ld = remove_false_negatives(df_s_low, array_s_ld, df_s_low_r, array_s_r_ld);

# ╔═╡ 9848ed54-530c-46c9-b47c-5214c401c6b2
begin
	r_squared_reprod_s_ld, rms_values_reprod_s_ld, fitted_line_reprod_s_ld, coefficient_reprod_s_ld = prepare_linear_regression(df_s_r_large_ld, df_s_r_medium_ld, df_s_r_small_ld, df_s_large_ld, df_s_medium_ld, df_s_small_ld)
	
	r_squared_reprod_s_ld = round.(r_squared_reprod_s_ld; digits=3)
	rms_values_reprod_s_ld = round.(rms_values_reprod_s_ld; digits=3)
	coefficient_reprod_s_ld = round.(coefficient_reprod_s_ld; digits=3)
end

# ╔═╡ 3d7de16f-5a40-4737-af20-7bdd535239d7
function reprod_ld()
    f = Figure()

    ##-- A --##
    ax = Axis(
		f[1, 1],
		xticks = [0, 5, 10, 15, 20, 25],
		yticks = [0, 15, 30],
	    xlabel = "Mass 1 (mg)",
	    ylabel = "Mass 2 (mg)",
	    title = "Integrated (Low Density)",
		)
    scatter!(
		df_i_low_r[!, :calculated_mass_large], 
		df_i_low[!, :calculated_mass_large]
	)
    scatter!(
		df_i_low_r[!, :calculated_mass_medium], 
		df_i_low[!, :calculated_mass_medium]
	)
    scatter!(
        df_i_low_r[!, :calculated_mass_small],
        df_i_low[!, :calculated_mass_small];
        color=:red,
    )
    lines!([-1000, 1000], [-1000, 1000]; label="Unity")
    lines!(collect(1:1000), fitted_line_reprod_i_ld; linestyle=:dashdot)
	create_textbox(f[1, 1], coefficient_reprod_i_ld, r_squared_reprod_i_ld, rms_values_reprod_i_ld)

    xlims!(ax; low=0, high=25)
    ylims!(ax; low=-5, high=30)

	##-- B --##
    ax = Axis(
		f[2, 1],
		xticks = [0, 5, 10, 15, 20, 25],
		yticks = [0, 15, 30],
    	xlabel = "Mass 1 (mg)",
    	ylabel = "Mass 2 (mg)",
    	title = "Volume Fraction (Low Density)",
	)
    scatter!(
        df_vf_low_r[!, :calculated_mass_large],
        df_vf_low[!, :calculated_mass_large];
        label="Large Inserts",
    )
    scatter!(
        df_vf_low_r[!, :calculated_mass_medium],
        df_vf_low[!, :calculated_mass_medium];
        label="Medium Inserts",
    )
    scatter!(
        df_vf_low_r[!, :calculated_mass_small],
        df_vf_low[!, :calculated_mass_small];
        label="Small Inserts",
        color=:red,
    )
    lines!([-1000, 1000], [-1000, 1000]; label="Unity")
    lines!(collect(1:1000), fitted_line_reprod_vf_ld; linestyle=:dashdot)
	create_textbox(f[2, 1], coefficient_reprod_vf_ld, r_squared_reprod_vf_ld, rms_values_reprod_vf_ld)

    xlims!(ax; low=0, high=25)
    ylims!(ax; low=-5, high=30)


    ##-- C --##
    ax = Axis(
		f[1, 2],
		xticks = [0, 5, 10, 15, 20, 25],
		yticks = [0, 15, 30],
    	xlabel = "Mass 1 (mg)",
    	ylabel = "Mass 2 (mg)",
    	title = "Agatston (Low Density)",
	)
    scatter!(
        df_a_low_r[!, :calculated_mass_large],
        df_a_low[!, :calculated_mass_large];
        label="Large Inserts",
    )
    scatter!(
        df_a_low_r[!, :calculated_mass_medium],
        df_a_low[!, :calculated_mass_medium];
        label="Medium Inserts",
    )
    scatter!(
        df_a_low_r[!, :calculated_mass_small],
        df_a_low[!, :calculated_mass_small];
        label="Small Inserts",
        color=:red,
    )
    lines!([-1000, 1000], [-1000, 1000]; label="Unity")
    lines!(collect(1:1000), fitted_line_reprod_a_ld; linestyle=:dashdot)
	create_textbox(f[1, 2], coefficient_reprod_a_ld, r_squared_reprod_a_ld, rms_values_reprod_a_ld)
	
    xlims!(ax; low=0, high=25)
    ylims!(ax; low=-5, high=30)


    # ##-- D --##
    ax = Axis(
		f[2, 2],
		xticks = [0, 100, 200, 300],
    	yticks = [0, 100, 200, 300],
    	xlabel = "SWCS 1",
    	ylabel = "SWCS 2",
    	title = "Spatially Weighted (Low Density)",
	)
    scatter!(
        df_s_low_r[!, :calculated_swcs_large],
        df_s_low[!, :calculated_swcs_large];
        label="Large Inserts",
    )
    scatter!(
        df_s_low_r[!, :calculated_swcs_medium],
        df_s_low[!, :calculated_swcs_medium];
        label="Medium Inserts",
    )
    scatter!(
        df_s_low_r[!, :calculated_swcs_small],
        df_s_low[!, :calculated_swcs_small];
        label="Small Inserts",
        color=:red,
    )
    lines!([-1000, 1000], [-1000, 1000]; label="Unity")
    lines!(collect(1:1000), fitted_line_reprod_s_ld; linestyle=:dashdot, label="Fitted Line")
	create_textbox(f[2, 2], coefficient_reprod_s_ld, r_squared_reprod_s_ld, rms_values_reprod_s_ld)
	

    xlims!(ax; low=0, high=300)
    ylims!(ax; low=0, high=300)


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

    # save(joinpath(dirname(pwd()),"figures", FIGURE_PATH, "reproducibility.png"), f)
    return f
end

# ╔═╡ b3c5c174-24be-447f-90f7-31a6252a55ae
with_theme(reprod_ld, medphys_theme)

# ╔═╡ 87c5676f-3ed2-49fd-ae79-aee072b2bc80
md"""
### Patient Size
"""

# ╔═╡ 5b853ca7-f2e0-4d38-87b8-e48faa9102c2
begin 
	df_i_small_patient_r, df_i_medium_patient_r, df_i_large_patient_r = groupby(df_i_r, :SIZE)
	df_vf_small_patient_r, df_vf_medium_patient_r, df_vf_large_patient_r = groupby(df_vf_r, :SIZE)
	df_a_small_patient_r, df_a_medium_patient_r, df_a_large_patient_r = groupby(df_a_r, :SIZE)
	df_s_small_patient_r, df_s_medium_patient_r, df_s_large_patient_r = groupby(df_s_r, :SIZE)
end;

# ╔═╡ bd05fce5-a7d7-4e85-afa3-995224183406
md"""
#### Small Patient
"""

# ╔═╡ b3cc1fa7-c6d4-4002-bada-ea882f740e8b
md"""
##### Integrated
"""

# ╔═╡ a82a6cf4-9599-419b-8ffc-56dda90954fb
array_i_sp = hcat(
	df_i_small_patient[!, :calculated_mass_large], 
	df_i_small_patient[!, :calculated_mass_medium], 
	df_i_small_patient[!, :calculated_mass_small]
);

# ╔═╡ f8d840c1-ec0b-4165-ae75-79cd1e0a4fd6
array_i_r_sp = hcat(
    df_i_small_patient_r[!, :calculated_mass_large],
    df_i_small_patient_r[!, :calculated_mass_medium],
    df_i_small_patient_r[!, :calculated_mass_small],
);

# ╔═╡ 38605356-7165-44df-95c6-532b411ca708
df_i_large_sp, df_i_r_large_sp, df_i_medium_sp, df_i_r_medium_sp, df_i_small_sp, df_i_r_small_sp = remove_false_negatives(df_i_small_patient, array_i_sp, df_i_small_patient_r, array_i_r_sp);

# ╔═╡ 79c52fd3-8b91-4982-8c5e-2b212c0d4297
begin
	r_squared_reprod_i_sp, rms_values_reprod_i_sp, fitted_line_reprod_i_sp, coefficient_reprod_i_sp = prepare_linear_regression(df_i_r_large_sp, df_i_r_medium_sp, df_i_r_small_sp, df_i_large_sp, df_i_medium_sp, df_i_small_sp)
	
	r_squared_reprod_i_sp = round.(r_squared_reprod_i_sp; digits=3)
	rms_values_reprod_i_sp = round.(rms_values_reprod_i_sp; digits=3)
	coefficient_reprod_i_sp = round.(coefficient_reprod_i_sp; digits=3)
end

# ╔═╡ d36e7476-a540-46a0-86de-88ad0c78036d
md"""
##### Volume Fraction
"""

# ╔═╡ 76c1dc7e-e1fd-481a-bdff-f3b17b32f584
array_vf_sp = hcat(
	df_vf_small_patient[!, :calculated_mass_large], 
	df_vf_small_patient[!, :calculated_mass_medium], 
	df_vf_small_patient[!, :calculated_mass_small]
);

# ╔═╡ 96f9e011-c879-4f39-801e-ca7469c6ef0a
array_vf_r_sp = hcat(
    df_vf_small_patient_r[!, :calculated_mass_large],
    df_vf_small_patient_r[!, :calculated_mass_medium],
    df_vf_small_patient_r[!, :calculated_mass_small],
);

# ╔═╡ a9a49218-1239-4caf-9774-3a5614ab40be
df_vf_large_sp, df_vf_r_large_sp, df_vf_medium_sp, df_vf_r_medium_sp, df_vf_small_sp, df_vf_r_small_sp = remove_false_negatives(df_vf_small_patient, array_vf_sp, df_vf_small_patient_r, array_vf_r_sp);

# ╔═╡ 66ad4a90-27bd-4857-b0be-09b6fa25741a
begin
	r_squared_reprod_vf_sp, rms_values_reprod_vf_sp, fitted_line_reprod_vf_sp, coefficient_reprod_vf_sp = prepare_linear_regression(df_vf_r_large_sp, df_vf_r_medium_sp, df_vf_r_small_sp, df_vf_large_sp, df_vf_medium_sp, df_vf_small_sp)
	
	r_squared_reprod_vf_sp = round.(r_squared_reprod_vf_sp; digits=3)
	rms_values_reprod_vf_sp = round.(rms_values_reprod_vf_sp; digits=3)
	coefficient_reprod_vf_sp = round.(coefficient_reprod_vf_sp; digits=3)
end

# ╔═╡ 4217ba04-b8b9-425f-8991-10a058621478
md"""
##### Agatston
"""

# ╔═╡ b641aca4-af92-4110-8cd4-2ba4c122df04
array_a_sp = hcat(
	df_a_small_patient[!, :calculated_mass_large], 
	df_a_small_patient[!, :calculated_mass_medium], 
	df_a_small_patient[!, :calculated_mass_small]
);

# ╔═╡ eb4d0f50-65e4-44ba-9eb6-88f84ed4b070
array_a_r_sp = hcat(
    df_a_small_patient_r[!, :calculated_mass_large],
    df_a_small_patient_r[!, :calculated_mass_medium],
    df_a_small_patient_r[!, :calculated_mass_small],
);

# ╔═╡ ba357e75-bcc8-47db-9a11-6baf050ae925
df_a_large_sp, df_a_r_large_sp, df_a_medium_sp, df_a_r_medium_sp, df_a_small_sp, df_a_r_small_sp = remove_false_negatives(df_a_small_patient, array_a_sp, df_a_small_patient_r, array_a_r_sp);

# ╔═╡ c767a4a8-a521-4d68-990b-50c177788193
begin
	r_squared_reprod_a_sp, rms_values_reprod_a_sp, fitted_line_reprod_a_sp, coefficient_reprod_a_sp = prepare_linear_regression(df_a_r_large_sp, df_a_r_medium_sp, df_a_r_small_sp, df_a_large_sp, df_a_medium_sp, df_a_small_sp)
	
	r_squared_reprod_a_sp = round.(r_squared_reprod_a_sp; digits=3)
	rms_values_reprod_a_sp = round.(rms_values_reprod_a_sp; digits=3)
	coefficient_reprod_a_sp = round.(coefficient_reprod_a_sp; digits=3)
end

# ╔═╡ b4fe4133-256f-4e33-b4e5-6fa38a126c4f
md"""
##### Spatially Weighted
"""

# ╔═╡ 74ea613c-e4c4-4124-ae48-8ec6e1df68ff
array_s_sp = hcat(
	df_s_small_patient[!, :calculated_swcs_large], 
	df_s_small_patient[!, :calculated_swcs_medium], 
	df_s_small_patient[!, :calculated_swcs_small]
);

# ╔═╡ 6c7cdb74-cae8-4842-872c-65afc4ddd3a6
array_s_r_sp = hcat(
    df_s_small_patient_r[!, :calculated_swcs_large],
    df_s_small_patient_r[!, :calculated_swcs_medium],
    df_s_small_patient_r[!, :calculated_swcs_small],
);

# ╔═╡ d5b1182a-fbdb-4422-abda-6df0a5490003
df_s_large_sp, df_s_r_large_sp, df_s_medium_sp, df_s_r_medium_sp, df_s_small_sp, df_s_r_small_sp = remove_false_negatives(df_s_small_patient, array_s_sp, df_s_small_patient_r, array_s_r_sp);

# ╔═╡ 5f2a36aa-3b66-441d-a909-14f6d2820f2c
begin
	r_squared_reprod_s_sp, rms_values_reprod_s_sp, fitted_line_reprod_s_sp, coefficient_reprod_s_sp = prepare_linear_regression(df_s_r_large_sp, df_s_r_medium_sp, df_s_r_small_sp, df_s_large_sp, df_s_medium_sp, df_s_small_sp)
	
	r_squared_reprod_s_sp = round.(r_squared_reprod_s_sp; digits=3)
	rms_values_reprod_s_sp = round.(rms_values_reprod_s_sp; digits=3)
	coefficient_reprod_s_sp = round.(coefficient_reprod_s_sp; digits=3)
end

# ╔═╡ f46a6aac-eaee-4016-b628-b31c02be1798
function reprod_sp()
    f = Figure()

    ##-- A --##
    ax = Axis(
		f[1, 1],
	    xticks = [0, 50, 100, 150, 200],
	    yticks = [0, 50, 100, 150, 200],
	    xlabel = "Mass 1 (mg)",
	    ylabel = "Mass 2 (mg)",
	    title = "Integrated (Small Patient)",
		)
    scatter!(
		df_i_small_patient_r[!, :calculated_mass_large], 
		df_i_small_patient[!, :calculated_mass_large]
	)
    scatter!(
		df_i_small_patient_r[!, :calculated_mass_medium], 
		df_i_small_patient[!, :calculated_mass_medium]
	)
    scatter!(
        df_i_small_patient_r[!, :calculated_mass_small],
        df_i_small_patient[!, :calculated_mass_small];
        color=:red,
    )
    lines!([-1000, 1000], [-1000, 1000]; label="Unity")
    lines!(collect(1:1000), fitted_line_reprod_i_sp; linestyle=:dashdot)
	create_textbox(f[1, 1], coefficient_reprod_i_sp, r_squared_reprod_i_sp, rms_values_reprod_i_sp)

    xlims!(ax; low=0, high=200)
    ylims!(ax; low=0, high=200)

	##-- B --##
    ax = Axis(
		f[2, 1],
		xticks = [0, 50, 100, 150, 200],
    	yticks = [0, 50, 100, 150, 200],
    	xlabel = "Mass 1 (mg)",
    	ylabel = "Mass 2 (mg)",
    	title = "Volume Fraction (Small Patient)",
	)
    scatter!(
        df_vf_large_patient_r[!, :calculated_mass_large],
        df_vf_large_patient[!, :calculated_mass_large];
        label="Large Inserts",
    )
    scatter!(
        df_vf_medium_patient_r[!, :calculated_mass_medium],
        df_vf_medium_patient[!, :calculated_mass_medium];
        label="Medium Inserts",
    )
    scatter!(
        df_vf_small_patient_r[!, :calculated_mass_small],
        df_vf_small_patient[!, :calculated_mass_small];
        label="Small Inserts",
        color=:red,
    )
    lines!([-1000, 1000], [-1000, 1000]; label="Unity")
    lines!(collect(1:1000), fitted_line_reprod_vf_sp; linestyle=:dashdot)
	create_textbox(f[2, 1], coefficient_reprod_vf_sp, r_squared_reprod_vf_sp, rms_values_reprod_vf_sp)

    xlims!(ax; low=0, high=200)
    ylims!(ax; low=0, high=200)


    ##-- C --##
    ax = Axis(
		f[1, 2],
		xticks = [0, 50, 100, 150, 200],
    	yticks = [0, 50, 100, 150, 200],
    	xlabel = "Mass 1 (mg)",
    	ylabel = "Mass 2 (mg)",
    	title = "Agatston (Small Patient)",
	)
    scatter!(
        df_a_large_patient_r[!, :calculated_mass_large],
        df_a_large_patient[!, :calculated_mass_large];
        label="Large Inserts",
    )
    scatter!(
        df_a_medium_patient_r[!, :calculated_mass_medium],
        df_a_medium_patient[!, :calculated_mass_medium];
        label="Medium Inserts",
    )
    scatter!(
        df_a_small_patient_r[!, :calculated_mass_small],
        df_a_small_patient[!, :calculated_mass_small];
        label="Small Inserts",
        color=:red,
    )
    lines!([-1000, 1000], [-1000, 1000]; label="Unity")
    lines!(collect(1:1000), fitted_line_reprod_a_sp; linestyle=:dashdot)
	create_textbox(f[1, 2], coefficient_reprod_a_sp, r_squared_reprod_a_sp, rms_values_reprod_a_sp)
	
    xlims!(ax; low=0, high=200)
    ylims!(ax; low=0, high=200)


    # ##-- D --##
    ax = Axis(
		f[2, 2],
		xticks = [0, 125, 250, 375, 500],
    	yticks = [0, 125, 250, 375, 500],
    	xlabel = "SWCS 1",
    	ylabel = "SWCS 2",
    	title = "Spatially Weighted (Small Patient)",
	)
    scatter!(
        df_s_large_patient_r[!, :calculated_swcs_large],
        df_s_large_patient[!, :calculated_swcs_large];
        label="Large Inserts",
    )
    scatter!(
        df_s_medium_patient_r[!, :calculated_swcs_medium],
        df_s_medium_patient[!, :calculated_swcs_medium];
        label="Medium Inserts",
    )
    scatter!(
        df_s_small_patient_r[!, :calculated_swcs_small],
        df_s_small_patient[!, :calculated_swcs_small];
        label="Small Inserts",
        color=:red,
    )
    lines!([-1000, 1000], [-1000, 1000]; label="Unity")
    lines!(collect(1:1000), fitted_line_reprod_s_sp; linestyle=:dashdot, label="Fitted Line")
	create_textbox(f[2, 2], coefficient_reprod_s_sp, r_squared_reprod_s_sp, rms_values_reprod_s_sp)
	

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

    # save(joinpath(dirname(pwd()),"figures", FIGURE_PATH, "reproducibility.png"), f)
    return f
end

# ╔═╡ 12869712-c5f0-4f31-8088-e885b5d9ddb9
with_theme(reprod_sp, medphys_theme)

# ╔═╡ aa0487c3-2fe8-49f4-bf4d-927ce73de052
md"""
#### Medium Patient
"""

# ╔═╡ 03bf4424-01f7-4d50-aae2-fc4c6c3a19ba
md"""
##### Integrated
"""

# ╔═╡ d4e37208-0eab-4091-b360-94f74c9bfbd7
array_i_mp = hcat(
	df_i_medium_patient[!, :calculated_mass_large], 
	df_i_medium_patient[!, :calculated_mass_medium], 
	df_i_medium_patient[!, :calculated_mass_small]
);

# ╔═╡ 9bd4a8d3-c641-4fb5-b0bb-f6d699e944c4
array_i_r_mp = hcat(
    df_i_medium_patient_r[!, :calculated_mass_large],
    df_i_medium_patient_r[!, :calculated_mass_medium],
    df_i_medium_patient_r[!, :calculated_mass_small],
);

# ╔═╡ a09c2653-21b8-4a17-8b02-2451cda66ac4
df_i_large_mp, df_i_r_large_mp, df_i_medium_mp, df_i_r_medium_mp, df_i_small_mp, df_i_r_small_mp = remove_false_negatives(df_i_medium_patient, array_i_mp, df_i_medium_patient_r, array_i_r_mp);

# ╔═╡ cbb80103-63ac-40c0-8eae-d24c058695a5
begin
	r_squared_reprod_i_mp, rms_values_reprod_i_mp, fitted_line_reprod_i_mp, coefficient_reprod_i_mp = prepare_linear_regression(df_i_r_large_mp, df_i_r_medium_mp, df_i_r_small_mp, df_i_large_mp, df_i_medium_mp, df_i_small_mp)
	
	r_squared_reprod_i_mp = round.(r_squared_reprod_i_mp; digits=3)
	rms_values_reprod_i_mp = round.(rms_values_reprod_i_mp; digits=3)
	coefficient_reprod_i_mp = round.(coefficient_reprod_i_mp; digits=3)
end

# ╔═╡ 6af52148-3473-4660-b6cc-aa9879ee46cb
md"""
##### Volume Fraction
"""

# ╔═╡ fc88c296-38f3-4ccf-a1d9-73bd40e26b40
array_vf_mp = hcat(
	df_vf_medium_patient[!, :calculated_mass_large], 
	df_vf_medium_patient[!, :calculated_mass_medium], 
	df_vf_medium_patient[!, :calculated_mass_small]
);

# ╔═╡ 836f13b7-4aca-4050-b6ee-4522a45738ec
array_vf_r_mp = hcat(
    df_vf_medium_patient_r[!, :calculated_mass_large],
    df_vf_medium_patient_r[!, :calculated_mass_medium],
    df_vf_medium_patient_r[!, :calculated_mass_small],
);

# ╔═╡ 4181e1d4-503d-4dfe-ae62-be3e50af3d9f
df_vf_large_mp, df_vf_r_large_mp, df_vf_medium_mp, df_vf_r_medium_mp, df_vf_small_mp, df_vf_r_small_mp = remove_false_negatives(df_vf_medium_patient, array_vf_mp, df_vf_medium_patient_r, array_vf_r_mp);

# ╔═╡ e019181d-639b-474d-bd48-9e714502cd23
begin
	r_squared_reprod_vf_mp, rms_values_reprod_vf_mp, fitted_line_reprod_vf_mp, coefficient_reprod_vf_mp = prepare_linear_regression(df_vf_r_large_mp, df_vf_r_medium_mp, df_vf_r_small_mp, df_vf_large_mp, df_vf_medium_mp, df_vf_small_mp)
	
	r_squared_reprod_vf_mp = round.(r_squared_reprod_vf_mp; digits=3)
	rms_values_reprod_vf_mp = round.(rms_values_reprod_vf_mp; digits=3)
	coefficient_reprod_vf_mp = round.(coefficient_reprod_vf_mp; digits=3)
end

# ╔═╡ 889c5c63-146f-4b77-b067-3d4c6f977bc5
md"""
##### Agatston
"""

# ╔═╡ ed14543e-8444-4c00-aa32-e05d5f34bcf3
array_a_mp = hcat(
	df_a_medium_patient[!, :calculated_mass_large], 
	df_a_medium_patient[!, :calculated_mass_medium], 
	df_a_medium_patient[!, :calculated_mass_small]
);

# ╔═╡ 77ecb385-24c3-45b9-857b-3c3a1d5480a1
array_a_r_mp = hcat(
    df_a_medium_patient_r[!, :calculated_mass_large],
    df_a_medium_patient_r[!, :calculated_mass_medium],
    df_a_medium_patient_r[!, :calculated_mass_small],
);

# ╔═╡ 9358b9c7-5dfa-4dda-ab09-e522c5fcb202
df_a_large_mp, df_a_r_large_mp, df_a_medium_mp, df_a_r_medium_mp, df_a_small_mp, df_a_r_small_mp = remove_false_negatives(df_a_medium_patient, array_a_mp, df_a_medium_patient_r, array_a_r_mp);

# ╔═╡ 6eda4c2d-b56c-4366-a4f3-0e0dd710e061
begin
	r_squared_reprod_a_mp, rms_values_reprod_a_mp, fitted_line_reprod_a_mp, coefficient_reprod_a_mp = prepare_linear_regression(df_a_r_large_mp, df_a_r_medium_mp, df_a_r_small_mp, df_a_large_mp, df_a_medium_mp, df_a_small_mp)
	
	r_squared_reprod_a_mp = round.(r_squared_reprod_a_mp; digits=3)
	rms_values_reprod_a_mp = round.(rms_values_reprod_a_mp; digits=3)
	coefficient_reprod_a_mp = round.(coefficient_reprod_a_mp; digits=3)
end

# ╔═╡ a0da5d2f-835e-4250-8841-2d22b6531df0
md"""
##### Spatially Weighted
"""

# ╔═╡ 10dda852-aa33-49a6-adfc-08b95cfe70b5
array_s_mp = hcat(
	df_s_medium_patient[!, :calculated_swcs_large], 
	df_s_medium_patient[!, :calculated_swcs_medium], 
	df_s_medium_patient[!, :calculated_swcs_small]
);

# ╔═╡ 1e7b01d3-2077-4f22-9c14-02231a5bb8cb
array_s_r_mp = hcat(
    df_s_medium_patient_r[!, :calculated_swcs_large],
    df_s_medium_patient_r[!, :calculated_swcs_medium],
    df_s_medium_patient_r[!, :calculated_swcs_small],
);

# ╔═╡ 4e0b3ef7-5a75-46da-9797-94ad80005bec
df_s_large_mp, df_s_r_large_mp, df_s_medium_mp, df_s_r_medium_mp, df_s_small_mp, df_s_r_small_mp = remove_false_negatives(df_s_medium_patient, array_s_mp, df_s_medium_patient_r, array_s_r_mp);

# ╔═╡ 39260980-b87b-442e-82a9-532f3854b607
begin
	r_squared_reprod_s_mp, rms_values_reprod_s_mp, fitted_line_reprod_s_mp, coefficient_reprod_s_mp = prepare_linear_regression(df_s_r_large_mp, df_s_r_medium_mp, df_s_r_small_mp, df_s_large_mp, df_s_medium_mp, df_s_small_mp)
	
	r_squared_reprod_s_mp = round.(r_squared_reprod_s_mp; digits=3)
	rms_values_reprod_s_mp = round.(rms_values_reprod_s_mp; digits=3)
	coefficient_reprod_s_mp = round.(coefficient_reprod_s_mp; digits=3)
end

# ╔═╡ 457f79fe-4884-4bb8-b6ef-6da088cefc76
function reprod_mp()
    f = Figure()

    ##-- A --##
    ax = Axis(
		f[1, 1],
	    xticks = [0, 50, 100, 150, 200],
	    yticks = [0, 50, 100, 150, 200],
	    xlabel = "Mass 1 (mg)",
	    ylabel = "Mass 2 (mg)",
	    title = "Integrated (Med Patient)",
		)
    scatter!(
		df_i_medium_patient_r[!, :calculated_mass_large], 
		df_i_medium_patient[!, :calculated_mass_large]
	)
    scatter!(
		df_i_medium_patient_r[!, :calculated_mass_medium], 
		df_i_medium_patient[!, :calculated_mass_medium]
	)
    scatter!(
        df_i_medium_patient_r[!, :calculated_mass_small],
        df_i_medium_patient[!, :calculated_mass_small];
        color=:red,
    )
    lines!([-1000, 1000], [-1000, 1000]; label="Unity")
    lines!(collect(1:1000), fitted_line_reprod_i_mp; linestyle=:dashdot)
	create_textbox(f[1, 1], coefficient_reprod_i_mp, r_squared_reprod_i_mp, rms_values_reprod_i_mp)

    xlims!(ax; low=0, high=200)
    ylims!(ax; low=0, high=200)

	##-- B --##
    ax = Axis(
		f[2, 1],
		xticks = [0, 50, 100, 150, 200],
    	yticks = [0, 50, 100, 150, 200],
    	xlabel = "Mass 1 (mg)",
    	ylabel = "Mass 2 (mg)",
    	title = "Volume Fraction (Med Patient)",
	)
    scatter!(
        df_vf_large_patient_r[!, :calculated_mass_large],
        df_vf_large_patient[!, :calculated_mass_large];
        label="Large Inserts",
    )
    scatter!(
        df_vf_medium_patient_r[!, :calculated_mass_medium],
        df_vf_medium_patient[!, :calculated_mass_medium];
        label="Medium Inserts",
    )
    scatter!(
        df_vf_medium_patient_r[!, :calculated_mass_small],
        df_vf_medium_patient[!, :calculated_mass_small];
        label="Small Inserts",
        color=:red,
    )
    lines!([-1000, 1000], [-1000, 1000]; label="Unity")
    lines!(collect(1:1000), fitted_line_reprod_vf_mp; linestyle=:dashdot)
	create_textbox(f[2, 1], coefficient_reprod_vf_mp, r_squared_reprod_vf_mp, rms_values_reprod_vf_mp)

    xlims!(ax; low=0, high=200)
    ylims!(ax; low=0, high=200)


    ##-- C --##
    ax = Axis(
		f[1, 2],
		xticks = [0, 50, 100, 150, 200],
    	yticks = [0, 50, 100, 150, 200],
    	xlabel = "Mass 1 (mg)",
    	ylabel = "Mass 2 (mg)",
    	title = "Agatston (Med Patient)",
	)
    scatter!(
        df_a_large_patient_r[!, :calculated_mass_large],
        df_a_large_patient[!, :calculated_mass_large];
        label="Large Inserts",
    )
    scatter!(
        df_a_medium_patient_r[!, :calculated_mass_medium],
        df_a_medium_patient[!, :calculated_mass_medium];
        label="Medium Inserts",
    )
    scatter!(
        df_a_medium_patient_r[!, :calculated_mass_small],
        df_a_medium_patient[!, :calculated_mass_small];
        label="Small Inserts",
        color=:red,
    )
    lines!([-1000, 1000], [-1000, 1000]; label="Unity")
    lines!(collect(1:1000), fitted_line_reprod_a_mp; linestyle=:dashdot)
	create_textbox(f[1, 2], coefficient_reprod_a_mp, r_squared_reprod_a_mp, rms_values_reprod_a_mp)
	
    xlims!(ax; low=0, high=200)
    ylims!(ax; low=0, high=200)


    # ##-- D --##
    ax = Axis(
		f[2, 2],
		xticks = [0, 125, 250, 375, 500],
    	yticks = [0, 125, 250, 375, 500],
    	xlabel = "SWCS 1",
    	ylabel = "SWCS 2",
    	title = "Spatially Weighted (Med Patient)",
	)
    scatter!(
        df_s_large_patient_r[!, :calculated_swcs_large],
        df_s_large_patient[!, :calculated_swcs_large];
        label="Large Inserts",
    )
    scatter!(
        df_s_medium_patient_r[!, :calculated_swcs_medium],
        df_s_medium_patient[!, :calculated_swcs_medium];
        label="Medium Inserts",
    )
    scatter!(
        df_s_medium_patient_r[!, :calculated_swcs_small],
        df_s_medium_patient[!, :calculated_swcs_small];
        label="Small Inserts",
        color=:red,
    )
    lines!([-1000, 1000], [-1000, 1000]; label="Unity")
    lines!(collect(1:1000), fitted_line_reprod_s_mp; linestyle=:dashdot, label="Fitted Line")
	create_textbox(f[2, 2], coefficient_reprod_s_mp, r_squared_reprod_s_mp, rms_values_reprod_s_mp)
	

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

    # save(joinpath(dirname(pwd()),"figures", FIGURE_PATH, "reproducibility.png"), f)
    return f
end

# ╔═╡ c4c1495f-3366-4f39-be4c-667d5895bf17
with_theme(reprod_mp, medphys_theme)

# ╔═╡ 2afe38f3-f10a-47e3-a3e7-32b4523e17c5
md"""
#### Large Patient
"""

# ╔═╡ 83d13107-ad71-453e-b065-f009d7bfdc55
md"""
##### Integrated
"""

# ╔═╡ 4b4a25b0-2973-49cb-9be4-265e00b9417f
array_i_lp = hcat(
	df_i_large_patient[!, :calculated_mass_large], 
	df_i_large_patient[!, :calculated_mass_medium], 
	df_i_large_patient[!, :calculated_mass_small]
);

# ╔═╡ de7d4719-5337-43c1-a9d7-3d32f2560388
array_i_r_lp = hcat(
    df_i_large_patient_r[!, :calculated_mass_large],
    df_i_large_patient_r[!, :calculated_mass_medium],
    df_i_large_patient_r[!, :calculated_mass_small],
);

# ╔═╡ 227f49dc-e9fc-4560-b406-ebc959d9328f
df_i_large_lp, df_i_r_large_lp, df_i_medium_lp, df_i_r_medium_lp, df_i_small_lp, df_i_r_small_lp = remove_false_negatives(df_i_large_patient, array_i_lp, df_i_large_patient_r, array_i_r_lp);

# ╔═╡ 007e7280-0cf6-45c4-b6d6-240a901eac86
begin
	r_squared_reprod_i_lp, rms_values_reprod_i_lp, fitted_line_reprod_i_lp, coefficient_reprod_i_lp = prepare_linear_regression(df_i_r_large_lp, df_i_r_medium_lp, df_i_r_small_lp, df_i_large_lp, df_i_medium_lp, df_i_small_lp)
	
	r_squared_reprod_i_lp = round.(r_squared_reprod_i_lp; digits=3)
	rms_values_reprod_i_lp = round.(rms_values_reprod_i_lp; digits=3)
	coefficient_reprod_i_lp = round.(coefficient_reprod_i_lp; digits=3)
end

# ╔═╡ 3963efc0-8c06-4672-9c6e-3ae256e10d4d
md"""
##### Volume Fraction
"""

# ╔═╡ a832370e-af82-444a-9428-7114649ba203
array_vf_lp = hcat(
	df_vf_large_patient[!, :calculated_mass_large], 
	df_vf_large_patient[!, :calculated_mass_medium], 
	df_vf_large_patient[!, :calculated_mass_small]
);

# ╔═╡ 2969e17f-d3f5-4438-ac4a-05370df59a44
array_vf_r_lp = hcat(
    df_vf_large_patient_r[!, :calculated_mass_large],
    df_vf_large_patient_r[!, :calculated_mass_medium],
    df_vf_large_patient_r[!, :calculated_mass_small],
);

# ╔═╡ c4216708-b4b7-4625-831c-d7b9cce57ef9
df_vf_large_lp, df_vf_r_large_lp, df_vf_medium_lp, df_vf_r_medium_lp, df_vf_small_lp, df_vf_r_small_lp = remove_false_negatives(df_vf_large_patient, array_vf_lp, df_vf_large_patient_r, array_vf_r_lp);

# ╔═╡ 32ae56f6-2e86-4c6b-9343-8535a975da2b
begin
	r_squared_reprod_vf_lp, rms_values_reprod_vf_lp, fitted_line_reprod_vf_lp, coefficient_reprod_vf_lp = prepare_linear_regression(df_vf_r_large_lp, df_vf_r_medium_lp, df_vf_r_small_lp, df_vf_large_lp, df_vf_medium_lp, df_vf_small_lp)
	
	r_squared_reprod_vf_lp = round.(r_squared_reprod_vf_lp; digits=3)
	rms_values_reprod_vf_lp = round.(rms_values_reprod_vf_lp; digits=3)
	coefficient_reprod_vf_lp = round.(coefficient_reprod_vf_lp; digits=3)
end

# ╔═╡ ea5a016d-bed9-4f00-9960-5305e548fa28
md"""
##### Agatston
"""

# ╔═╡ 39132f1c-7a92-47b5-9cd7-2ee9b61cdf89
array_a_lp = hcat(
	df_a_large_patient[!, :calculated_mass_large], 
	df_a_large_patient[!, :calculated_mass_medium], 
	df_a_large_patient[!, :calculated_mass_small]
);

# ╔═╡ 13a516b5-6cb6-4384-80a6-78ffe7078cfc
array_a_r_lp = hcat(
    df_a_large_patient_r[!, :calculated_mass_large],
    df_a_large_patient_r[!, :calculated_mass_medium],
    df_a_large_patient_r[!, :calculated_mass_small],
);

# ╔═╡ 76194267-a83f-4aac-90a2-5e724ca1df60
df_a_large_lp, df_a_r_large_lp, df_a_medium_lp, df_a_r_medium_lp, df_a_small_lp, df_a_r_small_lp = remove_false_negatives(df_a_large_patient, array_a_lp, df_a_large_patient_r, array_a_r_lp);

# ╔═╡ 68a7b68b-0d92-4c7c-ad12-342d1c54efc8
begin
	r_squared_reprod_a_lp, rms_values_reprod_a_lp, fitted_line_reprod_a_lp, coefficient_reprod_a_lp = prepare_linear_regression(df_a_r_large_lp, df_a_r_medium_lp, df_a_r_small_lp, df_a_large_lp, df_a_medium_lp, df_a_small_lp)
	
	r_squared_reprod_a_lp = round.(r_squared_reprod_a_lp; digits=3)
	rms_values_reprod_a_lp = round.(rms_values_reprod_a_lp; digits=3)
	coefficient_reprod_a_lp = round.(coefficient_reprod_a_lp; digits=3)
end

# ╔═╡ 2dc008cf-5570-4a7e-a046-56ce1880b56e
md"""
##### Spatially Weighted
"""

# ╔═╡ 350a77af-9f11-4aea-aae5-31e30e1ccdc2
array_s_lp = hcat(
	df_s_large_patient[!, :calculated_swcs_large], 
	df_s_large_patient[!, :calculated_swcs_medium], 
	df_s_large_patient[!, :calculated_swcs_small]
);

# ╔═╡ 004bedd1-8237-4271-ab70-cd412d444027
array_s_r_lp = hcat(
    df_s_large_patient_r[!, :calculated_swcs_large],
    df_s_large_patient_r[!, :calculated_swcs_medium],
    df_s_large_patient_r[!, :calculated_swcs_small],
);

# ╔═╡ fb79376e-26a9-49e2-a6ba-301596d7e609
df_s_large_lp, df_s_r_large_lp, df_s_medium_lp, df_s_r_medium_lp, df_s_small_lp, df_s_r_small_lp = remove_false_negatives(df_s_large_patient, array_s_lp, df_s_large_patient_r, array_s_r_lp);

# ╔═╡ 221689f5-ac9a-4a69-a366-615aec7513b2
begin
	r_squared_reprod_s_lp, rms_values_reprod_s_lp, fitted_line_reprod_s_lp, coefficient_reprod_s_lp = prepare_linear_regression(df_s_r_large_lp, df_s_r_medium_lp, df_s_r_small_lp, df_s_large_lp, df_s_medium_lp, df_s_small_lp)
	
	r_squared_reprod_s_lp = round.(r_squared_reprod_s_lp; digits=3)
	rms_values_reprod_s_lp = round.(rms_values_reprod_s_lp; digits=3)
	coefficient_reprod_s_lp = round.(coefficient_reprod_s_lp; digits=3)
end

# ╔═╡ b65ed831-7a51-4035-b28a-24ae9ec64a01
function reprod_lp()
    f = Figure()

    ##-- A --##
    ax = Axis(
		f[1, 1],
	    xticks = [0, 50, 100, 150, 200],
	    yticks = [0, 50, 100, 150, 200],
	    xlabel = "Mass 1 (mg)",
	    ylabel = "Mass 2 (mg)",
	    title = "Integrated (Large Patient)",
		)
    scatter!(
		df_i_large_patient_r[!, :calculated_mass_large], 
		df_i_large_patient[!, :calculated_mass_large]
	)
    scatter!(
		df_i_large_patient_r[!, :calculated_mass_medium], 
		df_i_large_patient[!, :calculated_mass_medium]
	)
    scatter!(
        df_i_large_patient_r[!, :calculated_mass_small],
        df_i_large_patient[!, :calculated_mass_small];
        color=:red,
    )
    lines!([-1000, 1000], [-1000, 1000]; label="Unity")
    lines!(collect(1:1000), fitted_line_reprod_i_lp; linestyle=:dashdot)
	create_textbox(f[1, 1], coefficient_reprod_i_lp, r_squared_reprod_i_lp, rms_values_reprod_i_lp)

    xlims!(ax; low=0, high=200)
    ylims!(ax; low=0, high=200)

	##-- B --##
    ax = Axis(
		f[2, 1],
		xticks = [0, 50, 100, 150, 200],
    	yticks = [0, 50, 100, 150, 200],
    	xlabel = "Mass 1 (mg)",
    	ylabel = "Mass 2 (mg)",
    	title = "Volume Fraction (Large Patient)",
	)
    scatter!(
        df_vf_large_patient_r[!, :calculated_mass_large],
        df_vf_large_patient[!, :calculated_mass_large];
        label="Large Inserts",
    )
    scatter!(
        df_vf_large_patient_r[!, :calculated_mass_medium],
        df_vf_large_patient[!, :calculated_mass_medium];
        label="Medium Inserts",
    )
    scatter!(
        df_vf_large_patient_r[!, :calculated_mass_small],
        df_vf_large_patient[!, :calculated_mass_small];
        label="Small Inserts",
        color=:red,
    )
    lines!([-1000, 1000], [-1000, 1000]; label="Unity")
    lines!(collect(1:1000), fitted_line_reprod_vf_lp; linestyle=:dashdot)
	create_textbox(f[2, 1], coefficient_reprod_vf_lp, r_squared_reprod_vf_lp, rms_values_reprod_vf_lp)

    xlims!(ax; low=0, high=200)
    ylims!(ax; low=0, high=200)


    ##-- C --##
    ax = Axis(
		f[1, 2],
		xticks = [0, 50, 100, 150, 200],
    	yticks = [0, 50, 100, 150, 200],
    	xlabel = "Mass 1 (mg)",
    	ylabel = "Mass 2 (mg)",
    	title = "Agatston (Large Patient)",
	)
    scatter!(
        df_a_large_patient_r[!, :calculated_mass_large],
        df_a_large_patient[!, :calculated_mass_large];
        label="Large Inserts",
    )
    scatter!(
        df_a_large_patient_r[!, :calculated_mass_medium],
        df_a_large_patient[!, :calculated_mass_medium];
        label="Medium Inserts",
    )
    scatter!(
        df_a_large_patient_r[!, :calculated_mass_small],
        df_a_large_patient[!, :calculated_mass_small];
        label="Small Inserts",
        color=:red,
    )
    lines!([-1000, 1000], [-1000, 1000]; label="Unity")
    lines!(collect(1:1000), fitted_line_reprod_a_lp; linestyle=:dashdot)
	create_textbox(f[1, 2], coefficient_reprod_a_lp, r_squared_reprod_a_lp, rms_values_reprod_a_lp)
	
    xlims!(ax; low=0, high=200)
    ylims!(ax; low=0, high=200)


    # ##-- D --##
    ax = Axis(
		f[2, 2],
		xticks = [0, 125, 250, 375, 500],
    	yticks = [0, 125, 250, 375, 500],
    	xlabel = "SWCS 1",
    	ylabel = "SWCS 2",
    	title = "Spatially Weighted (Large Patient)",
	)
    scatter!(
        df_s_large_patient_r[!, :calculated_swcs_large],
        df_s_large_patient[!, :calculated_swcs_large];
        label="Large Inserts",
    )
    scatter!(
        df_s_large_patient_r[!, :calculated_swcs_medium],
        df_s_large_patient[!, :calculated_swcs_medium];
        label="Medium Inserts",
    )
    scatter!(
        df_s_large_patient_r[!, :calculated_swcs_small],
        df_s_large_patient[!, :calculated_swcs_small];
        label="Small Inserts",
        color=:red,
    )
    lines!([-1000, 1000], [-1000, 1000]; label="Unity")
    lines!(collect(1:1000), fitted_line_reprod_s_lp; linestyle=:dashdot, label="Fitted Line")
	create_textbox(f[2, 2], coefficient_reprod_s_lp, r_squared_reprod_s_lp, rms_values_reprod_s_lp)
	

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

    # save(joinpath(dirname(pwd()),"figures", FIGURE_PATH, "reproducibility.png"), f)
    return f
end

# ╔═╡ 7684a47b-a26d-4ab1-bc9d-54b4513c258a
with_theme(reprod_lp, medphys_theme)

# ╔═╡ c50117bb-bbb7-4347-b4bb-d446677f4509
md"""
## Sensitivity and Specificity
"""

# ╔═╡ 75e0f9d7-ef66-42bd-bf06-cebecb09a526
std_level = 1.5

# ╔═╡ 9ed516b4-232b-45c5-af90-2964ed2e6a36
md"""
#### False Negative
"""

# ╔═╡ 10725f84-8ce4-454a-8c6e-41cdaa25dd1a
md"""
##### SWCS
"""

# ╔═╡ 02dd53b4-3f0f-437f-885b-1c7f4c78d8e5
begin
	false_negative_s = []
	for i in 1:3:nrow(df_s)-2
		mean_s, std_s = mean(df_s[i:i+2, :swcs_bkg]), std(df_s[i:i+2, :swcs_bkg])*std_level
		array_s = hcat(df_s[i:i+2, :calculated_swcs_large], df_s[i:i+2, :calculated_swcs_medium], df_s[i:i+2, :calculated_swcs_small]);
		neg = length(findall(x -> x <= (mean_s + std_s), array_s))
		push!(false_negative_s, neg)
	end
end

# ╔═╡ 608e3046-5df2-4491-992c-46233e93f4c8
total_zero_s = sum(false_negative_s)

# ╔═╡ 65dc9755-c46f-49ea-8ecf-d68bcb42a628
md"""
##### Agatston
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

# ╔═╡ 0e5b9d7d-aa2f-44fa-9ca9-1d8b7f373f7b
r_squared_reprod_i, r_squared_reprod_vf, r_squared_reprod_a, r_squared_reprod_s

# ╔═╡ dfcd85cd-586c-49aa-96e2-2e265fb5f33d
coefficient_reprod_i, coefficient_reprod_vf, coefficient_reprod_a, coefficient_reprod_s

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
with_theme(reprod, medphys_theme)

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
##### Integrated
"""

# ╔═╡ 958cc5b3-1434-4eed-b7e7-7b3e36fe848c
begin
	false_negative_i = []
	for i in 1:3:nrow(df_i)-2
		mean_i, std_i = mean(df_i[i:i+2, :mass_bkg]), std(df_i[i:i+2, :mass_bkg])*std_level 
		array_i = hcat(df_i[i:i+2, :calculated_mass_large], df_i[i:i+2, :calculated_mass_medium], df_i[i:i+2, :calculated_mass_small]);
		neg = length(findall(x -> x <= mean_i + std_i, array_i))
		push!(false_negative_i, neg)
	end
end

# ╔═╡ f029a8c7-7cc1-49c1-8337-c76ac0d585cd
total_zero_i = sum(false_negative_i)

# ╔═╡ 73852846-3fd1-4246-8bd3-afabd7b8e480
md"""
##### Volume Fraction
"""

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

# ╔═╡ 4eef56ad-0ea4-422c-837b-3aa9846a2458


# ╔═╡ b708d2c9-f868-4814-b1b7-6908317aa5ba
total_zero_vf = sum(false_negative_vf)

# ╔═╡ 5e32f34c-2f32-457e-a9cf-59184170b6ed
total_zero_i, total_zero_vf, total_zero_s, num_zero_a, total_cac

# ╔═╡ 1d04d29f-401d-45c4-99ff-3bbdb0fb4f6b
md"""
#### False Positive
"""

# ╔═╡ 559f77d3-69c0-4155-8cc0-1f669f929bc3
md"""
##### SWCS
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

# ╔═╡ 802dc604-cd58-40db-ad6c-132ad38d7c3c
df_s

# ╔═╡ 045cb9a5-12ac-4016-a39c-10b2b15e90f9
total_zero_s_pos = sum(false_positive_s)

# ╔═╡ 1f0df2eb-a434-41fa-b690-4fcdae434a30
md"""
##### Agatston
"""

# ╔═╡ 4a0c4183-3ccc-468a-aa04-75719c21750d
array_a_pos = df_a[!, :mass_bkg];

# ╔═╡ 9d8911ab-05e3-4a07-8f2e-d9e1540447a5
total_cac_pos = length(array_a_pos)

# ╔═╡ f542ed82-d2f3-4033-8083-b2e7a3a54650
total_zero_a_pos = length(findall(x -> x > 0, array_a_pos))

# ╔═╡ 029cffc3-e313-46e6-9313-b50900b4dd1a
md"""
##### Integrated
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
##### Volume Fraction
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

# ╔═╡ ac33fc36-b7c2-47f1-b513-8a02f5fdd0ad
md"""
### Density
"""

# ╔═╡ 8aae5583-bf94-42f2-a7b7-b6b9a3a3c3eb
md"""
#### Normal Density
"""

# ╔═╡ a1f3606e-69c9-45c6-9602-4367ca63ae48
md"""
##### False Negative
"""

# ╔═╡ 32fb4eca-3ef6-4254-8fcf-79150c8ae39d
md"""
###### SWCS
"""

# ╔═╡ a894fc59-25b4-4536-bac6-eb2352c5baaf
begin
	false_negative_s_nd = []
	for i in 1:3:nrow(df_s_normal)-2
		mean_s, std_s = mean(df_s_normal[i:i+2, :swcs_bkg]), std(df_s_normal[i:i+2, :swcs_bkg])*std_level
		array_s = hcat(df_s_normal[i:i+2, :calculated_swcs_large], df_s_normal[i:i+2, :calculated_swcs_medium], df_s_normal[i:i+2, :calculated_swcs_small]);
		neg = length(findall(x -> x <= (mean_s + std_s), array_s))
		push!(false_negative_s_nd, neg)
	end
end

# ╔═╡ aeb44f0c-c77d-4f5e-8ad1-dbc3b7a0d7af
total_zero_s_nd = sum(false_negative_s_nd)

# ╔═╡ e0937eba-d445-4b73-a266-663e56c33246
md"""
###### Agatston
"""

# ╔═╡ e166f587-2a51-4135-86ce-648a5f92c0f1
num_zero_a_nd = length(findall(x -> x <= 0, array_a_nd))

# ╔═╡ 39f85f44-558d-4f55-9ed0-1db50e309847
total_cac_nd = length(array_a_nd)

# ╔═╡ 47d4c62d-6d01-4999-91ca-b0b9679a9bb5
md"""
###### Integrated
"""

# ╔═╡ 07a7060b-5490-4e1a-9739-eb9bf790a999
begin
	false_negative_i_nd = []
	for i in 1:3:nrow(df_i_normal)-2
		mean_i, std_i = mean(df_i_normal[i:i+2, :mass_bkg]), std(df_i_normal[i:i+2, :mass_bkg])*std_level 
		array_i = hcat(df_i_normal[i:i+2, :calculated_mass_large], df_i_normal[i:i+2, :calculated_mass_medium], df_i_normal[i:i+2, :calculated_mass_small]);
		neg = length(findall(x -> x <= mean_i + std_i, array_i))
		push!(false_negative_i_nd, neg)
	end
end

# ╔═╡ 7003f82d-d6a7-41b4-98e2-0cbefc4f1e24
total_zero_i_nd = sum(false_negative_i_nd)

# ╔═╡ 859ba958-96d5-4e26-9c98-963b1a3bcddf
md"""
###### Volume Fraction
"""

# ╔═╡ b780e9ea-4b9f-4ad9-b6dc-6c9b5070863d
begin
	false_negative_vf_nd = []
	for i in 1:3:nrow(df_vf_normal)-2
		mean_vf, std_vf = mean(df_vf_normal[i:i+2, :mass_bkg]), std(df_vf_normal[i:i+2, :mass_bkg])*std_level 
		array_vf = hcat(df_vf_normal[i:i+2, :calculated_mass_large], df_vf_normal[i:i+2, :calculated_mass_medium], df_vf_normal[i:i+2, :calculated_mass_small]);
		neg = length(findall(x -> x <= mean_vf + std_vf, array_vf))
		push!(false_negative_vf_nd, neg)
	end
end

# ╔═╡ 7734498b-ea86-47f7-8f17-c842a58e9c93
total_zero_vf_nd = sum(false_negative_vf_nd)

# ╔═╡ fa11f18f-3e05-435c-ab18-fa3faba49871
md"""
##### False Positive
"""

# ╔═╡ 2960b18c-f0df-40b6-8284-8d843ccd3783
md"""
###### SWCS
"""

# ╔═╡ 62d4be9f-046d-4b81-946a-8ddf359dc179
begin
	false_positive_s_nd = []
	for i in 1:3:nrow(df_s_normal)-2
		mean_s, std_s = mean(df_s_normal[i:i+2, :swcs_bkg]), std(df_s_normal[i:i+2, :swcs_bkg])*std_level
		array_s_pos = df_s_normal[i:i+2, :swcs_bkg]
		pos = length(findall(x -> x > (mean_s + std_s), array_s_pos))
		push!(false_positive_s_nd, pos)
	end
end

# ╔═╡ 295c146f-ffa0-4d62-9352-0b21c49584a1
total_zero_s_pos_nd = sum(false_positive_s_nd)

# ╔═╡ 7efe886e-76ca-4269-b01f-c87fd678032a
md"""
###### Agatston
"""

# ╔═╡ fbcdaf78-f366-428d-ad63-e0b17cfc1fa4
array_a_pos_nd = df_a_normal[!, :mass_bkg];

# ╔═╡ ad1303e7-a9eb-4655-a176-4a54b77db1b7
total_cac_pos_nd = length(array_a_pos_nd)

# ╔═╡ ad25959e-0cda-4292-8052-070ba4d73b59
total_zero_a_pos_nd = length(findall(x -> x > 0, array_a_pos_nd))

# ╔═╡ a94b4849-c2c9-4a22-b238-6212c9867f9f
md"""
###### Integrated
"""

# ╔═╡ c389a74e-4ab1-48f4-ac37-fcc0b794affb
begin
	false_positive_i_nd = []
	for i in 1:3:nrow(df_i_normal)-2
		mean_i, std_i = mean(df_i_normal[i:i+2, :mass_bkg]), std(df_i_normal[i:i+2, :mass_bkg])*std_level
		array_i_pos = df_i_normal[i:i+2, :mass_bkg]
		pos = length(findall(x -> x > (mean_i + std_level), array_i_pos))
		push!(false_positive_i_nd, pos)
	end
end

# ╔═╡ 36b574e8-6eba-4038-b889-d472c4fdc824
total_zero_i_pos_nd = sum(false_positive_i_nd)

# ╔═╡ a60d6f3e-7c31-47f9-901f-4b110d572b3b
md"""
###### Volume Fraction
"""

# ╔═╡ 0f60944f-e615-4117-bde6-9efaa9900806
begin
	false_positive_vf_nd = []
	for i in 1:3:nrow(df_vf_normal)-2
		mean_vf, std_vf = mean(df_vf_normal[i:i+2, :mass_bkg]), std(df_vf_normal[i:i+2, :mass_bkg])*std_level
		array_vf_pos = df_vf_normal[i:i+2, :mass_bkg]
		pos = length(findall(x -> x > (mean_vf + std_vf), array_vf_pos))
		push!(false_positive_vf_nd, pos)
	end
end

# ╔═╡ f1434cec-5655-46fe-af66-8eea8f639d94
total_zero_vf_pos_nd = sum(false_positive_vf_nd)

# ╔═╡ e4c021dc-2fee-401f-8bb8-4f72f0186231
function sensitivity_specificity_nd()
    f = Figure(resolution=(800, 800))
    colors = Makie.wong_colors()

    ##-- A --##
    ax = Axis(
		f[1, 1]; 
		xticks = (1:4, ["Integrated", "Volume Fraction", "Spatially Weighted", "Agatston"]),
		title = "False-Negative (Normal Density)",
		ylabel = "False-Negative (%)",
		yticks = [0, 25, 50, 75, 100]
	)

    table = [1, 2, 3, 4]
	h1 = (total_zero_i_nd / total_cac_nd) * 100
	h2 = (total_zero_vf_nd / total_cac_nd) * 100
    h3 = (total_zero_s_nd / total_cac_nd) * 100
    h4 = (num_zero_a_nd / total_cac_nd) * 100
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
		title = "False-Positive (Normal Density)",
		ylabel = "False-Positive (%)",
		yticks = [0, 25, 50, 75, 100]
	)

    table = [1, 2, 3, 4]
	h1 = (total_zero_i_pos_nd / total_cac_pos_nd) * 100
	h2 = (total_zero_vf_pos_nd / total_cac_pos_nd) * 100
    h3 = (total_zero_s_pos_nd / total_cac_pos_nd) * 100
    h4 = (total_zero_a_pos_nd / total_cac_pos_nd) * 100
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

    # save(joinpath(dirname(pwd()),"figures", FIGURE_PATH, "sensitivity_specificity.png"), f)


    return f
end

# ╔═╡ ce3956ea-9947-4bc9-91df-39079eec40ee
with_theme(sensitivity_specificity_nd, medphys_theme)

# ╔═╡ 0eff2135-b37f-413a-add6-80a07197d346
md"""
#### Low Density
"""

# ╔═╡ 8d88faa9-df59-4a4c-ad03-29612c2255a9
md"""
##### False Negative
"""

# ╔═╡ 2a2de10a-4792-45bd-ab03-9679a46ff6cd
md"""
###### SWCS
"""

# ╔═╡ fb2828fa-0c4f-409e-b739-cb1f2c135f8f
begin
	false_negative_s_ld = []
	for i in 1:3:nrow(df_s_low)-2
		mean_s, std_s = mean(df_s_low[i:i+2, :swcs_bkg]), std(df_s_low[i:i+2, :swcs_bkg])*std_level
		array_s = hcat(df_s_low[i:i+2, :calculated_swcs_large], df_s_low[i:i+2, :calculated_swcs_medium], df_s_low[i:i+2, :calculated_swcs_small]);
		neg = length(findall(x -> x <= (mean_s + std_s), array_s))
		push!(false_negative_s_ld, neg)
	end
end

# ╔═╡ d5ec6fbc-1577-4950-87e6-51a7e1b7da75
total_zero_s_ld = sum(false_negative_s_ld)

# ╔═╡ ed9e7420-0fb8-4882-8cd8-f38b1b254684
md"""
###### Agatston
"""

# ╔═╡ 0d1900e6-3681-4f20-afb6-63ccd92bd388
num_zero_a_ld = length(findall(x -> x <= 0, array_a_ld))

# ╔═╡ f266aca7-533e-48ac-a104-8f2940b97eb3
total_cac_ld = length(array_a_ld)

# ╔═╡ 9c2a5fda-d887-4101-9e37-1ef6ed71c4aa
md"""
###### Integrated
"""

# ╔═╡ 8f11feed-fa28-405d-a8c2-105e15b065d3
begin
	false_negative_i_ld = []
	for i in 1:3:nrow(df_i_low)-2
		mean_i, std_i = mean(df_i_low[i:i+2, :mass_bkg]), std(df_i_low[i:i+2, :mass_bkg])*std_level 
		array_i = hcat(df_i_low[i:i+2, :calculated_mass_large], df_i_low[i:i+2, :calculated_mass_medium], df_i_low[i:i+2, :calculated_mass_small]);
		neg = length(findall(x -> x <= mean_i + std_i, array_i))
		push!(false_negative_i_ld, neg)
	end
end

# ╔═╡ a77cc0a3-6eec-4bb3-aa25-299d8051dead
total_zero_i_ld = sum(false_negative_i_ld)

# ╔═╡ b8c525da-2805-4118-b116-10aeaad4cba5
md"""
###### Volume Fraction
"""

# ╔═╡ 310cb551-ca38-49d7-b570-38eb5221f2bb
begin
	false_negative_vf_ld = []
	for i in 1:3:nrow(df_vf_low)-2
		mean_vf, std_vf = mean(df_vf_low[i:i+2, :mass_bkg]), std(df_vf_low[i:i+2, :mass_bkg])*std_level 
		array_vf = hcat(df_vf_low[i:i+2, :calculated_mass_large], df_vf_low[i:i+2, :calculated_mass_medium], df_vf_low[i:i+2, :calculated_mass_small]);
		neg = length(findall(x -> x <= mean_vf + std_vf, array_vf))
		push!(false_negative_vf_ld, neg)
	end
end

# ╔═╡ 521a8a99-d5b1-4d38-8365-1cfbf64a1bb0
total_zero_vf_ld = sum(false_negative_vf_ld)

# ╔═╡ 1d65d6b1-16f4-4845-beb5-49d307be8a8b
md"""
##### False Positive
"""

# ╔═╡ 348bb217-e3b1-43ce-86f2-76629d1ca00e
md"""
###### SWCS
"""

# ╔═╡ 2a2492cd-8533-4149-bbe3-d1bdb70f22e2
begin
	false_positive_s_ld = []
	for i in 1:3:nrow(df_s_low)-2
		mean_s, std_s = mean(df_s_low[i:i+2, :swcs_bkg]), std(df_s_low[i:i+2, :swcs_bkg])*std_level
		array_s_pos = df_s_low[i:i+2, :swcs_bkg]
		pos = length(findall(x -> x > (mean_s + std_s), array_s_pos))
		push!(false_positive_s_ld, pos)
	end
end

# ╔═╡ 75a7f90b-1ccd-4770-a90b-ac7a5953443c
total_zero_s_pos_ld = sum(false_positive_s_ld)

# ╔═╡ 9eec9ed3-46a7-4c7e-9f91-9ea0b01afb9f
md"""
###### Agatston
"""

# ╔═╡ 08a34dcc-b4ca-434e-b9df-38666dfcbdc0
array_a_pos_ld = df_a_low[!, :mass_bkg];

# ╔═╡ 4fbfcb02-d325-4d75-86b5-bd944352bca6
total_cac_pos_ld = length(array_a_pos_ld)

# ╔═╡ 1b5074bf-b036-42f1-9853-42b1e53865c7
total_zero_a_pos_ld = length(findall(x -> x > 0, array_a_pos_ld))

# ╔═╡ 7f7ff260-ae32-4d78-9e85-c268ced39cba
md"""
###### Integrated
"""

# ╔═╡ bc9eb230-7ff3-4597-8b70-1a6e58918afb
begin
	false_positive_i_ld = []
	for i in 1:3:nrow(df_i_low)-2
		mean_i, std_i = mean(df_i_low[i:i+2, :mass_bkg]), std(df_i_low[i:i+2, :mass_bkg])*std_level
		array_i_pos = df_i_low[i:i+2, :mass_bkg]
		pos = length(findall(x -> x > (mean_i + std_level), array_i_pos))
		push!(false_positive_i_ld, pos)
	end
end

# ╔═╡ ac924022-6f03-4566-9c0a-23c3ae278938
total_zero_i_pos_ld = sum(false_positive_i_ld)

# ╔═╡ d6392619-9f11-4824-8168-c74fd08eece5
md"""
###### Volume Fraction
"""

# ╔═╡ 6950ed5a-fb5d-4b9e-a17e-314760933ed0
begin
	false_positive_vf_ld = []
	for i in 1:3:nrow(df_vf_low)-2
		mean_vf, std_vf = mean(df_vf_low[i:i+2, :mass_bkg]), std(df_vf_low[i:i+2, :mass_bkg])*std_level
		array_vf_pos = df_vf_low[i:i+2, :mass_bkg]
		pos = length(findall(x -> x > (mean_vf + std_vf), array_vf_pos))
		push!(false_positive_vf_ld, pos)
	end
end

# ╔═╡ 632e3310-7c7a-4d91-a2cd-d00fc2633492
total_zero_vf_pos_ld = sum(false_positive_vf_ld)

# ╔═╡ ab20b172-35c7-49d6-8b57-ae8c2a310a26
function sensitivity_specificity_ld()
    f = Figure(resolution=(800, 800))
    colors = Makie.wong_colors()

    ##-- A --##
    ax = Axis(
		f[1, 1]; 
		xticks = (1:4, ["Integrated", "Volume Fraction", "Spatially Weighted", "Agatston"]),
		title = "False-Negative (Low Density)",
		ylabel = "False-Negative (%)",
		yticks = [0, 25, 50, 75, 100]
	)

    table = [1, 2, 3, 4]
	h1 = (total_zero_i_ld / total_cac_ld) * 100
	h2 = (total_zero_vf_ld / total_cac_ld) * 100
    h3 = (total_zero_s_ld / total_cac_ld) * 100
    h4 = (num_zero_a_ld / total_cac_ld) * 100
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
		title = "False-Positive (Low Density)",
		ylabel = "False-Positive (%)",
		yticks = [0, 25, 50, 75, 100]
	)

    table = [1, 2, 3, 4]
	h1 = (total_zero_i_pos_ld / total_cac_pos_ld) * 100
	h2 = (total_zero_vf_pos_ld / total_cac_pos_ld) * 100
    h3 = (total_zero_s_pos_ld / total_cac_pos_ld) * 100
    h4 = (total_zero_a_pos_ld / total_cac_pos_ld) * 100
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

    save(joinpath(dirname(pwd()),"figures", FIGURE_PATH, "sensitivity_specificity_ld.png"), f)


    return f
end

# ╔═╡ 6e0c8aa0-89ac-4e5b-b5e1-f83ad87c246a
with_theme(sensitivity_specificity_ld, medphys_theme)

# ╔═╡ 91f9ec30-e6e8-4bc7-83c3-b49c98b0b32d
md"""
### Patient Size
"""

# ╔═╡ 15dea22a-662b-43c4-ac6b-33ea06a5f2c4
md"""
#### Small Patient
"""

# ╔═╡ a6b18ed3-080c-4691-bd53-83b878cbec74
md"""
##### False Negative
"""

# ╔═╡ 0270d409-cb8f-45b2-b30c-cdd398e7fe5f
md"""
###### SWCS
"""

# ╔═╡ 649ca894-ab4b-4a60-96f7-c91be28a092b
begin
	false_negative_s_sp = []
	for i in 1:3:nrow(df_s_small_patient)-2
		mean_s, std_s = mean(df_s_small_patient[i:i+2, :swcs_bkg]), std(df_s_small_patient[i:i+2, :swcs_bkg])*std_level
		array_s = hcat(df_s_small_patient[i:i+2, :calculated_swcs_large], df_s_small_patient[i:i+2, :calculated_swcs_medium], df_s_small_patient[i:i+2, :calculated_swcs_small]);
		neg = length(findall(x -> x <= (mean_s + std_s), array_s))
		push!(false_negative_s_sp, neg)
	end
end

# ╔═╡ d86480a0-041e-402b-a56d-f10a9903529f
total_zero_s_sp = sum(false_negative_s_sp)

# ╔═╡ 85dc1af9-033d-4371-b0e9-2888d1478bb8
md"""
###### Agatston
"""

# ╔═╡ 722824a9-e6b3-4a63-8131-f2043f7269f1
num_zero_a_sp = length(findall(x -> x <= 0, array_a_sp))

# ╔═╡ 04e8dfd9-2476-4a54-ad4c-cc381fad3087
total_cac_sp = length(array_a_sp)

# ╔═╡ e2a1b07c-1b14-464c-8589-2b3fd8d55d00
md"""
###### Integrated
"""

# ╔═╡ e35098ae-4147-41a4-b4ff-dfa6dcc7ff2f
begin
	false_negative_i_sp = []
	for i in 1:3:nrow(df_i_small_patient)-2
		mean_i, std_i = mean(df_i_small_patient[i:i+2, :mass_bkg]), std(df_i_small_patient[i:i+2, :mass_bkg])*std_level 
		array_i = hcat(df_i_small_patient[i:i+2, :calculated_mass_large], df_i_small_patient[i:i+2, :calculated_mass_medium], df_i_small_patient[i:i+2, :calculated_mass_small]);
		neg = length(findall(x -> x <= mean_i + std_i, array_i))
		push!(false_negative_i_sp, neg)
	end
end

# ╔═╡ 4bee3c44-f9c8-4428-b6d2-471ec732f3d1
total_zero_i_sp = sum(false_negative_i_sp)

# ╔═╡ 2c888b16-5385-4ea1-a5b1-e8a5a08c89c9
md"""
###### Volume Fraction
"""

# ╔═╡ 1ad040bb-e0a4-421b-b80b-90397db8121b
begin
	false_negative_vf_sp = []
	for i in 1:3:nrow(df_vf_small_patient)-2
		mean_vf, std_vf = mean(df_vf_small_patient[i:i+2, :mass_bkg]), std(df_vf_small_patient[i:i+2, :mass_bkg])*std_level 
		array_vf = hcat(df_vf_small_patient[i:i+2, :calculated_mass_large], df_vf_small_patient[i:i+2, :calculated_mass_medium], df_vf_small_patient[i:i+2, :calculated_mass_small]);
		neg = length(findall(x -> x <= mean_vf + std_vf, array_vf))
		push!(false_negative_vf_sp, neg)
	end
end

# ╔═╡ 0fb3dedd-782b-4f2f-add6-041c72b51d12
total_zero_vf_sp = sum(false_negative_vf_sp)

# ╔═╡ cfe68dce-c66c-4c98-9bc9-340859d34e82
md"""
##### False Positive
"""

# ╔═╡ 8e68a1fd-a9d7-4bd9-a863-0543029d9a90
md"""
###### SWCS
"""

# ╔═╡ 6d16a3ff-4de7-4fab-9c37-05497b2c895e
begin
	false_positive_s_sp = []
	for i in 1:3:nrow(df_s_small_patient)-2
		mean_s, std_s = mean(df_s_small_patient[i:i+2, :swcs_bkg]), std(df_s_small_patient[i:i+2, :swcs_bkg])*std_level
		array_s_pos = df_s_small_patient[i:i+2, :swcs_bkg]
		pos = length(findall(x -> x > (mean_s + std_s), array_s_pos))
		push!(false_positive_s_sp, pos)
	end
end

# ╔═╡ 6bee3d94-ce21-4e54-82a1-2fa926096d40
total_zero_s_pos_sp = sum(false_positive_s_sp)

# ╔═╡ 3fae6b4d-4ffe-47b2-92e9-a4e9f509759e
md"""
###### Agatston
"""

# ╔═╡ 8019a1c3-650b-40c2-9bf7-cc8d245f1173
array_a_pos_sp = df_a_small_patient[!, :mass_bkg];

# ╔═╡ 7f9119e0-a71d-49d1-a441-037dbb5be4a2
total_cac_pos_sp = length(array_a_pos_sp)

# ╔═╡ a4f44a8c-70f1-40ff-b87f-9ab1082bf435
total_zero_a_pos_sp = length(findall(x -> x > 0, array_a_pos_sp))

# ╔═╡ 4dfe6af4-21c3-4c13-8c0b-dc2ae437129b
md"""
###### Integrated
"""

# ╔═╡ 22d33469-000e-46d3-8693-f71a4682180f
begin
	false_positive_i_sp = []
	for i in 1:3:nrow(df_i_small_patient)-2
		mean_i, std_i = mean(df_i_small_patient[i:i+2, :mass_bkg]), std(df_i_small_patient[i:i+2, :mass_bkg])*std_level
		array_i_pos = df_i_small_patient[i:i+2, :mass_bkg]
		pos = length(findall(x -> x > (mean_i + std_level), array_i_pos))
		push!(false_positive_i_sp, pos)
	end
end

# ╔═╡ 176af30d-9078-43a6-a5f8-5b93d536349a
total_zero_i_pos_sp = sum(false_positive_i_sp)

# ╔═╡ 48684053-4258-4dd7-94d5-bbc73196ea10
md"""
###### Volume Fraction
"""

# ╔═╡ 5e96e6b5-07cd-4421-bdcc-ad1524a29b62
begin
	false_positive_vf_sp = []
	for i in 1:3:nrow(df_vf_small_patient)-2
		mean_vf, std_vf = mean(df_vf_small_patient[i:i+2, :mass_bkg]), std(df_vf_small_patient[i:i+2, :mass_bkg])*std_level
		array_vf_pos = df_vf_small_patient[i:i+2, :mass_bkg]
		pos = length(findall(x -> x > (mean_vf + std_vf), array_vf_pos))
		push!(false_positive_vf_sp, pos)
	end
end

# ╔═╡ 4db669a2-c259-4dfa-aa19-3862da5380d9
total_zero_vf_pos_sp = sum(false_positive_vf_sp)

# ╔═╡ 32cc0ed0-c025-43db-a107-83af5a8f1571
function sensitivity_specificity_sp()
    f = Figure(resolution=(800, 800))
    colors = Makie.wong_colors()

    ##-- A --##
    ax = Axis(
		f[1, 1]; 
		xticks = (1:4, ["Integrated", "Volume Fraction", "Spatially Weighted", "Agatston"]),
		title = "False-Negative (Small Patient)",
		ylabel = "False-Negative (%)",
		yticks = [0, 25, 50, 75, 100]
	)

    table = [1, 2, 3, 4]
	h1 = (total_zero_i_sp / total_cac_sp) * 100
	h2 = (total_zero_vf_sp / total_cac_sp) * 100
    h3 = (total_zero_s_sp / total_cac_sp) * 100
    h4 = (num_zero_a_sp / total_cac_sp) * 100
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
		title = "False-Positive (Small Patient)",
		ylabel = "False-Positive (%)",
		yticks = [0, 25, 50, 75, 100]
	)

    table = [1, 2, 3, 4]
	h1 = (total_zero_i_pos_sp / total_cac_pos_sp) * 100
	h2 = (total_zero_vf_pos_sp / total_cac_pos_sp) * 100
    h3 = (total_zero_s_pos_sp / total_cac_pos_sp) * 100
    h4 = (total_zero_a_pos_sp / total_cac_pos_sp) * 100
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

    # save(joinpath(dirname(pwd()),"figures", FIGURE_PATH, "sensitivity_specificity.png"), f)


    return f
end

# ╔═╡ 67dcee1b-f488-431e-8c1d-4211bf04696f
with_theme(sensitivity_specificity_sp, medphys_theme)

# ╔═╡ d55054cd-c63b-43c8-b0f9-5f4e60b6415d
md"""
#### Medium Patient
"""

# ╔═╡ 0ebdfc89-960e-45cc-a129-54a003d74d42
md"""
##### False Negative
"""

# ╔═╡ 9a37db8b-7e94-45ce-b459-9c10246270d7
md"""
###### SWCS
"""

# ╔═╡ d57d99ae-4aae-4081-87f5-31ff8d3e85c3
begin
	false_negative_s_mp = []
	for i in 1:3:nrow(df_s_medium_patient)-2
		mean_s, std_s = mean(df_s_medium_patient[i:i+2, :swcs_bkg]), std(df_s_medium_patient[i:i+2, :swcs_bkg])*std_level
		array_s = hcat(df_s_medium_patient[i:i+2, :calculated_swcs_large], df_s_medium_patient[i:i+2, :calculated_swcs_medium], df_s_medium_patient[i:i+2, :calculated_swcs_small]);
		neg = length(findall(x -> x <= (mean_s + std_s), array_s))
		push!(false_negative_s_mp, neg)
	end
end

# ╔═╡ 7c9c9b09-74e4-427e-af9f-f047b41d3ff9
total_zero_s_mp = sum(false_negative_s_mp)

# ╔═╡ 986bfd16-6520-43b9-912a-09f842e14724
md"""
###### Agatston
"""

# ╔═╡ ca0c55ac-6e17-4219-9fc6-0b2288fb7b1c
num_zero_a_mp = length(findall(x -> x <= 0, array_a_mp))

# ╔═╡ 2440eb1b-681c-47e0-b14b-1fdc201d3577
total_cac_mp = length(array_a_mp)

# ╔═╡ f08d6966-7132-4c98-866c-b1ea531acc9e
md"""
###### Integrated
"""

# ╔═╡ b3015cee-3cfa-4968-bd8d-8bcd0a961278
begin
	false_negative_i_mp = []
	for i in 1:3:nrow(df_i_medium_patient)-2
		mean_i, std_i = mean(df_i_medium_patient[i:i+2, :mass_bkg]), std(df_i_medium_patient[i:i+2, :mass_bkg])*std_level 
		array_i = hcat(df_i_medium_patient[i:i+2, :calculated_mass_large], df_i_medium_patient[i:i+2, :calculated_mass_medium], df_i_medium_patient[i:i+2, :calculated_mass_small]);
		neg = length(findall(x -> x <= mean_i + std_i, array_i))
		push!(false_negative_i_mp, neg)
	end
end

# ╔═╡ 34c20a4e-3ec2-45ae-b83a-b3f4f9b0ee46
total_zero_i_mp = sum(false_negative_i_mp)

# ╔═╡ bf25092a-f6ca-4605-994a-a13bbb5b41d3
md"""
###### Volume Fraction
"""

# ╔═╡ dbfc643d-81cd-4c99-8665-79c007115119
begin
	false_negative_vf_mp = []
	for i in 1:3:nrow(df_vf_medium_patient)-2
		mean_vf, std_vf = mean(df_vf_medium_patient[i:i+2, :mass_bkg]), std(df_vf_medium_patient[i:i+2, :mass_bkg])*std_level 
		array_vf = hcat(df_vf_medium_patient[i:i+2, :calculated_mass_large], df_vf_medium_patient[i:i+2, :calculated_mass_medium], df_vf_medium_patient[i:i+2, :calculated_mass_small]);
		neg = length(findall(x -> x <= mean_vf + std_vf, array_vf))
		push!(false_negative_vf_mp, neg)
	end
end

# ╔═╡ 5e94d82c-05ab-4bef-a5bc-d4cb5791a338
total_zero_vf_mp = sum(false_negative_vf_mp)

# ╔═╡ 5f6d3a21-96f1-4618-8de4-3add807572cb
md"""
##### False Positive
"""

# ╔═╡ 9b207e45-1ef3-4742-be7d-dd8f67abea1b
md"""
###### SWCS
"""

# ╔═╡ a07d7c66-c17c-4f17-8cc2-7db04d93a696
begin
	false_positive_s_mp = []
	for i in 1:3:nrow(df_s_medium_patient)-2
		mean_s, std_s = mean(df_s_medium_patient[i:i+2, :swcs_bkg]), std(df_s_medium_patient[i:i+2, :swcs_bkg])*std_level
		array_s_pos = df_s_medium_patient[i:i+2, :swcs_bkg]
		pos = length(findall(x -> x > (mean_s + std_s), array_s_pos))
		push!(false_positive_s_mp, pos)
	end
end

# ╔═╡ 3f0ee19a-55ea-4d75-be27-0a07bb33d909
total_zero_s_pos_mp = sum(false_positive_s_mp)

# ╔═╡ 2392ea54-a176-4d70-a144-06f31c28cda1
md"""
###### Agatston
"""

# ╔═╡ 93f8369f-f26f-4736-b86d-b306e5c53dce
array_a_pos_mp = df_a_medium_patient[!, :mass_bkg];

# ╔═╡ 6b0f51bf-b2be-4d5a-8ca1-733dc6a16628
total_cac_pos_mp = length(array_a_pos_mp)

# ╔═╡ 4b341998-336c-4ec4-a1b4-08e35c155f4c
total_zero_a_pos_mp = length(findall(x -> x > 0, array_a_pos_mp))

# ╔═╡ 8a0d81d1-863a-4899-ad19-fcbd5f87d857
md"""
###### Integrated
"""

# ╔═╡ dc9f382d-f05d-4ba6-8a96-cc963ed4722d
begin
	false_positive_i_mp = []
	for i in 1:3:nrow(df_i_medium_patient)-2
		mean_i, std_i = mean(df_i_medium_patient[i:i+2, :mass_bkg]), std(df_i_medium_patient[i:i+2, :mass_bkg])*std_level
		array_i_pos = df_i_medium_patient[i:i+2, :mass_bkg]
		pos = length(findall(x -> x > (mean_i + std_level), array_i_pos))
		push!(false_positive_i_mp, pos)
	end
end

# ╔═╡ f2cc4db1-65ce-4ebe-a10f-82cd9de5bef7
total_zero_i_pos_mp = sum(false_positive_i_mp)

# ╔═╡ 89e0e212-8ddd-45eb-8054-2045fa3e29f2
md"""
###### Volume Fraction
"""

# ╔═╡ 6cd86ec7-3e1a-4e7a-9f61-4077ca986fd2
begin
	false_positive_vf_mp = []
	for i in 1:3:nrow(df_vf_medium_patient)-2
		mean_vf, std_vf = mean(df_vf_medium_patient[i:i+2, :mass_bkg]), std(df_vf_medium_patient[i:i+2, :mass_bkg])*std_level
		array_vf_pos = df_vf_medium_patient[i:i+2, :mass_bkg]
		pos = length(findall(x -> x > (mean_vf + std_vf), array_vf_pos))
		push!(false_positive_vf_mp, pos)
	end
end

# ╔═╡ 4dceb05b-f463-4789-a79a-3cb9888b6c2e
total_zero_vf_pos_mp = sum(false_positive_vf_mp)

# ╔═╡ a2acf8b3-f56f-4144-8990-bf2058182cde
function sensitivity_mpecificity_mp()
    f = Figure(resolution=(800, 800))
    colors = Makie.wong_colors()

    ##-- A --##
    ax = Axis(
		f[1, 1]; 
		xticks = (1:4, ["Integrated", "Volume Fraction", "Spatially Weighted", "Agatston"]),
		title = "False-Negative (Medium Patient)",
		ylabel = "False-Negative (%)",
		yticks = [0, 25, 50, 75, 100]
	)

    table = [1, 2, 3, 4]
	h1 = (total_zero_i_mp / total_cac_mp) * 100
	h2 = (total_zero_vf_mp / total_cac_mp) * 100
    h3 = (total_zero_s_mp / total_cac_mp) * 100
    h4 = (num_zero_a_mp / total_cac_mp) * 100
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
		title = "False-Positive (Medium Patient)",
		ylabel = "False-Positive (%)",
		yticks = [0, 25, 50, 75, 100]
	)

    table = [1, 2, 3, 4]
	h1 = (total_zero_i_pos_mp / total_cac_pos_mp) * 100
	h2 = (total_zero_vf_pos_mp / total_cac_pos_mp) * 100
    h3 = (total_zero_s_pos_mp / total_cac_pos_mp) * 100
    h4 = (total_zero_a_pos_mp / total_cac_pos_mp) * 100
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

    # save(joinpath(dirname(pwd()),"figures", FIGURE_PATH, "sensitivity_mpecificity.png"), f)


    return f
end

# ╔═╡ cb2f16bc-c027-486e-a91f-7e730a05a71c
with_theme(sensitivity_mpecificity_mp, medphys_theme)

# ╔═╡ 8d1942eb-72bf-40a5-a92e-dfb42b058ac1
md"""
#### Large Patient
"""

# ╔═╡ a16bc287-e8ba-4dec-a67f-653452fd002c
md"""
##### False Negative
"""

# ╔═╡ b0264c59-f68c-4a59-a5bb-de5893d53d9d
md"""
###### SWCS
"""

# ╔═╡ 8e217351-0bc6-4783-85b6-2016b2b14a70
begin
	false_negative_s_lp = []
	for i in 1:3:nrow(df_s_large_patient)-2
		mean_s, std_s = mean(df_s_large_patient[i:i+2, :swcs_bkg]), std(df_s_large_patient[i:i+2, :swcs_bkg])*std_level
		array_s = hcat(df_s_large_patient[i:i+2, :calculated_swcs_large], df_s_large_patient[i:i+2, :calculated_swcs_large], df_s_large_patient[i:i+2, :calculated_swcs_small]);
		neg = length(findall(x -> x <= (mean_s + std_s), array_s))
		push!(false_negative_s_lp, neg)
	end
end

# ╔═╡ 6a738641-d188-4932-8156-ce43381e94b3
total_zero_s_lp = sum(false_negative_s_lp)

# ╔═╡ 627a7536-ef67-4b66-8c2d-75707a7c130d
md"""
###### Agatston
"""

# ╔═╡ 853d5979-fcb3-4061-8812-96ea3b07feb4
num_zero_a_lp = length(findall(x -> x <= 0, array_a_lp))

# ╔═╡ acea136f-7c9f-4df4-ad7a-61a8dd305c8b
total_cac_lp = length(array_a_lp)

# ╔═╡ 9f3c31b6-cfe2-4db3-84c1-0c6ba279396b
md"""
###### Integrated
"""

# ╔═╡ 214897cc-b61b-4846-8823-e2d854a6daf2
begin
	false_negative_i_lp = []
	for i in 1:3:nrow(df_i_large_patient)-2
		mean_i, std_i = mean(df_i_large_patient[i:i+2, :mass_bkg]), std(df_i_large_patient[i:i+2, :mass_bkg])*std_level 
		array_i = hcat(df_i_large_patient[i:i+2, :calculated_mass_large], df_i_large_patient[i:i+2, :calculated_mass_large], df_i_large_patient[i:i+2, :calculated_mass_small]);
		neg = length(findall(x -> x <= mean_i + std_i, array_i))
		push!(false_negative_i_lp, neg)
	end
end

# ╔═╡ 7f68a07c-f41b-4011-88ab-dcfb98d17e7d
total_zero_i_lp = sum(false_negative_i_lp)

# ╔═╡ a2c2ee38-16de-42c1-ad39-0c9683240eaa
md"""
###### Volume Fraction
"""

# ╔═╡ 21be7405-ec78-4d62-8e13-2f887dd137c8
begin
	false_negative_vf_lp = []
	for i in 1:3:nrow(df_vf_large_patient)-2
		mean_vf, std_vf = mean(df_vf_large_patient[i:i+2, :mass_bkg]), std(df_vf_large_patient[i:i+2, :mass_bkg])*std_level 
		array_vf = hcat(df_vf_large_patient[i:i+2, :calculated_mass_large], df_vf_large_patient[i:i+2, :calculated_mass_large], df_vf_large_patient[i:i+2, :calculated_mass_small]);
		neg = length(findall(x -> x <= mean_vf + std_vf, array_vf))
		push!(false_negative_vf_lp, neg)
	end
end

# ╔═╡ acd59c06-8b77-4d7e-8cd6-c12cc9930a6c
total_zero_vf_lp = sum(false_negative_vf_lp)

# ╔═╡ 89c3bc5e-fd37-40c6-866b-be9ff5be5018
md"""
##### False Positive
"""

# ╔═╡ 1bfb8f31-bddd-4d7c-96e0-40d9123b07a1
md"""
###### SWCS
"""

# ╔═╡ 31664e97-6f72-493d-833e-2eaf68910cc0
begin
	false_positive_s_lp = []
	for i in 1:3:nrow(df_s_large_patient)-2
		mean_s, std_s = mean(df_s_large_patient[i:i+2, :swcs_bkg]), std(df_s_large_patient[i:i+2, :swcs_bkg])*std_level
		array_s_pos = df_s_large_patient[i:i+2, :swcs_bkg]
		pos = length(findall(x -> x > (mean_s + std_s), array_s_pos))
		push!(false_positive_s_lp, pos)
	end
end

# ╔═╡ 8657d123-8ec3-4852-962e-fd2d3d842e61
total_zero_s_pos_lp = sum(false_positive_s_lp)

# ╔═╡ be90ad05-e241-4f0a-acbe-dd647a8a6730
md"""
###### Agatston
"""

# ╔═╡ 3726789b-c662-44c9-8a64-f5dbf8c88717
array_a_pos_lp = df_a_large_patient[!, :mass_bkg];

# ╔═╡ 025615c0-edeb-4f44-a705-2ae799c87c01
total_cac_pos_lp = length(array_a_pos_lp)

# ╔═╡ 8fc458f5-6e45-4944-827c-689ab17193ad
total_zero_a_pos_lp = length(findall(x -> x > 0, array_a_pos_lp))

# ╔═╡ 9c0373d6-53ce-4169-9a6f-942a690272d3
md"""
###### Integrated
"""

# ╔═╡ 66b4b908-3564-49e9-b556-708a9f23a80e
begin
	false_positive_i_lp = []
	for i in 1:3:nrow(df_i_large_patient)-2
		mean_i, std_i = mean(df_i_large_patient[i:i+2, :mass_bkg]), std(df_i_large_patient[i:i+2, :mass_bkg])*std_level
		array_i_pos = df_i_large_patient[i:i+2, :mass_bkg]
		pos = length(findall(x -> x > (mean_i + std_level), array_i_pos))
		push!(false_positive_i_lp, pos)
	end
end

# ╔═╡ 6eac382e-4a29-4b50-951e-4490fc52f97c
total_zero_i_pos_lp = sum(false_positive_i_lp)

# ╔═╡ faf26dba-bc46-4d48-aa3f-c27f6083706a
md"""
###### Volume Fraction
"""

# ╔═╡ 6eb07ce0-7f51-4cc1-8d69-50a13228e570
begin
	false_positive_vf_lp = []
	for i in 1:3:nrow(df_vf_large_patient)-2
		mean_vf, std_vf = mean(df_vf_large_patient[i:i+2, :mass_bkg]), std(df_vf_large_patient[i:i+2, :mass_bkg])*std_level
		array_vf_pos = df_vf_large_patient[i:i+2, :mass_bkg]
		pos = length(findall(x -> x > (mean_vf + std_vf), array_vf_pos))
		push!(false_positive_vf_lp, pos)
	end
end

# ╔═╡ abe62336-29db-43e6-a0bf-1c42c21a3483
total_zero_vf_pos_lp = sum(false_positive_vf_lp)

# ╔═╡ 83755b0f-3e5f-45a7-89cc-a232ac8ee6d6
function sensitivity_lpecificity_lp()
    f = Figure(resolution=(800, 800))
    colors = Makie.wong_colors()

    ##-- A --##
    ax = Axis(
		f[1, 1]; 
		xticks = (1:4, ["Integrated", "Volume Fraction", "Spatially Weighted", "Agatston"]),
		title = "False-Negative (Large Patient)",
		ylabel = "False-Negative (%)",
		yticks = [0, 25, 50, 75, 100]
	)

    table = [1, 2, 3, 4]
	h1 = (total_zero_i_lp / total_cac_lp) * 100
	h2 = (total_zero_vf_lp / total_cac_lp) * 100
    h3 = (total_zero_s_lp / total_cac_lp) * 100
    h4 = (num_zero_a_lp / total_cac_lp) * 100
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
		title = "False-Positive (Large Patient)",
		ylabel = "False-Positive (%)",
		yticks = [0, 25, 50, 75, 100]
	)

    table = [1, 2, 3, 4]
	h1 = (total_zero_i_pos_lp / total_cac_pos_lp) * 100
	h2 = (total_zero_vf_pos_lp / total_cac_pos_lp) * 100
    h3 = (total_zero_s_pos_lp / total_cac_pos_lp) * 100
    h4 = (total_zero_a_pos_lp / total_cac_pos_lp) * 100
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

    # save(joinpath(dirname(pwd()),"figures", FIGURE_PATH, "sensitivity_lpecificity.png"), f)


    return f
end

# ╔═╡ 8233cfc2-0cff-44e2-8503-aafb0008d700
with_theme(sensitivity_lpecificity_lp, medphys_theme)

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
# ╠═f8e0735d-7854-48d0-85a9-ccfb3f8fc898
# ╠═44341e09-1593-4d24-8c77-1e3fceb219ee
# ╠═766474b7-8158-48b5-a4fb-6ed4ab055338
# ╠═3cf4a3a1-0997-425c-9f4f-0785a45eea1c
# ╟─aaa022d9-2d23-4bc3-96fc-dc349848820e
# ╟─decd4320-112a-4c4c-9ed5-680a721b81c1
# ╟─0fa91a0b-c205-4b10-9876-9a9e163f3d7e
# ╟─bc9939d7-f3a8-48e5-ada4-8fcf0fc1137a
# ╟─3f02194d-9c69-45d1-80fc-4cc50b1f3d8f
# ╟─e3c5f720-adec-4b59-b8a7-a26385c002ff
# ╟─0822b8d6-6285-4b69-80e1-23538044e6cc
# ╟─9b3c6bda-5c21-4c07-be6f-939a1ac0496a
# ╠═b4e25c03-ac59-4104-acf5-2f049497d460
# ╠═074a6437-9745-4680-b9b5-baf9d767afc9
# ╠═f51b5a1e-a688-470d-8a83-b65f03fff875
# ╟─91a157fb-a515-4ab0-be1a-8d9148db6d67
# ╠═1478eef5-bcde-4513-8794-9b9cbe388ee2
# ╠═ffeb8498-d7a3-4122-b71f-56326490042b
# ╠═ad784be1-5252-4914-bdd2-085edf5a95e0
# ╟─77f50463-59bc-477d-8a26-97aa08b0e21f
# ╟─971ba7ff-c14a-4b24-8e3c-5df3f68a05fd
# ╟─89ace836-d183-4872-8037-c48a5db02b89
# ╟─2555e6e0-2712-4ba0-a1f8-38a8e89b025b
# ╠═086bdda7-cdf5-4510-95e9-c405a2acb9d2
# ╟─5b2b1e40-d75e-4281-9ddd-4445775b79ba
# ╟─c41d53f9-a74b-4096-8709-5a92c2bc7eaa
# ╠═8eac91f6-7872-4065-bc78-6ed611b2b192
# ╠═46a5b2bd-3f00-4e3a-b799-f1a9f1041928
# ╠═7dd8bbb5-80dc-490d-b518-d7b14452e305
# ╟─4cbeb9d2-4ed1-4190-8fd0-aeff30f3c81d
# ╟─d5a323d0-5ca7-46f5-b169-1e17c4ec0413
# ╠═398e4d28-d01b-4bd7-8145-7d78a3de672f
# ╠═a6f72994-a3e1-4bb0-963f-12f8c2788ad8
# ╠═134fc008-1eb8-496c-a079-0b6cf177a3d1
# ╟─7ae5193f-2d29-4152-b062-791629a5f1fa
# ╟─a28a1259-9161-4c18-82f2-9bfb26a50690
# ╠═6e2f2819-6421-416f-a097-df1faefd9748
# ╠═e6f08aec-06e2-4dd4-86c6-0e87db2d6f88
# ╠═384f8f09-6df0-4fd8-91a5-f901001afea2
# ╟─d4753644-29fb-41a9-b8d1-5f46bb87a613
# ╠═aef525b3-be1a-4a80-ba18-7ade91baea2d
# ╠═0e5b9d7d-aa2f-44fa-9ca9-1d8b7f373f7b
# ╠═dfcd85cd-586c-49aa-96e2-2e265fb5f33d
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
# ╠═598e1797-837f-4bef-8fac-716f0690142a
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
# ╟─473918d5-df19-439d-a125-e0241495d29d
# ╟─8af3d98c-f33c-4cb2-b388-387bf2385530
# ╟─b3c5c174-24be-447f-90f7-31a6252a55ae
# ╟─fc3c711f-dee1-4e88-9337-9d930a40e3ea
# ╟─b7c5e3f3-d0a5-4501-b260-73960c00f160
# ╟─a50a0ed7-5733-44ec-857b-710156ce40e3
# ╠═4877368a-2491-4f45-8d52-5050e335ad4d
# ╠═6720e3d7-6e6a-4a7e-b749-f4777cba587c
# ╠═a5409851-8bc1-434a-be9f-d6a5bd8e392a
# ╠═1309668f-d000-4062-a582-a90e23ef43f8
# ╟─768e6d0a-9c75-4b2f-b52a-bb68b834f68e
# ╠═ac366a6b-8384-4735-9aa5-51c7f7604877
# ╠═3802ae31-ae12-45ba-b4b5-2854729b1ff0
# ╠═ae19ce13-bba7-41cb-8be2-2e56ce7fa879
# ╠═2011184c-54ec-4243-a42b-af2a6eee5552
# ╟─8ac18600-6073-4f2f-8709-43101f264a22
# ╠═9c4fc2f5-3b58-4f8a-a681-e246ac757dc2
# ╠═cb402344-d304-4ec8-9fe4-2bda0a0d7bfd
# ╠═f40d6a0e-384d-4d61-af93-e7398046d364
# ╠═022674ae-64de-4e3a-b70d-89ef428cc780
# ╟─16812354-9c81-48c7-83c0-97c0c7c6963f
# ╠═dfe97497-f0b5-4cce-a6b8-b68382a1142a
# ╠═832051a2-8290-49eb-be5a-cfba36f19238
# ╠═9492ef06-4986-415b-9bec-cf3a84d686c6
# ╠═d3f7af8b-4a3a-4d1d-81ff-e529bd8994c5
# ╟─0c5c5a7b-3aa4-4193-b54a-7c753b819a60
# ╟─3d7de16f-5a40-4737-af20-7bdd535239d7
# ╟─1eea3eb6-8006-4a9f-adc2-00e19aec2af9
# ╠═a72aed45-966b-48fd-a1f5-ae39c1f1ca05
# ╠═6e5c06d3-7b5e-4aae-985d-e7458a3cb263
# ╠═8c68de3d-8126-47e4-bdbd-982ed0244b65
# ╠═3863e320-368c-4816-82d1-eeccc8ed004d
# ╟─d3278fe2-e57b-4983-a6eb-8dead9a70c94
# ╠═d205dd2f-b1a2-4285-ae27-31ac9cb5bfe2
# ╠═967baf16-5b6b-4611-9703-30990bf05852
# ╠═b83a626e-d45d-40ba-89ef-4787d44332ad
# ╠═7969e05f-5148-4d1b-ac77-f8b9dde09e64
# ╟─ef9255e1-3300-4354-8c9f-6226e3f32145
# ╠═c77c6d5c-a989-42fd-99ef-4adbf6aa6db5
# ╠═adcc69bc-c067-49cc-a017-5dc11ad92063
# ╠═ed7cad03-0e74-441f-b5d2-2426c0cf5767
# ╠═27abe3ed-eaa6-4acf-b83f-56a0b7838501
# ╟─8b39efb6-b0dc-418c-9e1f-7bf7231db159
# ╠═a2d4a55c-da21-43b8-aec3-d80b5e29bc04
# ╠═919efe63-7c77-45b8-b8d6-ed7f6e0bf694
# ╠═51bc03d9-6821-47cd-9435-1d79aa37f779
# ╠═9848ed54-530c-46c9-b47c-5214c401c6b2
# ╟─87c5676f-3ed2-49fd-ae79-aee072b2bc80
# ╟─12869712-c5f0-4f31-8088-e885b5d9ddb9
# ╟─c4c1495f-3366-4f39-be4c-667d5895bf17
# ╟─7684a47b-a26d-4ab1-bc9d-54b4513c258a
# ╠═5b853ca7-f2e0-4d38-87b8-e48faa9102c2
# ╟─bd05fce5-a7d7-4e85-afa3-995224183406
# ╟─f46a6aac-eaee-4016-b628-b31c02be1798
# ╟─b3cc1fa7-c6d4-4002-bada-ea882f740e8b
# ╠═a82a6cf4-9599-419b-8ffc-56dda90954fb
# ╠═f8d840c1-ec0b-4165-ae75-79cd1e0a4fd6
# ╠═38605356-7165-44df-95c6-532b411ca708
# ╠═79c52fd3-8b91-4982-8c5e-2b212c0d4297
# ╟─d36e7476-a540-46a0-86de-88ad0c78036d
# ╠═76c1dc7e-e1fd-481a-bdff-f3b17b32f584
# ╠═96f9e011-c879-4f39-801e-ca7469c6ef0a
# ╠═a9a49218-1239-4caf-9774-3a5614ab40be
# ╠═66ad4a90-27bd-4857-b0be-09b6fa25741a
# ╟─4217ba04-b8b9-425f-8991-10a058621478
# ╠═b641aca4-af92-4110-8cd4-2ba4c122df04
# ╠═eb4d0f50-65e4-44ba-9eb6-88f84ed4b070
# ╠═ba357e75-bcc8-47db-9a11-6baf050ae925
# ╠═c767a4a8-a521-4d68-990b-50c177788193
# ╟─b4fe4133-256f-4e33-b4e5-6fa38a126c4f
# ╠═74ea613c-e4c4-4124-ae48-8ec6e1df68ff
# ╠═6c7cdb74-cae8-4842-872c-65afc4ddd3a6
# ╠═d5b1182a-fbdb-4422-abda-6df0a5490003
# ╠═5f2a36aa-3b66-441d-a909-14f6d2820f2c
# ╟─aa0487c3-2fe8-49f4-bf4d-927ce73de052
# ╟─457f79fe-4884-4bb8-b6ef-6da088cefc76
# ╟─03bf4424-01f7-4d50-aae2-fc4c6c3a19ba
# ╠═d4e37208-0eab-4091-b360-94f74c9bfbd7
# ╠═9bd4a8d3-c641-4fb5-b0bb-f6d699e944c4
# ╠═a09c2653-21b8-4a17-8b02-2451cda66ac4
# ╠═cbb80103-63ac-40c0-8eae-d24c058695a5
# ╟─6af52148-3473-4660-b6cc-aa9879ee46cb
# ╠═fc88c296-38f3-4ccf-a1d9-73bd40e26b40
# ╠═836f13b7-4aca-4050-b6ee-4522a45738ec
# ╠═4181e1d4-503d-4dfe-ae62-be3e50af3d9f
# ╠═e019181d-639b-474d-bd48-9e714502cd23
# ╟─889c5c63-146f-4b77-b067-3d4c6f977bc5
# ╠═ed14543e-8444-4c00-aa32-e05d5f34bcf3
# ╠═77ecb385-24c3-45b9-857b-3c3a1d5480a1
# ╠═9358b9c7-5dfa-4dda-ab09-e522c5fcb202
# ╠═6eda4c2d-b56c-4366-a4f3-0e0dd710e061
# ╟─a0da5d2f-835e-4250-8841-2d22b6531df0
# ╠═10dda852-aa33-49a6-adfc-08b95cfe70b5
# ╠═1e7b01d3-2077-4f22-9c14-02231a5bb8cb
# ╠═4e0b3ef7-5a75-46da-9797-94ad80005bec
# ╠═39260980-b87b-442e-82a9-532f3854b607
# ╟─2afe38f3-f10a-47e3-a3e7-32b4523e17c5
# ╟─b65ed831-7a51-4035-b28a-24ae9ec64a01
# ╟─83d13107-ad71-453e-b065-f009d7bfdc55
# ╠═4b4a25b0-2973-49cb-9be4-265e00b9417f
# ╠═de7d4719-5337-43c1-a9d7-3d32f2560388
# ╠═227f49dc-e9fc-4560-b406-ebc959d9328f
# ╠═007e7280-0cf6-45c4-b6d6-240a901eac86
# ╟─3963efc0-8c06-4672-9c6e-3ae256e10d4d
# ╠═a832370e-af82-444a-9428-7114649ba203
# ╠═2969e17f-d3f5-4438-ac4a-05370df59a44
# ╠═c4216708-b4b7-4625-831c-d7b9cce57ef9
# ╠═32ae56f6-2e86-4c6b-9343-8535a975da2b
# ╟─ea5a016d-bed9-4f00-9960-5305e548fa28
# ╠═39132f1c-7a92-47b5-9cd7-2ee9b61cdf89
# ╠═13a516b5-6cb6-4384-80a6-78ffe7078cfc
# ╠═76194267-a83f-4aac-90a2-5e724ca1df60
# ╠═68a7b68b-0d92-4c7c-ad12-342d1c54efc8
# ╟─2dc008cf-5570-4a7e-a046-56ce1880b56e
# ╠═350a77af-9f11-4aea-aae5-31e30e1ccdc2
# ╠═004bedd1-8237-4271-ab70-cd412d444027
# ╠═fb79376e-26a9-49e2-a6ba-301596d7e609
# ╠═221689f5-ac9a-4a69-a366-615aec7513b2
# ╟─c50117bb-bbb7-4347-b4bb-d446677f4509
# ╟─34991466-181f-4036-8f47-12e20f2144bc
# ╟─bf35fb67-20b9-4b74-bd15-3315d5e37855
# ╠═75e0f9d7-ef66-42bd-bf06-cebecb09a526
# ╠═5e32f34c-2f32-457e-a9cf-59184170b6ed
# ╠═819c99a9-5a22-447b-9136-0fb3372af1e4
# ╟─9ed516b4-232b-45c5-af90-2964ed2e6a36
# ╟─10725f84-8ce4-454a-8c6e-41cdaa25dd1a
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
# ╠═beaa94a0-0de3-40a3-86a1-912cbc018775
# ╠═4eef56ad-0ea4-422c-837b-3aa9846a2458
# ╠═b708d2c9-f868-4814-b1b7-6908317aa5ba
# ╟─1d04d29f-401d-45c4-99ff-3bbdb0fb4f6b
# ╟─559f77d3-69c0-4155-8cc0-1f669f929bc3
# ╠═b540e5f0-d271-4d42-a167-bf79766a53ac
# ╠═802dc604-cd58-40db-ad6c-132ad38d7c3c
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
# ╟─ac33fc36-b7c2-47f1-b513-8a02f5fdd0ad
# ╟─ce3956ea-9947-4bc9-91df-39079eec40ee
# ╟─6e0c8aa0-89ac-4e5b-b5e1-f83ad87c246a
# ╟─8aae5583-bf94-42f2-a7b7-b6b9a3a3c3eb
# ╟─e4c021dc-2fee-401f-8bb8-4f72f0186231
# ╟─a1f3606e-69c9-45c6-9602-4367ca63ae48
# ╟─32fb4eca-3ef6-4254-8fcf-79150c8ae39d
# ╠═a894fc59-25b4-4536-bac6-eb2352c5baaf
# ╠═aeb44f0c-c77d-4f5e-8ad1-dbc3b7a0d7af
# ╟─e0937eba-d445-4b73-a266-663e56c33246
# ╠═e166f587-2a51-4135-86ce-648a5f92c0f1
# ╠═39f85f44-558d-4f55-9ed0-1db50e309847
# ╟─47d4c62d-6d01-4999-91ca-b0b9679a9bb5
# ╠═07a7060b-5490-4e1a-9739-eb9bf790a999
# ╠═7003f82d-d6a7-41b4-98e2-0cbefc4f1e24
# ╟─859ba958-96d5-4e26-9c98-963b1a3bcddf
# ╠═b780e9ea-4b9f-4ad9-b6dc-6c9b5070863d
# ╠═7734498b-ea86-47f7-8f17-c842a58e9c93
# ╟─fa11f18f-3e05-435c-ab18-fa3faba49871
# ╟─2960b18c-f0df-40b6-8284-8d843ccd3783
# ╠═62d4be9f-046d-4b81-946a-8ddf359dc179
# ╠═295c146f-ffa0-4d62-9352-0b21c49584a1
# ╟─7efe886e-76ca-4269-b01f-c87fd678032a
# ╠═fbcdaf78-f366-428d-ad63-e0b17cfc1fa4
# ╠═ad1303e7-a9eb-4655-a176-4a54b77db1b7
# ╠═ad25959e-0cda-4292-8052-070ba4d73b59
# ╟─a94b4849-c2c9-4a22-b238-6212c9867f9f
# ╠═c389a74e-4ab1-48f4-ac37-fcc0b794affb
# ╠═36b574e8-6eba-4038-b889-d472c4fdc824
# ╟─a60d6f3e-7c31-47f9-901f-4b110d572b3b
# ╠═0f60944f-e615-4117-bde6-9efaa9900806
# ╠═f1434cec-5655-46fe-af66-8eea8f639d94
# ╟─0eff2135-b37f-413a-add6-80a07197d346
# ╟─ab20b172-35c7-49d6-8b57-ae8c2a310a26
# ╟─8d88faa9-df59-4a4c-ad03-29612c2255a9
# ╟─2a2de10a-4792-45bd-ab03-9679a46ff6cd
# ╠═fb2828fa-0c4f-409e-b739-cb1f2c135f8f
# ╠═d5ec6fbc-1577-4950-87e6-51a7e1b7da75
# ╟─ed9e7420-0fb8-4882-8cd8-f38b1b254684
# ╠═0d1900e6-3681-4f20-afb6-63ccd92bd388
# ╠═f266aca7-533e-48ac-a104-8f2940b97eb3
# ╟─9c2a5fda-d887-4101-9e37-1ef6ed71c4aa
# ╠═8f11feed-fa28-405d-a8c2-105e15b065d3
# ╠═a77cc0a3-6eec-4bb3-aa25-299d8051dead
# ╟─b8c525da-2805-4118-b116-10aeaad4cba5
# ╠═310cb551-ca38-49d7-b570-38eb5221f2bb
# ╠═521a8a99-d5b1-4d38-8365-1cfbf64a1bb0
# ╟─1d65d6b1-16f4-4845-beb5-49d307be8a8b
# ╟─348bb217-e3b1-43ce-86f2-76629d1ca00e
# ╠═2a2492cd-8533-4149-bbe3-d1bdb70f22e2
# ╠═75a7f90b-1ccd-4770-a90b-ac7a5953443c
# ╟─9eec9ed3-46a7-4c7e-9f91-9ea0b01afb9f
# ╠═08a34dcc-b4ca-434e-b9df-38666dfcbdc0
# ╠═4fbfcb02-d325-4d75-86b5-bd944352bca6
# ╠═1b5074bf-b036-42f1-9853-42b1e53865c7
# ╟─7f7ff260-ae32-4d78-9e85-c268ced39cba
# ╠═bc9eb230-7ff3-4597-8b70-1a6e58918afb
# ╠═ac924022-6f03-4566-9c0a-23c3ae278938
# ╟─d6392619-9f11-4824-8168-c74fd08eece5
# ╠═6950ed5a-fb5d-4b9e-a17e-314760933ed0
# ╠═632e3310-7c7a-4d91-a2cd-d00fc2633492
# ╟─91f9ec30-e6e8-4bc7-83c3-b49c98b0b32d
# ╟─67dcee1b-f488-431e-8c1d-4211bf04696f
# ╠═cb2f16bc-c027-486e-a91f-7e730a05a71c
# ╟─8233cfc2-0cff-44e2-8503-aafb0008d700
# ╟─15dea22a-662b-43c4-ac6b-33ea06a5f2c4
# ╟─32cc0ed0-c025-43db-a107-83af5a8f1571
# ╟─a6b18ed3-080c-4691-bd53-83b878cbec74
# ╟─0270d409-cb8f-45b2-b30c-cdd398e7fe5f
# ╠═649ca894-ab4b-4a60-96f7-c91be28a092b
# ╠═d86480a0-041e-402b-a56d-f10a9903529f
# ╟─85dc1af9-033d-4371-b0e9-2888d1478bb8
# ╠═722824a9-e6b3-4a63-8131-f2043f7269f1
# ╠═04e8dfd9-2476-4a54-ad4c-cc381fad3087
# ╟─e2a1b07c-1b14-464c-8589-2b3fd8d55d00
# ╠═e35098ae-4147-41a4-b4ff-dfa6dcc7ff2f
# ╠═4bee3c44-f9c8-4428-b6d2-471ec732f3d1
# ╟─2c888b16-5385-4ea1-a5b1-e8a5a08c89c9
# ╠═1ad040bb-e0a4-421b-b80b-90397db8121b
# ╠═0fb3dedd-782b-4f2f-add6-041c72b51d12
# ╟─cfe68dce-c66c-4c98-9bc9-340859d34e82
# ╟─8e68a1fd-a9d7-4bd9-a863-0543029d9a90
# ╠═6d16a3ff-4de7-4fab-9c37-05497b2c895e
# ╠═6bee3d94-ce21-4e54-82a1-2fa926096d40
# ╟─3fae6b4d-4ffe-47b2-92e9-a4e9f509759e
# ╠═8019a1c3-650b-40c2-9bf7-cc8d245f1173
# ╠═7f9119e0-a71d-49d1-a441-037dbb5be4a2
# ╠═a4f44a8c-70f1-40ff-b87f-9ab1082bf435
# ╟─4dfe6af4-21c3-4c13-8c0b-dc2ae437129b
# ╠═22d33469-000e-46d3-8693-f71a4682180f
# ╠═176af30d-9078-43a6-a5f8-5b93d536349a
# ╟─48684053-4258-4dd7-94d5-bbc73196ea10
# ╠═5e96e6b5-07cd-4421-bdcc-ad1524a29b62
# ╠═4db669a2-c259-4dfa-aa19-3862da5380d9
# ╟─d55054cd-c63b-43c8-b0f9-5f4e60b6415d
# ╟─a2acf8b3-f56f-4144-8990-bf2058182cde
# ╟─0ebdfc89-960e-45cc-a129-54a003d74d42
# ╟─9a37db8b-7e94-45ce-b459-9c10246270d7
# ╠═d57d99ae-4aae-4081-87f5-31ff8d3e85c3
# ╠═7c9c9b09-74e4-427e-af9f-f047b41d3ff9
# ╟─986bfd16-6520-43b9-912a-09f842e14724
# ╠═ca0c55ac-6e17-4219-9fc6-0b2288fb7b1c
# ╠═2440eb1b-681c-47e0-b14b-1fdc201d3577
# ╟─f08d6966-7132-4c98-866c-b1ea531acc9e
# ╠═b3015cee-3cfa-4968-bd8d-8bcd0a961278
# ╠═34c20a4e-3ec2-45ae-b83a-b3f4f9b0ee46
# ╟─bf25092a-f6ca-4605-994a-a13bbb5b41d3
# ╠═dbfc643d-81cd-4c99-8665-79c007115119
# ╠═5e94d82c-05ab-4bef-a5bc-d4cb5791a338
# ╟─5f6d3a21-96f1-4618-8de4-3add807572cb
# ╟─9b207e45-1ef3-4742-be7d-dd8f67abea1b
# ╠═a07d7c66-c17c-4f17-8cc2-7db04d93a696
# ╠═3f0ee19a-55ea-4d75-be27-0a07bb33d909
# ╟─2392ea54-a176-4d70-a144-06f31c28cda1
# ╠═93f8369f-f26f-4736-b86d-b306e5c53dce
# ╠═6b0f51bf-b2be-4d5a-8ca1-733dc6a16628
# ╠═4b341998-336c-4ec4-a1b4-08e35c155f4c
# ╟─8a0d81d1-863a-4899-ad19-fcbd5f87d857
# ╠═dc9f382d-f05d-4ba6-8a96-cc963ed4722d
# ╠═f2cc4db1-65ce-4ebe-a10f-82cd9de5bef7
# ╟─89e0e212-8ddd-45eb-8054-2045fa3e29f2
# ╠═6cd86ec7-3e1a-4e7a-9f61-4077ca986fd2
# ╠═4dceb05b-f463-4789-a79a-3cb9888b6c2e
# ╟─8d1942eb-72bf-40a5-a92e-dfb42b058ac1
# ╟─83755b0f-3e5f-45a7-89cc-a232ac8ee6d6
# ╟─a16bc287-e8ba-4dec-a67f-653452fd002c
# ╟─b0264c59-f68c-4a59-a5bb-de5893d53d9d
# ╠═8e217351-0bc6-4783-85b6-2016b2b14a70
# ╠═6a738641-d188-4932-8156-ce43381e94b3
# ╟─627a7536-ef67-4b66-8c2d-75707a7c130d
# ╠═853d5979-fcb3-4061-8812-96ea3b07feb4
# ╠═acea136f-7c9f-4df4-ad7a-61a8dd305c8b
# ╟─9f3c31b6-cfe2-4db3-84c1-0c6ba279396b
# ╠═214897cc-b61b-4846-8823-e2d854a6daf2
# ╠═7f68a07c-f41b-4011-88ab-dcfb98d17e7d
# ╟─a2c2ee38-16de-42c1-ad39-0c9683240eaa
# ╠═21be7405-ec78-4d62-8e13-2f887dd137c8
# ╠═acd59c06-8b77-4d7e-8cd6-c12cc9930a6c
# ╟─89c3bc5e-fd37-40c6-866b-be9ff5be5018
# ╟─1bfb8f31-bddd-4d7c-96e0-40d9123b07a1
# ╠═31664e97-6f72-493d-833e-2eaf68910cc0
# ╠═8657d123-8ec3-4852-962e-fd2d3d842e61
# ╟─be90ad05-e241-4f0a-acbe-dd647a8a6730
# ╠═3726789b-c662-44c9-8a64-f5dbf8c88717
# ╠═025615c0-edeb-4f44-a705-2ae799c87c01
# ╠═8fc458f5-6e45-4944-827c-689ab17193ad
# ╟─9c0373d6-53ce-4169-9a6f-942a690272d3
# ╠═66b4b908-3564-49e9-b556-708a9f23a80e
# ╠═6eac382e-4a29-4b50-951e-4490fc52f97c
# ╟─faf26dba-bc46-4d48-aa3f-c27f6083706a
# ╠═6eb07ce0-7f51-4cc1-8d69-50a13228e570
# ╠═abe62336-29db-43e6-a0bf-1c42c21a3483
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

function weighted_volume_fraction(vol, mu2, sigma2, voxel_size, ρ_calcium)
	d = Distributions.Normal(μ, σ)

    scaled_array = zeros(size(vol))
    for i in axes(vol, 1)
        for j in axes(vol, 2)
            for z in axes(vol, 3)
                scaled_array[i, j, z] = Distributions.cdf(d, vol[i, j, z])
            end
        end
    end

    weighted_arr = zeros(size(vol))
    for z in axes(scaled_array, 3)
        kern = [0 0.2 0; 0.2 0.2 0.2; 0 0.2 0]
        weighted_arr[:, :, z] = DSP.conv(scaled_array[:, :, z], kern)[2:end-1, 2:end-1]
    end

    return sum(weighted_arr) * voxel_size * ρ_calcium
end
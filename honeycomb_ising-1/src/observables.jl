module Observables

using Statistics

# 平均磁化を計算する関数
function mean_magnetization(spins)
    return mean(spins)
end

# スピン相関を計算する関数
function spin_correlation(spins, max_distance)
    L = size(spins, 1)
    correlations = zeros(Float64, max_distance + 1)
    
    for r in 0:max_distance
        sum_corr = 0.0
        for i in 1:(L - r), j in 1:L
            sum_corr += spins[i, j] * spins[i + r, j]
        end
        correlations[r + 1] = sum_corr / ((L - r) * L)
    end
    
    return correlations
end

end # module Observables
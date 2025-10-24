module Mypackage

using Random
using Statistics

# 固定端境界：端を +1 で固定する
# 内部は LxL 格子
function init_spins_fixed(L; ordered=false, rng=Random.GLOBAL_RNG)
    s = fill(Int8(1), L+2, L+2)
    if !ordered
        # 内部だけランダムに
        s[2:end-1, 2:end-1] .= Int8.(rand(rng, Bool, L, L)) .* 2 .- 1
    end
    return s
end

# 境界を考慮したエネルギー差 ΔE 計算
function local_energy_change_fixed(s, i, j; J=1.0, h=0.0)
    si = s[i,j]
    nb = s[i+1,j] + s[i-1,j] + s[i,j+1] + s[i,j-1]
    return 2.0 * si * (J * nb + h)
end

# メトロポリス更新（固定境界では端を更新しない）
function metropolis_sweep_fixed!(s, β; J=1.0, h=0.0, rng=Random.GLOBAL_RNG)
    L = size(s,1) - 2
    for _ in 1:(L^2)
        i = rand(rng, 2:(L+1))
        j = rand(rng, 2:(L+1))
        ΔE = local_energy_change_fixed(s, i, j; J=J, h=h)
        if ΔE <= 0.0 || rand(rng) < exp(-β * ΔE)
            s[i,j] = -s[i,j]
        end
    end
end

# 内部スピンの平均磁化
function mag_per_spin(s)
    L = size(s,1) - 2
    return mean(Float64.(s[2:end-1, 2:end-1]))
end

# スピン相関（x方向のみで簡易版）
function spin_correlation_fixed(s, maxr)
    L = size(s,1) - 2
    C = zeros(Float64, maxr+1)
    for r in 0:maxr
        sumval = 0.0
        for i in 2:(L+1-r), j in 2:(L+1)
            sumval += s[i,j] * s[i+r,j]
        end
        C[r+1] = sumval / ((L - r + 1) * L)
    end
    return C
end

# 実行関数
function run_ising_fixed(L, T; J=1.0, h=0.0, n_eq=1000, n_mc=5000, measure_interval=10, rng_seed=1234)
    rng = MersenneTwister(rng_seed)
    s = init_spins_fixed(L; ordered=false, rng=rng)
    β = 1.0 / T
    maxr = Int(floor(L/2))

    # 熱化
    for t in 1:n_eq
        metropolis_sweep_fixed!(s, β; J=J, h=h, rng=rng)
    end

    # 測定
    mags = Float64[]
    C_accum = zeros(Float64, maxr+1)
    n_measure = 0
    for t in 1:n_mc
        metropolis_sweep_fixed!(s, β; J=J, h=h, rng=rng)
        if t % measure_interval == 0
            push!(mags, mag_per_spin(s))
            C_accum .+= spin_correlation_fixed(s, maxr)
            n_measure += 1
        end
    end

    return mean(abs.(mags)), mean(abs.(C_accum ./ n_measure))
end

end # module Mypackage

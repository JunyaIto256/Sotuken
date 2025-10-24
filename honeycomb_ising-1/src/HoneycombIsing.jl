module HoneycombIsing

using Random
using Statistics

# ハニカム格子のスピンの初期化
function init_spins_honeycomb(L; ordered=false, rng=Random.GLOBAL_RNG)
    s = fill(Int8(1), L, L)
    if !ordered
        s .= Int8.(rand(rng, Bool, L, L)) .* 2 .- 1
    end
    return s
end

# 隣接スピンの取得
function get_neighbors(i, j, L)
    neighbors = []
    if i > 1
        push!(neighbors, (i-1, j))
    end
    if i < L
        push!(neighbors, (i+1, j))
    end
    if j > 1
        push!(neighbors, (i, j-1))
    end
    if j < L
        push!(neighbors, (i, j+1))
    end
    if i % 2 == 1
        if j < L
            push!(neighbors, (i+1, j+1))
        end
        if j > 1
            push!(neighbors, (i+1, j-1))
        end
    else
        if j < L
            push!(neighbors, (i-1, j+1))
        end
        if j > 1
            push!(neighbors, (i-1, j-1))
        end
    end
    return neighbors
end

# エネルギー差の計算
function local_energy_change_honeycomb(s, i, j; J=1.0, h=0.0)
    si = s[i,j]
    nb = sum(s[n...] for n in get_neighbors(i, j, size(s, 1)))
    return 2.0 * si * (J * nb + h)
end

# メトロポリス更新
function metropolis_sweep_honeycomb!(s, β; J=1.0, h=0.0, rng=Random.GLOBAL_RNG)
    L = size(s, 1)
    for _ in 1:(L^2)
        i = rand(rng, 1:L)
        j = rand(rng, 1:L)
        ΔE = local_energy_change_honeycomb(s, i, j; J=J, h=h)
        if ΔE <= 0.0 || rand(rng) < exp(-β * ΔE)
            s[i,j] = -s[i,j]
        end
    end
end

# 実行関数
function run_honeycomb_ising(L, T; J=1.0, h=0.0, n_eq=1000, n_mc=5000, measure_interval=10, rng_seed=1234)
    rng = MersenneTwister(rng_seed)
    s = init_spins_honeycomb(L; ordered=false, rng=rng)
    β = 1.0 / T

    # 熱化
    for t in 1:n_eq
        metropolis_sweep_honeycomb!(s, β; J=J, h=h, rng=rng)
    end

    # 測定
    mags = Float64[]
    for t in 1:n_mc
        metropolis_sweep_honeycomb!(s, β; J=J, h=h, rng=rng)
        if t % measure_interval == 0
            push!(mags, mean(abs.(s)))
        end
    end

    return mean(mags)
end

end # module HoneycombIsing
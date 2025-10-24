# ...existing code...
module HoneycombIsing

using Random
using Statistics

# src 内の補助ファイルを読み込む
include(joinpath(@__DIR__, "..", "src", "lattice.jl"))
using .Lattice

# --- ハニカム用の補助関数（簡易実装） ---
# 固定端: 周囲の境界（i==1 || i==L || j==1 || j==L）を +1 に固定する
function init_spins_fixed(L; ordered=false, rng=Random.GLOBAL_RNG)
    s = fill(Int8(1), L, L)
    if !ordered
        if L > 2
            s[2:end-1, 2:end-1] .= Int8.(rand(rng, Bool, L-2, L-2)) .* 2 .- 1
        end
    end
    # 境界は既に +1 になっている
    return s
end

# エネルギー差 ΔE（近傍は Lattice.get_neighbors または事前作成した neighbor_list を使用）
function local_energy_change_fixed(s, i, j; J=1.0, h=0.0, neighbor_list=nothing)
    si = s[i,j]
    L = size(s,1)
    nb_sum = 0
    if neighbor_list === nothing
        for (ni,nj) in Lattice.get_neighbors(i, j, L)
            nb_sum += s[ni,nj]
        end
    else
        idx = (i-1)*L + j
        for (ni,nj) in neighbor_list[idx]
            nb_sum += s[ni,nj]
        end
    end
    return 2.0 * si * (J * nb_sum + h)
end

# メトロポリス掃引（境界は更新しない）、neighbor_list をオプションで受け取る
function metropolis_sweep_fixed!(s, β; J=1.0, h=0.0, rng=Random.GLOBAL_RNG, neighbor_list=nothing)
    L = size(s,1)
    if L <= 2
        return
    end
    for _ in 1:((L-2)^2)
        i = rand(rng, 2:(L-1))
        j = rand(rng, 2:(L-1))
        ΔE = local_energy_change_fixed(s, i, j; J=J, h=h, neighbor_list=neighbor_list)
        if ΔE <= 0.0 || rand(rng) < exp(-β * ΔE)
            s[i,j] = -s[i,j]
        end
    end
end

# 内部スピンの平均磁化（境界を除く）
function mag_per_spin(s)
    L = size(s,1)
    if L <= 2
        return mean(Float64.(s))
    end
    return mean(Float64.(s[2:end-1, 2:end-1]))
end

# 簡易スピン相関（x方向のみ、境界考慮）
function spin_correlation_fixed(s, maxr)
    L = size(s,1)
    innerL = max(0, L-2)
    if innerL == 0
        return zeros(Float64, maxr+1)
    end
    C = zeros(Float64, maxr+1)
    for r in 0:maxr
        sumval = 0.0
        count = 0
        for i in 2:(L-1-r), j in 2:(L-1)
            sumval += s[i,j] * s[i+r,j]
            count += 1
        end
        C[r+1] = count > 0 ? sumval / count : 0.0
    end
    return C
end

# 実行関数で neighbor_list を作って渡すように変更、innerL に基づいて maxr を決定、n_measure ガード追加
function run_simulation(L, T; J=1.0, h=0.0, n_eq=1000, n_mc=5000, measure_interval=10, rng_seed=1234)
    rng = MersenneTwister(rng_seed)
    s = init_spins_fixed(L; ordered=false, rng=rng)
    β = 1.0 / T
    innerL = max(0, L-2)
    maxr = innerL > 0 ? Int(floor(innerL/2)) : 0

    # 近傍リストを事前作成（高速化）
    neighbor_list = Lattice.build_neighbor_list(L)

    # 熱化
    for t in 1:n_eq
        metropolis_sweep_fixed!(s, β; J=J, h=h, rng=rng, neighbor_list=neighbor_list)
    end

    # 測定
    mags = Float64[]
    C_accum = zeros(Float64, maxr+1)
    n_measure = 0
    for t in 1:n_mc
        metropolis_sweep_fixed!(s, β; J=J, h=h, rng=rng, neighbor_list=neighbor_list)
        if t % measure_interval == 0
            push!(mags, mag_per_spin(s))
            C_accum .+= spin_correlation_fixed(s, maxr)
            n_measure += 1
        end
    end

    mean_mag = isempty(mags) ? 0.0 : mean(abs.(mags))
    mean_corr = n_measure == 0 ? 0.0 : mean(abs.(C_accum ./ n_measure))

    return mean_mag, mean_corr
end

end # module HoneycombIsing

# メインスクリプト

using Plots

# パラメータ設定（必要に応じて変更）
L = 30                    # 格子のサイズ
J = 1.0                   # 交換定数
h = 0.0                   # 外部磁場
n_eq = 1000               # 熱化ステップ数
n_mc = 5000               # モンテカルロステップ数
measure_interval = 10     # 測定間隔
rng_seed = 1234           # 乱数シード（固定するならそのまま）

# 温度スキャン
Tlist = 0.1:0.1:2.0
mags_T = Float64[]
corrs_T = Float64[]

println("=== Temperature sweep ===")
println("L=$L, J=$J, T ∈ [$(minimum(Tlist)), $(maximum(Tlist))]")

for T in Tlist
    mean_mag, mean_corr = HoneycombIsing.run_simulation(L, T; J=J, h=h, n_eq=n_eq, n_mc=n_mc, measure_interval=measure_interval, rng_seed=rng_seed)
    push!(mags_T, mean_mag)
    push!(corrs_T, mean_corr)
    println("T=$(round(T,digits=2)) -> ⟨|m|⟩=$(round(mean_mag,digits=4)), ⟨|corr|⟩=$(round(mean_corr,digits=4))")
end

# プロット作成（左右に並べる）
p1 = plot(Tlist, mags_T, lw=2, marker=:circle, xlabel="Temperature T (J=1)", ylabel="⟨|m|⟩", title="Magnetization vs Temperature", legend=false)
p2 = plot(Tlist, corrs_T, lw=2, marker=:circle, xlabel="Temperature T (J=1)", ylabel="⟨|corr|⟩", title="Correlation vs Temperature", legend=false)

plot(p1, p2, layout=(1,2), size=(1000,400))

# 保存
output_file = joinpath(pwd(), "results_temperature_scan.png")
savefig(output_file)
println("プロットを保存しました: ", output_file)
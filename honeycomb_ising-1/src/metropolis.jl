module Metropolis

using Random

# メトロポリスアルゴリズムによるスピンの更新
function update_spin!(s, i, j, β; J=1.0, h=0.0, rng=Random.GLOBAL_RNG)
    ΔE = local_energy_change(s, i, j; J=J, h=h)
    if ΔE <= 0.0 || rand(rng) < exp(-β * ΔE)
        s[i,j] = -s[i,j]
    end
end

# 境界を考慮したエネルギー差 ΔE 計算
function local_energy_change(s, i, j; J=1.0, h=0.0)
    si = s[i,j]
    nb = s[i+1,j] + s[i-1,j] + s[i,j+1] + s[i,j-1]
    return 2.0 * si * (J * nb + h)
end

# メトロポリススイープ
function metropolis_sweep!(s, β; J=1.0, h=0.0, rng=Random.GLOBAL_RNG)
    L = size(s, 1)
    for _ in 1:(L^2)
        i = rand(rng, 1:L)
        j = rand(rng, 1:L)
        update_spin!(s, i, j, β; J=J, h=h, rng=rng)
    end
end

end # module Metropolis
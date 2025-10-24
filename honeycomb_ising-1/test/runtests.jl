module Main

using Test
using HoneycombIsing

# テスト用のスピン初期化関数
function test_init_spins()
    L = 4
    s = init_spins_fixed(L; ordered=false)
    @test size(s) == (L + 2, L + 2)
    @test all(s[2:end-1, 2:end-1] .== 1) || all(s[2:end-1, 2:end-1] .== -1)
end

# テスト用のエネルギー計算関数
function test_local_energy_change()
    L = 4
    s = init_spins_fixed(L; ordered=true)
    ΔE = local_energy_change_fixed(s, 2, 2)
    @test ΔE == 0.0  # ordered state should have zero energy change
end

# テスト用のメトロポリス更新関数
function test_metropolis_sweep()
    L = 4
    s = init_spins_fixed(L; ordered=false)
    original_s = copy(s)
    metropolis_sweep_fixed!(s, 1.0)
    @test s != original_s  # スピンが変更されているべき
end

# テスト用の平均磁化計算関数
function test_mag_per_spin()
    L = 4
    s = init_spins_fixed(L; ordered=true)
    mag = mag_per_spin(s)
    @test mag == 1.0  # ordered state should have magnetization of 1
end

# テスト用のスピン相関関数
function test_spin_correlation()
    L = 4
    s = init_spins_fixed(L; ordered=true)
    C = spin_correlation_fixed(s, 2)
    @test length(C) == 3  # maxr = 2 の場合、Cは3要素であるべき
end

# 全てのテストを実行
function run_tests()
    test_init_spins()
    test_local_energy_change()
    test_metropolis_sweep()
    test_mag_per_spin()
    test_spin_correlation()
    println("全てのテストが成功しました！")
end

run_tests()

end # module Main
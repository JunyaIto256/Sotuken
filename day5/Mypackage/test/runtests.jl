using Test
using Mypackage
using ITensors
import Revise

@testset "Tensor Operations" begin
    # テンソルの規格化
    T = random_itensor(Index(2, "i"), Index(2, "j"))
    T_normalized = Mypackage.normalize_tensor(T)
    @test isapprox(norm(T_normalized), 1.0, atol=1e-6)

    # 原始テンソルの生成
    W, P, C = Mypackage.tensor_origin_fix(1.0)
    @test norm(W) > 0
    @test norm(P) > 0
    @test norm(C) > 0
end

@testset "Tensor Decomposition and Compression" begin
    # 固有値分解
    T = random_itensor(Index(4, "i"), Index(4, "j"))
    V = random_itensor(Index(4, "k"), Index(2, "l"), Index(4, "m"))
    D, U = Mypackage.diagonalize_tensor(T)
    @test size(D) == (4, 4)
    @test size(U) == (4, 4)

    # 圧縮
    chi = 2
    C_compressed, P_compressed = Mypackage.compress_tensor(T, V, chi)
    @test size(C_compressed) == (chi, chi)
end

@testset "Physical Measurements" begin
    W, P, C = Mypackage.tensor_origin_fix(0.0)
    G, O = Mypackage.environment_tensor(W, P, C)

    # スピン相関の測定
    spin_corr = Mypackage.measure_spin_correlation(W, G, O)
    @test spin_corr isa Number
    @test isfinite(spin_corr)
    @test isapprox(spin_corr, 0.0, atol=1e-6)

    # 自発磁化の測定
    magnetization = Mypackage.measure_magnetization(W, G, O)
    @test magnetization isa Number
    @test isfinite(magnetization)
    @test isapprox(magnetization, 0.0, atol=1e-6)
end

nothing
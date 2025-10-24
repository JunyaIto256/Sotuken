module Lattice

using Random

# ハニカム格子上を LxL の格子で表現する（各サイトはパリティで接続パターンが変わる）
function create_honeycomb_lattice(L)::Array{Int8,2}
    lattice = zeros(Int8, L, L)
    for i in 1:L, j in 1:L
        lattice[i, j] = (i + j) % 2 == 0 ? Int8(1) : Int8(-1)
    end
    return lattice
end

# ハニカム格子の近傍取得（各サイトは最大3近傍）
# 戻り型を明示: Vector{Tuple{Int,Int}}
function get_neighbors(i::Int, j::Int, L::Int)::Vector{Tuple{Int,Int}}
    neighbors = Vector{Tuple{Int,Int}}()
    # 左右は共通
    if j > 1
        push!(neighbors, (i, j-1))
    end
    if j < L
        push!(neighbors, (i, j+1))
    end
    # 縦方向はパリティで決める（簡易的なハニカム表現）
    if (i + j) % 2 == 0
        # parity 0: 上方向を接続
        if i > 1
            push!(neighbors, (i-1, j))
        end
    else
        # parity 1: 下方向を接続
        if i < L
            push!(neighbors, (i+1, j))
        end
    end
    return neighbors
end

# 追加: 全サイトの近傍リストを事前作成（高速化のため）
function build_neighbor_list(L::Int)::Vector{Vector{Tuple{Int,Int}}}
    list = Vector{Vector{Tuple{Int,Int}}}(undef, L * L)
    for i in 1:L, j in 1:L
        idx = (i-1)*L + j
        list[idx] = get_neighbors(i, j, L)
    end
    return list
end

end # module Lattice
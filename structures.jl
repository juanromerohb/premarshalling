struct CraneTimes
    he :: Matrix{Float64}
    hf :: Matrix{Float64}
    h_ini :: Vector{Float64}
    ve :: Vector{Float64}
    vf :: Vector{Float64}
    twist :: Float64
end

struct Instance
    H :: Int
    S :: Int
    P :: Int
    C :: Int
    Cp :: Vector{Int}
    stack :: Vector{Vector{Int}}
    lct :: Float64
    cts :: CraneTimes
end

mutable struct Container
    gv :: Int
    preGv :: Int
    relAcc :: Bool
end

mutable struct Stack
    const c :: Vector{Container}
    he :: Int
    fstRelAcc :: Int
    nClean :: Int
end

mutable struct Bay
    const inst :: Instance
    const stack :: Vector{Stack}
    const relAccGv :: Vector{Int}
    lastAccGv :: Int
    acc :: Int
    ordered :: Bool
end

mutable struct Solution
    const inst :: Instance
    const bay :: Bay
    const reloc :: Vector{Tuple{Int, Int}}
    nReloc :: Int
    ct :: Float64
end
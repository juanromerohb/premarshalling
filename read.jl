function equation1(d, a)
    # Copied from C code of Parreno-Torres et al. (2020)
    return (2 * sqrt(d * a)) / a
end

function equation2(d, a, v, dtomax)
    # Copied from C code of Parreno-Torres et al. (2020)
    return (2 * (v / a)) + ((d - dtomax) / v)
end

function MovEquations(v, d)
    # Copied from C code of Parreno-Torres et al. (2020)
    return (v^2) / (2 * d)
end

function computeTimes(n_stack::Int, s_height::Int)
    # Copied from C code of Parreno-Torres et al. (2020)

    v_he = 130 / 60.0
    v_hf = 70 / 60.0
    v_ve = 60 / 60.0
    v_vf = 30 / 60.0
    cont_wide = 2.438
    cont_high = 2.591
    cont_sep = 0.3
    cont_inisep = 1.0
    vert_sep = 2.0
    d_vh = 2.50
    d_vv = 2.650
    twist_lock = 5.0

    a_he = MovEquations(v_he, d_vh)
    a_hf = MovEquations(v_hf, d_vh)
    a_ve = MovEquations(v_ve, d_vv)
    a_vf = MovEquations(v_vf, d_vv)

    cost_he = zeros(Float64, n_stack, n_stack)
    cost_hf = zeros(Float64, n_stack, n_stack)
    cost_h_ini = zeros(Float64, n_stack)
    cost_ve = zeros(Float64, s_height)
    cost_vf = zeros(Float64, s_height)

    abs_ks = 0
    dist_h_ks = 0.0
    dist_v_h = 0.0

    for s in 1:n_stack
        for k in 1:n_stack
            if s != k
                abs_ks = abs(k - s)
                dist_h_ks = abs_ks * cont_wide + abs_ks * cont_sep
                if dist_h_ks < 2 * d_vh
                    cost_he[s, k] = equation1(dist_h_ks, a_he)
                    cost_hf[s, k] = equation1(dist_h_ks, a_hf)
                else
                    cost_he[s, k] = equation2(dist_h_ks, a_he, v_he, 2 * d_vh)
                    cost_hf[s, k] = equation2(dist_h_ks, a_hf, v_hf, 2 * d_vh)
                end
            end
        end
    end

    for s in 1:n_stack
        dist_h_ks = s * cont_wide + (s - 1) * cont_sep + cont_inisep
        if dist_h_ks < 2 * d_vh
            cost_h_ini[s] = equation1(dist_h_ks, a_he)
        else
            cost_h_ini[s] = equation2(dist_h_ks, a_he, v_he, 2 * d_vh)
        end
    end

    for h in 1:s_height
        dist_v_h = vert_sep + (s_height - (h - 1)) * cont_high
        if dist_v_h < 2 * d_vv
            cost_ve[h] = equation1(dist_v_h, a_ve)
            cost_vf[h] = equation1(dist_v_h, a_vf)
        else
            cost_ve[h] = equation2(dist_v_h, a_ve, v_ve, 2 * d_vv)
            cost_vf[h] = equation2(dist_v_h, a_vf, v_vf, 2 * d_vv)
        end
    end

    return CraneTimes(cost_he, cost_hf, cost_h_ini, cost_ve, cost_vf, twist_lock)
end

function readInst(path::String, lct::Float64)
    lines = open(path, "r") do file
        readlines(file)
    end

    lines = [parse.(Int, split(line)) for line in lines]

    H = length(lines)
    S = length(lines[1])
    stack = Vector{Vector{Int}}(undef, S)

    for s in 1:S
        stack[s] = Vector{Int}(undef, 0)
        for h in 1:H
            if lines[H - h + 1][s] == 0
                break
            end

            push!(stack[s], lines[H - h + 1][s])
        end
    end

    P = maximum(maximum.([isempty(stack[s]) ? [0] : stack[s] for s in 1:S]))
    Cp = [sum([sum([stack[s][h] == gv for h in eachindex(stack[s])]) for s in 1:S]) for gv in 1:P]
    C = sum(Cp)

    Instance(H, S, P, C, Cp, stack, lct, computeTimes(S, H))
end

readCV(h::Int, s::Int, n::Int; lct::Float64 = 0.0) = readInst("instances/CV/h$(h)s$(s)n$(n).txt", lct)

readBZ(p::Int, f::Int, h::Int, s::Int, n::Int; lct::Float64 = 0.0) = readInst("instances/BZ/p$(p)f$(f)h$(h)s$(s)n$(n).txt", lct)

readZJY(h::Int, s::Int, n::Int; lct::Float64 = 0.0) = readInst("instances/ZJY/h$(h)s$(s)n$(n).txt", lct)
function selectSrcHGFH(bay::Bay)
    inst = bay.inst
    stack = bay.stack

    p, src = findmax(
        [
        (stack[s].he == stack[s].nClean) ? 0 : maximum(
            [stack[s].c[h].gv for h in stack[s].nClean+1:stack[s].he]
        )
        for s in 1:inst.S
    ]
    )

    p, src
end

function selectDstHGFH(bay::Bay, p::Int, src::Int)
    inst = bay.inst
    stack = bay.stack

    dst = 0
    r = stack[src].he - findfirst([stack[src].c[h].gv < p for h in 1:stack[src].nClean]) + 1
    maxHe = maximum([stack[s].he for s in 1:inst.S if (stack[s].he < inst.H) && (s != src)])
    dists = [abs(s - src) for s in 1:inst.S]
    validDists = [(stack[s].he == maxHe) && (s != src) for s in 1:inst.S]
    v = view(dists, validDists)
    res = first(parentindices(v))[argmin(v)]
    e = sum([inst.H - stack[k].he for k in 1:inst.S if (k != src && k != res)])

    if r - 1 <= e
        dst = src
    end

    f = 0
    for h in stack[src].he:-1:stack[src].nClean+1
        if stack[src].c[h].gv == p
            break
        end
        f += 1
    end

    #f = findfirst([stack[src].c[h].gv == p for h in stack[src].he:-1:stack[src].nClean+1]) - 1
    g = 0

    for s in 1:inst.S
        if stack[s].he == inst.H || s == src
            continue
        end

        if stack[s].nClean == 0
            g2 = stack[s].he
        elseif (stack[s].nClean == stack[s].he && stack[s].c[stack[s].he].gv >= p)
            g2 = 0
        else
            g2 = stack[s].he - findfirst([stack[s].c[h].gv < p for h in 1:stack[s].nClean]) + 1
        end

        r2 = f + g2
        if dst != 0
            if r2 < r
                e = 0
                for k in 1:inst.S
                    if k != src && k != s
                        e += inst.H - stack[k].he
                    end
                end

                if r2 <= e
                    dst = s
                    r = r2
                    g = g2
                end
            end
        else
            e = sum([inst.H - stack[k].he for k in 1:inst.S if (k != src && k != s)])

            if r2 <= e
                dst = s
                r = r2
                g = g2
            end
        end
    end

    dst, r, f, g, res
end

function auxMoveHGFH(sol::Solution, p::Int, src::Int)
    bay = sol.bay
    inst = bay.inst
    stack = bay.stack

    dists = [abs(s - src) for s in 1:inst.S]
    validDists = [(stack[s].he == inst.H) && (stack[s].c[stack[s].he].gv < p) for s in 1:inst.S]

    if !any(validDists)
        return false
    end

    v = view(dists, validDists)
    from = first(parentindices(v))[argmin(v)]

    dists = [abs(s - from) for s in 1:inst.S]
    validDists = [stack[s].he <= inst.H - 2 for s in 1:inst.S]

    if !any(validDists)
        return false
    end

    v = view(dists, validDists)
    to = first(parentindices(v))[argmin(v)]

    move(sol, from, to)
end

function relocateContHGFH(sol::Solution, p::Int, src::Int, dst::Int, r::Int, f::Int, g::Int, res::Int)
    bay = sol.bay
    inst = bay.inst
    stack = bay.stack

    if src != dst
        for _ in 1:f
            toSt = findall(
                [
                    (s != src) && (s != dst) &&
                    (stack[s].he < inst.H) &&
                    (stack[s].he == stack[s].nClean) &&
                    (stack[s].c[stack[s].he].gv >= p)
                    for s in 1:inst.S
                ]
            )

            if isempty(toSt)
                toSt = findall(
                    [
                        (s != src) && (s != dst) && (stack[s].he < inst.H)
                        for s in 1:inst.S
                    ]
                )
            elseif length(toSt) == 1
                to = toSt[1]
                if !move(sol, src, to)
                    return false
                end
                continue
            end

            maxGvs = [
                (stack[s].he == 0) ? 0 :
                ((stack[s].he == stack[s].nClean) ? inst.P + 1 :
                maximum([stack[s].c[h].gv for h in stack[s].nClean+1:stack[s].he]))
                for s in toSt
            ]
            minMaxGv = minimum(maxGvs)
            toSt = toSt[findall(maxGvs .== minMaxGv)]

            if length(toSt) == 1
                to = toSt[1]
                if !move(sol, src, to)
                    return false
                end
                continue
            end

            hes = [stack[s].he for s in toSt]
            minHe = minimum(hes)
            toSt = toSt[findall(hes .== minHe)]
            
            to = toSt[1]
            if !move(sol, src, to)
                return false
            end
        end

        for _ in 1:g
            toSt = findall(
                [
                    (s != src) && (s != dst) &&
                    (stack[s].he < inst.H) &&
                    (stack[s].he == stack[s].nClean) &&
                    (stack[s].c[stack[s].he].gv >= p)
                    for s in 1:inst.S
                ]
            )

            if isempty(toSt)
                toSt = findall(
                    [
                        (s != src) && (s != dst) && (stack[s].he < inst.H)
                        for s in 1:inst.S
                    ]
                )
            elseif length(toSt) == 1
                to = toSt[1]
                if !move(sol, dst, to)
                    return false
                end
                continue
            end

            maxGvs = [
                (stack[s].he == 0) ? 0 :
                ((stack[s].he == stack[s].nClean) ? inst.P + 1 :
                maximum([stack[s].c[h].gv for h in stack[s].nClean+1:stack[s].he]))
                for s in toSt
            ]
            minMaxGv = minimum(maxGvs)
            toSt = toSt[findall(maxGvs .== minMaxGv)]

            if length(toSt) == 1
                to = toSt[1]
                if !move(sol, dst, to)
                    return false
                end
                continue
            end

            hes = [stack[s].he for s in toSt]
            minHe = minimum(hes)
            toSt = toSt[findall(hes .== minHe)]
            
            to = toSt[1]
            if !move(sol, dst, to)
                return false
            end
        end

        if !move(sol, src, dst)
            return false
        end
    else
        wasMoved = false

        for _ in 1:r
            if !wasMoved && (stack[src].c[stack[src].he].gv == p)
                if !move(sol, src, res)
                    return false
                end

                wasMoved = true
                continue
            end

            toSt = findall(
                [
                    (s != src) && (s != res) &&
                    (stack[s].he < inst.H) &&
                    (stack[s].he == stack[s].nClean) &&
                    (stack[s].c[stack[s].he].gv >= p)
                    for s in 1:inst.S
                ]
            )

            if isempty(toSt)
                toSt = findall(
                    [
                        (s != src) && (s != res) && (stack[s].he < inst.H)
                        for s in 1:inst.S
                    ]
                )
            elseif length(toSt) == 1
                to = toSt[1]
                if !move(sol, src, to)
                    return false
                end
                continue
            end

            maxGvs = [
                (stack[s].he == 0) ? 0 :
                ((stack[s].he == stack[s].nClean) ? inst.P + 1 :
                maximum([stack[s].c[h].gv for h in stack[s].nClean+1:stack[s].he]))
                for s in toSt
            ]
            minMaxGv = minimum(maxGvs)
            toSt = toSt[findall(maxGvs .== minMaxGv)]

            if length(toSt) == 1
                to = toSt[1]
                if !move(sol, src, to)
                    return false
                end
                continue
            end

            hes = [stack[s].he for s in toSt]
            minHe = minimum(hes)
            toSt = toSt[findall(hes .== minHe)]
            
            to = toSt[1]
            if !move(sol, src, to)
                return false
            end
        end

        if stack[res].nClean != stack[res].he
            if !move(sol, res, src)
                return false
            end
        end
    end

    true
end

function HGFH(sol::Solution)
    bay = sol.bay

    while !bay.ordered
        # SelectSRC
        p, src = selectSrcHGFH(bay)

        @label beforeSelectDST
        # SelectDST
        dst, r, f, g, res = selectDstHGFH(bay, p, src)

        if dst == 0
            if !auxMoveHGFH(sol, p, src)
                return sol
            else
                @goto beforeSelectDST
            end
        end

        if !relocateContHGFH(sol, p, src, dst, r, f, g, res)
            return sol
        end
    end

    return sol
end

function selectSrcHGFH_R(bay::Bay)
    inst = bay.inst
    stack = bay.stack

    maxPs = [
        (stack[s].he == stack[s].nClean) ? 0 : maximum(
            [stack[s].c[h].gv for h in stack[s].nClean+1:stack[s].he]
        )
        for s in 1:inst.S
    ]

    p = maximum(maxPs)

    src = rand(findall(maxPs .== p))

    p, src
end

function selectDstHGFH_R(sol::Solution, bay::Bay, p::Int, src::Int)
    inst = bay.inst
    stack = bay.stack

    dstSts = Vector{Tuple{Int, Int, Int, Int, Int}}(undef, 0)

    dst = 0
    r = stack[src].he - findfirst([stack[src].c[h].gv < p for h in 1:stack[src].nClean]) + 1
    maxHe = maximum([stack[s].he for s in 1:inst.S if (stack[s].he < inst.H) && (s != src)])
    dists = [abs(s - src) for s in 1:inst.S]
    validDists = [(stack[s].he == maxHe) && (s != src) for s in 1:inst.S]
    v = view(dists, validDists)
    res = first(parentindices(v))[argmin(v)]
    e = sum([inst.H - stack[k].he for k in 1:inst.S if (k != src && k != res)])

    if r - 1 <= e
        dst = src
        push!(dstSts, (dst, r, 0, 0, res))
    end

    f = findfirst([stack[src].c[h].gv == p for h in stack[src].he:-1:stack[src].nClean+1]) - 1
    g = 0

    for s in 1:inst.S
        if s == src || stack[s].he == inst.H
            continue
        end

        if stack[s].nClean == 0
            g2 = stack[s].he
        elseif (stack[s].nClean == stack[s].he && stack[s].c[stack[s].he].gv >= p)
            g2 = 0
        else
            g2 = stack[s].he - findfirst([stack[s].c[h].gv < p for h in 1:stack[s].nClean]) + 1
        end

        r2 = f + g2
        if dst != 0
            if r2 < r
                e = sum([inst.H - stack[k].he for k in 1:inst.S if (k != src && k != s)])

                if r2 <= e
                    dst = s
                    r = r2
                    g = g2
                    push!(dstSts, (dst, r, f, g, res))
                end
            end
        else
            e = sum([inst.H - stack[k].he for k in 1:inst.S if (k != src && k != s)])

            if r2 <= e
                dst = s
                r = r2
                g = g2
                push!(dstSts, (dst, r, f, g, res))
            end
        end
    end

    if isempty(dstSts)
        return 0, 0, 0, 0, 0
    end

    minR = minimum([dstSts[i][2] for i in eachindex(dstSts)])
    dstSts = dstSts[findall([dstSts[i][2] == minR for i in eachindex(dstSts)])]

    rand(dstSts)
end

function auxMoveHGFH_R(sol::Solution, p::Int, src::Int)
    bay = sol.bay
    inst = bay.inst
    stack = bay.stack

    validSts = [
        s
        for s in 1:inst.S
        if (stack[s].he == inst.H) && (stack[s].c[stack[s].he].gv < p)
    ]

    if isempty(validSts)
        return false
    end

    from = rand(validSts)

    validSts = [s for s in 1:inst.S if stack[s].he <= inst.H - 2]

    if isempty(validSts)
        return false
    end

    to = rand(validSts)

    move(sol, from, to)
end

function relocateContHGFH_R(sol::Solution, p::Int, src::Int, dst::Int, r::Int, f::Int, g::Int, res::Int)
    bay = sol.bay
    inst = bay.inst
    stack = bay.stack

    if src != dst

        while f > 0 || g > 0
            if f == 0
                from = dst
                g -= 1
            elseif g == 0
                from = src
                f -= 1
            else
                from = rand([src, dst])
                if from == src
                    f -= 1
                else
                    g -= 1
                end
            end

            toSt = findall(
                [
                    (s != src) && (s != dst) && (stack[s].he < inst.H)
                    for s in 1:inst.S
                ]
            )

            if isempty(toSt)
                return false
            end

            to = rand(toSt)

            if !move(sol, from, to)
                return false
            end
        end

        if !move(sol, src, dst)
            return false
        end
    else
        wasMoved = false

        for _ in 1:r
            if !wasMoved && (stack[src].c[stack[src].he].gv == p)
                if !move(sol, src, res)
                    return false
                end

                wasMoved = true
                continue
            end

            toSt = findall(
                [
                    (s != src) && (s != res) && (stack[s].he < inst.H)
                    for s in 1:inst.S
                ]
            )

            if isempty(toSt)
                return false
            end

            to = rand(toSt)

            if !move(sol, src, to)
                return false
            end
        end

        if stack[res].nClean != stack[res].he
            if !move(sol, res, src)
                return false
            end
        end
    end

    true
end


function HGFH_R(sol::Solution)
    bay = sol.bay

    while !bay.ordered
        # SelectSRC
        p, src = selectSrcHGFH_R(bay)

        @label beforeSelectDST
        # SelectDST
        dst, r, f, g, res = selectDstHGFH_R(sol, bay, p, src)

        if dst == 0
            if !auxMoveHGFH_R(sol, p, src)
                return sol
            else
                if bay.ordered
                    return sol
                end
                @goto beforeSelectDST
            end
        end

        if !relocateContHGFH_R(sol, p, src, dst, r, f, g, res)
            return sol
        end
    end

    return sol
end
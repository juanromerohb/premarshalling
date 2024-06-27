function selectStackAGH(bay::Bay, p::Int, blStacks::Vector{Int})
    stack = bay.stack
    nBlocking = Vector{Int}(undef, length(blStacks))

    for i in eachindex(blStacks)
        blockingH = 0

        for h in stack[blStacks[i]].fstRelAcc:-1:1
            if stack[blStacks[i]].c[h].gv > p
                blockingH = h
            elseif stack[blStacks[i]].c[h].gv == p
                break
            end
        end

        nBlocking[i] = stack[blStacks[i]].he - blockingH + 1
    end

    minBlocking = minimum(nBlocking)
    blStacks = [blStacks[i] for i in eachindex(blStacks) if nBlocking[i] == minBlocking]
    nBlocking = [nBlocking[i] for i in eachindex(nBlocking) if nBlocking[i] == minBlocking]

    if length(blStacks) == 1
        return blStacks[1], nBlocking[1]
    end

    minGvBlocking = [
        minimum(
            [
                stack[blStacks[i]].c[h].gv
                for h in stack[blStacks[i]].he:-1:stack[blStacks[i]].he-nBlocking[i]+1
            ]
        )
        for i in eachindex(blStacks)
    ]

    maxMinGvBlocking = maximum(minGvBlocking)
    blStacks = [blStacks[i] for i in eachindex(blStacks) if minGvBlocking[i] == maxMinGvBlocking]
    nBlocking = [nBlocking[i] for i in eachindex(nBlocking) if minGvBlocking[i] == maxMinGvBlocking]

    if length(blStacks) == 1
        return blStacks[1], nBlocking[1]
    end

    maxH = maximum([stack[s].he for s in blStacks])
    nBlocking = [nBlocking[i] for i in eachindex(nBlocking) if stack[blStacks[i]].he == maxH]
    blStacks = [s for s in blStacks if stack[s].he == maxH]

    if length(blStacks) == 1
        return blStacks[1], nBlocking[1]
    end

    center = (bay.inst.S + 1) / 2

    minDistC = minimum([abs(s - center) for s in blStacks])
    nBlocking = [nBlocking[i] for i in eachindex(nBlocking) if abs(blStacks[i] - center) == minDistC]
    blStacks = [s for s in blStacks if abs(s - center) == minDistC]

    blStacks[1], nBlocking[1]
end

function tryRelocate(sol::Solution, p::Int, st::Int)
    bay = sol.bay
    inst = bay.inst
    stack = bay.stack
    cTimeSK = 0.0

    Phi = Vector{Int}(undef, 0)

    for s in 1:inst.S
        if (s == st) || (stack[s].he == inst.H) || (sol.ct + calcCtime(sol, st, s) > inst.lct)
            continue
        end
        
        willBlock = false

        if stack[s].he > 0
            for h in stack[s].fstRelAcc:stack[s].he
                if (stack[s].c[h].gv <= p) && (stack[s].c[h].gv < stack[st].c[stack[st].he].gv)
                    willBlock = true
                    break
                end
            end
        end

        if !willBlock
            push!(Phi, s)
        end
    end

    if !isempty(Phi)
        if length(Phi) == 1
            to = Phi[1]

            if !move(sol, st, to)
                return false
            else
                return true
            end
        end

        toSt = [s for s in Phi if stack[s].nClean == stack[s].he]

        if isempty(toSt)
            toSt = Phi
        end

        if length(toSt) == 1
            to = toSt[1]

            if !move(sol, st, to)
                return false
            else
                return true
            end
        end

        minGv = [
            (stack[s].he == 0) ? inst.P + 1 : minimum([stack[s].c[h].gv for h in 1:stack[s].he])
            for s in toSt
        ]
        maxMinGv = maximum(minGv)
        toSt = [toSt[i] for i in eachindex(toSt) if minGv[i] == maxMinGv]

        if length(toSt) == 1
            to = toSt[1]

            if !move(sol, st, to)
                return false
            else
                return true
            end
        end

        cTimes = [calcCtime(sol, st, s) for s in toSt]
        minCTime = minimum(cTimes)
        toSt = [toSt[i] for i in eachindex(toSt) if cTimes[i] == minCTime]

        if length(toSt) == 1
            to = toSt[1]

            if !move(sol, st, to)
                return false
            else
                return true
            end
        end

        minHe = minimum([stack[s].he for s in toSt])
        toSt = [s for s in toSt if stack[s].he == minHe]

        to = toSt[1]

        if !move(sol, st, to)
            return false
        else
            return true
        end
    else
        Psi = Vector{Int}(undef, 0)

        for s in 1:inst.S
            if s == st || stack[s].he == 0
                continue
            end

            stack[s].he -= 1
            cTimeSK = calcCtime(sol, st, s)
            stack[s].he += 1
            if sol.ct + cTimeSK > inst.lct
                continue
            end

            willBlock = false

            for h in stack[s].fstRelAcc:stack[s].he-1
                if (stack[s].c[h].gv <= p) && (stack[s].c[h].gv < stack[st].c[stack[st].he].gv)
                    willBlock = true
                    break
                end
            end

            if willBlock
                continue
            end

            canBeMoved = false
            for q in 1:inst.S
                if (q == st) || (q == s) ||
                    (stack[q].he == inst.H) ||
                    (sol.ct + cTimeSK + calcCtime(sol, s, q) > inst.lct)
                    continue
                end

                willBlock = false

                for h in stack[q].fstRelAcc:stack[q].he
                    if (stack[q].c[h].gv <= p) && (stack[q].c[h].gv < stack[s].c[stack[s].he].gv)
                        willBlock = true
                        break
                    end
                end

                if !willBlock
                    canBeMoved = true
                    break
                end
            end

            if canBeMoved
                push!(Psi, s)
            end
        end

        if isempty(Psi)
            return false
        end

        toSt = Psi

        if length(toSt) > 1
            toSt = [s for s in toSt if stack[s].nClean >= stack[s].he - 1]

            if isempty(toSt)
                toSt = Psi
            end

            if length(toSt) > 1
                minGv = [
                    (stack[s].he == 1) ? inst.P + 1 : minimum([stack[s].c[h].gv for h in 1:stack[s].he-1])
                    for s in toSt
                ]
                maxMinGv = maximum(minGv)
                toSt = [toSt[i] for i in eachindex(toSt) if minGv[i] == maxMinGv]

                if length(toSt) > 1
                    topGv = [stack[s].c[stack[s].he].gv for s in toSt]
                    maxTopGv = maximum(topGv)
                    toSt = [toSt[i] for i in eachindex(toSt) if topGv[i] == maxTopGv]

                    if length(toSt) > 1
                        cTimes = Vector{Float64}(undef, length(toSt))

                        for i in eachindex(toSt)
                            stack[toSt[i]].he -= 1
                            cTimes[i] = calcCtime(sol, st, toSt[i])
                            stack[toSt[i]].he += 1
                        end

                        minCTime = minimum(cTimes)
                        toSt = [toSt[i] for i in eachindex(toSt) if cTimes[i] == minCTime]

                        if length(toSt) > 1
                            center = (inst.S + 1) / 2
                            minDistC = minimum([abs(s - center) for s in toSt])
                            toSt = [s for s in toSt if abs(s - center) == minDistC]

                            toSt = [toSt[1]]
                        end
                    end
                end
            end
        end

        to = toSt[1]

        Psi2 = Vector{Int}(undef, 0)

        for q in 1:inst.S
            if (q == st) || (q == to) ||
                (stack[q].he == inst.H) ||
                (sol.ct + cTimeSK + calcCtime(sol, to, q) > inst.lct)
                continue
            end

            willBlock = false

            for h in stack[q].fstRelAcc:stack[q].he
                if (stack[q].c[h].gv <= p) && (stack[q].c[h].gv < stack[to].c[stack[to].he].gv)
                    willBlock = true
                    break
                end
            end

            if !willBlock
                push!(Psi2, q)
            end
        end

        if length(Psi2) == 1
            to2 = Psi2[1]

            if !move(sol, to, to2)
                return false
            else
                if !move(sol, st, to)
                    return false
                else
                    return true
                end
            end
        end

        toSt = [s for s in Psi2 if stack[s].nClean == stack[s].he]

        if isempty(toSt)
            toSt = Psi2

            if isempty(toSt)
                return false
            end
        end

        if length(toSt) == 1
            to2 = toSt[1]

            if !move(sol, to, to2)
                return false
            else
                if !move(sol, st, to)
                    return false
                else
                    return true
                end
            end
        end

        minGv = [
            (stack[s].he == 0) ? inst.P + 1 : minimum([stack[s].c[h].gv for h in 1:stack[s].he])
            for s in toSt
        ]
        maxMinGv = maximum(minGv)
        toSt = [toSt[i] for i in eachindex(toSt) if minGv[i] == maxMinGv]

        if length(toSt) == 1
            to2 = toSt[1]

            if !move(sol, to, to2)
                return false
            else
                if !move(sol, st, to)
                    return false
                else
                    return true
                end
            end
        end

        cTimes = [calcCtime(sol, st, s) for s in toSt]
        minCTime = minimum(cTimes)
        toSt = [toSt[i] for i in eachindex(toSt) if cTimes[i] == minCTime]

        if length(toSt) == 1
            to2 = toSt[1]

            if !move(sol, to, to2)
                return false
            else
                if !move(sol, st, to)
                    return false
                else
                    return true
                end
            end
        end

        minHe = minimum([stack[s].he for s in toSt])
        toSt = [s for s in toSt if stack[s].he == minHe]

        to2 = toSt[1]

        if !move(sol, to, to2)
            return false
        else
            if !move(sol, st, to)
                return false
            else
                return true
            end
        end
    end
end

function tryReorder(sol::Solution, p::Int, st::Int)
    bay = sol.bay
    inst = bay.inst
    stack = bay.stack

    maxHe = maximum([stack[s].he for s in 1:inst.S if s != st && stack[s].he < inst.H])
    resSt = [s for s in 1:inst.S if stack[s].he == maxHe && s != st && stack[s].he < inst.H]
    dists = [abs(s - st) for s in resSt]
    res = resSt[argmin(dists)]

    movedTo = zeros(Int, inst.S)

    for _ in 1:inst.H
        toSt = [s for s in 1:inst.S if stack[s].he < inst.H && s != st && s != res]

        if isempty(toSt)
            return false
        end

        cTimes = [calcCtime(sol, st, s) for s in toSt]
        to = toSt[argmin(cTimes)]

        if !move(sol, st, to)
            return false
        end

        movedTo[to] += 1

        if stack[st].c[stack[st].he].gv == p
            break
        end
    end

    if !move(sol, st, res)
        return false    
    end

    fromSt = [s for s in 1:inst.S if movedTo[s] > 0]
    dists = [abs(s - st) for s in fromSt]
    from = fromSt[argmin(dists)]

    if !move(sol, from, st)
        return false
    end

    movedTo[from] -= 1

    for s in 1:inst.S
        for _ in 1:movedTo[s]
            if !move(sol, s, st)
                return false
            end
        end
    end

    isBlocking = false

    for h in 1:stack[res].he-1
        if stack[res].c[h].gv < p
            isBlocking = true
            break
        end
    end

    if isBlocking
        if !move(sol, res, st)
            return false
        end
    end

    true
end


function AGH(sol::Solution)
    bay = sol.bay
    inst = bay.inst
    stack = bay.stack

    while !bay.ordered
        p = bay.lastAccGv
        blStacks = [
            s
            for s in 1:inst.S
            if stack[s].he > 0 &&
                any(
                    [
                        stack[s].c[h].gv == p && !stack[s].c[h].relAcc
                        for h in 1:stack[s].he
                    ]
                )
        ]

        while !isempty(blStacks)
            st, t = selectStackAGH(bay, p, blStacks)

            success = true
            for _ in 1:t
                if !tryRelocate(sol, p, st)
                    deleteat!(blStacks, findfirst(isequal(st), blStacks))
                    success = false
                    break
                end
            end

            if success
                if !any(
                    [
                        stack[st].c[h].gv == p && !stack[st].c[h].relAcc
                        for h in 1:stack[st].he
                    ]
                )
                    deleteat!(blStacks, findfirst(isequal(st), blStacks))
                end
            end
        end

        blStacks = [
            s
            for s in 1:inst.S
            if stack[s].he > 0 &&
                any(
                    [
                        stack[s].c[h].gv == p && !stack[s].c[h].relAcc
                        for h in 1:stack[s].he
                    ]
                )
        ]

        if !isempty(blStacks)
            st = selectStackAGH(bay, p, blStacks)[1]
            solSave = copySol(sol)

            if !tryReorder(sol, p, st)
                sol = solSave
                return sol
            end
        end

        blStacks = [
            s
            for s in 1:inst.S
            if stack[s].he > 0 &&
                any(
                    [
                        stack[s].c[h].gv == p && !stack[s].c[h].relAcc
                        for h in 1:stack[s].he
                    ]
                )
        ]
    end

    sol
end

function reduceBay(bay::Bay, p::Int)
    stack = bay.stack
    inst = bay.inst

    for s in 1:inst.S
        for h in 1:stack[s].he
            if stack[s].c[h].gv > p
                inst.Cp[stack[s].c[h].gv] -= 1
                stack[s].c[h].gv = p + 1
                inst.Cp[p + 1] += 1
            end

            stack[s].c[h].relAcc = false
        end
    end

    for s in 1:inst.S
        if stack[s].he == 1
            stack[s].nClean = 1
        end

        if stack[s].he > 0
            stack[s].fstRelAcc = 1
        end

        for h in stack[s].he-1:-1:1
            if stack[s].c[h].gv < stack[s].c[h+1].gv
                stack[s].fstRelAcc = h + 1
                break
            end
        end

        stack[s].nClean = stack[s].he

        for h in 2:stack[s].he
            if stack[s].c[h].gv > stack[s].c[h-1].gv
                stack[s].nClean = h - 1
                break
            end
        end

        for h in 1:stack[s].he
            if all([stack[s].c[h].gv >= stack[s].c[e].gv for e in h+1:stack[s].he])
                stack[s].c[h].relAcc = true
            end
        end
    end

    for gv in 1:inst.P
        bay.relAccGv[gv] = sum(
            [
                sum([stack[s].c[h].relAcc for h in 1:stack[s].he if stack[s].c[h].gv == gv])
                for s in 1:inst.S
            ]
        )
    end

    bay.acc = 0

    for gv in 1:inst.P
        bay.acc += bay.relAccGv[gv]

        if bay.relAccGv[gv] < inst.Cp[gv]
            bay.lastAccGv = gv
            break
        end
    end

    bay.ordered = (bay.acc == inst.C)

    bay
end

function recoverBay(bay::Bay)
    stack = bay.stack
    inst = bay.inst

    for s in 1:inst.S
        for h in 1:stack[s].he
            if stack[s].c[h].gv != stack[s].c[h].preGv
                inst.Cp[stack[s].c[h].gv] -= 1
                stack[s].c[h].gv = copy(stack[s].c[h].preGv)
                inst.Cp[stack[s].c[h].gv] += 1
            end

            stack[s].c[h].relAcc = false
        end
    end

    for s in 1:inst.S
        if stack[s].he == 1
            stack[s].nClean = 1
        end

        if stack[s].he > 0
            stack[s].fstRelAcc = 1
        end

        for h in stack[s].he-1:-1:1
            if stack[s].c[h].gv < stack[s].c[h+1].gv
                stack[s].fstRelAcc = h + 1
                break
            end
        end

        stack[s].nClean = stack[s].he

        for h in 2:stack[s].he
            if stack[s].c[h].gv > stack[s].c[h-1].gv
                stack[s].nClean = h - 1
                break
            end
        end

        for h in 1:stack[s].he
            if all([stack[s].c[h].gv >= stack[s].c[e].gv for e in h+1:stack[s].he])
                stack[s].c[h].relAcc = true
            end
        end
    end

    for gv in 1:inst.P
        bay.relAccGv[gv] = sum(
            [
                sum([stack[s].c[h].relAcc for h in 1:stack[s].he if stack[s].c[h].gv == gv])
                for s in 1:inst.S
            ]
        )
    end

    bay.acc = 0

    for gv in 1:inst.P
        bay.acc += bay.relAccGv[gv]

        if bay.relAccGv[gv] < inst.Cp[gv]
            bay.lastAccGv = gv
            break
        end
    end

    bay.ordered = (bay.acc == inst.C)

    bay
end

function SM_HGFH_R(sol::Solution, iter::Int)
    inSol = copySol(sol)
    bestSol = copySol(inSol)

    normalCp = deepcopy(inSol.bay.inst.Cp)

    for p in sol.bay.inst.P:-1:1
        newSol = copySol(inSol)
        reduceBay(newSol.bay, p)
        reducedCp = deepcopy(newSol.bay.inst.Cp)

        tempSol = copySol(newSol)
        tempSol = HGFH(tempSol)

        if tempSol.bay.ordered
            recoverBay(tempSol.bay)
            tempSol = AGH(tempSol)

            if tempSol.bay.acc > bestSol.bay.acc ||
                (tempSol.bay.acc == bestSol.bay.acc && tempSol.ct < bestSol.ct)
                bestSol = copySol(tempSol)
            end

            [tempSol.inst.Cp[gv] = reducedCp[gv] for gv in 1:tempSol.inst.P]
        end

        for _ in 1:iter
            tempSol = copySol(newSol)
            tempSol = HGFH_R(tempSol)

            if tempSol.bay.ordered
                recoverBay(tempSol.bay)
                tempSol = AGH(tempSol)

                if tempSol.bay.acc > bestSol.bay.acc ||
                    (tempSol.bay.acc == bestSol.bay.acc && tempSol.ct < bestSol.ct)
                    bestSol = copySol(tempSol)
                end

                [tempSol.inst.Cp[gv] = reducedCp[gv] for gv in 1:tempSol.inst.P]
            end
        end

        [newSol.inst.Cp[gv] = normalCp[gv] for gv in 1:newSol.inst.P]
    end
    
    bestSol
end
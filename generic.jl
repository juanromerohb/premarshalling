function initialSol(inst::Instance)
    stack = [
        Stack(
            [
                (h > length(inst.stack[s])) ?
                Container(0, 0, false) :
                Container(inst.stack[s][h], inst.stack[s][h], false)
                for h in 1:inst.H
            ],
            0,
            0,
            0,
        )
        for s in 1:inst.S
    ]

    for s in 1:inst.S
        stack[s].he = length(inst.stack[s])
    end

    bay = Bay(inst, stack, zeros(Int, inst.P), inst.P, 0, false)

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

    for gv in 1:inst.P
        bay.acc += bay.relAccGv[gv]

        if bay.relAccGv[gv] < inst.Cp[gv]
            bay.lastAccGv = gv
            break
        end
    end

    bay.ordered = (bay.acc == inst.C)

    Solution(inst, bay, Vector{Tuple{Int, Int}}(undef, 0), 0, 0.0)
end

function calcCtime(sol::Solution, src::Int, dst::Int)
    bay = sol.bay
    stack = bay.stack
    inst = sol.inst
    cts = inst.cts
    reloc = sol.reloc
    newct = 0.0

    if sol.nReloc == 0
        newct += cts.h_ini[src]
    else
        newct += cts.he[reloc[end][2], src]
    end

    newct += cts.hf[src, dst] + cts.ve[stack[src].he] + cts.vf[stack[src].he] + (inst.H - (stack[src].he - 1))*cts.twist

    newct += cts.ve[stack[dst].he + 1] + cts.vf[stack[dst].he + 1]

    return newct
end

function move(sol::Solution, src::Int, dst::Int)
    bay = sol.bay
    stack = bay.stack
    reloc = sol.reloc
    inst = sol.inst

    if src == dst
        error("Source and destination stacks are the same")
    end

    if stack[src].he == 0
        error("Source stack is empty")
    end

    if stack[dst].he == bay.inst.H
        error("Destination stack is full")
    end

    newct = calcCtime(sol, src, dst)

    if sol.ct + newct > sol.inst.lct && sol.inst.lct > 0.0
        return false
    end

    stack[dst].c[stack[dst].he + 1], stack[src].c[stack[src].he] = stack[src].c[stack[src].he], stack[dst].c[stack[dst].he + 1]
    stack[src].he -= 1
    stack[dst].he += 1
    sol.nReloc += 1
    
    if sol.nReloc > 100000
        error("More than 1000 relocations")
    end

    push!(reloc, (src, dst))
    sol.ct += newct

    cont = stack[dst].c[stack[dst].he]

    # Changes in src

    if stack[src].fstRelAcc == stack[src].he + 1
        if stack[src].he == 0
            stack[src].fstRelAcc = 0
        else
            stack[src].fstRelAcc = 1

            for h in stack[src].he-1:-1:1
                if stack[src].c[h].gv < stack[src].c[h+1].gv
                    stack[src].fstRelAcc = h + 1
                    break
                end
            end

            for h in 1:stack[src].he
                if !stack[src].c[h].relAcc &&
                    cont.gv > stack[src].c[h].gv &&
                    all([stack[src].c[h].gv >= stack[src].c[e].gv for e in h+1:stack[src].he])
                    
                    stack[src].c[h].relAcc = true
                    bay.relAccGv[stack[src].c[h].gv] += 1

                    if stack[src].c[h].gv == bay.lastAccGv
                        bay.acc += 1

                        if bay.relAccGv[stack[src].c[h].gv] == bay.inst.Cp[stack[src].c[h].gv]
                            bay.lastAccGv = bay.inst.P + 1
                            for gv in stack[src].c[h].gv+1:inst.P
                                bay.acc += bay.relAccGv[gv]
                                if bay.relAccGv[gv] < bay.inst.Cp[gv]
                                    bay.lastAccGv = gv
                                    break
                                end
                            end
                        end

                        if bay.acc == bay.inst.C
                            bay.ordered = true
                        end
                    end
                end
            end
        end
    end



    if stack[src].nClean == stack[src].he + 1
        stack[src].nClean -= 1
    end

    # Changes in dst

    if stack[dst].he > 1 && cont.gv > stack[dst].c[stack[dst].he - 1].gv
        stack[dst].fstRelAcc = stack[dst].he
    elseif stack[dst].he == 1
        stack[dst].fstRelAcc = 1
    end

    if stack[dst].he == 1 || (stack[dst].nClean == stack[dst].he - 1 && cont.gv <= stack[dst].c[stack[dst].he - 1].gv)
        stack[dst].nClean += 1
    end

    for h in 1:stack[dst].he - 1
        if stack[dst].c[h].relAcc && cont.gv > stack[dst].c[h].gv
            stack[dst].c[h].relAcc = false
            bay.relAccGv[stack[dst].c[h].gv] -= 1

            if stack[dst].c[h].gv <= bay.lastAccGv
                bay.acc = 0

                for gv in 1:stack[dst].c[h].gv
                    bay.acc += bay.relAccGv[gv]
                    if bay.relAccGv[gv] < bay.inst.Cp[gv]
                        bay.lastAccGv = gv
                        break
                    end
                end

                if bay.acc < bay.inst.C
                    bay.ordered = false
                end
            end
        end
    end

    
    #print("Move from $src to $dst\n")
    #printBay(bay)

    true
end

function copyContainer(cont::Container)
    Container(cont.gv, cont.preGv, cont.relAcc)
end

function copyStack(st::Stack, H::Int)
    Stack([copyContainer(st.c[h]) for h in 1:H], st.he, st.fstRelAcc, st.nClean)
end

function copyBay(bay::Bay)
    Bay(bay.inst, [copyStack(bay.stack[s], bay.inst.H) for s in 1:bay.inst.S], deepcopy(bay.relAccGv), bay.lastAccGv, bay.acc, bay.ordered)
end

function copySol(sol::Solution)
    Solution(sol.inst, copyBay(sol.bay), deepcopy(sol.reloc), sol.nReloc, sol.ct)
end

function ctMin(bay::Bay, h::Int)
    inst = bay.inst
    cts = inst.cts
    
    cts.he[1, 2] + cts.ve[h] + cts.vf[h] + (bay.inst.H - h + 1)*cts.twist + cts.hf[1, 2] + cts.vf[inst.H] + cts.ve[inst.H]
end
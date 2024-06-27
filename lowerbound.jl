function lbCPMPLCT(bay::Bay, n::Int)
    inst = bay.inst
    stack = bay.stack
    lb = 0.0
    m = 0
    p0 = 0

    sumCp = 0
    for p in 1:inst.P
        sumCp += inst.Cp[p]

        if sumCp > n
            sumCp -= inst.Cp[p]
            p0 = p - 1
            m = n - sumCp
            break
        end
    end

    blHeights = [stack[s].he + 1 for s in 1:inst.S]

    for s in 1:inst.S
        h = 1

        while h < blHeights[s]
            if stack[s].c[h].gv <= p0 && !stack[s].c[h].relAcc
                for e in h+1 : blHeights[s]-1
                    if stack[s].c[e].gv > stack[s].c[h].gv
                        blHeights[s] = e
                        break
                    end
                end
            end

            h += 1
        end
    end

    for s in 1:inst.S
        for h in blHeights[s]:stack[s].he
            lb += ctMin(bay, h)
        end
    end

    RE = Vector{Int}(undef, 0)

    hMax = 1

    for s in 1:inst.S
        for h in stack[s].he:-1:1
            if stack[s].c[h].gv == p0 + 1
                if stack[s].c[h].relAcc
                    push!(RE, 0)
                    continue
                end

                val = 0

                for e in h+1 : blHeights[s]-1
                    if stack[s].c[e].gv > p0 + 1
                        val = blHeights[s] - e
                        blHeights[s] = e
                        break
                    end
                end

                push!(RE, val)

                if val > 0 && stack[s].he > hMax
                    hMax = stack[s].he
                end
            end
        end
    end

    sort!(RE)
    ctMinHMax = ctMin(bay, hMax)

    for i in 1:m
        lb += ctMinHMax*RE[i]
    end

    lb
end
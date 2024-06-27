function printBay(bay::Bay; f=nothing)
    inst = bay.inst
    stack = bay.stack

    if isnothing(f)
        for h in inst.H:-1:1
            if h < 10
                print("$(h)  ")
            else
                print("$(h) ")
            end
            for s in 1:inst.S
                if bay.stack[s].he < h
                    print("     ")
                else
                    if stack[s].c[h].gv < 10
                        print("[$(stack[s].c[h].gv) ] ")
                    else
                        print("[$(stack[s].c[h].gv)] ")
                    end
                end
            end
    
            print("\n")
        end
    
        print("    ")

        for s in 1:inst.S
            print("$(s)    ")
        end

        print("\n\n")
    else
        open("out/$(f)", "w") do file
            for h in inst.H:-1:1
                if h < 10
                    write(file, "$(h)  ")
                else
                    write(file, "$(h) ")
                end

                for s in 1:inst.S
                    if bay.stack[s].he < h
                        write(file, "     ")
                    else
                        if stack[s].c[h].gv < 10
                            write(file, "[$(stack[s].c[h].gv) ] ")
                        else
                            write(file, "[$(stack[s].c[h].gv)] ")
                        end
                    end
                end
        
                write(file, "\n")
            end
        
            write(file, "    ")

            for s in 1:inst.S
                write(file, "$(s)    ")
            end

            write(file, "\n\n")
        end
    end
end
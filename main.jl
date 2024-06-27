include("structures.jl");
include("read.jl");
include("print.jl");
include("generic.jl");
include("cpmp_heur.jl");
include("cpmplct_heur.jl");
include("lowerbound.jl");

inst = readCV(5, 6, 10; lct=900.0)
sol = initialSol(inst)
sol = SM_HGFH_R(sol, 100)

printBay(sol.bay)
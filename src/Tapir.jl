module Tapir

using JuMP
using Ipopt
using FastGaussQuadrature
using QuadGK

export Tapir_model, feasibility_optimize!, control_optimize!,
    add_pathcost!, add_boundarycost!,
    add_path_eq_constr!,add_path_ineq_constr!,
    add_boundary_eq_constr!,add_boundary_ineq_constr!,
    update_settings!, plot_solution,
    Tapir_settings,  Tapir_mesh,
    plot

include("Tapir_Settings.jl");
include("Tapir_Mesh.jl");
include("Tapir_Solution.jl");
include("Tapir_Model.jl");
include("Tapir_Solve.jl");
include("Tapir_plotting.jl");



end # module Tapir

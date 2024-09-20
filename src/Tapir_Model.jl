mutable struct Tapir_model
    nx::Int
    nu::Int
    ne::Int
    
    path_cost::Function
    terminal_cost::Function

    path_eq_const::Function
    path_ineq_const::Function

    boundary_eq_const::Function
    boundary_ineq_const::Function

    settings::Tapir_settings;
    mesh::Tapir_mesh;
    mesh_hist::Vector{Tapir_mesh}
    solution::Tapir_solution;
    solution_hist::Vector{Tapir_solution}

    Xguess::Vector{Float64}
    Uguess::Vector{Float64}
    Tguess::Vector{Float64}
end
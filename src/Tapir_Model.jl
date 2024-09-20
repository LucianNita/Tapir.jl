mutable struct Tapir_model
    path_cost::Function
    terminal_cost::Function

    path_eq_const::Function
    path_ineq_const::Function

    boundary_eq_const::Function
    boundary_ineq_const::Function

    problem_settings::Tapir_settings;
    problem_mesh::Tapir_mesh;
    mesh_hist::Vector{Tapir_mesh}
    solution::Tapir_solution;
    solution_hist::Vector{Tapir_solution}
end
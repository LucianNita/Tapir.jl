mutable struct Tapir_mesh
    N::Integer
    Px::Integer
    Pu::Integer
    Q::Integer

    model.mesh.st_len_x::Int
    model.mesh.st_len_u::Int

    outer_mesh::Vector{Float64}
    inner_x::Vector{Float64}
    inner_u::Vector{Float64}

    quad_pts::Vector{Float64}
    quad_weights::Vector{Float64}

    quad_pts_check::Vector{Float64}
    quad_weights_check::Vector{Float64}

    evalQX::Vector{Float64}
    evalQU::Vector{Float64}
    evalXU::Vector{Float64}
    evalUX::Vector{Float64}

    Dx::Matrix{Float64}
    Tapir_mesh(N,Px,Pu,Q) = new(N, Px, Pu, Q) 
end

function set_length!(model)
    model.mesh.st_len_x=model.nx*(model.mesh.Px+1);
    model.mesh.st_len_u=model.nu*(model.mesh.Pu+1);

    if model.settings.cont_type_x=="samevar"
        model.mesh.st_len_x=st_len_x-model.nx;
    end

    if model.settings.cont_type_u=="samevar"
        model.mesh.st_len_u=st_len_u-model.nu;
    end
end

function generate_quad!(mesh, quad_type; check_Q=2*mesh.Q)
    if quad_type=="gausslegendre"
        mesh.quad_pts,mesh.quad_weights=gausslegendre(mesh.Q)
        mesh.quad_pts_check,mesh.quad_weights_check=gausslegendre(check_Q)
    elseif quad_type=="gaussradau"
        mesh.quad_pts,mesh.quad_weights=gaussradau(mesh.Q)
        mesh.quad_pts_check,mesh.quad_weights_check=gausslegendre(check_Q)
    elseif quad_type=="gausslobatto"
        mesh.quad_pts,mesh.quad_weights=gausslobatto(mesh.Q)
        mesh.quad_pts_check,mesh.quad_weights_check=gausslegendre(check_Q)
    elseif quad_type=="gausskronrod"
        mesh.quad_pts_check,mesh.quad_weights_check,mesh.quad_weights=kronrod(mesh.Q)
        mesh.quad_pts=mesh.quad_pts_check[2:2:end];
    else
        error("Quadrature type not supported")
    end
end

function generate_outer!(mesh, outer_mesh_start)
    if outer_mesh_start=="equidistant"
        mesh.outer_mesh=collect(0:1/mesh.N:1);
    else
        error("Mesh type not supported")
    end
end

function generate_interpolation_mesh!(mesh,inner_mesh_type)
    if inner_mesh_type=="cheb2"
        mesh.inner_x=cos.((mesh.Px:-1:0).*π./mesh.Px);
        mesh.inner_u=cos.((mesh.Pu:-1:0).*π./mesh.Pu);
    elseif inner_mesh_type=="cheb1"
        mesh.inner_x=[cos((2 * k - 1) * π / (2 * (mesh.Px+1))) for k = mesh.Px+1:-1:1]
        mesh.inner_u=[cos((2 * k - 1) * π / (2 * (mesh.Pu+1))) for k = mesh.Pu+1:-1:1]
    elseif inner_mesh_type=="equi"
        mesh.inner_x=collect(0:1/mesh.Px:1);
        mesh.inner_u=collect(0:1/mesh.Pu:1);
    elseif inner_mesh_type=="legendre"
        mesh.inner_x,~=gausslegendre(mesh.Px+1);
        mesh.inner_u,~=gausslegendre(mesh.Pu+1);
    elseif inner_mesh_type=="lobatto"
        mesh.inner_x,~=gausslobatto(mesh.Px+1);
        mesh.inner_u,~=gausslobatto(mesh.Pu+1);
    else
        error("Interpolation point distribution not supported")
    end
end

function generate_interpolation!(mesh,interpolation_type)
    if interpolation_type=="lagrange_bary"
        mesh.evalQX=L(mesh.quad_pts,mesh.inner_x);
        mesh.evalQU=L(mesh.quad_pts,mesh.inner_u);
        mesh.evalUX=L(mesh.inner_u,mesh.inner_x);
        mesh.evalXU=L(mesh.inner_x,mesh.inner_u);

        mesh.Dx=D(mesh.inner_x);
    elseif interpolation_type=="lagrange"
        display("Not yet supported")
    elseif interpolation_type=="monomials"
        display("Not yet supported")
    elseif interpolation_type=="bernstein"
        display("Not yet supported")
    else
        error("Interpolation point distribution not supported")
    end
end


function L(x,ip)
    return [prod(x[ix] .- ip[ collect( Iterators.filter(j -> j != i, 1:length(ip)) ) ] ) /  prod(ip[i] .- ip[ collect( Iterators.filter(j -> j != i, 1:length(ip)) ) ] ) for i in eachindex(ip), ix in eachindex(x)]
end

function D(ip)
    return [i==j ? sum(1 ./ (ip[i] .- ip[ collect(Iterators.filter(k -> k != i, eachindex(ip)))] ) ) : prod(ip[j] .- ip[ collect(Iterators.filter(k -> (k != i && k != j), eachindex(ip)))]) / prod(ip[i] .- ip[ collect(Iterators.filter(k -> (k != i), eachindex(ip) ))]) for i in eachindex(ip), j in eachindex(ip)]
end
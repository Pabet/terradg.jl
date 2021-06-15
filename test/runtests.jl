using TerraDG
using Test
using LinearAlgebra

@testset "All tests" begin

@testset "Quadrature on the reference element" begin
    n = 5
    basis = TerraDG.Basis(n,1)
    p = basis.quadpoints
    wp = basis.quadweights
    func = (x,y) -> x*y*y
    sum = 0.0
    for i in CartesianIndices((n,n))
        sum += func(p[i[1]],p[i[2]]) * wp[i[1]]*wp[i[2]]
    end
    @test abs(sum - 1. / 6.) <= 10e-15
end

@testset "Maps reference coordinate to correct physical coord" begin
    cellcenter = [1.0, 1.0]
    cellsize = [0.5, 0.5]
    cell = TerraDG.Cell(cellcenter, 
            cellsize,
            [],
            [],
            -1)

    @test TerraDG.globalposition(cell, [0.5, 0.5]) == cellcenter
    @test TerraDG.globalposition(cell, [0.0, 0.0]) == [0.75, 0.75]
    @test TerraDG.globalposition(cell, [1.0, 1.0]) == [1.25, 1.25]
    @test_throws BoundsError TerraDG.globalposition(cell, [1.1, 0.0])
    @test_throws BoundsError TerraDG.globalposition(cell, [-0.1, 0.0])
end

@testset "Global to reference to global is same" begin
    cellcenter = [1.0, 1.0]
    cellsize = [0.5, 0.5]
    cell = TerraDG.Cell(cellcenter, 
            cellsize,
            [],
            [],
            -1)

    coords = [
        [0.5, 0.5],
        [0.0, 0.0],
        [1.0, 1.0]
    ]
    for coord ∈ coords
        global_coord = TerraDG.globalposition(cell, coord)
        local_coord = TerraDG.localposition(cell, global_coord)
        @test local_coord == coord
    end
end


@testset "Volume of 2D cells is correct" begin
    sizes = [
        [1.0, 1.0],
        [2.0, 2.0],
        [0.5, 0.5]
    ] 

    volumes = [
        1.0,
        4.0,
        0.25
    ]
    for (size, volume) ∈ zip(sizes, volumes)
        cell = TerraDG.Cell([0.0, 0.0], 
            size,
            [],
            [],
            -1)
        @test TerraDG.volume(cell) == volume
    end
end

@testset "Area of 2D cells is correct" begin
    sizes = [
        [1.0, 1.0],
        [2.0, 2.0],
        [0.5, 0.5]
    ] 

    areas = [
        1.0,
        2.0,
        0.5
    ]
    for (size, area) ∈ zip(sizes, areas)
        cell = TerraDG.Cell([0.0, 0.0], 
            size,
            [],
            [],
            -1)
        @test TerraDG.area(cell) == area
    end
end

@testset "Lagrange polynomials are correct" begin
    δ(a,b) = if (a == b) 1.0 else 0.0 end
    for n in range(1, length=6)
        basis = TerraDG.Basis(n, 1)
        points = basis.quadpoints
        for (i,p) in enumerate(points)
            for (j,q) in enumerate(points)
                @test TerraDG.lagrange_1d(points, i, q) == δ(i,j)
            end
        end
    end
end

@testset "Lagrange polynomials sum to one" begin
    for n in range(1, length=6)
        for x in [0.0, 0.5, 1.0]
            basis = TerraDG.Basis(n, 1)
            points = basis.quadpoints
            sum = 0.0
            for (i,p) in enumerate(points)
                sum += TerraDG.lagrange_1d(points, i, x)
            end
            @test isapprox(sum, 1, atol=10^-14)
        end
    end
end

@testset "Integral of derivatives of 1D-Lagrange polynomials is correct" begin
    for n = 1:6, i = 1:n, x = 0:0.1:1
        basis = TerraDG.Basis(n, 1)
        roots = basis.quadpoints
        integral_exact = TerraDG.lagrange_1d(roots, i, x) - TerraDG.lagrange_1d(roots, i, 0.0)
        integral_quad = sum(j ->  x * basis.quadweights[j] * TerraDG.lagrange_diff(roots, i, x * basis.quadpoints[j]), eachindex(basis.quadpoints))

        @test isapprox(integral_quad, integral_exact, atol=10e-14)
    end
end

@testset "Integral of derivatives of 2D-Lagrange polynomials is also correct" begin
    for n = 1:6
        basis = TerraDG.Basis(n, 2)
        qps = basis.quadpoints
        qws = basis.quadweights

        roots = qps

        cartesians = CartesianIndices(size(basis))
        n_phis = length(basis)
        for n_phi = 1:n_phis, L = 0:0.2:1
            i = cartesians[n_phi][1]
            j = cartesians[n_phi][2]

            Phi = function (x, y)
                coeffs = zeros(n_phis)
                coeffs[n_phi] = 1
                TerraDG.evaluate_basis(basis, coeffs, [x, y])
            end

            xderiv = (x,y) -> TerraDG.lagrange_diff(roots, i, x) * TerraDG.lagrange_1d(roots, j, y)
            yderiv = (x,y) -> TerraDG.lagrange_1d(roots, i, x) * TerraDG.lagrange_diff(roots, j, y)
            
            lhs_x = sum(n -> qws[n] * (Phi(L, qps[n]) - Phi(0, qps[n])), eachindex(qps))
            lhs_y = sum(n -> qws[n] * (Phi(qps[n], L) - Phi(qps[n], 0)), eachindex(qps))

            rhs_x = 0
            rhs_y = 0

            for n1 = eachindex(qps), n2 = eachindex(qps)
                rhs_x += L * qws[n1] * qws[n2] * xderiv(L*qps[n1], qps[n2])
                rhs_y += L * qws[n1] * qws[n2] * yderiv(qps[n1], L*qps[n2])
            end

            @test isapprox(rhs_x, lhs_x, atol=10e-14)
            @test isapprox(rhs_y, lhs_y, atol=10e-14)
        end
    end
end

@testset "Derivative of const = 0" begin
    for n=1:6
        basis = TerraDG.Basis(n, 1)
        for x in [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
            @test isapprox(
            sum([TerraDG.lagrange_diff(basis.quadpoints, i, x) for i=1:size(basis, 1)]), 0, atol=10e-14)
        end
    end
end

@testset "Mass matrix is correct for 2D" begin
    basis_order1 = TerraDG.Basis(1, 2)
    basis_order2 = TerraDG.Basis(2, 2)
    @test TerraDG.massmatrix(basis_order1, 2) == reshape([1], 1, 1)
    massmatrix_order2 = Diagonal([0.25,0.25,0.25,0.25])
    @test TerraDG.massmatrix(basis_order2, 2) == massmatrix_order2
end

@testset "Projection/evaluation works for polynomials" begin
    ns = [1,2,3,4,5,6]
    funs = [
        (x, y) -> [1],
        (x, y) -> [x + y],
        (x, y) -> [(x + y)^2],
        (x, y) -> [(x + y)^3],
        (x, y) -> [(x + y)^4],
        (x, y) -> [(x + y)^5]
    ]
    points = [
        [0.21, 0.23],
        [0.38, 0.93],
        [0.92, 0.23],
        [0.01, 0.01],
        [0.99, 0.99],
    ]
    proj_reshaped(basis, func) = reshape(
        TerraDG.project_to_reference_basis(func, basis, 1), length(basis.quadweights)^2 )

    for (n, fun) in zip(ns, funs)
        basis = TerraDG.Basis(n, 2)
        coeffs = proj_reshaped(basis, fun)
        for point in points
            evaluated = TerraDG.evaluate_basis(basis, coeffs, point)
            reference = fun(point[1], point[2])[1]
            @test isapprox(evaluated, reference, atol=10e-14)
        end
    end
end

@testset "Derivative matrix is correct" begin
    ns = [1,2,3,2,3,4,5,6]
    funs = [
        (x, y) -> [1]
        (x, y) -> [1]
        (x, y) -> [1]
        (x, y) -> [x + y]
        (x, y) -> [(x + y)^2]
        (x, y) -> [(x + y)^3]
        (x, y) -> [(x + y)^4]
        (x, y) -> [(x + y)^5]
    ]
    funs_deriv_x = [
        (x,y) -> [0]
        (x,y) -> [0]
        (x,y) -> [0]
        (x,y) -> [1]
        (x,y) -> [2 * (x + y)]
        (x,y) -> [3 * (x + y)^2]
        (x,y) -> [4 * (x + y)^3]
        (x,y) -> [5 * (x + y)^4]
    ]
    funs_deriv_y = [
        (x,y) -> [0]
        (x,y) -> [0]
        (x,y) -> [0]
        (x,y) -> [1]
        (x,y) -> [2 * (x + y)]
        (x,y) -> [3 * (x + y)^2]
        (x,y) -> [4 * (x + y)^3]
        (x,y) -> [5 * (x + y)^4]
    ]

    for (n, fun, deriv_x, deriv_y) in zip(ns, funs, funs_deriv_x, funs_deriv_y)
        proj_reshaped(basis, func) = reshape(
            TerraDG.project_to_reference_basis(func, basis, 1), length(basis.quadweights)^2 )
        basis = TerraDG.Basis(n, 2)
        fun_evaluated = proj_reshaped(basis, fun)
        reference_derivx_evaluated = proj_reshaped(basis, deriv_x)
        reference_derivy_evaluated = proj_reshaped(basis, deriv_y)
        ∇ = TerraDG.derivativematrix(basis)
        derivx_evaluated = transpose(∇[:, 1:length(basis)]) * fun_evaluated
        derivy_evaluated = transpose(∇[:, length(basis) + 1 : end]) * fun_evaluated

        @test all(isapprox.(derivx_evaluated, reference_derivx_evaluated, atol=10e-12))
        @test all(isapprox.(derivy_evaluated, reference_derivy_evaluated, atol=10e-12))
    end
end


@testset "Face projection matrix is correct" begin
    ns = [1,2,3,4,5,6]
    funs = [
        (x, y) -> [1,1,1]
        (x, y) -> [0,0,0]
        (x, y) -> [0,1,0]
    ]
    for n in ns
        basis = TerraDG.Basis(n, 2)
        for func in funs
            func_proj = TerraDG.project_to_reference_basis(func, basis, 3)
            for face in [TerraDG.top, TerraDG.bottom, TerraDG.left, TerraDG.right]
                fp = TerraDG.face_projection_matrix(basis, TerraDG.top)
                func_proj_face = fp * func_proj
                @test sum(func_proj)/length(func_proj) ≈ sum(func_proj_face)/length(func_proj_face)
            end
        end
    end
end

@testset "Linearization is correct" begin
    ns = [1,2,3,4,5,6]
    average = 3
    x_slope = -1.5
    y_slope = 4
    linearf = (x, y) -> [average + x_slope * (x - 0.5) + y_slope * (y - 0.5)]
    
    for n in ns
        basis = TerraDG.Basis(n, 2)
        globals = TerraDG.GlobalMatrices(basis, TerraDG.IdentityFilter(), 2)

        ref_proj = TerraDG.project_to_reference_basis(linearf, basis, 1)
        proj = similar(ref_proj)

        TerraDG.linearize(globals, average, x_slope, y_slope, proj)
        @test all(proj .≈ ref_proj)
    end
end

@testset "Projection to 'slope' basis is correct" begin
    ns = [2,3,4,5,6]
    fs = [
        (x, y) -> [1, x - 0.5, y - 0.5],
    ]

    averages = [
        [1.0, 0.0, 0.0]
    ]

    slopes_x = [
        [0, 1, 0]
    ]

    slopes_y = [
        [0, 0, 1]
    ]

    for n in ns
        basis = TerraDG.Basis(n, 2)
        globals = TerraDG.GlobalMatrices(basis, TerraDG.IdentityFilter(), 2)

        for (f, avg_ref, sx_ref, sy_ref) in zip(fs, averages, slopes_x, slopes_y)
            f_proj = TerraDG.project_to_reference_basis(f, basis, size(avg_ref, 1))

            avg = similar(avg_ref)
            TerraDG.project_to_basis(globals, avg, f_proj)
            @test all(isapprox.(avg, avg_ref; atol=10^-14))

            sx = similar(avg)
            TerraDG.project_to_basis(globals, sx, f_proj; f = ((x, _), ) -> x - 0.5, c = 12)
            @test all(isapprox.(sx, sx_ref; atol=10^-14))

            sy = similar(avg)
            TerraDG.project_to_basis(globals, sy, f_proj; f = ((_, y), ) -> y - 0.5, c = 12)
            @test all(isapprox.(sy, sy_ref; atol=10^-14))
        end
    end
end
end

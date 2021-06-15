using FastGaussQuadrature

"""
    lagrange_1d(points, i, x)

Evaluate the Lagrange interpolation polynomial defined over nodal `points` with index `i`
at point `x`.
"""
function lagrange_1d(points, i, x)
    result = 1.0
    for j = eachindex(points)
        if i != j
            result *= (x - points[j])/(points[i] - points[j])
        end
    end
    result
end

"""
    lagrange_diff(points, i, x)

Evaluate the derivative of the Lagrange interpolation polynomial defined over nodal `points` 
with index `i` at point `x`.
"""
function lagrange_diff(points, i, x)
    result = 0.0
    for j = eachindex(points)
        if i != j
            acc = 1 / (points[i] - points[j])
            for m = eachindex(points)
                if m != i && m != j
                    acc *= (x - points[m]) / (points[i] - points[m])
                end
            end
            result += acc
        end
    end

    result
end

"""
    get_quadpoints(n)

Compute quadrature points and weights for Gaussian quadrature
of order `n`.
The points (and thus the corresponding weights) are normalized to the range
``[0.0, 1.0]``.

Return a tuple of `points, weights`.
"""
function get_quadpoints(n)
    x, w = gausslegendre(n)
    # They are scaled in [-1.0, 1.0]
    # and we need [0.0, 1.0]
    (x .+ 1)./2, w ./ 2
end

"""
    Basis

A standard 1-dimensional basis of `order::Integer`
with `quadpoints::Array{Float64,1}` and 
corresponding `quadweights::Array{Float64,1}`
Is the basis (pun intended) for tensor-product
bases.

    Basis(order::Integer, dimensions)

Initialize a basis of `order::Integer` and
dimensions `dimensions`.
"""
struct Basis
    quadpoints::Array{Float64,1}
    quadweights::Array{Float64,1}
    order::Int64
    dimensions::Int64

    function Basis(order, dimensions)
        quadpoints, quadweights = get_quadpoints(order)
        new(quadpoints, quadweights, order, dimensions)
    end
end


"""
    Base.length(basis::Basis)

Return number of points for `basis` in n-dimensions.
"""
Base.length(basis::Basis) = basis.order ^ basis.dimensions

"""
    Base.size(basis::Basis)

Return number of points for `basis` for each dimensions as tuple.
"""
function Base.size(basis::Basis) 
    ntuple(basis.dimensions) do x 
       basis.order
    end
end

"""
    Base.size(basis::Basis, dim)

Return number of points for `basis` for dimensions `dim`.
"""
function Base.size(basis::Basis, dim::Integer)
    basis.order
end

"""
    evaluate_basis(basis::Basis, coeffs, x)

Evaluate the `basis` with coefficients
`coeffs` at point `x`.
"""
function evaluate_basis(basis::Basis, coeffs, x)
    cartesian = CartesianIndices(size(basis))
    
    sum(i -> coeffs[i] * evaluate_phi(basis, x, cartesian[i]), 1:length(basis))
end

function evaluate_phi(basis::Basis, x, cindex::CartesianIndex)
    prod(i -> lagrange_1d(basis.quadpoints, cindex[i], x[i]), 1:basis.dimensions)
end

"""
    project_to_reference_basis(fun, basis::Basis, ndofs)

Project the result of the function `fun` to coefficients
of the basis built of a tensor-product of `basis`.
The function `fun(x,y)`  takes in the ``x, y``-coordinates
and returns a vector with size `ndofs`.
The corresponding coefficients are returned.
"""
function project_to_reference_basis(fun, basis::Basis, ndofs)
    if basis.dimensions == 2
        cartesian = CartesianIndices(size(basis))
        coeffs = zeros((length(basis), ndofs))

        for i in 1:length(basis)
            xi = cartesian[i][1]
            yi = cartesian[i][2]
            coeffs[i,:] .= fun(basis.quadpoints[xi], basis.quadpoints[yi])
        end

        return coeffs
    else
        error(string("Number of dimensions not supported: ", basis.dimensions))
    end
end


"""
    massmatrix(basis, dimensions)

Return the mass-matrix for a `dimensions`-dimensional
tensor-product basis built up from the 1d-basis `basis`.
"""
function massmatrix(basis, dimensions)
    cartesians = CartesianIndices(ntuple(_ -> basis.order, dimensions))
    diag = zeros(basis.order ^ dimensions)

    for i = eachindex(diag)
        diag[i] = prod(j -> basis.quadweights[cartesians[i][j]], 1:dimensions)
    end
    Diagonal(diag)
end

"""
    derivativematrix(basis)

Returns the 2-dimensional derivative matrix for `basis`.
Multiplying this with flux-coefficients of shape
`(dimensions * basissize_2d, ndofs)` returns the
coefficients of the corresponding derivative.
"""
function derivativematrix(basis)
    D = zeros(length(basis), basis.dimensions * length(basis))
    cartesians = CartesianIndices(size(basis))
    for j = 1:length(basis), p = 1:length(basis), d = 1:basis.dimensions
        D[j, p + (d-1)*length(basis)] = evaluate_derivative(basis, basis.quadpoints, cartesians[j], cartesians[p], d)
    end
    D
end

function evaluate_derivative(basis, x, fun_index, qp_index, deriv_dim)
    deriv = lagrange_diff(basis.quadpoints, fun_index[deriv_dim], x[qp_index[deriv_dim]])

    for i = 1:basis.dimensions
        if i == deriv_dim
            continue
        else
            deriv *= lagrange_1d(basis.quadpoints, fun_index[i], x[qp_index[i]])
        end
    end
    deriv
end
"""
    get_face_quadpoints(basis::Basis, face)

Return the quadrature points at the face `face` for basis `basis`.
"""
function get_face_quadpoints(basis::Basis, face)
    to_face = Dict(
        left => x -> (0, x),
        right => x -> (1, x),
        top => x -> (x, 1),
        bottom => x -> (x, 0)
    )

    map(to_face[face], basis.quadpoints)
end

"""
    face_projection_matrix(basis, face)

Return the face projection matrix for `basis` and `face`.
Multiplying it with coefficient vector for the right basis
returns the coefficients of the solution evaluated at the 
quadrature nodes of the face.
"""
function face_projection_matrix(basis, face)
    P = ones((basis.order, length(basis)))

    qps = get_face_quadpoints(basis, face)

    coeffs = zeros(length(basis))
    for j = 1:length(basis)
        coeffs[j] = 1
        for i = 1:basis.order
            P[i, j] = evaluate_basis(basis, coeffs, qps[i])
        end
        coeffs[j] = 0
    end

    P
end

"""
    evaluate_m_to_n_vandermonde_basis(basis) 

Return the Vandermonde matrix that converts between the 2D-modal 
(normalized) Legendre-basis and the 2D-nodal tensor-product basis built
with `basis`.
"""
function evaluate_m_to_n_vandermonde_basis(basis)
    Array{Float64}(I, length(basis), length(basis))
end

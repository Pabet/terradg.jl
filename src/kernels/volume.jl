"""
    BuffersVolume(basis, ndofs)

This gives buffers that are used to store intermediate
results during `evaluate_volume`.
"""
struct BuffersVolume
    scaled_fluxcoeff::Array{Float64,2}
    χ::Diagonal{Float64, Array{Float64, 1}} 

    function BuffersVolume(basis, ndofs)
        basissize_nd = length(basis)

        scaled_fluxcoeff = zeros(basissize_nd * 2,ndofs)
        scaling = Diagonal(zeros(basissize_nd * 2, basissize_nd * 2))

        new(scaled_fluxcoeff, scaling)
    end

end

"""
    evaluate_volume(globals, buffers, flux_coeff, basis, inverse_jacobian, volume, celldu)

Evaluates the volume term of our pde.

# Arguments
- `globals`: struct containing global matrices
- `buffers`: buffers to store intermediate results
- `flux_coeff`: coefficients of the flux
- `basis`: 1d-basis
- `inverse_jacobian`: inverse of the Jacobian matrix. Assumed to be diagonal.
- `volume`: n-dimensional volume of the cell. E.g., area of 2d-cell
- `celldu`: update of dofs, return value is added to this
"""
function evaluate_volume(globals, buffers, flux_coeff, basis, inverse_jacobian, volume, celldu)
    quadweights = globals.quadweights_nd
    for i = 1:length(basis)
        buffers.χ[i,i] = volume * inverse_jacobian[1,1] * quadweights[i]
        buffers.χ[length(basis) + i, length(basis) + i] = volume * inverse_jacobian[2,2] * quadweights[i]
    end

    buffers.scaled_fluxcoeff .= buffers.χ * flux_coeff
    celldu .+= globals.reference_derivative_matrix * buffers.scaled_fluxcoeff
end

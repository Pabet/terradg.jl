struct MHD <: Equation end
@declare_dofs MHD [:rho, :rhou, :rhov, :rhow, :rhoE, :Bx, :By, :Bz, :psi] 

struct IdealRotor <: Scenario end
struct KelvinHelmholtz <: Scenario end
struct OrszagTang <: Scenario end
struct AlfenWave <: Scenario end

function is_periodic_boundary(eq::MHD, scenario::IdealRotor)
    false
end

function is_periodic_boundary(eq::MHD, scenario::KelvinHelmholtz)
    true
end

function is_periodic_boundary(eq::MHD, scenario::OrszagTang)
    true
end

function is_periodic_boundary(eq::MHD, scenario::GaussianWave)
    true
end

function is_periodic_boundary(eq::MHD, scenario::AlfenWave)
    true
end

function evaluate_boundary(eq::MHD, scenario::IdealRotor, face, normalidx, dofsface, dofsfaceneigh)
    dofsfaceneigh .= dofsface
end

function evaluate_boundary(eq::MHD, scenario::Scenario, face, normalidx, dofsface, dofsfaceneigh)
    #temporarily keep the reflecting boundaries of the euler equations
    s = MHDShortcuts()
    dofsfaceneigh .= dofsface
    if normalidx == 1
        dofsfaceneigh[:,s.rhou] .= -dofsface[:,s.rhou]
        #dofsfaceneigh[:,s.Bx] .= -dofsface[:,s.Bx]
    else
        dofsfaceneigh[:,s.rhov] .= -dofsface[:,s.rhov]
        #dofsfaceneigh[:,s.By] .= -dofsface[:,s.By]
    end
end

function get_γ(scenario::Scenario)
    1.4
end

function get_γ(scenario::OrszagTang)
    5/3
end

function get_γ(scenario::KelvinHelmholtz)
    5/3
end

function get_initial_values(eq::MHD, scenario::IdealRotor, global_pos; t=0.0)
    px, py = global_pos
    # the Dumbser paper uses domain [-0.5, 0.5]²
    x = px - 0.5
    y = py - 0.5
    # rotor has radius 0.1
    r = sqrt(x^2 + y^2)
    taper_width = 0.03 # 6 cells (at grid 200)
    taper_end = 0.1 + taper_width
    if r <= 0.1
        ρ = 10
        u, v, w = cross([x, y, 0], [0, 0, 10])
    # linear taper at border
    elseif 0.1 < r <= taper_end
        ρ = (10*(taper_end - r) + 1*(r-0.1))/(taper_end - 0.1)
        vel_at_radius = cross(normalize([x, y, 0]) * 0.1, [0, 0, 10])
        u, v, w = (vel_at_radius * (taper_end - r)) / (taper_end - 0.1)
    else
        ρ = 1
        u = v = w = 0
    end
    p = 1
    ρu = ρ*u
    ρv = ρ*v
    ρw = ρ*w
    Bx = 2.5
    By = 0
    Bz = 0
    γ = get_γ(scenario)
    E = evaluate_energy(eq, γ, p, ρu, ρv, ρw, Bx, By, Bz, inv(ρ))
    ψ = 0

    # TODO velocity "2 units" in Balsara paper, Dumbser leads to 1 unit at 0.1
    #      but in Balsara the magnetic field has value 5 instead of 2.5?
    [ρ, ρu, ρv, ρw, E, Bx, By, Bz, ψ]
end

function get_initial_values(eq::MHD, scenario::KelvinHelmholtz, global_pos; t=0.0)
    px, py = global_pos
    # the Dumbser paper uses domain [0,2]×[-1,1]
    x = 2*px
    y = 2*(py - 0.5)

    ρ = 1
    p = 0.6
    a = 1/25
    U₀ = 1
    δv = 0.01
    B₀ = 0.07
    χ = pi * (y - 0.5 + a) / (4*a)

    u = -0.5 * U₀ * tanh((abs(y) - 0.5)/a)
    v = δv * sinpi(2*x) * sinpi(abs(y))
    w = 0

    if 0.5 + a < abs(y) < 1
        Bx = B₀
        By = Bz = 0
    elseif 0.5 - a < abs(y) < 0.5 + a
        Bx = B₀ * sin(χ)
        By = 0
        Bz = cos(χ)
    else # 0 < abs(y) < 0.5 - a
        Bx = By = 0
        Bz = B₀
    end

    ρu = ρ*u
    ρv = ρ*v
    ρw = ρ*w

    γ = get_γ(scenario)
    E = evaluate_energy(eq, γ, p, ρu, ρv, ρw, Bx, By, Bz, inv(ρ))

    ψ = 0
    [ρ, ρu, ρv, ρw, E, Bx, By, Bz, ψ]
end

function get_initial_values(eq::MHD, scenario::OrszagTang, global_pos; t=0.0)
    x, y = global_pos
    γ = get_γ(scenario)
    ρ = γ^2
    p = γ

    ρu = ρ * -sin(y)
    ρv = ρ * sin(x)
    ρw = 0

    Bx = sqrt(4*pi)*-sin(y)
    By = sqrt(4*pi)*sin(2*x)
    Bz = 0

    E = evaluate_energy(eq, γ, p, ρu, ρv, ρw, Bx, By, Bz, inv(ρ))
    [ρ, ρu, ρv, ρw, E, Bx, By, Bz, 0]
end

function get_initial_values(eq::MHD, scenario::GaussianWave, global_pos; t=0.0)
    pxg, pyg = global_pos
    p = exp(-100 * (pxg - 0.5)^2 - 100 *(pyg - 0.5)^2) + 1
    rho = 1.0
    rhoE = evaluate_energy(eq, get_γ(scenario), p, 0, 0, 0, 1, 0, 0, inv(rho))
    [rho,0,0,0,rhoE,1,0,0,0]
end

function get_initial_values(eq::MHD, scenario::AlfenWave, global_pos; t=0.0)
    x, y = global_pos

    ρ = 1
    p = 1
    η = 1
    B₀ = 1

    γ = get_γ(scenario)

    hh = 1 + γ / (γ - 1) * p / ρ
    tempaa = ρ * hh + B₀^2 * (1 + η^2)
    tempab = 2 * η * B₀^2 / tempaa
    tempac = 0.5 * (1.0 + sqrt(1 - tempab^2))
    va2 = B₀^2 / (tempaa * tempac)
    vax = sqrt(va2)

    Bx = B₀
    By = η * B₀ * cospi(2 * (x - vax*t))
    Bz = η * B₀ * sinpi(2 * (x - vax*t))

    u = 0
    v = -vax * By / B₀
    w = -vax * Bz / B₀

    E = evaluate_energy(eq, γ, p, ρ*u, ρ*v, ρ*w, Bx, By, Bz, inv(ρ))

    [ρ, ρ*u, ρ*v, ρ*w, E, Bx, By, Bz, 0]
end

function is_analytical_solution(eq::MHD, scenario::AlfenWave)
    true
end

function is_analytical_solution(eq::MHD, scenario::Scenario)
    false
end

function evaluate_flux(eq::MHD, scenario::Scenario, celldofs, cellflux)
    s = MHDShortcuts()
    x_range = 1:size(celldofs, 1)
    y_range = x_range .+ size(celldofs, 1)

    eval_pressure = (E, ρu, ρv, ρw, Bx, By, Bz, inv_rho) -> evaluate_pressure(eq, get_γ(scenario), E, ρu, ρv, ρw, Bx, By, Bz, inv_rho)
    
    inv_rho = inv.(celldofs[:,s.rho])
    pi4_inv = inv(4*pi)
    p = eval_pressure.(celldofs[:,s.rhoE], celldofs[:,s.rhou], celldofs[:,s.rhov], celldofs[:,s.rhow], celldofs[:,s.Bx], celldofs[:,s.By], celldofs[:,s.Bz], inv_rho)

    # ρu
    cellflux[x_range, s.rho] .= celldofs[:,s.rhou]
    # ρu² + p + 1/2 (By² + Bz² - Bx²)
    cellflux[x_range, s.rhou] .= celldofs[:,s.rhou].^ 2 .* inv_rho + p + 0.5 * pi4_inv * (celldofs[:,s.By].^2 + celldofs[:,s.Bz].^2 - celldofs[:,s.Bx].^2)
    # ρu*v - Bx*By
    cellflux[x_range, s.rhov] .= celldofs[:,s.rhou] .* celldofs[:,s.rhov] .* inv_rho - celldofs[:,s.Bx] .* celldofs[:,s.By] * pi4_inv
    # ρu*w - Bx*Bz
    cellflux[x_range, s.rhow] .= celldofs[:,s.rhou] .* celldofs[:,s.rhow] .* inv_rho - celldofs[:,s.Bx] .* celldofs[:,s.Bz] * pi4_inv
    # ψ
    cellflux[x_range, s.Bx] .= celldofs[:,s.psi]
    # By*u - Bx*v
    cellflux[x_range, s.By] .= inv_rho .* (celldofs[:,s.By] .* celldofs[:,s.rhou] - celldofs[:,s.Bx] .* celldofs[:,s.rhov])
    # Bz*u - Bx*w
    cellflux[x_range, s.Bz] .= inv_rho .* (celldofs[:,s.Bz] .* celldofs[:,s.rhou] - celldofs[:,s.Bx] .* celldofs[:,s.rhow])
    # u*(rhoE + p + 0.5*B²) - Bx*(u*Bx + v*By + w*Bz)
    cellflux[x_range, s.rhoE] .= inv_rho .* (celldofs[:,s.rhou] .* (celldofs[:,s.rhoE] .+ p .+ 0.5 * pi4_inv * (celldofs[:,s.Bx].^2 + celldofs[:,s.By].^2 + celldofs[:,s.Bz].^2))
                                - pi4_inv * celldofs[:,s.Bx] .* (celldofs[:,s.rhou] .* celldofs[:,s.Bx] + celldofs[:,s.rhov] .* celldofs[:,s.By] + celldofs[:,s.rhow] .* celldofs[:,s.Bz]))
    # cₕ² * Bx
    cellflux[x_range, s.psi] .= celldofs[:,s.Bx]

    # ρv
    cellflux[y_range, s.rho] .= celldofs[:,s.rhov]
    # ρv*u - By*Bx
    cellflux[y_range, s.rhou] .= celldofs[:,s.rhov] .* celldofs[:,s.rhou] .* inv_rho - celldofs[:,s.By] .* celldofs[:,s.Bx] * pi4_inv
    # ρv² + p + 1/2 (Bx² + Bz² - By²)
    cellflux[y_range, s.rhov] .= celldofs[:,s.rhov].^2 .* inv_rho + p + 0.5 * pi4_inv * (celldofs[:,s.Bx].^2 + celldofs[:,s.Bz].^2 - celldofs[:,s.By].^2)
    # ρv*w - By*Bz
    cellflux[y_range, s.rhow] .= celldofs[:,s.rhov] .* celldofs[:,s.rhow] .* inv_rho - celldofs[:,s.By] .* celldofs[:,s.Bz] * pi4_inv
    # Bx*v - By*u
    cellflux[y_range, s.Bx] .= inv_rho .* (celldofs[:,s.Bx] .* celldofs[:,s.rhov] - celldofs[:,s.By] .* celldofs[:,s.rhou])
    # ψ
    cellflux[y_range, s.By] .= celldofs[:,s.psi]
    # Bz*v - By*w
    cellflux[y_range, s.Bz] .= inv_rho .* (celldofs[:,s.Bz] .* celldofs[:,s.rhov] - celldofs[:,s.By] .* celldofs[:,s.rhow])
    # v*(rhoE + p + 0.5*B²) - By*(u*Bx + v*By + w*Bz)
    cellflux[y_range, s.rhoE] .= inv_rho .* (celldofs[:,s.rhov] .* (celldofs[:,s.rhoE] .+ p .+ 0.5 * pi4_inv * (celldofs[:,s.Bx].^2 + celldofs[:,s.By].^2 + celldofs[:,s.Bz].^2))
                                - pi4_inv * celldofs[:,s.By] .* (celldofs[:,s.rhou] .* celldofs[:,s.Bx] + celldofs[:,s.rhov] .* celldofs[:,s.By] + celldofs[:,s.rhow] .* celldofs[:,s.Bz]))
    # cₕ² * Bx
    cellflux[y_range, s.psi] .= celldofs[:,s.By]
end

function max_eigenval(eq::MHD, scenario::Scenario, celldata, normalidx)
    max_eigenval(eq, celldata, normalidx; γ = get_γ(scenario))
end

function max_eigenval(eq::MHD, scenario::KelvinHelmholtz, celldata, normalidx)
    max_eigenval(eq, celldata, normalidx; γ = 5/3)
end

function max_eigenval(eq::MHD, celldata, normalidx; γ = 1.4)
    s = MHDShortcuts()
    v_idx = normalidx == 1 ? s.rhou : s.rhov
    max_eigenval = abs(celldata[1,v_idx]) / celldata[1,s.rho]
    for i = 1:(size(celldata,1))
        inv_rho = inv(celldata[i, s.rho])
        u = celldata[i, v_idx] * inv_rho
        ca = celldata[i, s.Bx]/sqrt(4*pi*celldata[i, s.rho])
        p = evaluate_pressure(eq::MHD, γ, celldata[i, s.rhoE], celldata[i, s.rhou], celldata[i, s.rhov], celldata[i, s.rhow], celldata[i, s.Bx], celldata[i, s.By], celldata[i, s.Bz], inv_rho)
        B² = celldata[i, s.Bx]^2 + celldata[i, s.By]^2 + celldata[i, s.Bz]^2

        c = sqrt((γ*p) * inv_rho)
        b = sqrt(B²/(4*pi*celldata[i, s.rho]))
        #cs = sqrt(0.5*(b + c - sqrt((b+c)^2 - 4*ca*c)))        
        cf = sqrt(0.5*(b + c + sqrt((b+c)^2 - 4*ca*c)))
        #max_eigenval = max(max_eigenval, u + cf, u - cf, u + ca, u - ca, u + cs, u - cs, u)
        max_eigenval = max(max_eigenval,abs(u)+cf)
    end
    return max_eigenval
end

function evaluate_energy(eq::MHD, γ, p, ρu, ρv, ρw, Bx, By, Bz, inv_rho)
    # p/(γ-1) + 0.5ρ|u|² + 0.5|B|²
    p / (γ-1) + 0.5 * (inv_rho * (ρu^2 + ρv^2 + ρw^2) + (Bx^2 + By^2 + Bz^2) * inv(4*pi))
end

function evaluate_pressure(eq::MHD, γ, E, ρu, ρv, ρw, Bx, By, Bz, inv_rho)
    # (γ - 1)(E - 0.5ρ|u|² - 0.5|B|²)
    (γ-1) * (E - 0.5 * (inv_rho * (ρu^2 + ρv^2 + ρw^2) + (Bx^2 + By^2 + Bz^2) * inv(4*pi)))
end
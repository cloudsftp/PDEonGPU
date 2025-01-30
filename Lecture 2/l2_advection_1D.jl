using Plots, Plots.Measures, Printf
default(size=(1200, 400), framestyle=:box, label=false, grid=false, margin=10mm, lw=6, labelfontsize=20, tickfontsize=20, titlefontsize=24)

@views function diffusion_1D()
    # physics
    lx   = 20.0
    vx   = 1.0

    # numerics
    nx   = 200
    nvis = 2

    # derived numerics
    dx   = lx / nx
    dt   = dx / abs(vx)
    nt   = nx^2 ÷ 100
    xc   = LinRange(dx / 2, lx - dx / 2, nx)

    # array initialisation
    C   = @. exp(-(xc - lx/4)^2)
    C_1 = copy(C)

    # time loop
    for it = 1:nt
        if vx > 0
            C[2:nx]   .-= dt.*vx.*diff(C)./dx
        else
            C[1:nx-1] .-= dt.*vx.*diff(C)./dx
        end

        if it % nvis == 0
            display(plot(xc, [C_1, C]; xlims=(0, lx), ylims=(-0.1, 1.1),
                         xlabel="lx", ylabel="Concentration",
                         title="time = $(round(it*dt,digits=1))"))
        end
    end
end

diffusion_1D()

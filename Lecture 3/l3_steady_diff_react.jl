using Plots, Plots.Measures, Printf
default(size=(1200, 400), framestyle=:box, label=false, grid=false, margin=10mm, lw=6, labelfontsize=20, tickfontsize=20, titlefontsize=24)

@views function diffusion_1D()
    # physics
    lx   = 20.0::Float64
    dc   = 1.0::Float64
    da   = 10.0::Float64
    re   = π + sqrt(π^2 + da)
    re   = 2π
    ρ    = (lx/(dc*re))^2
    C_eq = 0.1::Float64
    ξ    = lx^2/dc/da
    re   = 2π

    # numerics
    nx          = 200::Int64
    ϵtol        = 1e-8::Float64
    maxiter     = 20nx
    ncheck      = ceil(Int,0.25nx)

    # derived numerics
    dx   = lx / nx
    dt   = dx / √(1 / ρ)
    xc   = LinRange(dx / 2, lx - dx / 2, nx)

    # array initialisation
    C    = @. 1.0 + exp(-(xc-lx/4)^2) - xc/lx
    C_i  = copy(C)
    qx   = zeros(Float64, nx - 1)

    # iteration loop
    iter = 1; err = 2ϵtol; iter_evo = Float64[]; err_evo = Float64[]
    while err >= ϵtol && iter <= maxiter
        qx         .-= dt ./ (ρ * dc + dt) .* (qx + dc .* diff(C) ./ dx) # implicit
        #C[2:end-1] .-= dt  / dx .* diff(qx)
        C[2:end-1] .-= dt./(1 + dt/ξ) .*((C[2:end-1] .- C_eq)./ξ .+ diff(qx)./dx)

        if iter % ncheck == 0
            err = maximum(abs.(diff(dc.*diff(C)./dx)./dx .- ((C .- C_eq) ./ ξ)[2:nx-1]))
            push!(iter_evo,iter/nx); push!(err_evo,err)

            p1 = plot(xc,[C_i,C];xlims=(0,lx),ylims=(-0.1,2.0),
                    xlabel="lx",ylabel="Concentration",title="iter/nx=$(round(iter/nx,sigdigits=3))")

            p2 = plot(iter_evo,err_evo;xlabel="iter/nx",ylabel="err",
                    yscale=:log10,grid=true,markershape=:circle,markersize=10)

            display(plot(p1,p2;layout=(2,1)))
        end
        iter += 1
    end

    @show iter

    read(stdin, 1)
end

diffusion_1D()

using Plots, Plots.Measures, Printf
default(size=(1200, 400), framestyle=:box, label=false, grid=false, margin=10mm, lw=6, labelfontsize=20, tickfontsize=20, titlefontsize=24)

@views function diffusion_1D()
    # physics
    lx,ly = 20.0,20.0
    dc    = 1.0::Float64
    da    = 10.0::Float64
    re    = π + sqrt(π^2 + da)
    re    = 2π
    ρ     = (lx/(dc*re))^2
    C_eq  = 0.1::Float64
    ξ     = lx^2/dc/da
    re    = 2π

    # numerics
    nx,ny       = 100,100
    ϵtol        = 1e-8::Float64
    maxiter     = 20nx
    ncheck      = ceil(Int,0.25nx)

    # derived numerics
    dx,dy   = lx/nx,ly/ny
    xc,yc   = LinRange(dx/2,lx-dx/2,nx),LinRange(dy/2,ly-dy/2,ny)
    dt      = dx/sqrt(1/ρ)/sqrt(2)

    # array initialisation
    C       = @. 1.0 + exp(-(xc-lx/4)^2-(yc'-ly/4)^2) - xc/lx
    qx,qy   = zeros(nx-1,ny),zeros(nx,ny-1)

    # iteration loop
    iter = 1; err = 2ϵtol; iter_evo = Float64[]; err_evo = Float64[]
    while err >= ϵtol && iter <= maxiter

        qx                 .-= dt./(ρ + dt/dc).*(qx./dc .+ diff(C,dims=1)./dx)
        qy                 .-= dt./(ρ + dt/dc).*(qy./dc .+ diff(C,dims=2)./dy)
        C[2:end-1,2:end-1] .-= dt./(1 + dt/ξ) .*((C[2:end-1,2:end-1] .- C_eq)./ξ .+ diff(qx[:,2:end-1],dims=1)./dx .+
                                                                                    diff(qy[2:end-1,:],dims=2)./dy)

        if iter % ncheck == 0
            err = maximum(abs.(diff(dc.*diff(C, dims=1)./dx, dims=1)./dx .- ((C .- C_eq) ./ ξ)[2:nx-1]))
            push!(iter_evo,iter/nx); push!(err_evo,err)

            p1 = heatmap(xc,yc,C';xlims=(0,lx),ylims=(0,ly),clims=(0,1),aspect_ratio=1,
                        xlabel="lx",ylabel="ly",title="iter/nx=$(round(iter/nx,sigdigits=3))")

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

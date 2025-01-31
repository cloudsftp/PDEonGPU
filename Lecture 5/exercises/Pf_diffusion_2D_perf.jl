using Plots, Plots.Measures, Printf
default(size=(600, 500), framestyle=:box, label=false, grid=false, margin=10mm, lw=6, labelfontsize=11, tickfontsize=11, titlefontsize=11)

function Pf_diffusion_2D(; do_check=false)
    # physics
    lx, ly = 20.0, 20.0
    k_ηf   = 1.0

    # numerics
    nx, ny  = 500, 500
    ϵtol    = 1e-8
    maxiter = max(nx, ny)
    ncheck  = ceil(Int, 0.25max(nx, ny))
    cfl     = 1.0 / sqrt(2.1)
    re      = 2π

    # derived numerics
    dx, dy  = lx / nx, ly / ny
    xc, yc  = LinRange(dx / 2, lx - dx / 2, nx), LinRange(dy / 2, ly - dy / 2, ny)
    θ_dτ    = max(lx, ly) / re / cfl / min(dx, dy)
    β_dτ    = (re * k_ηf) / (cfl * min(dx, dy) * max(lx, ly))

    # array initialisation
    Pf      = @. exp(-(xc - lx / 2)^2 - (yc' - ly / 2)^2)
    qDx     = zeros(Float64, nx + 1, ny)
    qDy     = zeros(Float64, nx, ny + 1)
    r_Pf    = zeros(nx, ny)

    # precompute
    k_ηf_dx = k_ηf / dx
    k_ηf_dy = k_ηf / dy

    _1_θ_dτ = 1.0 / (1.0 + θ_dτ)
    _β_dτ = 1.0 / β_dτ

    _dx = 1.0 / dx
    _dy = 1.0 / dy
    # end precompute

    tic = time()
    niter = 0
    warmup = 11

    # iteration loop
    iter = 1; err_Pf = 2ϵtol
    while err_Pf >= ϵtol && iter <= maxiter
        qDx[2:end-1, :] .-= (qDx[2:end-1, :] .+ k_ηf_dx .* diff(Pf, dims=1)) .* _1_θ_dτ
        qDy[:, 2:end-1] .-= (qDy[:, 2:end-1] .+ k_ηf_dy .* diff(Pf, dims=2)) .* _1_θ_dτ
        Pf              .-= (diff(qDx, dims=1) .*_dx .+ diff(qDy, dims=2) .*_dy) .*_β_dτ

        if do_check && iter % ncheck == 0
            r_Pf .= diff(qDx, dims=1) ./ dx .+ diff(qDy, dims=2) ./ dy
            err_Pf = maximum(abs.(r_Pf))

            @printf("  iter/nx=%.1f, err_Pf=%1.3e\n", iter / nx, err_Pf)

            display(heatmap(xc, yc, Pf'; xlims=(xc[1], xc[end]), ylims=(yc[1], yc[end]), aspect_ratio=1, c=:turbo, clim=(0, 1)))
        end

        if iter == warmup
            tic = time()
            niter = 0
        end

        iter += 1
        niter += 1
    end

    toc = time() - tic
    tit = toc / niter
    A_eff = (3*2)/1e9*nx*ny*sizeof(Float64)
    T_eff = A_eff/tit

    @printf("Time = %1.3f sec, T_eff = %1.2f GB/s (niter = %d)\n", toc, round(T_eff, sigdigits=3), niter)

    return
end

Pf_diffusion_2D()

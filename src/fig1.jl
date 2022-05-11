function fig1()

    ################################################################################
    # Classical LL-type anisotropy
    ################################################################################
    N = 3
    D = -1.0
    θ = (4/5)*(π/2)
    ω = 2cos(θ) # angular frequency proportional to B only
    Δt = 0.5

    time = 10.
    nsteps = round(Int, time / Δt)

    sys = System(; spin_rescaling=1.0, N, D, L=1)
    set_state_from_dipole!(sys, 1, Vec3(sin(θ), 0, cos(θ)))

    t_hist = Δt * (0:nsteps)
    t2_hist = 0.1Δt * (0:10nsteps)

    sx_hist = Float64[]
    sy_hist = Float64[]
    e_hist = Float64[]

    for i=0:nsteps
        push!(sx_hist, sys.s[1][1])
        push!(sy_hist, sys.s[1][2])
        push!(e_hist, energy(sys))
        midpoint_step!(sys, Δt)
    end


    ## Plotting parameters
    color1 = 1
    color2 = 2
    plotparams = (;
        xlabelfontsize = 16,
        ylabelfontsize = 14,
        legendfontsize = 14,
        xtickfontsize = 12,
        ytickfontsize = 12,
        linewidth = 3.0,
        xlim=(0.0,10.0),
        ylim=(-1.0,1.0),
        xlabel=raw"Time",
        ylabel="Spin component",
        palette=:seaborn_colorblind,
        linealpha = 0.8,
    )
    plotparams_ref = (;
        label=false,
        linestyle=:dash,
        linecolor=:black,
        linealpha=0.85,
    )


    p1 = plot(;
        plotparams...,
        legend=:bottomright,
    )

    plot!(t_hist, sx_hist;
        plotparams...,
        label=raw"$s^x(t)$",
        color=color1,
        xticks = (0.0:2.0:10.0, ["" for _ = 1:6]))

    plot!(p1, t_hist, sy_hist;
        plotparams...,
        color=color2,
        label=raw"$s^y(t)$",
    )

    plot!(t2_hist, sin(θ) * cos.(ω*t2_hist);
        plotparams_ref...,
    )

    plot!(p1, t2_hist, -sin(θ) * sin.(ω*t2_hist);
        plotparams_ref...,
    )


    ################################################################################
    # Generalized spin dynamics (SU(3) anisotropy)
    ################################################################################
    Sz = spin_operators(N)[3]
    Λ = [-D*Sz^2]
    sys = System(; spin_rescaling=1.0, N, Λ, L=1)
    set_state_from_dipole!(sys, 1, Vec3(sin(θ), 0, cos(θ)))

    sx_hist_2 = Float64[]
    sy_hist_2 = Float64[]
    e_hist_2 = Float64[]

    for i = 0:nsteps
        push!(sx_hist_2, sys.s[1][1])
        push!(sy_hist_2, sys.s[1][2])
        push!(e_hist_2, energy(sys))
        midpoint_step!(sys, Δt)
    end
    t_hist = Δt * (0:nsteps)


    ## Plot results
    p2 = plot(;
        plotparams...,
        legend=:none,
    )


    plot!(p2, t_hist, sx_hist_2;
        plotparams...,
        label=raw"$s^x$",
        color=color1,
    )

    plot!(p2, t_hist, sy_hist_2;
        plotparams...,
        color=color2,
        label=raw"$s^y$",
    )

    plot!(p2, t2_hist, sin(θ) * cos.(D*t2_hist);
        plotparams_ref...,
    )

    plot!(p2, t2_hist, -cos(θ) * sin(θ) * sin.(D*t2_hist);
        plotparams_ref...,
    )

    ## Final layout
    p = plot(p1, p2, layout=(2, 1),
        title=["(a)" "(b)"],
        titleloc=:left,
    )

    return p
end

function fig4()

    ################################################################################
    # Generate data
    ################################################################################

    function collect_trajectory_gsd(; Δt)
        D = -1.0
        N = 3
        L = 100

        Sz = spin_operators(N)[3]
        Λ = fill(D*Sz^2, L)

        sys = System(; N, L, periodic=true, J=-1.0, Λ)

        for (i, x) in enumerate(range(0, 1, length=L))
            s = Vec3(cos(2π*x^2)*sin(2π*x^3),
                    sin(2π*x^2)*sin(2π*x^3),
                    cos(2π*x^3))
            set_state_from_dipole!(sys, i, s)
        end

        e0 = energy(sys)
        energies = Float64[]
        dipolemags = Float64[]

        t_max = 500.0
        nsteps = round(Int, t_max/Δt)
        for i = 1:nsteps
            push!(energies, abs(energy(sys) - e0))
            push!(dipolemags, norm(sys.s[60]))
            midpoint_step!(sys, Δt)
        end

        times = Δt * (0:nsteps-1)
        return (times, energies, dipolemags)
    end

    params = [
        (; Δt=0.08),
        (; Δt=0.04),
        (; Δt=0.02),
    ]
    data = map(params) do p
        collect_trajectory_gsd(; p...)
    end


    ################################################################################
    # Plot results
    ################################################################################
    colors = [1, 3, 2]
    linewidth=1.25
    plotparams = (;
        xlabelfontsize = 16,
        ylabelfontsize = 16,
        xtickfontsize = 12,
        ytickfontsize = 12,
        legendfontsize = 12,
        palette=:seaborn_colorblind,
    )

    plt1 = plot(;
        plotparams...,
        ylabel=L"$|\Delta H|$",
        legend=:topright,
    )

    plt2 = plot(;
        plotparams...,
        ylabel=L"$|\mathbf{s}_{60}|$",
        xlabel="Time",
        xrange=(0, 30),
        legend=false,
    )

    for (p, (ts, es, ds), c) = zip(params, data, colors)
        plot!(plt1, ts, es; linecolor=c, label="Δt=$(p.Δt)", linewidth)
        plot!(plt2, ts, ds; linecolor=c, label="Δt=$(p.Δt)", linewidth)
    end

    plt = plot(plt1, plt2;
        layout=grid(2, 1, heights=[0.5, 0.5]),
        margin=1.5mm,
        title=["(a)" "(b)"],
        titleloc=:left,
    )

    return plt
end

function fig2()

    ################################################################################
    # Generate data
    ################################################################################

    function collect_trajectory(; N, D, Δt, step!fn, label)
        L = 100
        spin_rescaling = N==2 ? 2.0 : 1.0

        sys = System(; N, spin_rescaling, periodic=true, L, J=-1.0, D)

        for (i, x) in enumerate(range(0, 1, length=L))
            s = Vec3(cos(2π*x^2)*sin(2π*x^3),
                    sin(2π*x^2)*sin(2π*x^3),
                    cos(2π*x^3))
            set_state_from_dipole!(sys, i, s)
        end

        e0 = energy(sys)
        energies = Float64[]
        mags = Float64[]

        t_max = 500.0
        nsteps = round(Int, t_max/Δt)
        for i = 1:nsteps
            push!(energies, abs(energy(sys) - e0))
            push!(mags, sum(n[1] for n=sys.s))
            step!fn(sys, Δt)
        end

        times = Δt * (0:nsteps-1)
        return (times, energies, mags)
    end


    params = [
        (; N = 2, D = 0., Δt = 0.1, step!fn = heun_step!, label="HeunP (non-symplectic)"),
        (; N = 3, D = 0., Δt = 0.1, step!fn = midpoint_step!, label="Schrodinger midpoint, spin-1"),
        (; N = 2, D = 0., Δt = 0.1, step!fn = midpoint_step!, label=raw"Schrodinger midpoint, spin-$\frac{1}{2}$"),
        (; N = 0, D = 0., Δt = 0.1, step!fn = spherical_midpoint_step!, label="Spherical midpoint")
    ]

    data = map(params) do p
        collect_trajectory(; p...)
    end


    ################################################################################
    # Plot results
    ################################################################################

    colors = [5, 1, 3, 2]
    linewidth = 1.25
    plotparams = (;
        palette=:seaborn_colorblind,
        xlabelfontsize = 16,
        ylabelfontsize = 16,
        xtickfontsize = 12,
        ytickfontsize = 12,
        linewidth = linewidth,
        legendfontsize = 11,
    )

    plt1 = plot(;
        plotparams...,
        ylabel=raw"$|\Delta H|$",
        yrange=(0, 5),
        legend=:topleft,
        xformatter=_->"",
    )

    yticks2 = [1e-5, 3e-5, 5e-5]
    plt2 = plot(;
        plotparams...,
        ylabel=raw"$|\Delta H|$",
        yrange=(0, 5.5e-5),
        yticks=yticks2,
        legend=false,
        xformatter=_->"",
    )

    plt3 = plot(;
        plotparams...,
        ylabel=raw"$|\Delta H|$",
        yrange=(0, 9e-7),
        xlabel=raw"Time",
        legend=false,
    )

    for (p, (ts, es), c) = zip(params, data, colors)
        plot!(plt1, ts, es; label=p.label, linecolor=c, linewidth)
        plot!(plt2, ts, es; linecolor=c, linewidth)
    end

    (ts, es) = data[end]
    plot!(plt3, ts, es; linecolor=colors[end], linewidth)

    p = plot(plt1, plt2, plt3, layout=grid(3, 1, heights=[0.45, 0.275, 0.275]))

    return p
end

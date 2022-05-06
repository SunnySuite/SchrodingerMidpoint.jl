function fig3()
    
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
        (; N = 2, D = -1., Δt = 0.02, step!fn = heun_step!, label="HeunP (non-symplectic)"),
        (; N = 3, D = -1., Δt = 0.02, step!fn = midpoint_step!, label="Schrodinger midpoint, spin-1"),
        (; N = 2, D = -1., Δt = 0.02, step!fn = midpoint_step!, label=raw"Schrodinger midpoint, spin-$\frac{1}{2}$"),
        (; N = 0, D = -1., Δt = 0.02, step!fn = spherical_midpoint_step!, label="Spherical midpoint")
    ]

    data = map(params) do p
        collect_trajectory(; p...)
    end


    ################################################################################
    # Plot results
    ################################################################################
    xlabelfontsize = 16
    ylabelfontsize = 16
    xtickfontsize = 12
    ytickfontsize = 12
    legendfontsize = 12
    linewidth = 1.1
    palette=:seaborn_colorblind
    colors = [5, 1, 3, 2]


    plt = plot(;
            xlabel=raw"Time",
            ylabel=raw"$|\Delta H|$",
            palette,
            xlabelfontsize,
            ylabelfontsize,
            xtickfontsize,
            ytickfontsize,
            legendfontsize,
            yrange=(0, 0.025),
            legend=:topright,
    )

    for (p, (ts, es, mags), c) = zip(params, data, colors)
        plot!(plt, ts, es; label=p.label, linecolor=c, linewidth)
    end

    # savefig(plt, "fig3.pdf")
    return plt

end

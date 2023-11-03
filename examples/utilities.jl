using Plots
using DiffEqDevTools: test_convergence
using PrettyTables: pretty_table

function convergence_tab_plot(prob, algs; dts = 0.5 .^ (1:10))
    for i in eachindex(algs)
        #convergence order
        sim = test_convergence(dts, prob, algs[i])
        err = sim.errors[:l∞]
        p = -log2.(err[2:end] ./ err[1:(end - 1)])
        #table
        #algname = string(Base.typename(typeof(algs[i])).wrapper)
        algname = string(algs[i])
        pretty_table([dts err [NaN; p]]; header = (["dt", "err", "p"]),
                     title = string("\n\n", string(algs[i])), title_alignment = :c,
                     title_autowrap = true, title_same_width_as_table = true)
        #plot
        pop!(sim.errors, :final)
        label = algname .* " " .* ["l∞" "l2"]
        if i == 1
            plot(sim; label = label)
        else
            plot!(sim; label = label)
        end
        plot!(legend = :bottomright)
    end
    display(current())
end

function myplot(sol, name = "", analytic = false)
    N = length(sol.u[1])
    if analytic == true
        plot(sol, color = palette(:default)[1:(2 * N)]', legend = :right,
             plot_analytic = true)
    else
        plot(sol, color = palette(:default)[1:N]', legend = :right, plot_analytic = false)
    end
    p = plot!(sol, color = palette(:default)[1:N]', denseplot = false,
              markershape = :circle, markerstrokecolor = palette(:default)[1:N]',
              linecolor = invisible(), label = "")
    title!(name)
    return p
end

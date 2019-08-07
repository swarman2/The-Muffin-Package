using Plots

p = plot([sin,cos],0,π)
plot!(xlims=(-π,2*π),ylims=(-2,2))
display(p)

x = plot([sin],0,π)
plot!(xlims=(-π,2*π),ylims=(-2,2))
display(x)
display(plot(p,x))

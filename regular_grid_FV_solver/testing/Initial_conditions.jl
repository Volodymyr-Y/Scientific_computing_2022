using CairoMakie


xs = range(0, 10, length = 25)
ys = range(0, 15, length = 25)
zs = [cos(x) * sin(y) for x in xs, y in ys]
println(size(zs))
heatmap(xs, ys, zs)
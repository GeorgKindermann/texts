using Plots
c1 = 0.75
c3 = 3
c2 = exp(-8)
c0 = 8

t = 0:.5:50
h = c0.*log1p.(c2.*t.^c3).^c1

#ih = c0.*c1.*c2.*c3.*t.^(c3.-1).*log1p.(c2.*t.^c3).^(c1.-1)./(c2.*t.^c3.+1)
#ih[1] = 0
ih = [0; diff(h) ./ diff(t)]
it = [0; t[2:end] .- diff(t) ./ 2.]
hi = [0; h[2:end] .- diff(h) ./ 2.]

mih = h ./ t
mih[1] = 0

plot(t, h, labels=nothing, xlabel="Alter [Jahre]", ylabel="Höhe [m]", color=:black, lw=3)
annotate!([t[end-5]], [h[end-5]+.8], Plots.text("Höhe", :right, rotation = 20))

axis2 = twinx()
plot!(axis2, it, ih, ylabel="Höhenzuwachs [m/Jahr]", labels=nothing, color=:gray, linestyle=:dot, lw=2, yguidefontcolor=:gray)
annotate!(axis2, [it[5]-1.3], [ih[5]], Plots.text("Laufend", :left, rotation = 73, color=:gray))
plot!(axis2, t, mih, labels=nothing, color=:gray, linestyle=:dash, lw=2)
annotate!(axis2, [t[end]], [mih[end]+.03], Plots.text("Durchschnitt", :right, rotation = -10, color=:gray))

axis3 = twiny()
plot!(axis3, ih, hi, labels=nothing, xlabel="Höhenzuwachs [m/Jahr]", linestyle=:dot, color=:blue, xguidefontcolor=:blue)
annotate!(axis3, [ih[15]], [hi[15]+.9], Plots.text("Laufend", :left, rotation = 9, color=:blue))
plot!(axis3, mih, hi, labels=nothing, xlabel="Höhenzuwachs [m/Jahr]", linestyle=:dash, color=:blue, xguidefontcolor=:blue)
annotate!(axis3, [mih[19]], [hi[19]+.9], Plots.text("Durchschnitt", :left, rotation = 30, color=:blue))

i = findmax(ih)[2]
vline!([it[i]], label=nothing, color=:gray)
hline!([hi[i]], label=nothing, color=:gray)
scatter!(axis2, [it[i]], [ih[i]], label=nothing, color=:gray)
scatter!([it[i]], [hi[i]], label=nothing, color=:white)
scatter!(axis3, [ih[i]], [hi[i]], label=nothing, color=:blue)
annotate!([it[i]], [hi[i]], Plots.text(" max. L.", :left, rotation = -45))

i = findmax(mih)[2]
plot!([0, 43], [0, (h[i]) / t[i] * 43], label=nothing, color=:black, linestyle=:dashdotdot)
vline!([t[i]], label=nothing, color=:gray)
hline!([h[i]], label=nothing, color=:gray)
scatter!(axis2, [t[i]], [mih[i]], label=nothing, color=:gray)
scatter!([t[i]], [h[i]], label=nothing, color=:white)
scatter!(axis3, [mih[i]], [h[i]], label=nothing, color=:blue)
annotate!([t[i]], [h[i]], Plots.text(" max. D.", :left, rotation = -45))

savefig("hoehenzuwachs.pdf")


using Optim
x = 0:.025:.6
f(x) = x - (x^1.8)
y = f.(x)
plot(x./x[end], y./y[end], xlab="Natürlicher Bestockungsgrad", ylab="Relativer Zuwachs", label="Mit Optimum", color=:black, lw=3)
hline!([.97], label=nothing, color=:gray, linestyle=:dash)
#i = findmax(y)[2]
#scatter!([x[i]/maximum(x)], [y[i]/y[end]], label=nothing, color=:white)
#xi = optimize(y->(f(y) - f(x[end]))^2, 0.1, 0.59).minimizer
#scatter!([xi/x[end]], [f(xi)/y[end]], label=nothing, color=:white)
xi = optimize(y->(f(y) - .97*f(x[end]))^2, 0.1, 0.59).minimizer
scatter!([xi/x[end]], [f(xi)/y[end]], label=nothing, color=:gray)
vline!([[xi/x[end]]], label=nothing, color=:gray, linestyle=:dash, z_order=:back)
xi = optimize(x->-f(x), 0.1, 0.59).minimizer
scatter!([xi/x[end]], [f(xi)/y[end]], label=nothing, color=:black)
xil = optimize(x->(f(x)-.97*f(xi))^2, 0.1, xi).minimizer
xir = optimize(x->(f(x)-.97*f(xi))^2, xi, 0.59).minimizer
plot!([xil/x[end], xir/x[end]], [f(xil)/y[end], f(xir)/y[end]], label=nothing, color=:gray, linesyle=:dash)
scatter!([xil/x[end]], [f(xil)/y[end]], label=nothing, color=:white)
scatter!([xir/x[end]], [f(xir)/y[end]], label=nothing, color=:white)

f(x) = x - (x^3)
y = f.(x)
plot!(x./maximum(x), y./y[end], label="Ohne Optimum", color=:black, lw=3, linestyle=:dot)
xi = optimize(y->(f(y) - .97*f(x[end]))^2, 0.1, 0.59).minimizer
scatter!([xi/x[end]], [f(xi)/y[end]], label=nothing, color=:gray)
vline!([[xi/x[end]]], label=nothing, color=:gray, linestyle=:dash, z_order=:back)

savefig("zuwachsBestockungsgrad.pdf")


x = 0:.025:.6
f(x) = x - (x^1.8)
y = f.(x)
plot(x./x[end], y./y[end], xlab="Natürlicher Bestockungsgrad", ylab="Zuwachs", label="Ausgangssituation", color=:black, lw=3, legend=:bottomright)
scatter!([x[12]/x[end]], [y[12]/y[end]], label=nothing, color=:white)
scatter!([x[20]/x[end]], [y[20]/y[end]], label=nothing, color=:gray)
f(x) = x - (x^1.4)
y0 = f.(x)
plot!(x[1:15]./x[end], 1.2 .* y0[1:15]./y0[end], label="Zukunft niedere Dichte", color=:black, lw=3, linestyle=:dash)
scatter!([x[15]/x[end]], [1.2 * y0[15]/y0[end]], label=nothing, color=:white)
f(x) = x - (x^2.5)
y1 = f.(x)
plot!(x[1:23]./x[end], 1.15 .* y1[1:23]./y1[end], label="Zukunft hohe Dichte", color=:black, lw=3, linestyle=:dot)
scatter!([x[23]/x[end]], [1.15 * y1[23]/y1[end]], label=nothing, color=:gray)
plot!([x[12]/x[end], x[15]/x[end]], [y[12]/y[end], 1.2*y0[15]/y0[end]], label=nothing, color=:gray, lw=2, arrow=true)
plot!([x[20]/x[end], x[23]/x[end]], [y[20]/y[end], 1.15*y1[23]/y1[end]], label=nothing, color=:gray, lw=2, arrow=true)
savefig("zuwachsveraenderungBestockungsgrad.pdf")


x = .6 ./ (1 .+ [0:1/6:1.5; 2])
f(x) = x - (x^1.8)
y = f.(x)
p1 = plot([0], [0], ylab="Relativer Zuwachs", label=nothing, legend=:bottomright, xformatter=_->"")
plot!([1], [1], label=nothing)
scatter!(x[2:end]./x[begin], y[2:end]./y[begin], label=nothing, color=:black)
f(x) = x - (x^3)
y = f.(x)
scatter!(x[2:end]./x[begin], y[2:end]./y[begin], label=nothing, color=:white)

x = [.6; .8 * .6 ./ (1 .+ [0:1/6:1.5; 2])]
p2 = plot([0], [0], xlab="Natürlicher Bestockungsgrad", ylab="Relativer Zuwachs", label=nothing, legend=:bottomright)
plot!([1], [1], label=nothing)
f(x) = x - (x^1.8)
y = f.(x)
scatter!(x[3:end]./x[begin], y[3:end]./y[begin], label=nothing, color=:black)
f(x) = x - (x^3)
y = f.(x)
scatter!(x[3:end]./x[begin], y[3:end]./y[begin], label=nothing, color=:white)
vline!([.8], label=nothing, color=:gray)

x = .6 ./ ((1.5 .+ (0:1/2:3/2)) ./ 1.5)
f(x) = x - (x^1.8)
y = f.(x)
p3 = plot([0], [0], label=nothing, xformatter=_->"", yformatter=_->"")
plot!([1], [1], label=nothing)
scatter!(x[2:end]./x[begin], y[2:end]./y[begin], label=nothing, color=:black)
f(x) = x - (x^3)
y = f.(x)
scatter!(x[2:end]./x[begin], y[2:end]./y[begin], label=nothing, color=:white)

x = [.6; .8 * .6 ./ ((1.5 .+ (0:1/2:3/2)) ./ 1.5)]
p4 = plot([0], [0], xlab="Natürlicher Bestockungsgrad", label=nothing, yformatter=_->"", legend=:bottomright)
plot!([1], [1], label=nothing)
f(x) = x - (x^1.8)
y = f.(x)
scatter!(x[3:end]./x[begin], y[3:end]./y[begin], label="Mit Optimum", color=:black)
f(x) = x - (x^3)
y = f.(x)
scatter!(x[3:end]./x[begin], y[3:end]./y[begin], label="Ohne Optimum", color=:white)
vline!([.8], label=nothing, color=:gray)

tt = stack([[collect(xlims(x)); collect(ylims(x))] for x = [p1, p2, p3, p4]])
tt = [minimum(tt[1,:]), maximum(tt[2,:]), minimum(tt[3,:]), maximum(tt[4,:])]
for x = [p1, p2, p3, p4]
  xlims!(x, tt[1], tt[2])
  ylims!(x, tt[3], tt[4])
end

plot(p1, p3, p2, p4, layout = (2, 2))
savefig("bestockungsgradstufen.pdf")


x = range(0, step = 0.01, stop = 2pi)
y1 = sin.(x)
y2 = cos.(x)

plt1 = plot(x, y1)
plt2 = plot(x, y2)
plt1 = plot(x, y1, xformatter=_->"")
plot(plt1, plt2, layout = (2, 1), framestyle=:box)



using VoronoiCells
using GeometryBasics
rect = Rectangle(Point2(0., 0.), Point2(5., 2*sqrt(3)))
points = [Point2(.5+x, y*sqrt(3)/2) for y in 0:4 for x in -y%2/2:5 if x < 5]
i = trues(length(points))
tess = voronoicells(points[i], rect)
p1 = plot(tess, label=nothing, aspect_ratio=:equal, color=:black, axis=([], :off), xlims=(0.7,4.3), ylims=(0.1,3.4), lw=3, size=(360, 340))
scatter!(points[i], label=nothing, color=:black)
scatter!(points[14], label=nothing, color=:black, markersize=10)

i[20] = false
p2 = plot(tess, label=nothing, aspect_ratio=:equal, color=:gray, axis=([], :off), xlims=(0.7,4.3), ylims=(0.1,3.4), size=(360, 340))
tess2 = voronoicells(points[i], rect)
plot!(tess2, label=nothing, xlims=(0.7,4.3), ylims=(0.1,3.4), color=:black, lw=3)
scatter!(points[i], label=nothing, color=:black)
scatter!(points[14], label=nothing, color=:black, markersize=10)

i[[2, 5, 6, 9, 13, 16, 17, 24, 27]] .= false
p3 = plot(tess, label=nothing, aspect_ratio=:equal, color=:gray, axis=([], :off), xlims=(0.7,4.3), ylims=(0.1,3.4), size=(360, 340))
tess2 = voronoicells(points[i], rect)
plot!(tess2, label=nothing, xlims=(0.7,4.3), ylims=(0.1,3.4), color=:black, lw=3)
scatter!(points[i], label=nothing, color=:black)
scatter!(points[14], label=nothing, color=:black, markersize=10)

i[[1, 4, 8, 11, 12, 15, 19, 22, 23, 26]] .= false
p4 = plot(tess, label=nothing, aspect_ratio=:equal, color=:gray, axis=([], :off), xlims=(0.7,4.3), ylims=(0.1,3.4), size=(360, 340))
tess2 = voronoicells(points[i], rect)
plot!(tess2, label=nothing, xlims=(0.7,4.3), ylims=(0.1,3.4), color=:black, lw=3)
scatter!(points[i], label=nothing, color=:black)
scatter!(points[14], label=nothing, color=:black, markersize=10)

plot(p1, p2, p3, p4, layout = (2, 2), size=(360*1.5, 340*1.5))
savefig("bestandesdichteEinzelbaum.pdf")

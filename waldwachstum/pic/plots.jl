using Plots

c1 = 0.75
c3 = 3
c2 = exp(-8)
c0 = 8
f(c0, t) = c0.*log1p.(c2.*t.^c3).^c1
#t = [x:2:x+20 for x in 0:10:140]
#h = [f(c0, x) .+ [0, 1, -1][x[1] ÷ 10 % 3 + 1] ./ (x[end] / 50)^.4 for x in t]
t = [x:2:x+20 for x in 0:16:140]
h = [f(c0, x) .+ [0, 1][x[1] ÷ 16 % 2 + 1] ./ ((20 .+ x) ./ 25).^.5 for x in t]
plot(t, h, label=nothing, xlabel="Alter [Jahre]", ylabel="Höhe [m]", color=:black, lw=2)
savefig("wuchsreiheHoehe.pdf")


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
t = 0:.5:150
f(c0, t) = c0.*log1p.(c2.*t.^c3).^c1
c1 = 0.75
c3 = 3.
c2 = exp(-6)
c0 = Optim.minimizer(optimize(c0 -> (30-f(c0[1], 100))^2, [8.], BFGS()))[1]
h = f(c0, t)
ih = [0; diff(h) ./ diff(t)]
it = [0; t[2:end] .- diff(t) ./ 2.]
plot(it, ih, ylabel="Höhenzuwachs [m/Jahr]", labels=nothing, color=:gray, lw=2, yguidefontcolor=:gray, linestyle=:dash, z_order=:front)
axis2 = twinx()
plot!(axis2, t, h, labels="Früh", xlabel="Alter [Jahre]", ylabel="Höhe [m]", color=:black, lw=2, linestyle=:dash, z_order=:front, legend=:right)
i = findmax(ih)[2]
scatter!([it[i]], [ih[i]], label=nothing, mc=:white, msc=:gray)
scatter!(axis2, [t[i]], [h[i]], label=nothing, mc=:white, msc=:black)
c2 = exp(-9)
c0 = Optim.minimizer(optimize(c0 -> (30-f(c0[1], 100))^2, [8.], BFGS()))[1]
h = f(c0, t)
ih = [0; diff(h) ./ diff(t)]
it = [0; t[2:end] .- diff(t) ./ 2.]
plot!(it, ih, labels=nothing, color=:gray, lw=2, z_order=:back)
plot!(axis2, t, h, labels="Mittel", color=:black, lw=2)
i = findmax(ih)[2]
scatter!([it[i]], [ih[i]], label=nothing, mc=:white, msc=:gray)
scatter!(axis2, [t[i]], [h[i]], label=nothing, mc=:white, msc=:black)
c2 = exp(-11.5)
c0 = Optim.minimizer(optimize(c0 -> (30-f(c0[1], 100))^2, [8.], BFGS()))[1]
h = f(c0, t)
ih = [0; diff(h) ./ diff(t)]
it = [0; t[2:end] .- diff(t) ./ 2.]
plot!(it, ih, labels=nothing, color=:gray, linestyle=:dot, lw=2, z_order=:back)
plot!(axis2, t, h, labels="Spät", color=:black, lw=2, linestyle=:dot)
i = findmax(ih)[2]
scatter!([it[i]], [ih[i]], label=nothing, mc=:white, msc=:gray)
scatter!(axis2, [t[i]], [h[i]], label=nothing, mc=:white, msc=:black)
savefig("hoehenzuwachsFrueMittelSpaet.pdf")


t = 0:.5:150
f(c0, t) = c0.*log1p.(c2.*t.^c3).^c1
c1 = 0.75
c3 = 3.
c2 = exp(-7)
c0 = 9
h = f(c0, t)
ih = [0; diff(h) ./ diff(t)]
it = [0; t[2:end] .- diff(t) ./ 2.]
plot(it, ih, ylabel="Höhenzuwachs [m/Jahr]", labels=nothing, color=:gray, lw=2, yguidefontcolor=:gray, linestyle=:dash, z_order=:front)
axis2 = twinx()
plot!(axis2, t, h, labels="Hoch", xlabel="Alter [Jahre]", ylabel="Höhe [m]", color=:black, lw=2, linestyle=:dash, z_order=:front, legend=(.8,.25))
i = findmax(ih)[2]
scatter!([it[i]], [ih[i]], label=nothing, mc=:white, msc=:gray)
scatter!(axis2, [t[i]], [h[i]], label=nothing, mc=:white, msc=:black)
c2 = exp(-8.4)
c0 = 7.6
h = f(c0, t)
ih = [0; diff(h) ./ diff(t)]
it = [0; t[2:end] .- diff(t) ./ 2.]
plot!(it, ih, labels=nothing, color=:gray, lw=2, z_order=:back)
plot!(axis2, t, h, labels="Mittel", color=:black, lw=2)
i = findmax(ih)[2]
scatter!([it[i]], [ih[i]], label=nothing, mc=:white, msc=:gray)
scatter!(axis2, [t[i]], [h[i]], label=nothing, mc=:white, msc=:black)
c2 = exp(-10)
c0 = 6
h = f(c0, t)
ih = [0; diff(h) ./ diff(t)]
it = [0; t[2:end] .- diff(t) ./ 2.]
plot!(it, ih, labels=nothing, color=:gray, linestyle=:dot, lw=2, z_order=:back)
plot!(axis2, t, h, labels="Nieder", color=:black, lw=2, linestyle=:dot)
i = findmax(ih)[2]
scatter!([it[i]], [ih[i]], label=nothing, mc=:white, msc=:gray)
scatter!(axis2, [t[i]], [h[i]], label=nothing, mc=:white, msc=:black)
savefig("hoehenzuwachsBonitaet.pdf")


using Optim
x = 0:.025:.6
f(x) = x - (x^1.8)
y = f.(x)
plot(x./x[end], y./y[end], xlab="Natürlicher Bestockungsgrad [1]", ylab="Relativer Zuwachs [1/ha/Jahr]", label="Mit Scheitelpunkt", color=:black, lw=3)
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
plot!(x./maximum(x), y./y[end], label="Ohne Scheitelpunkt", color=:black, lw=3, linestyle=:dot)
xi = optimize(y->(f(y) - .97*f(x[end]))^2, 0.1, 0.59).minimizer
scatter!([xi/x[end]], [f(xi)/y[end]], label=nothing, color=:gray)
vline!([[xi/x[end]]], label=nothing, color=:gray, linestyle=:dash, z_order=:back)

savefig("zuwachsBestockungsgrad.pdf")


pythonplot()  # arrow groese veraendern geht derzeit mit gr (noch) nicht
x = 0:.025:.6
f(x) = x - (x^1.8)
y = f.(x)
plot(x./x[end], y./y[end], xlab="Natürlicher Bestockungsgrad [1]", ylab="Zuwachs [m²/ha/Jahr]", label="Ausgangssituation", color=:black, lw=3, legend=:bottomright)
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
plot!([x[12]/x[end], x[15]/x[end]], [y[12]/y[end], 1.2*y0[15]/y0[end]], label=nothing, color=:gray, lw=2, arrow=(:closed, 1, 0.7))
plot!([x[20]/x[end], x[23]/x[end]], [y[20]/y[end], 1.15*y1[23]/y1[end]], label=nothing, color=:gray, lw=2, arrow=(:closed, 1, 0.7))
savefig("zuwachsveraenderungBestockungsgrad.pdf")
gr()


x = .6 ./ (1 .+ [0:1/6:1.5; 2])
f(x) = x - (x^1.8)
y = f.(x)
p1 = plot([0], [0], ylab="Relativer Zuwachs/ha", label=nothing, legend=:bottomright, xformatter=_->"")
plot!([1], [1], label=nothing)
scatter!(x[2:end]./x[begin], y[2:end]./y[begin], label=nothing, color=:black)
f(x) = x - (x^3)
y = f.(x)
scatter!(x[2:end]./x[begin], y[2:end]./y[begin], label=nothing, color=:white)

x = [.6; .8 * .6 ./ (1 .+ [0:1/6:1.5; 2])]
p2 = plot([0], [0], xlab="Natürlicher Bestockungsgrad", ylab="Relativer Zuwachs/ha", label=nothing, legend=:bottomright)
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
for (x, s) = zip([p1, p2, p3, p4], ["A", "C", "B", "D"])
  xlims!(x, tt[1], tt[2])
  ylims!(x, tt[3], tt[4])
  annotate!(x, tt[2], tt[4], (s, :right, :top))
end

plot(p1, p3, p2, p4, layout = (2, 2))
savefig("bestockungsgradstufen.pdf")


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

for (x, s) = zip([p1, p2, p3, p4], ["A", "B", "C", "D"])
#  annotate!(x, 4, 3, (s, :right, :top))
#  annotate!(x, 3.7, 3.4, (s, :right, :top))
  annotate!(x, 0.8, 3.5, (s, :right, :top))
end
plot(p1, p2, p3, p4, layout = (2, 2), size=(360*1.5, 340*1.5))
savefig("bestandesdichteEinzelbaum.pdf")


e0 = 1.6
n0 = 200000

d = [1, 100]
n = n0 .* d.^(-e0)
plot(d, n, axis=:log, xlims=extrema(d), ylims=extrema(n), label=nothing, color=:gray, linestyle=:dot, xlabel="Durchmesser [cm]", ylabel="Stammzahl [n/ha]", minorticks=9, lw=3)

#o = (2/3) .^ (0:16)
o = ((1/3)^.5) .^ (0:16)
#o = 1
#o = [o; vec(stack([global o = o[end] .* [2/3, 1/3] for i in 1:8]))]

for p in [[50000, 2.2, :solid, :black], [10000, 2, :dot, :black], [1500, 1.8, :solid, :gray]]
  n1 = p[1]
  e1 = p[2]
  print(n1)
  n = n1 .* o
  d = (n0 ./ n[1]).^(1/e0)
  n1 = d^e1 * n[1]
  d = (n1 ./ n).^(1/e1)
  d[2:2:end] .*= 1 .- 1 ./ (1 .+ n[2:2:end].^ .4)
  x = [1; repeat(d[1:end-1], inner=2)]
  y = [repeat(n[1:end-1], inner=2); n[end]]
  display(plot!(x, y, lw=2, color=p[4], label=nothing, linestyle=p[3]))
end
savefig("bestandesdichteDN.pdf")

using Interpolations
x = [550., 600., 500., 400] .* .8
y = [0.15, 0.8, 4., 10.]
jr = y[1]:.1:y[end]
f = interpolate(y, x, FiniteDifferenceMonotonicInterpolation())
plot(jr, f.(jr), label="Nadelholz", xlabel="Jahrringbreite [mm]", ylabel="Darrdichte ρ₀ [kg/m³]", legend=:bottomleft, lw=3, color=:black)
#x = [600., 680., 700.]
#y = [0.5, 4., 10.]
x = [550., 630., 680., 700.]
y = [0.5, 2., 4., 10.]
jr = y[1]:.1:y[end]
f = interpolate(y, x, FiniteDifferenceMonotonicInterpolation())
plot!(jr, f.(jr), label="Ringporig", lw=3, linestyle=:dot, color=:black)
x = [510., 530., 550, 550.]
y = [0.25, 0.5, 1., 10.]
jr = y[1]:.1:y[end]
f = interpolate(y, x, SteffenMonotonicInterpolation())
plot!(jr, f.(jr), label="Zerstreutporig", lw=3, linestyle=:dash, color=:black)
savefig("holzdichteJahrringbreite.pdf")


h = [ 1.3, 4.3, 8.3,12.3,16.3,19.3,21.3,23.3,25.3,27.3,29.3]
d = [1.06 0.97 0.91 0.98 1.05 1.07 1.19 1.18 1.30 1.36 2.02
     1.09 1.00 0.81 0.84 0.95 0.83 0.82 0.80 0.76 0.69 0.74
     1.49 1.28 1.15 1.09 1.05 0.95 0.94 0.86 0.81 0.70 0.64]
hi = h[1]:.25:h[end]
PF = plot(interpolate(h, d[1,:], SteffenMonotonicInterpolation()).(hi), hi, label=nothing, xlabel="Duchmesserzuwachs [cm/Jahr]", ylabel="Höhe [m]", legend=:outertop, xlim=(0,2.1), ylim=(0, 30), lw=2, color=:black)
plot!(interpolate(h, d[2,:], SteffenMonotonicInterpolation()).(hi), hi, label=nothing, lw=2, color=:black, linestyle=:dot)
plot!(interpolate(h, d[3,:], SteffenMonotonicInterpolation()).(hi), hi, label=nothing, lw=2, color=:black, linestyle=:dash)
plot!([], label="Jahrzehnt vor der Freistellung", color=:black)
plot!([], label="1. Jahrzehnt nach der Freistellung", color=:black, linestyle=:dot)
plot!([], label="2. Jahrzehnt nach der Freistellung", color=:black, linestyle=:dash)

h = [ 0.3, 1.3, 3.3, 5.3, 8.3,12.3,16.3,19.3,21.3]
d = [4.26 3.18 2.82 2.94 2.70 2.48 2.48 2.43 2.46
     3.77 3.10 2.76 2.95 2.74 2.58 2.38 2.32 2.18
     3.35 2.59 2.41 2.57 2.32 2.34 2.17 1.80 2.00] ./ 10
hi = h[1]:.25:h[end]
PS = plot(interpolate(h, d[1,:], SteffenMonotonicInterpolation()).(hi), hi, leg_title="Alter", label=nothing, xlabel="Duchmesserzuwachs [cm/Jahr]", ylabel="Höhe [m]", xlim=(0,maximum(d)), ylim=(0, maximum(h)), lw=2, color=:black)
plot!(interpolate(h, d[2,:], SteffenMonotonicInterpolation()).(hi), hi, label=nothing, lw=2, color=:black, linestyle=:dot)
plot!(interpolate(h, d[3,:], SteffenMonotonicInterpolation()).(hi), hi, label=nothing, lw=2, color=:black, linestyle=:dash)
plot!([], label="190-199", color=:black)
plot!([], label="200-209", color=:black, linestyle=:dot)
plot!([], label="209-219", color=:black, linestyle=:dash)

plot(PF, PS, layout = (1, 2))
savefig("freistellung.pdf")
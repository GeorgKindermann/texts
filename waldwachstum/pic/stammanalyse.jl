# Stammanalyse Guttenberg 1915
#h = [0.4, 1.3, 4.3, 8.3, 12.3, 16.3, 20.3, 23.3, 25.3, 27.3, 29.3, 31.3]
#d = [[3.31, 1.44],
#  [7.95, 6.96, 2.16],
#  [12.68, 11.54, 8.85, 3.50],
#  [17.32, 16.08, 14.63, 11.38, 6.22],
#  [21.95, 19.79, 18.68, 16.65, 13.16, 8.07, 1.07],
#  [27.08, 23.83, 22.65, 21.00, 18.57, 13.82, 6.99, 0.50],
#  [31.99, 27.14, 25.95, 24.42, 22.32, 18.49, 12.71, 7.11, 2.75],
#  [36.05, 30.36, 28.85, 27.20, 25.25, 22.10, 17.71, 13.10, 9.08, 5.08, 0.30],
#  [39.81, 33.28, 31.37, 29.71, 27.58, 24.80, 21.41, 17.52, 13.92, 10.00, 4.68],
#  [42.71, 35.81, 33.57, 31.75, 29.52, 26.92, 23.83, 20.38, 17.27, 13.90, 8.31, 3.24]]

using DelimitedFiles, Interpolations, Plots, Optim
SA = readdlm("kon_10_2.txt", '\t', Float32, '\n')
hs = [0.; SA[1,:]; 35.35]
SA[1,:] .= 0
t = 0:size(SA, 1)-1
i = findfirst(==(1.3f0), hs)-1
# d
d = SA[:,i]
PD = plot(t, d, labels=nothing, ylabel="BHD [cm]", color=:black, lw=3)
#annotate!(PD, [t[end-5]], [d[end-5]+1.5], Plots.text("Durchmesser", :right, rotation = 11))
di = [0; diff(d) ./ diff(t)]
ti = [0; t[2:end] .- diff(t) ./ 2.]
tp = twinx(PD)
plot!(tp, ti, di.*10, ylabel="Δd [mm/Jahr]", labels=nothing, color=:gray, linestyle=:dot, lw=2, yguidefontcolor=:gray, legend=(.8,.75))
#annotate!(tp, [ti[5]+2], [di[5]+.7], Plots.text("Laufend", :left, rotation = 85, color=:gray))
da = [0.; d[2:end]./t[2:end]]
plot!(tp, t, da.*10, labels=nothing, color=:gray, linestyle=:dash, lw=2)
#annotate!(tp, [t[end]], da[end]+0.04, Plots.text("Durchschnitt", :right, rotation = -12, color=:gray))
i = findmax(da)[2]
scatter!(tp, [t[i]], [da[i]*10], label=nothing, color=:gray)
vline!([t[i]], label=nothing, color=:gray)

plot!(PD, [], labels="Laufend", color=:gray, linestyle=:dot, legend=(.55,.7), fg_legend = :false, bg_legend = :false)
plot!(PD, [], labels="Durchschnitt", color=:gray, linestyle=:dash)
scatter!(PD, [], labels="Max ∅", color=:gray)

# g
g = d.^2 .*(pi/4.)
PG = plot(t, g./100, labels=nothing, xlabel="Alter [Jahre]", ylabel="Kreisfläche [dm²]", color=:black, lw=3)
#annotate!(PG, [t[end-5]], [g[end-5]+40], Plots.text("Kreisfläche", :right, rotation = 20))
gi = [0; diff(g) ./ diff(t)]
tp = twinx(PG)
plot!(tp, ti, gi, ylabel="Δg [cm²/Jahr]", labels=nothing, color=:gray, linestyle=:dot, lw=2, yguidefontcolor=:gray)
#annotate!(tp, [ti[5]+5], [gi[5]+7], Plots.text("Laufend", :left, rotation = 80, color=:gray))
ga = [0.; g[2:end]./t[2:end]]
plot!(tp, t, ga, labels=nothing, color=:gray, linestyle=:dash, lw=2)
#annotate!(tp, [t[end]], ga[end]+.8, Plots.text("Durchschnitt", :right, rotation = -10, color=:gray))
i = findmax([0.; g[2:end]./t[2:end]])[2]
scatter!(tp, [t[i]], [ga[i]], label=nothing, color=:gray)
vline!([t[i]], label=nothing, color=:gray)

# h
th = [0; vec(size(SA, 1)-0.5 .- sum(SA .> 0, dims=1)); 105]
f = interpolate(th, hs, SteffenMonotonicInterpolation())
h = f.(t)
hi = [0; diff(h) ./ diff(t)]
using Random
Random.seed!(42)
hi = hi .* (rand(length(h)) .* .4 .+ .8)
x = cumsum(hi)
i = 1 .+ ceil.(Int, th)
f = interpolate(i, h[i] - x[i], SteffenMonotonicInterpolation())
h = x .+ f.(1:length(x))
PH = plot(t, h, labels=nothing, ylabel="Höhe [m]", color=:black, lw=3)
#annotate!(PH, [t[end-5]], [h[end-5]+1], Plots.text("Höhe", :right, rotation = 20))
hi = [0; diff(h) ./ diff(t)]
tp = twinx(PH)
plot!(tp, ti, hi.*10, ylabel="Δh [dm/Jahr]", labels=nothing, color=:gray, linestyle=:dot, lw=2, yguidefontcolor=:gray)
#annotate!(tp, [ti[5]+5], [hi[5]+0.4], Plots.text("Laufend", :left, rotation = 80, color=:gray))
ha = [0.; h[2:end]./t[2:end]]
plot!(tp, t, ha.*10, labels=nothing, color=:gray, linestyle=:dash, lw=2)
#annotate!(tp, [t[end]], ha[end]+.05, Plots.text("Durchschnitt", :right, rotation = -10, color=:gray))
i = findmax(ha)[2]
scatter!(tp, [t[i]], [ha[i]*10], label=nothing, color=:gray)
vline!([t[i]], label=nothing, color=:gray)

# v
using Integrals
d0 = Vector{Float64}(undef, length(d))
for i in 1:length(d0)
  x = SA[i,:]
  x = x[x .> 0]
  if(length(x) > 1)
    d0[i] = extrapolate(interpolate(hs[2:length(x)+1], x, LinearMonotonicInterpolation()), Linear())(0)
  else
    d0[i] = -1.;
  end
end
i = findfirst(>(0), d0)
d0[1:i-1] .= interpolate([1; i:length(d0)], [0; d0[i:length(d0)]], SteffenMonotonicInterpolation()).(1:i-1)
v = Vector{Float64}(undef, length(d))
v[1] = 0.;
for i in 2:length(d0)
  x = SA[i,:]
  x = x[x .> 0]
  f = interpolate([hs[1:length(x)+1]; h[i]], [d0[i]; x; 0]./100., SteffenMonotonicInterpolation())
  v[i] = solve(IntegralProblem((x, p)->f(x)^2*pi/4, 0, h[i]), QuadGKJL()).u
end
PV = plot(t, v, labels=nothing, xlabel="Alter [Jahre]", ylabel="Volumen [m³]", color=:black, lw=3)
#annotate!(PV, [t[end-5]], [v[end-5]+1], Plots.text("Volumen", :right, rotation = 20))
vi = [0; diff(v) ./ diff(t)] .* 1000.
tp = twinx(PV)
plot!(tp, ti, vi, ylabel="Δv [dm³/Jahr]", labels=nothing, color=:gray, linestyle=:dot, lw=2, yguidefontcolor=:gray)
#annotate!(tp, [ti[5]+5], [vi[5]+0.4], Plots.text("Laufend", :left, rotation = 80, color=:gray))
va = [0.; v[2:end]./t[2:end]] .* 1000.
plot!(tp, t, va, labels=nothing, color=:gray, linestyle=:dash, lw=2)
#annotate!(tp, [t[end]], va[end]+.05, Plots.text("Durchschnitt", :right, rotation = -10, color=:gray))
i = findmax(va)[2]
scatter!(tp, [t[i]], [va[i]], label=nothing, color=:gray)
vline!([t[i]], label=nothing, color=:gray)

plot(PD, PH, PG, PV, layout = grid(2, 2), link=:x)
savefig("stammanalyseDGHV.pdf")


idx = 11:10:size(SA, 1)
S1 = [filter(>(0), SA[i,:]) for i in idx]

#f = [extrapolate(interpolate([0.; hs[2:length(S1[i])+1]; h[idx[i]]], [d0[idx[i]]; S1[i]; 0.], SteffenMonotonicInterpolation()), 0.) for i in 1:length(S1)]
f = [interpolate([0.; hs[2:length(S1[i])+1]; h[idx[i]]], [d0[idx[i]]; S1[i]; 0.], SteffenMonotonicInterpolation()) for i in 1:length(S1)]
f = [extrapolate(x, 0.) for x in f]
f = [x->0; f]

p1 = plot([-d0[101]; -S1[end]; 0; reverse(S1[end]); d0[101]], [hs[1:end-1]; h[101]; reverse(hs[1:end-1])], color=:black, fillrange=0, aspect_ratio=200, label=nothing, ylabel="Höhe [m]", xticks=false, xaxis=false, grid=false, yticks=([0; 1.3; 10:10:hs[end-1]]), lw=0.1)

p2 = plot(aspect_ratio=5, xlabel="d [cm]", showaxis=:x, left_margin=-12*Plots.mm, xticks=(-20:10:20, abs.(-20:10:20)), grid=false)
for i in 1:length(S1)
  x = [range(0, min(5, h[idx[i]]), 15); range(min(5, h[idx[i]]), h[idx[i]], 50)[2:end]]
  y = f[i+1](x)
  display(plot!([-y./2.; reverse(y./2.)], [x; reverse(x)], label=nothing, color=palette(:darkrainbow, length(S1))[i]))
end

p3 = plot(xlabel="   Δd [cm/10 Jahre]", showaxis=:x, yaxis=nothing, left_margin=-5Plots.mm)
for i in 1:length(S1)
  x = sort([range(0, h[idx[i]], 55); hs[2:i+1]])
  y = f[i+1](x) .- f[i](x)
  display(plot!(y, x, label=nothing, color=palette(:darkrainbow, length(S1))[i], lw=2))
end
xlims!(p3, -1, xlims(p3)[2])

p4 = plot(xlabel="Δg [cm²/10 Jahre]", showaxis=:x, yaxis=nothing, left_margin=-5Plots.mm)
axis2 = twinx()
plot!(axis2, yaxis="Jahrringe", yticks=([0; h[idx]], 100:-10:0), legend=:topright, leg_title="Alter")
for i in length(S1):-1:1
  x = range(0, h[idx[i]], 55)
  y = (f[i+1](x).^2 .- f[i](x).^2) .* (pi / 4.)
  display(plot!(axis2, y, x, label=i*10, color=palette(:darkrainbow, length(S1))[i], lw=2))
end
xlims!(p4, -20, xlims(p4)[2])

for p in [p1, p2, p3, p4]
  for i in 1:length(S1)
    hline!(p, [h[idx[i]]], color=palette(:darkrainbow, length(S1))[i], label=nothing, alpha=.2, z_order=:back)
  end
  hline!(p, [1.3], color=:gray, label=nothing, linestyle=:dash, z_order=:back)
  hline!(p, [0], color=:gray, label=nothing, z_order=:front)
end

for p in [p3, p4]
  vline!(p, [0], color=:gray, label=nothing, linestyle=:dash, z_order=:back)
end

plot(p1, p2, p3, p4, layout = grid(1, 4, widths=[0.07 ,0.23, .25, .45]), link=:y)
savefig("stammanalyse.pdf")


# Schaftform
using LambertW, Statistics, Interpolations, Optim, Plots, QuadGK, Integrals, LaTeXStrings

hm = hs
dm = [d0[end]; SA[end,:]; 0.0]
H  = hm[end]

f_messung = interpolate(hm, dm, SteffenMonotonicInterpolation())
hr = range(0, H, length=200)
dr_messung = f_messung.(hr)

ff(c0, c1, h, x) = c0 * (h - x)^c1
fob = function(c1, h, H, V)
    c3 = 2 * c1
    c2 = optimize(c2->(c2[1] * H^(c3+1) / (c3 + 1) - V)^2, [.001]).minimizer
    c0 = sqrt(c2[1] * 40000 / pi)
    (ff.(c0, c1, H, h), c0)
end

V = solve(IntegralProblem((x, p)->f_messung(x)^2*pi/40000., 0, hm[end]), QuadGKJL()).u

function schaftkurve_komplex(x, p)
    c0, c1, c2, c3, c4, h_wp = p
    if x >= h_wp
        return c0 * (H - x)^c1
    else
        return c2 + c3 * x^c4
    end
end
schaftkurve_komplex_vec(hr_vec, p) = [schaftkurve_komplex(x, p) for x in hr_vec]

function objective(p)
    c0, c1, c2, c3, c4, h_wp = p
    error = 0.0
    for (h_i, d_i) in zip(hr, dr_messung)
        d_pred = schaftkurve_komplex(h_i, p)
        error += (d_pred - d_i)^2
    end
    wert_oben = c0 * (H - h_wp)^c1
    wert_unten = c2 + c3 * h_wp^c4
    stetigkeit_penalty = 5000.0 * (wert_oben - wert_unten)^2
    
    abl_oben = -c0 * c1 * (H - h_wp)^(c1 - 1)
    abl_unten = (h_wp > 1e-3 && c4 > 0) ? (c3 * c4 * h_wp^(c4 - 1)) : 0.0
    knick_penalty = 5000.0 * (abl_oben - abl_unten)^2
    
    return error + stetigkeit_penalty + knick_penalty
end

lower = [0.001, 0.05, 0.0, -100.0, 0.01, 1.5]
upper = [100.0, 5.0, 200.0, 0.0, 10.0, 0.5 * H]
p_start = [1.0, 0.5, dm[1], -0.05, 2.0, 8.0]

res = optimize(objective, lower, upper, p_start, Fminbox(LBFGS()))
p_opt = res.minimizer
h_wp_opt = p_opt[6]
dr_komplex = schaftkurve_komplex_vec(hr, p_opt)

i = findfirst(==(1.3f0), hm)
G = (dm[i]/200)^2*pi
b = lambertw(G*(H-1.3)*log(1-1.3/H)/(V)) / log(1 - 1.3/H) - 1
gM = G / (1-1.3/H)^b
dr_bhd_ident = sqrt(gM/pi)*200/H^(b/2) .* (H .- hr) .^(b/2)

c1a = optimize(c1->(fob(c1[1], .3*H, H, V)[1] - f_messung(.3*H))^2, [1.5]).minimizer[1]
dr_d03_a = fob(c1a, hr, H, V)[1]
c1b = optimize(c1->(fob(c1[1], .3*H, H, V)[1] - f_messung(.3*H))^2, [.4]).minimizer[1]
dr_d03_b = fob(c1b, hr, H, V)[1]

c1c = optimize(c1->sum(sort((fob(c1[1], hr, H, V)[1] .- dr_messung).^2) .* range(1, 0, length(hr)).^.5), [0.5]).minimizer[1]
dr_nahe_v_ident = fob(c1c, hr, H, V)[1]

c1d = optimize(x->sum(sort((ff.(x[1], x[2], H, hr) .- dr_messung).^2) .* range(1, 0, length(hr)).^.5), [1., 1.]).minimizer
dr_nahe_messung = ff.(c1d[1], c1d[2], H, hr)

PS = plot(dr_messung, hr, label="Messung (Spline)", color=:lightgray, lw=5, linestyle=:solid,
          xlab="Durchmesser [cm]", ylab="Höhe [m]", legend=:topright, fontfamily="sans-serif")

plot!(dr_bhd_ident, hr, label="Durch BHD, V ident", lw=1.5, color=:red, linealpha=.7)
plot!(dr_d03_a, hr, label="Durch d03, V ident", lw=1.5, color=:green, linealpha=.7)
plot!(dr_d03_b, hr, label=nothing, lw=1.5, color=:green, linealpha=.7)
plot!(dr_nahe_v_ident, hr, label="Nahe Messung, V ident", lw=1.5, color=:blue, linealpha=.7)
plot!(dr_nahe_messung, hr, label="Nahe Messung", lw=1.5, color=:black, linealpha=.7)

plot!(dr_komplex, hr, label="Zweiteilige Funktion", color=:purple, lw=3, linestyle=:dot)

hline!([h_wp_opt], label=L"\mathsf{h}_{\mathsf{wp}}", color=:orange, linestyle=:dot, lw=1.5)

PD = plot([0, 0], hr[[1, end]], color=:lightgray, lw=5, linestyle=:solid, 
          xlab="Messung - Kurve [cm]", label=nothing, fontfamily="sans-serif")

plot!(dr_messung .- dr_bhd_ident, hr, label=nothing, lw=1.5, color=:red, linealpha=.7)
plot!(dr_messung .- dr_d03_a, hr, label=nothing, lw=1.5, color=:green, linealpha=.4)
plot!(dr_messung .- dr_d03_b, hr, label=nothing, lw=1.5, color=:green, linealpha=.4)
plot!(dr_messung .- dr_nahe_v_ident, hr, label=nothing, lw=1.5, color=:blue, linealpha=.7)
plot!(dr_messung .- dr_nahe_messung, hr, label=nothing, lw=1.5, color=:black, linealpha=.7)

plot!(dr_messung .- dr_komplex, hr, label=nothing, lw=3, color=:purple, linestyle=:dot)

hline!([h_wp_opt], label=nothing, color=:orange, linestyle=:dot, lw=1.5)
xlims!((-2.5, 2.5))

plot(PS, PD, layout = grid(1, 2), link=:y, guidefontsize=10, tickfontsize=9, legendfontsize=7, legend=:topright)
savefig("schaftkurven_vergleich.pdf")



# Verschiedene Durchmesser
dr = Vector{Float64}(undef, length(d))
for i in 1:length(d)
    H_total = h[i]
    if H_total >= 1.3
        dt = 0.8
        dr[i] = d[i] + 1.3
    else
        dt = 0.1 + 0.7 * (1.0 - ((1.3 - H_total) / 1.3)^2)
        dr[i] = dt + H_total
    end
end
df = sqrt.(v ./ h .* (4/pi)) .* 100
df[1] = 0
hdf = Vector{Float64}(undef, length(d))
hdf[1] = 0
d03 = Vector{Float64}(undef, length(d))
d03[1] = 0
for i in 2:length(d03)
  x = SA[i,:]
  x = x[x .> 0]
  f = interpolate([hs[1:length(x)+1]; h[i]], [d0[i]; x; 0], SteffenMonotonicInterpolation())
  hdf[i] = optimize(x->(df[i]-f(x))^2, 0, h[i]).minimizer
  d03[i] = f(.3 * h[i])
end
using LaTeXStrings
PD = plot(t, d, label="BHD", xlabel="", ylabel="Durchmesser [cm]", lw=2, color=:black, legend=(.5, .32), legend_columns=4, fg_legend = :false, background_color_legend = nothing)
plot!(t, d03, label=L"\mathrm{d_{03}}", lw=2, color=:red)
plot!(t, df, label=L"\mathrm{d_{f}}", lw=2, color=:blue)
#plot!(t, d00e, label=L"\mathrm{d_{00e}}", lw=2, color=:green)
plot!(t, dr, label=L"\mathrm{d_{r}}", lw=2, color=:green)
tp = twinx()
plot!(tp, t, 1.3 ./ h, label=nothing, ylabel="Messhöhe d / Höhe", lw=2, color=:black, linestyle=:dash)
plot!(tp, t[[1, end]], [.3, .3], label=nothing, lw=2, color=:red, linestyle=:dash)
plot!(tp, t, hdf ./ h, label=nothing, lw=2, color=:blue, linestyle=:dash)
plot!(tp, t[[1, end]], [0, 0], label=nothing, lw=2, color=:green, linestyle=:dash)
ylims!(tp, (-0.015, .55))
plot!(tp, [], label="Durchmesser", lw=2, color=:gray, legend=(.5, .25), legend_columns=2, fg_legend = :false, background_color_legend = nothing)
plot!(tp, [], label="Messhöhe", lw=2, color=:gray, linestyle=:dash)

# h/d-Wert
PHD = plot(t, h ./ d .* 100, label="BHD", xlabel="Alter [Jahre]", ylabel="h/d-Wert [1]", lw=2, color=:black, legend=(.35, .12), legend_columns=4, fg_legend = :false, background_color_legend = nothing)
plot!(t, h ./ d03 .* 100, label=L"\mathrm{d_{03}}", lw=2, color=:red)
plot!(t, h ./ df .* 100, label=L"\mathrm{d_{f}}", lw=2, color=:blue)
#plot!(t, h ./ d00e .* 100, label=L"\mathrm{d_{00e}}", lw=2, color=:green)
plot!(t, h ./ dr .* 100, label=L"\mathrm{d_{r}}", lw=2, color=:green)
ylims!((ylims()[1], 150))
# Formzahl
tp = twinx()
plot!(tp, t, v ./ (d.^2 .* pi ./ 40000 .* h), label=nothing, ylabel="Formzahl", lw=2, color=:black, linestyle=:dash)
plot!(tp, t, v ./ (d03.^2 .* pi ./ 40000 .* h), label=nothing, lw=2, color=:red, linestyle=:dash)
plot!(tp, t[[1, end]], [1, 1], label=nothing, lw=2, color=:blue, linestyle=:dash)
#plot!(tp, t, v ./ (d00e.^2 .* pi ./ 40000 .* h), label=nothing, lw=2, color=:green, linestyle=:dash)
plot!(tp, t, v ./ (dr.^2 .* pi ./ 40000 .* h), label=nothing, lw=2, color=:green, linestyle=:dash)
ylims!(tp, (0.15, 1.15))
plot!(tp, [], label="h/d-Wert", lw=2, color=:gray, legend=(.35, .05), legend_columns=2, fg_legend = :false, background_color_legend = nothing)
plot!(tp, [], label="Formzahl", lw=2, color=:gray, linestyle=:dash)

plot!(PD, xticks=(xticks(PHD)[1][1], []))

plot(PD, PHD, layout = grid(2, 1), link=:x)
savefig("stammanalyseDfz.pdf")


#Entwicklung der Stammform
using LambertW, Statistics, Interpolations, Optim, Plots, QuadGK, Integrals, LaTeXStrings
using Plots.Measures

b0 = Vector{Float64}(undef, length(d)); b0 .= NaN
b1 = Vector{Float64}(undef, length(d)); b1 .= NaN
b2 = Vector{Float64}(undef, length(d)); b2 .= NaN

h_wp_vec = Vector{Float64}(undef, length(d)); h_wp_vec .= NaN
d_wp_vec = Vector{Float64}(undef, length(d)); d_wp_vec .= NaN
c1_upper = Vector{Float64}(undef, length(d)); c1_upper .= NaN
c4_lower = Vector{Float64}(undef, length(d)); c4_lower .= NaN

last_stable_rel_wp = 0.15 

for i in 2:length(d)
    x_data = SA[i,:]
    x_data = x_data[x_data .> 0]
    hm = [hs[1:length(x_data)+1]; h[i]]
    dm = [d0[i]; x_data; 0.0]
    H  = hm[end]
    G  = (d[i]/200)^2*pi
    V  = v[i]

    f = interpolate(hm, dm, SteffenMonotonicInterpolation())
    hr = range(0, H, length=100)
    dr = f.(hr)

    if (G > .001 && H > 1.3 && V > .01)
        b0[i] = (lambertw(G*(H-1.3)*log(1-1.3/H)/(V)) / log(1 - 1.3/H) - 1) / 2.
    else
        b0[i] = NaN
    end

    b1[i] = optimize(c1->sum(sort((fob(c1[1], hr, H, V)[1] .- dr).^2) .* range(1, 0, length(hr)).^.5), [0.7]).minimizer[1]
    b2[i] = optimize(x->sum(sort((ff.(x[1], x[2], H, hr) .- dr).^2) .* range(1, 0, length(hr)).^.5), [1., .7]).minimizer[2]
    
    if H > 3.0
        X = [ones(length(hm)) hm hm.^2 hm.^3]
        beta = X \ dm  
        c_koeff = beta[3]
        d_koeff = beta[4]
        
        h_wp_analytisch = -c_koeff / (3.0 * d_koeff)
        rel_wp_aktuell = h_wp_analytisch / H
        
        if 0.10 <= rel_wp_aktuell <= 0.30 && abs(rel_wp_aktuell - last_stable_rel_wp) < 0.05
            h_wp_opt = h_wp_analytisch
            last_stable_rel_wp = rel_wp_aktuell 
        else
            h_wp_opt = last_stable_rel_wp * H
        end
        
        d_wp_opt = f(h_wp_opt)
        
        mask_up = (hr .>= h_wp_opt) .& (hr .< H)
        hr_up   = hr[mask_up]
        dr_up   = dr[mask_up]
        
        loss_upper(c1) = sum((( (d_wp_opt / (H - h_wp_opt)^c1) .* (H .- hr_up).^c1 ) .- dr_up).^2)
        res_up = optimize(loss_upper, 0.1, 3.0)
        c1_o = res_up.minimizer
        
        d0_eff = d0[i]
        mask_low = hr .< h_wp_opt
        hr_low   = hr[mask_low]
        dr_low   = dr[mask_low]
        
        loss_lower(c4) = sum(((d0_eff .+ ((d_wp_opt - d0_eff) / (h_wp_opt^c4)) .* hr_low.^c4) .- dr_low).^2)
        res_low = optimize(loss_lower, 0.5, 5.0)
        c4_o = res_low.minimizer
        
        h_wp_vec[i] = h_wp_opt
        d_wp_vec[i] = d_wp_opt
        c1_upper[i] = c1_o
        c4_lower[i] = c4_o
    end
end

PS = plot(t, b0, label="Durch BHD, V ident", xlabel="Alter [Jahre]", 
          ylabel="Formfaktor " * L"  \mathrm{c_1 : d_x = c_0 * (h - x)^{c_1}}", 
          lw=2, color=:red, alpha=.6, legend=:topright, legendfontsize=6.0, fontfamily="sans-serif")
plot!(t, b1, label="Nahe Messung, V ident", lw=2, color=:blue, alpha=.6)
plot!(t, b2, label="Nahe Messung", lw=2, color=:black, alpha=.6)

plot!(t, c1_upper, label="Zweiteilige Schaftfunktion (oberer Teil)", lw=3, color=:purple, linestyle=:dot)

hr_geo = 0:.01:1
dr_geo = ff.(1., 2.0, hr_geo[end], hr_geo)
P2 = plot([-dr_geo; reverse(dr_geo)], [hr_geo; reverse(hr_geo)], lw=2, label="c1=2.0", 
          legend=:bottom, legendfontsize=6, xlabel="Durchmesser", ylabel="\nHöhe", ticks=false, fontfamily="sans-serif")
dr_geo = ff.(1., 1., hr_geo[end], hr_geo)
plot!([-dr_geo; reverse(dr_geo)], [hr_geo; reverse(hr_geo)], lw=2, label="c1=1.0")
dr_geo = ff.(1., .5, hr_geo[end], hr_geo)
plot!([-dr_geo; reverse(dr_geo)], [hr_geo; reverse(hr_geo)], lw=2, label="c1=0.5")
hline!([0], label=nothing, color=:gray, lw=2)

plot(PS, P2, layout = grid(1, 2))
savefig("stammanalyse_kompakt.pdf")

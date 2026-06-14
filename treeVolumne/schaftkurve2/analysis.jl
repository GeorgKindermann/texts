#download("https://www.envidat.ch/dataset/b23a6bfa-3a55-4677-91e2-d32f91945574/resource/8ba7cb02-3d90-409f-9113-4e9aa5e7d33c/download/stembranchchwood_2024-04-25.txt", "stembranchchwood_2024-04-25.txt")
# https://doi.org/10.16904/envidat.486
# http://dx.doi.org/10.1038/s41597-024-03336-7

using DataFrames, CSV
df = CSV.read("../stembranchchwood_2024-04-25.txt", DataFrame, missingstring="NA")
filter!(:TreeSpecies => sp -> !ismissing(sp) && sp == 21, df)  # Fichte

using StatsBase, Plots
describe(df[:, [:H_total, :DBH]])

Plots.scalefontsizes(1.4)
scatter(df.DBH ./ 10.0, df.H_total ./ 10, legend=false, xlabel="DBH [cm]", ylabel="Height [m]", color=:black, ylim=[-1, 50], ms=1, ma=0.2)
savefig("dhScatter.png")
Plots.scalefontsizes()

Plots.scalefontsizes(1.4)
i = round.(Int, range(1, size(df, 1), 100))
plot(sort(df.DBH ./ 10.0)[i], i, label="DBH", xlabel="DBH [cm]", ylabel="Observations", lw=2, formatter=:plain, color=:black, xlim=[0, 100], widen=true, xticks=0:20:100)
plot!([0], [0], label="Height", lw=2, color=:red)
plot!(twiny(), sort(df.H_total ./ 10.0)[i], i, label=nothing, xlabel="Height [m]", lw=2, color=:red, xlim=[0, 50], widen=true, xticks=0:10:50, grid=true, foreground_color_guide=:red, foreground_color_text=:red, foreground_color_axis=:red, foreground_color_border=:red)
savefig("dhCumObs.pdf")
Plots.scalefontsizes()

h = [-1.0; 0.0; 0.65; 1; 1.3; 3:2:53]  # h-1, h0, h065, h1, h1.3, h3, h5, ..., h51, h53
ht = df[:, [:H_total, :L_coarsestem, :L_coarsestemfinal, :L_top]]
ht.hc = ht.L_coarsestem .- ht.L_coarsestemfinal ./ 2.0
ht.ht = ht.L_coarsestem .+ ht.L_top ./ 2.0
ht = ht[:, [:hc, :L_coarsestem, :ht, :H_total]] ./ 10.0  # hd>7, hd=7, hd<7, h

d = df[:, Cols(:DBH, r"^DM", :D_coarsestemfinal, :D_top)] ./ 10.0
replace!(d.D_coarsestemfinal, 0.0 => missing)
replace!(d.D_top, 0.0 => missing)
insertcols!(d, :D_top, :d7 => 7.0)
d.dw .= 0.8 # Annahme
d.d0 = Vector{Union{Missing,Float64}}(missing, nrow(d))
d.du1 = Vector{Union{Missing,Float64}}(missing, nrow(d))
select!(d, :du1, :d0, :DM065, :DM1, :) # d-1, d0, d065, d1, d1.3, d3, d5, ..., d51, d53, d>7, d=7, d<7, dw

#d bei h=0
using LsqFit
@. model(h, p) = p[1] + p[2] * h^p[3]
@. model1(h, p) = p[1] + p[2] * h
x = Matrix(coalesce.(d[:, [:DM065, :DM1, :DBH, :DM3, :DM5]], NaN))
y = 5.0 .- [0.65, 1.0, 1.3, 3.0, 5.0]
d0 = map(eachrow(x)) do x
  if (all(.!isnan.(x)) && length(unique(x)) > 2 && issorted(x, rev=true))
    fit = curve_fit(model1, y, x, [x[end], 1.0])
    fit = curve_fit(model, y, x, [coef(fit); 1.0])
    model([6.0, 5.0, 3.70], coef(fit))
  else
    [missing, missing, missing]
  end
end |> stack
sum(.!ismissing.(d0), dims=2)  # 2828
describe(d.DBH ./ d0[3, :])
#i = findall(coalesce.(d0[1, :] ./ d.DBH .> 3.0, false))  # 412
i = findall(coalesce.(d0[1, :] ./ d0[3, :] .> 3.0, false))  # 410
d0[1, i] .= missing
i = findall(coalesce.(d0[2, :] ./ d0[3, :] .> 2.0, false))  # 47
d0[2, i] .= missing
d00 = copy(d0)


#Hilfsfunktion um Residuen zu plotten
using Loess
function rp(y, x, sp=0.1; h=0., xlab="", ylab="", ms=2, ma=0.2, col=:red, ax=:xy, trim=0.02, trx=:identity)
  scatter(x, y, color=:black, ms=ms, ma=ma, xlabel=xlab, ylabel=ylab, showaxis=ax, xscale=trx)
  i = @. !ismissing(x) && isfinite(x) && !ismissing(y) && isfinite(y)
  xi, yi = Float64.(x[i]), Float64.(y[i])
  le = loess(xi, yi; span=sp)
  x_min, x_max = quantile(filter(!isnan, x), [trim, 1.0 - trim])
  us = range(x_min, x_max; length=200) 
  vs = predict(le, us)
  hline!([h], lw=2, color=:grey, linestyle=:dash)
  plot!(us, vs, legend=false, lw=3, color=col)
end

#Hilfsfunktion für LsqFit Ergebnisse
using DataFrames, Distributions
function quick_summary(fit)
    cc = coef(fit)
    se = stderror(fit)
    tv = cc ./ se
    # p-Werte basierend auf der t-Verteilung
    pv = 2 * (1 .- cdf(TDist(dof(fit)), abs.(tv)))
    ci = confint(fit)
    DataFrame(
        Parameter = ["p[$i]" for i in 1:length(cc)],
        Estimate = cc, StdError = se,
        t_stat = tv, p_val = pv,
        Lower95 = [c[1] for c in ci], Upper95 = [c[2] for c in ci]
    )
end

# Duchmesser bei h=0 aus dAux (bhd + 1.3/2) und h schaetzen
using GLM
using LaTeXStrings
#D = DataFrame(dbh=d.DBH, d0=d00[2, :], h=ht.H_total, dAux = d.DBH .+ 1.3 ./ 2)
D = DataFrame(dbh=d00[3, :], d0=d00[2, :], h=ht.H_total, dAux = d.DBH .+ 1.3 ./ 2, dD = d00[3, :])
dropmissing!(D)
size(D)  # 2781
D.s = (D.d0 .- D.dbh) ./ 130 .- .015
rp(D.s, D.h, col=:red)
rp(D.s, D.dAux, col=:red)
rp(D.s, D.dD, col=:red)
sum(D.s .< 0)  # 223
sum(D.s .< 0.001)  # 238
D.s[D.s .< 0.001] .= 0.001
reg = glm(@formula(s ~ inv(h) + inv(dD)), D, Gamma(), LogLink())
#                  Coef.  Std. Error       z  Pr(>|z|)   Lower 95%   Upper 95%
#(Intercept)   -0.300035   0.0433191   -6.93    <1e-11   -0.384939   -0.215132
#inv(h)       -39.0406     2.14614    -18.19    <1e-73  -43.247     -34.8343
#inv(dD)      -31.4588     2.21926    -14.18    <1e-44  -35.8085    -27.1092
i = fill(true, size(d)[1])
x = predict(reg, DataFrame(dD=d[i, :DBH], h=ht[i, :H_total])) .+ .015
d0[2, i] .= d[i, :DBH] .+ x .* 130
Plots.scalefontsizes(1.4)
rp(response(reg) .- predict(reg), predict(reg), col=:red, ylab=L"s_{0} - \hat{s}_{0}", xlab=L"\hat{s}_{0}")
savefig("d0ResidPred.png")
rp(response(reg) .- predict(reg), D.h, col=:red, ylab=L"s_{0} - \hat{s}_{0}", xlab="Height [m]")
savefig("d0ResidHeight.png")
rp(response(reg) .- predict(reg), D.dAux, col=:red, ylab=L"s_{0} - \hat{s}_{0}", xlab=L"d_a [cm]")
savefig("d0ResidDAux.png")
Plots.scalefontsizes()
predict(reg, DataFrame(h=[0, 1.3, 30, 70], dD=[0, 1.45, 45, 200])) .+ .015
plot((predict(reg, DataFrame(h=90, dD=0:200)) .+ .015) .* 130)
# d bei h-1
D = DataFrame(dbh=d00[3, :], d0=d00[2, :], dm1=d00[1, :], h=ht.H_total, dAux = d00[3, :] .+ 1.3 ./ 2, dD=d00[3, :])
dropmissing!(D)
size(D)  # 2418
D.s = (D.dm1 .- D.dbh) ./ 230 .- (D.d0 .- D.dbh) ./ 130
rp(D.s, D.h, col=:red)
rp(D.s, D.dD, col=:red)
sum(D.s .< 0)  # 86
sum(D.s .< 0.001)  # 250
D.s[D.s .< 0.001] .= 0.001
reg = glm(@formula(s ~ inv(h) + inv(dD)), D, Gamma(), LogLink())
reg = glm(@formula(s ~ inv(h)), D, Gamma(), LogLink())
#                 Coef.  Std. Error       z  Pr(>|z|)  Lower 95%   Upper 95%
#(Intercept)   -1.04765   0.0586381  -17.87    <1e-70   -1.16258   -0.932721
#inv(h)       -59.9838    1.43177    -41.89    <1e-99  -62.79     -57.1775
i = fill(true, size(d)[1])
x = predict(reg, DataFrame(h=ht[i, :H_total]))
d0[1, i] .= d[i, :DBH] .+ ((d0[2, i] .- d[i, :DBH]) ./ 130 .+ x) .* 200
Plots.scalefontsizes(1.4)
rp(response(reg) .- predict(reg), predict(reg), col=:red, ylab=L"ds_{-1m} - \hat{ds}_{-1m}", xlab=L"\hat{ds}_{-1m}")
savefig("dN1ResidPred.png")
rp(response(reg) .- predict(reg), D.h, col=:red, ylab=L"ds_{-1m} - \hat{ds}_{-1m}", xlab="Height [m]")
savefig("dN1ResidHeight.png")
rp(response(reg) .- predict(reg), D.dAux, col=:red, ylab=L"ds_{-1m} - \hat{ds}_{-1m}", xlab=L"d_a [cm]")
savefig("dN1ResidDAux.png")
Plots.scalefontsizes()
predict(reg, DataFrame(h=[0, 1.3, 30, 70], dAux=[0, 1.45, 45, 200])) .* 200
plot(predict(reg, DataFrame(h=0:90, dAux=500)) .* 200)
#
d0 = round.(d0, digits=1)
d00 = round.(d00, digits=1)
d.du1 = d0[1, :]
i = .!ismissing.(d00[1, :])  # optional
d.du1[i] .= d00[1, i]
d.d0 = d0[2, :]
i = .!ismissing.(d00[2, :])
d.d0[i] .= d00[2, i]
d.d0 = vec(maximum(coalesce.(Matrix(d[:, 2:end]), 0.0), dims=2))
d.du1 = vec(maximum(coalesce.(Matrix(d), 0.0), dims=2))
#
D = d0 = d00 = x = i = y = y2 = reg = nothing

# Anstiege bestimmen
x = coalesce.(collect(Matrix(d[:, 3:end])'), NaN)
y = collect(Matrix(ht[:, 1:end])')
y = vcat(repeat(h[3:end], 1, size(y, 2)), y)
y[isnan.(x)] .= NaN

Plots.scalefontsizes(1.4)
scatter(vec((x[1:end-1,:]) ./ (1.3 .+ x[2,:])'), vec(y[1:end-1,:] ./ y[end,:]'), xlabel=L"d / d_r", ylabel=L"h / h_{max}", color=:black, ms=1, ma=.1, legend=false, fmt = :png)
savefig("relativeTaper.png")
Plots.scalefontsizes()

res = map(eachcol(x), eachcol(y)) do xd, yh
  i = findall(.!isnan.(xd))
  D = combine(groupby(DataFrame(h=yh[i], d=xd[i]), :h, sort=true), :d => mean => :d)
  dh = diff(D.h)
  s = diff(D.d) ./ 50.0 ./ dh  # Neigung bezogen auf Markroehre
  rhs = (D.h[1:end-1] + dh ./ 2) ./ D.h[end]
#  ds = diff(s)
  ds = diff(s) ./ diff(rhs .* D.h[end])
#  rhd = D.h[2:end-1] ./ D.h[end]
  rhd = rhs[1:end-1] + diff(rhs) ./ 2
  mdh = minimum([dh[1:end-1] dh[2:end]], dims=2)
  [s dh rhs [ds mdh rhd; NaN NaN NaN]]
end |> x -> reduce(vcat, x)

using Interpolations
ra = Matrix{Any}(undef, 2, 2)
t1 = DataFrame(h=round.(res[:, 3] .* 2.5, digits=1) ./ 2.5, s=res[:, 1], w=res[:, 2])
t1 = groupby(t1, :h, sort=true)
t1 = combine(t1, [:s, :w] => ((s, w) -> [[quantile(s, weights(w), [0.025, 0.5, 0.975]); length(s)]]) => [:l, :m, :u, :n])
filter!(:n => >(300), t1)
#
ra[1, 1] = linear_interpolation(t1.h, t1.l, extrapolation_bc=Line())
ra[1, 2] = linear_interpolation(t1.h, t1.u, extrapolation_bc=Line())
#
Plots.scalefontsizes(1.4)
scatter(res[:, 1], res[:, 3], xlabel="Slope", ylabel="Relative Height", color=:black, ms=1, ma=0.1, label=nothing, xlim=[-0.2, 0.04], ylim=[0, 1], widen=true, fmt=:png)
vline!([0], color=:gray, lw=2, label=nothing)
plot!(Matrix(t1[:, [:m, :l, :u]]), repeat(t1.h, 1, 3), color=[:green :red :red], lw=2, label=nothing)
plot!([0 0], [0 0], label=["Median" "Q.025 Q.975"], lw=2, color=[:green :red])
savefig("slope.png")
Plots.scalefontsizes()

t1 = DataFrame(h=round.(res[:, 6] .* 2.5, digits=1) ./ 2.5, s=res[:, 4], w=res[:, 5])
subset!(t1, All() .=> ByRow(!isnan))
t1 = groupby(t1, :h, sort=true)
t1 = combine(t1, [:s, :w] => ((s, w) -> [[quantile(s, weights(w), [0.025, 0.5, 0.975]); length(s)]]) => [:l, :m, :u, :n])
filter!(:n => >(300), t1)
#
ra[2, 1] = linear_interpolation(t1.h, t1.l, extrapolation_bc=Line())
ra[2, 2] = linear_interpolation(t1.h, t1.u, extrapolation_bc=Line())
#
Plots.scalefontsizes(1.4)
scatter(res[:, 4], res[:, 6], xlabel="Slope Change", ylabel="Relative Height", color=:black, ms=1, ma=0.1, label=nothing, xlim=[-0.05, 0.1], ylim=[0, 1], widen=true, fmt=:png)
vline!([0], color=:gray, lw=2, label=nothing)
plot!(Matrix(t1[:, [:m, :l, :u]]), repeat(t1.h, 1, 3), color=[:green :red :red], lw=2, label=nothing)
plot!([0 0], [0 0], label=["Median" "Q.025 Q.975"], lw=2, color=[:green :red])
savefig("dSlope.png")
Plots.scalefontsizes()

s2 = function (D, i)
  if (length(i) > 0)
    ip = i .+ 1
    cnt = countmap([i; ip])
    cnt = Dict(keys(cnt) .=> 1.0 ./ values(cnt))
    w = get.(Ref(cnt), i, 0)
    wp = get.(Ref(cnt), ip, 0)
    D = vcat(D[Not(unique([i; ip])), :],
      DataFrame(h=(D.h[i] .* w .+ D.h[ip] .* wp) ./ (w .+ wp),
        d=(D.d[i] .* w .+ D.d[i.+1] .* wp) ./ (w .+ wp)))
    sort!(D, :h)
  end
  D
end

Plots.scalefontsizes(1.4)
D = DataFrame(h=[5, 10, 11, 15, 20], d=[18, 13, 14, 9, 5])
i = 2
p1 = plot(D.d, D.h, marker=true, xlabel="", ylabel="Height [m]", label=nothing, color=:black)
plot!(p1, D.d[i:i+1], D.h[i:i+1], marker=true, lw=2, color=:red, label=nothing)
D2 = s2(D, i)
p2 = plot(D.d, D.h, marker=true, xlabel="Diameter [cm]", ylabel="Height [m]", label=nothing, color=:gray)
plot!(p2, D2.d, D2.h, marker=true, label=nothing, color=:black)
plot!(p2, D2.d[i-1:i+1], D2.h[i-1:i+1], lw=2, color=:green, label=nothing)
scatter!(p2, [D2.d[i]], [D2.h[i]], label=nothing, color=:green, ms=5)
#
i = 3
p3 = plot(D.d, D.h, marker=true, xlabel="", ylabel="", label=nothing, color=:black)
plot!(p3, D.d[i-1:i+1], D.h[i-1:i+1], marker=true, lw=2, color=:red, label=nothing)
D2 = s2(D, [i - 1, i])
p4 = plot(D.d, D.h, marker=true, xlabel="Diameter [cm]", ylabel="", label=nothing, color=:gray)
plot!(p4, D2.d, D2.h, marker=true, label=nothing, color=:black)
plot!(p4, D2.d[i-2:i+1], D2.h[i-2:i+1], lw=2, color=:green, label=nothing)
scatter!(p4, D2.d[i-1:i], D2.h[i-1:i], label=nothing, color=:green, ms=5)
plot(p1, p3, p2, p4)
savefig("sampleSmooth.pdf")
Plots.scalefontsizes()

res = map(eachcol(x), eachcol(y)) do xd, yh
  i = findall(.!isnan.(xd))
  D = combine(groupby(DataFrame(h=yh[i], d=xd[i]), :h, sort=true), :d => mean => :d)
  dh = diff(D.h)
  s = diff(D.d) ./ 50.0 ./ dh
  rhs = (D.h[1:end-1] + dh ./ 2) ./ D.h[end]
  i = findall(.!(ra[1, 1].(rhs) .<= s .<= ra[1, 2].(rhs)))
  D = s2(D, i)
  dh = diff(D.h)
  s = diff(D.d) ./ 50.0 ./ dh
#  ds = diff(s)
#  rhd = D.h[2:end-1] ./ D.h[end]
  rhs = (D.h[1:end-1] + dh ./ 2) ./ D.h[end]
  ds = diff(s) ./ diff(rhs .* D.h[end])
  rhd = rhs[1:end-1] + diff(rhs) ./ 2
  i = findall(.!(ra[2, 1].(rhd) .<= ds .<= ra[2, 2].(rhd))) .+ 1
  s2(D, unique([i .- 1; i]))
end
#
#d0 und d-1 dazuhaengen
res = map(res, eachrow(Matrix(d[:, 1:2]))) do x, y
  combine(groupby(vcat(DataFrame(h=[-1.0, 0.0], d=y), x), :h, sort=true), :d => mean => :d)
end

D = res[1]
#D = res[4]
dh = diff(D.h)
s = diff(D.d) ./ 50.0 ./ dh  # Neigung bezogen auf Markroehre
rhs = (D.h[1:end-1] + dh ./ 2) ./ D.h[end]
scatter(rhs, s)
ds = diff(s) ./ diff(rhs .* D.h[end])
rhd = rhs[1:end-1] + diff(rhs) ./ 2
scatter(rhd, ds)
idx = findfirst(i -> ds[i] > 0 && ds[i+1] < 0, 1:length(ds)-1)
if idx !== nothing && idx < length(ds) && all(ds[1:idx] .> 0)
    y1, y2 = ds[idx], ds[idx+1]
    x1, x2 = rhd[idx], rhd[idx+1]
    wp = x1 - y1 * (x2 - x1) / (y2 - y1)
else
    wp = 0.0
end
scatter!([wp], [0.])

wp = map(res) do D
  dh = diff(D.h)
  s = diff(D.d) ./ 50.0 ./ dh
  rhs = (D.h[1:end-1] + dh ./ 2) ./ D.h[end]
  ds = diff(s) ./ diff(rhs .* D.h[end])
  rhd = rhs[1:end-1] + diff(rhs) ./ 2
  idx = findfirst(i -> ds[i] > 0 && ds[i+1] < 0, 1:length(ds)-1)
  if idx !== nothing && idx < length(ds) && all(ds[1:idx] .> 0)
    y1, y2 = ds[idx], ds[idx+1]
    x1, x2 = rhd[idx], rhd[idx+1]
    wp = x1 - y1 * (x2 - x1) / (y2 - y1)
  else
    wp = missing
  end
  wp
end
wp[coalesce.(wp .< 0.0, false)] .= missing

i = round.(Int, range(1, size(wp, 1), 100))
plot(sort(wp)[i], i)

rp(wp, d.DBH, .3)
rp(wp, ht.H_total, .3)
rp(wp .* ht.H_total, ht.H_total, .3)

@. model(x, p) = p[1] * x^p[2]
@. model(x, p) = p[1] * (1 .+ x)^p[2] - p[1]
i = .!ismissing.(wp) .& .!isnan.(wp)
fit = curve_fit(model, ht.H_total[i], wp[i] .* ht.H_total[i], [1., 1.])
quick_summary(fit)
# Row │ Parameter  Estimate  StdError   t_stat   p_val        Lower95   Upper95 
#   1 │ p[1]       1.8104    0.240531    7.5267  5.57332e-14  1.33892   2.28188
#   2 │ p[2]       0.358966  0.0280861  12.7809  0.0          0.303913  0.41402
Plots.scalefontsizes(1.4)
rp(-fit.resid, wp[i] .* ht.H_total[i] .+ fit.resid, ylab=L"H_{ip} - \hat{H}_{ip} [m]", xlab=L"\hat{H}_{ip} [m]")
savefig("hIpResidPred.png")
rp(-fit.resid, ht.H_total[i], ylab=L"H_{ip} - \hat{H}_{ip} [m]", xlab="Total Tree Height [m]")
savefig("hIpResidHoehe.png")
rp(-fit.resid, d.DBH[i], ylab=L"H_{ip} - \hat{H}_{ip} [m]", xlab=L"d_r [cm]")  # BHD bzw d_r
savefig("hIpResidDiameter.png")
Plots.scalefontsizes()
rp(wp .* ht.H_total, ht.H_total, .2)
plot!(0:50, model(0:50., coef(fit)))
rp(wp, ht.H_total, .2)
plot!(1:50, model(1:50, coef(fit)) ./ (1:50))
model(0.01, coef(fit))

using Integrals, Optim
v = map(res) do x
  s = interpolate(x.h, x.d, SteffenMonotonicInterpolation())
  v = solve(IntegralProblem((h, p) -> s(h)^2 * pi / 40000.0, (0.0, max(0.0, x.h[end]))), QuadGKJL()).u
  hWp = model(x.h[end], coef(fit))
  dWp = s(hWp)
  vWu = solve(IntegralProblem((h, p) -> s(h)^2 * pi / 40000.0, (0.0, hWp)), QuadGKJL()).u
  vWi = solve(IntegralProblem((h, p) -> s(h)^2 * pi / 40000.0, (hWp, x.h[end])), QuadGKJL()).u
  HT = hWp + (x.h[end] - hWp) * .6
  [hWp; x.h[end]; v; vWu; vWi; s(0.); s(1.3); s(1.3) + 1.3; s(hWp); x.d[end]; s.(min.((0.:.05:1)*x.h[end], x.h[end])); HT; s(HT)]
end |> stack

D = DataFrame(d = v[7,:], d0 = v[6,:], da = v[8,:], dWp = v[9,:], dt = v[10,:], h = v[2,:], hWp = v[1,:], vWu = v[4,:], vWi = v[5,:])

using DataFrames, Optim

struct TaperProfile
    H::Float64     # Gesamthöhe
    hw::Float64    # Übergangshöhe (hw)
    d0::Float64    # Durchmesser am Boden (m)
    dw::Float64    # Durchmesser an der Schnittstelle (m)
    dt::Float64    # Terminaler Durchmesser an der Spitze (m)
    a::Float64     # Skalierung oben
    b::Float64     # Exponent oben
    c1::Float64    # Skalierung basal
    c2::Float64    # Exponent basal
end

function volume_upper_calc(b, dw, dt, H, hw)
    a = (dw - dt) / (H - hw)^b
    return (pi/4) * ((a^2 * (H - hw)^(2b + 1)) / (2b + 1) + (2 * a * dt * (H - hw)^(b + 1)) / (b + 1) + (dt^2 * (H - hw)))
end

function volume_lower_calc(d0, dw, sw, hw)
    if d0 <= dw 
        return -1.0 
    end
    
    c2 = (abs(sw) * hw) / (d0 - dw)
    if c2 <= 0.001 
        return -1.0 
    end 
    
    c1 = (d0 - dw) / (hw^c2)
    
    return (pi / 4.0) * (
        d0^2 * hw - 
        (2.0 * d0 * c1 * hw^(c2 + 1.0)) / (c2 + 1.0) + 
        (c1^2 * hw^(2.0 * c2 + 1.0)) / (2.0 * c2 + 1.0)
    )
end

function getVolume(p::TaperProfile, h1::Float64, h2::Float64)
    hs = max(0.0, min(p.H, h1))
    he = max(0.0, min(p.H, h2))
    if hs >= he return 0.0 end
    
    if hs < p.hw && he > p.hw
        return getVolume(p, hs, p.hw) + getVolume(p, p.hw, he)
    end
    
    v_raw = 0.0
    if he <= p.hw
        F_low(h) = p.d0^2 * h + (2.0 * p.d0 * p.c1 * h^(p.c2 + 1.0)) / (p.c2 + 1.0) + (p.c1^2 * h^(2.0 * p.c2 + 1.0)) / (2.0 * p.c2 + 1.0)
        v_raw = (pi / 4.0) * (F_low(he) - F_low(hs))
    else
        z1 = p.H - hs  
        z2 = p.H - he  
        F_up(z) = p.dt^2 * z + (2.0 * p.a * p.dt * z^(p.b + 1.0)) / (p.b + 1.0) + (p.a^2 * z^(2.0 * p.b + 1.0)) / (2.0 * p.b + 1.0)
        v_raw = (pi / 4.0) * (F_up(z1) - F_up(z2))
    end
    return v_raw
end

function getParam(bhd_cm::Float64, H::Float64, b_val::Float64, c2_val::Float64)
    dbh = bhd_cm / 100.0
    dt = H >= 1.3 ? 0.008 : (0.001 + 0.007 * (1.0 - ((1.3 - H)/1.3)^2))
    hw = 1.8104 * (1.0 + H)^0.3590 - 1.8104
    
    local_c2 = max(0.02, c2_val)
    local_d0, local_c1, local_a, local_dip = 0.0, 0.0, 0.0, 0.0

    if 1.3 >= hw
        local_a = (dbh - dt) / (H - 1.3)^b_val
        local_dip = local_a * (H - hw)^b_val + dt
        sip = -local_a * b_val * (H - hw)^(b_val - 1.0)
        local_c1 = sip / (local_c2 * hw^(local_c2 - 1.0))
        local_d0 = local_dip - local_c1 * hw^local_c2
        
    else
        X = (H - hw)^b_val
        Y = -b_val * (H - hw)^(b_val - 1.0)

        Klammerterm = X - ((Y * hw) / local_c2) * (1.0 - (1.3 / hw)^local_c2)
        
        local_a = (dbh - dt) / max(1e-6, Klammerterm)
        local_dip = local_a * X + dt
        sip = local_a * Y
        
        local_c1 = sip / (local_c2 * hw^(local_c2 - 1.0))
        local_d0 = local_dip - local_c1 * hw^local_c2
    end
    
    return TaperProfile(H, hw, local_d0, local_dip, dt, local_a, b_val, local_c1, local_c2)
end

function solve_tree_global_clean(row)
    H, hw = row.h, row.hWp
    dw, dt = row.dWp / 100.0, row.dt / 100.0
    v_total_target = row.vWu + row.vWi
    d0_meas = row.d0 / 100.0
    
    bhd_cm = row.d 

    function cost_stage1(p::Vector{Float64})
        b = p[1]
        d0 = p[2]
        
        if b < 0.01 || b > 1.0 || d0 <= dw || d0 > dw * 4.0
            return 1e15
        end
        v_up = volume_upper_calc(b, dw, dt, H, hw)
        a = (dw - dt) / (H - hw)^b
        sw = -a * b * (H - hw)^(b - 1)
        
        c2 = (sw * hw) / (dw - d0)
        if c2 < 0.01 return 1e15 end 

        v_low = volume_lower_calc(d0, dw, sw, hw)
        if v_low < 0.0 return 1e15 end 
        
        v_total_pred = v_up + v_low
        return (log(v_total_pred) - log(v_total_target))^2 * 1e9 + (d0 - d0_meas)^2 * 1e2 + (b - 1.0)^2 * 0.01
    end

    start_p1 = [0.6, d0_meas]
    res1 = optimize(cost_stage1, start_p1, NelderMead(), Optim.Options(iterations=2000))
    b_stage1, d0_stage1 = Optim.minimizer(res1)[1], Optim.minimizer(res1)[2]
    
    a_s1 = (dw - dt) / (H - hw)^b_stage1
    sw_s1 = -a_s1 * b_stage1 * (H - hw)^(b_stage1 - 1)
    c2_stage1 = (sw_s1 * hw) / (dw - d0_stage1)

    start_p2 = [1.0, 1.0] 
    
    function cost_stage2(p::Vector{Float64})
      f_b = p[1]
      f_c2 = p[2]

      ofr = 0.
        
        if f_b < 0.10 || f_b > 2.0 || f_c2 < 0.10 || f_c2 > 3.0
          ofr += 1e15
        end
        
        b_tuned = b_stage1 * f_b
        c2_tuned = c2_stage1 * f_c2
        
        if b_tuned < 0.01 || b_tuned > 1.0 || c2_tuned < 0.01 || c2_tuned > 1.0
          ofr += 1e15
        end
      
        p_profile = getParam(bhd_cm, H, b_tuned, c2_tuned)
        v_total_pred = getVolume(p_profile, 0.0, H)
        
        v_err = (log(v_total_pred) - log(v_total_target))^2 * 1e14
        tuning_penalty = (f_b - 1.0)^2 * 1000.0 + (f_c2 - 1.0)^2 * 1000.0
        
      return v_err + tuning_penalty + ofr
    end
    
    res2 = optimize(cost_stage2, start_p2, NelderMead(), Optim.Options(iterations=2000, g_tol=1e-12))
    f_b_fin, f_c2_fin = Optim.minimizer(res2)[1], Optim.minimizer(res2)[2]
    
    b_fin = b_stage1 * f_b_fin
    c2_fin = c2_stage1 * f_c2_fin
    
    p_fin = getParam(bhd_cm, H, b_fin, c2_fin)
    sip_fin = -p_fin.a * b_fin * (H - p_fin.hw)^(b_fin - 1.0)
    
  return (p_fin.a, b_fin, sip_fin, p_fin.d0 * 100.0, p_fin.c1, c2_fin, Optim.minimum(res2), b_stage1, c2_stage1)
end

D = DataFrame(d = v[7,:], d0 = v[6,:], da = v[8,:], dWp = v[9,:], dt = v[10,:], h = v[2,:], hWp = v[1,:], vWu = v[4,:], vWi = v[5,:])

resForm = map(solve_tree_global_clean, eachrow(D)) |> stack

D = DataFrame(h = ht.H_total, dRef = d.DBH .+ 1.3, b = resForm[2,:], c2 = resForm[6,:], bhd = d.DBH, v = v[3,:], bS1 = resForm[8,:], c2S1 = resForm[9,:], d1 = d.DBH, d2 = v[17,:], h1 = 1.3, h2 = 0.3 .* ht.H_total, d3 = v[24,:], h3 = 0.65 .* ht.H_total, d4 = v[13,:], h4 = 0.1 .* ht.H_total)

D.bO = copy(D.b)
D.c2O = copy(D.c2)

D.b[D.b .> .99] .= .99
D.c2[D.c2 .> .99] .= .99

p = getParam(30., 20., .5, .5)
getVolume(p, 0., 20.)
D.vEst = map(eachrow(D)) do x
  p = getParam(x.bhd, x.h, x.bO, x.c2O)
  getVolume(p, 0., x.h)
end
rp(D.vEst ./ D.v, D.vEst)
rp(D.v ./ D.vEst, D.vEst)
std(D.v ./ D.vEst)
std(D.vEst ./ D.v)  # 0.03748412187743803
D.r = D.v ./ D.vEst
ql, qh = quantile(D.r, [.005, .995])
i = D.r .>= ql .&& D.r .<= qh
std(D.r[i])  # 0.0273961716410582
D.w = D.r .>= 0.99 .&& D.r .<= 1.01
sum(D.w)  # 11946

rp(D.b, D.h)
rp(D.b, D.dRef)
rp(D.b, D.h ./ D.dRef)
rp(D.b, D.dRef ./ D.h)
reg = glm(@formula(b ~ dRef / h), D, Normal(), LogitLink())
reg = glm(@formula(b ~ dRef / h + h + log(dRef)), D, Normal(), LogitLink())
reg = glm(@formula(b ~ inv(h) + log1p(dRef)), D, Normal(), LogitLink())
reg = glm(@formula(b ~ inv(h) + log1p(dRef)), D, Normal(), LogitLink(), wts=D.w)
#                Coef.  Std. Error       z  Pr(>|z|)  Lower 95%  Upper 95%
#(Intercept)  -4.31616   0.0578332  -74.63    <1e-99   -4.42951   -4.20281
#inv(h)       21.7497    0.297931    73.00    <1e-99   21.1658    22.3337
#log1p(dRef)   1.0813    0.0139944   77.27    <1e-99    1.05388    1.10873
D.bEst = predict(reg)
Mb1d = reg
Plots.scalefontsizes(1.4)
rp(response(reg)[D.w] .- predict(reg)[D.w], predict(reg)[D.w], ylab=L"b - \hat{b}", xlab=L"\hat{b}")
savefig("bResidPred.png")
rp(response(reg)[D.w] .- predict(reg)[D.w], D.h[D.w], ylab=L"b - \hat{b}", xlab="Total Tree Height [m]")
savefig("bResidH.png")
rp(response(reg)[D.w] .- predict(reg)[D.w], D.dRef[D.w], ylab=L"b - \hat{b}", xlab=L"d_{r} [cm]")
savefig("bResidDref.png")
Plots.scalefontsizes()
rp(D.b, D.h)
plot!(0:50, predict(reg, DataFrame(h=0:50, dRef=0:50)), lw=2)
plot!(0:50, predict(reg, DataFrame(h=0:50, dRef=1)), lw=2)
plot!(0:50, predict(reg, DataFrame(h=0:50, dRef=20)), lw=2)
plot!(0:50, predict(reg, DataFrame(h=0:50, dRef=30)), lw=2)
plot!(0:50, predict(reg, DataFrame(h=0:50, dRef=80)), lw=2)
plot!(0:50, predict(reg, DataFrame(h=0:50, dRef=200)), lw=2)
rp(D.b, D.dRef)
plot!(0:100, predict(reg, DataFrame(h=0:100, dRef=0:100)), lw=2)
plot!(0:50, predict(reg, DataFrame(h=0:50, dRef=1)), lw=2)
plot!(0:50, predict(reg, DataFrame(h=0:50, dRef=30)), lw=2)
plot!(0:50, predict(reg, DataFrame(h=0:50, dRef=80)), lw=2)

rp(D.c2, D.h)
rp(D.c2, D.dRef)
rp(D.c2, D.h ./ D.dRef)
rp(D.c2, D.dRef ./ D.h)
reg = glm(@formula(c2 ~ log(h) + dRef + inv(dRef) + log1p(dRef)), D, Normal(), LogitLink())
reg = glm(@formula(c2 ~ log1p(h) + log(inv(1+dRef))), D, Normal(), LogitLink())
reg = glm(@formula(c2 ~ log1p(h) + inv(dRef)), D, Normal(), LogitLink())
reg = glm(@formula(c2 ~ log1p(h) + inv(3. + dRef)), D, Normal(), LogitLink())
reg = glm(@formula(c2 ~ log1p(h) + log1p(dRef)), D, Normal(), LogitLink())
reg = glm(@formula(c2 ~ log1p(h) + log1p(dRef)), D, Normal(), LogitLink(), wts=D.w)
#                 Coef.  Std. Error       z  Pr(>|z|)  Lower 95%  Upper 95%
#(Intercept)   7.0974     0.0520578  136.34    <1e-99   6.99536    7.19943
#log1p(h)     -3.17164    0.0370744  -85.55    <1e-99  -3.2443    -3.09897
#log1p(dRef)   0.631264   0.0281408   22.43    <1e-99   0.576109   0.686419
D.c2Est = predict(reg)
Mc21d = reg
Plots.scalefontsizes(1.4)
rp(response(reg)[D.w] .- predict(reg)[D.w], predict(reg)[D.w], ylab=L"c_2 - \hat{c}_2", xlab=L"\hat{c}_2")
savefig("c2ResidPred.png")
rp(response(reg)[D.w] .- predict(reg)[D.w], D.h[D.w], ylab=L"c_2 - \hat{c}_2", xlab="Total Tree Height [m]")
savefig("c2ResidH.png")
rp(response(reg)[D.w] .- predict(reg)[D.w], D.dRef[D.w], ylab=L"c_2 - \hat{c}_2", xlab=L"d_{r} [cm]")
savefig("c2ResidDref.png")
Plots.scalefontsizes()

rp(response(reg) .- predict(reg), predict(reg))
rp(response(reg) .- predict(reg), D.h)
rp(response(reg) .- predict(reg), D.dRef)
rp(D.c2, D.h)
plot!(0:50, predict(reg, DataFrame(h=0:50, dRef=0:50)), lw=2)
plot!(0:50, predict(reg, DataFrame(h=0:50, dRef=1)), lw=2)
plot!(0:50, predict(reg, DataFrame(h=0:50, dRef=20)), lw=2)
plot!(0:50, predict(reg, DataFrame(h=0:50, dRef=30)), lw=2)
plot!(0:50, predict(reg, DataFrame(h=0:50, dRef=80)), lw=2)
plot!(0:50, predict(reg, DataFrame(h=0:50, dRef=200)), lw=2)
plot!(0:50, predict(reg, DataFrame(h=0, dRef=0:50)), lw=2)
predict(reg, DataFrame(h=0.1, dRef=0.1))

D.vEst = map(eachrow(D)) do x
  p = getParam(x.bhd, x.h, x.bEst, x.c2Est)
  getVolume(p, 0., x.h)
end
D.r = D.v ./ D.vEst
ql, qh = quantile(D.r, [.005, .995])
i = D.r .>= ql .&& D.r .<= qh
std(D.r[i])  # 0.07171019787999719
Plots.scalefontsizes(1.4)
rp(D.r[i], D.vEst[i], h=1, trx=:log10, ylab=L"V / \hat{V}", xlab=L"\hat{V} [m^3]")
savefig("VResidPred.png")
rp(D.r[i], D.h[i], h=1, ylab=L"V / \hat{V}", xlab="Total Tree Height [m]")
savefig("VResidH.png")
j = i .&& D.dRef .< 70
rp(D.r[j], D.dRef[j], h=1, ylab=L"V / \hat{V}", xlab=L"d_r [cm]")
savefig("VResidDRef.png")
Plots.scalefontsizes()

function p2cm(p::TaperProfile)
    return TaperProfile(
        p.H,          # H bleibt gleich
        p.hw,         # hw bleibt gleich
        p.d0 * 100.0, # d0 in cm
        p.dw * 100.0, # dw in cm
        p.dt * 100.0, # dt in cm
        p.a * 100.,          # a bleibt gleich
        p.b,          # b bleibt gleich
        p.c1 * 100.,         # c1 bleibt gleich
        p.c2          # c2 bleibt gleich
    )
end

function getVolume2(p::TaperProfile, h1::Float64, h2::Float64)
    hs = max(0.0, min(p.H, h1))
    he = max(0.0, min(p.H, h2))
    if hs >= he return 0.0 end
    
    if hs < p.hw && he > p.hw
        return getVolume2(p, hs, p.hw) + getVolume2(p, p.hw, he)
    end
    
    v_raw = 0.0
    if he <= p.hw
      F_low(h) = p.d0^2 * h + (2.0 * p.d0 * p.c1 * h^(p.c2 + 1.0)) / (p.c2 + 1.0) + (p.c1^2 * h^(2.0 * p.c2 + 1.0)) / (2.0 * p.c2 + 1.0)
      v_raw = (pi / 40000.0) * (F_low(he) - F_low(hs))
    else
      z1 = p.H - hs  
      z2 = p.H - he  
      F_up(z) = p.dt^2 * z + (2.0 * p.a * p.dt * z^(p.b + 1.0)) / (p.b + 1.0) + (p.a^2 * z^(2.0 * p.b + 1.0)) / (2.0 * p.b + 1.0)
      v_raw = (pi / 40000.0) * (F_up(z1) - F_up(z2))
    end
    return v_raw
end

function get_diameter_at_height(p::TaperProfile, h::Float64)
  if h < 0.0 || h > p.H return 0.0 end
  return h < p.hw ? p.d0 + p.c1 * h^p.c2 : p.a * (p.H - h)^p.b + p.dt
end

D.d03Est = map(eachrow(D)) do x
  p = getParam(x.bhd, x.h, x.bEst, x.c2Est)
  get_diameter_at_height(p2cm(p), x.h2)
end
rp(D.d2 .- D.d03Est, D.d03Est, ylab=L"s_{0} - \hat{s}_{0}", xlab=L"\hat{s}_{0}")
D.r = D.d2 ./ D.d03Est
ql, qh = quantile(D.r, [.005, .995])
i = D.r .>= ql .&& D.r .<= qh
std(D.r[i])  # 0.044504231473490884
Plots.scalefontsizes(1.4)
rp(D.r[i], D.d03Est[i], h=1, ylab=L"d_{0.3h} / \hat{d}_{0.3h}", xlab=L"\hat{d}_{0.3h} [cm]")
savefig("d03Est.png")
Plots.scalefontsizes()

D.d065Est = map(eachrow(D)) do x
  p = getParam(x.bhd, x.h, x.bEst, x.c2Est)
  get_diameter_at_height(p2cm(p), x.h3)
end
rp(D.d3 .- D.d065Est, D.d065Est)
D.r = D.d3 ./ D.d065Est
ql, qh = quantile(D.r, [.005, .995])
i = D.r .>= ql .&& D.r .<= qh
std(D.r[i])  # 0.0777538250714001
Plots.scalefontsizes(1.4)
rp(D.r[i], D.d065Est[i], h=1, ylab=L"d_{0.65h} / \hat{d}_{0.65h}", xlab=L"\hat{d}_{0.65h} [cm]")
savefig("d065Est.png")
Plots.scalefontsizes()

function plot_tree(p::TaperProfile; color=:black)
  h_steps = range(0.0, p.H, length=150)
  d_ob_vals = [get_diameter_at_height(p, h) for h in h_steps]
  r_ob = d_ob_vals ./ 2.0
  x_wood = vcat(r_ob, reverse(-r_ob))
  y_wood = vcat(h_steps, reverse(h_steps))
  plot!(x_wood, y_wood, color=color, label=nothing)
end

function getDiameterAtHeight(p::TaperProfile, h::Float64)
    if h <= p.hw
        return p.d0 + p.c1 * h^p.c2
    else
        return p.a * (p.H - h)^p.b + p.dt
    end
end

using Roots
using Optim

function getParam2(bhd_cm::Float64, H::Float64, b_glm::Float64, c2_glm::Float64, d1_cm::Float64, d2_cm::Float64, h1_m::Float64, h2_m::Float64)
    # =================================================================
    # 1. INTERNE SKALIERUNG & AUTOMATISCHE SORTIERUNG NACH HÖHE
    # =================================================================
    hw = 1.8104 * (1.0 + H)^0.3590 - 1.8104
    
    # Sortierung, damit h_oben immer der höhere und h_unten der tiefere Punkt ist
    if h1_m >= h2_m
        h_oben, h_unten = h1_m, h2_m
        d_oben, d_unten = d1_cm, d2_cm
    else
        h_oben, h_unten = h2_m, h1_m
        d_oben, d_unten = d2_cm, d1_cm
    end

    dt_cm = H >= 1.3 ? 0.8 : 0.1 + 0.7 * (1.0 - ((1.3 - H)/1.3)^2)
    
    d_ref_cm = abs(H - 1.3) < 1e-3 ? bhd_cm + 1.3 : H * 100.0
    limit_c0 = d_ref_cm * 1.8

    local_c2 = c2_glm
    local_b = b_glm

    local_a = 0.0
    local_c1 = 0.0
    local_c0 = 0.0
    local_dip = 0.0  # entspricht dw

    # Fall-Entscheidung basierend auf den sortierten Höhen
    if h_unten >= hw
        # Fall A: Beide Punkte liegen über hw
        term_num = max(1e-5, d_unten - dt_cm) / max(1e-5, d_oben - dt_cm)
        term_den = (H - h_unten) / (H - h_oben)
        
        local_b = clamp(log(term_num) / log(term_den), 0.1, 0.95)
        
        local_a = (d_unten - dt_cm) / max(1e-6, H - h_unten)^local_b
        local_dip = local_a * max(1e-6, H - hw)^local_b + dt_cm
        sip = -local_a * local_b * max(1e-6, H - hw)^(local_b - 1.0)
        
        local_c1 = (sip * hw^(1.0 - local_c2)) / local_c2
        local_c0 = local_dip - local_c1 * hw^local_c2
        
        if local_c0 > limit_c0 || local_c0 < 0.0
            local_c2 = 2.0
            local_c1 = (sip * hw^(1.0 - local_c2)) / local_c2
            local_c0 = local_dip - local_c1 * hw^local_c2
        end
        
    elseif h_unten < hw && h_oben >= hw
        # Fall B: Ein Punkt unten, ein Punkt oben
        # Wir suchen das b, das den unteren Punkt (d_unten) exakt trifft
        find_b_trajectory_error = (b_trial) -> begin
            a_t = (d_oben - dt_cm) / max(1e-6, H - h_oben)^b_trial
            dw_t = a_t * max(1e-6, H - hw)^b_trial + dt_cm
            sw_t = -a_t * b_trial * max(1e-6, H - hw)^(b_trial - 1.0)
            c1_t = sw_t / (local_c2 * hw^(local_c2 - 1.0))
            c0_t = dw_t - (sw_t * hw) / local_c2
            return (c0_t + c1_t * h_unten^local_c2) - d_unten
        end

        err_low  = find_b_trajectory_error(0.05)
        err_high = find_b_trajectory_error(1.95)

        if sign(err_low) != sign(err_high)
            local_b = find_zero(find_b_trajectory_error, (0.05, 1.95), Bisection())
        else
            local_b = b_glm
        end

        local_a = (d_oben - dt_cm) / max(1e-6, H - h_oben)^local_b
        local_dip = local_a * max(1e-6, H - hw)^local_b + dt_cm
        sip = -local_a * local_b * max(1e-6, H - hw)^(local_b - 1.0)
        local_c1 = sip / (local_c2 * hw^(local_c2 - 1.0))
        local_c0 = local_dip - (sip * hw) / local_c2
        
        # Abfangjäger gegen unphysikalische Verläufe bei Kleinstbäumen
        if local_c0 > limit_c0 || local_c0 < 0.0
            local_c0 = clamp(local_c0, 0.0, limit_c0)
            local_c1 = (local_dip - local_c0) / max(1e-6, hw^local_c2)
        end

    else
        # Fall C: Beide Punkte unten
        f_opt_c2 = (c2_trial) -> begin
            denom = max(1e-6, h_oben^c2_trial - h_unten^c2_trial)
            c1_t = (d_oben - d_unten) / denom
            c0_t = d_unten - c1_t * h_unten^c2_trial
            
            dw_t = c0_t + c1_t * hw^c2_trial
            sw_t = c1_t * c2_trial * hw^(c2_trial - 1.0)
            a_t = (dw_t - dt_cm) / max(1e-6, H - hw)^local_b
            sw_expected = -a_t * local_b * max(1e-6, H - hw)^(local_b - 1.0)
            
            return (sw_t - sw_expected)^2
        end

        res = optimize(f_opt_c2, 0.1, 1.5, Brent())
        local_c2 = Optim.minimizer(res)
        
        denom = max(1e-6, h_oben^local_c2 - h_unten^local_c2)
        local_c1 = (d_oben - d_unten) / denom
        local_c0 = d_unten - local_c1 * h_unten^local_c2
        local_dip = local_c0 + local_c1 * hw^local_c2
        local_a = (local_dip - dt_cm) / max(1e-6, H - hw)^local_b
    end

    # =================================================================
    # 2. DIREKTE RÜCKGABE
    # =================================================================
    return TaperProfile(
        H, 
        hw, 
        local_c0,   
        local_dip,  
        dt_cm,      
        local_a,    
        local_b,    
        local_c1,   
        local_c2    
    )
end
  
function getParam3(bhd_cm::Float64, H::Float64, b_glm::Float64, c2_glm::Float64, 
                          d1_cm::Float64, d2_cm::Float64, d3_cm::Float64, 
                          h1_m::Float64, h2_m::Float64, h3_m::Float64)
    # =================================================================
    # 1. INTERNE SKALIERUNG & STRIKTE SORTIERUNG NACH HÖHE
    # =================================================================
    hw = 1.8104 * (1.0 + H)^0.3590 - 1.8104
    
    pts = sort([(h1_m, d1_cm), (h2_m, d2_cm), (h3_m, d3_cm)], by = x -> x)
    
    h_meas = [pts[1][1], pts[2][1], pts[3][1]]
    d_meas = [pts[1][2], pts[2][2], pts[3][2]]
    d_meas_sq = d_meas .^ 2
    
    dt_cm = H >= 1.3 ? 0.8 : 0.1 + 0.7 * (1.0 - ((1.3 - H)/1.3)^2)
    
    d_ref_cm = abs(H - 1.3) < 1e-3 ? bhd_cm + 1.3 : H * 100.0
    limit_c0 = d_ref_cm * 1.8
    
    pts_below = sum(h_meas .<= hw)
    pts_above = sum(h_meas .> hw)
    
    local_b  = b_glm
    local_c2 = c2_glm
    local_a  = 0.0
    local_c1 = 0.0
    local_c0 = 0.0

    # =================================================================
    # SONDERFALL A: Alle 3 Messpunkte liegen unten
    # =================================================================
    if pts_above == 0
        lower_bounds = [0.1, 0.01, 1.0]           
        upper_bounds = [1.2,  15.0, limit_c0]      
        initial_guess = [b_glm, c2_glm, min(limit_c0, d_meas[1])]
        
        f_opt_all_below_3d = (x) -> begin
            b_t, c2_t, c0_t = x[1], x[2], x[3]
            
            A = hw^c2_t
            B = (c2_t * hw^(c2_t - 1.0) * max(1e-6, H - hw)) / b_t
            Nenner = A + B
            
            c1_t = (dt_cm - c0_t) / max(1e-6, Nenner)
            
            total_error = 0.0
            for i in 1:3
                d_mod_sq = (c0_t + c1_t * h_meas[i]^c2_t)^2
                total_error += (d_mod_sq - d_meas_sq[i])^2
            end
            return total_error
        end
        
        res = optimize(f_opt_all_below_3d, lower_bounds, upper_bounds, initial_guess, Fminbox(NelderMead()))
        best = Optim.minimizer(res)
        local_b, local_c2, local_c0 = best[1], best[2], best[3]
        
        A_f = hw^local_c2
        B_f = (local_c2 * hw^(local_c2 - 1.0) * max(1e-6, H - hw)) / local_b
        local_c1 = (dt_cm - local_c0) / max(1e-6, A_f + B_f)
        
        dw_final = local_c0 + local_c1 * hw^local_c2
        local_a = (dw_final - dt_cm) / max(1e-6, H - hw)^local_b

    # =================================================================
    # SONDERFALL B: Alle 3 Messpunkte liegen oben
    # =================================================================
    elseif pts_below == 0
        a_init = (d_meas[1] - dt_cm) / max(1e-6, H - h_meas[1])^b_glm
        a_init = clamp(a_init, 0.1, 500.0)
        
        lower_bounds = [0.1, 1e-4, 0.01]
        upper_bounds = [1.2, 2000.0, 15.0] 
        initial_guess = [b_glm, a_init, c2_glm]
        
        f_opt_all_above_3d = (x) -> begin
            b_t, a_t, c2_t = x[1], x[2], x[3]
            
            total_error = 0.0
            for i in 1:3
                d_mod = a_t * max(1e-6, H - h_meas[i])^b_t + dt_cm
                total_error += (d_mod^2 - d_meas_sq[i])^2
            end
            return total_error
        end
        
        res = optimize(f_opt_all_above_3d, lower_bounds, upper_bounds, initial_guess, Fminbox(NelderMead()))
        best = Optim.minimizer(res)
        local_b, local_a, local_c2 = best[1], best[2], best[3]
        
        dw_final = local_a * max(1e-6, H - hw)^local_b + dt_cm
        sw_final = -local_a * local_b * max(1e-6, H - hw)^(local_b - 1.0)
        
        local_c1 = sw_final / (local_c2 * hw^(local_c2 - 1.0))
        local_c0 = dw_final - local_c1 * hw^local_c2
        
        if local_c0 > limit_c0 || local_c0 < 0.0
            local_c0 = clamp(local_c0, 0.0, limit_c0)
            local_c1 = (dw_final - local_c0) / hw^local_c2
        end
    # =================================================================
    # NORMALFALL: ABSOLUTE QUADRATE MIT VOLLER MATHEMATISCHER FREIHEIT
    # =================================================================
    else
        # Weite, stabile Grenzen erlauben dem Löser den exakten mathematischen Fit
        lower_bounds = [0.05, 1e-4,   0.01]  
        upper_bounds = [1.5,  2000.0, 15.0]   
        
        # Ein zentrierter Startwert für b sichert die stabile Konvergenz
        initial_guess = [clamp(b_glm, 0.1, 1.2), 10.0, clamp(c2_glm, 0.1, 1.5)] 
        
        f_opt_3d = (x) -> begin
            b_t, a_t, c2_t = x[1], x[2], x[3]
            
            dw_t = a_t * max(1e-6, H - hw)^b_t + dt_cm
            sw_t = -a_t * b_t * max(1e-6, H - hw)^(b_t - 1.0)
            c1_t = sw_t / (c2_t * hw^(c2_t - 1.0))
            c0_t = dw_t - c1_t * hw^c2_t
            
            if c0_t > limit_c0
                c0_t = limit_c0
                c1_t = (dw_t - limit_c0) / hw^c2_t
            elseif c0_t < 0.0
                c0_t = 0.0
                c1_t = dw_t / hw^c2_t
            end
            
            total_error = 0.0
            for i in 1:3
                d_mod = (h_meas[i] <= hw) ? (c0_t + c1_t * h_meas[i]^c2_t) : (a_t * max(1e-6, H - h_meas[i])^b_t + dt_cm)
                
                # Strikt die absoluten Quadrate minimieren (Volumen-äquivalent)
                total_error += (d_mod^2 - d_meas_sq[i])^2
            end
            return total_error
        end
        
        res = optimize(f_opt_3d, lower_bounds, upper_bounds, initial_guess, Fminbox(NelderMead()))
        best = Optim.minimizer(res)
        
        local_b, local_a, local_c2 = best[1], best[2], best[3]
        
        dw_final = local_a * max(1e-6, H - hw)^local_b + dt_cm
        sw_final = -local_a * local_b * max(1e-6, H - hw)^(local_b - 1.0)
        local_c1 = sw_final / (local_c2 * hw^(local_c2 - 1.0))
        local_c0 = dw_final - local_c1 * hw^local_c2
        
        if local_c0 > limit_c0 || local_c0 < 0.0
            local_c0 = clamp(local_c0, 0.0, limit_c0)
            local_c1 = (dw_final - local_c0) / hw^local_c2
        end
    end

    dw_val_cm = local_a * max(1e-6, H - hw)^local_b + dt_cm

    return TaperProfile(
        H, 
        hw, 
        local_c0,   
        dw_val_cm,  
        dt_cm,      
        local_a,    
        local_b,    
        local_c1,   
        local_c2    
    )
end


D.vEst2 = map(eachrow(D)) do x
  p = getParam2.(x.bhd, x.h, x.bEst, x.c2Est, x.d1, x.d2, x.h1, x.h2)
  getVolume2(p, 0., x.h)
end
D.r = D.v ./ D.vEst2
ql, qh = quantile(D.r, [.005, .995])
i = D.r .>= ql .&& D.r .<= qh
std(D.r[i])  # 0.06632768872295262
Plots.scalefontsizes(1.4)
rp(D.r[i], D.vEst2[i], h=1, trx=:log10, ylab=L"V / \hat{V}", xlab=L"\hat{V} [m^3]")
savefig("VResidPred2.png")
Plots.scalefontsizes()

D.vEst2b = map(eachrow(D)) do x
  p = getParam2.(x.bhd, x.h, x.bEst, x.c2Est, x.d1, x.d3, x.h1, x.h3)
  getVolume2(p, 0., x.h)
end
D.r = D.v ./ D.vEst2b
ql, qh = quantile(D.r, [.005, .995])
i = D.r .>= ql .&& D.r .<= qh
std(D.r[i])  # 0.04243298816973128
Plots.scalefontsizes(1.4)
rp(D.r[i], D.vEst2b[i], h=1, trx=:log10, ylab=L"V / \hat{V}", xlab=L"\hat{V} [m^3]")
savefig("VResidPred2b.png")
Plots.scalefontsizes()

D.vEst3 = map(eachrow(D)) do x
  p = getParam3.(x.bhd, x.h, x.bEst, x.c2Est, x.d1, x.d2, x.d3, x.h1, x.h2, x.h3)
  getVolume2(p, 0., x.h)
end
D.r = D.v ./ D.vEst3
ql, qh = quantile(D.r, [.005, .995])
i = D.r .>= ql .&& D.r .<= qh
std(D.r[i])  # 0.027925046999304325
Plots.scalefontsizes(1.4)
rp(D.r[i], D.vEst3[i], h=1, trx=:log10, ylab=L"V / \hat{V}", xlab=L"\hat{V} [m^3]")
savefig("VResidPred3.png")
Plots.scalefontsizes()


ND = DataFrame(h=25., bhd=25., b=[.5, .3, 1.2], c2=[.25, .3, .7])
p = p2cm.(getParam.(ND.bhd, ND.h, ND.b, ND.c2))
p1 = plot()
plot_tree.(p)
scatter!(p1, ND.bhd[1] * [-.5, .5], [1.3, 1.3], markersize=5, label=nothing, xlabel="", ylabel="Height [m]", xticks=(-10:10:10, abs.(-10:10:10)), title="1 Diameter", grid=:x, color=:black)
ND.v1 = getVolume2.(p, 0., ND.h)

ND = DataFrame(h=25., bhd=25., b=0., c2=[.3, .6, 1.], d03=20.7, h03=7.5)
p = getParam2.(ND.bhd, ND.h, ND.b, ND.c2, ND.bhd, ND.d03, 1.3, ND.h03)
p2 = plot()
plot_tree.(p)
#scatter!(p2, [ND.bhd[1], ND.d03[1]] * [-.5, .5,]', [1.3, ND.h03[1]], markersize=5, label=nothing, xlabel="Radius [cm]", xticks=(-10:10:10, abs.(-10:10:10)), showaxis=:x, yaxis=nothing, title="2 Diameters", color=:black)
scatter!(p2, [ND.bhd[1], ND.d03[1]] * [-.5, .5,]', [1.3, ND.h03[1]], markersize=5, label=nothing, xlabel="", xticks=(-10:10:10, abs.(-10:10:10)), showaxis=:x, yaxis=nothing, title="2 Diameters", color=:black)
ND.v2 = getVolume2.(p, 0., ND.h)

ND = DataFrame(h=25., bhd=25., b=0., c2=[.3, .6, 1.], d03=14.2, h03=16.25)
p = getParam2.(ND.bhd, ND.h, ND.b, ND.c2, ND.bhd, ND.d03, 1.3, ND.h03)
p2b = plot()
plot_tree.(p)
scatter!(p2b, [ND.bhd[1], ND.d03[1]] * [-.5, .5,]', [1.3, ND.h03[1]], markersize=5, label=nothing, xlabel="", xticks=(-10:10:10, abs.(-10:10:10)), showaxis=:x, yaxis=nothing, title="2 Diameters", color=:black)
ND.v2b = getVolume2.(p, 0., ND.h)

ND = DataFrame(h=25., bhd=25., b=0., c2=0., d03=20.7, h03=7.5, d065=14.2, h065=16.25)
p = getParam3.(ND.bhd, ND.h, ND.b, ND.c2, ND.bhd, ND.d03, ND.d065, 1.3, ND.h03, ND.h065)
p3 = plot()
plot_tree.(p)
scatter!(p3, [ND.bhd[1], ND.d03[1], ND.d065[1]] * [-.5, .5,]', [1.3, ND.h03[1], ND.h065[1]], markersize=5, label=nothing, xlabel="", xticks=(-10:10:10, abs.(-10:10:10)), showaxis=:x, yaxis=nothing, legend=(-.2,.85), fg_legend = :false, bg_legend = :false, title="3 Diameters", color=:black)
ND.v3 = getVolume2.(p, 0., ND.h)

for p in [p1, p2, p2b, p3]
  hline!(p, 0:5:25, color=:gray, label=nothing, alpha=.2, z_order=:back)
end
plot(p1, p2, p2b, p3, layout = grid(1, 4, link=:y))
#savefig("nMesspunkte.pdf")

using Plots.PlotMeasures
meinLayout = @layout [grid(1, 4); label_zeile{0.0001h}]
pLabel = plot(framestyle=:none, ticks=false, annotation=(0.5, 0.5, text("Radius [cm]", :center, 11)))
plot(p1, p2, p2b, p3, pLabel, layout = meinLayout, link = :y)
savefig("nMesspunkte.pdf")


ND = DataFrame(h=20., bhd=30., b=.5, c2=.25, d2=7., h2=18.5)
ND.dRef = ND.bhd .+ 1.3
ND.c2 = predict(Mc21d, ND)
ND.b = predict(Mb1d, ND)
p = p2cm.(getParam.(ND.bhd, ND.h, ND.b, ND.c2))
getVolume2.(p, 0., ND.h)
p = getParam2.(ND.bhd, ND.h, ND.b, ND.c2, ND.bhd, ND.d2, 1.3, ND.h2)


ND = DataFrame(h=25., bhd=24:26., b=.5, c2=.25, h03=7.5, h065=16.25, h13=1.3)
ND.dRef = ND.bhd .+ 1.3
ND.c2 = predict(Mc21d, ND)
ND.b = predict(Mb1d, ND)
p = p2cm.(getParam.(ND.bhd, ND.h, ND.b, ND.c2))
ND.d03 .= getDiameterAtHeight(p[2], ND.h03[1])
ND.d065 .= getDiameterAtHeight(p[2], ND.h065[1])
p1 = plot()
plot_tree.(p, color=1)
scatter!(p1, ND.bhd * [-.5, .5]', [1.3, 1.3], markersize=3, label=nothing, xlabel="", ylabel="Height [m]", xticks=(-10:10:10, abs.(-10:10:10)), title="1 Diameter", grid=:x, color=:black)
ND.v1 = getVolume2.(p, 0., ND.h)

p = getParam2.(ND.bhd, ND.h, ND.b, ND.c2, ND.bhd, ND.d03, 1.3, ND.h03)
p2 = plot()
plot_tree.(p, color=2)
#scatter!(p2, [ND.bhd; ND.d03[1]] * [-.5, .5,]', [ND.h13; ND.h03[1]], markersize=3, label=nothing, xlabel="Radius [cm]", xticks=(-10:10:10, abs.(-10:10:10)), showaxis=:x, yaxis=nothing, title="2 Diameters", color=:black)
scatter!(p2, [ND.bhd; ND.d03[1]] * [-.5, .5,]', [ND.h13; ND.h03[1]], markersize=3, label=nothing, xlabel="", xticks=(-10:10:10, abs.(-10:10:10)), showaxis=:x, yaxis=nothing, title="2 Diameters", color=:black)
ND.v2 = getVolume2.(p, 0., ND.h)

p = getParam2.(ND.bhd, ND.h, ND.b, ND.c2, ND.bhd, ND.d065, 1.3, ND.h065)
p2b = plot()
plot_tree.(p, color=4)
scatter!(p2b, [ND.bhd; ND.d065[1]] * [-.5, .5,]', [ND.h13; ND.h065[1]], markersize=3, label=nothing, xlabel="", xticks=(-10:10:10, abs.(-10:10:10)), showaxis=:x, yaxis=nothing, title="2 Diameters", color=:black)
ND.v2b = getVolume2.(p, 0., ND.h)

p = getParam3.(ND.bhd, ND.h, ND.b, ND.c2, ND.bhd, ND.d03, ND.d065, 1.3, ND.h03, ND.h065)
p3 = plot()
plot_tree.(p, color=3)
scatter!(p3, [ND.bhd; ND.d03[1]; ND.d065[1]] * [-.5, .5,]', [ND.h13; ND.h03[1]; ND.h065[1]], markersize=3, label=nothing, xlabel="", xticks=(-10:10:10, abs.(-10:10:10)), showaxis=:x, yaxis=nothing, legend=(-.2,.85), fg_legend = :false, bg_legend = :false, title="3 Diameters", color=:black)
ND.v3 = getVolume2.(p, 0., ND.h)

xlims!(p3, xlims(p1))
for p in [p1, p2, p2b, p3]
  hline!(p, 0:5:25, color=:gray, label=nothing, alpha=.2, z_order=:back)
end
plot(p1, p2, p3, layout = grid(1, 3, link=:y))
#savefig("dbhMeasurementError.pdf")

meinLayout = @layout [grid(1, 4); label_zeile{0.0001h}]
pLabel = plot(framestyle=:none, ticks=false, annotation=(0.5, 0.5, text("Radius [cm]", :center, 11)))
plot(p1, p2, p2b, p3, pLabel, layout = meinLayout, link = :y)
savefig("dbhMeasurementError.pdf")

ND = DataFrame(h=25., bhd=24:.1:26., b=.5, c2=.25, h03=7.5, h065=16.25, h13=1.3)
ND.dRef = ND.bhd .+ 1.3
ND.c2 = predict(Mc21d, ND)
ND.b = predict(Mb1d, ND)
p = p2cm.(getParam.(ND.bhd, ND.h, ND.b, ND.c2))
i = findfirst(ND.bhd .>= 25.)
ND.d03 .= getDiameterAtHeight(p[i], ND.h03[1])
ND.d065 .= getDiameterAtHeight(p[i], ND.h065[1])
ND.v1 = getVolume2.(p, 0., ND.h)
p = getParam2.(ND.bhd, ND.h, ND.b, ND.c2, ND.bhd, ND.d03, 1.3, ND.h03)
ND.v2 = getVolume2.(p, 0., ND.h)
p = getParam2.(ND.bhd, ND.h, ND.b, ND.c2, ND.bhd, ND.d065, 1.3, ND.h065)
ND.v2b = getVolume2.(p, 0., ND.h)
p = getParam3.(ND.bhd, ND.h, ND.b, ND.c2, ND.bhd, ND.d03, ND.d065, 1.3, ND.h03, ND.h065)
ND.v3 = getVolume2.(p, 0., ND.h)

plot(ND.bhd, ND.v1, color=1, label=nothing, xlabel="DBH [cm]", ylabel="Volume [m³]")
plot!(ND.bhd, ND.v2, color=2, label=nothing)
plot!(ND.bhd, ND.v2b, color=4, label=nothing)
plot!(ND.bhd, ND.v3, color=3, label=nothing)

plot!(Shape(ND.bhd, ND.v1), fillalpha=0.15, linecolor=:transparent, label=nothing, fillcolor=1)
plot!(Shape(ND.bhd, ND.v2), fillalpha=0.15, linecolor=:transparent, label=nothing, fillcolor=2)
plot!(Shape(ND.bhd, ND.v2b), fillalpha=0.15, linecolor=:transparent, label=nothing, fillcolor=4)
plot!(Shape(ND.bhd, ND.v3), fillalpha=0.15, linecolor=:transparent, label=nothing, fillcolor=3)
ylims!((0.538, 0.743))
savefig("dbhMeasurementErrorV.pdf")
println(ND[1:4:end, [:bhd, :v1, :v2, :v2b, :v3]])
println(ND[:, [:bhd, :v3]])


ND = DataFrame(h=25., bhd=25., b=.5, c2=.25, h03=7.5, h065=16.25, h13=1.3, d03=1:3.)
ND.dRef = ND.bhd .+ 1.3
ND.c2 = predict(Mc21d, ND)
ND.b = predict(Mb1d, ND)
p = p2cm.(getParam.(ND.bhd, ND.h, ND.b, ND.c2))
ND.d03 = getDiameterAtHeight(p[2], ND.h03[1]) .+ [-1., 0., 1.]
ND.d065 .= getDiameterAtHeight(p[2], ND.h065[1])
ND.d065b .= ND.d065 .+ [-1., 0., 1.]

p1 = plot()
plot_tree.(p, color=1)
scatter!(p1, ND.bhd[1] * [-.5, .5], [1.3, 1.3], markersize=3, label=nothing, xlabel="", ylabel="Height [m]", xticks=(-10:10:10, abs.(-10:10:10)), title="1 Diameter", grid=:x, color=:black)
ND.v1 = getVolume2.(p, 0., ND.h)

p = getParam2.(ND.bhd, ND.h, ND.b, ND.c2, ND.bhd, ND.d03, 1.3, ND.h03)
p2 = plot()
plot_tree.(p, color=2)
scatter!(p2, [ND.bhd[1]; ND.d03] * [-.5, .5,]', [ND.h13[1]; ND.h03], markersize=3, label=nothing, xlabel="", xticks=(-10:10:10, abs.(-10:10:10)), showaxis=:x, yaxis=nothing, title="2 Diameters", color=:black)
ND.v2 = getVolume2.(p, 0., ND.h)

p = getParam2.(ND.bhd, ND.h, ND.b, ND.c2, ND.bhd, ND.d065b, 1.3, ND.h065)
p2b = plot()
plot_tree.(p, color=4)
scatter!(p2b, [ND.bhd[1]; ND.d065b] * [-.5, .5,]', [ND.h13[1]; ND.h065], markersize=3, label=nothing, xlabel="", xticks=(-10:10:10, abs.(-10:10:10)), showaxis=:x, yaxis=nothing, title="2 Diameters", color=:black)
ND.v2b = getVolume2.(p, 0., ND.h)

p = getParam3.(ND.bhd, ND.h, ND.b, ND.c2, ND.bhd, ND.d03, ND.d065, 1.3, ND.h03, ND.h065)
p3 = plot()
plot_tree.(p, color=3)
scatter!(p3, [ND.bhd[1]; ND.d03; ND.d065[1]] * [-.5, .5,]', [ND.h13[1]; ND.h03; ND.h065[1]], markersize=3, label=nothing, xlabel="", xticks=(-10:10:10, abs.(-10:10:10)), showaxis=:x, yaxis=nothing, legend=(-.2,.85), fg_legend = :false, bg_legend = :false, title="3 Diameters", color=:black)
ND.v3 = getVolume2.(p, 0., ND.h)

xlims!(p3, xlims(p1))
for p in [p1, p2, p2b, p3]
  hline!(p, 0:5:25, color=:gray, label=nothing, alpha=.2, z_order=:back)
end
plot(p1, p2, p3, layout = grid(1, 3, link=:y))
#savefig("d03MeasurementError.pdf")

meinLayout = @layout [grid(1, 4); label_zeile{0.0001h}]
pLabel = plot(framestyle=:none, ticks=false, annotation=(0.5, 0.5, text("Radius [cm]", :center, 11)))
plot(p1, p2, p2b, p3, pLabel, layout = meinLayout, link = :y)
savefig("d03MeasurementError.pdf")

ND = DataFrame(h=25., bhd=25., b=.5, c2=.25, h03=7.5, h065=16.25, h13=1.3, d03=1:.1:3.)
ND.dRef = ND.bhd .+ 1.3
ND.c2 = predict(Mc21d, ND)
ND.b = predict(Mb1d, ND)
p = p2cm.(getParam.(ND.bhd, ND.h, ND.b, ND.c2))
i = findfirst(ND.bhd .>= 25.)
ND.d03 = getDiameterAtHeight(p[i], ND.h03[1]) .+ (-1.:.1:1.)
ND.d065 .= getDiameterAtHeight(p[i], ND.h065[1])
ND.d065b .= ND.d065 .+ (-1.:.1:1.)
ND.v1 = getVolume2.(p, 0., ND.h)
p = getParam2.(ND.bhd, ND.h, ND.b, ND.c2, ND.bhd, ND.d03, 1.3, ND.h03)
ND.v2 = getVolume2.(p, 0., ND.h)
p = getParam2.(ND.bhd, ND.h, ND.b, ND.c2, ND.bhd, ND.d065b, 1.3, ND.h065)
ND.v2b = getVolume2.(p, 0., ND.h)
p = getParam3.(ND.bhd, ND.h, ND.b, ND.c2, ND.bhd, ND.d03, ND.d065, 1.3, ND.h03, ND.h065)
ND.v3 = getVolume2.(p, 0., ND.h)

plot(ND.d03, ND.v1, color=1, label=nothing, xlabel=L"d_{0.3h} [cm]", ylabel="Volume [m³]")
plot!(ND.d03, ND.v2, color=2, label=nothing)
plot!(ND.d03, ND.v2b, color=4, label=nothing)
plot!(ND.d03, ND.v3, color=3, label=nothing)

plot!(Shape(ND.d03, ND.v1), fillalpha=0.15, linecolor=:transparent, label=nothing, fillcolor=1)
plot!(Shape(ND.d03, ND.v2), fillalpha=0.15, linecolor=:transparent, label=nothing, fillcolor=2)
plot!(Shape(ND.d03, ND.v2b), fillalpha=0.15, linecolor=:transparent, label=nothing, fillcolor=4)
plot!(Shape(ND.d03, ND.v3), fillalpha=0.15, linecolor=:transparent, label=nothing, fillcolor=3)
ylims!((0.538, 0.743))
savefig("d03MeasurementErrorV.pdf")
println(ND[1:2:end, [:bhd, :d03, :d065, :d065b, :v1, :v2, :v2b, :v3]])
println(ND[:, [:bhd, :d03, :d065, :v3]])


ND = DataFrame(h=25., bhd=24:26., b=.5, c2=.25, h03=7.5, h065=16.25, h13=1.3, d03=1.)
ND.dRef = ND.bhd .+ 1.3
ND.c2 = predict(Mc21d, ND)
ND.b = predict(Mb1d, ND)
p = p2cm.(getParam.(ND.bhd, ND.h, ND.b, ND.c2))
ND.d03 = getDiameterAtHeight(p[2], ND.h03[2]) .+ [-1., 0., 1.]
ND.d065 .= getDiameterAtHeight(p[2], ND.h065[2])
ND.d065b .= ND.d065 .+ [-1., 0., 1.]

p1 = plot()
plot_tree.(p, color=1)
scatter!(p1, ND.bhd * [-.5, .5]', [1.3, 1.3], markersize=3, label=nothing, xlabel="", ylabel="Height [m]", xticks=(-10:10:10, abs.(-10:10:10)), title="1 Diameter", grid=:x, color=:black)
ND.v1 = getVolume2.(p, 0., ND.h)

p = getParam2.(ND.bhd, ND.h, ND.b, ND.c2, ND.bhd, ND.d03, 1.3, ND.h03)
p2 = plot()
plot_tree.(p, color=2)
scatter!(p2, [ND.bhd; ND.d03] * [-.5, .5,]', [ND.h13; ND.h03], markersize=3, label=nothing, xlabel="", xticks=(-10:10:10, abs.(-10:10:10)), showaxis=:x, yaxis=nothing, title="2 Diameters", color=:black)
ND.v2 = getVolume2.(p, 0., ND.h)

p = getParam2.(ND.bhd, ND.h, ND.b, ND.c2, ND.bhd, ND.d065b, 1.3, ND.h065)
p2b = plot()
plot_tree.(p, color=4)
scatter!(p2b, [ND.bhd; ND.d065b] * [-.5, .5,]', [ND.h13; ND.h065], markersize=3, label=nothing, xlabel="", xticks=(-10:10:10, abs.(-10:10:10)), showaxis=:x, yaxis=nothing, title="2 Diameters", color=:black)
ND.v2b = getVolume2.(p, 0., ND.h)

p = getParam3.(ND.bhd, ND.h, ND.b, ND.c2, ND.bhd, ND.d03, ND.d065, 1.3, ND.h03, ND.h065)
p3 = plot()
plot_tree.(p, color=3)
scatter!(p3, [ND.bhd; ND.d03; ND.d065[1]] * [-.5, .5,]', [ND.h13; ND.h03; ND.h065[1]], markersize=3, label=nothing, xlabel="", xticks=(-10:10:10, abs.(-10:10:10)), showaxis=:x, yaxis=nothing, legend=(-.2,.85), fg_legend = :false, bg_legend = :false, title="3 Diameters", color=:black)
ND.v3 = getVolume2.(p, 0., ND.h)

xlims!(p3, xlims(p1))
for p in [p1, p2, p2b, p3]
  hline!(p, 0:5:25, color=:gray, label=nothing, alpha=.2, z_order=:back)
end
plot(p1, p2, p3, layout = grid(1, 3, link=:y))
#savefig("bhdUd03MeasurementErrorG.pdf")

meinLayout = @layout [grid(1, 4); label_zeile{0.0001h}]
pLabel = plot(framestyle=:none, ticks=false, annotation=(0.5, 0.5, text("Radius [cm]", :center, 11)))
plot(p1, p2, p2b, p3, pLabel, layout = meinLayout, link = :y)
savefig("bhdUd03MeasurementErrorG.pdf")

ND = DataFrame(h=25., bhd=24:.1:26., b=.5, c2=.25, h03=7.5, h065=16.25, h13=1.3, d03=1.)
ND.dRef = ND.bhd .+ 1.3
ND.c2 = predict(Mc21d, ND)
ND.b = predict(Mb1d, ND)
p = p2cm.(getParam.(ND.bhd, ND.h, ND.b, ND.c2))
i = findfirst(ND.bhd .>= 25.)
ND.d03 = getDiameterAtHeight(p[i], ND.h03[i]) .+ (-1.:.1:1.)
ND.d065 .= getDiameterAtHeight(p[i], ND.h065[i])
ND.d065b .= ND.d065 .+ (-1.:.1:1.)
ND.v1 = getVolume2.(p, 0., ND.h)
p = getParam2.(ND.bhd, ND.h, ND.b, ND.c2, ND.bhd, ND.d03, 1.3, ND.h03)
ND.v2 = getVolume2.(p, 0., ND.h)
p = getParam2.(ND.bhd, ND.h, ND.b, ND.c2, ND.bhd, ND.d065b, 1.3, ND.h065)
ND.v2b = getVolume2.(p, 0., ND.h)
p = getParam3.(ND.bhd, ND.h, ND.b, ND.c2, ND.bhd, ND.d03, ND.d065, 1.3, ND.h03, ND.h065)
ND.v3 = getVolume2.(p, 0., ND.h)

p = plot(ND.bhd, ND.v1, color=1, label=nothing, xlabel="dbh [cm]", ylabel="Volume [m³]")
vline!(xticks(p)[1][1], color=:black, lw=0.5, label=nothing, alpha=0.15)
plot!(ND.bhd, ND.v2, color=2, label=nothing)
plot!(ND.bhd, ND.v2b, color=4, label=nothing)
plot!(ND.bhd, ND.v3, color=3, label=nothing)

plot!(Shape(ND.bhd, ND.v1), fillalpha=0.15, linecolor=:transparent, label=nothing, fillcolor=1)
plot!(Shape(ND.bhd, ND.v2), fillalpha=0.15, linecolor=:transparent, label=nothing, fillcolor=2)
plot!(Shape(ND.bhd, ND.v2b), fillalpha=0.15, linecolor=:transparent, label=nothing, fillcolor=4)
plot!(Shape(ND.bhd, ND.v3), fillalpha=0.15, linecolor=:transparent, label=nothing, fillcolor=3)

plot!(twiny(), [], [], xlabel=L"d_{0.3h} [cm]", xlims=extrema(ND.d03), widen=true, showaxis=:x, label=nothing)
ylims!((0.538, 0.743))
savefig("bhdd03MeasurementErrorGV.pdf")
println(ND[1:2:end, [:bhd, :d03, :v1, :v2, :v2b, :v3]])
println(ND[:, [:bhd, :d03, :d065, :v3]])


ND = DataFrame(h=25., bhd=26:-1:24., b=.5, c2=.25, h03=7.5, h065=16.25, h13=1.3, d03=1.)
ND.dRef = ND.bhd .+ 1.3
ND.c2 = predict(Mc21d, ND)
ND.b = predict(Mb1d, ND)
p = p2cm.(getParam.(ND.bhd, ND.h, ND.b, ND.c2))
ND.d03 = getDiameterAtHeight(p[2], ND.h03[2]) .+ [-1., 0., 1.]
ND.d065 .= getDiameterAtHeight(p[2], ND.h065[2])
ND.d065b .= ND.d065 .+ [-1., 0., 1.]

p1 = plot()
plot_tree.(p, color=1)
scatter!(p1, ND.bhd * [-.5, .5]', [1.3, 1.3], markersize=3, label=nothing, xlabel="", ylabel="Height [m]", xticks=(-10:10:10, abs.(-10:10:10)), title="1 Diameter", grid=:x, color=:black)
ND.v1 = getVolume2.(p, 0., ND.h)

p = getParam2.(ND.bhd, ND.h, ND.b, ND.c2, ND.bhd, ND.d03, 1.3, ND.h03)
p2 = plot()
plot_tree.(p, color=2)
scatter!(p2, [ND.bhd; ND.d03] * [-.5, .5,]', [ND.h13; ND.h03], markersize=3, label=nothing, xlabel="", xticks=(-10:10:10, abs.(-10:10:10)), showaxis=:x, yaxis=nothing, title="2 Diameters", color=:black)
ND.v2 = getVolume2.(p, 0., ND.h)

p = getParam2.(ND.bhd, ND.h, ND.b, ND.c2, ND.bhd, ND.d065b, 1.3, ND.h065)
p2b = plot()
plot_tree.(p, color=4)
scatter!(p2b, [ND.bhd; ND.d065b] * [-.5, .5,]', [ND.h13; ND.h065], markersize=3, label=nothing, xlabel="", xticks=(-10:10:10, abs.(-10:10:10)), showaxis=:x, yaxis=nothing, title="2 Diameters", color=:black)
ND.v2b = getVolume2.(p, 0., ND.h)

p = getParam3.(ND.bhd, ND.h, ND.b, ND.c2, ND.bhd, ND.d03, ND.d065, 1.3, ND.h03, ND.h065)
p3 = plot()
plot_tree.(p, color=3)
scatter!(p3, [ND.bhd; ND.d03; ND.d065[1]] * [-.5, .5,]', [ND.h13; ND.h03; ND.h065[1]], markersize=3, label=nothing, xlabel="", xticks=(-10:10:10, abs.(-10:10:10)), showaxis=:x, yaxis=nothing, legend=(-.2,.85), fg_legend = :false, bg_legend = :false, title="3 Diameters", color=:black)
ND.v3 = getVolume2.(p, 0., ND.h)

xlims!(p3, xlims(p1))
for p in [p1, p2, p2b, p3]
  hline!(p, 0:5:25, color=:gray, label=nothing, alpha=.2, z_order=:back)
end
plot(p1, p2, p3, layout = grid(1, 3, link=:y))
#savefig("bhdUd03MeasurementError.pdf")

meinLayout = @layout [grid(1, 4); label_zeile{0.0001h}]
pLabel = plot(framestyle=:none, ticks=false, annotation=(0.5, 0.5, text("Radius [cm]", :center, 11)))
plot(p1, p2, p2b, p3, pLabel, layout = meinLayout, link = :y)
savefig("bhdUd03MeasurementError.pdf")

ND = DataFrame(h=25., bhd=26:-.1:24., b=.5, c2=.25, h03=7.5, h065=16.25, h13=1.3, d03=1.)
ND.dRef = ND.bhd .+ 1.3
ND.c2 = predict(Mc21d, ND)
ND.b = predict(Mb1d, ND)
p = p2cm.(getParam.(ND.bhd, ND.h, ND.b, ND.c2))
i = findfirst(ND.bhd .<= 25.)
ND.d03 = getDiameterAtHeight(p[i], ND.h03[i]) .+ (-1.:.1:1.)
ND.d065 .= getDiameterAtHeight(p[i], ND.h065[i])
ND.d065b .= ND.d065 .+ (-1.:.1:1.)
ND.v1 = getVolume2.(p, 0., ND.h)
p = getParam2.(ND.bhd, ND.h, ND.b, ND.c2, ND.bhd, ND.d03, 1.3, ND.h03)
ND.v2 = getVolume2.(p, 0., ND.h)
p = getParam2.(ND.bhd, ND.h, ND.b, ND.c2, ND.bhd, ND.d065b, 1.3, ND.h065)
ND.v2b = getVolume2.(p, 0., ND.h)
p = getParam3.(ND.bhd, ND.h, ND.b, ND.c2, ND.bhd, ND.d03, ND.d065, 1.3, ND.h03, ND.h065)
ND.v3 = getVolume2.(p, 0., ND.h)

p = plot(ND.bhd, ND.v1, color=1, label=nothing, xlabel="dbh [cm]", ylabel="Volume [m³]")
vline!(xticks(p)[1][1], color=:black, lw=0.5, label=nothing, alpha=0.15)
plot!(ND.bhd, ND.v2, color=2, label=nothing)
plot!(ND.bhd, ND.v2b, color=4, label=nothing)
plot!(ND.bhd, ND.v3, color=3, label=nothing)

plot!(Shape(ND.bhd, ND.v1), fillalpha=0.15, linecolor=:transparent, label=nothing, fillcolor=1)
plot!(Shape(ND.bhd, ND.v2), fillalpha=0.15, linecolor=:transparent, label=nothing, fillcolor=2)
plot!(Shape(ND.bhd, ND.v2b), fillalpha=0.15, linecolor=:transparent, label=nothing, fillcolor=4)
plot!(Shape(ND.bhd, ND.v3), fillalpha=0.15, linecolor=:transparent, label=nothing, fillcolor=3)

plot!(twiny(), [], [], xlabel=L"d_{0.3h} [cm]", xlims=extrema(ND.d03), widen=true, xflip=true, showaxis=:x, label=nothing)
ylims!((0.538, 0.743))
savefig("bhdd03MeasurementErrorV.pdf")
println(ND[1:2:end, [:bhd, :d03, :v1, :v2, :v2b, :v3]])
println(ND[:, [:bhd, :d03, :d065, :v2b, :v3]])


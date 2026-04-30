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

#Stockhoehe - Varianten
t = [0:10:100; 200; 500]
0.08 .+ t ./ 400.
0.08 .+ 0.04 .* log.(1 .+ t)
0.08 .+ 0.02 .* sqrt.(t)
0.08 .+ 0.2 .* (1 .- exp.(-.03 .* t))
0.08 .+ 0.3 .* t ./ (30 .+ t)

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
    hSt = 5. - (0.08 + 0.3 * (x[3] + 0.65) / (30 + (x[3] + 0.65)))
    model([6.0, 5.0, hSt, 3.70], coef(fit))
  else
    [missing, missing, missing, missing]
  end
end |> stack
sum(.!ismissing.(d0), dims=2)  # 2828
describe(d.DBH ./ d0[4, :])
#i = findall(coalesce.(d0[1, :] ./ d.DBH .> 3.0, false))  # 412
i = findall(coalesce.(d0[1, :] ./ d0[4, :] .> 3.0, false))  # 410
d0[1, i] .= missing
i = findall(coalesce.(d0[2, :] ./ d0[4, :] .> 2.0, false))  # 47
d0[2, i] .= missing
scatter(x[:,3], d0[3, :])
i = findall(coalesce.(d0[3, :] ./ d0[4, :] .> 1.8, false))  # 11
d0[3, i] .= missing
d00 = copy(d0)


#Hilfsfunktion um Residuen zu plotten
using Loess
function rp(y, x, sp=0.1; h=0., xlab="", ylab="", ms=2, ma=0.2, col=:black, ax=:xy, trim=0.02, trx=:identity)
  scatter(x, y, color=:black, ms=ms, ma=ma, xlabel=xlab, ylabel=ylab, showaxis=ax, xscale=trx)
  le = loess(x, y; span=sp)
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
D = DataFrame(dbh=d00[4, :], d0=d00[2, :], h=ht.H_total, dAux = d.DBH .+ 1.3 ./ 2)
dropmissing!(D)
size(D)  # 2781
D.s = (D.d0 .- D.dbh) ./ 130 .- .015
rp(D.s, D.h, col=:red)
rp(D.s, D.dAux, col=:red)
sum(D.s .< 0)  # 223
sum(D.s .< 0.001)  # 238
D.s[D.s .< 0.001] .= 0.001
reg = glm(@formula(s ~ inv(h) + inv(dAux)), D, Gamma(), LogLink())
#                  Coef.  Std. Error       z  Pr(>|z|)   Lower 95%   Upper 95%
#(Intercept)   -0.274453   0.0431239   -6.36    <1e-09   -0.358974   -0.189931
#inv(h)       -37.6638     2.14314    -17.57    <1e-68  -41.8643    -33.4633
#inv(dAux)    -34.7816     2.31706    -15.01    <1e-50  -39.323     -30.2403
#i = ismissing.(d0[2, :])
i = fill(true, size(d)[1])
x = predict(reg, DataFrame(dAux=d[i, :DBH] .+ 1.3 ./ 2, h=ht[i, :H_total])) .+ .015
d0[2, i] .= d[i, :DBH] .+ x .* 130
Plots.scalefontsizes(1.4)
rp(response(reg) .- predict(reg), predict(reg), col=:red, ylab=L"s_{0} - \hat{s}_{0}", xlab=L"\hat{s}_{0}")
savefig("d0ResidPred.png")
rp(response(reg) .- predict(reg), D.h, col=:red, ylab=L"s_{0} - \hat{s}_{0}", xlab="Height [m]")
savefig("d0ResidHeight.png")
rp(response(reg) .- predict(reg), D.dAux, col=:red, ylab=L"s_{0} - \hat{s}_{0}", xlab=L"d_a [cm]")
savefig("d0ResidDAux.png")
Plots.scalefontsizes()
predict(reg, DataFrame(h=[0, 1.3, 30, 70], dAux=[0, 1.45, 45, 200])) .+ .015
plot((predict(reg, DataFrame(h=90, dAux=0:200)) .+ .015) .* 130)
# d bei h-1
D = DataFrame(dbh=d00[4, :], d0=d00[2, :], dm1=d00[1, :], h=ht.H_total, dAux = d00[4, :] .+ 1.3 ./ 2)
dropmissing!(D)
size(D)  # 2418
D.s = (D.dm1 .- D.dbh) ./ 230 .- (D.d0 .- D.dbh) ./ 130
rp(D.s, D.h, col=:red)
rp(D.s, D.dAux, col=:red)
sum(D.s .< 0)  # 86
sum(D.s .< 0.001)  # 250
D.s[D.s .< 0.001] .= 0.001
reg = glm(@formula(s ~ inv(h) + inv(dAux)), D, Gamma(), LogLink())
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
# d am Stock
D = DataFrame(dbh=d00[4, :], d0=d00[2, :], dSt=d00[3, :], h=ht.H_total, dAux = d00[4, :] .+ 1.3 ./ 2)
dropmissing!(D)
size(D)  # 2781
D.s = (D.d0 .- D.dbh) ./ 130
D.hSt = 0.08 .+ 0.3 .* D.dAux ./ (30 .+ D.dAux)
D.rHSt = D.hSt ./ 1.3
D.r = (D.d0 .- D.dSt) ./ (D.d0 - D.dbh)
D.e = log.(D.r) ./ log.(D.rHSt)
rp(D.e, D.s, col=:red)
reg = glm(@formula(e ~ log(0.01 + s)), D, Gamma(), LogLink())
rp(response(reg) .- predict(reg), predict(reg), col=:red)
rp(response(reg) .- predict(reg), D.s, col=:red)
@. model(x, p) = exp(p[1] + p[2] * log(p[3] + x))
fit = curve_fit(model, D.s, D.e, [-.57, -.15, 0.01])
#i = ismissing.(d0[3, :])
i = fill(true, size(d)[1])
x = model((d0[2,i] - d[i,:DBH]) ./ 130 , coef(fit))
y = d[i,:DBH] .+ 1.3 ./ 2
y = 0.08 .+ 0.3 .* y ./ (30 .+ y)
y ./= 1.3
d0[3,i] .= d0[2,i] .- y.^x .* (d0[2,i] .- d[i,:DBH])
quick_summary(fit)
# Row │ Parameter  Estimate    StdError    t_stat     p_val        Lower95    Upper95    
#   1 │ p[1]       -0.629332   0.0138572   -45.4156   0.0          -0.656503  -0.60216
#   2 │ p[2]       -0.184223   0.00940268  -19.5926   0.0          -0.20266   -0.165786
#   3 │ p[3]        0.0259615  0.00390674    6.64529  3.62823e-11   0.018301   0.0336219
coef(fit)
stderror(fit)
confint(fit; level=0.95)
Plots.scalefontsizes(1.4)
rp(-fit.resid, D.e + fit.resid, col=:red, ylab=L"e_{st} - \hat{e}_{st}", xlab=L"\hat{e}_{st}")
savefig("dEStResidPred.png")
rp(-fit.resid, D.s, col=:red, ylab=L"e_{st} - \hat{e}_{st}", xlab=L"d_0 - d_{1.3}/130")
savefig("dEStResidSlp.png")
Plots.scalefontsizes()
model(0:0.1:1, coef(fit))
plot((0:0.1:1).^.5)
plot((0:0.1:1).^1.05)
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
scatter(vec((x[1:end-1,:]) ./ (0.65 .+ x[2,:])'), vec(y[1:end-1,:] ./ y[end,:]'), xlabel=L"d / d_a", ylabel=L"h / h_{max}", color=:black, ms=1, ma=.1, legend=false, fmt = :png)
savefig("relativeTaper.png")
Plots.scalefontsizes()

res = map(eachcol(x), eachcol(y)) do xd, yh
  i = findall(.!isnan.(xd))
  D = combine(groupby(DataFrame(h=yh[i], d=xd[i]), :h, sort=true), :d => mean => :d)
  dh = diff(D.h)
  s = diff(D.d) ./ 50.0 ./ dh  # Neigung bezogen auf Markroehre
  rhs = (D.h[1:end-1] + dh ./ 2) ./ D.h[end]
  ds = diff(s)
  mdh = minimum([dh[1:end-1] dh[2:end]], dims=2)
  rhd = D.h[2:end-1] ./ D.h[end]
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
Plots.scalefontsizes(1.1)
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
scatter(res[:, 4], res[:, 6], xlabel="Slope Change", ylabel="Relative Height", color=:black, ms=1, ma=0.1, label=nothing, xlim=[-0.1, 0.25], ylim=[0, 1], widen=true, fmt=:png)
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
  ds = diff(s)
  rhd = D.h[2:end-1] ./ D.h[end]
  i = findall(.!(ra[2, 1].(rhd) .<= ds .<= ra[2, 2].(rhd))) .+ 1
  s2(D, unique([i .- 1; i]))
end
#
#d0 und d-1 dazuhaengen
res = map(res, eachrow(Matrix(d[:, 1:2]))) do x, y
  combine(groupby(vcat(DataFrame(h=[-1.0, 0.0], d=y), x), :h, sort=true), :d => mean => :d)
end

using Integrals, Optim
v = map(res) do x
  s = interpolate(x.h, x.d, SteffenMonotonicInterpolation())
  v = solve(IntegralProblem((h, p) -> s(h)^2 * pi / 40000.0, (0.0, max(0.0, x.h[end]))), QuadGKJL()).u
  da = s(1.3) + 0.65 
  hSt = min(x.h[end], 0.08 + 0.3 * da / (30 + da))
  vSt = solve(IntegralProblem((h, p) -> s(h)^2 * pi / 40000.0, (0.0, hSt)), QuadGKJL()).u
  hDh = 0.
  vWi = if s(hSt) > 7.0
    hDh = optimize(h -> abs(7.0 - s(h[1])), hSt, x.h[end]).minimizer[1]
    solve(IntegralProblem((h, p) -> s(h)^2 * pi / 40000.0, (hDh, x.h[end])), QuadGKJL()).u
  else
    v - vSt
  end
  [s(1.3); x.h[end]; v; vSt; vWi; hSt; s(0.); s(hSt); hDh; da; s.(min.((0.:.05:1)*x.h[end], x.h[end]))]
end |> stack

#Hoehe Derbholzgrenze
D = DataFrame(h=v[2,:], hSt=v[6,:], hD7=v[9,:], dhSt=v[8,:], da=v[9,:])
filter!(:hD7 => >(0.), D)
filter!(:hD7 => >(1.), D)  #
D.hs = D.h .- D.hSt
D.r = (D.hD7 .- D.hSt) ./ D.hs
D.dhStD .= D.dhSt .- 7
filter!(:dhStD => >(0.), D)
filter!(:dhStD => >(1.), D)  #
D.rs = 1 .- (7 .- .8) ./ (D.dhSt .- .8)
scatter(D.r, D.rs)
D.e = log.(D.r) ./ log.(D.rs)
rp(D.e, D.dhStD, col=:red)
rp(D.e, D.hs ./ D.dhSt, col=:red)
@. model(x, p) = exp(p[1] + p[2] * log(p[3] * x[:,1])^2 + p[4] * log(x[:,2]))
fit = curve_fit(model, [D.dhStD (D.hs ./ D.dhSt)], D.e, [-.6, 0.1, .04, 0.3])
quick_summary(fit)
# Row │ Parameter  Estimate    StdError     t_stat     p_val    Lower95     Upper95    
#   1 │ p[1]       -0.580281   0.00285244   -203.433       0.0  -0.585872   -0.57469
#   2 │ p[2]        0.0909597  0.001664       54.6632      0.0   0.0876981   0.0942213
#   3 │ p[3]        0.039944   0.000747524    53.4351      0.0   0.0384788   0.0414092
#   4 │ p[4]       -0.328553   0.00778766    -42.189       0.0  -0.343818   -0.313288
Plots.scalefontsizes(1.4)
#rp(-fit.resid, D.e .+ fit.resid, col=:red, trim=.015, ylab=L"e_{d7} - \hat{e_{d7}}", xlab=L"\hat{e_{d7}}")
rp(-fit.resid, D.e .+ fit.resid, col=:red, ylab=L"Residuals \ e_{d7} - \hat{e}_{d7}", xlab=L"Predicted \ \hat{e}_{d7}")
savefig("hd7ResidPred.png")
#rp(-fit.resid, D.dhStD, col=:red, ylab=L"e - \hat{e}", xlab=L"d_{merc,hst}")
rp(-fit.resid, D.dhStD, col=:red, ylab=L"Residuals \ e_{d7} - \hat{e}_{d7}", xlab=L"d_{st} - 7 \ [cm]")
savefig("hd7ResidDhstD.png")
#rp(-fit.resid, D.hs ./ D.dhSt, col=:red, ylab=L"e - \hat{e}", xlab=L"(h_{Total} - h_{St}) / d_{hst}")
rp(-fit.resid, D.hs ./ D.dhSt, col=:red, ylab=L"e_{d7} - \hat{e}_{d7}", xlab=L"Slenderness \ (H_{total} - h_{st}) / d_{st}")
savefig("hd7ResidHD.png")
Plots.scalefontsizes()

#Volumen Schaftderbholz
D = DataFrame(vSD = v[3,:] .- v[4,:] .- v[5,:], dhSt=v[8,:], hSt=v[6,:], da=v[1,:], h=v[2,:], hD7=v[9,:])
filter!(:vSD => >(0.), D)
D.lSD = D.hD7 .- D.hSt  # Laenge Schaftderbholz
D.dnd = min.(7., D.dhSt)  # Nicht Derbholzduchmesser
D.vKeSD = D.lSD .* pi ./ 3. .* ((D.dhSt./200.).^2 .+ (D.dhSt ./ 200.) .* (D.dnd ./ 200.) .+ (D.dnd ./ 200.).^2)  # Kegelstumpf
filter!(:vKeSD => >(0.), D)
D.fz = D.vKeSD ./ D.vSD
rp(D.fz, D.h ./ D.da, col=:red)
rp(D.fz, D.h, col=:red)
rp(D.fz, D.lSD, col=:red)
rp(D.fz, D.dhSt, col=:red)
rp(D.fz, D.da, col=:red)
@. model(x, p) = p[1] + p[2] / (1. + exp(p[3] + p[4] / max(1, x)))
fit = curve_fit(model, D.lSD, D.fz, [1., 1., -2., 10.])
quick_summary(fit)
# Row │ Parameter  Estimate    StdError    t_stat    p_val    Lower95     Upper95
#   1 │ p[1]         1.05708   0.00165828  637.457       0.0    1.05383     1.06033
#   2 │ p[2]         0.443337  0.0248524    17.8388      0.0    0.394623    0.492051
#   3 │ p[3]        -4.98345   0.369294    -13.4945      0.0   -5.70731    -4.25959
#   4 │ p[4]       138.223     7.16239      19.2984      0.0  124.184     152.262
plot(model(0:50, coef(fit)))
Plots.scalefontsizes(1.4)
#rp(-fit.resid, D.fz .+ fit.resid, col=:red, ylab=L"\zeta_d - \hat{\zeta_d}", xlab=L"\hat{\zeta_d}")
rp(-fit.resid, D.fz .+ fit.resid, col=:red, ylab=L"Residuals \ \zeta_d - \hat{\zeta}_d", xlab=L"Predicted \ \hat{\zeta}_d")
savefig("fzDResidPred.png")
#rp(-fit.resid, D.lSD, col=:red, ylab=L"\zeta_d - \hat{\zeta_d}", xlab=L"l_{SD}")
rp(-fit.resid, D.lSD, col=:red, ylab=L"Residuals \ \zeta_d - \hat{\zeta}_d", xlab=L"Merchantable \ length \ l_{SD} \ [m]")
savefig("fzDResidLsd.png")
Plots.scalefontsizes()

#Volumen Wipfel
D = DataFrame(vWi = v[5,:], dhSt=v[8,:], hSt=v[6,:], da=v[1,:], h=v[2,:], hD7=v[9,:])
D.dwi .= .8 #Wipfeldurchmesser
D.dnd = min.(7., D.dhSt)  # Nicht Derbholzduchmesser
D.lWi = D.h .- D.hD7
D.vKeWi = D.lWi .* pi ./ 3. .* ((D.dwi./200.).^2 .+ (D.dwi ./ 200.) .* (D.dnd ./ 200.) .+ (D.dnd ./ 200.).^2)  # Kegelstumpf
D.fz = D.vKeWi ./ D.vWi
rp(D.fz, D.lWi, col=:red)
rp(D.fz, D.da, col=:red)
rp(D.fz, D.h, col=:red)
rp(D.fz, D.lWi ./ D.h, col=:red)
rp(D.fz, D.da ./ D.h, col=:red)
rp(D.fz, D.hD7, col=:red)

@. model(x, p) = p[1] + p[2] / (1. + exp(p[3] + p[4] * x))
fit = curve_fit(model, D.lWi, D.fz, [1., 1., 0., .1])
fit = curve_fit(model, D.da ./ D.h, D.fz, [1., 1., 0., .1])
@. model(x, p) = p[1] + p[2] / (1. + exp(p[3] + p[4] * x[:,1] + p[5] * x[:,2]))
fit = curve_fit(model, [D.lWi (D.da ./ D.h)], D.fz, [.5, 1., .4, 1.2, -8.])
@. model(x, p) = p[1] + p[2] / (1. + exp(p[3] * x[:,1] + p[4] * x[:,2]))
fit = curve_fit(model, [D.lWi (D.da ./ D.h)], D.fz, [.5, 1., 1.2, -8.])
quick_summary(fit)
#Row │ Parameter  Estimate   StdError    t_stat    p_val    Lower95    Upper95  
#  1 │ p[1]        0.866187  0.00496459  174.473       0.0   0.856455   0.875918
#  2 │ p[2]        0.114083  0.00533969   21.365       0.0   0.103616   0.124549
#  3 │ p[3]        1.27105   0.116922     10.871       0.0   1.04187    1.50023
#  4 │ p[4]       -7.76699   0.622445    -12.4782      0.0  -8.98706   -6.54693
std(fit.resid)
plot(model([0:10 repeat([1.], 11)], coef(fit)))
plot(model([repeat([8.], 21) 0:.1:2.], coef(fit)))
Plots.scalefontsizes(1.4)
#rp(-fit.resid, D.fz .+ fit.resid, .5, col=:red, ylab=L"fz_t - \hat{fz_t}", xlab=L"\hat{fz_t}")
rp(-fit.resid, D.fz .+ fit.resid, .4, col=:red, ylab=L"Residuals \ \zeta_t - \hat{\zeta}_t", xlab=L"Predicted \ \hat{\zeta}_t")
savefig("fzTResidPred.png")
#rp(-fit.resid, D.lWi, col=:red, ylab=L"fz_t - \hat{fz_t}", xlab=L"l_{Top} [m]")
rp(-fit.resid, D.lWi, col=:red, ylab=L"Residuals \ \zeta_t - \hat{\zeta}_t", xlab=L"Top \ length \ l_{top} \ [m]")
savefig("fzTResidLTop.png")
#rp(-fit.resid, D.da ./ D.h, col=:red, ylab=L"fz_t - \hat{fz_t}", xlab=L"d_{Aux} / h_{Total}")
rp(-fit.resid, D.da ./ D.h, col=:red, ylab=L"\zeta_t - \hat{\zeta_t}", xlab=L"Relative \ diameter \ d_a / H_{total}")
savefig("fzTResidDDH.png")

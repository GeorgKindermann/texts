#Typical the following is done only once to download the data so I commented it out
#download("https://www.envidat.ch/dataset/b23a6bfa-3a55-4677-91e2-d32f91945574/resource/8ba7cb02-3d90-409f-9113-4e9aa5e7d33c/download/stembranchchwood_2024-04-25.txt", "stembranchchwood_2024-04-25.txt")
# https://doi.org/10.16904/envidat.486
# http://dx.doi.org/10.1038/s41597-024-03336-7

using DataFrames, CSV
df = CSV.read("stembranchchwood_2024-04-25.txt", DataFrame, missingstring="NA")
filter!(:TreeSpecies => sp -> !ismissing(sp) && sp == 21, df)  # Fichte

using StatsBase, Plots
summarystats(df.DBH)
describe(df.DBH)

describe(df[:,[:H_total, :DBH, :DM065, :DM1, :D_coarsestemfinal, :D_top]])
describe(df[ismissing.(df.DM1),:DBH])

Plots.scalefontsizes(1.4)
scatter(df.DBH ./ 10., df.H_total ./ 10, legend=false, xlabel="DBH [cm]", ylabel="Height [m]", color=:black, ylim=[-1, 50], ms=1, ma=0.2)
savefig("dhScatter.png")
Plots.scalefontsizes()

Plots.scalefontsizes(1.4)
i = round.(Int, range(1, size(df, 1), 100))
plot(sort(df.DBH ./ 10.)[i], i, label="DBH", xlabel="DBH [cm]", ylabel="Observations", lw=2, formatter=:plain, color=:black, xlim=[0, 100], widen=true, xticks=0:20:100)
plot!([0], [0], label="Height", lw=2, color=:black, linestyle=:dot)
plot!(twiny(), sort(df.H_total ./ 10.)[i], i, label=nothing, xlabel="Height [m]", lw=2, color=:black, linestyle=:dot, xlim=[0, 50], widen=true, xticks=0:10:50, grid=true)
savefig("dhCumObs.pdf")
Plots.scalefontsizes()

using Interpolations
#df.TreeID .== 272, 117
d = [122.9, 68.3, 50.9, 45.1, 41.7, 35.4, 34.8, 33.3, 31.2, 30.3, 30.2, 27.3, 25.7, 24.2, 22.4, 18.7, 15.5, 12.1, 8.5, 7.0, 4.4, 0.7]
h = [-1.0, 0.0, 0.65, 1.0, 1.3, 3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0, 17.0, 19.0, 21.0, 23.0, 25.0, 27.0, 28.5, 29.0, 30.5, 32.0]
#
i = [4; 6:19; 21; 22]
p1 = scatter(d[i]./2., h[i], markersize=3, label=nothing, xlabel="", ylabel="Height [m]", xticks=(-20:10:20, abs.(-20:10:20)), title="Steps", grid=:x, color=:black)
xh = [0.; repeat([h[i[1:end-3]] .+ 1; h[i[end-2]] - (h[i[end-3]] + 1) + h[i[end-2]]], inner=2); h[i[end]]]
yd = repeat(d[i[1:end-1]], inner=2) ./ 2
plot!([yd; -reverse(yd)], [xh; reverse(xh)], fillrange = [0.], fillalpha = 0.3, label=nothing, color=:black)
#
i = [2; 4; 6:19; 21; 22]
p2 = scatter(d[i]./2., h[i], markersize=3, label=nothing, xlabel="Radius [cm]", xticks=(-20:10:20, abs.(-20:10:20)), showaxis=:x, yaxis=nothing, title="Linear", color=:black)
fi = interpolate(h[i], d[i], LinearMonotonicInterpolation())
xh = h[i]
yd = fi.(xh) ./ 2.
plot!([yd; -reverse(yd)], [xh; reverse(xh)], fillrange = [0.], fillalpha = 0.3, label=nothing, color=:black)
#
i = [2; 4; 6:19; 21; 22]
p3 = scatter(d[i]./2., h[i], markersize=3, label="Measurements", xlabel="", xticks=(-20:10:20, abs.(-20:10:20)), showaxis=:x, yaxis=nothing, legend=(-.2,.85), fg_legend = :false, bg_legend = :false, title="Curve", color=:black)
fi = interpolate(h, d, SteffenMonotonicInterpolation())
xh = unique(sort([range(0, maximum(h), 100); h[i]]))
yd = fi.(xh) ./ 2.
plot!([yd; -reverse(yd)], [xh; reverse(xh)], fillrange = [0.], fillalpha = 0.3, label="Volume", color=:black)
#
xlims!(p1, xlims(p2))
for p in [p1, p2, p3]
  hline!(p, 0:10:30, color=:gray, label=nothing, alpha=.2, z_order=:back)
end
plot(p1, p2, p3, layout = grid(1, 3, link=:y))
savefig("vMethods.pdf")

h = [-1.; 0.; 0.65; 1; 1.3; 3:2:53]  # h-1, h0, h065, h1, h1.3, h3, h5, ..., h51, h53
ht = df[:, [:H_total, :L_coarsestem, :L_coarsestemfinal, :L_top]]
ht.hc = ht.L_coarsestem .- ht.L_coarsestemfinal ./ 2.0
ht.ht = ht.L_coarsestem .+ ht.L_top ./ 2.0
ht = ht[:, [:hc, :L_coarsestem, :ht, :H_total]] ./ 10.0  # hd>7, hd=7, hd<7, h

d = df[:, Cols(:DBH, r"^DM", :D_coarsestemfinal, :D_top)] ./ 10.0
replace!(d.D_coarsestemfinal, 0.0 => missing)
replace!(d.D_top, 0.0 => missing)
insertcols!(d, :D_top, :d7 => 7.0)
d.dw .= 0.7 # Annahme
d.d0 = Vector{Union{Missing, Float64}}(missing, nrow(d))
d.du1 = Vector{Union{Missing, Float64}}(missing, nrow(d))
select!(d, :du1, :d0, :DM065, :DM1, :) # d-1, d0, d065, d1, d1.3, d3, d5, ..., d51, d53, d>7, d=7, d<7, dw

print(mapcols(x -> count(!ismissing, x), d[:,3:end]))
mapcols(x -> count(>(0), x), ht)

#Durchmesser des Wipfels in Abhaengigkeit von der Baumhoehe
x = 0:0.01:0.5  # Hoehe in m
y = 0.1 .+ tanh.(x .^ 1.7 .* 15.0) .* 0.6 # Durchmesser in cm
using Plots
Plots.scalefontsizes(1.4)
plot(x, y, label=nothing, xlabel="Height [m]", ylabel="Diameter [cm]", lw=2, color=:black)
savefig("dAtTop.pdf")
Plots.scalefontsizes()

#d bei h=0
using LsqFit
@. model(h, p) = p[1] + p[2] * h^p[3]
@. model1(h, p) = p[1] + p[2] * h
x = Matrix(coalesce.(d[:,[:DM065, :DM1, :DBH, :DM3, :DM5]], NaN))
y = 5. .- [.65, 1., 1.3, 3., 5.]
d0 = map(eachrow(x)) do x
  if(all(.!isnan.(x)) && length(unique(x)) > 2 && issorted(x, rev=true))
    fit = curve_fit(model1, y, x, [x[end], 1.])
    fit = curve_fit(model, y, x, [coef(fit); 1.])
    model([6., 5.], coef(fit))
  else
    [missing, missing]
  end
end |> stack
sum(.!ismissing.(d0), dims=2)
i = findall(coalesce.(d0[1,:] ./ d.DBH .> 3., false))
d0[1,i] .= missing
i = findall(coalesce.(d0[2,:] ./ d.DBH .> 2., false))
d0[2,i] .= missing
scatter(d.DBH, d0[1,:] ./ d.DBH)
scatter(d.DBH, d0[2,:] ./ d.DBH)
using Statistics
mean(skipmissing(d0[1,:] ./ d.DBH))
median(skipmissing(d0[1,:] ./ d.DBH))
mean(skipmissing(d0[2,:] ./ d.DBH))
median(skipmissing(d0[2,:] ./ d.DBH))
#
#Hilfsfunktion um Residuen zu plotten
using Loess
#function rp(y, x, sp=0.1)
function rp(y, x, sp=0.1; h=0., xlab="", ylab="", ms=2, ma=.2, col=:black, ax=:xy)
  scatter(x, y, color=:black, ms=ms, ma=ma, xlabel=xlab, ylabel=ylab, showaxis=ax)
  le = loess(x, y; span=sp)
  us = range(extrema(x)...; step=0.01)
  vs = predict(le, us)
  hline!([h], lw=2, color=:grey,linestyle=:dash)
  plot!(us, vs, legend=false, lw=3, color=col)
end
# Duchmesser bei h=0 aus bhd und h schaetzen
using GLM
D = DataFrame(d=d.DBH .+ 1.3, d0=d0[2,:], h=ht.H_total)
dropmissing!(D)
D.d0r = D.d0 ./ D.d .- 1
filter!(:d0r => >(0.0), D)
reg = glm(@formula(d0r ~ log(h) + d), D, Normal(), LogLink())
#Coef.  Std. Error       z  Pr(>|z|)   Lower 95%    Upper 95%
#───────────────────────────────────────────────────────────────────────────────
#(Intercept)  -5.87166     0.306584    -19.15    <1e-81  -6.47255    -5.27077
#log(h)        1.48349     0.100722     14.73    <1e-48   1.28608     1.6809
#d            -0.00899125  0.00136781   -6.57    <1e-10  -0.0116721  -0.00631039
(cor(response(reg), predict(reg)), std(response(reg) - predict(reg)))
x = 0:400
plot(x, predict(reg, DataFrame(d=x, h=30)))
x = 0:60
plot(x, predict(reg, DataFrame(d=30, h=x)))
using LaTeXStrings
p1 = rp(predict(reg) - response(reg), predict(reg), ylab="Residual", xlab="Predicted "*L"d_0 / d_a - 1")
p2 = rp(predict(reg) - response(reg), D.d, ylab="Residual", xlab=L"d_a"*" [cm]")
p3 = rp(predict(reg) - response(reg), D.h, ylab="Residual", xlab="Tree Height [m]")
p4 = rp(predict(reg) - response(reg), 100 .* D.h ./ D.d, ylab="Residual", xlab=L"h/d_a")
plot(p1, p2, p3, p4, layout = grid(2, 2, link=:y))
savefig("d0Resid.png")

i = ismissing.(d0[2,:])
x = predict(reg, DataFrame(d=d[i, :DBH] .+ 1.3, h=ht[i, :H_total]))
d0[2,i] .= (d[i, :DBH] .+ 1.3) .* (1. .+ x)
#
D = DataFrame(d=d.DBH .+ 1.3, d0=d0[1,:], h=ht.H_total)
dropmissing!(D)
D.d0r = D.d0 ./ D.d .- 1
filter!(:d0r => >(0.0), D)
reg = glm(@formula(d0r ~ log(h) + d), D, Normal(), LogLink())
#Coef.  Std. Error       z  Pr(>|z|)   Lower 95%    Upper 95%
#──────────────────────────────────────────────────────────────────────────────
#(Intercept)  -5.39029     0.321219   -16.78    <1e-62  -6.01986    -4.76071
#log(h)        1.64376     0.106548    15.43    <1e-52   1.43493     1.85259
#d            -0.0125484   0.0015311   -8.20    <1e-15  -0.0155493  -0.00954749
p1 = rp(predict(reg) - response(reg), predict(reg), ylab="Residual", xlab="Predicted "*L"d_{-1m} / d_a - 1")
p2 = rp(predict(reg) - response(reg), D.d, ylab="Residual", xlab=L"d_a"*" [cm]")
p3 = rp(predict(reg) - response(reg), D.h, ylab="Residual", xlab="Tree Height [m]")
p4 = rp(predict(reg) - response(reg), 100 .* D.h ./ D.d, ylab="Residual", xlab=L"h/d_a")
plot(p1, p2, p3, p4, layout = grid(2, 2, link=:y))
savefig("dm1Resid.png")
#
i = ismissing.(d0[1,:])
x = predict(reg, DataFrame(d=d[i, :DBH] .+ 1.3, h=ht[i, :H_total]))
d0[1,i] .= (d[i, :DBH] .+ 1.3) .* (1. .+ x)
#
d0 = round.(d0, digits=1)
d.du1 = d0[1,:]
d.d0 = d0[2,:]
d.d0 = vec(maximum(coalesce.(Matrix(d[:,2:end]), 0.), dims=2))
d.du1 = vec(maximum(coalesce.(Matrix(d), 0.), dims=2))
D = d0 = x = i = y = y2 = reg = nothing

findall(.!ismissing.(df.DM065))[2]
x = Vector(d[117,:])
i = .!ismissing.(x)
print(x[i])
print([h; Vector(ht[117,:])][i])

# Vergleich von Interpolationsverfahren
x = coalesce.(collect(Matrix(d)'), NaN)
y = collect(Matrix(ht)')
y = vcat(repeat(h, 1, size(y, 2)), y)
y[isnan.(x)] .= NaN
#d = coalesce.(collect(Matrix(d)'), NaN)
#ht = collect(Matrix(ht)')
using Interpolations, Dierckx, StatsBase, Statistics, PPInterpolation
fu = [LinearMonotonicInterpolation();
  2:4;
  FiniteDifferenceMonotonicInterpolation();
  CardinalMonotonicInterpolation(0.1);
  CardinalMonotonicInterpolation(0.5);
  CardinalMonotonicInterpolation(0.9);
  CardinalMonotonicInterpolation(1.0);
  FritschCarlsonMonotonicInterpolation();
  FritschButlandMonotonicInterpolation();
  SteffenMonotonicInterpolation(); #12 Steffen, M. (1990) A Simple Method for Monotonic Interpolation in One Dimen sion. Astronomy and Astrophysics, 239, 443.
  AkimaMonotonicInterpolation();
  ("pp", C2());
  ("pp", C2Hyman89());  #15 Dougherty, R.L., Edelman, A. and Hyman, J.M. (1989) Nonnegativity-, Monotonic  ity-, or Convexity-Preserving Cubic and Quintic Hermite Interpolation. Mathemat  ics of Computation (United States), 52, 471-494.    https://doi.org/10.1090/S0025-5718-1989-0962209-1
  ("pp", C2MP());
  ("pp", PPInterpolation.C2MP2());
  ("pp", HuynRational());
  ("pp", VanAlbada());
  ("pp", VanLeer());
#  ("pp", FritschButland());
  ("pp", Brodlie());
  ("pp", C2HymanNonNegative());
  ("pp", Bessel());
  ("pp", Fukasawa());
  ("lagrange", "")]
function f(ipf, d, h, pos; ri = [-1, 0, 1], iso=false)
  map(ipf) do k
    res = map(eachcol(d), eachcol(h)) do x, yh
      i = .!isnan.(x)
      if (issorted(yh[i]))  # Bei kleinen Baeumen
        map(pos) do p  # Position(en) wo Wert abgefragt wird
          j = copy(i)
          #r = [p-1, p, p+1]  # Position(en) wo Messungen noetig sind
          r = p .+ ri
          l = [p]  # Position(en) dessen Wert fuer Modell nicht verwendet wird
          if (all(j[r]))
            j[l] .= false
            D = DataFrame(d = x[j], h = yh[j])
            if (iso)
              D.d = -isotonic(D.h, -D.d)
              D = combine(groupby(D, :d), [:h] .=> mean .=> [:h])
            end
            #x[p] .-
            if (isa(k, Int))
              if (sum(j) > k)
                Spline1D(D.h, D.d, k=k).(yh[p])
              else
                NaN
              end
            elseif (isa(k, Tuple))
              if (k[1] == "pp")
                makeCubicPP(D.h, D.d, PPInterpolation.SECOND_DERIVATIVE, 0.0, PPInterpolation.SECOND_DERIVATIVE, 0.0, k[2]).(yh[p])
              else
                QuadraticLagrangePP(D.h, D.d).(yh[p])
              end
            else
              interpolate(D.h, D.d, k).(yh[p])
            end
          else
            NaN
          end
        end
      else
        fill(NaN, length(pos))
      end
    end
    stack(res)
  end
end

res1 = stack(f(fu, x, y, 3:size(x, 1)-4)) .- x[3:size(x, 1)-4,:]
res2 = stack(f(fu, x, y, [size(x, 1)-3], ri = [0, 1])) .- x[end-3,:]'
res3 = stack(f(fu, x, y, collect(size(x, 1)-2 .+ (0:1)))) .- x[end-2 .+ (0:1),:]
res = vcat(res1, res2, res3)  # DM065, 1, 1.3, 3, 5, ..., 51, 53, hd>7, hd=7, hd<7

y = map(eachslice(res; dims=1)) do x
  i = map(eachslice(x; dims=1)) do y
    all(isfinite.(y))
  end
  n = sum(i)
  if (n > 0)
    map(eachslice(view(x, i, :); dims=2)) do y
      #std(y)
      #iqr(y)
      #diff(quantile(y, [0.025, 0.975]))[1]
      #diff(quantile(y, [0.05, 0.95]))[1]
      mm = quantile(y, [0.025, 0.975])
      z = filter(a -> a>=mm[1] && a<=mm[2], y)
      [n
      median(y)
      mean(z)
      std(z)
      mm[2] - mm[1]
      iqr(y)]
    end |> stack
  else
    fill(NaN, 6, size(x, 2))
  end
end |> stack
y[:,:, 3] #bhd
y[:,:, 4] #3m
y[:,:, 32] #hd<7
i = [1:16; 30:32]
plot(y[2,1,i])
plot!(y[2,3,i])
plot!(y[2,12,i])
plot!(y[2,15,i])
y[:,12,:]
join(names(d)[3:end][i], " & ")
join(round.(Int, y[1,1,i]), " & ")  #n

i = [1:16; 30:32]
t1 = ["A"; "B2"; "B3"; "B4"; "C"; "D0.1"; "D0.5"; "D0.9"; "D1"; string.('E':'T')]
using Format
f2 = function(x, y)
  t2 = abs.(round.(y, digits=2))
  t3 = stack([x .== minimum(x) for x in eachcol(t2)])
  t2 = cfmt.("%.2f", y)
  t2[t3] = "\\B " .* t2[t3]
  for i in eachindex(t1)
    print(join([x[i]; t2[i,:]], " & "), " \\\\")
    if(i % 5 == 0) print("[.35em]") end
    print("\n")
  end
  sum(t3, dims=2)'
end
f2(t1, y[2,:,i])  #median
f2(t1, y[3,:,i])  #mean
f2(t1, y[4,:,i])  #std
f2(t1, y[5,:,i])  #2.5% - 97.5%
f2(t1, y[6,:,i])  #iqr
round.(maximum(abs.(y[4,:,i]), dims=2), digits=2)'


using Integrals
x = coalesce.(collect(Matrix(d)'), NaN)
y = collect(Matrix(ht)')
y = vcat(repeat(h, 1, size(y, 2)), y)
y[isnan.(x)] .= NaN

res = map(eachcol(x), eachcol(y)) do xd, yh
  i = findall(.!isnan.(xd))
  i = unique(i -> yh[i], i)
  i = i[sortperm(yh[i])]
  if(length(i) > 1)
    map([interpolate(yh[i], xd[i], LinearMonotonicInterpolation()),
    Spline1D(yh[i], xd[i], k=min(3, length(i)-1)),
    #makeCubicPP(yh[i], xd[i], PPInterpolation.SECOND_DERIVATIVE, 0.0, PPInterpolation.SECOND_DERIVATIVE, 0.0, C2Hyman89()),
    interpolate(yh[i], xd[i], SteffenMonotonicInterpolation())]) do s
      solve(IntegralProblem((x, p)->s.(x).^2 .*pi./40000., 0., yh[i[end]]), QuadGKJL()).u
    end
  end
end |> stack

V = x[[32,34],:]'.^2 .*pi./40000 .*
[(y[33,:] .- y[32,:]) .* 2;;
(y[35,:] .- y[33,:])]
V = sum.(eachrow(replace(V, NaN => 0.)))
V += [sum(replace(x[[4;6:31],i], NaN => 0.) .^2 .*pi./40000 .* 2.) for i in 1:size(x, 2)]
res = [V'; res]
mean(res, dims=2)

using Plots.PlotMeasures
Plots.scalefontsizes(1.4)
i = res[1,:] .> 0.
p1 = rp(min.(1.5, res[1,i] ./ res[4,i]), ht.H_total[i], .2, h=1., ylab="Section", col=:red, ax=:y)
p2 = rp(res[2,:] ./ res[4,:], ht.H_total, .2, h=1., ylab="Linear", col=:red, ax=:y)
p3 = rp(min.(1.5, res[3,:] ./ res[4,:]), ht.H_total, .2, h=1., ylab="Spline", xlab="Tree height [m]", col=:red)
plot!(p1, top_margin= -4px, bottom_margin = -10px)
plot!(p2, top_margin= -10px, bottom_margin = -10px)
plot!(p3, top_margin= -10px)
plot(p1, p2, p3, layout = grid(3, 1, link=:x), size=(600, 600))
savefig("vResid.png")
Plots.scalefontsizes()

findfirst(res[2,:] ./ res[4,:] .> 1.15) # 6804
findfirst(res[2,:] ./ res[4,:] .< 0.85) # 2759
findall(d.DBH .== minimum(d.DBH))
findall(d.DBH .< 5)
Plots.scalefontsizes(1.4)
p = map([1, 22, 6804, 2759]) do j
  xd = x[:,j]
yh = y[:,j]
i = findall(.!isnan.(xd))
i = unique(i -> yh[i], i)
i = i[sortperm(yh[i])]
scatter(xd[i[2:end]], yh[i[2:end]], xlab="Diameter [cm]", ylab="Height [m]", label="Measurements", xlims=[0, maximum(xd[i])], widen=true)
annotate!(0, 0, text("TreeID: "*string(df.TreeID[j]), 14, :left, :bottom))
ih = unique(sort([range(yh[i[2]], yh[i[end]], 30); yh[i[2:end]]]))
plot!(xd[i[2:end]], yh[i[2:end]], label="Linear", lw=2)
plot!(Spline1D(yh[i], xd[i], k=2).(ih), ih, label="Quadratic Spline", lw=2)
#plot!(Spline1D(yh[i], xd[i], k=3).(ih), ih)
#plot!(makeCubicPP(yh[i], xd[i], PPInterpolation.SECOND_DERIVATIVE, 0.0, PPInterpolation.SECOND_DERIVATIVE, 0.0, C2Hyman89()).(ih), ih)
plot!(interpolate(yh[i], xd[i], SteffenMonotonicInterpolation()).(ih), ih, label="Steffen", lw=2)
end
plot!(p[2], p[3], p[4], legend=false)
plot!(p[3], yticks = 0:2:10)
#plot(p[1], p[2], p[3], p[4], layout = grid(2, 2))
plot(p[1], p[2], p[3], p[4], layout = grid(4, 1, link=:x), size=(600,3*400), leftmargin=20pt)
savefig("interpolExamples.pdf")
Plots.scalefontsizes()
D = 25.2  # DBH
H = 27.5  # Height
d = [26.1, 25, 24.2, 24.5, 23.7, 22.8, 22, 21.8, 21.3, 20.6, 20.3, 19.7, 19,
    18.2, 17.6, 16.5, 15.7, 14.9, 14, 12.6, 11.6, 11.2, 9.2, 8.2, 6.3, 4.6, 2.6]
h = .5:1:26.5

using(Plots)

V = sum((d ./ 200).^2 .* pi)
using(LambertW)

G = (D/200)^2*pi
b = lambertw(G*(H-1.3)*log(1-1.3/H)/(V)) / log(1 - 1.3/H) - 1
gM = G / (1-1.3/H)^b
plot(sqrt(gM/pi)*200/H^(b/2) .* [H;H.-h;0].^(b/2), [0;h;H],
label="through DBH", color=1, linewidth=2)

gM = ((D+1.3)/200)^2*pi
b = gM*H / V - 1
plot!(sqrt(gM/pi)*200/H^(b/2) .* [H;H.-h;0].^(b/2), [0;h;H],
label="through DBH+1.3", color=3, linewidth=2)

v = V - G * 1.3
b = G*(H-1.3) / v - 1
plot!([D;sqrt(G/pi)*200/(H-1.3)^(b/2) .*
[H-1.3;H.-h[2:end];0].^(b/2)], [0;1.3;h[2:end];H], label="Cylinder",
linestyle=:dash, color=2, linewidth=2)

scatter!(d, h, xlabel="Diameter [cm]", ylabel="Height [m]", label="Measured")
scatter!([D], [1.3], label="DBH", color=:black)
savefig("/tmp/trunk.pdf")

savefig("trunk.pdf")

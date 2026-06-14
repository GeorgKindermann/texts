using Plots
using Printf
using Plots.Measures
using Optim, Roots
using Statistics

# =========================================================================
# 1. FORSTLICHE DATENSTRUKTUR (D in cm, H in m)
# =========================================================================
struct TaperProfile
    H::Float64   # H_total: Gesamthöhe (m)
    hw::Float64  # H_ip: Inflexionspunkt (m)
    c0::Float64  # d_0: Durchmesser am Boden (cm)
    dw::Float64  # d_ip: Durchmesser am Inflexionspunkt (cm)
    dt::Float64  # d_t: Terminaler Durchmesser an der Spitze (cm)
    a::Float64   # Parameter a (oberes Segment)
    b::Float64   # Parameter b (oberes Segment)
    c1::Float64  # Parameter c1 (unteres Segment, negativ)
    c2::Float64  # Parameter c2 (unteres Segment, positiv)
end

function get_bark_coefficients()
    alpha = 0.85149 / 10.0
    beta  = 0.60934 / 10.0
    gamma = -0.00228 / 10.0
    
    A0 = alpha^2
    A1 = -2.0 * alpha * (1.0 - beta)
    A2 = (1.0 - beta)^2 + 2.0 * alpha * gamma
    A3 = -2.0 * gamma * (1.0 - beta)
    A4 = gamma^2
    
    return (A0, A1, A2, A3, A4)
end

inline_plogis(x) = 1.0 / (1.0 + exp(-x))

function binomial_coeff(n::Int, k::Int)
    if k < 0 || k > n return 0.0 end
    res = 1.0
    for i in 1:k; res *= (n - i + 1) / i; end
    return res
end

# =========================================================================
# 2. UNIVERSAL CALIBRATION SOLVER (Rein in cm und m)
# =========================================================================

# INTERNE KERN-RECHENFUNKTION (Erwartet ein fest vorgegebenes dr)
function _solve_taper_with_dr(h_meas::Vector{Float64}, d_meas::Vector{Float64}, H::Float64, dt::Float64, hw::Float64, dr::Float64; 
                              max_c0_factor::Float64=2.0, h_ref::Union{Float64, Nothing}=nothing)
    n_pts = length(h_meas)
    limit_c0 = max_c0_factor * dr
    
    # Biologische Startwerte aus den GLMs basierend auf dem übergebenen dr
    b_glm  = clamp(inline_plogis(-4.31616 + 21.7497 * (1.0 / max(0.1, H)) + 1.0813 * log1p(dr)), 0.01, 1.0)
    c2_glm = clamp(inline_plogis(7.0974 - 3.17164 * log1p(max(0.1, H)) + 0.631264 * log1p(dr)), 0.01, 1.0)
    
    local_b   = b_glm
    local_c2  = c2_glm
    local_a   = 0.0
    local_dip = 0.0
    sip       = 0.0
    local_c1  = 0.0
    local_c0  = 0.0

    # ---------------------------------------------------------------------
    # N = 0 PUNKTE: REIN HÖHENBASIERTE SCHÄTZUNG
    # ---------------------------------------------------------------------
      if n_pts == 0
        local_c0 = min(dt + H, limit_c0)
        f_zero_n0 = (a_trial) -> begin
            dw_t = a_trial * (H - hw)^local_b + dt
            sw_t = -a_trial * local_b * (H - hw)^(local_b - 1.0)
            c1_t = sw_t / (local_c2 * hw^(local_c2 - 1.0))
            return (dw_t - c1_t * hw^local_c2) - local_c0
        end
        #local_a = find_zero(f_zero_n0, (0.0001, local_c0 * 2.0), Bisection())
        prob_n0 = ZeroProblem(f_zero_n0, (0.0001, local_c0 * 2.0))
        local_a = solve(prob_n0, Bisection(), verbose=false)
        if isnan(local_a); local_a = (local_c0 - dt) / ((H - hw)^local_b); end
        local_dip = local_a * (H - hw)^local_b + dt
        sip = -local_a * local_b * (H - hw)^(local_b - 1.0)
        local_c1 = sip / (local_c2 * hw^(local_c2 - 1.0))

    # ---------------------------------------------------------------------
    # N = 1 PUNKT: EINZELMESSUNG
    # ---------------------------------------------------------------------
    elseif n_pts == 1
        h1, d1 = h_meas[1], d_meas[1]
        if h1 <= hw
            Omega = (local_b * hw) / (local_c2 * (H - hw)) * (1.0 - (h1 / hw)^local_c2)
            trial_dip = (d1 + dt * Omega) / (1.0 + Omega)
            trial_a = (trial_dip - dt) / (H - hw)^local_b
            trial_sip = -trial_a * local_b * (H - hw)^(local_b - 1.0)
            trial_c1 = trial_sip / (local_c2 * hw^(local_c2 - 1.0))
            trial_c0 = trial_dip - trial_c1 * hw^local_c2
            
            if trial_c0 > limit_c0
                local_c0 = limit_c0
                f_zero_c2 = (c2_trial) -> begin
                    Omega_loop = (local_b * hw) / (c2_trial * (H - hw)) * (1.0 - (h1 / hw)^c2_trial)
                    tdip = (d1 + dt * Omega_loop) / (1.0 + Omega_loop)
                    ta = (tdip - dt) / (H - hw)^local_b
                    tsip = -ta * local_b * (H - hw)^(local_b - 1.0)
                    return tdip - (tsip / (c2_trial * hw^(c2_trial - 1.0))) * hw^c2_trial - limit_c0
                end
                local_c2 = find_zero(f_zero_c2, (0.01, 1.0), Bisection())
                Omega = (local_b * hw) / (local_c2 * (H - hw)) * (1.0 - (h1 / hw)^local_c2)
                local_dip = (d1 + dt * Omega) / (1.0 + Omega)
                local_a = (local_dip - dt) / (H - hw)^local_b
                sip = -local_a * local_b * (H - hw)^(local_b - 1.0)
                local_c1 = sip / (local_c2 * hw^(local_c2 - 1.0))
            else
                local_dip, local_a, sip, local_c1, local_c0 = trial_dip, trial_a, trial_sip, trial_c1, trial_c0
            end
        else
            local_a = (d1 - dt) / (H - h1)^local_b
            local_dip = local_a * (H - hw)^local_b + dt
            sip = -local_a * local_b * (H - hw)^(local_b - 1.0)
            local_c1 = sip / (local_c2 * hw^(local_c2 - 1.0))
            local_c0 = local_dip - local_c1 * hw^local_c2
            if local_c0 > limit_c0
                local_c0 = limit_c0
                local_c2 = clamp((sip * hw) / (local_dip - local_c0), 0.01, 1.0)
                local_c1 = (local_dip - local_c0) / (hw^local_c2)
            end
        end

    # ---------------------------------------------------------------------
    # N = 2 PUNKTE
    # ---------------------------------------------------------------------
    elseif n_pts == 2
        h1, h2 = h_meas[1], h_meas[2]
        d1, d2 = d_meas[1], d_meas[2]

              if h1 >= hw && h2 >= hw
            term_num = max(1e-6, d2 - dt) / max(1e-6, d1 - dt)
            term_den = max(1e-6, H - h2) / max(1e-6, H - h1)
            local_b = clamp(log(term_num) / log(term_den), 0.2, 0.9) # Grenzen geweitet gegen Flachheit
            local_a = (d1 - dt) / max(1e-6, H - h1)^local_b
            local_dip = local_a * (H - hw)^local_b + dt
            sip = -local_a * local_b * (H - hw)^(local_b - 1.0)
            local_c2 = c2_glm # Schützt die untere Krümmung vor dem Kollaps
            local_c1 = sip / (local_c2 * hw^(local_c2 - 1.0))
            local_c0 = min(local_dip - local_c1 * hw^local_c2, limit_c0)
            
        elseif h1 < hw && h2 >= hw
            find_b_trajectory_error = (b_trial) -> begin
                a_t = (d2 - dt) / max(1e-6, H - h2)^b_trial
                dw_t = a_t * (H - hw)^b_trial + dt
                sw_t = -a_t * b_trial * (H - hw)^(b_trial - 1.0)
                c1_t = sw_t / (local_c2 * hw^(local_c2 - 1.0))
                c0_t = dw_t - c1_t * hw^local_c2
                return (c0_t + c1_t * h1^local_c2) - d1
            end
            # 1. Berechne die Fehler an den Intervallgrenzen
            err_low  = find_b_trajectory_error(0.01)
            err_high = find_b_trajectory_error(1.0)

            # 2. Prüfe, ob ein Vorzeichenwechsel vorliegt
            if sign(err_low) != sign(err_high)
            # Sicherer Aufruf, da mathematisch lösbar
              local_b = find_zero(find_b_trajectory_error, (0.01, 1.0), Bisection())
            else
            # FALLBACK: Wenn das System unlösbar ist, weiche auf das biologische GLM-b aus
            # und passe das obere Segment an, ohne abzustürzen.
              local_b = b_glm
            end
            local_a = (d2 - dt) / max(1e-6, H - h2)^local_b
            local_dip = local_a * (H - hw)^local_b + dt
            sip = -local_a * local_b * (H - hw)^(local_b - 1.0)
            local_c1 = sip / (local_c2 * hw^(local_c2 - 1.0))
                local_c0 = min(local_dip - local_c1 * hw^local_c2, limit_c0)
        else
            # BEIDE UNTEN: Fixieren Sie b und passen Sie die Krümmung unten an
            local_b = b_glm
            f_opt_c2 = (c2_trial) -> begin
                # Berechne c0 und c1 direkt aus den zwei Punkten für dieses c2
                # Lineares Gleichungssystem für c0 und c1:
                # c0 + c1 * h1^c2 = d1
                # c0 + c1 * h2^c2 = d2
                denom = max(1e-6, h2^c2_trial - h1^c2_trial)
                c1_t = (d2 - d1) / denom
                c0_t = d1 - c1_t * h1^c2_trial
                
                # Wie gut passt das zum biologischen Wendepunkt-GLM?
                # Errechnung des theoretischen Übergangs
                dw_t = min(limit_c0, c0_t + c1_t * hw^c2_trial)
                sw_t = c1_t * c2_trial * hw^(c2_trial - 1.0)
                a_t = (dw_t - dt) / ((H - hw)^local_b)
                sw_expected = -a_t * local_b * (H - hw)^(local_b - 1.0)
                
                return (sw_t - sw_expected)^2 # Minimiert den Knick am Übergang
            end
            res = optimize(f_opt_c2, 0.1, 1.5, Brent())
            local_c2 = Optim.minimizer(res)
            denom = max(1e-6, h2^local_c2 - h1^local_c2)
            local_c1 = (d2 - d1) / denom
            local_c0 = d1 - local_c1 * h1^local_c2
            local_dip = local_c0 + local_c1 * hw^local_c2
            sip = local_c1 * local_c2 * hw^(local_c2 - 1.0)
            local_a = (local_dip - dt) / ((H - hw)^local_b)
        end

    # ---------------------------------------------------------------------
    # N >= 3 PUNKTE
    # ---------------------------------------------------------------------
    elseif n_pts >= 3
        d_meas_sq = d_meas .^ 2
        
        # Prüfen, wie die Punkte relativ zu hw liegen
        pts_below = sum(h_meas .<= hw)
        pts_above = sum(h_meas .> hw)
        
        # =================================================================
        # SONDERFALL A: Alle Messpunkte liegen im unteren Stammfußbereich
        # =================================================================
        if pts_above == 0
            if !isnothing(h_ref)
                # --- SONDERFALL A MIT ANKERPUNKT (Dein Originalcode) ---
                lower_bounds = [0.1]
                upper_bounds = [1.0]
                initial_guess = [c2_glm]
                
                idx_anchor = argmin(abs.(h_meas .- h_ref))
                h_anchor = h_meas[idx_anchor]
                d_anchor = d_meas[idx_anchor]
                
                f_opt_all_below = (x) -> begin
                    c2_t = x[1]
                    b_fixed = b_glm
                    
                    k_slope = -(b_fixed * max(1e-6, H - hw)^(b_fixed - 1.0)) / (c2_t * hw^(c2_t - 1.0))
                    Nenner = max(1e-6, H - hw)^b_fixed - k_slope * (hw^c2_t - h_anchor^c2_t)
                    
                    if abs(Nenner) < 1e-6; Nenner = sign(Nenner) * 1e-6; end
                    
                    a_t = (d_anchor - dt) / Nenner
                    c1_t = a_t * k_slope
                    c0_t = d_anchor - c1_t * h_anchor^c2_t
                    
                    if c0_t > limit_c0
                        c0_t = limit_c0
                        c1_t = (d_anchor - limit_c0) / max(1e-6, h_anchor^c2_t)
                    end
                    
                    total_error = 0.0
                    for i in 1:n_pts
                        d_mod_sq = (c0_t + c1_t * h_meas[i]^c2_t)^2
                        total_error += (d_mod_sq - d_meas_sq[i])^2
                    end
                    return total_error
                end
                
                res = optimize(f_opt_all_below, lower_bounds, upper_bounds, initial_guess, Fminbox(NelderMead()))
                local_c2 = Optim.minimizer(res)[1]
                local_b = b_glm
                
                k_slope = -(local_b * max(1e-6, H - hw)^(local_b - 1.0)) / (local_c2 * hw^(local_c2 - 1.0))
                Nenner = max(1e-6, H - hw)^local_b - k_slope * (hw^local_c2 - h_anchor^local_c2)
                if abs(Nenner) < 1e-6; Nenner = sign(Nenner) * 1e-6; end
                
                local_a = (d_anchor - dt) / Nenner
                local_c1 = local_a * k_slope
                local_c0 = d_anchor - local_c1 * h_anchor^local_c2
                
                if local_c0 > limit_c0
                    local_c0 = limit_c0
                    local_c1 = (d_anchor - limit_c0) / max(1e-6, h_anchor^local_c2)
                end
            else
                # --- NEU: SONDERFALL A OHNE ANKERPUNKT (Freier Fit im unteren Bereich) ---
                # Wir optimieren c2 (Krümmung) und c0 (theoretischer Basisdurchmesser) frei.
                lower_bounds = [0.1, 1.0]           # [c2, c0]
                upper_bounds = [1.0, limit_c0]      # c0 wird hier direkt über die Grenze gedeckelt
                initial_guess = [c2_glm, min(limit_c0, sqrt(d_meas_sq[1]))] # Startwert für c0 aus dem ersten Messpunkt schätzen
                
                local_b = b_glm
                
                f_opt_all_below_free = (x) -> begin
                    c2_t, c0_t = x[1], x[2]
                    
                    # Da wir keine Messpunkte oben haben, schätzen wir c1 so, 
                    # dass der Übergang zu einem plausiblen oberen Baumteil (b_glm, dt) passt.
                    # Mathematische Bedingung für Stetigkeit und Ableitung bei hw:
                    # c1_t * c2_t * hw^(c2_t - 1) = -a_t * b_fixed * (H - hw)^(b_fixed - 1)
                    # c0_t + c1_t * hw^c2_t = a_t * (H - hw)^b_fixed + dt
                    
                    A = hw^c2_t
                    B = (c2_t * hw^(c2_t - 1.0) * max(1e-6, H - hw)) / local_b
                    Nenner = A + B
                    
                    # c1 berechnen aus der Differenz zum asymptotischen Zopf-Durchmesser dt
                    c1_t = (dt - c0_t) / max(1e-6, Nenner)
                    
                    total_error = 0.0
                    for i in 1:n_pts
                        d_mod_sq = (c0_t + c1_t * h_meas[i]^c2_t)^2
                        total_error += (d_mod_sq - d_meas_sq[i])^2
                    end
                    return total_error
                end
                
                res = optimize(f_opt_all_below_free, lower_bounds, upper_bounds, initial_guess, Fminbox(NelderMead()))
                best = Optim.minimizer(res)
                local_c2, local_c0 = best[1], best[2]
                
                # Finale Parameter-Ableitung für den oberen (ungemessenen) Teil
                A_f = hw^local_c2
                B_f = (local_c2 * hw^(local_c2 - 1.0) * max(1e-6, H - hw)) / local_b
                local_c1 = (dt - local_c0) / max(1e-6, A_f + B_f)
                
                # a berechnen, damit die obere Kurve stetig anschließt
                dw_final = local_c0 + local_c1 * hw^local_c2
                local_a = (dw_final - dt) / max(1e-6, H - hw)^local_b
            end

        # =================================================================
        # SONDERFALL B: Alle Messpunkte liegen im oberen Schaftbereich
        # =================================================================
        elseif pts_below == 0
            if !isnothing(h_ref)
                # --- SONDERFALL B MIT ANKERPUNKT (Dein Originalcode) ---
                lower_bounds = [0.2]
                upper_bounds = [1.0]
                initial_guess = [b_glm]
                
                idx_anchor = argmin(abs.(h_meas .- h_ref))
                h_anchor = h_meas[idx_anchor]
                d_anchor = d_meas[idx_anchor]
                
                f_opt_all_above = (x) -> begin
                    b_t = x[1]
                    c2_fixed = c2_glm
                    
                    a_t = (d_anchor - dt) / max(1e-6, H - h_anchor)^b_t
                    dw_t = a_t * (H - hw)^b_t + dt
                    sw_t = -a_t * b_t * max(1e-6, H - hw)^(b_t - 1.0)
                    
                    c1_t = sw_t / (c2_fixed * hw^(c2_fixed - 1.0))
                    c0_t = dw_t - c1_t * hw^c2_fixed
                    
                    if c0_t > limit_c0
                        c0_t = limit_c0
                        c1_t = (dw_t - limit_c0) / hw^c2_fixed
                    end
                    
                    total_error = 0.0
                    for i in 1:n_pts
                        d_mod_sq = (a_t * (H - h_meas[i])^b_t + dt)^2
                        total_error += (d_mod_sq - d_meas_sq[i])^2
                    end
                    return total_error
                end
                
                res = optimize(f_opt_all_above, lower_bounds, upper_bounds, initial_guess, Fminbox(NelderMead()))
                local_b = Optim.minimizer(res)[1]
                local_c2 = c2_glm
                
                local_a = (d_anchor - dt) / max(1e-6, H - h_anchor)^local_b
                dw_final = local_a * (H - hw)^local_b + dt
                sw_final = -local_a * local_b * max(1e-6, H - hw)^(local_b - 1.0)
                
                local_c1 = sw_final / (local_c2 * hw^(local_c2 - 1.0))
                local_c0 = dw_final - local_c1 * hw^local_c2
                
                if local_c0 > limit_c0
                    local_c0 = limit_c0
                    local_c1 = (dw_final - limit_c0) / hw^local_c2
                end
              else
                # --- NEU: SONDERFALL B OHNE ANKERPUNKT (Freier Fit im oberen Bereich) ---
                # Wir schätzen einen sinnvollen Startwert für 'a' aus dem ersten oberen Messpunkt
                # d_meas[1] = a * (H - h_meas[1])^b_glm + dt  => nach 'a' umgeformt:
                a_init = (d_meas[1] - dt) / max(1e-6, H - h_meas[1])^b_glm
                a_init = clamp(a_init, 0.01, 45.0) # Absicherung für den Startwert
                
                lower_bounds = [0.2, 0.001]
                upper_bounds = [1.0, 50.0] 
                initial_guess = [b_glm, a_init] # <--- DYNAMISCHER STARTWERT STATT 10.0
                
                local_c2 = c2_glm
                
                f_opt_all_above_free = (x) -> begin
                    b_t, a_t = x[1], x[2]
                    
                    total_error = 0.0
                    for i in 1:n_pts
                        # max(1e-6, ...) verhindert negative Basen bei der Potenzierung am Wipfel
                        d_mod = a_t * max(1e-6, H - h_meas[i])^b_t + dt
                        total_error += (d_mod^2 - d_meas_sq[i])^2
                    end
                    return total_error
                end
                
                res = optimize(f_opt_all_above_free, lower_bounds, upper_bounds, initial_guess, Fminbox(NelderMead()))
                best = Optim.minimizer(res)
                local_b, local_a = best[1], best[2]
                
                # Untere (ungemessene) Parameter analytisch über den Übergang hw herleiten
                dw_final = local_a * max(1e-6, H - hw)^local_b + dt
                sw_final = -local_a * local_b * max(1e-6, H - hw)^(local_b - 1.0)
                
                local_c1 = sw_final / (local_c2 * hw^(local_c2 - 1.0))
                local_c0 = dw_final - local_c1 * hw^local_c2
                
                # Auch hier die c0 Limit-Prüfung analog anwenden
                if local_c0 > limit_c0
                    local_c0 = limit_c0
                    local_c1 = (dw_final - limit_c0) / hw^local_c2
                end
            end
        # =================================================================
        # NORMALFALL: Punkte sind sowohl über als auch unter hw verteilt
        # =================================================================
        else
          # --- STRATEGIE NORMAL 1: Genau 3 Punkte -> Stabilisierte 2D-Interpolation ---
    if n_pts == 3
        # Hilfsfunktion zur Parameterberechnung aus dem Ankerpunkt
        function calc_parameters_3pts(b, c2, h_anc, d_anc)
            k_slope = -(b * max(1e-6, H - hw)^(b - 1.0)) / (c2 * hw^(c2 - 1.0))
            if h_anc <= hw
                Nenner = max(1e-6, H - hw)^b - k_slope * (hw^c2 - h_anc^c2)
                if abs(Nenner) < 1e-6; Nenner = sign(Nenner) * 1e-6; end
                a = (d_anc - dt) / Nenner
                c1 = a * k_slope
                c0 = d_anc - c1 * h_anc^c2
            else
                a = (d_anc - dt) / max(1e-6, H - h_anc)^b
                dw = a * max(1e-6, H - hw)^b + dt
                sw = -a * b * max(1e-6, H - hw)^(b - 1.0)
                c1 = sw / (c2 * hw^(c2 - 1.0))
                c0 = dw - c1 * hw^c2
            end
            return a, c0, c1
        end

        # Überprüfen, ob ein Ankerpunkt gewünscht ist
        if !isnothing(h_ref)
            # --- VARIANTE 1: MIT ANKERPUNKT (Fixiert bei h_ref) ---
            idx_anchor = argmin(abs.(h_meas .- h_ref))
            h_anchor = h_meas[idx_anchor]
            d_anchor = d_meas[idx_anchor]
            
            lower_bounds = [0.01, 0.01]
            upper_bounds = [1.0, 1.0]
            initial_guess = [b_glm, max(0.6, c2_glm)]
            
            f_opt_3pts = (x) -> begin
                b_t, c2_t = x[1], x[2]
                a_t, c0_t, c1_t = calc_parameters_3pts(b_t, c2_t, h_anchor, d_anchor)
                
                if c0_t > limit_c0
                    c0_t = limit_c0
                    if h_anchor <= hw
                        c1_t = (d_anchor - limit_c0) / max(1e-6, h_anchor^c2_t)
                        dw_t = c0_t + c1_t * hw^c2_t
                        a_t = (dw_t - dt) / max(1e-6, H - hw)^b_t
                    else
                        dw_t = a_t * max(1e-6, H - hw)^b_t + dt
                        c1_t = (dw_t - limit_c0) / hw^c2_t
                    end
                end
                
                total_error = 0.0
                for i in 1:3
                    d_mod = (h_meas[i] <= hw) ? (c0_t + c1_t * h_meas[i]^c2_t) : (a_t * max(1e-6, H - h_meas[i])^b_t + dt)
                    total_error += (d_mod^2 - d_meas_sq[i])^2
                end
                return total_error
            end
            
            res = optimize(f_opt_3pts, lower_bounds, upper_bounds, initial_guess, Fminbox(NelderMead()))
            best = Optim.minimizer(res)
            
            local_b, local_c2 = best[1], best[2]
            local_a, local_c0, local_c1 = calc_parameters_3pts(local_b, local_c2, h_anchor, d_anchor)

        else
            # --- VARIANTE 2: OHNE ANKERPUNKT (Freier Fit über alle Punkte) ---
            lower_bounds = [0.01, 0.01]  # x[1] = b, x[2] = a
            upper_bounds = [1.0, 50.0]   
            initial_guess = [b_glm, 10.0] 
            
            c2_fixed = max(0.6, c2_glm)  # c2 bleibt fest, da wir ohne Ankerpunkt sonst zu viele Freiheitsgrade für 3 Punkte haben
            
            f_opt_free = (x) -> begin
                b_t, a_t = x[1], x[2]
                
                # Stetigkeit und Ableitung direkt am Übergangspunkt hw berechnen
                dw_t = a_t * max(1e-6, H - hw)^b_t + dt
                sw_t = -a_t * b_t * max(1e-6, H - hw)^(b_t - 1.0)
                c1_t = sw_t / (c2_fixed * hw^(c2_fixed - 1.0))
                c0_t = dw_t - c1_t * hw^c2_fixed
                
                if c0_t > limit_c0
                    c0_t = limit_c0
                    c1_t = (dw_t - limit_c0) / hw^c2_fixed
                end
                
                total_error = 0.0
                for i in 1:3
                    d_mod = (h_meas[i] <= hw) ? (c0_t + c1_t * h_meas[i]^c2_fixed) : (a_t * max(1e-6, H - h_meas[i])^b_t + dt)
                    total_error += (d_mod^2 - d_meas_sq[i])^2
                end
                return total_error
            end
            
            res = optimize(f_opt_free, lower_bounds, upper_bounds, initial_guess, Fminbox(NelderMead()))
            best = Optim.minimizer(res)
            
            local_b = best[1]
            local_a = best[2]
            local_c2 = c2_fixed
            
            dw_final = local_a * max(1e-6, H - hw)^local_b + dt
            sw_final = -local_a * local_b * max(1e-6, H - hw)^(local_b - 1.0)
            local_c1 = sw_final / (local_c2 * hw^(local_c2 - 1.0))
            local_c0 = dw_final - local_c1 * hw^local_c2
            if local_c0 > limit_c0
                local_c0 = limit_c0
                local_c1 = (dw_final - limit_c0) / hw^local_c2
            end
        end



            # --- STRATEGIE NORMAL 2: Mehr als 3 Punkte -> 3D-Ausgleich (Regression) ---
               else
                if !isnothing(h_ref)
                    # --- NEU: Mehr als 3 Punkte MIT ANKERPUNKT ---
                    # Hier bestimmen wir den Ankerpunkt aus den Daten
                    idx_anchor = argmin(abs.(h_meas .- h_ref))
                    h_anchor = h_meas[idx_anchor]
                    d_anchor = d_meas[idx_anchor]

                    # Da der Ankerpunkt die Kurve fixiert, reduzieren wir das Problem auf 2D (b und c2)
                    lower_bounds = [0.2, 0.1]
                    upper_bounds = [1.0, 1.0]
                    initial_guess = [b_glm, c2_glm]

                    # Nutzen der algebraischen Verknüpfung über den Ankerpunkt (wie im 3-Punkt-Fall)
                    function calc_parameters_multi_anchor(b, c2, h_anc, d_anc)
                        k_slope = -(b * max(1e-6, H - hw)^(b - 1.0)) / (c2 * hw^(c2 - 1.0))
                        if h_anc <= hw
                            Nenner = max(1e-6, H - hw)^b - k_slope * (hw^c2 - h_anc^c2)
                            if abs(Nenner) < 1e-6; Nenner = sign(Nenner) * 1e-6; end
                            a = (d_anc - dt) / Nenner
                            c1 = a * k_slope
                            c0 = d_anc - c1 * h_anc^c2
                        else
                            a = (d_anc - dt) / max(1e-6, H - h_anc)^b
                            dw = a * max(1e-6, H - hw)^b + dt
                            sw = -a * b * max(1e-6, H - hw)^(b - 1.0)
                            c1 = sw / (c2 * hw^(c2 - 1.0))
                            c0 = dw - c1 * hw^c2
                        end
                        return a, c0, c1
                    end

                    f_opt_multi_anchor = (x) -> begin
                        b_t, c2_t = x[1], x[2]
                        a_t, c0_t, c1_t = calc_parameters_multi_anchor(b_t, c2_t, h_anchor, d_anchor)
                        
                        if c0_t > limit_c0
                            c0_t = limit_c0
                            if h_anchor <= hw
                                c1_t = (d_anchor - limit_c0) / max(1e-6, h_anchor^c2_t)
                                dw_t = c0_t + c1_t * hw^c2_t
                                a_t = (dw_t - dt) / max(1e-6, H - hw)^b_t
                            else
                                dw_t = a_t * max(1e-6, H - hw)^b_t + dt
                                c1_t = (dw_t - limit_c0) / hw^c2_t
                            end
                        end
                        
                        d_mod(h) = (h <= hw) ? (c0_t + c1_t * h^c2_t) : (a_t * max(1e-6, H - h)^b_t + dt)
                        
                        total_error = 0.0
                        for i in 1:n_pts
                            total_error += (d_mod(h_meas[i])^2 - d_meas_sq[i])^2
                        end
                        return total_error
                    end

                    res = optimize(f_opt_multi_anchor, lower_bounds, upper_bounds, initial_guess, Fminbox(NelderMead()))
                    best = Optim.minimizer(res)
                    local_b, local_c2 = best[1], best[2]
                    local_a, local_c0, local_c1 = calc_parameters_multi_anchor(local_b, local_c2, h_anchor, d_anchor)

                else
                    # --- DEIN ORIGINALCODE: Mehr als 3 Punkte OHNE ANKERPUNKT ---
                    lower_bounds = [0.2, 0.1, 0.001]
                    upper_bounds = [1.0, 1.0, 5.0]
                    initial_guess = [b_glm, c2_glm, 1.0]
                    
                    function calc_lower_parameters(b, c2, a)
                        dw = a * max(1e-6, H - hw)^b + dt
                        sw = -a * b * max(1e-6, H - hw)^(b - 1.0)
                        c1 = sw / (c2 * hw^(c2 - 1.0))
                        c0 = dw - c1 * hw^c2
                        return c0, c1
                    end

                    f_opt_multi = (x) -> begin
                        b_t, c2_t, a_t = x[1], x[2], x[3]
                        c0_t, c1_t = calc_lower_parameters(b_t, c2_t, a_t)
                        
                        if c0_t > limit_c0
                            c0_t = limit_c0
                            dw_t = a_t * max(1e-6, H - hw)^b_t + dt
                            c1_t = (dw_t - limit_c0) / hw^c2_t
                        end
                        
                        d_mod(h) = (h <= hw) ? (c0_t + c1_t * h^c2_t) : (a_t * max(1e-6, H - h)^b_t + dt)
                        
                        total_error = 0.0
                        for i in 1:n_pts
                            total_error += (d_mod(h_meas[i])^2 - d_meas_sq[i])^2
                        end
                        return total_error
                    end
                    
                    res = optimize(f_opt_multi, lower_bounds, upper_bounds, initial_guess, Fminbox(NelderMead()))
                    best = Optim.minimizer(res)
                    local_b, local_c2, local_a = best[1], best[2], best[3]
                    local_c0, local_c1 = calc_lower_parameters(local_b, local_c2, local_a)
                end
            end # Schließt die "if n_pts == 3"-Unterscheidung
         
            # Gemeinsame, abschließende Limit-Prüfung für den Normalfall
            if local_c0 > limit_c0
                local_c0 = limit_c0
                
                if !isnothing(h_ref)
                    if h_anchor <= hw
                        local_c1 = (d_anchor - limit_c0) / max(1e-6, h_anchor^local_c2)
                        dw_final = local_c0 + local_c1 * hw^local_c2
                        local_a = (dw_final - dt) / max(1e-6, H - hw)^local_b
                    else
                        dw_final = local_a * max(1e-6, H - hw)^local_b + dt
                        local_c1 = (dw_final - limit_c0) / hw^local_c2
                    end
                else
                    dw_final = local_a * max(1e-6, H - hw)^local_b + dt
                    local_c1 = (dw_final - limit_c0) / hw^local_c2
                end
            end
        end
        
        local_dip = local_a * (H - hw)^local_b + dt
    end

    return TaperProfile(H, hw, local_c0, local_dip, dt, local_a, local_b, local_c1, local_c2)
end

# DIE HAUPTFUNKTION: REGELT DIE IMPLIZITE BHD-BESTIMMUNG ZIRKELFREI
function solve_taper(h_meas::Vector{Float64}, d_meas::Vector{Float64}, H::Float64; max_c0_factor::Float64=2.0, h_ref::Union{Float64, Nothing}=nothing)
    n_pts = length(h_meas)
    if n_pts != length(d_meas) error("Vektoren müssen gleich lang sein.") end
    
    dt = H >= 1.3 ? 0.8 : 0.1 + 0.7 * (1.0 - ((1.3 - H) / 1.3)^2)
    hw = 1.8104 * (1.0 + H)^0.358966 - 1.8104

    # 1. AUSNAHME FÜR KLEINE BÄUME (H < 1.3m): Hier gibt es keinen BHD!
    if H < 1.3
        dr_small = dt + H
        return _solve_taper_with_dr(h_meas, d_meas, H, dt, hw, dr_small; max_c0_factor=max_c0_factor, h_ref=h_ref)
    end

    # 2. PRÜFEN, OB EIN ECHTER BHD VORLIEGT (Messung exakt bei h = 1.3)
    dbh_idx = findfirst(h -> abs(h - 1.3) < 1e-4, h_meas)
    
    if dbh_idx !== nothing
        # Fall A: BHD direkt gemessen. Keine iterative Schleife nötig!
        dbh_measured = d_meas[dbh_idx]
        dr = dbh_measured + 1.3
      return _solve_taper_with_dr(h_meas, d_meas, H, dt, hw, dr; max_c0_factor=max_c0_factor, h_ref=h_ref)
      else
        # Fall B: Kein BHD gemessen (und H >= 1.3m). Wir suchen den impliziten BHD!
        
        # 1. ZUERST DEN SPEZIALFALL OHNE MESSUNG ABFANGEN (Zurück zu Ihrer Originalversion)
        if n_pts == 0
            start_dbh = H # Ihr originaler Startwert für n_pts == 0
            
            f_outer_root = (dbh_trial) -> begin
                dr_trial = dbh_trial + 1.3
                prof = _solve_taper_with_dr(h_meas, d_meas, H, dt, hw, dr_trial; max_c0_factor=max_c0_factor, h_ref=h_ref)
                
                d_model_13 = if 1.3 <= prof.hw
                    prof.c0 + prof.c1 * 1.3^prof.c2
                else
                    prof.a * (H - 1.3)^prof.b + prof.dt
                end
                return d_model_13 - dbh_trial
            end
            
            final_dbh = find_zero(f_outer_root, start_dbh, Order0())
            dr_final = final_dbh + 1.3
            return _solve_taper_with_dr(h_meas, d_meas, H, dt, hw, dr_final; max_c0_factor=max_c0_factor, h_ref=h_ref)
        end

        # 2. DER REPARIERTE ZWEIG FÜR N >= 1 MESSREIHEN
        pts_above = d_meas[h_meas .> 1.3]
        pts_below = d_meas[h_meas .<= 1.3]
        
        lower_bound_from_above = isempty(pts_above) ? 0.0 : maximum(pts_above)
        upper_bound_from_below = isempty(pts_below) ? Inf : minimum(pts_below)
        
        min_dbh = max(0.1, lower_bound_from_above)
        max_dbh = isinf(upper_bound_from_below) ? (min_dbh * 3.0 + 10.0) : upper_bound_from_below

        f_outer_root_with_pts = (dbh_trial) -> begin
            dr_trial = dbh_trial + 1.3
            prof = try
                _solve_taper_with_dr(h_meas, d_meas, H, dt, hw, dr_trial; max_c0_factor=max_c0_factor, h_ref=h_ref)
            catch e
                return Inf 
            end
            
            d_model_13 = if 1.3 <= prof.hw
                prof.c0 + prof.c1 * 1.3^prof.c2
            else
                prof.a * (H - 1.3)^prof.b + prof.dt
            end
            return d_model_13 - dbh_trial
        end
        
        if sign(f_outer_root_with_pts(min_dbh)) == sign(f_outer_root_with_pts(max_dbh))
            final_dbh = (min_dbh + max_dbh) / 2.0
        else
            final_dbh = find_zero(f_outer_root_with_pts, (min_dbh, max_dbh), Bisection())
        end
        
        dr_final = final_dbh + 1.3
        return _solve_taper_with_dr(h_meas, d_meas, H, dt, hw, dr_final; max_c0_factor=max_c0_factor, h_ref=h_ref)
    end
end

# =========================================================================
# 3. APPLICATIONS (Durchmesser in cm, Integrale mit Faktor 40000)
# =========================================================================

function get_diameter_at_height(p::TaperProfile, h::Float64)
    if h < 0.0 || h > p.H return 0.0 end
    return h <= p.hw ? p.c0 + p.c1 * h^p.c2 : p.a * (p.H - h)^p.b + p.dt
end

function get_diameter_under_bark(p::TaperProfile, h::Float64)
    d_ob = get_diameter_at_height(p, h)
    if d_ob <= 0.0 return 0.0 end
    dbt = (0.85149 + 0.60934 * d_ob - 0.00228 * d_ob^2) / 10.0
    return max(0.0, d_ob - dbt)
end

function get_height_at_diameter(p::TaperProfile, d_target::Float64)
    if d_target >= p.c0 return 0.0 end
    if d_target <= p.dt return p.H end
    if d_target >= p.dw
        return max(0.0, ((d_target - p.c0) / p.c1)^(1.0 / p.c2))
    else
        return min(p.H, p.H - ((d_target - p.dt) / p.a)^(1.0 / p.b))
    end
end

function get_segment_volume(p::TaperProfile, h1::Float64, h2::Float64, type::String="total")
    hs = max(0.0, min(p.H, h1))
    he = max(0.0, min(p.H, h2))
    if hs >= he return 0.0 end
    
    # Segmentgrenze (hw) sauber im Metersystem überbrücken
    if hs < p.hw && he > p.hw
        return get_segment_volume(p, hs, p.hw, type) + get_segment_volume(p, p.hw, he, type)
    end
    
    is_basal = he <= p.hw
    
    if type == "total"
        if is_basal
            # Unteres Segment: d(h) = c0 + c1 * h^c2
          F_low(h) = p.c0^2 * h + (2.0 * p.c0 * p.c1 * h^(p.c2 + 1.0)) / (p.c2 + 1.0) + (p.c1^2 * h^(2.0 * p.c2 + 1.0)) / (2.0 * p.c2 + 1.0)
            return (pi / 40000.0) * (F_low(he) - F_low(hs))
        else
            # Oberes Segment: Direkt nach h integriert (Minuszeichen beachten!)
            # F_up(h) für d(h) = a * (H - h)^b + dt
            F_up(h) = p.dt^2 * h - (2.0 * p.a * p.dt * (p.H - h)^(p.b + 1.0)) / (p.b + 1.0) - (p.a^2 * (p.H - h)^(2.0 * p.b + 1.0)) / (2.0 * p.b + 1.0)
            return (pi / 40000.0) * (F_up(he) - F_up(hs))
        end
    elseif type == "wood"
        A = get_bark_coefficients()
        if is_basal
            sum_val_e = 0.0; sum_val_s = 0.0
            for n in 0:4
                for k in 0:n
                    bc = binomial_coeff(n, k)
                    sum_val_e += A[n + 1] * bc * p.c0^(n - k) * p.c1^k * (he^(k * p.c2 + 1.0)) / (k * p.c2 + 1.0)
                    sum_val_s += A[n + 1] * bc * p.c0^(n - k) * p.c1^k * (hs^(k * p.c2 + 1.0)) / (k * p.c2 + 1.0)
                end
            end
            return (pi / 40000.0) * (sum_val_e - sum_val_s)
        else
            # Oberes Segment für Holz: Integration über h (führt zu einem Minuszeichen vor dem inneren Integral)
            sum_val_e = 0.0; sum_val_s = 0.0
            for n in 0:4
                for k in 0:n
                    bc = binomial_coeff(n, k)
                    # Weil dh = -dz, wird der Term bei der Stammfunktion nach h negativiert!
                    term_e = p.dt^(n - k) * p.a^k * ((p.H - he)^(k * p.b + 1.0)) / (k * p.b + 1.0)
                    term_s = p.dt^(n - k) * p.a^k * ((p.H - hs)^(k * p.b + 1.0)) / (k * p.b + 1.0)
                    sum_val_e -= A[n + 1] * bc * term_e
                    sum_val_s -= A[n + 1] * bc * term_s
                end
            end
            return (pi / 40000.0) * (sum_val_e - sum_val_s)
        end
    end
    return 0.0
end

# =========================================================================
# 4. INTERFACES & OUTPUT
# =========================================================================

function calculate_assortments(h_meas::Vector{Float64}, d_meas::Vector{Float64}, H::Float64; h_stump=0.3, d_top=7.0)
    # BUGFIX: n_pts lokal definieren, damit es in dieser Funktion verfügbar ist!
    n_pts = length(h_meas)
    
    p = solve_taper(h_meas, d_meas, H)
    h_derb = get_height_at_diameter(p, d_top)
    if h_derb < h_stump h_derb = h_stump end
    
    # Brutto-Volumina (m3)
    v_stump_tot = get_segment_volume(p, 0.0, h_stump, "total")
    v_merc_tot  = get_segment_volume(p, h_stump, h_derb, "total")
    v_tip_tot   = get_segment_volume(p, h_derb, H, "total")
    v_total_tot = v_stump_tot + v_merc_tot + v_tip_tot
    
    # Netto-Volumina (Holz m3)
    v_stump_wood = get_segment_volume(p, 0.0, h_stump, "wood")
    v_merc_wood  = get_segment_volume(p, h_stump, h_derb, "wood")
    v_tip_wood   = get_segment_volume(p, h_derb, H, "wood")
    v_total_wood = v_stump_wood + v_merc_wood + v_tip_wood
    
    # BUGFIX: Sichere Textausgabe für die Konsole generieren
    dbh_print = n_pts > 0 ? "$(d_meas) cm" : "Kein Durchmesser angegeben (H/D=100 Fallback)"

    println("\n=======================================================")
    println("ERGEBNISSE JULIA-EINZELBAUM (H = $(H)m, Basis-DBH = $(dbh_print))")
    println("Derbholzgrenze ($(d_top)cm) erreicht bei h = $(round(h_derb, digits=2))m")
    println("=======================================================")
    
    arrangements = [
        ("Stock (Stump)", 0.0, h_stump, v_stump_tot, v_stump_wood, v_stump_tot-v_stump_wood),
        ("Derbholz (Merc)", h_stump, round(h_derb, digits=2), v_merc_tot, v_merc_wood, v_merc_tot-v_merc_wood),
        ("Wipfel (Tip)", round(h_derb, digits=2), H, v_tip_tot, v_tip_wood, v_tip_tot-v_tip_wood),
        ("Gesamt (Total)", 0.0, H, v_total_tot, v_total_wood, v_total_tot-v_total_wood)
    ]
    
    @printf("%-18s | %-7s | %-7s | %-12s | %-11s | %-11s\n", "Sektion", "Start_m", "Ende_m", "Vol_Gross_m3", "Vol_Wood_m3", "Vol_Bark_m3")
    println("--------------------------------------------------------------------------------")
    for r in arrangements
        @printf("%-18s | %-7.2f | %-7.2f | %-12.4f | %-11.4f | %-11.4f\n", r[1], r[2], r[3], r[4], r[5], r[6])
    end
    return p
end

function plot_tree(p::TaperProfile; aspect_ratio=:auto)
h_steps = range(0.0, p.H, length=300)
r_ob = [get_diameter_at_height(p, h) / 2.0 for h in h_steps]
r_ub = [get_diameter_under_bark(p, h) / 2.0 for h in h_steps]
x_wood = vcat(r_ub, reverse(-r_ub))
y_wood = vcat(h_steps, reverse(h_steps))
x_bark_r = vcat(r_ob, reverse(r_ub))
y_bark_r = vcat(h_steps, reverse(h_steps))
x_bark_l = vcat(-r_ob, reverse(-r_ub))
y_bark_l = vcat(h_steps, reverse(h_steps))
plt = plot(title = "Analytical Taper Profile (H = $(p.H)m)",
xaxis = "Radius [cm]", yaxis = "Height [m]",
size = (500, 850), left_margin = 6mm, bottom_margin = 4mm,
aspect_ratio = aspect_ratio)
plot!(plt, x_wood, y_wood, fill = (0, 0.6, :yellow), color = :orange, lw = 1, label = "Wood under Bark")
plot!(plt, x_bark_r, y_bark_r, fill = (0, 0.4, :brown), color = :brown, lw = 1, label = "Bark Layer")
plot!(plt, x_bark_l, y_bark_l, fill = (0, 0.4, :brown), color = :brown, lw = 1, label = "")
vline!(plt, [0.0], color = :black, linestyle = :dash, lw = 0.5, label = "")
display(plt)
return plt
end

# =========================================================================
# 5. TEST UMGEBUNG (JULIA) - EXAKT IDENTISCH ZU C++26
# =========================================================================

println("\n=======================================================")
println("[TEST 1] Pfad: N = 0 (Rein höhenbasierte Schätzung)")
println("=======================================================")
h0 = Float64[]
d0 = Float64[]
TP = calculate_assortments(h0, d0, 15.0, h_stump=0.3, d_top=7.0)
plot_tree(TP)
get_diameter_at_height(TP, 0.)
TP = solve_taper(h0, d0, 15.0, max_c0_factor = 1.)
plot_tree(TP)

println("\n=======================================================")
println("[TEST 2] Pfad: N = 1 (Punkt UNTER hw; z.B. BHD bei 1.3m)")
println("=======================================================")
h1_under = [1.3]
d1_under = [28.0]
TP = calculate_assortments(h1_under, d1_under, 20.0, h_stump=0.3, d_top=7.0)
plot_tree(TP)
get_diameter_at_height(TP, h1_under[1])
get_height_at_diameter(TP, d1_under[1])
TP = solve_taper(h1_under, d1_under, 20.0, max_c0_factor = 4.)
h3_std = [1.3, 7.0, 15.0]
map(x -> get_diameter_at_height(TP, x), h3_std)

println("\n=======================================================")
println("[TEST 3] Pfad: N = 1 (Punkt OBER hw; z.B. Kronenansatz)")
println("=======================================================")
h1_over = [12.0]
d1_over = [14.0]
TP = calculate_assortments(h1_over, d1_over, 18.0, h_stump=0.3, d_top=7.0)
plot_tree(TP)
get_diameter_at_height(TP, h1_over[1])
get_height_at_diameter(TP, d1_over[1])

println("\n=======================================================")
println("[TEST 4] Pfad: N = 2 (Beide Punkte OBER hw)")
println("=======================================================")
h2_over = [6.0, 14.0]
d2_over = [22.0, 15.0]
TP = calculate_assortments(h2_over, d2_over, 22.0, h_stump=0.3, d_top=7.0)
plot_tree(TP)
TP = solve_taper(h2_over, d2_over, 22.0, max_c0_factor = 4.)
TP = solve_taper(h2_over, d2_over, 22.0, max_c0_factor = 1.1)
map(x -> get_diameter_at_height(TP, x), h2_over)
map(x -> get_height_at_diameter(TP, x), d2_over)

println("\n=======================================================")
println("[TEST 5] Pfad: N = 2 (Gemischt: h1 UNTER hw, h2 OBER hw)")
println("=======================================================")
h2_mix = [1.3, 10.0]
d2_mix = [35.0, 27.0]
TP = calculate_assortments(h2_mix, d2_mix, 25.0, h_stump=0.3, d_top=7.0)
plot_tree(TP)
map(x -> get_diameter_at_height(TP, x), h2_mix)
map(x -> get_height_at_diameter(TP, x), d2_mix)

println("\n=======================================================")
println("[TEST 6] Pfad: N = 2 (Beide Punkte UNTER hw) -> Optimierung")
println("=======================================================")
h2_under = [0.5, 2.0]
d2_under = [48.0, 42.0]
TP = calculate_assortments(h2_under, d2_under, 30.0, h_stump=0.3, d_top=7.0)
plot_tree(TP)
map(x -> get_diameter_at_height(TP, x), h2_under)
map(x -> get_height_at_diameter(TP, x), d2_under)

println("\n=======================================================")
println("[TEST 7] Pfad: N >= 3 (Sektionierung mit plausiblen Daten)")
println("=======================================================")
h3_std = [1.3, 7.0, 15.0]
d3_std = [28.0, 22.0, 12.6]
#d3_std = [28.0, 22.0, 10.6]
TP = calculate_assortments(h3_std, d3_std, 20.0, h_stump=0.3, d_top=7.0)
plot_tree(TP)
TP = solve_taper(h3_std, d3_std, 20.0, max_c0_factor = 4.)
TP = solve_taper(h3_std, d3_std, 20.0, max_c0_factor = 1.1)
map(x -> get_diameter_at_height(TP, x), h3_std)
map(x -> get_height_at_diameter(TP, x), d3_std)

println("\n=======================================================")
println("[TEST 8] Pfad: N >= 3 (Anomalie: d3 >= d2 -> Fallback auf GLM)")
println("=======================================================")
h3_anomaly = [1.3, 8.0, 12.0]
d3_anomaly = [25.0, 16.0, 16.5] # d3 ist größer als d2
TP = calculate_assortments(h3_anomaly, d3_anomaly, 20.0, h_stump=0.3, d_top=7.0)
plot_tree(TP)
TP = solve_taper(h3_anomaly, d3_anomaly, 20.0, max_c0_factor = 4.)
map(x -> get_diameter_at_height(TP, x), h3_anomaly)   ## FEHLER
map(x -> get_height_at_diameter(TP, x), d3_anomaly)   ## FEHLER

# Alle Messhöhen liegen unterhalb von hw (2.0m)
println("[TEST 2] Pfad: SONDERFALL A (Alle Punkte UNTER hw)")
h = [0.3, 0.8, 1.5]
d = [35.0, 31.0, 28.0]
TP = calculate_assortments(h, d, 20.0, h_stump=0.3, d_top=7.0)
TP = solve_taper(h, d, 20.0, max_c0_factor = 1.1)
plot_tree(TP)
map(x -> get_diameter_at_height(TP, x), h)

# Alle Messhöhen liegen oberhalb von hw (2.0m)
println("[TEST 3] Pfad: SONDERFALL B (Alle Punkte OBERHALB hw)")
h = [6.0, 10.0, 16.0]
d = [21.0, 18.0, 11.0]
TP = calculate_assortments(h, d, 20.0, h_stump=0.3, d_top=7.0)
plot_tree(TP)
map(x -> get_diameter_at_height(TP, x), h)

# Extrem dicker unterer Punkt, der den mathematischen c0 weit über das Limit (30.0) treiben würde
println("[TEST 4] Pfad: LIMIT_C0 AKTIV (Erzwungene Stammfuß-Kappung)")
h = [0.5, 5.0, 12.0]
d = [45.0, 20.0, 14.0]
TP = calculate_assortments(h, d, 20.0, h_stump=0.3, d_top=7.0)
TP = solve_taper(h, d, 20.0, max_c0_factor = 4.)

plot_tree(TP)
map(x -> get_diameter_at_height(TP, x), h)

h = [1.3, 7.0, 15.0]
d = [28.0, 22.0, 12.6]
TP = calculate_assortments(h, d, 20.0, h_stump=0.3, d_top=7.0)
plot_tree(TP)
map(x -> get_diameter_at_height(TP, x), h)


h = [1.3]
d = [25.]
H = 25.
TP = solve_taper(h, d, H)
h = [1.3, 7.5, 16.25]
map(x -> get_diameter_at_height(TP, x), h)

d = [25., 20.8, 14.6]
TP = solve_taper(h, d, H)
plot_tree(TP)
scatter!(d ./ 2., h)

d = [25., 22.8, 14.6]
TP = solve_taper(h, d, H, h_ref=14.)
TP = solve_taper(h, d, H)
plot_tree(TP)
scatter!(d ./ 2., h)


h = [1.3, 7.0, 15.0]
d = [28.0, 22.0, 12.6]
TP = solve_taper(h, d, 20.)

d = [28.0, 22.0, 10.6]

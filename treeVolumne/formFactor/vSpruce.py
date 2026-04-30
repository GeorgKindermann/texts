import math

# Koeffizienten für die Fichte (Norway Spruce)
SPRUCE_PARAMS = {
    "basal_taper": {"b0": -0.2745, "b1": -37.6638, "b2": -34.7816},
    "stump_curv":  {"g0": -0.6293, "g1": -0.1842, "g2": 0.0260},
    "taper_curv":  {"b0": -0.5803, "b1": 0.0910,  "b2": 0.0399, "b3": -0.3286},
    "zeta_merc":   {"g0": 1.0571,  "g1": 0.4433,  "g2": -4.9835, "g3": 138.2230},
    "zeta_top":    {"d0": 0.8662,  "d1": 0.1141,  "d2": 1.2711,  "d3": -7.7670}
}

def calc_tree_volume(BHD, H_total, p):
    # --- 1. Referenzparameter ---
    d_t = 0.8 if H_total >= 1.3 else 0.1 + 0.7 * (1 - ((1.3 - H_total) / 1.3)**2)
    d_a = BHD + 0.65 if H_total >= 1.3 else d_t + H_total / 2
    h_ref = min(1.3, H_total)
    d_ref = BHD if H_total >= 1.3 else d_t
    
    # --- 2. Basale Geometrie ---
    h_st = min(H_total, 0.08 + 0.3 * (d_a / (30 + d_a)))
    s_0 = math.exp(p["basal_taper"]["b0"] + (p["basal_taper"]["b1"] / H_total) + (p["basal_taper"]["b2"] / d_a)) + 0.015
    d_0 = d_ref + s_0 * (h_ref * 100)
    e_st = math.exp(p["stump_curv"]["g0"] + p["stump_curv"]["g1"] * math.log(p["stump_curv"]["g2"] + s_0))
    d_st = d_0 - (h_st / h_ref)**e_st * (d_0 - d_ref)
    
    # --- 3. Stumpfvolumen (Analytisch) ---
    V_st = (math.pi * h_st / 40000) * (
        d_0**2 - (2 * d_0 * (d_0 - d_ref) / (e_st + 1)) * (h_st / h_ref)**e_st + 
        ((d_0 - d_ref)**2 / (2 * e_st + 1)) * (h_st / h_ref)**(2 * e_st)
    )
    
    # --- 4. Sektionslängen ---
    delta_d_st = max(0.01, d_st - 7)
    e_d7 = math.exp(p["taper_curv"]["b0"] + p["taper_curv"]["b1"] * (math.log(p["taper_curv"]["b2"] * delta_d_st))**2 + 
                    p["taper_curv"]["b3"] * math.log((H_total - h_st) / d_st))
    
    if d_st > 7:
        r_s = 1 - (7 - d_t) / (d_st - d_t)
        l_merc = (r_s**e_d7) * (H_total - h_st)
    else:
        l_merc = 0.0
    
    l_top = H_total - h_st - l_merc
    
    # --- 5. Dynamische Korrekturfaktoren (Zeta) ---
    zeta_m = p["zeta_merc"]["g0"] + (p["zeta_merc"]["g1"] / (1 + math.exp(p["zeta_merc"]["g2"] + p["zeta_merc"]["g3"] / max(0.1, l_merc))))
    zeta_t = p["zeta_top"]["d0"] + (p["zeta_top"]["d1"] / (1 + math.exp(p["zeta_top"]["d2"] * l_top + p["zeta_top"]["d3"] * (d_a / H_total))))
    
    # --- 6. Sektionsvolumen ---
    V_merc = 0.0
    if l_merc > 0:
        V_frust_m = (math.pi * l_merc / 120000) * (d_st**2 + d_st * 7 + 49)
        V_merc = V_frust_m / zeta_m
        
    d_up = 7.0 if l_merc > 0 else d_st
    V_frust_t = (math.pi * l_top / 120000) * (d_up**2 + d_up * d_t + d_t**2)
    V_top = V_frust_t / zeta_t
    
    return {
        "V_total": V_st + V_merc + V_top,
        "V_st": V_st,
        "V_merc": V_merc,
        "V_top": V_top,
        "h_st": h_st,
        "l_merc": l_merc,
        "l_top": l_top
    }

# Beispielaufruf


result = calc_tree_volume(30, 20, SPRUCE_PARAMS)

for key, value in result.items():
    print(f"{key:8}: {value:.6f}")

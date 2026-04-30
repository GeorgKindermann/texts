# Koeffizienten-Objekt definieren
spruce_params <- list(
  basal_taper = c(b0 = -0.2745, b1 = -37.6638, b2 = -34.7816),
  stump_curv  = c(g0 = -0.6293, g1 = -0.1842, g2 = 0.0260),
  taper_curv  = c(b0 = -0.5803, b1 = 0.0910,  b2 = 0.0399, b3 = -0.3286),
  zeta_merc   = c(g0 = 1.0571,  g1 = 0.4433,  g2 = -4.9835, g3 = 138.2230),
  zeta_top    = c(d0 = 0.8662,  d1 = 0.1141,  d2 = 1.2711, d3 = -7.7670)
)

calc_tree_volume <- function(BHD, H_total, p) {
  # 1. Referenzparameter
  d_t   <- ifelse(H_total >= 1.3, 0.8, 0.1 + 0.7 * (1 - ((1.3 - H_total) / 1.3)^2))
  d_a   <- ifelse(H_total >= 1.3, BHD + 0.65, d_t + H_total / 2)
  h_ref <- min(1.3, H_total)
  d_ref <- ifelse(H_total >= 1.3, BHD, d_t)
  
  # 2. Basale Geometrie
  h_st <- min(H_total, 0.08 + 0.3 * (d_a / (30 + d_a)))
  s_0  <- exp(p$basal_taper["b0"] + (p$basal_taper["b1"] / H_total) + (p$basal_taper["b2"] / d_a)) + 0.015
  d_0  <- d_ref + s_0 * (h_ref * 100)
  e_st <- exp(p$stump_curv["g0"] + p$stump_curv["g1"] * log(p$stump_curv["g2"] + s_0))
  d_st <- d_0 - (h_st / h_ref)^e_st * (d_0 - d_ref)
  
  # 3. Stumpfvolumen
  V_st <- (pi * h_st / 40000) * (
    d_0^2 - (2 * d_0 * (d_0 - d_ref) / (e_st + 1)) * (h_st / h_ref)^e_st + 
    ((d_0 - d_ref)^2 / (2 * e_st + 1)) * (h_st / h_ref)^(2 * e_st)
  )
  
  # 4. Sektionslängen
  delta_d_st <- max(0.01, d_st - 7)
  e_d7 <- exp(p$taper_curv["b0"] + p$taper_curv["b1"] * (log(p$taper_curv["b2"] * delta_d_st))^2 + 
              p$taper_curv["b3"] * log((H_total - h_st) / d_st))
  
  l_merc <- if(d_st > 7) ( (1 - (7 - d_t) / (d_st - d_t))^e_d7 ) * (H_total - h_st) else 0
  l_top  <- H_total - h_st - l_merc
  
  # 5. Zeta & Volumen
  z_m <- p$zeta_merc["g0"] + (p$zeta_merc["g1"] / (1 + exp(p$zeta_merc["g2"] + p$zeta_merc["g3"] / max(0.1, l_merc))))
  z_t <- p$zeta_top["d0"] + (p$zeta_top["d1"] / (1 + exp(p$zeta_top["d2"] * l_top + p$zeta_top["d3"] * (d_a / H_total))))
  
  V_m <- if(l_merc > 0) ((pi * l_merc / 120000) * (d_st^2 + d_st * 7 + 49)) / z_m else 0
  d_up <- if(l_merc > 0) 7 else d_st
  V_t <- ((pi * l_top / 120000) * (d_up^2 + d_up * d_t + d_t^2)) / z_t
  
  #return(list(V_total = V_st + V_m + V_t, V_st = V_st, V_merc = V_m, V_top = V_t))
  return(list(V_total = V_st + V_m + V_t, V_st = V_st, V_merc = V_m, V_top = V_t, h_st = h_st, l_merc = l_merc, l_top = l_top))
}

# Aufruf:
calc_tree_volume(30, 20, spruce_params)

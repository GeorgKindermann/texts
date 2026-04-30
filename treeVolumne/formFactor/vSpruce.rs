use std::f64::consts::PI;

// Struktur für die Parameter (ähnlich einer Liste in R)
struct ModelParams {
    basal_taper: [f64; 3], // b0, b1, b2
    stump_curv:  [f64; 3], // g0, g1, g2
    taper_curv:  [f64; 4], // b0, b1, b2, b3
    zeta_merc:   [f64; 4], // g0, g1, g2, g3
    zeta_top:    [f64; 4], // d0, d1, d2, d3
}

// Die Ergebnisse als Struktur
#[derive(Debug)]
struct TreeResult {
    v_total: f64,
    v_st: f64,
    v_merc: f64,
    v_top: f64,
    h_st: f64,
    l_merc: f64,
    l_top: f64,
}

fn calc_tree_volume(bhd: f64, h_total: f64, p: &ModelParams) -> TreeResult {
    // --- 1. Referenzparameter ---
    let d_t = if h_total >= 1.3 { 0.8 } else { 0.1 + 0.7 * (1.0 - ((1.3 - h_total) / 1.3).powi(2)) };
    let d_a = if h_total >= 1.3 { bhd + 0.65 } else { d_t + h_total / 2.0 };
    let h_ref = h_total.min(1.3);
    let d_ref = if h_total >= 1.3 { bhd } else { d_t };

    // --- 2. Basale Geometrie ---
    let h_st = h_total.min(0.08 + 0.3 * (d_a / (30.0 + d_a)));
    let s_0 = (p.basal_taper[0] + (p.basal_taper[1] / h_total) + (p.basal_taper[2] / d_a)).exp() + 0.015;
    let d_0 = d_ref + s_0 * (h_ref * 100.0);
    let e_st = (p.stump_curv[0] + p.stump_curv[1] * (p.stump_curv[2] + s_0).ln()).exp();
    let d_st = d_0 - (h_st / h_ref).powf(e_st) * (d_0 - d_ref);

    // --- 3. Stumpfvolumen (Analytisch) ---
    let v_st = (PI * h_st / 40000.0) * (
        d_0.powi(2) - 
        (2.0 * d_0 * (d_0 - d_ref) / (e_st + 1.0)) * (h_st / h_ref).powf(e_st) + 
        ((d_0 - d_ref).powi(2) / (2.0 * e_st + 1.0)) * (h_st / h_ref).powf(2.0 * e_st)
    );

    // --- 4. Sektionslängen ---
    let delta_d_st = 0.01_f64.max(d_st - 7.0);
    let e_d7 = (p.taper_curv[0] + p.taper_curv[1] * (p.taper_curv[2] * delta_d_st).ln().powi(2) + 
                p.taper_curv[3] * ((h_total - h_st) / d_st).ln()).exp();

    let l_merc = if d_st > 7.0 {
        let r_s = 1.0 - (7.0 - d_t) / (d_st - d_t);
        r_s.powf(e_d7) * (h_total - h_st)
    } else {
        0.0
    };
    let l_top = h_total - h_st - l_merc;

    // --- 5. Dynamische Korrekturfaktoren (Zeta) ---
    let zeta_m = p.zeta_merc[0] + (p.zeta_merc[1] / (1.0 + (p.zeta_merc[2] + p.zeta_merc[3] / 0.1_f64.max(l_merc)).exp()));
    let zeta_t = p.zeta_top[0] + (p.zeta_top[1] / (1.0 + (p.zeta_top[2] * l_top + p.zeta_top[3] * (d_a / h_total)).exp()));

    // --- 6. Sektionsvolumen ---
    let mut v_merc = 0.0;
    if l_merc > 0.0 {
        let v_frust_m = (PI * l_merc / 120000.0) * (d_st.powi(2) + d_st * 7.0 + 49.0);
        v_merc = v_frust_m / zeta_m;
    }

    let d_up = if l_merc > 0.0 { 7.0 } else { d_st };
    let v_frust_t = (PI * l_top / 120000.0) * (d_up.powi(2) + d_up * d_t + d_t.powi(2));
    let v_top = v_frust_t / zeta_t;

    TreeResult {
        v_total: v_st + v_merc + v_top,
        v_st,
        v_merc,
        v_top,
        h_st,
        l_merc,
        l_top,
    }
}

fn main() {
    let spruce = ModelParams {
        basal_taper: [-0.2745, -37.6638, -34.7816],
        stump_curv:  [-0.6293, -0.1842, 0.0260],
        taper_curv:  [-0.5803, 0.0910, 0.0399, -0.3286],
        zeta_merc:   [1.0571, 0.4433, -4.9835, 138.2230],
        zeta_top:    [0.8662, 0.1141, 1.2711, -7.7670],
    };

    let result = calc_tree_volume(30.0, 20.0, &spruce);
    
    println!("{:#?}", result);
}

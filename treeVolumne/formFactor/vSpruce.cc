#include <iostream>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <vector>

// Struktur für die Koeffizienten (Wartbarkeit)
struct ModelParams {
    std::vector<double> basal_taper; // b0, b1, b2
    std::vector<double> stump_curv;  // g0, g1, g2
    std::vector<double> taper_curv;  // b0, b1, b2, b3
    std::vector<double> zeta_merc;   // g0, g1, g2, g3
    std::vector<double> zeta_top;    // d0, d1, d2, d3
};

// Struktur für das detaillierte Ergebnis (entspricht deiner R-Liste)
struct TreeResult {
    double V_total, V_st, V_merc, V_top;
    double h_st, l_merc, l_top;
};

// Die globale Parameter-Definition für Fichte
const ModelParams SPRUCE_PARAMS = {
    {-0.2745, -37.6638, -34.7816},          // basal_taper
    {-0.6293, -0.1842, 0.0260},             // stump_curv
    {-0.5803, 0.0910, 0.0399, -0.3286},     // taper_curv
    {1.0571, 0.4433, -4.9835, 138.2230},    // zeta_merc
    {0.8662, 0.1141, 1.2711, -7.7670}       // zeta_top
};

TreeResult calc_tree_volume(double BHD, double H_total, const ModelParams& p) {
    const double PI = std::acos(-1.0);
    TreeResult res;

    // --- 1. Referenzparameter ---
    double d_t = (H_total >= 1.3) ? 0.8 : 0.1 + 0.7 * (1.0 - std::pow((1.3 - H_total) / 1.3, 2));
    double d_a = (H_total >= 1.3) ? BHD + 0.65 : d_t + H_total / 2.0;
    double h_ref = std::min(1.3, H_total);
    double d_ref = (H_total >= 1.3) ? BHD : d_t;

    // --- 2. Basale Geometrie ---
    res.h_st = std::min(H_total, 0.08 + 0.3 * (d_a / (30.0 + d_a)));
    double s_0 = std::exp(p.basal_taper[0] + (p.basal_taper[1] / H_total) + (p.basal_taper[2] / d_a)) + 0.015;
    double d_0 = d_ref + s_0 * (h_ref * 100.0);
    double e_st = std::exp(p.stump_curv[0] + p.stump_curv[1] * std::log(p.stump_curv[2] + s_0));
    double d_st = d_0 - std::pow(res.h_st / h_ref, e_st) * (d_0 - d_ref);

    // --- 3. Stumpfvolumen (V_st) ---
    res.V_st = (PI * res.h_st / 40000.0) * (
        std::pow(d_0, 2) - 
        (2.0 * d_0 * (d_0 - d_ref) / (e_st + 1.0)) * std::pow(res.h_st / h_ref, e_st) + 
        (std::pow(d_0 - d_ref, 2) / (2.0 * e_st + 1.0)) * std::pow(res.h_st / h_ref, 2.0 * e_st)
    );

    // --- 4. Sektionslängen (l_merc, l_top) ---
    double delta_d_st = std::max(0.01, d_st - 7.0);
    double e_d7 = std::exp(p.taper_curv[0] + p.taper_curv[1] * std::pow(std::log(p.taper_curv[2] * delta_d_st), 2) + 
                           p.taper_curv[3] * std::log((H_total - res.h_st) / d_st));

    if (d_st > 7.0) {
        double r_s = 1.0 - (7.0 - d_t) / (d_st - d_t);
        res.l_merc = std::pow(r_s, e_d7) * (H_total - res.h_st);
    } else {
        res.l_merc = 0.0;
    }
    res.l_top = H_total - res.h_st - res.l_merc;

    // --- 5. Zeta & Volumen-Korrektur ---
    double z_m = p.zeta_merc[0] + (p.zeta_merc[1] / (1.0 + std::exp(p.zeta_merc[2] + p.zeta_merc[3] / std::max(0.1, res.l_merc))));
    double z_t = p.zeta_top[0] + (p.zeta_top[1] / (1.0 + std::exp(p.zeta_top[2] * res.l_top + p.zeta_top[3] * (d_a / H_total))));

    res.V_merc = 0.0;
    if (res.l_merc > 0) {
        double V_frust_m = (PI * res.l_merc / 120000.0) * (std::pow(d_st, 2) + d_st * 7.0 + 49.0);
        res.V_merc = V_frust_m / z_m;
    }

    double d_up = (res.l_merc > 0) ? 7.0 : d_st;
    double V_frust_t = (PI * res.l_top / 120000.0) * (std::pow(d_up, 2) + d_up * d_t + std::pow(d_t, 2));
    res.V_top = V_frust_t / z_t;

    res.V_total = res.V_st + res.V_merc + res.V_top;
    return res;
}

int main() {
    TreeResult r = calc_tree_volume(30.0, 20.0, SPRUCE_PARAMS);

    std::cout << std::fixed << std::setprecision(5);
    std::cout << "V_total: " << r.V_total << "\nV_st:    " << r.V_st << "\nV_merc:  " << r.V_merc << "\nV_top:   " << r.V_top << std::endl;
    std::cout << "h_st:    " << r.h_st << "\nl_merc:  " << r.l_merc << "\nl_top:   " << r.l_top << std::endl;

    return 0;
}

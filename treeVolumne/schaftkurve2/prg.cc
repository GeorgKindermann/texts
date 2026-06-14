// /opt/gcc-16/bin/g++-16 prg.cc -std=c++26 -Wall -Wextra -Wconversion -Wl,-rpath,/opt/gcc-16/lib64

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include <string_view>
#include <format>
#include <stdexcept>
#include <print>
#include <optional>

// =========================================================================
// 1. FORSTLICHE DATENSTRUKTUR (D in cm, H in m)
// =========================================================================
struct TaperProfile {
    double H;   // H_total: Gesamthöhe (m)
    double hw;  // H_ip: Inflexionspunkt (m)
    double c0;  // d_0: Durchmesser am Boden (cm)
    double dw;  // d_ip: Durchmesser am Inflexionspunkt (cm)
    double dt;  // d_t: Terminaler Durchmesser an der Spitze (cm)
    double a;   // Parameter a (oberes Segment)
    double b;   // Parameter b (oberes Segment)
    double c1;  // Parameter c1 (unteres Segment, negativ)
    double c2;  // Parameter c2 (unteres Segment, positiv)
};

// =========================================================================
// Typ-Aliase, um Markdown-Darstellungsfehler mit spitzen Klammern zu umgehen
// =========================================================================
using DoubleVector = std::vector<double>;

struct BarkCoefficients {
    double A0, A1, A2, A3, A4;
};

inline BarkCoefficients get_bark_coefficients() noexcept {
    constexpr double alpha = 0.85149 / 10.0;
    constexpr double beta  = 0.60934 / 10.0;
    constexpr double gamma = -0.00228 / 10.0;
    return BarkCoefficients{
        alpha * alpha,
        -2.0 * alpha * (1.0 - beta),
        (1.0 - beta) * (1.0 - beta) + 2.0 * alpha * gamma,
        -2.0 * gamma * (1.0 - beta),
        gamma * gamma
    };
}

inline double inline_plogis(double x) noexcept {
    return 1.0 / (1.0 + std::exp(-x));
}

inline double binomial_coeff(int n, int k) noexcept {
    if (k < 0 || k > n) return 0.0;
    double res = 1.0;
    for (int i = 1; i <= k; ++i) {
        res *= static_cast<double>(n - i + 1) / static_cast<double>(i);
    }
    return res;
}

// Bisektion für Nullstellensuche (Äquivalent zu Roots.jl Bisection)
template <typename F>
double find_zero_bisection(F&& f, double lower, double upper, double tol = 1e-7, int max_iter = 100) {
    double f_lower = f(lower);
    double f_upper = f(upper);
    if (f_lower * f_upper > 0.0) return 0.5 * (lower + upper);
    for (int i = 0; i < max_iter; ++i) {
        double mid = lower + 0.5 * (upper - lower);
        double f_mid = f(mid);
        if (std::abs(f_mid) < tol || (0.5 * (upper - lower)) < tol) return mid;
        if ((f_lower > 0.0) == (f_mid > 0.0)) { lower = mid; f_lower = f_mid; } else { upper = mid; }
    }
    return lower + 0.5 * (upper - lower);
}

// Sekantenverfahren für offene Nullstellensuche (Äquivalent zu Roots.jl Order0 / Secant)
template <typename F>
double find_zero_secant(F&& f, double x0, double tol = 1e-7, int max_iter = 100) {
    double x1 = x0 * 1.05 + 0.1; // Zweiter Startpunkt nahe dem ersten
    double y0 = f(x0);
    double y1 = f(x1);
    
    for (int i = 0; i < max_iter; ++i) {
        if (std::abs(y1 - y0) < 1e-12) break;
        double x_next = x1 - y1 * (x1 - x0) / (y1 - y0);
        if (std::abs(x_next - x1) < tol) return x_next;
        
        x0 = x1; y0 = y1;
        x1 = x_next; y1 = f(x1);
    }
    return x1;
}

// Brent-Verfahren für Minimierung (Äquivalent zu Optim.jl)
template <typename F>
double optimize_brent(F&& f, double ax, double cx, double tol = 1e-8, int max_iter = 100) {
    constexpr double cgold = 0.3819660112501051;
    
    // Mathematisch korrekte, asymmetrische Initialisierung des mittleren Punktes bx
    double bx = ax + cgold * (cx - ax);
    
    double a = std::min(ax, cx), b = std::max(ax, cx);
    double v = bx, w = bx, x = bx, d = 0.0, e = 0.0;
    double fx = f(x), fv = fx, fw = fx;
    
    for (int iter = 0; iter < max_iter; ++iter) {
        double xm = 0.5 * (a + b);
        double tol1 = tol * std::abs(x) + 1e-10;
        double tol2 = 2.0 * tol1;
        
        if (std::abs(x - xm) <= (tol2 - 0.5 * (b - a))) return x;
        
        if (std::abs(e) > tol1) {
            double r = (x - w) * (fx - fv);
            double q = (x - v) * (fx - fw);
            double p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if (q > 0.0) p = -p;
            q = std::abs(q);
            double etemp = e; e = d;
            
            if (std::abs(p) >= std::abs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x)) {
                e = (x >= xm ? a - x : b - x);
                d = cgold * e;
            } else {
                d = p / q;
                double u = x + d;
                if (u - a < tol2 || b - u < tol2) {
                    d = (xm - x >= 0.0 ? tol1 : -tol1);
                }
            }
        } else {
            e = (x >= xm ? a - x : b - x);
            d = cgold * e;
        }
        
        double u = (std::abs(d) >= tol1 ? x + d : x + (d >= 0.0 ? tol1 : -tol1));
        double fu = f(u);
        
        if (fu <= fx) {
            if (u >= x) a = x; else b = x;
            v = w; fv = fw; w = x; fw = fx; x = u; fx = fu;
        } else {
            if (u < x) a = u; else b = u;
            if (fu <= fw || w == x) { v = w; fv = fw; w = u; fw = fu; }
            else if (fu <= fv || v == x || v == w) { v = u; fv = fu; }
        }
    }
    return x;
}

template <typename F>
double optimize_brentOLD(F&& f, double ax, double cx, double tol = 1e-7, int max_iter = 100) {
    constexpr double cgold = 0.3819660;
    double bx = 0.5 * (ax + cx);
    double a = std::min(ax, cx), b = std::max(ax, cx);
    double v = bx, w = bx, x = bx, d = 0.0, e = 0.0;
    double fx = f(x), fv = fx, fw = fx;
    for (int iter = 0; iter < max_iter; ++iter) {
        double xm = 0.5 * (a + b);
        double tol1 = tol * std::abs(x) + 1e-10;
        double tol2 = 2.0 * tol1;
        if (std::abs(x - xm) <= (tol2 - 0.5 * (b - a))) return x;
        if (std::abs(e) > tol1) {
            double r = (x - w) * (fx - fv), q = (x - v) * (fx - fw);
            double p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r); if (q > 0.0) p = -p; q = std::abs(q);
            double etemp = e; e = d;
            if (std::abs(p) >= std::abs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x)) {
                e = (x >= xm ? a - x : b - x); d = cgold * e;
            } else {
                d = p / q; double u = x + d;
                if (u - a < tol2 || b - u < tol2) d = (xm - x >= 0.0 ? std::abs(tol1) : -std::abs(tol1));
            }
        } else { e = (x >= xm ? a - x : b - x); d = cgold * e; }
        double u = (std::abs(d) >= tol1 ? x + d : x + (d >= 0.0 ? std::abs(tol1) : -std::abs(tol1)));
        double fu = f(u);
        if (fu <= fx) {
            if (u >= x) a = x; else b = x;
            v = w; fv = fw; w = x; fw = fx; x = u; fx = fu;
        } else {
            if (u < x) a = u; else b = u;
            if (fu <= fw || w == x) { v = w; fv = fw; w = u; fw = fu; }
            else if (fu <= fv || v == x || v == w) { v = u; fv = fu; }
        }
    }
    return x;
}

// Hilfsfunktion zum Sortieren zweier Vektoren nach der Höhe h_meas
inline void sort_by_height(DoubleVector& h, DoubleVector& d) noexcept {
    if (h.empty()) return;
    std::vector<size_t> p(h.size());
    std::iota(p.begin(), p.end(), 0);
    std::sort(p.begin(), p.end(), [&](size_t i, size_t j) { return h[i] < h[j]; });
    DoubleVector h_sorted(h.size()), d_sorted(d.size());
    for (size_t i = 0; i < h.size(); ++i) {
        h_sorted[i] = h[p[i]];
        d_sorted[i] = d[p[i]];
    }
    h = std::move(h_sorted);
    d = std::move(d_sorted);
}

// Kompakter 2D Nelder-Mead Optimierer mit Box-Constraints (Schranken) für den 3-Punkt-Fall
// 2D Nelder-Mead Hilfsoperatoren für std::vector
/*
inline std::vector<double> operator+(const std::vector<double>& a, const std::vector<double>& b) {
    return { a[0] + b[0], a[1] + b[1] };
}
inline std::vector<double> operator-(const std::vector<double>& a, const std::vector<double>& b) {
    return { a[0] - b[0], a[1] - b[1] };
}
inline std::vector<double> operator*(double scalar, const std::vector<double>& v) {
    return { scalar * v[0], scalar * v[1] };
}
*/

template <typename F>
std::vector<double> optimize_nelder_mead_2d(F&& f, const std::vector<double>& lb, const std::vector<double>& ub, const std::vector<double>& init, int max_iter = 500) {
    struct Point2D {
        double b, c2;
        Point2D operator+(const Point2D& o) const { return {b + o.b, c2 + o.c2}; }
        Point2D operator-(const Point2D& o) const { return {b - o.b, c2 - o.c2}; }
        Point2D operator*(double s) const { return {b * s, c2 * s}; }
        std::vector<double> to_vec() const { return {b, c2}; }
    };

    // Weiche Straffunktion: Erlaubt dem Simplex das Entlanggleiten an den Rändern (wie Julias Fminbox)
    auto penalty_f = [&](const Point2D& p) -> double {
        double penalty = 0.0;
        if (p.b < lb[0])  penalty += std::pow(lb[0] - p.b, 2) * 1e6;
        if (p.b > ub[0])  penalty += std::pow(p.b - ub[0], 2) * 1e6;
        if (p.c2 < lb[1]) penalty += std::pow(lb[1] - p.c2, 2) * 1e6;
        if (p.c2 > ub[1]) penalty += std::pow(p.c2 - ub[1], 2) * 1e6;
        
        // Parameter clambern für die eigentliche Funktionsauswertung
        double clamped_b = std::clamp(p.b, lb[0], ub[0]);
        double clamped_c2 = std::clamp(p.c2, lb[1], ub[1]);
        return f({clamped_b, clamped_c2}) + penalty;
    };

    // Initialisierung mit an den Suchraum angepassten, feinen Schrittweiten
    Point2D p0 = {init[0], init[1]};
    Point2D p1 = {init[0] + 0.05, init[1]};
    Point2D p2 = {init[0], init[1] + 0.05};
    //Point2D p1 = {init[0] + 0.005, init[1]};
    //Point2D p2 = {init[0], init[1] + 0.005};
    
    std::vector<Point2D> s = {p0, p1, p2};
    std::vector<double> f_val = {penalty_f(s[0]), penalty_f(s[1]), penalty_f(s[2])};

    for (int iter = 0; iter < max_iter; ++iter) {
        std::vector<int> idx = {0, 1, 2};
        std::sort(idx.begin(), idx.end(), [&](int i, int j) { return f_val[i] < f_val[j]; });
        int b_i = idx[0], w_i = idx[2], g_i = idx[1];
        
        Point2D cen = (s[b_i] + s[g_i]) * 0.5;
        
        // 1. Reflektion
        Point2D r = cen + (cen - s[w_i]);
        double fr = penalty_f(r);
        if (fr < f_val[g_i] && fr >= f_val[b_i]) { s[w_i] = r; f_val[w_i] = fr; continue; }
        
        // 2. Expansion
        if (fr < f_val[b_i]) {
            Point2D e = cen + (r - cen) * 2.0;
            double fe = penalty_f(e);
            if (fe < fr) { s[w_i] = e; f_val[w_i] = fe; } else { s[w_i] = r; f_val[w_i] = fr; }
            continue;
        }
        
        // 3. Kontraktion
        Point2D c = cen + (s[w_i] - cen) * 0.5;
        double fc = penalty_f(c);
        if (fc < f_val[w_i]) { s[w_i] = c; f_val[w_i] = fc; continue; }
        
        // 4. Shrink
        for (int i = 0; i < 3; ++i) {
            if (i != b_i) {
                s[i] = s[b_i] + (s[i] - s[b_i]) * 0.5;
                f_val[i] = penalty_f(s[i]);
            }
        }
    }
    
    std::vector<int> final_idx = {0, 1, 2};
    std::sort(final_idx.begin(), final_idx.end(), [&](int i, int j) { return f_val[i] < f_val[j]; });
    
    // Ergebnis innerhalb der harten Schranken zurückgeben
    double final_b = std::clamp(s[final_idx[0]].b, lb[0], ub[0]);
    double final_c2 = std::clamp(s[final_idx[0]].c2, lb[1], ub[1]);
    return {final_b, final_c2};
}

template <typename F>
std::vector<double> optimize_nelder_mead_2dOLD(F&& f, const std::vector<double>& lb, const std::vector<double>& ub, const std::vector<double>& init, int max_iter = 300) {
    // Wrapper-Funktion, die harte Schranken über eine unendliche Straffunktion erzwingt
    auto penalty_f = [&](const std::vector<double>& x) -> double {
        if (x[0] < lb[0] || x[0] > ub[0] || x[1] < lb[1] || x[1] > ub[1]) {
            return std::numeric_limits<double>::max(); // Verhindert das Festfahren an den Ecken
        }
        return f(x);
    };

    // Erzeugung eines stabilen, größeren Start-Simplex
    std::vector<std::vector<double>> s(3, std::vector<double>(2));
    s[0] = init;
    s[1] = {init[0] + 0.10, init[1]};
    s[2] = {init[0], init[1] + 0.10};
    
    std::vector<double> f_val = {penalty_f(s[0]), penalty_f(s[1]), penalty_f(s[2])};

    for (int iter = 0; iter < max_iter; ++iter) {
        std::vector<int> idx = {0, 1, 2};
        std::sort(idx.begin(), idx.end(), [&](int i, int j) { return f_val[i] < f_val[j]; });
        int b_i = idx[0], w_i = idx[2], g_i = idx[1];
        
        std::vector<double> cen = {0.5 * (s[b_i][0] + s[g_i][0]), 0.5 * (s[b_i][1] + s[g_i][1])};
        
        // 1. Reflektion
        std::vector<double> r = {2.0 * cen[0] - s[w_i][0], 2.0 * cen[1] - s[w_i][1]};
        double fr = penalty_f(r);
        if (fr < f_val[g_i] && fr >= f_val[b_i]) { s[w_i] = r; f_val[w_i] = fr; continue; }
        
        // 2. Expansion
        if (fr < f_val[b_i]) {
            std::vector<double> e = {cen[0] + 2.0 * (r[0] - cen[0]), cen[1] + 2.0 * (r[1] - cen[1])};
            double fe = penalty_f(e); s[w_i] = (fe < fr) ? e : r; f_val[w_i] = (fe < fr) ? fe : fr; continue;
        }
        
        // 3. Kontraktion
        std::vector<double> c = {cen[0] + 0.5 * (s[w_i][0] - cen[0]), cen[1] + 0.5 * (s[w_i][1] - cen[1])};
        double fc = penalty_f(c);
        if (fc < f_val[w_i]) { s[w_i] = c; f_val[w_i] = fc; continue; }
        
        // 4. Shrink
        for (int i = 0; i < 3; ++i) {
            if (i != b_i) {
                s[i][0] = 0.5 * (s[i][0] + s[b_i][0]);
                s[i][1] = 0.5 * (s[i][1] + s[b_i][1]);
                f_val[i] = penalty_f(s[i]);
            }
        }
    }
    std::vector<int> final_idx = {0, 1, 2};
    std::sort(final_idx.begin(), final_idx.end(), [&](int i, int j) { return f_val[i] < f_val[j]; });
    return s[final_idx[0]];
}


// =========================================================================
// 2. INTERNE PROFILBERECHNUNG BEI FESTEM dr
// =========================================================================
TaperProfile _solve_taper_with_dr(DoubleVector h_meas, DoubleVector d_meas, double H, double dt, double hw, double dr, double max_c0_factor, bool& out_is_valid, std::optional<double> h_ref = std::nullopt) {
    size_t n_pts = h_meas.size();
    sort_by_height(h_meas, d_meas);
    
    double limit_c0 = max_c0_factor * dr;
    double b_glm  = std::clamp(inline_plogis(-4.31616 + 21.7497 * (1.0 / std::max(0.1, H)) + 1.0813 * std::log1p(dr)), 0.01, 1.0);
    double c2_glm = std::clamp(inline_plogis(7.0974 - 3.17164 * std::log1p(std::max(0.1, H)) + 0.631264 * std::log1p(dr)), 0.01, 1.0);
    
    double local_b   = b_glm;
    double local_c2  = c2_glm;
    double local_a   = 0.0, local_dip = 0.0, sip = 0.0, local_c1 = 0.0, local_c0 = 0.0;

    if (n_pts == 0) {
        local_c0 = std::min(dt + H, limit_c0);
        auto f_zero_n0 = [&](double a_trial) noexcept -> double {
            double dw_t = a_trial * std::pow(H - hw, local_b) + dt;
            double sw_t = -a_trial * local_b * std::pow(H - hw, local_b - 1.0);
            double c1_t = sw_t / (local_c2 * std::pow(hw, local_c2 - 1.0));
            return (dw_t - c1_t * std::pow(hw, local_c2)) - local_c0;
        };
        local_a = find_zero_bisection(f_zero_n0, 0.0001, local_c0 * 2.0);
        local_dip = local_a * std::pow(H - hw, local_b) + dt;
        sip = -local_a * local_b * std::pow(H - hw, local_b - 1.0);
        local_c1 = sip / (local_c2 * std::pow(hw, local_c2 - 1.0));

    } else if (n_pts == 1) {
        double h1 = h_meas[0], d1 = d_meas[0];
        if (h1 <= hw) {
            double Omega = (local_b * hw) / (local_c2 * (H - hw)) * (1.0 - std::pow(h1 / hw, local_c2));
            double trial_dip = (d1 + dt * Omega) / (1.0 + Omega);
            double trial_a = (trial_dip - dt) / std::pow(H - hw, local_b);
            double trial_sip = -trial_a * local_b * std::pow(H - hw, local_b - 1.0);
            double trial_c1 = trial_sip / (local_c2 * std::pow(hw, local_c2 - 1.0));
            double trial_c0 = trial_dip - trial_c1 * std::pow(hw, local_c2);
            
            if (trial_c0 > limit_c0) {
                local_c0 = limit_c0;
                auto f_zero_c2 = [&](double c2_trial) noexcept -> double {
                    double Omega_loop = (local_b * hw) / (c2_trial * (H - hw)) * (1.0 - std::pow(h1 / hw, c2_trial));
                    double tdip = (d1 + dt * Omega_loop) / (1.0 + Omega_loop);
                    double ta = (tdip - dt) / std::pow(H - hw, local_b);
                    double tsip = -ta * local_b * std::pow(H - hw, local_b - 1.0);
                    return tdip - (tsip / (c2_trial * std::pow(hw, c2_trial - 1.0))) * std::pow(hw, c2_trial) - limit_c0;
                };
                local_c2 = find_zero_bisection(f_zero_c2, 0.01, 1.0);
                Omega = (local_b * hw) / (local_c2 * (H - hw)) * (1.0 - std::pow(h1 / hw, local_c2));
                local_dip = (d1 + dt * Omega) / (1.0 + Omega);
                local_a = (local_dip - dt) / std::pow(H - hw, local_b);
                sip = -local_a * local_b * std::pow(H - hw, local_b - 1.0);
                local_c1 = sip / (local_c2 * std::pow(hw, local_c2 - 1.0));
            } else {
                local_dip = trial_dip; local_a = trial_a; sip = trial_sip; local_c1 = trial_c1; local_c0 = trial_c0;
            }
        } else {
            local_a = (d1 - dt) / std::pow(H - h1, local_b);
            local_dip = local_a * std::pow(H - hw, local_b) + dt;
            sip = -local_a * local_b * std::pow(H - hw, local_b - 1.0);
            local_c1 = sip / (local_c2 * std::pow(hw, local_c2 - 1.0));
            local_c0 = local_dip - local_c1 * std::pow(hw, local_c2);
            if (local_c0 > limit_c0) {
                local_c0 = limit_c0;
                local_c2 = std::clamp((sip * hw) / (local_dip - local_c0), 0.01, 1.0);
                local_c1 = (local_dip - local_c0) / std::pow(hw, local_c2);
            }
        }

    } else if (n_pts == 2) {
        double h1 = h_meas[0], h2 = h_meas[1];
        double d1 = d_meas[0], d2 = d_meas[1];
        
        if (h1 >= hw && h2 >= hw) {
            double term_num = std::max(1e-6, d2 - dt) / std::max(1e-6, d1 - dt);
            double term_den = std::max(1e-6, H - h2) / std::max(1e-6, H - h1);
            local_b = std::clamp(std::log(term_num) / std::log(term_den), 0.2, 1.0);
            local_a = (d1 - dt) / std::pow(std::max(1e-6, H - h1), local_b);
            local_dip = local_a * std::pow(H - hw, local_b) + dt;
            sip = -local_a * local_b * std::pow(H - hw, local_b - 1.0);
            local_c1 = sip / (local_c2 * std::pow(hw, local_c2 - 1.0));
            local_c0 = std::min(local_dip - local_c1 * std::pow(hw, local_c2), limit_c0);
        } else if (h1 < hw && h2 >= hw) {
            auto find_b_trajectory_error = [&](double b_trial) noexcept -> double {
                double a_t = (d2 - dt) / std::pow(std::max(1e-6, H - h2), b_trial);
                double dw_t = a_t * std::pow(H - hw, b_trial) + dt;
                double sw_t = -a_t * b_trial * std::pow(H - hw, b_trial - 1.0);
                double c1_t = sw_t / (local_c2 * std::pow(hw, local_c2 - 1.0));
                double c0_t = dw_t - c1_t * std::pow(hw, local_c2);
                return (c0_t + c1_t * std::pow(h1, local_c2)) - d1;
            };
            if ((find_b_trajectory_error(0.01) > 0.0) != (find_b_trajectory_error(1.0) > 0.0)) {
              local_b = find_zero_bisection(find_b_trajectory_error, 0.01, 1.0);
            } else {
              local_b = b_glm;
            }
            local_a = (d2 - dt) / std::pow(std::max(1e-6, H - h2), local_b);
            local_dip = local_a * std::pow(H - hw, local_b) + dt;
            sip = -local_a * local_b * std::pow(H - hw, local_b - 1.0);
            local_c1 = sip / (local_c2 * std::pow(hw, local_c2 - 1.0));
            local_c0 = std::min(local_dip - local_c1 * std::pow(hw, local_c2), limit_c0);
        } else {
          local_b = b_glm;
            auto f_opt_c2 = [&](double c2_trial) noexcept -> double {
                double denom = std::max(1e-6, std::pow(h2, c2_trial) - std::pow(h1, c2_trial));
                double c1_t = (d2 - d1) / denom;
                double c0_t = d1 - c1_t * std::pow(h1, c2_trial);
                double dw_t = std::min(limit_c0, c0_t + c1_t * std::pow(hw, c2_trial));
                double sw_t = c1_t * c2_trial * std::pow(hw, c2_trial - 1.0);
                double a_t = (dw_t - dt) / std::pow(H - hw, local_b);
                double sw_expected = -a_t * local_b * std::pow(H - hw, local_b - 1.0);
                return (sw_t - sw_expected) * (sw_t - sw_expected);
            };
            
            // Löst das Optimierungsproblem genau wie Julia via Brent-Verfahren im Intervall [0.1, 1.0]
            local_c2 = optimize_brent(f_opt_c2, 0.1, 1.0); 
            
            double denom = std::max(1e-6, std::pow(h2, local_c2) - std::pow(h1, local_c2));
            local_c1 = (d2 - d1) / denom;
            local_c0 = d1 - local_c1 * std::pow(h1, local_c2);
            local_dip = local_c0 + local_c1 * std::pow(hw, local_c2);
            sip = local_c1 * local_c2 * std::pow(hw, local_c2 - 1.0);
            //sip = local_c1 * local_c2 * std::pow(hw, 0.0001 + local_c2 - 1.0);
            local_a = (local_dip - dt) / std::pow(H - hw, local_b);
        }
    } else {
        // N >= 3 PUNKTE:
        // Vorab-Berechnung der quadrierten Messdurchmesser für Performance
        std::vector<double> d_meas_sq(n_pts);
        size_t pts_below = 0;
        size_t pts_above = 0;
        for (size_t i = 0; i < n_pts; ++i) {
            d_meas_sq[i] = d_meas[i] * d_meas[i];
            if (h_meas[i] <= hw) pts_below++;
            else pts_above++;
        }

        // Variablen für Ergebnisse deklarieren (Kein double davor! Verhindert Shadowing)
        local_b = b_glm;
        local_c2 = c2_glm;
        local_a = 0.0; 
        local_dip = 0.0; 
        local_c1 = 0.0; 
        local_c0 = 0.0;

        // =================================================================
        // SONDERFALL A: Alle Messpunkte liegen im unteren Stammfußbereich
        // =================================================================
        if (pts_above == 0) {
            if (h_ref.has_value()) {
                // --- Variante mit Ankerpunkt ---
                size_t idx_anchor = 0;
                double min_dist = std::abs(h_meas[0] - h_ref.value());
                for (size_t i = 1; i < n_pts; ++i) {
                    double dist = std::abs(h_meas[i] - h_ref.value());
                    if (dist < min_dist) { min_dist = dist; idx_anchor = i; }
                }
                double h_anchor = h_meas[idx_anchor];
                double d_anchor = d_meas[idx_anchor];

                auto f_opt_all_below = [&](const std::vector<double>& x) noexcept -> double {
                    double c2_t = x[1];
                    double b_fixed = b_glm;
                    double k_slope = -(b_fixed * std::pow(std::max(1e-6, H - hw), b_fixed - 1.0)) / (c2_t * std::pow(hw, c2_t - 1.0));
                    double Nenner = std::pow(std::max(1e-6, H - hw), b_fixed) - k_slope * (std::pow(hw, c2_t) - std::pow(h_anchor, c2_t));
                    if (std::abs(Nenner) < 1e-6) Nenner = (Nenner >= 0.0 ? 1e-6 : -1e-6);

                    double a_t = (d_anchor - dt) / Nenner;
                    double c1_t = a_t * k_slope;
                    double c0_t = d_anchor - c1_t * std::pow(h_anchor, c2_t);

                    if (c0_t > limit_c0) {
                        c0_t = limit_c0;
                        c1_t = (d_anchor - limit_c0) / std::max(1e-6, std::pow(h_anchor, c2_t));
                    }

                    double total_error = 0.0;
                    for (size_t i = 0; i < n_pts; ++i) {
                        double d_mod = c0_t + c1_t * std::pow(h_meas[i], c2_t);
                        double err = (d_mod * d_mod) - d_meas_sq[i];
                        total_error += err * err;
                    }
                    return total_error;
                };

                std::vector<double> lb = {b_glm, 0.1};
                std::vector<double> ub = {b_glm, 1.0};
                std::vector<double> init = {b_glm, c2_glm};
                std::vector<double> best_params = optimize_nelder_mead_2d(f_opt_all_below, lb, ub, init);
                
                local_b = b_glm;
                local_c2 = best_params[1];

                double k_slope = -(local_b * std::pow(std::max(1e-6, H - hw), local_b - 1.0)) / (local_c2 * std::pow(hw, local_c2 - 1.0));
                double Nenner = std::pow(std::max(1e-6, H - hw), local_b) - k_slope * (std::pow(hw, local_c2) - std::pow(h_anchor, local_c2));
                if (std::abs(Nenner) < 1e-6) Nenner = (Nenner >= 0.0 ? 1e-6 : -1e-6);

                local_a = (d_anchor - dt) / Nenner;
                local_c1 = local_a * k_slope;
                local_c0 = d_anchor - local_c1 * std::pow(h_anchor, local_c2);

                if (local_c0 > limit_c0) {
                    local_c0 = limit_c0;
                    local_c1 = (d_anchor - limit_c0) / std::max(1e-6, std::pow(h_anchor, local_c2));
                    local_dip = local_c0 + local_c1 * std::pow(hw, local_c2);
                    local_a = (local_dip - dt) / std::pow(std::max(1e-6, H - hw), local_b);
                } else {
                    local_dip = local_a * std::pow(std::max(1e-6, H - hw), local_b) + dt;
                }
            } else {
                // --- Variante OHNE Ankerpunkt ---
                auto f_opt_all_below_free = [&](const std::vector<double>& x) noexcept -> double {
                    double c0_t = x[0];
                    double c2_t = x[1];
                    double b_fixed = b_glm;

                    double A = std::pow(hw, c2_t);
                    double B = (c2_t * std::pow(hw, c2_t - 1.0) * std::max(1e-6, H - hw)) / b_fixed;
                    double c1_t = (dt - c0_t) / std::max(1e-6, A + B);

                    double total_error = 0.0;
                    for (size_t i = 0; i < n_pts; ++i) {
                        double d_mod = c0_t + c1_t * std::pow(h_meas[i], c2_t);
                        double err = (d_mod * d_mod) - d_meas_sq[i];
                        total_error += err * err;
                    }
                    return total_error;
                };

                std::vector<double> lb = {1.0, 0.1};
                std::vector<double> ub = {limit_c0, 1.0};
                std::vector<double> init = {std::sqrt(d_meas_sq[0]), c2_glm};
                std::vector<double> best_params = optimize_nelder_mead_2d(f_opt_all_below_free, lb, ub, init);

                local_b = b_glm;
                local_c0 = best_params[0];
                local_c2 = best_params[1];

                double A_f = std::pow(hw, local_c2);
                double B_f = (local_c2 * std::pow(hw, local_c2 - 1.0) * std::max(1e-6, H - hw)) / local_b;
                local_c1 = (dt - local_c0) / std::max(1e-6, A_f + B_f);
                
                local_dip = local_c0 + local_c1 * std::pow(hw, local_c2);
                local_a = (local_dip - dt) / std::pow(std::max(1e-6, H - hw), local_b);
            }
        }
        // =================================================================
        // SONDERFALL B: Alle Messpunkte liegen im oberen Schaftbereich
        // =================================================================
        else if (pts_below == 0) {
            if (h_ref.has_value()) {
                // --- Variante mit Ankerpunkt ---
                size_t idx_anchor = 0;
                double min_dist = std::abs(h_meas[0] - h_ref.value());
                for (size_t i = 1; i < n_pts; ++i) {
                    double dist = std::abs(h_meas[i] - h_ref.value());
                    if (dist < min_dist) { min_dist = dist; idx_anchor = i; }
                }
                double h_anchor = h_meas[idx_anchor];
                double d_anchor = d_meas[idx_anchor];

                auto f_opt_all_above = [&](const std::vector<double>& x) noexcept -> double {
                    double b_t = x[0];
                    double c2_fixed = c2_glm;
                    double a_t = (d_anchor - dt) / std::pow(std::max(1e-6, H - h_anchor), b_t);
                    double dw_t = a_t * std::pow(std::max(1e-6, H - hw), b_t) + dt;
                    double sw_t = -a_t * b_t * std::pow(std::max(1e-6, H - hw), b_t - 1.0);
                    double c1_t = sw_t / (c2_fixed * std::pow(hw, c2_fixed - 1.0));
                    double c0_t = dw_t - c1_t * std::pow(hw, c2_fixed);

                    if (c0_t > limit_c0) {
                        c0_t = limit_c0;
                        c1_t = (dw_t - limit_c0) / std::pow(hw, c2_fixed);
                    }

                    double total_error = 0.0;
                    for (size_t i = 0; i < n_pts; ++i) {
                        double d_mod = a_t * std::pow(std::max(1e-6, H - h_meas[i]), b_t) + dt;
                        double err = (d_mod * d_mod) - d_meas_sq[i];
                        total_error += err * err;
                    }
                    return total_error;
                };

                std::vector<double> lb = {0.2, c2_glm};
                std::vector<double> ub = {1.0, c2_glm};
                std::vector<double> init = {b_glm, c2_glm};
                std::vector<double> best_params = optimize_nelder_mead_2d(f_opt_all_above, lb, ub, init);

                local_b = best_params[0];
                local_c2 = c2_glm;
                local_a = (d_anchor - dt) / std::pow(std::max(1e-6, H - h_anchor), local_b);
                
                local_dip = local_a * std::pow(std::max(1e-6, H - hw), local_b) + dt;
                double sw_final = -local_a * local_b * std::pow(std::max(1e-6, H - hw), local_b - 1.0);
                local_c1 = sw_final / (local_c2 * std::pow(hw, local_c2 - 1.0));
                local_c0 = local_dip - local_c1 * std::pow(hw, local_c2);

                if (local_c0 > limit_c0) {
                    local_c0 = limit_c0;
                    local_c1 = (local_dip - limit_c0) / std::pow(hw, local_c2);
                }
            } else {
                // --- Variante OHNE Ankerpunkt ---
                auto f_opt_all_above_free = [&](const std::vector<double>& x) noexcept -> double {
                    double b_t = x[0];
                    double a_t = x[1];

                    double total_error = 0.0;
                    for (size_t i = 0; i < n_pts; ++i) {
                        double d_mod = a_t * std::pow(std::max(1e-6, H - h_meas[i]), b_t) + dt;
                        double err = (d_mod * d_mod) - d_meas_sq[i];
                        total_error += err * err;
                    }
                    return total_error;
                };

                std::vector<double> lb = {0.2, 0.01};
                std::vector<double> ub = {1.0, 50.0};
                std::vector<double> init = {b_glm, 10.0};
                std::vector<double> best_params = optimize_nelder_mead_2d(f_opt_all_above_free, lb, ub, init);

                local_b = best_params[0];
                local_a = best_params[1];
                local_c2 = c2_glm;
local_dip = local_a * std::pow(std::max(1e-6, H - hw), local_b) + dt;
double sw_final = -local_a * local_b * std::pow(std::max(1e-6, H - hw), local_b - 1.0);
local_c1 = sw_final / (local_c2 * std::pow(hw, local_c2 - 1.0));
local_c0 = local_dip - local_c1 * std::pow(hw, local_c2);
if (local_c0 > limit_c0) {
local_c0 = limit_c0;
local_c1 = (local_dip - limit_c0) / std::pow(hw, local_c2);
}
}
}
       // =================================================================
        // NORMALFALL: Punkte sind sowohl über als auch unter hw verteilt
        // =================================================================
        else {
            if (h_ref.has_value()) {
                // --- Normalfall MIT Ankerpunkt ---
                size_t idx_anchor = 0;
                double min_dist = std::abs(h_meas[0] - h_ref.value());
                for (size_t i = 1; i < n_pts; ++i) {
                    double dist = std::abs(h_meas[i] - h_ref.value());
                    if (dist < min_dist) { min_dist = dist; idx_anchor = i; }
                }
                double h_anchor = h_meas[idx_anchor];
                double d_anchor = d_meas[idx_anchor];

                if (n_pts == 3) {
                    auto calc_parameters_3pts = [&](double b, double c2, double& a, double& c0, double& c1) {
                        double k_slope = -(b * std::pow(std::max(1e-6, H - hw), b - 1.0)) / (c2 * std::pow(hw, c2 - 1.0));
                        if (h_anchor <= hw) {
                            double Nenner = std::pow(std::max(1e-6, H - hw), b) - k_slope * (std::pow(hw, c2) - std::pow(h_anchor, c2));
                            if (std::abs(Nenner) < 1e-6) Nenner = (Nenner >= 0.0 ? 1e-6 : -1e-6);
                            a = (d_anchor - dt) / Nenner;
                            c1 = a * k_slope;
                            c0 = d_anchor - c1 * std::pow(h_anchor, c2);
                        } else {
                            a = (d_anchor - dt) / std::pow(std::max(1e-6, H - h_anchor), b);
                            double dw = a * std::pow(std::max(1e-6, H - hw), b) + dt;
                            double sw = -a * b * std::pow(std::max(1e-6, H - hw), b - 1.0);
                            c1 = sw / (c2 * std::pow(hw, c2 - 1.0));
                            c0 = dw - c1 * std::pow(hw, c2);
                        }
                    };

                    // HIER FIX: <double> explizit hinzugefügt
                    auto f_opt_3pts = [&](const std::vector<double>& x) noexcept -> double {
                        double b_t = x[0];   
                        double c2_t = x[1];  
                        double a_t, c0_t, c1_t;
                        calc_parameters_3pts(b_t, c2_t, a_t, c0_t, c1_t);
                        if (c0_t > limit_c0) {
                            c0_t = limit_c0;
                            if (h_anchor <= hw) {
                                c1_t = (d_anchor - limit_c0) / std::max(1e-6, std::pow(h_anchor, c2_t));
                                double dw_t = c0_t + c1_t * std::pow(hw, c2_t);
                                a_t = (dw_t - dt) / std::pow(std::max(1e-6, H - hw), b_t);
                            } else {
                                double dw_t = a_t * std::pow(std::max(1e-6, H - hw), b_t) + dt;
                                c1_t = (dw_t - limit_c0) / std::pow(hw, c2_t);
                            }
                        }

                        double total_error = 0.0;
                        for (size_t i = 0; i < 3; ++i) {
                            double d_mod = (h_meas[i] <= hw) ? (c0_t + c1_t * std::pow(h_meas[i], c2_t)) 
                                                             : (a_t * std::pow(std::max(1e-6, H - h_meas[i]), b_t) + dt);
                            double err = (d_mod * d_mod) - d_meas_sq[i];
                            total_error += err * err;
                        }
                        return total_error;
                    };

                    std::vector<double> lb = {0.2, 0.1};
                    std::vector<double> ub = {1.0, 1.0};
                    std::vector<double> init = {b_glm, c2_glm};
                    std::vector<double> best_params = optimize_nelder_mead_2d(f_opt_3pts, lb, ub, init);
                    local_b = best_params[0];
                    local_c2 = best_params[1]; 
                    calc_parameters_3pts(local_b, local_c2, local_a, local_c0, local_c1);

                    if (local_c0 > limit_c0) {
                        local_c0 = limit_c0;
                        if (h_anchor <= hw) {
                            local_c1 = (d_anchor - limit_c0) / std::max(1e-6, std::pow(h_anchor, local_c2));
                            local_dip = local_c0 + local_c1 * std::pow(hw, local_c2);
                            local_a = (local_dip - dt) / std::pow(std::max(1e-6, H - hw), local_b);
                        } else {
                            local_dip = local_a * std::pow(std::max(1e-6, H - hw), local_b) + dt;
                            local_c1 = (local_dip - limit_c0) / std::pow(hw, local_c2);
                        }
                    } else {
                        local_dip = local_a * std::pow(std::max(1e-6, H - hw), local_b) + dt;
                    }
                } else {
                    // Mehr als 3 Punkte mit Anker
                    auto calc_lower_parameters = [&](double b, double c2, double a, double& c0, double& c1) {
                        double dw = a * std::pow(std::max(1e-6, H - hw), b) + dt;
                        double sw = -a * b * std::pow(std::max(1e-6, H - hw), b - 1.0);
                        c1 = sw / (c2 * std::pow(hw, c2 - 1.0));
                        c0 = dw - c1 * std::pow(hw, c2);
                    };

                    // HIER FIX: <double> explizit hinzugefügt
                    auto f_opt_multi = [&](const std::vector<double>& x) noexcept -> double {
                        double b_t = x[0];
                        double c2_t = x[1];
                        double a_t = (d_anchor - dt) / std::pow(std::max(1e-6, H - h_anchor), b_t); 
                        double c0_t, c1_t;
                        calc_lower_parameters(b_t, c2_t, a_t, c0_t, c1_t);

                        if (c0_t > limit_c0) {
                            c0_t = limit_c0;
                            double dw_t = a_t * std::pow(std::max(1e-6, H - hw), b_t) + dt;
                            c1_t = (dw_t - limit_c0) / std::pow(hw, c2_t);
                        }

                        double total_error = 0.0;
                        for (size_t i = 0; i < n_pts; ++i) {
                            double d_mod = (h_meas[i] <= hw) ? (c0_t + c1_t * std::pow(h_meas[i], c2_t)) 
                                                             : (a_t * std::pow(std::max(1e-6, H - h_meas[i]), b_t) + dt);
                            double err = (d_mod * d_mod) - d_meas_sq[i];
                            total_error += err * err;
                        }
                        return total_error;
                    };

                    std::vector<double> lb = {0.2, 0.1};
                    std::vector<double> ub = {1.0, 1.0};
                    std::vector<double> init = {b_glm, c2_glm};
                    std::vector<double> best_params = optimize_nelder_mead_2d(f_opt_multi, lb, ub, init);
                    local_b = best_params[0];
                    local_c2 = best_params[1];
                    local_a = (d_anchor - dt) / std::pow(std::max(1e-6, H - h_anchor), local_b);
                    calc_lower_parameters(local_b, local_c2, local_a, local_c0, local_c1);

                    if (local_c0 > limit_c0) {
                        local_c0 = limit_c0;
                        local_dip = local_a * std::pow(std::max(1e-6, H - hw), local_b) + dt;
                        local_c1 = (local_dip - limit_c0) / std::pow(hw, local_c2);
                    } else {
                        local_dip = local_a * std::pow(std::max(1e-6, H - hw), local_b) + dt;
                    }
                }
            } else {
                // --- Normalfall OHNE Ankerpunkt (Freier Fit) ---
                auto calc_lower_parameters_free = [&](double b, double c2, double a, double& c0, double& c1) {
                    double dw = a * std::pow(std::max(1e-6, H - hw), b) + dt;
                    double sw = -a * b * std::pow(std::max(1e-6, H - hw), b - 1.0);
                    c1 = sw / (c2 * std::pow(hw, c2 - 1.0));
                    c0 = dw - c1 * std::pow(hw, c2);
                };

                // HIER FIX: <double> explizit hinzugefügt
                auto f_opt_multi_free = [&](const std::vector<double>& x) noexcept -> double {
                    double b_t = x[0]; 
                    double a_t = x[1]; 
                    double c2_fixed = c2_glm;
                    double c0_t, c1_t;
                    calc_lower_parameters_free(b_t, c2_fixed, a_t, c0_t, c1_t);

                    if (c0_t > limit_c0) {
                        c0_t = limit_c0;
                        double dw_t = a_t * std::pow(std::max(1e-6, H - hw), b_t) + dt;
                        c1_t = (dw_t - limit_c0) / std::pow(hw, c2_fixed);
                    }

                    double total_error = 0.0;
                    for (size_t i = 0; i < n_pts; ++i) {
                        double d_mod = (h_meas[i] <= hw) ? (c0_t + c1_t * std::pow(h_meas[i], c2_fixed)) 
                                                         : (a_t * std::pow(std::max(1e-6, H - h_meas[i]), b_t) + dt);
                        double err = (d_mod * d_mod) - d_meas_sq[i];
                        total_error += err * err;
                    }
                    return total_error;
                };

                std::vector<double> lb = {0.2, 0.01};
                std::vector<double> ub = {1.0, 50.0};
                std::vector<double> init = {b_glm, 1.0}; 
                
                std::vector<double> best_params = optimize_nelder_mead_2d(f_opt_multi_free, lb, ub, init);
                local_b = best_params[0]; 
                local_a = best_params[1]; 
                local_c2 = c2_glm;
                
                calc_lower_parameters_free(local_b, local_c2, local_a, local_c0, local_c1);

                if (local_c0 > limit_c0) {
local_c0 = limit_c0;
local_dip = local_a * std::pow(std::max(1e-6, H - hw), local_b) + dt;
local_c1 = (local_dip - limit_c0) / std::pow(hw, local_c2);
} else {
local_dip = local_a * std::pow(std::max(1e-6, H - hw), local_b) + dt;
}
}
}
        
}
out_is_valid = !(std::isnan(local_c0) || std::isinf(local_c0) || local_c0 > limit_c0);
return TaperProfile{H, hw, local_c0, local_dip, dt, local_a, local_b, local_c1, local_c2};
}
      
// =========================================================================
// 3. UNIVERSAL HAUPTFUNKTION (REGELT DIE BHD-SUCHE)
// =========================================================================

// 4. HAUPTFUNKTION (Löst den Verzweigungs- und Argument-Fehler)
TaperProfile solve_taper(DoubleVector h_meas, DoubleVector d_meas, double H, double max_c0_factor = 2.0, std::optional<double> h_ref = std::nullopt) {
    size_t n_pts = h_meas.size();
    if (n_pts != d_meas.size()) {
        throw std::invalid_argument("Vektoren muessen gleich lang sein.");
    }
    sort_by_height(h_meas, d_meas); // Absolute strukturelle Segment-Sicherheit

    double dt = (H >= 1.3) ? 0.8 : 0.1 + 0.7 * (1.0 - std::pow((1.3 - H) / 1.3, 2.0));
    double hw = 1.8104 * std::pow(1.0 + H, 0.358966) - 1.8104;
    bool unused_valid = true;

    // 1. SPEZIALFALL: KLEINSTBÄUME ABFANGEN
    if (H < 1.3) {
        double dr_small = dt + H;
        return _solve_taper_with_dr(h_meas, d_meas, H, dt, hw, dr_small, max_c0_factor, unused_valid, h_ref);
    }

    // 2. PRÜFEN, OB EIN ECHTER BHD VORLIEGT (h == 1.3)
    ptrdiff_t dbh_idx = -1;
    for (size_t i = 0; i < n_pts; ++i) {
        if (std::abs(h_meas[i] - 1.3) < 1e-4) {
            dbh_idx = static_cast<ptrdiff_t>(i);
            break;
        }
    }

    if (dbh_idx != -1) {
        // Fall A: BHD direkt bekannt (Keine Iteration nötig)
        double dr = d_meas[static_cast<size_t>(dbh_idx)] + 1.3;
        return _solve_taper_with_dr(h_meas, d_meas, H, dt, hw, dr, max_c0_factor, unused_valid, h_ref);
    } else {
        // Fall B: Kein BHD gemessen -> SYNCHRONISATION AUF find_zero_secant (Order0)
        // Dynamisch optimierter, biologischer Startwert für die Secant-Suche
        double start_dbh = H;
        if (n_pts > 0) {
            // Falls Messungen vorliegen, nutzen wir den dicksten gemessenen Durchmesser
            // und schätzen ihn leicht nach oben, falls die Messung höher am Stamm lag
            double max_measured_d = 0.0;
            double lowest_measured_h = H;
            for (size_t i = 0; i < n_pts; ++i) {
                if (d_meas[i] > max_measured_d) {
                    max_measured_d = d_meas[i];
                    lowest_measured_h = h_meas[i];
                }
            }
            start_dbh = (lowest_measured_h > 1.3) ? max_measured_d * (H - 1.3) / (H - lowest_measured_h) : max_measured_d;
        }

        auto f_outer_root_secant = [&](double dbh_trial) noexcept -> double {
            double dr_trial = dbh_trial + 1.3;
            bool inner_valid = true;
            TaperProfile prof = _solve_taper_with_dr(h_meas, d_meas, H, dt, hw, dr_trial, max_c0_factor, inner_valid, h_ref);
            
            // Falls die innere Modellierung mathematisch unlösbar wird, signalisieren wir dies
            if (!inner_valid) return std::numeric_limits<double>::max();
            
            double d_model_13 = (1.3 <= prof.hw) ? (prof.c0 + prof.c1 * std::pow(1.3, prof.c2)) 
                                                 : (prof.a * std::pow(H - 1.3, prof.b) + prof.dt);
            return d_model_13 - dbh_trial;
        };

        // Aufruf des implementierten Sekantenverfahrens (äquivalent zu Julias Order0)
        double final_dbh = find_zero_secant(f_outer_root_secant, start_dbh);
        
        // Letzte Absicherung falls die Secant-Suche wegen extremer max_c0_factor-Werte ins Leere lief
        if (std::isnan(final_dbh) || std::isinf(final_dbh) || final_dbh <= 0.0) {
            final_dbh = (n_pts > 0) ? d_meas[0] : H; // Robuster biologischer Fallback
        }

        double dr_final = final_dbh + 1.3;
        return _solve_taper_with_dr(h_meas, d_meas, H, dt, hw, dr_final, max_c0_factor, unused_valid, h_ref);
    }
}
TaperProfile solve_taperOLD(DoubleVector h_meas, DoubleVector d_meas, double H, double max_c0_factor = 2.0) {
  size_t n_pts = h_meas.size();
  if (n_pts != d_meas.size()) {
    throw std::invalid_argument("Vektoren muessen gleich lang sein.");
  }
  sort_by_height(h_meas, d_meas);
  double dt = (H >= 1.3) ? 0.8 : 0.1 + 0.7 * (1.0 - std::pow((1.3 - H) / 1.3, 2.0));
  double hw = 1.8104 * std::pow(1.0 + H, 0.358966) - 1.8104;
  bool unused_valid = true;
  if (H < 1.3) {
    double dr_small = dt + H;
    return _solve_taper_with_dr(h_meas, d_meas, H, dt, hw, dr_small, max_c0_factor, unused_valid);
  }
  ptrdiff_t dbh_idx = -1;
  for (size_t i = 0; i < n_pts; ++i) {
    if (std::abs(h_meas[i] - 1.3) < 1e-4) {
      dbh_idx = static_cast<ptrdiff_t>(i);
      break;
    }
  }
  if (dbh_idx != -1) {
    double dr = d_meas[static_cast<size_t>(dbh_idx)] + 1.3;
    return _solve_taper_with_dr(h_meas, d_meas, H, dt, hw, dr, max_c0_factor, unused_valid);
  } else {
    if (n_pts == 0) {
      double lower_dbh = 0.1;
            double upper_dbh = H * 3.0; // Biologische Maximalschranke ohne Messpunkte
            
            auto f_outer_root_n0 = [&](double dbh_trial) -> double {
                double dr_trial = dbh_trial + 1.3;
                bool inner_valid = true;
                TaperProfile prof = _solve_taper_with_dr(h_meas, d_meas, H, dt, hw, dr_trial, max_c0_factor, inner_valid);
                double d_model_13 = (1.3 <= prof.hw) ? (prof.c0 + prof.c1 * std::pow(1.3, prof.c2)) : (prof.a * std::pow(H - 1.3, prof.b) + prof.dt);
                return d_model_13 - dbh_trial;
            };
            
            double final_dbh = find_zero_bisection(f_outer_root_n0, lower_dbh, upper_dbh);
            return _solve_taper_with_dr(h_meas, d_meas, H, dt, hw, final_dbh + 1.3, max_c0_factor, unused_valid);
    }
    double lower_bound_from_above = 0.0;
    double upper_bound_from_below = std::numeric_limits<double>::max();
    for (size_t i = 0; i < n_pts; ++i) {
      if (h_meas[i] > 1.3) lower_bound_from_above = std::max(lower_bound_from_above, d_meas[i]);
      if (h_meas[i] <= 1.3) upper_bound_from_below = std::min(upper_bound_from_below, d_meas[i]);
    }
    double min_dbh = std::max(0.1, lower_bound_from_above);
    double max_dbh = (upper_bound_from_below == std::numeric_limits<double>::max()) ? (min_dbh * 3.0 + 10.0) : upper_bound_from_below;
    auto f_outer_root_with_pts = [&](double dbh_trial) -> double {
      double dr_trial = dbh_trial + 1.3;
      bool inner_valid = true;
      TaperProfile prof = _solve_taper_with_dr(h_meas, d_meas, H, dt, hw, dr_trial, max_c0_factor, inner_valid);
      if (!inner_valid) return std::numeric_limits<double>::max();
      double d_model_13 = (1.3 <= prof.hw) ? (prof.c0 + prof.c1 * std::pow(1.3, prof.c2)) : (prof.a * std::pow(H - 1.3, prof.b) + prof.dt);
      return d_model_13 - dbh_trial;
    };
    double final_dbh = 0.0;
    double f_min = f_outer_root_with_pts(min_dbh);
    double f_max = f_outer_root_with_pts(max_dbh);
    if ((f_min > 0.0) == (f_max > 0.0) || f_min == std::numeric_limits<double>::max() || f_max == std::numeric_limits<double>::max()) {
      final_dbh = (min_dbh + max_dbh) / 2.0;
    } else {
      final_dbh = find_zero_bisection(f_outer_root_with_pts, min_dbh, max_dbh);
    }
    double dr_final = final_dbh + 1.3;
    return _solve_taper_with_dr(h_meas, d_meas, H, dt, hw, dr_final, max_c0_factor, unused_valid);
}
}
                


// =========================================================================
// 3. APPLICATIONS & ANALYTICAL INTEGRATION
// =========================================================================

inline double get_diameter_at_height(const TaperProfile& p, double h) noexcept {
    if (h < 0.0 || h > p.H) return 0.0;
    return (h <= p.hw) ? (p.c0 + p.c1 * std::pow(h, p.c2)) 
                       : (p.a * std::pow(p.H - h, p.b) + p.dt);
}

inline double get_diameter_under_bark(const TaperProfile& p, double h) noexcept {
    double d_ob = get_diameter_at_height(p, h);
    if (d_ob <= 0.0) return 0.0;
    double dbt = (0.85149 + 0.60934 * d_ob - 0.00228 * d_ob * d_ob) / 10.0;
    return std::max(0.0, d_ob - dbt);
}

inline double get_height_at_diameterOLD(const TaperProfile& p, double d_target) noexcept {
    if (d_target >= p.c0) return 0.0;
    if (d_target <= p.dt) return p.H;
    if (d_target >= p.dw) {
        return std::max(0.0, std::pow((d_target - p.c0) / p.c1, 1.0 / p.c2));
    } else {
        return std::min(p.H, p.H - std::pow((d_target - p.dt) / p.a, 1.0 / p.b));
    }
}

double get_height_at_diameter(const TaperProfile& prof, double d_target) {
    if (d_target <= prof.dt) return prof.H; // Spitze erreicht
    
    // 1. Prüfen, ob der Grenzdurchmesser im oberen Segment liegt
    double d_at_hw = prof.a * std::pow(prof.H - prof.hw, prof.b) + prof.dt;
    
    if (d_target <= d_at_hw) {
        // Analytische Umkehrfunktion des oberen Segments (Wipfel)
        double ratio = (d_target - prof.dt) / prof.a;
        double h_target = prof.H - std::pow(ratio, 1.0 / prof.b);
        return h_target;
    } else {
        // Numerische Umkehrfunktion für das untere Segment (Stammfuß)
        auto f_height = [&](double h_trial) {
            return (prof.c0 + prof.c1 * std::pow(h_trial, prof.c2)) - d_target;
        };
        // Nutzen Sie Ihre implementierte Bisektion für die exakte Höhe unten
        return find_zero_bisection(f_height, 0.0, prof.hw);
    }
}

// Analytische Integration des Stammvolumens (Rinde / Holz)
double get_segment_volume(const TaperProfile& p, double h1, double h2, std::string_view type = "total") {
    double hs = std::max(0.0, std::min(p.H, h1));
    double he = std::max(0.0, std::min(p.H, h2));
    if (hs >= he) return 0.0;
    
    // Falls das Segment den Inflexionspunkt (hw) schneidet, splitten wir die Integration auf
    if (hs < p.hw && he > p.hw) {
        return get_segment_volume(p, hs, p.hw, type) + get_segment_volume(p, p.hw, he, type);
    }
    
    const bool is_basal = (he <= p.hw);
    constexpr double factor = std::numbers::pi / 40000.0; // Umrechnung von cm² in m²

    if (type == "total") {
        if (is_basal) {
            // Unteres Segment: Integral über (c0 + c1 * h^c2)^2
            auto F_low = [&](double h) noexcept -> double {
                return p.c0 * p.c0 * h 
                     + (2.0 * p.c0 * p.c1 * std::pow(h, p.c2 + 1.0)) / (p.c2 + 1.0) 
                     + (p.c1 * p.c1 * std::pow(h, 2.0 * p.c2 + 1.0)) / (2.0 * p.c2 + 1.0);
            };
            return factor * (F_low(he) - F_low(hs));
    /*
    if (type == "total") {
    if (is_basal) {
        // Exakt an Julias fehlerhaften Nenner (p.c2 + 1.0) beim letzten Term angepasst
        auto F_low = [&](double h) noexcept -> double {
            return p.c0 * p.c0 * h 
                 + (2.0 * p.c0 * p.c1 * std::pow(h, p.c2 + 1.0)) / (p.c2 + 1.0) 
                 + (p.c1 * p.c1 * std::pow(h, 2.0 * p.c2 + 1.0)) / (p.c2 + 1.0); // <-- HIER zurückgesetzt!
        };
        return factor * (F_low(he) - F_low(hs));
    */
        } else {
            // Oberes Segment: Integral über (a * (H - h)^b + dt)^2 (Minuszeichen durch innere Ableitung!)
            auto F_up = [&](double h) noexcept -> double {
                return p.dt * p.dt * h 
                     - (2.0 * p.a * p.dt * std::pow(p.H - h, p.b + 1.0)) / (p.b + 1.0) 
                     - (p.a * p.a * std::pow(p.H - h, 2.0 * p.b + 1.0)) / (2.0 * p.b + 1.0);
            };
            return factor * (F_up(he) - F_up(hs));
        }
    } else if (type == "wood") {
        const BarkCoefficients A = get_bark_coefficients();
        const double bark_array[5] = {A.A0, A.A1, A.A2, A.A3, A.A4};
        
        if (is_basal) {
            auto compute_sum = [&](double h) noexcept -> double {
                double total_sum = 0.0;
                for (int n = 0; n <= 4; ++n) {
                    for (int k = 0; k <= n; ++k) {
                        double bc = binomial_coeff(n, k);
                        total_sum += bark_array[n] * bc 
                                   * std::pow(p.c0, n - k) 
                                   * std::pow(p.c1, k) 
                                   * std::pow(h, static_cast<double>(k) * p.c2 + 1.0) 
                                   / (static_cast<double>(k) * p.c2 + 1.0);
                    }
                }
                return total_sum;
            };
            return factor * (compute_sum(he) - compute_sum(hs));
        } else {
            // Oberes Segment für Holz unter Rinde (Minuszeichen bei Substitution beachten)
            auto compute_sum_up = [&](double h) noexcept -> double {
                double total_sum = 0.0;
                for (int n = 0; n <= 4; ++n) {
                    for (int k = 0; k <= n; ++k) {
                        double bc = binomial_coeff(n, k);
                        double term = std::pow(p.dt, n - k) 
                                    * std::pow(p.a, k) 
                                    * std::pow(p.H - h, static_cast<double>(k) * p.b + 1.0) 
                                    / (static_cast<double>(k) * p.b + 1.0);
                        total_sum -= bark_array[n] * bc * term; // Minus durch -dh Substitution
                    }
                }
                return total_sum;
            };
            return factor * (compute_sum_up(he) - compute_sum_up(hs));
        }
    }
    return 0.0;
}

// =========================================================================
// 4. SORTIMENTIERUNG & AUSGABE
// =========================================================================

/**
 * Berechnet Holzsortimente (Stock, Derbholz, Wipfel) und gibt eine 
 * strukturierte forstliche Bilanzierung in der Konsole aus.
 */
TaperProfile calculate_assortments(const std::vector<double>& h_meas, const std::vector<double>& d_meas, 
                                   double H, double h_stump = 0.3, double d_top = 7.0) {
    size_t n_pts = h_meas.size();
    
    // Modell kalibrieren
    TaperProfile p = solve_taper(h_meas, d_meas, H);
    
    // Derbholzgrenze bestimmen
    double h_derb = get_height_at_diameter(p, d_top);
    if (h_derb < h_stump) {
        h_derb = h_stump;
    }
    
    // Brutto-Volumina (m³) mit Rinde
    double v_stump_tot = get_segment_volume(p, 0.0, h_stump, "total");
    double v_merc_tot  = get_segment_volume(p, h_stump, h_derb, "total");
    double v_tip_tot   = get_segment_volume(p, h_derb, H, "total");
    double v_total_tot = v_stump_tot + v_merc_tot + v_tip_tot;
    
    // Netto-Volumina (m³) ohne Rinde (Holz)
    double v_stump_wood = get_segment_volume(p, 0.0, h_stump, "wood");
    double v_merc_wood  = get_segment_volume(p, h_stump, h_derb, "wood");
    double v_tip_wood   = get_segment_volume(p, h_derb, H, "wood");
    double v_total_wood = v_stump_wood + v_merc_wood + v_tip_wood;
    
    // Text für gemessene Durchmesser aufbereiten (Sichere Alternative zu Julias Vektor-Print)
    std::string dbh_print;
    if (n_pts > 0) {
        dbh_print = std::format("{} Messpunkte vorhanden (Basis-D: {} cm)", n_pts, d_meas[0]);
    } else {
        dbh_print = "Kein Durchmesser angegeben (H/D=100 Fallback)";
    }

    // Ausgabe der Kopfdaten
    std::println("\n=======================================================");
    std::println("ERGEBNISSE C++26-EINZELBAUM (H = {:.2f}m, {})", H, dbh_print);
    std::println("Derbholzgrenze ({:.1f}cm) erreicht bei h = {:.2f}m", d_top, h_derb);
    std::println("=======================================================");
    
    // Tabellenkopf formatiert ausgeben
    std::println("{:<18} | {:<7} | {:<7} | {:<12} | {:<11} | {:<11}", 
                 "Sektion", "Start_m", "Ende_m", "Vol_Gross_m3", "Vol_Wood_m3", "Vol_Bark_m3");
    std::println("--------------------------------------------------------------------------------");
    
    // Lambda-Hilfsfunktion für das Zeilen-Formatting (C++20/C++23 kompatibel innerhalb der Funktion)
    auto print_row = [](std::string_view name, double start, double end, double gross, double wood) {
        std::println("{:<18} | {:<7.2f} | {:<7.2f} | {:<12.4f} | {:<11.4f} | {:<11.4f}", 
                     name, start, end, gross, wood, gross - wood);
    };

    // Zeilen ausgeben
    print_row("Stock (Stump)", 0.0, h_stump, v_stump_tot, v_stump_wood);
    print_row("Derbholz (Merc)", h_stump, h_derb, v_merc_tot, v_merc_wood);
    print_row("Wipfel (Tip)", h_derb, H, v_tip_tot, v_tip_wood);
    print_row("Gesamt (Total)", 0.0, H, v_total_tot, v_total_wood);
    
    return p;
}

// =========================================================================
// 5. TEST UMGEBUNG (MAIN)
// =========================================================================
int main() {
    try {
        // =========================================================================
        // 1. Pfad: N = 0 (Keine Punkte)
        // =========================================================================
        std::println("\n=======================================================");
        std::println("[TEST 1] Pfad: N = 0 (Rein höhenbasierte Schätzung)");
        std::println("=======================================================");
        std::vector<double> h0, d0;
        calculate_assortments(h0, d0, 15.0, 0.3, 7.0);

        // =========================================================================
        // 2. Pfad: N = 1 (Punkt UNTER hw)
        // =========================================================================
        std::println("\n=======================================================");
        std::println("[TEST 2] Pfad: N = 1 (Punkt UNTER hw; z.B. BHD bei 1.3m)");
        std::println("=======================================================");
        std::vector<double> h1_under = {1.3};
        std::vector<double> d1_under = {28.0};
        calculate_assortments(h1_under, d1_under, 20.0, 0.3, 7.0);

        // =========================================================================
        // 3. Pfad: N = 1 (Punkt OBER hw)
        // =========================================================================
        std::println("\n=======================================================");
        std::println("[TEST 3] Pfad: N = 1 (Punkt OBER hw; z.B. Kronenansatz)");
        std::println("=======================================================");
        std::vector<double> h1_over = {12.0}; // hw liegt bei H=18m ca. bei ~3.2m
        std::vector<double> d1_over = {14.0};
        calculate_assortments(h1_over, d1_over, 18.0, 0.3, 7.0);

        // =========================================================================
        // 4. Pfad: N = 2 (Beide Punkte OBER hw)
        // =========================================================================
        std::println("\n=======================================================");
        std::println("[TEST 4] Pfad: N = 2 (Beide Punkte OBER hw)");
        std::println("=======================================================");
        std::vector<double> h2_over = {6.0, 14.0}; // hw liegt bei H=22m ca. bei ~3.9m
        std::vector<double> d2_over = {22.0, 15.0};
        calculate_assortments(h2_over, d2_over, 22.0, 0.3, 7.0);

        // =========================================================================
        // 5. Pfad: N = 2 (Gemischt: h1 UNTER hw, h2 OBER hw)
        // =========================================================================
        std::println("\n=======================================================");
        std::println("[TEST 5] Pfad: N = 2 (Gemischt: h1 UNTER hw, h2 OBER hw)");
        std::println("=======================================================");
        std::vector<double> h2_mix = {1.3, 10.0}; // hw liegt bei H=25m ca. bei ~4.2m
        std::vector<double> d2_mix = {35.0, 27.0};
        calculate_assortments(h2_mix, d2_mix, 25.0, 0.3, 7.0);

        // =========================================================================
        // 6. Pfad: N = 2 (Beide Punkte UNTER hw) -> Löst Brent-Optimierung aus
        // =========================================================================
        std::println("\n=======================================================");
        std::println("[TEST 6] Pfad: N = 2 (Beide Punkte UNTER hw) -> Erfordert Optimierung");
        std::println("=======================================================");
        std::vector<double> h2_under = {0.5, 2.0}; // hw liegt bei H=30m ca. bei ~4.7m
        std::vector<double> d2_under = {48.0, 42.0};
        calculate_assortments(h2_under, d2_under, 30.0, 0.3, 7.0);

        // =========================================================================
        // 7. Pfad: N >= 3 (Standardfall mit stetigem Abfall)
        // =========================================================================
        std::println("\n=======================================================");
        std::println("[TEST 7] Pfad: N >= 3 (Sektionierung mit plausiblen Daten)");
        std::println("=======================================================");
        std::vector<double> h3_std = {1.3, 7.0, 15.0};
        std::vector<double> d3_std = {28.0, 22.0, 12.6};
        calculate_assortments(h3_std, d3_std, 20.0, 0.3, 7.0);

        // =========================================================================
        // 8. Pfad: N >= 3 (Gleichbleibender oder steigender Durchmesser) 
        // Löst den Fallback-Zweig 'local_b = b_glm' aus.
        // =========================================================================
        std::println("\n=======================================================");
        std::println("[TEST 8] Pfad: N >= 3 (Anomalie: d3 >= d2 -> Fallback auf GLM)");
        std::println("=======================================================");
        std::vector<double> h3_anomaly = {1.3, 8.0, 12.0};
        std::vector<double> d3_anomaly = {25.0, 16.0, 16.5}; // d3 ist größer als d2
        calculate_assortments(h3_anomaly, d3_anomaly, 20.0, 0.3, 7.0);

        {
          std::println("[TEST 2] Pfad: SONDERFALL A (Alle Punkte UNTER hw)"); 
          std::vector<double> h = {0.3, 0.8, 1.5};
          std::vector<double> d = {35.0, 31.0, 28.0};
          calculate_assortments(h, d, 20.0, 0.3, 7.0);
        }
        {
          std::println("[TEST 3] Pfad: SONDERFALL B (Alle Punkte OBERHALB hw)"); 
          std::vector<double> h = {6.0, 10.0, 16.0};
          std::vector<double> d = {21.0, 18.0, 11.0};
          calculate_assortments(h, d, 20.0, 0.3, 7.0);
        }
        {
          std::println("[TEST 4] Pfad: LIMIT_C0 AKTIV (Erzwungene Stammfuß-Kappung)"); 
          std::vector<double> h = {0.5, 5.0, 12.0};
          std::vector<double> d = {45.0, 20.0, 14.0};
          auto TP = calculate_assortments(h, d, 20.0, 0.3, 7.0);
          for (auto x : h) {std::cout << " " << get_diameter_at_height(TP, x);}
        }
        {
          std::println("[TEST 4] Pfad: LIMIT_C0 AKTIV (Erzwungene Stammfuß-Kappung)"); 
          std::vector<double> h = {1.3, 7.0, 15.0};
          std::vector<double> d = {28.0, 22.0, 12.6};
          auto TP = calculate_assortments(h, d, 20.0, 0.3, 7.0);
          for (auto x : h) {std::cout << " " << get_diameter_at_height(TP, x);}
        }
        {
          std::vector<double> h = {1.3, 7.0, 15.0};
          std::vector<double> d = {28.0, 22.0, 12.6};
          auto TP = solve_taper(h, d, 20.);
          std::cout << TP.H
                    << ", " << TP.hw
                    << ", " << TP.c0
                    << ", " << TP.dw
                    << ", " << TP.dt
                    << ", " << TP.a
                    << ", " << TP.b
                    << ", " << TP.c1
                    << ", " << TP.c2 << "\n";
          for (auto x : h) {std::cout << " " << get_diameter_at_height(TP, x);}
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Fehler bei den Testläufen: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}


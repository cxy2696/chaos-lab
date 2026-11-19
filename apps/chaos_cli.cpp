// chaos_cli: a small, dependency-free front end for experiments and data export.
// SPDX-License-Identifier: MIT
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "chaos/control.hpp"
#include "chaos/embedding.hpp"
#include "chaos/integrators.hpp"
#include "chaos/lyapunov.hpp"
#include "chaos/sindy.hpp"
#include "chaos/systems.hpp"
#include "chaos/timeseries.hpp"
#include "chaos/tools.hpp"
#include "chaos/version.hpp"

using namespace chaos;

namespace {

struct Options {
    std::map<std::string, std::string> values;
    std::vector<std::string> positional;

    bool has(const std::string& key) const { return values.count(key) != 0; }
    std::string get(const std::string& key, const std::string& fallback = "") const {
        const auto it = values.find(key);
        return it == values.end() ? fallback : it->second;
    }
    double get_double(const std::string& key, double fallback) const {
        const auto it = values.find(key);
        if (it == values.end()) return fallback;
        return std::stod(it->second);
    }
    int get_int(const std::string& key, int fallback) const {
        const auto it = values.find(key);
        if (it == values.end()) return fallback;
        return std::stoi(it->second);
    }
    std::size_t get_size(const std::string& key, std::size_t fallback) const {
        const auto it = values.find(key);
        if (it == values.end()) return fallback;
        return static_cast<std::size_t>(std::stoull(it->second));
    }
};

Options parse_options(int argc, char** argv) {
    Options opt;
    for (int i = 2; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg.rfind("--", 0) == 0) {
            const std::string key = arg.substr(2);
            if (i + 1 < argc && std::string(argv[i + 1]).rfind("--", 0) != 0) {
                opt.values[key] = argv[++i];
            } else {
                opt.values[key] = "1";
            }
        } else {
            opt.positional.push_back(arg);
        }
    }
    return opt;
}

void usage() {
    std::cout <<
        "chaos_cli " CHAOS_VERSION_STRING " - C/C++ chaos toolkit\n\n"
        "Usage: chaos_cli <command> [--key value ...]\n\n"
        "Commands:\n"
        "  models                          list ODEs and maps\n"
        "  simulate                        integrate/iterate and export CSV\n"
        "  lyapunov                        Benettin spectrum (ODE) or map exponents\n"
        "  lyapunov-scan                   lambda(r) scan for the logistic map\n"
        "  bifurcation                     Bifurcation diagram of the logistic map\n"
        "  poincare                        Poincare section of a flow\n"
        "  ami | fnn | dimension | rosenstein | zeroone | rqa | spectrum\n"
        "                                  scalar time-series diagnostics\n"
        "  sindy                           sparse regression from trajectory data\n"
        "  sync                            Lorenz drive-response synchronization\n"
        "  version                         print version\n\n"
        "Common options:\n"
        "  --system NAME  --map NAME  --params p1,p2,...  --x0 a,b,...\n"
        "  --h DT         --n N       --transient T       --out FILE\n"
        "  --input FILE   --column C  --stride S          --header\n";
}

std::vector<double> parse_list(const std::string& text) {
    std::vector<double> out;
    if (text.empty()) return out;
    std::stringstream ss(text);
    std::string item;
    while (std::getline(ss, item, ',')) {
        if (item.empty()) continue;
        out.push_back(std::stod(item));
    }
    return out;
}

std::vector<double> read_series(const std::string& path, int column, bool /*header*/) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    std::vector<double> out;
    std::string line;
    bool first = true;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string cell;
        int col = 0;
        std::string selected;
        while (std::getline(ss, cell, ',')) {
            if (col == column) selected = cell;
            ++col;
        }
        if (selected.empty()) {
            if (first) {
                first = false;
                continue;
            }
            throw std::runtime_error("column " + std::to_string(column) + " missing in " + path);
        }
        try {
            const double v = std::stod(selected);
            out.push_back(v);
            first = false;
        } catch (const std::exception&) {
            if (first) {
                first = false;  // header row
                continue;
            }
            throw std::runtime_error("non-numeric value in " + path + ": " + selected);
        }
    }
    if (out.size() < 10) throw std::runtime_error("series too short: " + path);
    return out;
}

State choose_state(const Options& opt) {
    const std::string text = opt.get("x0");
    if (!text.empty()) return parse_state(text);
    const std::string name = opt.has("system") ? opt.get("system") : opt.get("map");
    return default_initial_state(name);
}

OdeSystem choose_ode(const Options& opt) {
    const std::string name = opt.get("system");
    if (name.empty()) throw std::runtime_error("--system is required");
    return make_ode_system(name, parse_list(opt.get("params")));
}
MapSystem choose_map(const Options& opt) {
    const std::string name = opt.get("map");
    if (name.empty()) throw std::runtime_error("--map is required");
    return make_map_system(name, parse_list(opt.get("params")));
}

void cmd_models() {
    std::cout << "ODE systems:\n";
    for (const auto& s : ode_catalog()) {
        std::cout << "  " << std::left << std::setw(16) << s.name << " " << s.summary << "\n";
        std::cout << "      params:";
        for (std::size_t i = 0; i < s.param_names.size(); ++i) {
            std::cout << " " << s.param_names[i] << "=" << s.defaults[i];
        }
        std::cout << "\n";
    }
    std::cout << "\nMaps:\n";
    for (const auto& s : map_catalog()) {
        std::cout << "  " << std::left << std::setw(16) << s.name << " " << s.summary << "\n";
        std::cout << "      params:";
        for (std::size_t i = 0; i < s.param_names.size(); ++i) {
            std::cout << " " << s.param_names[i] << "=" << s.defaults[i];
        }
        std::cout << "\n";
    }
}

void cmd_simulate(const Options& opt) {
    const std::string out_path = opt.get("out");
    const std::size_t stride = std::max<std::size_t>(1, opt.get_size("stride", 1));
    const double h = opt.get_double("h", 0.01);
    const double transient = opt.get_double("transient", 0.0);

    if (opt.has("map") || (!opt.has("system") && opt.has("map"))) {
        MapSystem map = choose_map(opt);
        State x = opt.has("x0") ? parse_state(opt.get("x0")) : default_initial_state(map.name);
        const std::size_t n = opt.get_size("n", 1000);
        const std::size_t skip = static_cast<std::size_t>(std::llround(transient));
        for (std::size_t i = 0; i < skip; ++i) x = map.step(x);
        std::ofstream out;
        if (!out_path.empty()) {
            out.open(out_path);
            if (!out) throw std::runtime_error("cannot write " + out_path);
            out << "n";
            for (int j = 0; j < map.dim; ++j) out << ",x" << j;
            out << "\n";
            out << std::setprecision(12);
        }
        for (std::size_t i = 0; i < n; ++i) {
            x = map.step(x);
            if (i % stride != 0) continue;
            if (out) {
                out << i;
                for (double v : x) out << "," << v;
                out << "\n";
            }
        }
        if (!out) {
            std::cout << "map " << map.name << " done; last state:";
            for (double v : x) std::cout << " " << v;
            std::cout << "\n";
        }
        return;
    }

    OdeSystem sys = choose_ode(opt);
    State x = choose_state(opt);
    const std::size_t n = opt.get_size("n", 5000);
    const std::size_t skip = static_cast<std::size_t>(std::llround(transient / h));
    double t = 0.0;
    for (std::size_t i = 0; i < skip; ++i) {
        x = rk4_step(sys, x, t, h);
        t += h;
    }
    std::ofstream out;
    if (!out_path.empty()) {
        out.open(out_path);
        if (!out) throw std::runtime_error("cannot write " + out_path);
        out << "t";
        for (int j = 0; j < sys.dim; ++j) out << ",x" << j;
        out << "\n";
        out << std::setprecision(12);
    }
    for (std::size_t i = 0; i < n; ++i) {
        if (i % stride == 0 && out) {
            out << t;
            for (double v : x) out << "," << v;
            out << "\n";
        }
        x = rk4_step(sys, x, t, h);
        t += h;
    }
    if (!out) {
        std::cout << sys.name << " done; last state:";
        for (double v : x) std::cout << " " << v;
        std::cout << "\n";
    }
}

void cmd_lyapunov(const Options& opt) {
    if (opt.has("map")) {
        MapSystem map = choose_map(opt);
        State x = opt.has("x0") ? parse_state(opt.get("x0")) : default_initial_state(map.name);
        const std::size_t transient = opt.get_size("transient", 1000);
        const std::size_t n = opt.get_size("n", 200000);
        const LyapunovSpectrum sp = lyapunov_spectrum_map(map, x, transient, n);
        std::cout << "map " << map.name << " exponents (per iteration):";
        double sum = 0.0;
        for (double v : sp.exponents) {
            std::cout << " " << std::fixed << std::setprecision(6) << v;
            sum += v;
        }
        std::cout << "  sum=" << sum << "\n";
        if (!opt.get("out").empty()) {
            std::ofstream out(opt.get("out"));
            out << "index,exponent\n" << std::setprecision(12);
            for (std::size_t i = 0; i < sp.exponents.size(); ++i) out << i << "," << sp.exponents[i] << "\n";
        }
        return;
    }
    OdeSystem sys = choose_ode(opt);
    State x = choose_state(opt);
    LyapunovOptions lo;
    lo.h = opt.get_double("h", 0.01);
    lo.transient = opt.get_double("transient", 50.0);
    lo.duration = opt.get_double("duration", 500.0);
    const LyapunovSpectrum sp = lyapunov_spectrum(sys, x, lo);
    std::cout << "system " << sys.name << " exponents:";
    double sum = 0.0;
    for (double v : sp.exponents) {
        std::cout << " " << std::fixed << std::setprecision(6) << v;
        sum += v;
    }
    std::cout << "  sum=" << sum << "  T=" << sp.integration_time << "\n";
    if (!opt.get("out").empty()) {
        std::ofstream out(opt.get("out"));
        out << "index,exponent\n" << std::setprecision(12);
        for (std::size_t i = 0; i < sp.exponents.size(); ++i) out << i << "," << sp.exponents[i] << "\n";
    }
}

void cmd_lyapunov_scan(const Options& opt) {
    const std::string name = opt.get("map", "logistic");
    if (name != "logistic") throw std::runtime_error("lyapunov-scan supports --map logistic");
    const double pmin = opt.get_double("pmin", 2.5);
    const double pmax = opt.get_double("pmax", 4.0);
    const int nparams = opt.get_int("nparams", 1200);
    const std::size_t transient = opt.get_size("transient", 1000);
    const std::size_t n = opt.get_size("n", 20000);
    const double x0 = opt.get_double("x0", 0.123456789);
    const std::string path = opt.get("out");
    const bool to_file = !path.empty();
    std::ofstream out;
    if (to_file) {
        out.open(path);
        if (!out) throw std::runtime_error("cannot write " + path);
        out << "r,lambda\n" << std::setprecision(12);
    }
    for (int k = 0; k < nparams; ++k) {
        const double r = pmin + (pmax - pmin) * static_cast<double>(k) /
                                     static_cast<double>(nparams - 1);
        const double lambda = Logistic{r}.lyapunov(x0, transient, n);
        if (to_file) out << r << "," << lambda << "\n";
        else std::cout << r << "," << lambda << "\n";
    }
}

void cmd_bifurcation(const Options& opt) {
    const std::string name = opt.get("map", "logistic");
    if (name != "logistic") throw std::runtime_error("bifurcation currently supports --map logistic");
    const double pmin = opt.get_double("pmin", 2.5);
    const double pmax = opt.get_double("pmax", 4.0);
    const int nparams = opt.get_int("nparams", 1500);
    const std::size_t transient = opt.get_size("transient", 1000);
    const std::size_t keep = opt.get_size("keep", 200);
    const std::size_t stride = opt.get_size("stride", 1);
    const double x0 = opt.get_double("x0", 0.123456789);
    auto f = [](double x, double p) { return p * x * (1.0 - x); };
    const auto pts = bifurcation_scan_1d(f, x0, pmin, pmax, nparams, transient, keep, stride);
    const std::string path = opt.get("out");
    std::ofstream out;
    if (!path.empty()) {
        out.open(path);
        if (!out) throw std::runtime_error("cannot write " + path);
        out << "r,x\n" << std::setprecision(12);
    }
    for (const auto& p : pts) {
        if (out) {
            out << p.parameter << "," << p.value << "\n";
        }
    }
    std::cout << "bifurcation points: " << pts.size() << "\n";
}

void cmd_poincare(const Options& opt) {
    OdeSystem sys = choose_ode(opt);
    State x = choose_state(opt);
    PoincareOptions po;
    po.coordinate = opt.get_int("coord", 2);
    po.value = opt.get_double("value", 0.0);
    po.direction = opt.get_int("direction", 0);
    po.h = opt.get_double("h", 0.005);
    po.t_transient = opt.get_double("transient", 20.0);
    po.t_max = opt.get_double("t", 200.0);
    po.max_points = opt.get_size("max-points", 20000);
    const auto pts = poincare_section(sys, x, po);
    const std::string path = opt.get("out");
    std::ofstream out;
    if (!path.empty()) {
        out.open(path);
        if (!out) throw std::runtime_error("cannot write " + path);
        for (int j = 0; j < sys.dim; ++j) out << (j ? ",x" : "x") << j;
        out << "\n" << std::setprecision(12);
        for (const auto& s : pts) {
            for (int j = 0; j < sys.dim; ++j) out << (j ? "," : "") << s[static_cast<std::size_t>(j)];
            out << "\n";
        }
    }
    std::cout << "Poincare points: " << pts.size() << "\n";
}

std::vector<double> load_series_opts(const Options& opt) {
    const std::string path = opt.get("input");
    if (path.empty()) throw std::runtime_error("--input is required");
    std::vector<double> s = read_series(path, opt.get_int("column", 1), opt.has("header"));
    const std::size_t stride = std::max<std::size_t>(1, opt.get_size("stride", 1));
    if (stride > 1) {
        std::vector<double> sub;
        for (std::size_t i = 0; i < s.size(); i += stride) sub.push_back(s[i]);
        s.swap(sub);
    }
    return s;
}

void cmd_ami(const Options& opt) {
    const auto s = load_series_opts(opt);
    const int maxtau = opt.get_int("max-tau", 60);
    const int bins = opt.get_int("bins", 32);
    const auto ami = average_mutual_information(s, maxtau, bins);
    const int first = first_ami_minimum(ami);
    std::cout << "first AMI minimum: " << first << "\n";
    std::cout << "tau,ami\n" << std::setprecision(12);
    for (std::size_t i = 0; i < ami.size(); ++i) std::cout << i << "," << ami[i] << "\n";
}

void cmd_fnn(const Options& opt) {
    const auto s = load_series_opts(opt);
    const int maxm = opt.get_int("max-m", 10);
    const int tau = opt.get_int("tau", 1);
    const int theiler = opt.get_int("theiler", 0);
    const double rtol = opt.get_double("rtol", 15.0);
    const auto r = false_nearest_neighbors(s, maxm, tau, theiler, rtol);
    std::cout << "suggested embedding dimension: " << r.suggested_m << "\n";
    std::cout << "m,false_ratio\n" << std::setprecision(12);
    for (int m = 1; m <= maxm; ++m) std::cout << m << "," << r.false_ratio[static_cast<std::size_t>(m)] << "\n";
}

void cmd_dimension(const Options& opt) {
    const auto s = load_series_opts(opt);
    const int m = opt.get_int("m", 3);
    const int tau = opt.get_int("tau", 10);
    const int theiler = opt.get_int("theiler", 50);
    const int nr = opt.get_int("nr", 60);
    const double rmin = opt.get_double("rmin", 0.0);
    const double rmax = opt.get_double("rmax", 0.0);
    const double flo = opt.get_double("fit-lo", 0.2);
    const double fhi = opt.get_double("fit-hi", 0.8);
    const auto r = correlation_dimension(s, m, tau, rmin, rmax, nr, theiler, flo, fhi);
    std::cout << "D2=" << r.d2 << "  R2=" << r.r2 << "  r=[" << r.r_min << "," << r.r_max
              << "]  m=" << m << " tau=" << tau << "\n";
    const std::string path = opt.get("out");
    if (!path.empty()) {
        std::ofstream out(path);
        if (!out) throw std::runtime_error("cannot write " + path);
        out << "r,C,log_r,log_C,local_slope\n" << std::setprecision(12);
        for (std::size_t i = 0; i < r.radii.size(); ++i) {
            out << r.radii[i] << "," << std::exp(r.log_c[i]) << "," << r.log_r[i] << ","
                << r.log_c[i];
            if (i < r.local_slope.size() && std::isfinite(r.local_slope[i])) {
                out << "," << r.local_slope[i];
            }
            out << "\n";
        }
    }
}

void cmd_rosenstein(const Options& opt) {
    const auto s = load_series_opts(opt);
    const int m = opt.get_int("m", 3);
    const int tau = opt.get_int("tau", 10);
    const int theiler = opt.get_int("theiler", 100);
    const int maxk = opt.get_int("max-k", 300);
    const double dt = opt.get_double("dt", 0.01);
    const int fs = opt.get_int("fit-start", -1);
    const int fe = opt.get_int("fit-end", -1);
    const auto r = rosenstein_lyapunov(s, m, tau, theiler, maxk, dt, fs, fe);
    std::cout << "lambda_max=" << r.lambda << "  R2=" << r.r2 << "  fit k=["
              << r.k[r.fit_begin] << "," << r.k[r.fit_end] << "]  dt=" << dt << "\n";
    const std::string path = opt.get("out");
    if (!path.empty()) {
        std::ofstream out(path);
        if (!out) throw std::runtime_error("cannot write " + path);
        out << "k,mean_log_divergence\n" << std::setprecision(12);
        for (std::size_t i = 0; i < r.k.size(); ++i) out << r.k[i] << "," << r.divergence[i] << "\n";
    }
}

void cmd_zeroone(const Options& opt) {
    const auto s = load_series_opts(opt);
    const int nc = opt.get_int("nc", 64);
    const auto r = zero_one_test(s, nc);
    std::cout << "K=" << r.K << "\n";
    const std::string path = opt.get("out");
    if (!path.empty()) {
        std::ofstream out(path);
        out << "trial,K\n" << std::setprecision(12);
        for (std::size_t i = 0; i < r.K_values.size(); ++i) out << i << "," << r.K_values[i] << "\n";
    }
}

void cmd_rqa(const Options& opt) {
    const auto s = load_series_opts(opt);
    const int m = opt.get_int("m", 3);
    const int tau = opt.get_int("tau", 10);
    const int theiler = opt.get_int("theiler", 10);
    const Embedding emb = delay_embed(s, m, tau);
    double eps = opt.get_double("epsilon", 0.0);
    if (!(eps > 0.0)) {
        eps = epsilon_for_recurrence_rate(emb, opt.get_double("rate", 0.05), theiler);
    }
    const RqaResult r = rqa(emb, eps, theiler, opt.get_int("min-line", 2));
    std::cout << "epsilon=" << eps << " RR=" << r.recurrence_rate << " DET=" << r.determinism
              << " LAM=" << r.laminarity << " ENTR=" << r.entropy << " TT=" << r.trapping_time
              << " Lmax=" << r.longest_diagonal << "\n";
    const std::string path = opt.get("matrix-out");
    if (!path.empty()) {
        const auto R = recurrence_matrix(emb, eps, theiler);
        std::ofstream out(path);
        if (!out) throw std::runtime_error("cannot write " + path);
        for (std::size_t i = 0; i < emb.n; ++i) {
            for (std::size_t j = 0; j < emb.n; ++j) out << (j ? "," : "") << static_cast<int>(R[i * emb.n + j]);
            out << "\n";
        }
    }
}

void cmd_spectrum(const Options& opt) {
    const auto s = load_series_opts(opt);
    const double dt = opt.get_double("dt", 0.01);
    const Spectrum sp = power_spectrum(s, dt);
    std::cout << "df=" << sp.df << " nfft=" << sp.nfft << "\n";
    std::vector<std::size_t> idx(sp.power.size());
    for (std::size_t i = 0; i < idx.size(); ++i) idx[i] = i;
    std::sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) {
        return sp.power[a] > sp.power[b];
    });
    std::cout << "top peaks (frequency,power):\n";
    for (std::size_t i = 0; i < std::min<std::size_t>(10, idx.size()); ++i) {
        std::cout << "  " << sp.frequency[idx[i]] << " " << sp.power[idx[i]] << "\n";
    }
    const std::string path = opt.get("out");
    if (!path.empty()) {
        std::ofstream out(path);
        out << "frequency,power\n" << std::setprecision(12);
        for (std::size_t i = 0; i < sp.frequency.size(); ++i) out << sp.frequency[i] << "," << sp.power[i] << "\n";
    }
}

void cmd_sindy(const Options& opt) {
    OdeSystem sys = choose_ode(opt);
    State x = choose_state(opt);
    const double h = opt.get_double("h", 0.002);
    const std::size_t n = opt.get_size("n", 5000);
    const std::size_t skip = opt.get_size("transient", 5000);
    for (std::size_t i = 0; i < skip; ++i) x = rk4_step(sys, x, static_cast<double>(i) * h, h);
    std::vector<State> X;
    std::vector<State> dX;
    double t = static_cast<double>(skip) * h;
    for (std::size_t i = 0; i < n; ++i) {
        X.push_back(x);
        dX.push_back(sys.rhs(x, t));
        x = rk4_step(sys, x, t, h);
        t += h;
    }
    SindyOptions so;
    so.degree = opt.get_int("degree", 2);
    so.threshold = opt.get_double("threshold", 0.1);
    so.include_constant = !opt.has("no-constant");
    const SindyModel model = sindy_fit(X, dX, so);
    std::cout << model.pretty();
    const std::string path = opt.get("out");
    if (!path.empty()) {
        std::ofstream out(path);
        out << "term";
        for (std::size_t j = 0; j < model.coefficients.cols; ++j) out << ",d" << j;
        out << "\n" << std::setprecision(12);
        for (std::size_t i = 0; i < model.terms.size(); ++i) {
            out << model.terms[i].name;
            for (std::size_t j = 0; j < model.coefficients.cols; ++j) out << "," << model.coefficients(i, j);
            out << "\n";
        }
    }
}

void cmd_sync(const Options& opt) {
    const double h = opt.get_double("h", 0.005);
    const double tt = opt.get_double("t", 40.0);
    const double trans = opt.get_double("transient", 5.0);
    const State d0 = opt.has("x0") ? parse_state(opt.get("x0")) : State{1.0, 1.0, 1.0};
    const State r0 = opt.has("response0") ? parse_state(opt.get("response0")) : State{0.0, 0.0};
    const SyncResult r = pecora_carroll_lorenz(Lorenz{}, d0, r0, h, trans, tt);
    std::cout << "synchronized=" << (r.synchronized ? "true" : "false")
              << " final_error=" << r.final_error
              << " max_error_after_transient=" << r.max_error_after_transient
              << " settling_time=" << r.settling_time << "\n";
    const std::string path = opt.get("out");
    if (!path.empty()) {
        std::ofstream out(path);
        out << "t,error\n" << std::setprecision(12);
        for (std::size_t i = 0; i < r.t.size(); ++i) out << r.t[i] << "," << r.error[i] << "\n";
    }
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        usage();
        return 0;
    }
    const std::string command = argv[1];
    try {
        const Options opt = parse_options(argc, argv);
        if (command == "models") {
            cmd_models();
        } else if (command == "simulate") {
            cmd_simulate(opt);
        } else if (command == "lyapunov") {
            cmd_lyapunov(opt);
        } else if (command == "lyapunov-scan") {
            cmd_lyapunov_scan(opt);
        } else if (command == "bifurcation") {
            cmd_bifurcation(opt);
        } else if (command == "poincare") {
            cmd_poincare(opt);
        } else if (command == "ami") {
            cmd_ami(opt);
        } else if (command == "fnn") {
            cmd_fnn(opt);
        } else if (command == "dimension") {
            cmd_dimension(opt);
        } else if (command == "rosenstein") {
            cmd_rosenstein(opt);
        } else if (command == "zeroone") {
            cmd_zeroone(opt);
        } else if (command == "rqa") {
            cmd_rqa(opt);
        } else if (command == "spectrum") {
            cmd_spectrum(opt);
        } else if (command == "sindy") {
            cmd_sindy(opt);
        } else if (command == "sync") {
            cmd_sync(opt);
        } else if (command == "version" || command == "--version") {
            std::cout << CHAOS_VERSION_STRING << "\n";
        } else if (command == "help" || command == "--help" || command == "-h") {
            usage();
        } else {
            std::cerr << "unknown command: " << command << "\n";
            usage();
            return 2;
        }
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}

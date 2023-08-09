#include <algorithm>
#include <chrono>
#include <functional>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <random>
#include <cmath>
#include <vector>
//---------------------------------------------------------------
constexpr auto deg = 2;
constexpr auto N   = 1000;
constexpr auto dt  = 5e-3;
//---------------------------------------------------------------
constexpr auto rho = 1.0;
constexpr auto Ndof = deg*N;

constexpr auto margin = 0.4;
constexpr auto SKIN2  = (margin*0.5) * (margin*0.5);
//---------------------------------------------------------------
constexpr auto eta  = 0.2;      // non-additivity
constexpr auto smax = 1.62;     // sigma_max
constexpr auto smin = 0.73;     // sigma_min
constexpr auto Ap   = 1.0/(0.5/(smin*smin) - 0.5/(smax*smax)); // normalization for polydispersity

constexpr auto xc   = 1.25;// cut-off length
constexpr auto xc2  = xc*xc;
constexpr auto xc12 = xc2 * xc2 * xc2 * xc2 * xc2 * xc2;
constexpr auto c4  = -21.0/(xc12 * xc2 * xc2);
constexpr auto c2  = 48.0/(xc12 * xc2);
constexpr auto c0  = -28.0/xc12;
//---------------------------------------------------------------
const double Lbox = std::pow(N/rho, 1.0/deg);
const double Linv = 1.0/Lbox;
const auto Ncell = int(Lbox/(xc*smax+margin));

double conf[N][deg], velo[N][deg], NL_config[N][deg], force[N][deg], sigma[N];
std::vector<int> point(N), list(N*80), head(Ncell*Ncell), linklist(N);
int pairs[N*80][2];
double vxi1 = 0.0;

enum {X, Y};
//---------------------------------------------------------------
void init_lattice() {
    const auto ln   = std::ceil(std::pow(N, 1.0/deg));
    const auto haba = Lbox/ln;

    for (int i=0; i<N; i++) {
        int iy = std::floor(i/ln);
        int ix = i - iy*ln;

        conf[i][X] = haba*0.5 + haba * ix;
        conf[i][Y] = haba*0.5 + haba * iy;

        for (int d=0; d<deg; d++) {
            conf[i][d] -= Lbox * std::round(conf[i][d] * Linv);
        }
    }
}
void init_sigma(std::mt19937 &mt) {
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    for (int i=0; i<N; i++) {
        sigma[i] = 1.0/std::sqrt(1.0/(smin*smin) - 2.0*dist(mt)/Ap);
    }
}
constexpr auto Ninv  = 1.0/N;
inline void remove_drift() {
    double vel1 = 0.0, vel2 = 0.0;
    for (int i=0; i<N; i++) {
        vel1 += velo[i][X];
        vel2 += velo[i][Y];
    }
    vel1 *= Ninv;
    vel2 *= Ninv;
    for (int i=0; i<N; i++) {
        velo[i][X] -= vel1;
        velo[i][Y] -= vel2;
    }
}
void init_vel_MB(const double T_targ, std::mt19937 &mt) {
    std::normal_distribution<double> dist_trans(0.0, std::sqrt(T_targ));
    for (int i=0; i<N; i++) {
        velo[i][X] = dist_trans(mt);
        velo[i][Y] = dist_trans(mt);
    }
    remove_drift();
}
//---------------------------------------------------------------
const auto Lcell = Lbox/Ncell;
const auto Lcellinv = 1.0/Lcell;
const auto halfLbox = Lbox*0.5;
void new_llist() {
    std::fill(head.begin(), head.end(), -1);
    for (int i=0; i<N; i++) {
        const double xi = conf[i][X];
        const double yi = conf[i][Y];
        const int icell = int((xi + halfLbox)*Lcellinv) + Ncell*int((yi + halfLbox)*Lcellinv);
        linklist[i] = head[icell];
        head[icell] = i;
    }
}
constexpr auto Nend = N-1;
void counting_sort(const int count) {
    // count
    std::fill(linklist.begin(), linklist.end(), 0);
    for (int i=0; i<count; i++) {
        linklist[pairs[i][0]]++;
    }
    // cumsum
    point[0] = 0;
    int sum = 0;
    for (int i=0; i<Nend; i++) {
        sum += linklist[i];
        point[i+1] = sum;
    }
    // build list
    std::fill(linklist.begin(), linklist.end(), 0);
    for (int i=0; i<count; i++) {
        int pos = point[pairs[i][0]] + linklist[pairs[i][0]];
        list[pos] = pairs[i][1];
        linklist[pairs[i][0]]++;
    }
    // sort within pair particles j of i
    // Use this if you want to match the order of list with O(N^2) generate_NL.
    // for (int i=0; i<Nend; i++) {
    //     std::sort(&list[point[i]], &list[point[i+1]]);
    // }
}
//---------------------------------------------------------------
void search_other_cell(int ix, int iy, const int icell, int& nlist) {
    if (ix < 0) ix += Ncell;
    if (ix >= Ncell) ix -= Ncell;
    if (iy < 0) iy += Ncell;
    if (iy >= Ncell) iy -= Ncell;
    const int id2 = ix + Ncell*iy;

    int i = head[icell];
    while (i >= 0) {
        const double xi = conf[i][X];
        const double yi = conf[i][Y];
        const double si = sigma[i];

        int j = head[id2];
        while (j >= 0) {
            double dx = xi - conf[j][X];
            double dy = yi - conf[j][Y];
            dx -= Lbox * floor(dx * Linv + 0.5);
            dy -= Lbox * floor(dy * Linv + 0.5);
            const double rij2 = dx*dx + dy*dy;

            const double sj = sigma[j];
            const double sij1 = 0.5 * (si + sj) * (1.0 - eta * std::abs(si - sj));
            double rlist1 = 1.25 * sij1 + margin;
            double rlist2 = rlist1 * rlist1;
            if (rij2 < rlist2) {
                if (i < j) {
                    pairs[nlist][0] = i;
                    pairs[nlist][1] = j;
                } else {
                    pairs[nlist][0] = j;
                    pairs[nlist][1] = i;
                }
                nlist++;
            }
            j = linklist[j];
        }
        i = linklist[i];
    }
}
void generate_NL_LL() {
    new_llist();
    int count = 0;
    for (int icellx=0; icellx<Ncell; icellx++) {
        for (int icelly=0; icelly<Ncell; icelly++) {
            const int icell = icellx + Ncell*icelly;

            search_other_cell(icellx-1, icelly+1, icell, count);
            search_other_cell(icellx,   icelly+1, icell, count);
            search_other_cell(icellx+1, icelly+1, icell, count);

            search_other_cell(icellx+1, icelly,   icell, count);

            // search the same cell
            int i = head[icell];
            while (i >= 0) {
                const double xi = conf[i][X];
                const double yi = conf[i][Y];
                const double si = sigma[i];

                int j = head[icell];
                while (j >= 0) {
                    if (i >= j) {
                        j = linklist[j];
                        continue;
                    }
                    double dx = xi - conf[j][X];
                    double dy = yi - conf[j][Y];
                    dx -= Lbox * floor(dx * Linv + 0.5);
                    dy -= Lbox * floor(dy * Linv + 0.5);
                    const double rij2 = dx*dx + dy*dy;

                    const double sj = sigma[j];
                    const double sij1 = 0.5 * (si + sj) * (1.0 - eta * std::abs(si - sj));
                    double rlist1 = 1.25 * sij1 + margin;
                    double rlist2 = rlist1 * rlist1;
                    if (rij2 < rlist2) {
                        pairs[count][0] = i;
                        pairs[count][1] = j;
                        count++;
                    }
                    j = linklist[j];
                }
                i = linklist[i];
            }
        }
    }
    counting_sort(count);
    std::copy(*conf, *conf+Ndof, *NL_config);
}
//---------------------------------------------------------------
void calc_force() {
    std::fill(*force, *force+Ndof, 0.0);
    for (int i=0; i<Nend; i++) {
        const int pstart = point[i];
        const int pend = point[i+1];
        if (pstart == pend) continue;

        const double xi = conf[i][X];
        const double yi = conf[i][Y];
        const double si = sigma[i];
        double fix = force[i][X];
        double fiy = force[i][Y];
        for (int p=pstart; p<pend; p++) {
            int j = list[p];
            double dx = xi - conf[j][X];
            double dy = yi - conf[j][Y];
            dx -= Lbox * floor(dx * Linv + 0.5);
            dy -= Lbox * floor(dy * Linv + 0.5);
            const double rij2 = dx*dx + dy*dy;

            const double sj = sigma[j];
            const double sij1 = 0.5 * (si + sj) * (1.0 - eta * std::abs(si - sj));
            const double sij2 = sij1 * sij1;
            if (rij2 < xc2 * sij2) {
                const double rij6 = rij2 * rij2 * rij2;
                const double rij12 = rij6 * rij6;
                const double rij14 = rij12 * rij2;
                const double sij4 = sij2 * sij2;
                const double sij12 = sij4 * sij4 * sij4;

                const double temp = (-12.0 * sij12 * sij4 + 2.0 * c2 * rij14 * sij2 + 4.0 * c4 * rij14 * rij2) / (rij14 * sij4);
                fix -= temp * dx;
                fiy -= temp * dy;
                force[j][X] += temp * dx;
                force[j][Y] += temp * dy;
            }
        }
        force[i][X] = fix;
        force[i][Y] = fiy;
    }
}
//---------------------------------------------------------------
double calc_potential() {
    double ans = 0.0;
    for (int i=0; i<Nend; i++) {
        const int pend = point[i+1];
        const double xi = conf[i][X];
        const double yi = conf[i][Y];
        const double si = sigma[i];
        for (int p=point[i]; p<pend; p++) {
            const int j  = list[p];
            double dx = xi - conf[j][X];
            double dy = yi - conf[j][Y];
            dx -= Lbox * floor(dx * Linv + 0.5);
            dy -= Lbox * floor(dy * Linv + 0.5);
            const double rij2 = dx*dx + dy*dy;

            const double sj = sigma[j];
            const double sij1 = 0.5 * (si + sj) * (1.0 - eta * std::abs(si - sj));
            const double sij2 = sij1 * sij1;
            if (rij2 < xc2 * sij2) {
                const double rij6 = rij2 * rij2 * rij2;
                const double rij12 = rij6 * rij6;
                const double sij4 = sij2 * sij2;
                const double sij12 = sij4 * sij4 * sij4;
                ans += c0 + rij2 * (c2 * sij2 + c4 * rij2) / sij4 + sij12/rij12;
            }
        }
    }
    return ans;
}
double calc_potential_N2() {
    double ans = 0.0;
    for (int i=0; i<N; i++) {
        const double xi = conf[i][X];
        const double yi = conf[i][Y];
        const double si = sigma[i];
        for (int j=i+1; j<N; j++) {
            double dx = xi - conf[j][X];
            double dy = yi - conf[j][Y];
            dx -= Lbox * floor(dx * Linv + 0.5);
            dy -= Lbox * floor(dy * Linv + 0.5);
            const double rij2 = dx*dx + dy*dy;

            const double sj = sigma[j];
            const double sij1 = 0.5 * (si + sj) * (1.0 - eta * std::abs(si - sj));
            const double sij2 = sij1 * sij1;
            if (rij2 < xc2 * sij2) {
                const double rij6 = rij2 * rij2 * rij2;
                const double rij12 = rij6 * rij6;
                const double sij4 = sij2 * sij2;
                const double sij12 = sij4 * sij4 * sij4;
                ans += c0 + rij2 * (c2 * sij2 + c4 * rij2) / sij4 + sij12/rij12;
            }
        }
    }
    return ans;
}
//---------------------------------------------------------------
constexpr auto dt2 = dt*0.5;
constexpr auto dt4 = dt*0.25;
inline void velocity_update() {
    for (int i=0; i<N; i++) {
        velo[i][X] += dt2*force[i][X];
        velo[i][Y] += dt2*force[i][Y];
    }
}
inline void position_update() {
    for (int i=0; i<N; i++) {
        conf[i][X] += dt*velo[i][X];
        conf[i][Y] += dt*velo[i][Y];
    }
}
inline void PBC() {
    for (int i=0; i<N; i++) {
        conf[i][X] -= Lbox * floor(conf[i][X] * Linv + 0.5);
        conf[i][Y] -= Lbox * floor(conf[i][Y] * Linv + 0.5);
    }
}
inline void NL_check() {
    double dev_max = 0.0;
    for (int i=0; i<N; i++) {
        double xij = conf[i][X] - NL_config[i][X];
        double yij = conf[i][Y] - NL_config[i][Y];
        xij -= Lbox * floor(xij * Linv + 0.5);
        yij -= Lbox * floor(yij * Linv + 0.5);
        double rij2 = xij*xij + yij*yij;
        if (rij2 > dev_max) dev_max = rij2;
    }
    if (dev_max > SKIN2) {// renew neighbor list
        generate_NL_LL();
    }
}
//---------------------------------------------------------------
void print_log(long t) {
    double K = 0.5*std::inner_product(*velo, *velo+Ndof, *velo, 0.0);
    double U = calc_potential();
    std::cout << std::setprecision(6) << std::scientific
              << dt*t << ","
              << K*Ninv << ","
              << U*Ninv << ","
              << (K+U)*Ninv << std::endl;
}
void NVE(const double tsim) {
    calc_force();
    // for logging ///////////////////////////////////
    const auto logbin = std::pow(10.0, 1.0/9);
    int counter = 5;
    auto checker = 1e-3 * std::pow(logbin, counter);
    //////////////////////////////////////////////////

    long t = 0;
    const long steps = tsim/dt;
    while (t < steps) {
        velocity_update();
        position_update();
        PBC();
        NL_check();
        calc_force();
        velocity_update();

        t++;
        if (dt*t > checker) {
            checker *= logbin;
            print_log(t);
        }
    }
}
//---------------------------------------------------------------
void NVT(const double T_targ, const double tsim) {
    calc_force();
    // Nose-Hoover variable
    const auto gkBT = Ndof*T_targ;

    long t = 0;
    const long steps = tsim/dt;
    while (t < steps) {
        // Nose-Hoover chain (QMASS = 1.0)
        double uk = std::inner_product(*velo, *velo+Ndof, *velo, 0.0);
        vxi1 += dt4 * (uk - gkBT);
        double temp = std::exp(-vxi1 * dt2);
        std::transform(*velo, *velo+Ndof, *velo,
                       std::bind(std::multiplies<double>(), std::placeholders::_1, temp));
        vxi1 += dt4 * (uk*temp*temp - gkBT);

        velocity_update();
        position_update();
        PBC();
        NL_check();
        calc_force();
        velocity_update();

        // Nose-Hoover chain (QMASS = 1.0)
        uk    = std::inner_product(*velo, *velo+Ndof, *velo, 0.0);
        vxi1 += dt4 * (uk - gkBT);
        temp  = std::exp(-vxi1 * dt2);
        std::transform(*velo, *velo+Ndof, *velo,
                       std::bind(std::multiplies<double>(), std::placeholders::_1, temp));
        vxi1 += dt4 * (uk*temp*temp - gkBT);

        t++;
        if (!(t & 127)) remove_drift();
    }
}
//---------------------------------------------------------------
void initialize_simulation(int seed) {
    std::fill(*conf, *conf+Ndof, 0.0);
    std::fill(*velo, *velo+Ndof, 0.0);

    std::mt19937 mt(seed);
    init_lattice();
    init_vel_MB(0.2, mt);
    init_sigma(mt);
    generate_NL_LL();
    NVT(0.2, 1e3);
}
int main() {
    constexpr int random_seed = 123456789;
    constexpr double target_temperature = 0.2;
    constexpr double simulation_duration = 1e3;

    initialize_simulation(random_seed);
    auto start = std::chrono::system_clock::now();
    NVT(target_temperature, simulation_duration);
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
    std::cout << elapsed << std::endl;

    // for debug
    // NVE(1e5);
}

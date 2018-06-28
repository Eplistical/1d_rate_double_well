#ifndef _BCME_PARTICLE_1D_HPP
#define _BCME_PARTICLE_1D_HPP

#include <algorithm>
#include "cme_particle_1d.hpp"
#include "misc/randomer.hpp"
#include "misc/fermi.hpp"

namespace {

    template <typename PotentialType, typename InttableMgrType>
        struct BCME_Particle_1D final : public CME_Particle_1D<PotentialType>
    {
        public:
            using potential_t = PotentialType;
            using inttable_mgr_t = InttableMgrType;

        public:
            BCME_Particle_1D(double X, double V, double MASS, double KT,
                    int SURFACE,
                    double NUCLEAR_FRIC,
                    const potential_t& POTENTIAL,
                    const inttable_mgr_t& INTTABLE_MGR) noexcept :
                CME_Particle_1D<PotentialType>(X, V, MASS, KT, SURFACE, NUCLEAR_FRIC, POTENTIAL),
                inttable_mgr(INTTABLE_MGR), 
                int_gamma_t(0.0)
                {
                }

            ~BCME_Particle_1D() noexcept = default;

        public:
            double get_N() const noexcept {
                double f(misc::fermi(this->potential.cal_h(this->x) * this->kT_inv));
                double n(inttable_mgr.retrieve("n", this->x));
                return this->surface + (n - f) * (1.0 - exp(-int_gamma_t));
            }

        private:
            inline double cal_force(double x, int surf) const {
                double f(misc::fermi(this->potential.cal_h(x) * this->kT_inv));
                double dhdx(this->potential.cal_dhdx(x));
                double fBCME(dhdx * f + inttable_mgr.retrieve("force", x));
                return this->potential.cal_force(x, surf) + fBCME;
            }

            inline double cal_fric(double x) const {
                // unbroadened friction
                double f(misc::fermi(this->potential.cal_h(x) * this->kT_inv));
                double dhdx(this->potential.cal_dhdx(x));

                double unbroadened_elefric = this->kT_inv / this->potential.cal_gamma(x) * f * (1.0 - f) * pow(dhdx, 2);
                double broadened_elefric = inttable_mgr.retrieve("fric", x);
                double elefric = (broadened_elefric - unbroadened_elefric);

                return std::max(this->nuclear_fric + elefric, 0.0);
            }

            inline double rk4_f(double x, double v, int surf, double fric, double noise) const {
                return v;
            }

            inline double rk4_g(double x, double v, int surf, double fric, double noise) const {
                return (cal_force(x, surf) - fric * v + noise) * this->mass_inv;
            }

            void do_evolve(double dt, const std::string& alg) override {
                // hopping
                this->hopper(dt);
                // integral {gamma(x(t)) * dt}
                int_gamma_t += this->potential.cal_gamma(this->x) * dt;

                if (alg == "verlet") {
                    const double half_dt_mass_inv(0.5 * dt / this->mass);
                    double fric, noise_force;

                    fric = cal_fric(this->x);
                    noise_force = randomer::normal(0.0, sqrt(2.0 * fric * this->kT / dt));
                    assert(fric >= 0.0);

                    this->v += (cal_force(this->x, this->surface) - fric * this->v + noise_force) * half_dt_mass_inv;
                    this->x += this->v * dt;
                    this->v += (cal_force(this->x, this->surface) - fric * this->v + noise_force) * half_dt_mass_inv;
                }
                else if (alg == "rk4") {
                    const double half_dt(0.5 * dt);
                    double k1, k2, k3, k4;
                    double l1, l2, l3, l4;
                    double fric, noise_force;

                    fric = cal_fric(this->x);
                    noise_force = randomer::normal(0.0, sqrt(2.0 * fric * this->kT / dt));

                    k1 = rk4_f(this->x, this->v, this->surface, fric, noise_force);
                    l1 = rk4_g(this->x, this->v, this->surface, fric, noise_force);

                    k2 = rk4_f(this->x + half_dt * k1, this->v + half_dt * l1, this->surface, fric, noise_force);
                    l2 = rk4_g(this->x + half_dt * k1, this->v + half_dt * l1, this->surface, fric, noise_force);

                    k3 = rk4_f(this->x + half_dt * k2, this->v + half_dt * l2, this->surface, fric, noise_force);
                    l3 = rk4_g(this->x + half_dt * k2, this->v + half_dt * l2, this->surface, fric, noise_force);

                    k4 = rk4_f(this->x + dt * k3, this->v + dt * l3, this->surface, fric, noise_force);
                    l4 = rk4_g(this->x + dt * k3, this->v + dt * l3, this->surface, fric, noise_force);

                    this->x += dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4);
                    this->v += dt / 6.0 * (l1 + 2 * l2 + 2 * l3 + l4);
                }
            }

        private:
            const inttable_mgr_t& inttable_mgr;
            double int_gamma_t;
    };

};

#endif // _BCME_PARTICLE_1D_CPP

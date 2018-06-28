#ifndef _BFP_PARTICLE_1D_HPP
#define _BFP_PARTICLE_1D_HPP

#include "particle_1d_base.hpp"
#include "misc/randomer.hpp"
#include "misc/fermi.hpp"

namespace {

    template <typename PotentialType, typename InttableMgrType>
        struct BFP_Particle_1D final : public Particle_1D
    {
        public:
            using potential_t = PotentialType;
            using inttable_mgr_t = InttableMgrType;

        public:
            BFP_Particle_1D(double X, double V, double MASS, double KT,
                    double NUCLEAR_FRIC,
                    const potential_t& POTENTIAL,
                    const inttable_mgr_t& INTTABLE_MGR) noexcept :
                Particle_1D(X, V, MASS, KT),
                nuclear_fric(NUCLEAR_FRIC), 
                potential(POTENTIAL), 
                inttable_mgr(INTTABLE_MGR)
            {
            }

            ~BFP_Particle_1D() noexcept = default;

        public:
            // population on impurity
            double get_N() const {
                return inttable_mgr.retrieve("n", this->x);
            }

        private:
            inline double cal_force(double x) const {
                return potential.cal_force(x, 0) + inttable_mgr.retrieve("force", x);
            }

            inline double cal_fric(double x) const {
                // unbroadened friction
                /*
                double f(misc::fermi(potential.cal_h(x) * this->kT_inv));
                double dhdx(potential.cal_dhdx(x));
                return nuclear_fric + this->kT_inv / potential.cal_gamma(x) * f * (1.0 - f) * pow(dhdx, 2);
				*/
                return nuclear_fric + inttable_mgr.retrieve("fric", x);
            }

            inline double rk4_f(double x, double v, double fric, double noise) const {
                return v;
            }

            inline double rk4_g(double x, double v, double fric, double noise) const {
                return (cal_force(x) - fric * v + noise) * this->mass_inv;
            }

            void do_evolve(double dt, const std::string& alg) override {
                if (alg == "verlet") {
                    const double half_dt_mass_inv(0.5 * dt / this->mass);
                    double fric, noise;

                    fric = cal_fric(this->x);
                    noise = randomer::normal(0.0, sqrt(2.0 * fric * this->kT / dt));
                    assert(fric >= 0.0);

                    this->v += (cal_force(this->x) - fric * this->v + noise) * half_dt_mass_inv;
                    this->x += this->v * dt;
                    this->v += (cal_force(this->x) - fric * this->v + noise) * half_dt_mass_inv;
                }
                else if (alg == "rk4") {
                    const double half_dt(0.5 * dt);
                    double k1, k2, k3, k4;
                    double l1, l2, l3, l4;
                    double fric, noise;

                    fric = cal_fric(this->x);
                    noise = randomer::normal(0.0, sqrt(2.0 * fric * this->kT / dt));

                    k1 = rk4_f(this->x, this->v, fric, noise);
                    l1 = rk4_g(this->x, this->v, fric, noise);

                    k2 = rk4_f(this->x + half_dt * k1, this->v + half_dt * l1, fric, noise);
                    l2 = rk4_g(this->x + half_dt * k1, this->v + half_dt * l1, fric, noise);

                    k3 = rk4_f(this->x + half_dt * k2, this->v + half_dt * l2, fric, noise);
                    l3 = rk4_g(this->x + half_dt * k2, this->v + half_dt * l2, fric, noise);

                    k4 = rk4_f(this->x + dt * k3, this->v + dt * l3, fric, noise);
                    l4 = rk4_g(this->x + dt * k3, this->v + dt * l3, fric, noise);

                    this->x += dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4);
                    this->v += dt / 6.0 * (l1 + 2 * l2 + 2 * l3 + l4);
                }
            }

        public:
            double nuclear_fric;
            const potential_t& potential;
            const inttable_mgr_t& inttable_mgr;
    };

};

#endif // _BFP_PARTICLE_CPP

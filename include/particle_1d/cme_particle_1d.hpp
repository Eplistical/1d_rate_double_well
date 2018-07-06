#ifndef _CME_PARTICLE_1D_HPP
#define _CME_PARTICLE_1D_HPP

#include <string>
#include "particle_1d_base.hpp"
#include "misc/randomer.hpp"
#include "misc/fermi.hpp"

namespace {
    using std::mt19937;

    template <typename PotentialType>
        struct CME_Particle_1D : public Particle_1D
    {
        public:
            using potential_t = PotentialType;

        public:
            CME_Particle_1D(double X, double V, 
                    double MASS, double KT,
                    int SURFACE,
                    double NUCLEAR_FRIC,
                    const potential_t& POTENTIAL) noexcept :
                Particle_1D(X, V, MASS, KT),
                surface(SURFACE), nuclear_fric(NUCLEAR_FRIC),
                potential(POTENTIAL)
            {
            }

            virtual ~CME_Particle_1D() noexcept = default;

        public:
            // population on impurity
            double get_N() const noexcept {
                return surface;
            }

        private:
            inline double cal_force(double x, int surf) const {
                return potential.cal_force(x, surf);
            }

            inline double cal_fric(double x) const {
                return nuclear_fric;
            }

            inline double rk4_f(double x, double v, int surf, double fric, double noise) const {
                return v;
            }

            inline double rk4_g(double x, double v, int surf, double fric, double noise) const {
                return (cal_force(x, surf) - fric * v + noise) * this->mass_inv;
            }

            virtual void do_evolve(double dt, const std::string& alg) override {
                // hopping
                hopper(dt);

                if (alg == "verlet") {
                    const double half_dt_mass_inv(0.5 * dt * this->mass_inv);
                    double fric, noise;

                    fric = cal_fric(this->x);
                    noise = randomer::normal(0.0, sqrt(2.0 * fric * this->kT / dt));
                    assert(fric >= 0.0);

                    this->v += (cal_force(this->x, surface) - fric * this->v + noise) * half_dt_mass_inv;
                    this->x += this->v * dt;
                    this->v += (cal_force(this->x, surface) - fric * this->v + noise) * half_dt_mass_inv;
                }
                else if (alg == "rk4") {
                    const double half_dt(0.5 * dt);
                    double k1, k2, k3, k4;
                    double l1, l2, l3, l4;
                    double fric, noise;

                    fric = cal_fric(this->x);
                    noise = randomer::normal(0.0, sqrt(2.0 * fric * this->kT / dt));

                    k1 = rk4_f(this->x, this->v, surface, fric, noise);
                    l1 = rk4_g(this->x, this->v, surface, fric, noise);

                    k2 = rk4_f(this->x + half_dt * k1, this->v + half_dt * l1, surface, fric, noise);
                    l2 = rk4_g(this->x + half_dt * k1, this->v + half_dt * l1, surface, fric, noise);

                    k3 = rk4_f(this->x + half_dt * k2, this->v + half_dt * l2, surface, fric, noise);
                    l3 = rk4_g(this->x + half_dt * k2, this->v + half_dt * l2, surface, fric, noise);

                    k4 = rk4_f(this->x + dt * k3, this->v + dt * l3, surface, fric, noise);
                    l4 = rk4_g(this->x + dt * k3, this->v + dt * l3, surface, fric, noise);

                    this->x += dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4);
                    this->v += dt / 6.0 * (l1 + 2 * l2 + 2 * l3 + l4);
                }
            }

        protected:
            void hopper(double dt) {
                const double gamma(potential.cal_gamma(this->x));
                const double h(potential.cal_h(this->x));
                double prob;
                if (surface == 0) {
                    prob = gamma * misc::fermi(h * this->kT_inv) * dt;
                    if (randomer::rand() < prob) {
                        surface = 1;
                    }
                }
                else if (surface == 1) {
                    prob = gamma * (1.0 - misc::fermi(h * this->kT_inv)) * dt;
                    if (randomer::rand() < prob) {
                        surface = 0;
                    }
                }
            }
            
        public:
            int surface;
            double nuclear_fric;
            const potential_t& potential;
    };

};

#endif // _CME_PARTICLE_1D_HPP

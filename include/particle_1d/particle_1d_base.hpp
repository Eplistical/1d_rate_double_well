#ifndef _PARTICLE_1D_BASE_HPP
#define _PARTICLE_1D_BASE_HPP
// module for 1d particle base class (abstract class)

#include <string>

namespace {

    struct Particle_1D
    {
        public:
            Particle_1D(   double X, double V, double MASS, double KT) noexcept :
                x(X), v(V), mass(MASS), kT(KT), mass_inv(1.0 / MASS), kT_inv(1.0 / KT)
                {
                }

            virtual ~Particle_1D() noexcept = default;

        public:
            double get_Ek() const noexcept { 
                return 0.5 * mass * v * v; 
            }

        public:
            void evolve(double dt, const std::string& alg = "verlet") {
                do_evolve(dt, alg); 
            }

        public:
            const double mass, mass_inv;
            const double kT, kT_inv;
            double x, v;

        private:
            virtual void do_evolve(double dt, const std::string& alg) = 0;
    };
};

#endif // _PARTICLE_1D_BASE_HPP

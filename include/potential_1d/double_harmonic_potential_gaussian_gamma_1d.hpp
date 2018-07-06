#ifndef _DOUBLE_HARMONIC_POTENTIAL_GAUSSIAN_GAMMA_1D_HPP
#define _DOUBLE_HARMONIC_POTENTIAL_GAUSSIAN_GAMMA_1D_HPP

#include "potential_1d_base.hpp"

namespace {

    struct DoubleHarmonicPotentialGaussianGamma_1D final 
        : public Potential_1D
    {
        public:
            DoubleHarmonicPotentialGaussianGamma_1D(double MASS, double OMEGA, 
                    double HSHIFT, double VSHIFT, double GAMMA, double ALPHA) noexcept :
                mass(MASS), omega(OMEGA), mw2(MASS * OMEGA * OMEGA),
                hshift(HSHIFT), vshift(VSHIFT),
                gamma(GAMMA), alpha(ALPHA)
                {
                }

            ~DoubleHarmonicPotentialGaussianGamma_1D() noexcept = default;

        private:
            double cal_potential_impl(double x, int i) const noexcept override {
                if (i == 0)
                    return 0.5 * mw2 * x * x;
                else
                    return 0.5 * mw2 * (x - hshift) * (x - hshift) + vshift;
            }

            double cal_force_impl(double x, int i) const noexcept override {
                if (i == 0) 
                    return -mw2 *  x;
                else
                    return -mw2 * (x - hshift);
            }

            double cal_gamma_impl(double x) const noexcept override {
                static const double amw(alpha * mass * omega);
                return gamma * (1.0 + exp(-amw * x * x));
            }

            double cal_nabla_gamma_impl(double x) const noexcept override {
                static const double amw(alpha * mass * omega);
                return -2 * amw * x * gamma * exp(-amw * x * x);
            }

        public:
            double cal_h(double x) const noexcept {
                return 0.5 * mw2 * hshift * (hshift - 2 * x) + vshift;
            }

            double cal_dhdx(double x) const noexcept {
                static const double dhdx(-mw2 * hshift);
                return dhdx;
            }

        private:
            const double mass, omega, hshift, vshift;
            const double alpha;
            const double mw2;
            const double gamma;
    };

};
#endif // _DOUBLE_HARMONIC_POTENTIAL_GAUSSIAN_GAMMA_1D_HPP

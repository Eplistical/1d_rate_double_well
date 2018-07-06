#ifndef _DOUBLE_HARMONIC_POTENTIAL_CONST_GAMMA_1D_HPP
#define _DOUBLE_HARMONIC_POTENTIAL_CONST_GAMMA_1D_HPP

#include "potential_1d_base.hpp"

namespace {

    struct DoubleHarmonicPotentialConstGamma_1D final 
        : public Potential_1D
    {
        public:
            DoubleHarmonicPotentialConstGamma_1D(double MASS, double OMEGA, 
                    double HSHIFT, double VSHIFT, double GAMMA) noexcept :
                mass(MASS), omega(OMEGA), mw2(MASS * OMEGA * OMEGA),
                hshift(HSHIFT), vshift(VSHIFT),
                gamma(GAMMA)
                {
                }

            ~DoubleHarmonicPotentialConstGamma_1D() noexcept = default;

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
                return gamma;
            }

            double cal_nabla_gamma_impl(double x) const noexcept override {
                return 0.0;
            }

        public:
            double cal_h(double x) const noexcept {
                return cal_potential_impl(x, 1) - cal_potential_impl(x, 0);
                return 0.5 * mw2 * hshift * (hshift - 2 * x) + vshift;
            }

            double cal_dhdx(double x) const noexcept {
                return cal_force_impl(x, 0) - cal_force_impl(x, 1);
                static const double dhdx(-mw2 * hshift);
                return dhdx;
            }

        private:
            const double mass, omega, hshift, vshift;
            const double mw2;
            const double gamma;
    };

};
#endif // _DOUBLE_HARMONIC_POTENTIAL_CONST_GAMMA_1D_HPP

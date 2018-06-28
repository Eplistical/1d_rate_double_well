#ifndef _POTENTIAL_BASE_1D_HPP
#define _POTENTIAL_BASE_1D_HPP
// interface module for 1D potential

namespace {

    struct Potential_1D
    {
        public:
            Potential_1D() = default;
            virtual ~Potential_1D() = default;

        public:
            double cal_potential(double x, int i) const noexcept {
                return cal_potential_impl(x, i); 
            }

            double cal_force(double x, int i) const noexcept { 
                return cal_force_impl(x, i); 
            }

            double cal_gamma(double x) const noexcept {
                return cal_gamma_impl(x); 
            }

            double cal_nabla_gamma(double x) const noexcept { 
                return cal_nabla_gamma_impl(x); 
            }

        private:
            virtual double cal_potential_impl(double x, int i) const noexcept = 0;
            virtual double cal_force_impl(double x, int i) const noexcept = 0;
            virtual double cal_gamma_impl(double x) const noexcept = 0;
            virtual double cal_nabla_gamma_impl(double x) const noexcept = 0;
    };

};

#endif // _POTENTIAL_1D_BASE_HPP

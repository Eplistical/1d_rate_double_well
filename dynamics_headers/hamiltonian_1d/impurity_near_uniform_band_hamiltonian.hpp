#ifndef _IMPURITY_NEAR_UNIFORM_BAND_HAMILTONIAN_1D_HPP
#define _IMPURITY_NEAR_UNIFORM_BAND_HAMILTONIAN_1D_HPP

#include <vector>
#include "hamiltonian_base_1d.hpp"
#include "matrixop.hpp"

namespace 
{
	using std::vector;

	template <typename PotentialType> 
		struct ImpurityNearUniformBandHamiltonian_1D final : public Hamiltonian_1D
	{
		public:
			using potential_t = PotentialType;

		public:
			ImpurityNearUniformBandHamiltonian_1D(const potential_t& POTENTIAL, double BANDWIDTH, int NTOTSTATE, int NBANDSTATE, double GAMMA) noexcept :
				Hamiltonian_1D(NTOTSTATE), potential(POTENTIAL)
				{
					band_discretize(BANDWIDTH, NBANDSTATE, GAMMA);
				}

			~ImpurityNearUniformBandHamiltonian_1D() = default;

		// potential functions
		public:
			double cal_potential(double x, int i) const {
				return potential.cal_potential(x, i);
			}

			double cal_force(double x, int i) const {
				return potential.cal_force(x, i);
			}

			double cal_gamma(double x) const {
				return potential.cal_gamma(x);
			}

			double cal_nabla_gamma(double x) const {
				return potential.cal_nabla_gamma(x);
			}

		private:
			void cal_H(double x) override {
				this->H[0] = potential.cal_potential(x, 1) - potential.cal_potential(x, 0);

				for (int k(1); k < this->dim; ++k) {
					this->H[k + k * this->dim] = band_e.at(k - 1);
					this->H[0 + k * this->dim] = band_V;
					this->H[k + 0 * this->dim] = band_V;
				}
			}

			void cal_nabla_H(double x) override {
				this->nabla_H.resize(this->dim * this->dim, 0.0);
				this->nabla_H[0] = potential.cal_force(x, 0) - potential.cal_force(x, 1);
			}

		private:
			void band_discretize(double bandwidth, int Nbandstate, double gamma) {
				// band_V
				const double dos(static_cast<double>(Nbandstate) / bandwidth);
				band_V = sqrt(gamma * 0.5 / M_PI / dos);

				// band_e
				const double de(bandwidth / Nbandstate);
				const double spread(bandwidth * 0.5);
				band_e.clear();
				for (double ei(-spread + de / 2); ei < spread; ei += de) {
					band_e.push_back(ei);
				}
			}

		public:
			const potential_t& potential;
			vector<double> band_e;
			double band_V;
	};

};

#endif // _IMPURITY_NEAR_UNIFORM_BAND_HAMILTONIAN_1D_HPP

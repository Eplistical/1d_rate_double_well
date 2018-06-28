#ifndef _HAMILTONIAN_BASE_1D_HPP
#define _HAMILTONIAN_BASE_1D_HPP
// base module for 1D hamiltonian, abstract class

#include <vector>
#include "matrixop.hpp"

namespace 
{

	struct Hamiltonian_1D
	{
		public:
			explicit Hamiltonian_1D(int DIM) noexcept :
				dim(DIM)
				{
				}

			virtual ~Hamiltonian_1D() noexcept = default;

		public:
			/**
			 * update H, eva, evt
			 */
			void update_H(double x) {
				update_H_impl(x);
			}

			/**
			 * update H, eva
			 */
			void update_H_evaonly(double x) {
				update_H_evaonly(x);
			}

			/**
			 * update derivative coupling, force
			 */
			void update_dc(double x) {
				update_dc_impl(x);
			}

		private:
			virtual void update_H_impl(double x) {
				H.resize(dim * dim);
				eva.resize(dim);
				evt.resize(dim * dim);

				vector<double> last_evt(std::move(evt));

				cal_H(x);
				matrixop::hdiag(H, eva, evt);

				// maintain the wave-function to make sure it changes slowly for each update,
				// this assumes new x is very close to old x
				if (not last_evt.empty()) {
					for (int k(0); k < dim; ++k) {
						if (matrixop::innerproduct(last_evt, evt, k * dim, dim) < 0.0) {
							for (int i(0); i < dim; ++i) {
								evt[i + k * dim] *= -1;
							}
						}
					}
				}
			}

			virtual void update_H_evaonly_impl(double x) {
				H.resize(dim * dim);
				eva.resize(dim);

				cal_H(x);
				matrixop::hdiag(H, eva);
			}

			virtual void update_dc_impl(double x) {
				// calculate nabla_H
				nabla_H.resize(dim * dim);
				cal_nabla_H(x);

            	// calculate F & dc
            	dc = matrixop::matCmat(evt, matrixop::matmat(nabla_H, evt, dim), dim);
				F.resize(dim);
				for (int k(0); k < dim; ++k) {
					F[k] = -dc[k + k * dim];
					dc[k + k * dim] = 0.0;
					for (int j(0); j < k; ++j) {
						dc[k + j * dim] /= (eva[j] - eva[k]);
						dc[j + k * dim] = -dc[k + j * dim];
					}
				}
			}

		private:
			virtual void cal_H(double x) = 0;
			virtual void cal_nabla_H(double x) = 0;

		public:
			const int dim;
			vector<double> H;
			vector<double> eva;
			vector<double> evt;
			vector<double> nabla_H;
			vector<double> dc;
			vector<double> F;
	}; // struct Hamiltonian_1D

}; // namespace

#endif // _HAMILTONIAN_BASE_1D_HPP

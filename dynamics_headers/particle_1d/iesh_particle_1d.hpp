#ifndef _IESH_PARTICLE_HPP
#define _IESH_PARTICLE_HPP

#include <complex>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <random>
#include "misc/vector.hpp"
#include "misc/randomer.hpp"
#include "misc/matrixop_mkl.hpp"
#include "misc/sgn.hpp"
#include "particle_1d_base.hpp"

namespace {
	using std::complex;
	using std::vector;
	using std::mt19937;
	using std::normal_distribution;
	using std::sqrt;
	using std::replace;
    using misc::sgn;

	template <typename HamiltonianType>
		struct IESH_Particle_1D final : public Particle_1D
		{
			public:
				using hamiltonian_t = HamiltonianType;

			public:
				IESH_Particle_1D(  double X, double V, double MASS, double KT, 
						double NUCLEAR_FRIC, int NELE, int NHOLE,
						double NDTQ, double THERMAL_TAU,
						const hamiltonian_t& HAMILTONIAN) noexcept :
                	Particle_1D(X, V, MASS, KT),
					nuclear_fric(NUCLEAR_FRIC),
					Nele(NELE), Nhole(NHOLE),
					Ndtq(NDTQ), 
                    thermal_tau(THERMAL_TAU), thermal_tau_inv(1.0 / THERMAL_TAU),
					hamiltonian(HAMILTONIAN)
					{
						Norb = Nele + Nhole;

						// init occ & uocc
						occ_vec.clear();
						occ_vec.reserve(Nele);
						uocc_vec.clear();
						uocc_vec.reserve(Nhole);
						for (int i(0); i < Norb; ++i) {
							(i < Nele) ? occ_vec.push_back(i) : uocc_vec.push_back(i);
						}

						// init hamiltonian
						hamiltonian.update_H(x);
						hamiltonian.update_dc(x);

						// init psi
						psi.assign(Norb * Nele, matrixop::ZERO_z);
						for (int i(0); i < Nele; ++i) {
							psi[occ_vec[i] + i * Norb].real(1.0);
						}

                        // hop statistics
                        hop_count = 0;
					}

				~IESH_Particle_1D() noexcept = default;

			public:
				double cal_E(const vector<int>& orb) const
				{
					double rst(0.0);
					for (const auto& i : orb)
						rst += hamiltonian.eva.at(i);
					return rst;
				}

				double cal_force(double x) const
				{
					double rst(hamiltonian.cal_force(x, 0));
					for (const auto& i : occ_vec) 
						rst += hamiltonian.F.at(i);
					return rst;
				}

				double get_N(int idiab = 0) const noexcept
				{
					complex<double> r_psi;
					double rst(0.0);

					for (int k(0); k < Nele; ++k) {
						r_psi = complex<double>(0.0, 0.0);
						for (int j(0); j < Norb; ++j) {
							r_psi += hamiltonian.evt[idiab + j * Norb] * psi[j + k * Norb];
						}
						rst += norm(r_psi);
					}
					return rst;
				}

                double get_nu_Ep() const 
                {
                    return hamiltonian.cal_potential(this->x, 0);
                }

                double get_el_Ep() const
                {
                    return cal_E(occ_vec);
                }

                double get_Ep() const
                {
                    return cal_E(occ_vec) 
                        + hamiltonian.cal_potential(this->x, 0);
                }

			private:
				void do_evolve(double dt, const std::string& alg) override
				{
					// nuclear part, velocity verlet
					const double fric(nuclear_fric);
					const double noise_force(randomer::normal(0.0, sqrt(2.0 * fric * this->kT / dt)));
                    assert(fric >= 0.0);

					this->v += (cal_force(this->x) - fric * this->v + noise_force) * 0.5 * dt * this->mass_inv;
					this->x += this->v * dt;
					this->v += (cal_force(this->x) - fric * this->v + noise_force) * 0.5 * dt * this->mass_inv;

					// update hamiltonian
					vector<double> last_eva(hamiltonian.eva);
					hamiltonian.update_H(this->x);
					hamiltonian.update_dc(this->x);

					// electronic part
					bool hopped(false);
					const double dtq(dt / Ndtq);
					const vector<double> deva((hamiltonian.eva - last_eva) / Ndtq);
					vector< complex<double> > T(complex<double>(-this->v, 0.0) * hamiltonian.dc);
					vector< complex<double> > k1, k2, k3, k4;

					for (int idt(0); idt < Ndtq; ++idt) {
                        // evolve
						for (int k(0); k < Norb; ++k) {
							T[k + k * Norb] = -matrixop::I_z * (last_eva[k] + idt * deva[k]);
						}
						k1 = matrixop::matmat(T, psi, Norb, dtq);
						k2 = matrixop::matmat(T, psi + 0.5 * k1, Norb, dtq);
						k3 = matrixop::matmat(T, psi + 0.5 * k2, Norb, dtq);
						k4 = matrixop::matmat(T, psi + k3, Norb, dtq);
						psi = psi + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
						// hopping
						if (not hopped) {
							//hopped = hopper(dtq); 
							hopper(dtq); 
						}
						// el_thermal
						if (thermal_tau > 0) {
							el_thermal(dtq);
						}
					}
				}

				bool hopper(double dt)
				{
					const static int Nhop(Nele * Nhole);
					static vector< complex<double> > S_inv(Nele * Nele);
					static vector<double> probability(Nhop + 1);
					static vector<double> from_to_indices((Nhop + 1) * 2);
					static vector< complex<double> > tmp(Nele);

					// all occ rows, construct S maxtrix and get inverse
					for (int k(0); k < Nele; ++k) {
						matrixop::vcopy(Nele, &psi[occ_vec[k]], Norb, &S_inv[k], Nele);
					}
					matrixop::inv(S_inv);

					// calculate g_kj for all possible jvec
					double Phop(0.0);
					int count(0);
					const complex<double> neg_2_dt_v(-2.0 * dt * this->v, 0.0);

					for (int to_idx(0); to_idx < Nhole; ++to_idx) {
						const int a(uocc_vec[to_idx]);
						zgemv("T", &Nele, &Nele, &neg_2_dt_v, &S_inv[0], &Nele,
								&psi[a], &Norb, &matrixop::ZERO_z, 
								&tmp[0], &matrixop::ONE_i);
						for (int from_idx(0); from_idx < Nele; ++from_idx) {
							const int i(occ_vec[from_idx]);
							if (sgn(tmp[from_idx].real()) == sgn(hamiltonian.dc[a + i * Norb])) {
								probability[count] = tmp[from_idx].real() * hamiltonian.dc[a + i * Norb];
								Phop += probability[count];
								from_to_indices[count * 2] = from_idx;
								from_to_indices[1 + count * 2] = to_idx;
								++count;
							}
						}
					}

                    // if Phop is too large 
                    int Nattempt(1);
                    if (Phop > 1.0) {
                        Nattempt = static_cast<int>(Phop / 0.1) + 1;
                    }
                    probability = probability / Nattempt;
                    Phop = Phop / Nattempt;
                    
					// append non-hop probability
					probability[count] = 1.0 - Phop; 
					from_to_indices[count * 2] = -1;
					from_to_indices[1 + count * 2] = -1;

					// random number
                    for (int iattempt(0); iattempt < Nattempt; ++iattempt) {
                        const int attempt_hop(randomer::discrete(probability.begin(), probability.begin() + count + 1));
                        if (attempt_hop != count) {
                            const int attempt_from_idx(from_to_indices[attempt_hop * 2]);
                            const int attempt_to_idx(from_to_indices[1 + attempt_hop * 2]);
                            const int i(occ_vec[attempt_from_idx]);
                            const int a(uocc_vec[attempt_to_idx]);
                            const double Ek(this->get_Ek());
                            const double deltaE(hamiltonian.eva[a] - hamiltonian.eva[i]);
                            // not frustrated
                            if (deltaE <= Ek) {
                                this->v *= sqrt((Ek - deltaE) / Ek);
                                occ_vec[attempt_from_idx] = a;
                                uocc_vec[attempt_to_idx] = i;
                                hop_count += 1;
                                return true;
                            }
                        }
                    }
                    return false;
				}

                void el_thermal(double dt)
                {
                    if (randomer::rand() < dt * thermal_tau_inv) {
                        const int i(randomer::choice(occ_vec));
                        const int a(randomer::choice(uocc_vec));
                        const double Ea(hamiltonian.eva.at(a));
                        const double Ei(hamiltonian.eva.at(i));

                        if ((Ei >= Ea) or randomer::rand() < exp(-this->kT_inv * (Ea - Ei))) {
                            replace(occ_vec.begin(), occ_vec.end(), i, a);
                            replace(uocc_vec.begin(), uocc_vec.end(), a, i);
                        }
                    }
                }

			public:
				int Norb, Nele, Nhole;
				vector< complex<double> > psi;
				vector<int> occ_vec, uocc_vec;
				hamiltonian_t hamiltonian;
				double nuclear_fric;
				int Ndtq;
                const double thermal_tau, thermal_tau_inv;
                uint64_t hop_count;
		};

};

#endif // _IESH_PARTICLE_HPP

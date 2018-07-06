#ifndef _INTTABLE_MGR_HPP
#define _INTTABLE_MGR_HPP

#include <map>
#include <string>
#include <vector>
#include <string>
#include <cmath>
#include "inttable/inttable_mgr_base.hpp"
#include "misc/ioer.hpp"

namespace {
    using std::string;
    using std::map;
    using std::vector;

    template <typename PotentialType>
        struct InttableMgr final : InttableMgr_Base {
            public:
                using potential_t = PotentialType;

            public:
                InttableMgr(const string& FNAME, const potential_t& POTENTIAL) noexcept :
                    InttableMgr_Base(FNAME),
                    potential(POTENTIAL)
                    {
                    }

                ~InttableMgr() noexcept = default;

            private:
                void load_inttable_impl() override {
                    ioer::h5file_t h5f(this->fname, std::ios::in);
                    // read header
                    h5f.read_attr("para",
                            "kT", kT,
                            "gamma0", gamma0,
                            "hmin", this->xmin,
                            "hmax", this->xmax,
                            "dh", this->dx,
                            "Nh", this->Nx
                            );
                    this->dx_inv = 1.0 / this->dx;
                    // read data
                    this->inttable_dict.insert(make_pair("Af", vector<double>(this->Nx, 0.0)));
                    this->inttable_dict.insert(make_pair("A2dfde", vector<double>(this->Nx, 0.0)));
                    h5f.read_dataset(
                            "Af", this->inttable_dict.at("Af"),
                            "A2dfde", this->inttable_dict.at("A2dfde")
                            );
                    // close file
                    h5f.close();
                }

                double retrieve_impl(const string& key, double x) const override {
                    const static double pi_inv(1.0 / M_PI);
                    double h(potential.cal_h(x));
                    if (key == "n") {
                        return 0.5 * pi_inv * this->retrieve_raw("Af", h);
                    }
                    else if (key == "force") {
                        double dhdx(potential.cal_dhdx(x));
                        return -dhdx * 0.5 * pi_inv * this->retrieve_raw("Af", h);
                    }
                    else if (key == "fric") {
                        double dhdx(potential.cal_dhdx(x));
                        return -0.25 * pi_inv * pow(dhdx, 2) * this->retrieve_raw("A2dfde", h);
                    }
                    else {
                        return this->retrieve_raw(key, h);
                    }
                }

            public:
                const potential_t& potential;
                double kT, gamma0;
        };
};

#endif // _INTTABLE_MGR_HPP

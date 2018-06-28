#ifndef _PARA_HPP
#define _PARA_HPP

#include <cstdlib>
#include <string>
#include "ioer.hpp"

namespace
{
    using std::string;
    using ioer::keyval;
    using ioer::input_t;
    using ioer::output_t;
    using ioer::STDOUT;

    struct Para
    {

        int M;
        int Ne;
        double W;

        Para() = default;

        Para(const string& fname) {
            loadpara(fname);
        }

        void loadpara(const string& fname) {
            input_t inp(fname);

            inp.extract_para("M", M);
            inp.extract_para("Ne", Ne);
            inp.extract_para("W", W);

            inp.read_text();
            inp.close();
        }

        void showpara(output_t& outp=STDOUT)
        {
            outp.tabout("#", " M", M);
            outp.tabout("#", " Ne", Ne);
            outp.tabout("#", " W", W);
        }
    };

#endif // _PARA_HPP

#!/usr/bin/env python3
import sys
import re

hppfile ='''
#ifndef _PARA_HPP
#define _PARA_HPP

#include <cstdlib>
#include <string>
#include "misc/ioer.hpp"

namespace
{
    using std::string;
    using ioer::keyval;
    using ioer::input_t;
    using ioer::output_t;
    using ioer::STDOUT;

    struct Para
    {

%s
        Para() = default;

        Para(const string& fname) {
            loadpara(fname);
        }

        void loadpara(const string& fname) {
            input_t inp(fname);
            inp.read_text();

%s

            inp.close();
        }

        void showpara(output_t& outp=STDOUT)
        {
%s
        }
    };

};

#endif // _PARA_HPP
'''

def extract_varlist(fname):
    reg = re.compile('struct Para.*{(.*?)}', re.DOTALL)
    with open(fname, 'r') as f:
        paras = re.search(reg, f.read()).group(1)
    reg2 = re.compile('.* ([^ ]*);')
    varlist = []
    varpart = ""
    for line in paras.split('\n'):
        if line.strip():
            v = re.search(reg2, line).group(1)
            varlist.append(v)
            varpart += "        " + line.strip() + '\n'
    return varlist, varpart


def get_load_para_content(varlist):
    rst = ''''''
    tmp = '''
            inp.extract_para("%s", %s);
'''
    for v in varlist:
        rst += tmp[1:] % (v, v,)
    return rst[:-1]


def get_show_para_content(varlist):
    rst = ''''''
    tmp = '''
            outp.tabout("#", " %s", %s);
'''
    for v in varlist:
        rst += tmp[1:] % (v, v,)
    return rst[:-1]

def get_example_input(varlist):
    rst = ''''''
    tmp = '''
0        # %s
'''
    for v in varlist:
        rst += tmp[1:] % v
    return rst[:-1]

if __name__ == '__main__':
    fname = sys.argv[1]
    varlist, varpart = extract_varlist(fname)
    # write para.hpp
    with open('para.hpp', 'w') as f:
        f.write(hppfile[1:] % (varpart,
                get_load_para_content(varlist),
                get_show_para_content(varlist))
                )

    # generate example input file
    with open('sample.in', 'w') as f:
        f.write(get_example_input(varlist))

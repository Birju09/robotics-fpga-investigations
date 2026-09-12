#ifndef TB_COMMON_HPP
#define TB_COMMON_HPP

//! Same sources build against ap_fixed (Vitis csim/cosim) and against the
//! float reference (host g++, hls/tb/Makefile), so regressions run without
//! a Vitis install.
//
//! Vector files (model/gen_vectors.py) hold signed decimal Q16.16 words,
//! matching what the PS driver writes into the IP.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

#include "ik_config.hpp"
#include "ik_types.hpp"

//! HLS C-sim runs from a generated dir several levels down whose depth
//! varies by tool version, so try plausible roots in order.
//! $IK_VECTOR_DIR overrides everything.
static inline std::ifstream tb_open(const char* name) {
    static const char* roots[] = {
#ifdef TB_VECTOR_DIR
        TB_VECTOR_DIR,
#endif
        "vectors",
        "../vectors",
        "../tb/vectors",
        "../../tb/vectors",
        "../../../tb/vectors",
        "../../../../tb/vectors",
        "../../../../../tb/vectors",
        "../../../../../../tb/vectors",
    };

    std::vector<std::string> tried;

    const char* env = std::getenv("IK_VECTOR_DIR");
    if (env) {
        std::string p = std::string(env) + "/" + name;
        std::ifstream f(p.c_str());
        if (f.is_open())
            return f;
        tried.push_back(p);
    }

    for (size_t i = 0; i < sizeof(roots) / sizeof(roots[0]); i++) {
        std::string p = std::string(roots[i]) + "/" + name;
        std::ifstream f(p.c_str());
        if (f.is_open())
            return f;
        tried.push_back(p);
    }

    std::fprintf(stderr, "FATAL: cannot locate vector file '%s'. Tried:\n",
                 name);
    for (size_t i = 0; i < tried.size(); i++)
        std::fprintf(stderr, "         %s\n", tried[i].c_str());
    std::fprintf(stderr,
                 "       Generate them with 'python3 model/gen_vectors.py',\n"
                 "       or point $IK_VECTOR_DIR at the directory.\n");
    std::exit(2);
}

static inline double tb_dbl(ik_word_t w) {
    return (double)w / 65536.0;
}

//! Angular difference wrapped to [-pi, pi]: the kernels return joint angles
//! wrapped, the golden model does not, and they are equivalent mod 2*pi.
static inline double tb_angdiff(double a, double b) {
    double d = a - b;
    while (d > M_PI)
        d -= 2.0 * M_PI;
    while (d < -M_PI)
        d += 2.0 * M_PI;
    return d;
}

struct TbStats {
    const char* name;
    int n;
    int fail;
    double worst;

    TbStats(const char* nm) : name(nm), n(0), fail(0), worst(0.0) {}

    void note(double err, double tol) {
        n++;
        if (err > worst)
            worst = err;
        if (!(err <= tol))
            fail++;
    }

    int report(double tol) const {
        std::printf("  %-28s vectors=%-5d worst=%.3e tol=%.1e  %s\n", name, n,
                    worst, tol, fail ? "FAIL" : "pass");
        if (fail)
            std::printf("      %d/%d vectors exceeded tolerance\n", fail, n);
        return fail ? 1 : 0;
    }
};

#endif  //! TB_COMMON_HPP

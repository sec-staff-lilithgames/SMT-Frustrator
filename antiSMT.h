#pragma once
#include <string>
#include <tuple>

struct QF {
    int coeff_uu;
    int coeff_uv;
    int coeff_vv;
    int rhs;
};

struct PellRes {
    std::string expression;
    QF qf;
    int D;
    std::tuple<int, int, int, int> mat;
    int timeout;
    int noise;

    // operator bool
    explicit operator bool() const { return D >= 0; }
};

namespace PellOP {
    // 推荐C++风格主接口
    PellRes getPellOP(
        int D_low = 1, int D_high = 100,
        int coeff_range = 10,
        int timeout = 60000,
        double noise_lower = -0.02,
        double noise_upper = 0.02
    );

    // 若需要C导出接口，加转换和说明
    // extern "C" bool getPellOP(char pattern, float timeout, int u, int v);
}
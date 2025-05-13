#include "antiSMT.h"
#include <random>
#include <chrono>
#include <z3++.h>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <tuple>
#include <sstream>




// =============================生成器===============================//
// 生成一个D（完全平方数），返回-1说明失败
int generate_PellSeeds(int low, int high) {
    if (low > high || high <= 0) return -1; // 边界检查
    int minSqrt = (int)ceil(sqrt((double)low));
    int maxSqrt = (int)floor(sqrt((double)high));
    if (minSqrt > maxSqrt) return -1;
    int range = maxSqrt - minSqrt + 1;
    int r = rand() % range;
    int sq = (minSqrt + r) * (minSqrt + r);
    return sq;
}


// ============================转换器===============================//
// 过渡矩阵生成器
std::tuple<int,int,int,int> random_invertible_matrix(int range) {
    if (range <= 0) throw std::invalid_argument("Matrix range must be positive.");
    static std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<int> dist(-range, range);

    for(int trial=0; trial<100; ++trial) { // 最多尝试100次，避免死循环
        int a = dist(rng), b = dist(rng), c = dist(rng), d = dist(rng);
        if(a * d - b * c != 0)
            return std::make_tuple(a, b, c, d);
    }
    throw std::runtime_error("Failed to generate invertible matrix.");
}

// 二次型变换器
QF pell_transformer(int D, const std::tuple<int,int,int,int>& mat) {
    int a = std::get<0>(mat);
    int b = std::get<1>(mat);
    int c = std::get<2>(mat);
    int d = std::get<3>(mat);

    QF res;

    // https://li.feishu.cn/docx/YTfgdMT1qoTDhCx2n0NcHTMinQd
    // 这里直接完全展开是因为括号非常难处理
    int delta = a * d - b * c;
    res.coeff_uu = d * d - D * c * c;
    res.coeff_uv = -2 * d * b + 2 * D * c * a;
    res.coeff_vv = b * b - D * a * a;
    res.rhs = delta * delta;
    return res;
}

// ===========================验证器=================================

QF pell_verifer(const QF& qf, int timeout = 1000, int noise = 100){
    using namespace z3;
    if(timeout <= 0) throw std::invalid_argument("timeout must be positive");
    if(noise < 0) throw std::invalid_argument("noise must be non-negative");

    context ctx;
    solver s(ctx);

    expr u = ctx.int_const("u");
    expr v = ctx.int_const("v");
    expr eq = qf.coeff_uu * u * u + qf.coeff_uv * u * v + qf.coeff_vv * v * v == qf.rhs;

    s.add(eq);

    params p(ctx);
    p.set("timeout", static_cast<unsigned>(timeout));
    s.set(p);

    auto t0 = std::chrono::steady_clock::now();
    check_result result;
    try {
        result = s.check();
    } catch (const z3::exception& e) {
        // Z3库崩溃也不能炸进程
        return qf;
    }
    auto t1 = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double, std::milli>(t1 - t0).count();

    if(result == unknown ||
      (result == unsat && (elapsed > timeout - noise && elapsed < timeout + noise)))
    {
        return qf;
    }
    return qf;
}

// ==================================格式转换器===========================
// 变量表达式格式化
std::string QF_to_expr(const QF& qf, const std::string& u = "u", const std::string& v = "v") {
    std::ostringstream oss;
    bool first = true;

    if (qf.coeff_uu != 0) {
        if (qf.coeff_uu < 0) oss << "-";
        if (std::abs(qf.coeff_uu) != 1) oss << std::abs(qf.coeff_uu);
        oss << u << "^2";
        first = false;
    }
    if (qf.coeff_uv != 0) {
        if (!first) {
            if (qf.coeff_uv > 0) oss << " + ";
            else oss << " - ";
        } else {
            if (qf.coeff_uv < 0) oss << "-";
        }
        if (std::abs(qf.coeff_uv) != 1) oss << std::abs(qf.coeff_uv);
        oss << u << v;
        first = false;
    }
    if (qf.coeff_vv != 0) {
        if (!first) {
            if (qf.coeff_vv > 0) oss << " + ";
            else oss << " - ";
        } else {
            if (qf.coeff_vv < 0) oss << "-";
        }
        if (std::abs(qf.coeff_vv) != 1) oss << std::abs(qf.coeff_vv);
        oss << v << "^2";
        first = false;
    }
    if (first) { // 全为0
        oss << "0";
    }
    oss << " = " << qf.rhs;
    return oss.str();
}



namespace PellOP {

PellRes getPellOP(
    int D_low, int D_high,
    int coeff_range,
    int timeout,
    double noise_lower,
    double noise_upper
)
{
    PellRes ret;
    try {
        // ==== 检查参数 =====
        if(D_low > D_high || D_high <= 0)
            throw std::invalid_argument("Invalid D range.");
        if(coeff_range <= 0)
            throw std::invalid_argument("Coefficient range must be positive.");
        if(timeout <= 0)
            throw std::invalid_argument("Timeout must be positive.");
        if(noise_lower >= noise_upper)
            throw std::invalid_argument("Noise range invalid: lower must be < upper.");

        // ==== 步骤 1 - 2 ====
        int D = generate_PellSeeds(D_low, D_high);
        if(D == -1)
            throw std::runtime_error("Unable to generate valid D.");

        ret.D = D;
        ret.mat = random_invertible_matrix(coeff_range);

        // ==== 步骤 3 ====
        QF qf = pell_transformer(D, ret.mat);

        // ==== 步骤 4 ====
        static std::mt19937 rng(std::random_device{}());
        std::uniform_real_distribution<> dist(noise_lower, noise_upper);
        int noise = static_cast<int>(timeout * dist(rng));
        if(noise < 0) noise = -noise; // 确保噪声为正
        ret.timeout = timeout;
        ret.noise = noise;

        // ==== 步骤 5 ====
        QF checked_qf = pell_verifer(qf, timeout, noise);
        ret.qf = checked_qf;

        // ==== 步骤 6 ====
        std::ostringstream oss;
        oss << checked_qf.coeff_uu << "*u^2";
        if(checked_qf.coeff_uv != 0) {
            if(checked_qf.coeff_uv > 0) oss << " + ";
            else oss << " - ";
            oss << std::abs(checked_qf.coeff_uv) << "*uv";
        }
        if(checked_qf.coeff_vv != 0) {
            if(checked_qf.coeff_vv > 0) oss << " + ";
            else oss << " - ";
            oss << std::abs(checked_qf.coeff_vv) << "*v^2";
        }
        oss << " = " << checked_qf.rhs;
        ret.expression = oss.str();

    } catch (const std::exception& e) {
        // 错误填空表达式，其余置默认
        ret.expression = std::string("Error: ") + e.what();
        ret.D = -1; ret.timeout = ret.noise = 0;
        ret.mat = std::make_tuple(0,0,0,0);
        ret.qf = QF{0,0,0,0};
    }
    return ret;
}

} // namespace PellOP



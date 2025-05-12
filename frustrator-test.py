#########################################################################
# 这是一个非常简单的对抗**符号执行**的约束生成器
# 利用的是目前常见SMT求解器在对非线性约束支持上的缺陷
# 我们采用 Pell 方程作为约束种子，经过一组可逆的线性变换变成同构不同形的约束形式
# 这个约束生成器会生成任意**无解**的Pell-like约束，注意一定是整型变量
# 使用的时候只需要把结果替换为工程项目里面的逻辑False即可


import random
from z3 import *
import time

# ----------- 生成器 -----------
# 注意，这个生成器输出的约束一定是无解约束
# 因为经过测试SMT求解存在性问题和任意性问题上表现差异巨大

# 两种模式，只要是完全平方数情况都是无解的，效果恒定
def pell_generator(mode='simple'):
    constraints_list = []
    if mode == 'simple':
        # 默认仅取D=1
        D_list = [1]
    elif mode == 'complete':
        # 取D为[1^2, 2^2, ..., 7^2]（最大不超过62），可根据需求调整上界
        # 若需要更大，可以修改range(1, N)
        D_list = [k**2 for k in range(1, 8)]  # D=[1,4,9,16,25,36,49]
    else:
        raise RuntimeError('Unknown generator mode. Use "simple" or "complete"')
    for D in D_list:
        x = Int('x')
        y = Int('y')
        x_domain = x > 0
        y_domain = y > 0
        pell_constraint = x*x - D*y*y == 1
        constraints = [x_domain, y_domain, pell_constraint]
        constraints_list.append({'D': D, 'constraints': constraints, 'vars': (x, y),
                                'description': f"x > 0, y > 0, x^2 - {D}*y^2 == 1"})
    return constraints_list

# ----------- 变换器 -----------


# 这里可以控制过度矩阵的扰动范围
def random_invertible_matrix(entry_range=(-10, 10)):
    while True:
        a, b, c, d = [random.randint(*entry_range) for _ in range(4)]
        if a*d - b*c != 0:
            return a, b, c, d


# 作的是可逆的线性变换，因此变换后仍然无解
def pell_transformer(x, y, D, mat=None):
    if mat is None:
        a, b, c, d = random_invertible_matrix()
    else:
        a, b, c, d = mat
    u = Int('u')
    v = Int('v')
    new_constraint = (a*u + b*v) * (a*u + b*v) - D*(c*u + d*v)*(c*u + d*v) == 1
    desc = f"u,v ∈ Z, new constraint: ({a}*u+{b}*v)^2 - {D}*({c}*u+{d}*v)^2 == 1"
    constraints = [new_constraint]
    return {'constraints': constraints, 'vars': (u, v), 'matrix': (a, b, c, d),
            'description': desc}



# ----------- 验证器 -----------
# 只验证是否超时，目前只支持z3作为验证器
def timeout_verifier(constraints, timeout_limit=10.0):
    set_param('verbose', 0)
    
    # 指定z3专门用于求解非线性整型约束的策略，不指定也可以
    t = Tactic('qfnia')
    s = t.solver()
    s.set(timeout=int(timeout_limit * 1000))
    for c in constraints:
        s.add(c)
    start_time = time.time()
    result = s.check()
    elapsed_time = time.time() - start_time
    is_timeout = (result == unknown) and (elapsed_time >= timeout_limit)
    return is_timeout, elapsed_time, result



# ----------- 主体控制 -----------
if __name__ == '__main__':
    #mode = 'simple'  # 'simple' 或 'complete'
    mode = 'complete'
    timeout_limit = 60.0  # 超时门限（秒）
    gens = pell_generator(mode=mode)
    output_count = 0
    for item in gens:
        D = item['D']
        x, y = item['vars']
        trans_out = pell_transformer(x, y, D)
        new_constraints = trans_out['constraints']
        desc_str = trans_out['description']
        # <<< 修正：在此处解包 a, b, c, d >>>
        a, b, c, d = trans_out['matrix']
        is_timeout, elapsed, res = timeout_verifier(new_constraints, timeout_limit=timeout_limit)
        if is_timeout:
            output_count += 1
            print("="*30)
            print(f"合格约束#{output_count}:")
            print(f"Pell seed: D={D}")
            print(f"新变量及定义域: u, v ∈ Z")
            print(f"所用可逆线性变换矩阵: [[{a}, {b}], [{c}, {d}]]")
            print(f"新约束结构: {desc_str}")
            
            # 新增完全展开式输出
            coef_u2 = a**2 - D * c**2
            coef_uv = 2*a*b - 2*D*c*d
            coef_v2 = b**2 - D * d**2
            # 可选美化
            def coefstr(c, mono):
                if c == 1:
                    return f"{mono}"
                elif c == -1:
                    return f"-{mono}"
                else:
                    return f"{c}*{mono}"
            terms = []
            if coef_u2 != 0:
                terms.append(coefstr(coef_u2, "u^2"))
            if coef_uv != 0:
                terms.append(coefstr(coef_uv, "u*v"))
            if coef_v2 != 0:
                terms.append(coefstr(coef_v2, "v^2"))
            lin_comb = " + ".join(terms)
            print(f"完全展开式: {lin_comb} == 1")
            
            print(f"超时用时: {elapsed:.3f}s (门限: {timeout_limit}s)")
    if output_count == 0:
        print(f"\n模式[{mode}]下没有超时（合格）约束。")
    else:
        print(f"\n共输出{output_count}个合格约束。")

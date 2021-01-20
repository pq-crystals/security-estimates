from math import *
from model_BKZ import *
from proba_util import gaussian_center_weight

log_infinity = 9999

STEPS_b = 1
STEPS_m = 5


class MSISParameterSet:
    def __init__(self, n, w, h, B, q, norm=""):
        self.n = n      # Ring Dimension
        self.w = w      # MSIS dimension
        self.h = h      # Number of equations
        self.B = B      # norm bound
        self.q = q      # Modulus
        self.norm = norm


def SIS_l2_cost(q, w, h, B, b, cost_svp=svp_classical, verbose=False):
    """ Return the cost of finding a vector shorter than B with BKZ-b if it works.
    The equation is Ax = 0 mod q, where A has h rows, and w collumns (h equations in dim w).
    """
    if B>=q:
        if verbose:
            print ("Do not know how to handle B>=q in l2 norm. Concluding 0 bits of security.")
        return 0
    l = BKZ_first_length(q, h, w-h, b)
    if l > B:
        return log_infinity 
    if verbose:
        print ("Attack uses block-size %d and %d equations"%(b, h))
        print ("shortest vector used has length l=%.2f, q=%d, `l<q'= %d"%(l, q, l<q))
    return cost_svp(b)


def SIS_linf_cost(q, w, h, B, b, cost_svp=svp_classical, verbose=False):
    """ Return the cost of finding a vector shorter than B in infinity norm, using BKZ-b, if it works.
    The equation is Ax = 0 mod q, where A has h rows, and w columns (h equations in dim w).
    """
    (i, j, L) = construct_BKZ_shape_randomized(q, h, w-h, b)
    #(i, j, L) = construct_BKZ_shape(h * log(q), 1, w-1, b)

    
    l = exp(L[i])
    d = j - i + 1
    sigma = l / sqrt(j - i + 1)
    p_middle = gaussian_center_weight(sigma, B)
    p_head = 2.*B / q

    log2_eps = d * log(p_middle, 2) + i * log(p_head, 2)
    log2_R = max(0, - log2_eps - nvec_sieve(b))

    if verbose:
        print ("Attack uses block-size %d and %d dimensions, with %d q-vectors"%(b, w, i))
        print ("log2(epsilon) = %.2f, log2 nvector per run %.2f"%(log2_eps, nvec_sieve(b)))
        print ("shortest vector used has length l=%.2f, q=%d, `l<q'= %d"%(l, q, l<q))
    return cost_svp(b) + log2_R


def SIS_optimize_attack(q, max_w, h, B, cost_attack=SIS_linf_cost, cost_svp=svp_classical, verbose=False):
    """ Find optimal parameters for a given attack
    """
    best_cost = log_infinity

    for b in range(50, max_w, STEPS_b):
        if cost_svp(b) > best_cost:
            break
        for w in [max_w]:    # No need to exhaust w here as the attack will auto-adjust anyway  range(max(h+1, b+1), max_w, STEPS_m):
            cost = cost_attack(q, w, h, B, b, cost_svp)
            if cost<=best_cost:
                best_cost = cost
                best_w = w
                best_b = b

    if verbose:
        cost_attack(q, best_w, h, B, best_b, cost_svp=cost_svp, verbose=verbose)

    return (best_w, best_b, best_cost)


def check_eq(m_pc, m_pq, m_pp):
    if (m_pc != m_pq):
        print("m and b not equals among the three models")
    if (m_pq != m_pp):
        print("m and b not equals among the three models")


def MSIS_summarize_attacks(ps):
    """ Create a report on the best primal and dual BKZ attacks on an l_oo - MSIS instance
    """
    q = ps.q
    h = ps.n * ps.h
    max_w = ps.n * ps.w
    B = ps.B

    if ps.norm == "linf":
        attack = SIS_linf_cost
    elif ps.norm == "l2":
        attack = SIS_l2_cost
    else:
        raise ValueError("Unknown norm: "+ps.norm)

    (m_pc, b_pc, c_pc) = SIS_optimize_attack(q, max_w, h, B, cost_attack=attack, cost_svp=svp_classical, verbose=True)
    (m_pq, b_pq, c_pq) = SIS_optimize_attack(q, max_w, h, B, cost_attack=attack, cost_svp=svp_quantum, verbose=False)
    (m_pp, b_pp, c_pp) = SIS_optimize_attack(q, max_w, h, B, cost_attack=attack, cost_svp=svp_plausible, verbose=False)

    check_eq(m_pc, m_pq, m_pp)
    check_eq(b_pc, b_pq, b_pp)

    print("SIS & %d & %d & %d & %d & %d"%(m_pq, b_pq, int(floor(c_pc)), int(floor(c_pq)), int(floor(c_pp))))

    return (b_pq, int(floor(c_pc)), int(floor(c_pq)), int(floor(c_pp)))

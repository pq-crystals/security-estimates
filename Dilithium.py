from MSIS_security import MSIS_summarize_attacks, MSISParameterSet
from MLWE_security import MLWE_summarize_attacks, MLWEParameterSet
from math import sqrt


class UniformDilithiumParameterSet(object):
    def __init__(self, n, k, l, gamma1, gamma2, tau, q, eta, pkdrop=0):
        self.n = n
        self.k = k
        self.l = l
        self.gamma1 = gamma1
        self.gamma2 = gamma2
        self.q = q
        self.eta = eta
        # SIS l_oo bound for unforgeability
        self.zeta = max(gamma1, 2*gamma2 + 1 + 2**(pkdrop-1)*tau)                # Define in section 6.2
        # SIS l_oo bound for strong unforgeability
        self.zeta_prime = max(2*gamma1, 4*gamma2 + 1)                # Define in section 6.2
        self.pkdrop = pkdrop


# class GaussianDilithiumParameterSet(object):
#     def __init__(self, n, k, l, sigma, q, eta, pkdrop=0):
#         self.n = n
#         self.k = k
#         self.l = l
#         self.sigma = sigma
#         self.q = q
#         self.eta = eta
#         self.pkdrop = pkdrop
#         self.B = 2*equation5(self)


# def equation5(dps):
#     B2 = ((1.05 * dps.sigma * sqrt((dps.k + dps.l)*dps.n))**2 
#           +(2**(dps.pkdrop-1) * sqrt(60*dps.n*dps.k))**2)
#     return sqrt(B2)

n = 256
q = 8380417

# UnifWeakDilithium           = UniformDilithiumParameterSet(n, 3, 2, gamma, q, 7, pkdrop=14)

UnifMediumDilithium         = UniformDilithiumParameterSet(n, 4, 4, 2**17, (q-1)/88, 39, q, 2, pkdrop=13)
UnifRecommendedDilithium    = UniformDilithiumParameterSet(n, 6, 5, 2**19, (q-1)/32, 49, q, 4, pkdrop=13)
UnifVeryHighDilithium       = UniformDilithiumParameterSet(n, 8, 7, 2**19, (q-1)/32, 60, q, 2, pkdrop=13)


all_params_unif = [#("Uniform Dilithium Weak", UnifWeakDilithium),
                   ("Uniform Dilithium Medium", UnifMediumDilithium),
                   ("Uniform Dilithium Recommended", UnifRecommendedDilithium),
                   ("Uniform Dilithium Very High", UnifVeryHighDilithium)]

all_params = all_params_unif


def Dilithium_to_MSIS(dps, strong_uf = False):
    if strong_uf:
        return MSISParameterSet(dps.n, dps.k + dps.l + 1, dps.k, dps.zeta_prime, dps.q, norm="linf")
    else:
        return MSISParameterSet(dps.n, dps.k + dps.l + 1, dps.k, dps.zeta, dps.q, norm="linf")


def Dilithium_to_MLWE(dps):
    return MLWEParameterSet(dps.n, dps.l, dps.k, dps.eta, dps.q, distr="uniform")

text_SIS = ["BKZ block-size $b$ to break SIS","Best Known Classical bit-cost","Best Known Quantum bit-cost","Best Plausible bit-cost"]
text_LWE = ["BKZ block-size $b$ to break LWE","Best Known Classical bit-cost","Best Known Quantum bit-cost","Best Plausible bit-cost"]



table_SIS = [4*[0] for i in range(4)]
table_LWE = [4*[0] for i in range(4)]
j = 0

for (scheme, param) in all_params_unif:
    print("\n"+scheme)
    print(param.__dict__)
    print("")
    print("=== WEAK UF")
    v = MSIS_summarize_attacks(Dilithium_to_MSIS(param))
    print("=== STRONG UF")
    v = MSIS_summarize_attacks(Dilithium_to_MSIS(param, strong_uf=True))
    for i in range(4):
        table_SIS[i][j] = v[i]
    print("=== SECRET KEY RECOVERY")
    v = MLWE_summarize_attacks(Dilithium_to_MLWE(param))
    for i in range(4):
        table_LWE[i][j] = v[i]
    j+=1

print("UNIFORM DILITHIUM TABLE")
print("========================")

print("\\hline")
for j in range(4):
    print(text_SIS[j]+" & "),
    for i in range(4):
        print(table_SIS[j][i]),
        if i<3:
            print(" & "),
    print("\\\\")
print("\\hline")
for j in range(4):
    print(text_LWE[j]+" & "),
    for i in range(4):
        print(table_LWE[j][i]),
        if i<3:
            print(" & "),
    print("\\\\")
print("\\hline")

print("========================")

table_SIS = [4*[0] for i in range(4)]
table_LWE = [4*[0] for i in range(4)]
j = 0

from MSIS_security import MSIS_summarize_attacks, MSISParameterSet
from MLWE_security import MLWE_summarize_attacks, MLWEParameterSet
from math import sqrt


class UniformDilithiumParameterSet(object):
    def __init__(self, n, k, l, gamma, q, eta, pkdrop=0):
        self.n = n
        self.k = k
        self.l = l
        self.gamma = gamma
        self.q = q
        self.eta = eta
        self.B = max(2*gamma, 2**(pkdrop+1))
        self.pkdrop = pkdrop


class GaussianDilithiumParameterSet(object):
    def __init__(self, n, k, l, sigma, q, eta, pkdrop=0):
        self.n = n
        self.k = k
        self.l = l
        self.sigma = sigma
        self.q = q
        self.eta = eta
        self.pkdrop = pkdrop
        self.B = 2*equation5(self)


def equation5(dps):
    B2 = ((1.05 * dps.sigma * sqrt((dps.k + dps.l)*dps.n))**2 
          +(2**(dps.pkdrop-1) * sqrt(60*dps.n*dps.k))**2)
    return sqrt(B2)

n = 256
q = 8380417
gamma = (q-1)/16

UnifWeakDilithium           = UniformDilithiumParameterSet(n, 3, 2, gamma, q, 7, pkdrop=14)
UnifMediumDilithium         = UniformDilithiumParameterSet(n, 4, 3, gamma, q, 6, pkdrop=14)
UnifRecommendedDilithium    = UniformDilithiumParameterSet(n, 5, 4, gamma, q, 5, pkdrop=14)
UnifVeryHighDilithium       = UniformDilithiumParameterSet(n, 6, 5, gamma, q, 3, pkdrop=14)


all_params_unif = [("Uniform Dilithium Weak", UnifWeakDilithium),
                   ("Uniform Dilithium Medium", UnifMediumDilithium),
                   ("Uniform Dilithium Recommended", UnifRecommendedDilithium),
                   ("Uniform Dilithium Very High", UnifVeryHighDilithium)]

all_params = all_params_unif


def Dilithium_to_MSIS(dps):
    if type(dps)==UniformDilithiumParameterSet:
        return MSISParameterSet(dps.n, dps.k + dps.l + 1, dps.k, dps.B, dps.q, norm="linf")
    if type(dps)==GaussianDilithiumParameterSet:
        return MSISParameterSet(dps.n, dps.k + dps.l + 1, dps.k, dps.B, dps.q, norm="l2")
    else:
        raise ValueError("Unrecognized Dilithium Parameter Type")


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
    v = MSIS_summarize_attacks(Dilithium_to_MSIS(param))
    for i in range(4):
        table_SIS[i][j] = v[i]
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

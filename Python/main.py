import random
import math
import copy
import time


def afficher_poly(P):
    print(f"{P[0]}", end="")
    for i in range(1, len(P)):
        if P[i] >= 0:
            print(f" + {P[i]}x^{i}", end="")
        else:
            print(f" - {-P[i]}x^{i}", end="")
    print()


def gen_poly_ent(n, maxi):
    return [random.randint(0, maxi) for i in range(n + 1)]


def gen_poly_reel(n, maxi):
    return [random.random() * maxi for i in range(n + 1)]


def prod_poly_naif(P, Q):
    res = [0 for _ in range(len(P) + len(Q) - 1)]
    for i in range(len(P)):
        for j in range(len(Q)):
            res[i + j] += P[i] * Q[j]
    return res


def somme_poly(P, Q, deb):
    # on suppose len(P) >= len(Q)
    res = copy.deepcopy(P)
    for i in range(len(Q)):
        res[deb + i] += Q[i]
    return res


def oppose_poly(P):
    return [-P[i] for i in range(len(P))]


def prod_poly_karatsuba(P, Q):
    if len(P) == 1:
        return [P[0] * Q[0]]
    k = math.ceil(len(P) / 2)
    P1 = P[:k]
    P2 = P[k:]
    Q1 = Q[:k]
    Q2 = Q[k:]
    E1 = prod_poly_karatsuba(P1, Q1)
    E2 = prod_poly_karatsuba(P2, Q2)
    E3 = prod_poly_karatsuba(somme_poly(P1, P2, 0), somme_poly(Q1, Q2, 0))

    res = [0 for _ in range(len(P) + len(Q) - 1)]
    res = somme_poly(res, E1, 0)
    res = somme_poly(res, E2, 2 * int((len(P) + 1) / 2))
    res = somme_poly(res, somme_poly(E3, oppose_poly(somme_poly(E1, E2, 0)), 0), k)
    return res


def main():
    P = gen_poly_ent(8, 10)
    Q = gen_poly_ent(8, 10)
    afficher_poly(prod_poly_naif(P, Q))
    afficher_poly(prod_poly_karatsuba(P, Q))


def test_temps():
    for i in range(100, 100000, 100):
        P = gen_poly_ent(i, 100)
        Q = gen_poly_ent(i, 100)
        debut = time.time()
        prod_poly_naif(P, Q)
        fin = time.time()
        print("naif",i,fin-debut)
        debut = time.time()
        prod_poly_karatsuba(P, Q)
        fin = time.time()
        print("Kara",i,fin-debut)

# main()
test_temps()

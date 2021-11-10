# import csv
# attach worksheet.sage

# L_table = []
# bad_numbers = []
# for D in primes(20,24):
#     F.<a> = NumberField(x^2-D)
#     for p in primes(7,13):
#         if p != D:
#             k0 = 2
#             m = 10
#             psi = DirichletGroup(p, F)[0]
#             # F.dirichlet_group()[0]
#             print("Computing with data:\n", F,"\n p = ",  p, "\n precision m = ", m, "\n")
#             try:
#                 P, Q = diagonal_restriction_Lp(F, p, k0, psi, m)
#                 L_table.append([D, p, P])
#             except:
#                 print("Error computing diagonal restriction\n")
#                 bad_numbers.append([D, p])


# table(L_table)

# For making nice latex tables:
# n = 3

# L3_table = [L[1:3] for L in L_table[:6]]
# latex(table([["$p$",  "$L_p$"]] + L3_table, frame=True, header_row=True))
# this gives [ p | L_p]

# n = 5

# L5_table = [L[1:3] for L in L_table[6:11]]
# latex(table([["$p$",  "$L_p$"]] + L5_table, frame=True, header_row=True))

F.<a> = NumberField(x^2-12)
p = 5
k0 = 2
m = 6
r1 = F.signature()[0]       # nr of real embeddings = 2 
mfrak = F.modulus(F.ideal(4), range(r1))
H = HeckeCharacterGroup(mfrak)
G = F.ray_class_group(mfrak)

def character_table(G,H):
    r"""
    Returns table
    0    | g1       | g2      | ...
    chi1 | chi1(g1) | chi1(g2)| ...
    chi2 | chi2(g1) | chi2(g2)| ...
    ...
    
    """
    gs = [g.ideal().gens() for g in G.gens()]
    chis = list(H.gens())
    T = [[0] + gs]
    # print(chis)
    for chi in chis:
        T.append([chi] + [chi(g) for g in gs])
    # print(table(T))
    return(T)

# try:
#     psi = G.gens()[0]
# except:
#     psi = G.one()
#     # F.dirichlet_group()[0]
# print("F = ", F, "\nG = ", G.group(), "\np = ", p, "\nk0 = ", k0,
#       "\nm = ", m, "\npsi = ", psi)
# P, Q = diagonal_restriction_Lp(F, p, k0, psi, m,verbose=True);
# # print("P(s) = ", P)
# # print("Q(T) = ", Q)

# with open('pAdicLvalues.csv', newline='') as csvfile:

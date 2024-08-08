import numpy as np
import cvxpy as cp
import mosek

def Psatz_data_reader(filename):
    with open(filename, 'r') as file:
        content = file.read()
    Hleft, Hright, Hshapes, F, Fshapes, g = content.split("\n")
    Hleft = eval(Hleft)
    Hright = eval(Hright)
    Hshapes = eval(Hshapes)
    F = eval(F)
    Fshapes = eval(Fshapes)
    g = np.array(eval(g))
    return Hleft, Hright, Hshapes, F, Fshapes, g
path = "NCdata_p.txt"
Hleft, Hright, Hshapes, F, Fshapes, g = Psatz_data_reader(path)

# H_i in H <-> q_i, F_i in F <-> X_i

entry_num = len(g)
qleftlist = []
qrightlist = []
Xlist = []
vleftlist = []
vrightlist = []

for i in range(len(Hleft)):
    m = Hshapes[i][0]
    n = Hshapes[i][1]
    qleft = cp.Variable((n, 1))
    vleft = cp.Variable((m, 1))
    vleftlist.append(vleft)
    qleftlist.append(qleft)

for i in range(len(Hright)):
    m = Hshapes[i][0]
    n = Hshapes[i][1]
    qright = cp.Variable((n, 1))
    vright = cp.Variable((m, 1))
    vrightlist.append(vright)
    qrightlist.append(qright)

# recoverF = None

for j in range(len(F)):
    dim = Fshapes[j]
    recoverF = np.zeros([dim, dim])
    X = cp.Variable((dim, dim), symmetric=True)
    Xlist.append(X)

# for i in range(entry_num):
#     for a in range(len(F[0][i])):
#         recoverF[ F[0][i][a][0] - 1, F[0][i][a][1] - 1 ] = F[0][i][a][2]


constr = []

for i in range(len(Xlist)):
    constr += [Xlist[i] >> 0]
# H consist of list of linear transformation for each row of v, Hq = v
for i in range(len(vleftlist)):
    for j in range(entry_num):
        constr += [vleftlist[i][j] == 
                sum([Hleft[i][j][k][1] * qleftlist[i][Hleft[i][j][k][0] - 1] for k in range(len(Hleft[i][j]))])
                ]
        
for i in range(len(vrightlist)):
    for j in range(entry_num):
        constr += [vrightlist[i][j] == 
                sum([Hright[i][j][k][1] * qrightlist[i][Hright[i][j][k][0] - 1] for k in range(len(Hright[i][j]))])
                ]

print(len(vleftlist), len(vrightlist), len(Xlist))

# -------- For general case ----------------- #
# F[k][j] has entry form [x, y, f[x, y]] where f[x, y] is the coefficient
for j in range(entry_num):
    constr += [
        sum([vleftlist[i][j] for i in range(len(Hleft))]) + sum([vrightlist[i][j] for i in range(len(Hright))]) + sum([sum([Xlist[k][F[k][j][a][0] - 1, F[k][j][a][1] - 1] * F[k][j][a][2] for a in range(len(F[k][j])) if F[k][j][a][2] != 0]) for k in range(len(F)) ]) == g[j][0]
    ]



nullObj = cp.Maximize(0)
prob = cp.Problem(nullObj, constr)
prob.solve(solver=cp.MOSEK, verbose=True, mosek_params={"MSK_DPAR_OPTIMIZER_MAX_TIME": 300000.0, "MSK_IPAR_INTPNT_MAX_ITERATIONS": 100000, "MSK_DPAR_INTPNT_CO_TOL_NEAR_REL": 50000})

if prob.status == 'optimal':
    print("Problem state: optimal, find -1 in SOS + ideal(h)\n",
        "\bNo perfect quantum strategy exists for CSP"
        )
    # for X in Xlist:
    #     print(X.value)
    # for v in vlist:
    #     print(v.value)
    with open("result_p_dup.txt", 'w') as result:
        result.write("Imperfect\n")
        result.write("----- Refutation -----\n")
        result.write("X = \n")
        for X in Xlist:
            result.write(str(X.value.tolist()) + "\n")
        result.write("List of qleft: \n")
        result.write(str([qleftlist[i].value.tolist() for i in range(len(qleftlist))]) + "\n")
        result.write("List of qright: \n")
        result.write(str([qrightlist[i].value.tolist() for i in range(len(qrightlist))]) + "\n")
        

elif prob.status == 'infeasible':
    print("Problem state: infeasible, no -1 in SOS + ideal(h) found\n",
        "\bCSP is perfectly quantum satisfiable")
    with open("result_p_dup.txt", 'w') as result:
        result.write("Perfect")

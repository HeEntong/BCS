import numpy as np
import cvxpy as cp
import mosek

def Psatz_data_reader(filename):
    with open(filename, 'r') as file:
        content = file.read()
    H, Hshapes, F, Fshapes, g = content.split("\n")
    H = eval(H)
    Hshapes = eval(Hshapes)
    F = eval(F)
    Fshapes = eval(Fshapes)
    g = np.array(eval(g))
    return H, Hshapes, F, Fshapes, g
path = "mathematica_codes/NCdata.txt"
H, Hshapes, F, Fshapes, g = Psatz_data_reader(path)

# H_i in H <-> q_i, F_i in F <-> X_i

entry_num = len(g)
qlist = []
Xlist = []
vlist = []

for i in range(len(H)):
    m = Hshapes[i][0]
    n = Hshapes[i][1]
    q = cp.Variable((n, 1))
    v = cp.Variable((m, 1))
    vlist.append(v)
    qlist.append(q)

for j in range(len(F)):
    dim = Fshapes[j]
    X = cp.Variable((dim, dim), symmetric=True)
    Xlist.append(X)

constr = []

for i in range(len(Xlist)):
    constr += [Xlist[i] >> 0]
# H consist of list of linear transformation for each row of v, Hq = v
for i in range(len(vlist)):
    for j in range(entry_num):
        constr += [vlist[i][j] == 
                sum([H[i][j][k][1] * qlist[i][H[i][j][k][0] - 1] for k in range(len(H[i][j]))])
                ]

# F[k][j] has entry form [x, y, f[x, y]] where f[x, y] is the coefficient
for j in range(entry_num):
    constr += [
        sum([vlist[i][j] for i in range(len(H))]) + sum([
            sum([Xlist[k][F[k][j][a][0] - 1, F[k][j][a][1] - 1] * F[k][j][a][2] for a in range(len(F[k][j]))]) for k in range(len(F))
        ]) == g[j]
    ]


nullObj = cp.Maximize(0)
prob = cp.Problem(nullObj, constr)
prob.solve(solver=cp.MOSEK, verbose=True, mosek_params={"MSK_DPAR_OPTIMIZER_MAX_TIME": 30000000.0, "MSK_IPAR_INTPNT_MAX_ITERATIONS": 10000000, "MSK_DPAR_INTPNT_CO_TOL_NEAR_REL": 5000})

if prob.status == 'optimal':
    print("Problem state: optimal, find -1 in SOS + ideal(h)\n",
        "\bNo perfect quantum strategy exists for CSP"
        )
    # for X in Xlist:
    #     print(X.value)

else:
    print("Problem state: infeasible, no -1 in SOS + ideal(h) found\n",
        "\bCSP is perfectly quantum satisfiable")
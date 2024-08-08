import numpy as np
import cvxpy as cp
import mosek

def Psatz_data_reader(filename):
    with open(filename, 'r') as file:
        content = file.read()
    H, F, g = content.split("\n")
    H = eval(H)
    F = eval(F)
    H = [np.array(H[i]) for i in range(len(H))]
    F = [np.array(F[i]) for i in range(len(F))]
    g = np.array(eval(g))
    return H, F, g
path = "mathematica_codes/data.txt"
H, F, g = Psatz_data_reader(path)

# H_i in H <-> q_i, F_i in F <-> X_i

entry_num = len(g)
qlist = []
Xlist = []
vlist = []

for i in range(len(H)):
    m = len(H[i])
    n = len(H[i][0])
    q = cp.Variable((n, 1))
    v = cp.Variable((m, 1))
    vlist.append(v)
    qlist.append(q)

for j in range(len(F)):
    dim = F[j].shape[1]
    X = cp.Variable((dim, dim), symmetric=True)
    Xlist.append(X)

constr = []

for i in range(len(Xlist)):
    constr += [Xlist[i] >> 0]
for i in range(len(vlist)):
    constr += [vlist[i] == H[i] @ qlist[i]]

for j in range(entry_num):
    constr += [
        sum([vlist[i][j] for i in range(len(H))]) + sum([
            cp.trace(F[k][j] @ Xlist[k]) for k in range(len(F))
        ]) == g[j]
    ]

nullObj = cp.Maximize(0)
prob = cp.Problem(nullObj, constr)
prob.solve(solver=cp.MOSEK, verbose=True)

for i in range(len(Xlist)):
    print(Xlist[i].value)

for i in range(len(qlist)):
    print(qlist[i].value)
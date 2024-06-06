LoadPackage("GBNP");

edges := [[1,2], [2,3], [3,1], [4,5], [5,6], [6,4], [1,4], [2,5], [3,6]];
colors := 3;

size := Maximum(Concatenation(edges)); # Number of vertices
n := size*colors; # Number of variables

A := FreeAssociativeAlgebraWithOne(Rationals, n , "v");
x := GeneratorsOfAlgebra(A){[2..n+1]};
I := One(A);

GBNP.ConfigPrint(A);

id := function (v, c)
    return (v-1)*colors + c;
end;;

commutator := function (x, y)
    return x*y-y*x;
end;;

constraints := function (edges)
    local ip, i, j, k, s, e;
    ip := [];
    
    for i in [1..size] do
        s := -I;
        for j in [1..colors] do
            s := s + x[id(i,j)];
            Add(ip, -x[id(i,j)]+x[id(i,j)]^2);
        #    for k in [j+1..colors] do
        #        Add(ip, commutator(x[id(i,j)], x[id(i,k)]));
        #    od;
        od;
        Add(ip, s);
    od;
    
    for e in edges do
        i := e[1]; j := e[2];
        for k in [1..colors] do
            Add(ip, x[id(i,k)]*x[id(j,k)]);
            Add(ip, x[id(j,k)]*x[id(i,k)]);
        od;
    od;
    return GP2NPList(ip);
end;;

commutativity := function (i, j, gb)
    local f, g;
    f := GP2NP(commutator(x[i], x[j]));
    g := StrongNormalFormNP(f, gb);
    return Size(Concatenation(g)) = 0;
end;;

proof := function (f, ip)
    local g, gbt;
    gbt := SGrobnerTrace(ip);
    g := StrongNormalFormTraceDiff(f, gbt);
    return g;
end;;

Print("Edges: ", Size(edges), "\n", edges, "\n\n");

ip := constraints(edges);
Print("Constraints: ", Size(ip), "\n");
PrintNPList(ip);
Print("\n");

gb := SGrobner(ip);
Print("Grobner basis: ", Size(gb), ", ", IsGrobnerBasis(gb), "\n");
PrintNPList(gb);
Print("\n");

Print("Commutativity: ", commutativity(id(1,1), id(5,1), gb), "\n\n");

Print("Commutativity: ", commutativity(id(1,1), id(5,2), gb), "\n\n");

# Print("Proof:\n");
# f := GP2NP(commutator(x[id(1,1)], x[id(5,2)]));
# g := proof(f, ip);
# PrintTraceList([g]);

(* File group.m (C) 1995  William Paulsen.   Used by permission. *)

InitGroup::usage = "InitGroup[e] initializes the group, using e as the
                    identity element.  This clears any previously defined
                    groups.  The dot is used to separate letters in the 
                    words.";

InitGroup[a_] := Module[{x, y, z, n},
   Unprotect[Dot,Power];
   ClearAll[Dot,Power];
   z_ . a := z;
   a . z_ := z;
   z_ . (z_^(-1)) := a;
   (z_^(-1)) . z_ := a;
   a ^(-1) := a;
   Dot[x_ ,Dot[y_ , z_]] := Dot[Dot[x , y] , z];
   Dot[x__, y_, z_] := Dot[Dot[x , y] , z];
   x_ . y_ . (y_^(-1)) := x;
   x_ . (y_^(-1)) . y_ := x;
   (x_ . y_)^(-1) := (y^(-1)) . (x^(-1));
   x_^n_Integer := x . (x^(n-1)) /; n > 1;
   x_^n_Integer := (x^(-1)) . (x^(n+1)) /; n < -1;
   g_List . h_List := Module[{i,j}, 
   Union[Flatten[Table[g[[i]] . h[[j]],{i,Length[g]},{j,Length[h]}],1]]];
   g_List . x_ := Module[{i}, Union[Table[g[[i]] . x,{i,Length[g]}]]]; 
   x_ . g_List := Module[{i}, Union[Table[x . g[[i]],{i,Length[g]}]]];
   g_List^(-1) := Module[{i}, Union[Table[g[[i]]^(-1),{i,Length[g] }]]];
   Ident$ = a ]

Define::usage = "Define[u, v], where u and v are words, defines the inference
                 rule u -> v.";

Define[a_ . b_ , c_] := Module[{},
   a . b := c;
   Dot[Dot[z$_ , a] , b] := z$ . c ]

DefInvert::usage = "DefInvert[u, v], where u is a generator, defines the
                    inverse of u to be the word v.";

DefInvert[a_ , b_] := a^(-1) := b;

ResetDot::usage = "ResetDot[] resets the dot product, erasing all groups.";

ResetDot[] := Module[{},
   ClearAll[Dot,Power];
   Protect[Dot,Power] ]

Group::usage = "Group[X_List] gives a list of elements of the group
                symbolically presented by the set of generators X and the 
                inference rules previously defined.";

Group[G_List] := Module[{m, n, i, j, g, Repeat},
   g = Union[G];
   m = Length[g];
   Repeat = True;
   While[Repeat,
      Repeat = False;
      n = Length[g];
      Do[ Do[ If[ MemberQ[g, g[[i]] . g[[j]] ], ,
                  Repeat = True;
                  g = Append[g, g[[i]] . g[[j]] ]
                ], {i,1,n} ], {j,1,m} ] ];
   Return[g] ]

LeftCoset::usage = "LeftCoset[G_list, H_list] gives a list of all of the left
                    cosets of the subgroup H of G.";

LeftCoset[G_List, H_List] := Module[{g, i},
    g = { };
    Do[ If[MemberQ[Flatten[g,1],G[[i]] ], ,g = Append[g, G[[i]] . H] ],
        {i,Length[G]}];
    Return[g] ]

RightCoset::usage = "RightCoset[G_list, H_list] gives a list of all of the 
                     right cosets of the subgroup H of G.";

RightCoset[G_List, H_List] := Module[{g, i},
    g = { };
    Do[ If[MemberQ[Flatten[g,1],G[[i]] ], ,g = Append[g,H . G[[i]]] ],
        {i,Length[G]}];
    Return[g] ]

NormalGroup::usage = "NormalGroup[G_List, X_List] gives the smallest normal
                      subgroup of G which contains the elements of X.";

NormalGroup[G_List, H_List] := Module[{g, m, n, i, j, Repeat},
   g = Union[Flatten[
     Table[G[[i]] . H[[j]] . (G[[i]]^(-1)),{i,Length[G]},{j,Length[H]}],1]];
   m = Length[g];
   Repeat = True;
   While[Repeat,
      Repeat = False;
      n = Length[g];
      Do[ Do[ If[ MemberQ[g,g[[i]].g[[j]] ], ,
                  Repeat = True;
                  g = Append[g,g[[i]].g[[j]] ]
                ], {j,1,m}];
          If[2*Length[g] > Length[G], 
             g = G; 
             Repeat = False; 
             Break[ ]
          ], {i,1,n} ] ];
   Return[g] ]

ConjugacyClass::usage = "ConjugacyClass[G_List] gives the list of the 
                         conjugacy classes for the group G.";

ConjugacyClass[G_List] := Module[{i, j, g},
    g = {};
    Do[ If[MemberQ[Flatten[g,1],G[[i]] ], , g = Append[g,Union[Table[
           G[[j]] . G[[i]] . (G[[j]]^(-1)),{j,Length[G]}]]]
          ], {i,Length[G]}];
    Return[g] ]
    
Commutator::usage = "Commutator[G_List, X_List] gives the commutator subgroup
                     [G, H], where H is the subgroup generated by the elements
                     of X.  If X generates G, this gives the derived group
                     G'.";

Commutator[G_List, H_List] := Module[{g, m, n, i, j, Repeat},
   g = Union[Flatten[
     Table[G[[i]] . H[[j]] . (G[[i]]^(-1)) . (H[[j]]^(-1)),
           {i,Length[G]},{j,Length[H]}],1]];
   m = Length[g];
   Repeat = True;
   While[Repeat,
         Repeat = False;
         n = Length[g];
         Do[ Do[ If[ MemberQ[g,g[[i]].g[[j]] ], ,
                     Repeat = True;
                     g = Append[g,g[[i]].g[[j]] ]
                   ], {j,1,m}];
             If[2*Length[g] > Length[G], 
                g = G; 
                Repeat = False; 
                Break[ ]
             ], {i,1,n} ] ];
   Return[g] ]

Normalizer::usage = "Normalizer[G_List, {a}] gives the normalizer of a in the
                     group G.  If H is a subgroup of G, then 
                     Normalizer[G_List, H_List] gives the largest subgroup of
                     G in which H is normal.";

Normalizer[G_List,H_List] := Module[{g, i, j, Include},
   g = {};
   Do[ Include = True;
      Do[ If[MemberQ[H, G[[i]] . H[[j]] . (G[[i]]^(-1))], ,
             Include = False;
             Break[ ] 
            ], {j,Length[H] }];
      If[Include, g = Append[g, G[[i]] ] ], {i, Length[G] }];
   Return[g] ]

GroupCenter::usage = "GroupCenter[G_List] gives the center of the group G,
                      which is a normal subgroup.";

GroupCenter[G_List] := Module[{g,i,j,Include},
   g = {};
   Do[ Include = True;
      Do[ If[G[[i]] . G[[j]] === G[[j]] . G[[i]] , ,
             Include = False;
             Break[ ] 
            ], {j, Length[G] }];
      If[Include, g = Append[g, G[[i]] ] ], {i, Length[G] }];
   Return[g] ]

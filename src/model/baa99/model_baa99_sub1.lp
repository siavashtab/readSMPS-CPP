\ENCODING=ISO-8859-1
\Problem name: IloCplex

Minimize
 obj: - 8 w11 - 4 w12 - 4 w22 + 0.2 v1 + 0.2 v2 + 10 u1 + 10 u2
Subject To
 d1: w11 + u1  = 100
 d2: w12 + w22 + u2  = 100
 s1: w11 + w12 + v1  = 0
 s2: w22 + v2  = 0
End

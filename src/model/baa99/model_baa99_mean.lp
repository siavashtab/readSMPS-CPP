\ENCODING=ISO-8859-1
\Problem name: .\Input\baa99\baa99.cor

Minimize
 obj: 4 x1 + 2 x2 - 8 w11 - 4 w12 - 4 w22 + 0.2 v1 + 0.2 v2 + 10 u1 + 10 u2
Subject To
 d1: w11 + u1  = 100
 d2: w12 + w22 + u2  = 100
 s1: - x1 + w11 + w12 + v1  = 0
 s2: - x2 + w22 + v2  = 0
Bounds
 0 <= x1 <= 217
 0 <= x2 <= 217
End

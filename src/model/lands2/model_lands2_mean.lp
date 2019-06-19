\ENCODING=ISO-8859-1
\Problem name: .\Input\lands2\lands2.cor

Minimize
 OBJ: 10 X1 + 7 X2 + 16 X3 + 6 X4 + 40 Y11 + 45 Y21 + 32 Y31 + 55 Y41 + 24 Y12
      + 27 Y22 + 19.2 Y32 + 33 Y42 + 4 Y13 + 4.5 Y23 + 3.2 Y33 + 5.5 Y43
Subject To
 S1C1: X1 + X2 + X3 + X4 >= 12
 S1C2: 10 X1 + 7 X2 + 16 X3 + 6 X4 <= 120
 S2C1: - X1 + Y11 + Y12 + Y13 <= 0
 S2C2: - X2 + Y21 + Y22 + Y23 <= 0
 S2C3: - X3 + Y31 + Y32 + Y33 <= 0
 S2C4: - X4 + Y41 + Y42 + Y43 <= 0
 S2C5: Y11 + Y21 + Y31 + Y41 >= 1.98
 S2C6: Y12 + Y22 + Y32 + Y42 >= 1.98
 S2C7: Y13 + Y23 + Y33 + Y43 >= 1.98
End

# Unconventioanl Materials （Unconv-Mat）

# 1. Standardize the POSCAR (e.g. Ca2N)

$$ phonopy  --tolerance 0.01 --symmetry -c POSCAR > phonopy.out

$$ pos2aBR > aBR.out


-----aBR.out-----

  1          166 R-3m
  
  2  The matrix Kc2p
  
  3     0.666667   -0.333333   -0.333333
  
  4     0.333333    0.333333   -0.666667
  
  5     0.333333    0.333333    0.333333
  
  6  The matrix p2cR
  
  7     1.000000    0.000000    1.000000
  
  8    -1.000000    1.000000    1.000000
  
  9     0.000000   -1.000000    1.000000
  
 10  Origin shift
 
 11     0.000000    0.000000    0.000000
 
 12  \\
 
 13   SN  Mult. Wyck. Atom  s    p    d  Wyck. Name
 
 14     1    2    7   20    2    6    0    6c   Ca
 
 15     2    1    9    7    2    3    0    3a   N
 
 16  \\
 
 17   SN  Orb. @ Site     Symm.
 
 18    1  Ca-s @ 6c( 7)    3m(19) >>>   BCS  Int. Sch.      Basis
 
 19                                  1  GM1 ;GM1 ; A1 ;     z;x2+y2;z2
 
 20    1  Ca-p @ 6c( 7)    3m(19) >>>   BCS  Int. Sch.      Basis
 
 21                                  1  GM1 ;GM1 ; A1 ;     z;x2+y2;z2
 
 22                                  3  GM3 ;GM3 ; E  ;     x,y;xz,yz;x2-y2,xy;Jx,Jy
 
 23    2   N-s @ 3a( 9)   -3m(20) >>>   BCS  Int. Sch.      Basis
 
 24                                  1  GM1+;GM1+; A1g;     x2+y2;z2
 
 25    2   N-p @ 3a( 9)   -3m(20) >>>   BCS  Int. Sch.      Basis
 
 26                                  4  GM2-;GM2-; A2u;     z
 
 27                                  6  GM3-;GM3-; Eu ;     x,y
 
-----aBR.out-----


# 2. Run scf and band calculations in VASP for Maximial High-symmetry k-points (POTCAR: PAW_PBE Ca_sv 06Sep2000; PAW_PBE N 08Apr2002) 
$$ irvsp2 -sg 166 -nb 9 13 > outir2   !! generating tqc.txt 

-----tqc.txt-----

  1 Computed bands:  9 - 13
  
  2 GM: GM1+(1); GM2-(1); GM1+(1); GM3-(2); [5]
  
  3 T : T1+ (1); T2- (1); T3- (2); T2- (1); [5]
  
  4 F : F1+ (1); F2- (1); F2- (1); F1- (1); F1+ (1); [5]
  
  5 L : L1+ (1); L2- (1); L2- (1); L1- (1); L2- (1); [5]
  
----------------

# 3. Solving the BR decomposition for the set of energy bands.
$$ Python BR_decomp.py 
  A2u@3a + Eu@3a + A1g@3b

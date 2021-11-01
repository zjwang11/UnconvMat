# Unconventioanl Materials （Unconv. Mat.）

* Website: http://tm.iphy.ac.cn/UnconvMat.html <br> 
Refs: J. Gao, et al. https://arxiv.org/abs/2106.08035. <br>
S. Nie, et al. Phys. Rev. B 103, 205133 (2021). <br>
S. Nie, et al. "Six-fold Excitations in Electrides",  Phys. Rev. Research 3, L012028 (2021).

### 1. Standardize the PPOSCAR from phonopy and generate the aBRs (e.g. Ca2N)

$$ phonopy  --tolerance 0.01 --symmetry -c POSCAR > phonopy.out

$$ pos2aBR > aBR.out


--aBR.out--

       166 R-3m
    The matrix Kc2p
      0.666667   -0.333333   -0.333333
      0.333333    0.333333   -0.666667
      0.333333    0.333333    0.333333
    The matrix p2cR
       1.000000    0.000000    1.000000
      -1.000000    1.000000    1.000000
       0.000000   -1.000000    1.000000
    Origin shift
    0.000000    0.000000    0.000000
    \\
     SN  Mult. Wyck. Atom  s    p    d  Wyck. Name
       1    2    7   20    2    6    0    6c   Ca
       2    1    9    7    2    3    0    3a   N
    \\
     SN  Orb. @ Site     Symm.
      1  Ca-s @ 6c( 7)    3m(19) >>>   BCS  Int. Sch.      Basis
                                    1  GM1 ;GM1 ; A1 ;     z;x2+y2;z2
      1  Ca-p @ 6c( 7)    3m(19) >>>   BCS  Int. Sch.      Basis
                                    1  GM1 ;GM1 ; A1 ;     z;x2+y2;z2
                                    3  GM3 ;GM3 ; E  ;     x,y;xz,yz;x2-y2,xy;Jx,Jy
      2   N-s @ 3a( 9)   -3m(20) >>>   BCS  Int. Sch.      Basis
                                    1  GM1+;GM1+; A1g;     x2+y2;z2
      2   N-p @ 3a( 9)   -3m(20) >>>   BCS  Int. Sch.      Basis
                                    4  GM2-;GM2-; A2u;     z
                                    6  GM3-;GM3-; Eu ;     x,y


### 2. Run scf and band calculations in VASP for maximial High-symmetry k-points (HSKPs)

(POTCAR: PAW_PBE Ca_sv 06Sep2000; PAW_PBE N 08Apr2002) 

$$ irvsp2 -sg 166 -nb 9 13 > outir2    (* generating tqc.txt and tqc.data \*)

--tqc.txt--

    Computed bands:  9 - 13
    GM: GM1+(1); GM2-(1); GM1+(1); GM3-(2); [5]
    T : T1+ (1); T2- (1); T3- (2); T2- (1); [5]
    F : F1+ (1); F2- (1); F2- (1); F1- (1); F1+ (1); [5]
    L : L1+ (1); L2- (1); L2- (1); L1- (1); L2- (1); [5]


--tqc.data--

    166    4    5       # space group number, number of k-points, number of energy bands (nb)
    1  1  4  1  6       # k index, irrep index [nb]
    2  1  4  6  4
    4  1  4  4  2  1
    5  1  4  4  2  4

### 3. Solving the aBR/eBR decomposition for the set of energy bands online 

Website: http://tm.iphy.ac.cn/UnconvMat.html <br> Ref: J. Gao, et al. https://arxiv.org/abs/2106.08035.

$$ Paste "PPOSCAR" and "tqc.data" the website.

$$ Press the button "BR" to solve the eBR and aBR decomposition.


    A1g@3a (N-s) + A2u@3a (N-Pz) + Eu@3a (N-Px,Py) + A1g@3b (an empty site)
  

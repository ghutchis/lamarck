%nproc=1
%mem=1GB
%Chk=3.chk
#T PM6 OPT

c(c(nsn1)c12)ccc2c(c(nsn1)c12)ccc2c(c(nsn1)c12)ccc2c(c(nsn1)c12)ccc2

0  1
C          -0.44026        -1.90416        -0.31940
C          -1.46006        -1.28170         0.42968
N          -2.41701        -1.93396         1.11124
S          -3.35914        -0.77528         1.82503
N          -2.59775         0.60291         1.31324
C          -1.56672         0.23097         0.53231
C           0.44737        -1.04076        -0.96928
C           0.34406         0.37137        -0.89369
C          -0.63393         1.08618        -0.14581
C          -0.69962         2.55964        -0.06909
C           0.46507         3.39204        -0.11678
N           1.74034         2.97166        -0.20225
S           2.70811         4.30512        -0.26602
N           1.57980         5.50216        -0.16873
C           0.36840         4.92025        -0.06052
C          -1.89987         3.29657         0.05326
C          -1.98507         4.71211         0.14449
C          -0.88347         5.59750         0.07884
C          -0.99569         7.06903         0.14212
C          -2.10234         7.79238        -0.40226
N          -3.17768         7.25841        -1.01472
S          -4.17557         8.49576        -1.45117
N          -3.30071         9.79102        -0.92725
C          -2.17674         9.32417        -0.35244
C          -0.02738         7.91038         0.73904
C          -0.09230         9.32859         0.78201
C          -1.13709        10.11194         0.24043
C          -1.16773        11.58839         0.27917
C           0.02214        12.39044         0.20773
N           1.29179        11.96043         0.09056
S           2.27200        13.29443         0.07156
N           1.15023        14.50563         0.19061
C          -0.05180        13.90785         0.25277
C          -2.35902        12.35849         0.39343
C          -2.41978        13.77511         0.41799
C          -1.28313        14.58747         0.35212
H          -0.35603        -2.98195        -0.38531
H           1.25588        -1.47061        -1.56152
H           1.08702         0.93221        -1.46443
H          -2.85173         2.76142         0.08084
H          -2.98897         5.12090         0.27238
H           0.85156         7.46142         1.20625
H           0.74496         9.82680         1.27500
H          -3.31747        11.84144         0.47138
H          -3.39839        14.25079         0.49211
H          -1.33893        15.66886         0.37621


--Link1--
%nproc=1
%mem=1GB
%Chk=3.chk
%NoSave
# Geom=AllCheck ZINDO(NStates=15,Singlets)

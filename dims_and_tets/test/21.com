%nproc=1
%mem=1GB
%Chk=21.chk
#T PM6 OPT

c(s1)cc(c12)c(=O)c(c2s3)cc3c(s1)cc(c12)c(=O)c(c2s3)cc3c(s1)cc(c12)c(=O)c(c2s3)cc3c(s1)cc(c12)c(=O)c(c2s3)cc3

0  1
C           0.84930        -3.87004        -2.17279
S           1.19255        -2.78784        -0.86452
C           1.88964        -3.94988        -3.08631
C           2.95677        -3.10739        -2.69956
C           2.68510        -2.44489        -1.52089
C           4.26783        -2.76366        -3.22162
O           4.83682        -3.16168        -4.21746
C           4.73628        -1.79902        -2.19544
C           3.71400        -1.67338        -1.18983
S           4.16185        -0.62508         0.05805
C           5.83943        -1.07241        -2.02106
C           5.70125        -0.31623        -0.75426
C           6.62765         0.61557        -0.19145
S           7.65948         1.52066        -1.24874
C           6.80555         0.97984         1.14011
C           7.70938         2.05828         1.25590
C           8.20406         2.43958         0.02755
C           8.11844         2.97163         2.30710
O           7.86663         2.95344         3.49871
C           8.75647         4.03433         1.49780
C           8.91696         3.55649         0.14679
S           9.29412         4.82573        -0.91446
C           8.82343         5.35276         1.69991
C           8.88779         5.98370         0.35724
C           8.41226         7.30249         0.05265
S           7.80867         7.66003        -1.53712
C           8.28424         8.39897         0.90369
C           7.79156         9.53028         0.20660
C           7.50354         9.22768        -1.10432
C           7.59572        10.94454         0.48140
O           7.74288        11.56212         1.51591
C           7.28026        11.44539        -0.88278
C           7.11759        10.31544        -1.75786
S           6.90355        10.76765        -3.37619
C           7.30852        12.64412        -1.46426
C           7.22508        12.44110        -2.93141
C           7.54074        13.39054        -3.94931
S           8.80484        14.54115        -3.65764
C           6.99876        13.52589        -5.22204
C           7.60261        14.59919        -5.91686
C           8.58530        15.20501        -5.16553
C           7.45464        15.22796        -7.21785
O           6.68767        14.96661        -8.12162
C           8.49503        16.28623        -7.14067
C           9.13063        16.20218        -5.85151
S          10.34606        17.36166        -5.65375
C           8.94843        17.23249        -7.96234
C          10.03588        17.96480        -7.27807
H          -0.09858        -4.39062        -2.19871
H           1.87301        -4.57408        -3.97156
H           6.69969        -1.04135        -2.66936
H           6.24845         0.55619         1.96770
H           8.53714         5.86863         2.60216
H           8.58831         8.40411         1.94401
H           7.51901        13.58840        -0.99047
H           6.20284        12.90303        -5.61401
H           8.61068        17.44487        -8.96324
H          10.56346        18.78542        -7.78005


--Link1--
%nproc=1
%mem=1GB
%Chk=21.chk
%NoSave
# Geom=AllCheck ZINDO(NStates=15,Singlets)

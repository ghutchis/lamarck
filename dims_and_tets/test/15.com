%nproc=1
%mem=1GB
%Chk=15.chk
#T PM6 OPT

c(s1)cc(c12)[nH]c(c2s3)cc3c(s1)cc(c12)[nH]c(c2s3)cc3c(s1)cc(c12)[nH]c(c2s3)cc3c(s1)cc(c12)[nH]c(c2s3)cc3

0  1
C           1.66474         1.28938        -0.04059
S           0.08098         2.02003        -0.03843
C           1.66631        -0.10203        -0.02188
C           0.32711        -0.56642        -0.00573
C          -0.53821         0.47122        -0.01275
N          -0.44902        -1.69141         0.01459
H          -0.09608        -2.63905         0.02373
C          -1.77212        -1.34777         0.01943
C          -1.79956        -0.00105         0.00288
S          -3.28105         0.75672         0.00330
C          -3.08717        -1.87485         0.03424
C          -4.02754        -0.83603         0.02585
C          -5.48088        -0.96223         0.03194
S          -6.22672        -2.55559         0.05354
C          -6.42201         0.07605         0.02054
C          -7.73679        -0.45217         0.02379
C          -7.70837        -1.79881         0.04040
N          -9.06043        -0.11124         0.01407
H          -9.41534         0.83607        -0.00070
C          -9.83361        -1.23842         0.01884
C          -8.96894        -2.27237         0.03884
S          -9.58699        -3.81874         0.05132
C         -11.17109        -1.70600         0.00769
C         -11.18975        -3.10543         0.01562
C         -12.36357        -3.96128        -0.01433
S         -13.82839        -3.40236        -0.80301
C         -12.48279        -5.24788         0.52364
C         -13.78132        -5.75956         0.27450
C         -14.52295        -4.86234        -0.40564
N         -14.59489        -6.84259         0.46159
H         -14.32901        -7.69133         0.94385
C         -15.81941        -6.61178        -0.10153
C         -15.74679        -5.37565        -0.63435
S         -17.08705        -4.77274        -1.41747
C         -17.09976        -7.17774        -0.32631
C         -17.90269        -6.28023        -1.03989
C         -19.27847        -6.48363        -1.46248
S         -20.34266        -7.51053        -0.51742
C         -19.90213        -5.92991        -2.58649
C         -21.25153        -6.35834        -2.65399
C         -21.54497        -7.17136        -1.61915
N         -22.40228        -6.26107        -3.38583
H         -22.50654        -5.71826        -4.23295
C         -23.39227        -7.00297        -2.80453
C         -22.83243        -7.55873        -1.70714
S         -23.77840        -8.54274        -0.74770
C         -24.75039        -7.39152        -2.92624
C         -25.08047        -8.22941        -1.86542
H           2.52533         1.94423        -0.05593
H           2.55908        -0.71313        -0.02039
H          -3.32860        -2.92956         0.04802
H          -6.18085         1.13115         0.00731
H         -12.05264        -1.07801        -0.00667
H         -11.69551        -5.76931         1.05324
H         -17.41427        -8.16173        -0.00269
H         -19.41810        -5.27455        -3.29924
H         -25.43699        -7.09541        -3.70828
H         -26.03856        -8.69249        -1.67153


--Link1--
%nproc=1
%mem=1GB
%Chk=15.chk
%NoSave
# Geom=AllCheck ZINDO(NStates=15,Singlets)

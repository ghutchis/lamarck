%nproc=1
%mem=1GB
%Chk=13.chk
#T PM6 OPT

c(o1)ccc1c(o1)ccc1c(o1)ccc1c(o1)ccc1

0  1
C          -0.31318         1.18340        -0.10121
O           0.98230         0.78154        -0.07482
C          -1.15315         0.09763        -0.12616
C          -0.31457        -1.04285        -0.11429
C           0.99432        -0.59087        -0.08147
C           2.20043        -1.36200        -0.05171
O           2.21096        -2.73200        -0.05360
C           3.50596        -0.90950        -0.01039
C           4.33989        -2.05215         0.02125
C           3.51212        -3.15839        -0.00501
C           3.87760        -4.54276         0.02294
O           5.17725        -4.97086         0.09344
C           3.04796        -5.64763        -0.00255
C           3.87821        -6.79102         0.06400
C           5.18373        -6.34072         0.12547
C           6.38459        -7.11518         0.21695
O           6.38728        -8.48628         0.27976
C           7.69501        -6.66931         0.26438
C           8.52479        -7.81229         0.36115
C           7.67832        -8.89327         0.36823
H          -0.47028         2.25304        -0.09765
H          -2.23336         0.12394        -0.14890
H          -0.61814        -2.08077        -0.12536
H           3.81095         0.12799         0.00166
H           5.42031        -2.07550         0.06356
H           1.96799        -5.62244        -0.05560
H           3.57047        -7.82781         0.07246
H           8.00533        -5.63378         0.23389
H           9.60350        -7.84335         0.42008
H           7.82761        -9.96235         0.42843


--Link1--
%nproc=1
%mem=1GB
%Chk=13.chk
%NoSave
# Geom=AllCheck ZINDO(NStates=15,Singlets)
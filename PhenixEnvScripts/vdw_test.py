#%%

from __future__ import print_function

models = [
'''
CRYST1   55.738  117.016   69.874  90.00  92.43  90.00 P 1 21 1
SCALE1      0.017941  0.000000  0.000761        0.00000
SCALE2      0.000000  0.008546  0.000000        0.00000
SCALE3      0.000000  0.000000  0.014324        0.00000
ATOM      1  O   LYS A 566      57.339 -28.468  36.618  1.00 15.17           O
ANISOU    1  O   LYS A 566     2141   1693   1929     62     57    373       O
ATOM      2  CB  LYS A 566      57.923 -29.004  33.451  1.00 13.89           C
ANISOU    2  CB  LYS A 566     1975   1410   1891     95     53    217       C
ATOM      3  CG  LYS A 566      58.715 -28.298  32.334  1.00 16.05           C
ANISOU    3  CG  LYS A 566     2252   1696   2152    120     48    145       C
ATOM      4  CD  LYS A 566      58.238 -28.682  30.927  1.00 13.76           C
ANISOU    4  CD  LYS A 566     1957   1381   1892     97     58     71       C
ATOM      5  CE  LYS A 566      56.720 -28.302  30.684  1.00 15.09           C
ANISOU    5  CE  LYS A 566     2116   1591   2027     60     61     61       C
ATOM      6  NZ  LYS A 566      56.537 -26.846  30.270  1.00 16.27           N
ANISOU    6  NZ  LYS A 566     2262   1810   2109     70     50     39       N
ATOM      7  CB  VAL A 569      59.871 -26.460  38.172  1.00 17.22           C
ANISOU    7  CB  VAL A 569     2415   2082   2044    113     13    395       C
ATOM      8  CG1 VAL A 569      59.537 -24.942  38.059  1.00 22.77           C
ANISOU    8  CG1 VAL A 569     3118   2847   2687    117     33    320       C
ATOM      9  CG  LYS A 570      54.873 -26.375  36.871  1.00 15.92           C
ANISOU    9  CG  LYS A 570     2205   1926   1919     25    116    287       C
ATOM     10  CD  LYS A 570      53.914 -27.000  35.821  1.00 19.04           C
ANISOU   10  CD  LYS A 570     2585   2274   2377      8    120    260       C
ATOM     11  CE  LYS A 570      53.950 -26.155  34.574  1.00 19.96           C
ANISOU   11  CE  LYS A 570     2695   2390   2500     33    106    192       C
ATOM     12  NZ  LYS A 570      53.230 -26.881  33.455  1.00 17.67           N
ANISOU   12  NZ  LYS A 570     2392   2060   2262      9     99    165       N
ATOM     13  CG2 VAL A 791      62.776 -24.286  35.819  1.00 20.32           C
ANISOU   13  CG2 VAL A 791     2805   2456   2459    213    -18    233       C
ATOM     14  CG  ARG A 792      61.559 -25.134  31.781  1.00 20.02           C
ANISOU   14  CG  ARG A 792     2765   2310   2531    199     22     81       C
ATOM     15  CD  ARG A 792      60.589 -24.195  31.092  1.00 27.55           C
ANISOU   15  CD  ARG A 792     3719   3303   3446    181     25     49       C
ATOM     16  NE  ARG A 792      60.714 -23.959  29.670  1.00 18.50           N
ANISOU   16  NE  ARG A 792     2572   2162   2296    170     25     -2       N
ATOM     17  CZ  ARG A 792      59.667 -23.661  28.912  1.00 18.64           C
ANISOU   17  CZ  ARG A 792     2585   2201   2298    144     21    -19       C
ATOM     18  NH1 ARG A 792      58.447 -23.694  29.478  1.00 17.73           N
ANISOU   18  NH1 ARG A 792     2460   2091   2185    131     22      5       N
ATOM     19  CD1 LEU A 808      59.912 -20.380  30.793  1.00 14.84           C
ANISOU   19  CD1 LEU A 808     2096   1794   1749    191     10     24       C
ATOM     20  NH2 ARG A 812      52.559 -21.925  30.586  1.00 18.04           N
ANISOU   20  NH2 ARG A 812     2371   2222   2261     80     31     52       N
ATOM     21  C   MET A 856      55.537 -21.766  37.868  1.00 16.00           C
ANISOU   21  C   MET A 856     2194   2099   1785     84    148    137       C
ATOM     22  O   MET A 856      55.905 -22.676  37.098  1.00 16.66           O
ANISOU   22  O   MET A 856     2292   2131   1906     91    121    164       O
ATOM     23  CB  MET A 856      57.849 -21.805  38.929  1.00 16.40           C
ANISOU   23  CB  MET A 856     2285   2184   1763     94    109    181       C
ATOM     24  N   TYR A 857      54.329 -21.248  37.829  1.00 14.55           N
ANISOU   24  N   TYR A 857     1978   1930   1620     76    178    102       N
ATOM     25  CA  TYR A 857      53.285 -21.719  36.905  1.00 15.53           C
ANISOU   25  CA  TYR A 857     2080   2019   1802     71    175    102       C
ATOM     26  C   TYR A 857      53.657 -21.594  35.444  1.00 16.83           C
ANISOU   26  C   TYR A 857     2252   2136   2008     97    136     89       C
ATOM     27  O   TYR A 857      53.033 -22.207  34.568  1.00 16.94           O
ANISOU   27  O   TYR A 857     2255   2123   2059     85    123     93       O
ATOM     28  N   PHE A 858      54.626 -20.704  35.169  1.00 14.65           N
ANISOU   28  N   PHE A 858     1988   1857   1720    126    119     68       N
ATOM     29  CA  PHE A 858      55.043 -20.398  33.802  1.00 15.33           C
ANISOU   29  CA  PHE A 858     2080   1913   1832    144     86     55       C
ATOM     30  C   PHE A 858      56.060 -21.363  33.236  1.00 17.20           C
ANISOU   30  C   PHE A 858     2348   2120   2066    144     64     71       C
ATOM     31  O   PHE A 858      56.433 -21.239  32.064  1.00 18.02           O
ANISOU   31  O   PHE A 858     2457   2208   2182    150     43     55       O
ATOM     32  CB  PHE A 858      55.626 -18.978  33.792  1.00 15.30           C
ANISOU   32  CB  PHE A 858     2073   1917   1825    170     81     28       C
ATOM     33  CG  PHE A 858      56.687 -18.779  34.842  1.00 14.67           C
ANISOU   33  CG  PHE A 858     2014   1858   1701    175     91     26       C
ATOM     34  CD1 PHE A 858      58.016 -19.129  34.574  1.00 16.90           C
ANISOU   34  CD1 PHE A 858     2325   2128   1967    185     68     39       C
ATOM     35  CD2 PHE A 858      56.354 -18.204  36.101  1.00 15.68           C
ANISOU   35  CD2 PHE A 858     2129   2025   1804    166    125      3       C
ATOM     36  CE1 PHE A 858      58.971 -18.948  35.533  1.00 16.58           C
ANISOU   36  CE1 PHE A 858     2298   2114   1888    186     69     42       C
ATOM     37  N   ALA A 859      56.500 -22.317  34.041  1.00 16.14           N
ANISOU   37  N   ALA A 859     2231   1981   1922    136     71    102       N
ATOM     38  CA  ALA A 859      57.528 -23.252  33.600  1.00 21.06           C
ANISOU   38  CA  ALA A 859     2874   2564   2563    143     54    114       C
ATOM     39  C   ALA A 859      56.960 -24.272  32.645  1.00 24.43           C
ANISOU   39  C   ALA A 859     3298   2951   3035    123     52    104       C
ATOM     40  O   ALA A 859      55.744 -24.471  32.554  1.00 25.15           O
ANISOU   40  O   ALA A 859     3372   3047   3138     98     61    104       O
ATOM     41  CB  ALA A 859      58.133 -23.957  34.791  1.00 22.83           C
ANISOU   41  CB  ALA A 859     3110   2789   2776    141     55    165       C
ATOM     42  OXT ALA A 859      57.727 -24.901  31.935  1.00 18.85           O
ANISOU   42  OXT ALA A 859     2601   2207   2354    128     43     89       O
TER
ATOM     43  P     U B   1      52.419 -26.766  29.880  1.00 15.68           P
ANISOU   43  P     U B   1     1871   1958   2130    -53      5   -321       P
ATOM     44  OP2   U B   1      52.377 -25.310  29.540  1.00 15.13           O
ANISOU   44  OP2   U B   1     1799   1878   2071     -7      1   -367       O
ATOM     45  OP3   U B   1      53.835 -27.127  30.208  1.00 17.49           O
ANISOU   45  OP3   U B   1     2121   2176   2350    -61     -0   -238       O
ATOM     46  OP1   C B   3      57.650 -26.667  27.771  1.00 15.18           O
ANISOU   46  OP1   C B   3     1903   1758   2105     39    -43   -161       O
TER
HETATM   47  O   HOH A1113      53.094 -24.204  32.121  1.00 19.70           O
HETATM   48  O   HOH A1158      55.980 -21.792  29.414  1.00 19.04           O
HETATM   49  O   HOH B 209      54.279 -23.804  28.528  1.00 25.28           O
''',
'''
CRYST1   55.738  117.016   69.874  90.00  92.43  90.00 P 1 21 1
SCALE1      0.017941  0.000000  0.000761        0.00000
SCALE2      0.000000  0.008546  0.000000        0.00000
SCALE3      0.000000  0.000000  0.014324        0.00000
ATOM      1  CA  ARG A 196      55.932 -14.165  -3.214  1.00 38.55           C
ANISOU    1  CA  ARG A 196     3836   7492   3318   -451   -381  -1294       C
ATOM      2  CB  ARG A 196      55.596 -14.417  -1.760  1.00 39.37           C
ANISOU    2  CB  ARG A 196     3897   7511   3552   -494   -358  -1257       C
ATOM      3  CG  ARG A 196      54.317 -13.730  -1.361  1.00 43.23           C
ANISOU    3  CG  ARG A 196     4273   8150   4002   -433   -394  -1178       C
ATOM      4  CD  ARG A 196      53.727 -14.270  -0.082  1.00 45.10           C
ANISOU    4  CD  ARG A 196     4429   8364   4344   -508   -388  -1163       C
ATOM      5  NE  ARG A 196      52.432 -13.643   0.165  1.00 42.06           N
ANISOU    5  NE  ARG A 196     3915   8167   3897   -441   -427  -1098       N
ATOM      6  CZ  ARG A 196      51.759 -13.729   1.304  1.00 48.30           C
ANISOU    6  CZ  ARG A 196     4612   8995   4744   -458   -418  -1052       C
ATOM      7  NH1 ARG A 196      52.259 -14.420   2.322  1.00 53.35           N
ANISOU    7  NH1 ARG A 196     5279   9487   5506   -548   -371  -1057       N
ATOM      8  NH2 ARG A 196      50.589 -13.109   1.427  1.00 47.30           N
ANISOU    8  NH2 ARG A 196     4361   9067   4544   -377   -455   -998       N
ATOM      9  CA  ALA A 223      59.935 -12.202   1.469  1.00 28.37           C
ANISOU    9  CA  ALA A 223     2864   5473   2444   -291    -61  -1005       C
ATOM     10  C   ALA A 223      59.529 -10.899   0.782  1.00 28.81           C
ANISOU   10  C   ALA A 223     2930   5628   2388   -180    -90   -923       C
ATOM     11  O   ALA A 223      59.279  -9.896   1.446  1.00 29.84           O
ANISOU   11  O   ALA A 223     3073   5737   2527    -89    -87   -837       O
ATOM     12  CB  ALA A 223      58.683 -12.926   1.957  1.00 35.45           C
ANISOU   12  CB  ALA A 223     3650   6449   3372   -352   -105  -1037       C
ATOM     13  N   PHE A 224      59.458 -10.942  -0.551  1.00 28.97           N
ANISOU   13  N   PHE A 224     2950   5757   2301   -181   -122   -953       N
ATOM     14  CA  PHE A 224      59.011  -9.841  -1.407  1.00 29.77           C
ANISOU   14  CA  PHE A 224     3059   5973   2281    -82   -161   -878       C
ATOM     15  C   PHE A 224      57.941 -10.348  -2.380  1.00 34.60           C
ANISOU   15  C   PHE A 224     3579   6771   2795   -105   -231   -941       C
ATOM     16  O   PHE A 224      57.885 -11.542  -2.689  1.00 37.84           O
ANISOU   16  O   PHE A 224     3951   7205   3220   -208   -245  -1053       O
ATOM     17  N   TYR A 225      57.097  -9.448  -2.862  1.00 33.76           N
ANISOU   17  N   TYR A 225     3442   6794   2592     -8   -282   -873       N
ATOM     18  CA  TYR A 225      56.197  -9.835  -3.953  1.00 35.55           C
ANISOU   18  CA  TYR A 225     3586   7216   2706    -24   -351   -928       C
ATOM     19  CB  TYR A 225      54.988  -8.889  -4.012  1.00 36.92           C
ANISOU   19  CB  TYR A 225     3695   7531   2803     96   -411   -848       C
ATOM     20  CG  TYR A 225      54.099  -9.030  -2.811  1.00 37.63           C
ANISOU   20  CG  TYR A 225     3692   7633   2973    101   -419   -845       C
ATOM     21  CD1 TYR A 225      53.144 -10.039  -2.754  1.00 45.48           C
ANISOU   21  CD1 TYR A 225     4554   8751   3974      4   -462   -926       C
ATOM     22  CD2 TYR A 225      54.215  -8.170  -1.729  1.00 39.43           C
ANISOU   22  CD2 TYR A 225     3961   7756   3265    196   -387   -761       C
ATOM     23  CE1 TYR A 225      52.328 -10.191  -1.639  1.00 44.61           C
ANISOU   23  CE1 TYR A 225     4345   8676   3928     -1   -465   -912       C
ATOM     24  CE2 TYR A 225      53.402  -8.311  -0.601  1.00 37.92           C
ANISOU   24  CE2 TYR A 225     3676   7599   3134    208   -390   -759       C
ATOM     25  CZ  TYR A 225      52.455  -9.324  -0.574  1.00 45.50           C
ANISOU   25  CZ  TYR A 225     4494   8701   4094    107   -426   -829       C
ATOM     26  OH  TYR A 225      51.644  -9.473   0.528  1.00 49.62           O
ANISOU   26  OH  TYR A 225     4908   9280   4664    110   -425   -817       O
ATOM     27  C   GLN A 350      57.579  -5.583  -1.278  1.00 37.05           C
ANISOU   27  C   GLN A 350     4078   6957   3043    337   -260   -533       C
ATOM     28  O   GLN A 350      57.408  -6.733  -1.675  1.00 34.59           O
ANISOU   28  O   GLN A 350     3690   6729   2724    235   -259   -629       O
ATOM     29  N   ARG A 351      57.940  -5.296  -0.036  1.00 31.85           N
ANISOU   29  N   ARG A 351     3461   6153   2489    357   -225   -512       N
ATOM     30  CA  ARG A 351      58.188  -6.334   0.931  1.00 33.16           C
ANISOU   30  CA  ARG A 351     3580   6254   2765    259   -178   -589       C
ATOM     31  C   ARG A 351      56.848  -6.939   1.322  1.00 37.96           C
ANISOU   31  C   ARG A 351     4046   7005   3371    258   -213   -639       C
ATOM     32  O   ARG A 351      55.831  -6.228   1.360  1.00 38.80           O
ANISOU   32  O   ARG A 351     4105   7217   3419    374   -263   -596       O
ATOM     33  CB  ARG A 351      58.919  -5.746   2.132  1.00 35.32           C
ANISOU   33  CB  ARG A 351     3936   6349   3135    293   -136   -544       C
ATOM     34  CG  ARG A 351      58.981  -6.613   3.353  1.00 43.21           C
ANISOU   34  CG  ARG A 351     4883   7292   4244    226    -95   -602       C
ATOM     35  CD  ARG A 351      59.009  -5.678   4.558  1.00 48.63           C
ANISOU   35  CD  ARG A 351     5617   7879   4982    328    -89   -548       C
ATOM     36  NE  ARG A 351      59.634  -6.257   5.726  1.00 51.71           N
ANISOU   36  NE  ARG A 351     6009   8159   5480    263    -34   -579       N
ATOM     37  N  ACYS A 352      56.835  -8.238   1.590  0.63 39.35           N
ANISOU   37  N  ACYS A 352     4155   7190   3607    128   -193   -724       N
ATOM     38  CA ACYS A 352      55.621  -8.871   2.070  0.63 41.28           C
ANISOU   38  CA ACYS A 352     4259   7564   3860     98   -223   -761       C
ATOM     39  C  ACYS A 352      55.586  -8.745   3.582  0.63 39.91           C
ANISOU   39  C  ACYS A 352     4072   7314   3779    126   -186   -733       C
ATOM     40  O  ACYS A 352      56.482  -9.225   4.267  0.63 40.81           O
ANISOU   40  O  ACYS A 352     4238   7281   3988     55   -133   -751       O
ATOM     41  CB ACYS A 352      55.552 -10.333   1.636  0.63 39.29           C
ANISOU   41  CB ACYS A 352     3947   7350   3630    -65   -232   -860       C
ATOM     42  SG ACYS A 352      54.052 -11.124   2.200  0.63 43.05           S
ANISOU   42  SG ACYS A 352     4246   7994   4117   -133   -276   -890       S
ATOM     43  N  BCYS A 352      56.832  -8.253   1.560  0.37 39.64           N
ANISOU   43  N  BCYS A 352     4190   7228   3642    126   -193   -726       N
ATOM     44  CA BCYS A 352      55.653  -8.916   2.120  0.37 41.01           C
ANISOU   44  CA BCYS A 352     4226   7523   3833     93   -220   -763       C
ATOM     45  C  BCYS A 352      55.640  -8.654   3.602  0.37 40.09           C
ANISOU   45  C  BCYS A 352     4103   7328   3802    135   -185   -728       C
ATOM     46  O  BCYS A 352      56.599  -8.973   4.297  0.37 40.97           O
ANISOU   46  O  BCYS A 352     4279   7282   4005     80   -130   -737       O
ATOM     47  CB BCYS A 352      55.656 -10.430   1.881  0.37 39.94           C
ANISOU   47  CB BCYS A 352     4033   7405   3737    -76   -221   -862       C
ATOM     48  SG BCYS A 352      56.014 -10.964   0.226  0.37 48.76           S
ANISOU   48  SG BCYS A 352     5181   8578   4768   -142   -251   -935       S
ATOM     49  N   ILE A 353      54.549  -8.086   4.090  1.00 41.84           N
ANISOU   49  N   ILE A 353     4241   7671   3985    240   -217   -691       N
ATOM     50  CA  ILE A 353      54.408  -7.801   5.517  1.00 48.04           C
ANISOU   50  CA  ILE A 353     5006   8416   4832    298   -187   -662       C
ATOM     51  C   ILE A 353      53.514  -8.846   6.193  1.00 46.98           C
ANISOU   51  C   ILE A 353     4715   8408   4727    202   -187   -695       C
ATOM     52  O   ILE A 353      53.639  -9.130   7.382  1.00 50.04           O
ANISOU   52  O   ILE A 353     5081   8746   5187    178   -146   -687       O
ATOM     53  CB  ILE A 353      53.826  -6.374   5.716  1.00 50.71           C
ANISOU   53  CB  ILE A 353     5361   8804   5102    505   -223   -599       C
ATOM     54  CG1 ILE A 353      54.876  -5.333   5.337  1.00 54.49           C
ANISOU   54  CG1 ILE A 353     6016   9106   5580    580   -219   -552       C
ATOM     55  CG2 ILE A 353      53.332  -6.147   7.142  1.00 52.31           C
ANISOU   55  CG2 ILE A 353     5502   9035   5338    585   -205   -586       C
ATOM     56  CD1 ILE A 353      56.036  -5.301   6.286  1.00 51.67           C
ANISOU   56  CD1 ILE A 353     5759   8547   5326    543   -160   -549       C
ATOM     57  N   LYS A 354      52.602  -9.420   5.422  1.00 47.21           N
ANISOU   57  N   LYS A 354     4632   8609   4695    138   -236   -726       N
ATOM     58  CA  LYS A 354      51.720 -10.450   5.951  1.00 45.71           C
ANISOU   58  CA  LYS A 354     4288   8549   4530     19   -245   -748       C
ATOM     59  C   LYS A 354      52.510 -11.703   6.324  1.00 47.90           C
ANISOU   59  C   LYS A 354     4604   8674   4920   -164   -205   -793       C
ATOM     60  O   LYS A 354      53.547 -11.993   5.715  1.00 45.77           O
ANISOU   60  O   LYS A 354     4453   8254   4682   -217   -190   -833       O
ATOM     61  CB  LYS A 354      50.638 -10.791   4.936  1.00 48.41           C
ANISOU   61  CB  LYS A 354     4510   9104   4781    -24   -316   -777       C
ATOM     62  CG  LYS A 354      49.852  -9.591   4.420  1.00 50.74           C
ANISOU   62  CG  LYS A 354     4767   9557   4954    166   -365   -734       C
ATOM     63  CD2 LEU A 356      55.803 -15.133   3.561  1.00 38.20           C
ANISOU   63  CD2 LEU A 356     3696   6982   3836   -598   -200  -1090       C
TER
HETATM   64  O   HOH A1036      54.988  -3.825   0.933  1.00 40.41           O
END
''',
]

from iotbx.data_manager import DataManager
import mmtbx

def test_each_model(i_model):
  file_name = 'tst_vdw_and_h_bond_%02d.pdb' % (i_model+1)
  f=open(file_name, 'w')
  f.write(models[i_model])
  del f
  dm = DataManager()
  dm.process_model_file(file_name)
  model = dm.get_model(file_name)

  useNeutronDistances = False
  p = mmtbx.model.manager.get_default_pdb_interpretation_params()
  p.pdb_interpretation.use_neutron_distances = useNeutronDistances
  model.process(make_restraints=True)

  # Get information on VdW radii.
  mon_lib_srv = model.get_mon_lib_srv()
  ener_lib = mmtbx.monomer_library.server.ener_lib(use_neutron_distances = useNeutronDistances)
  ph = model.get_hierarchy()
  for atom in ph.atoms():
    ag = atom.parent()
    md, ani = mon_lib_srv.get_comp_comp_id_and_atom_name_interpretation(
      residue_name=ag.resname, atom_names=ag.atoms().extract_name())
    atom_dict = md.atom_dict()
    vdw = model.get_specific_vdw_radius(atom.i_seq)
    vdwh = model.get_specific_vdw_radius(atom.i_seq, vdw_radius_without_H=True)
    h_bond = model.get_specific_h_bond_type(atom.i_seq)
    ion_radius = model.get_specific_ion_radius(atom.i_seq)
    if ion_radius is None:
      ion_radius_str = '    -'
    else:
      ion_radius_str = '%7.2f' % ion_radius
    print('  %s  %7.2f %7.2f %s %s' % (atom.quote(),
                                       vdw,
                                       vdwh,
                                       h_bond,
                                       ion_radius_str))
    name = atom.name.strip()
    try:
      te = atom_dict[name].type_energy
    except:
      continue
    vdw_radius = ener_lib.lib_atom[te].vdw_radius
    assert vdw_radius==vdw
  print('OK')

def main():
  for i in range(len(models)):
    test_each_model(i)

if __name__ == '__main__':
  main()

# %%

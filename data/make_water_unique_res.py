#%%
waters = """ATOM   1897  O  AHOH S   1       5.713   3.830  -0.595  0.50  5.34           O
ATOM   1898  O  BHOH S   1       6.007   4.060  -0.356  0.50  4.89           O
ATOM   1899  O  AHOH S   2       9.440   1.580  -6.297  0.50  6.34           O
ATOM   1900  O  BHOH S   2       9.297   1.443  -6.046  0.50  6.32           O
ATOM   1901  O  AHOH S   3      12.691   2.957   0.295  0.50  6.86           O
ATOM   1902  O  BHOH S   3      12.693   3.391   0.125  0.50  7.32           O
ATOM   1903  O  AHOH S   4      15.296   4.501   1.823  0.50 10.79           O
ATOM   1904  O  BHOH S   4      15.572   4.285   2.041  0.50  8.95           O
ATOM   1905  O  AHOH S   5     -10.020   5.405  -1.444  0.45  9.14           O
ATOM   1906  O  BHOH S   5     -10.177   4.905  -1.386  0.45  9.10           O
ATOM   1907  O  AHOH S   6      -6.147   9.873  -3.215  0.49  8.09           O
ATOM   1908  O  BHOH S   6      -5.969  10.337  -2.680  0.49  8.85           O
ATOM   1909  O  AHOH S   7       5.178   1.838  -7.042  0.50  7.06           O
ATOM   1910  O  BHOH S   7       5.495   1.938  -7.739  0.50 11.16           O
ATOM   1911  O  AHOH S   8      11.421   7.239  -5.060  0.48  6.89           O
ATOM   1912  O  BHOH S   8      12.385   6.326  -5.078  0.48 15.33           O
ATOM   1913  O  AHOH S   9      11.931  10.348  -6.664  0.41  9.50           O
ATOM   1914  O  BHOH S   9      11.662  10.853  -6.750  0.41  8.15           O
ATOM   1915  O  AHOH S  10      15.500   5.894  -0.193  0.47 10.13           O
ATOM   1916  O  BHOH S  10      15.010   5.925  -0.539  0.47 15.34           O
ATOM   1917  O  AHOH S  11      10.046   5.726  -3.511  0.47 15.56           O
ATOM   1918  O  BHOH S  11       9.677   5.059  -3.575  0.47  7.25           O
ATOM   1919  O  AHOH S  12      -3.898  10.993  -9.913  0.50 12.82           O
ATOM   1920  O  BHOH S  12      -4.083  11.464  -9.391  0.50 10.30           O
ATOM   1921  O  AHOH S  13       6.097   4.883   7.278  0.46  8.31           O
ATOM   1922  O  BHOH S  13       5.877   5.537   7.452  0.46  9.09           O
ATOM   1923  O  AHOH S  14      12.414   2.753  -3.686  0.50  8.66           O
ATOM   1924  O  BHOH S  14      12.037   1.726  -4.004  0.50 23.46           O
ATOM   1925  O  AHOH S  15       0.218  -2.490  -3.378  0.50 10.34           O
ATOM   1926  O  BHOH S  15      -0.554  -2.666  -3.441  0.50 14.62           O
ATOM   1927  O  AHOH S  16      13.568   5.149  -2.394  0.47  9.94           O
ATOM   1928  O  BHOH S  16      13.657   5.397  -3.032  0.47 14.19           O
ATOM   1929  O  AHOH S  17       2.527   6.722 -10.734  0.50 17.80           O
ATOM   1930  O  BHOH S  17       2.150   7.171 -10.554  0.50 13.50           O
ATOM   1931  O  AHOH S  18      10.610  14.477  -6.707  0.40  8.45           O
ATOM   1932  O  BHOH S  18      10.814  14.077  -7.295  0.40  8.94           O
ATOM   1933  O  AHOH S  19       6.566   3.124   9.199  0.50  8.59           O
ATOM   1934  O  BHOH S  19       6.847   1.889  10.166  0.50 16.33           O
ATOM   1935  O  AHOH S  20      -6.952  10.420   1.410  0.49 10.37           O
ATOM   1936  O  BHOH S  20      -6.282  10.227   0.982  0.49 11.40           O
ATOM   1937  O  AHOH S  21      -3.932   5.354  11.557  0.49 19.24           O
ATOM   1938  O  BHOH S  21      -4.563   5.276  11.091  0.49 13.07           O
ATOM   1939  O  AHOH S  22      17.054   7.644  -1.462  0.50 20.06           O
ATOM   1940  O  BHOH S  22      16.778   7.594  -2.157  0.50 12.47           O
ATOM   1941  O  AHOH S  23      -5.558  -3.221  11.955  0.46 12.35           O
ATOM   1942  O  BHOH S  23      -5.528  -4.007  11.656  0.46 13.35           O
ATOM   1943  O  AHOH S  24      -0.943  10.400   8.822  0.47 11.01           O
ATOM   1944  O  BHOH S  24      -0.382   9.983   9.358  0.47 13.45           O
ATOM   1945  O  AHOH S  25      21.110   1.327   5.061  0.50 13.49           O
ATOM   1946  O  BHOH S  25      20.948   1.626   5.711  0.50 19.10           O
ATOM   1947  O  AHOH S  26       5.929  -5.743   2.774  0.41 15.04           O
ATOM   1948  O  BHOH S  26       5.410  -5.716   2.399  0.41 15.43           O
ATOM   1949  O  AHOH S  27      13.430  -3.924  -3.820  0.34  8.97           O
ATOM   1950  O  BHOH S  27      12.994  -4.507  -3.939  0.34 13.32           O
ATOM   1951  O  AHOH S  28      -1.145  13.462  -8.649  0.42 10.81           O
ATOM   1952  O  BHOH S  28      -1.973  13.098  -9.032  0.42 15.02           O
ATOM   1953  O  AHOH S  29       3.441  -2.156   2.612  0.28 10.42           O
ATOM   1954  O  BHOH S  29       3.560  -2.814   2.445  0.28  7.47           O
ATOM   1955  O  AHOH S  30      -2.875  -6.229   3.399  0.50 16.18           O
ATOM   1956  O  BHOH S  30      -3.687  -6.332   3.764  0.50 14.04           O
ATOM   1957  O  AHOH S  31      13.666   5.596  10.420  0.50 19.49           O
ATOM   1958  O  BHOH S  31      13.295   4.988  10.016  0.50 14.78           O
ATOM   1959  O  AHOH S  32      17.686   6.244   6.006  0.50 12.13           O
ATOM   1960  O  BHOH S  32      16.820   5.712   4.032  0.50 17.41           O
ATOM   1961  O  AHOH S  33      19.202   0.696  -2.412  0.50 10.49           O
ATOM   1962  O  BHOH S  33      18.482  -0.211  -2.859  0.50 13.33           O
ATOM   1963  O  AHOH S  34      -6.594  -1.894  14.184  0.43 20.59           O
ATOM   1964  O  BHOH S  34      -7.096  -2.104  13.537  0.43 16.13           O
ATOM   1965  O  AHOH S  35      -0.646   5.897  12.407  0.50 16.06           O
ATOM   1966  O  BHOH S  35       0.145   5.440  12.382  0.50 16.76           O
ATOM   1967  O  AHOH S  36      10.428  -2.129   7.605  0.50 13.09           O
ATOM   1968  O  BHOH S  36       7.085  -7.471   4.583  0.50 12.64           O
ATOM   1969  O  AHOH S  37      -6.249  10.835  -5.645  0.47 20.20           O
ATOM   1970  O  BHOH S  37      -6.573   9.731  -5.831  0.47 12.13           O
ATOM   1971  O  AHOH S  38      21.111   2.123  -1.969  0.41 18.97           O
ATOM   1972  O  BHOH S  38      20.407   2.414  -2.702  0.41 10.58           O
ATOM   1973  O  AHOH S  39      16.826  -3.846  -2.644  0.50 25.88           O
ATOM   1974  O  BHOH S  39      15.968  -4.111  -3.321  0.50 14.54           O
ATOM   1975  O  AHOH S  40       3.331  15.787  -2.406  0.50 13.69           O
ATOM   1976  O  BHOH S  40       3.234  15.444  -3.450  0.50 16.45           O
ATOM   1977  O  AHOH S  41      16.580  -1.117   4.154  0.39 12.32           O
ATOM   1978  O  BHOH S  41      16.259  -1.796   4.724  0.39 13.13           O
ATOM   1979  O  AHOH S  42      -7.537   9.945   5.754  0.50 18.16           O
ATOM   1980  O  BHOH S  42      -7.712  10.203   4.954  0.50 18.28           O
ATOM   1981  O  AHOH S  43      16.757  -3.929   1.971  0.50 15.23           O
ATOM   1982  O  BHOH S  43      16.523  -3.738   0.975  0.50 19.18           O
ATOM   1983  O  AHOH S  44       3.457  -6.367   4.799  0.40 17.54           O
ATOM   1984  O  BHOH S  44       3.444  -5.761   4.052  0.40 11.14           O
ATOM   1985  O  AHOH S  45      22.052   4.064   5.335  0.42 13.57           O
ATOM   1986  O  BHOH S  45      22.328   4.182   4.348  0.42 18.70           O
ATOM   1987  O  AHOH S  46      -5.761   7.861   9.287  0.35 10.41           O
ATOM   1988  O  BHOH S  46      -5.442   8.809   9.567  0.35 16.28           O
ATOM   1989  O  AHOH S  47      -2.826  -0.168  -9.056  0.50 25.99           O
ATOM   1990  O  BHOH S  47      -2.083   0.387  -9.807  0.50 14.63           O
ATOM   1991  O  AHOH S  48       5.306  -8.063   6.557  0.50 14.06           O
ATOM   1992  O  BHOH S  48       4.044  -7.601   6.358  0.50 18.84           O
ATOM   1993  O  AHOH S  49       0.907  14.962  -3.509  0.50 13.01           O
ATOM   1994  O  BHOH S  49       1.115  14.809  -4.830  0.50 20.12           O
ATOM   1995  O  AHOH S  50       1.024  14.272  -5.863  0.39 19.04           O
ATOM   1996  O  BHOH S  50       1.040  13.967  -7.007  0.39 11.50           O
ATOM   1997  O  AHOH S  51      10.271  17.085  -6.371  0.38 14.93           O
ATOM   1998  O  BHOH S  51      10.564  17.022  -7.343  0.38 14.73           O
ATOM   1999  O  AHOH S  52       3.098   8.482   9.520  0.49 24.45           O
ATOM   2000  O  BHOH S  52       3.662   9.493   9.281  0.49 18.40           O
ATOM   2001  O  AHOH S  53       7.351  17.753  -1.047  0.37 16.60           O
ATOM   2002  O  BHOH S  53       6.600  17.765  -1.734  0.37 15.29           O
ATOM   2003  O  AHOH S  54      -7.343   5.802  -7.702  0.41 15.98           O
ATOM   2004  O  BHOH S  54      -7.694   5.244  -7.022  0.41 17.20           O
ATOM   2005  O  AHOH S  55      -5.978  12.426  -8.136  0.34 19.07           O
ATOM   2006  O  BHOH S  55      -6.298  11.879  -7.340  0.34 14.73           O
ATOM   2007  O  AHOH S  56      -5.159   1.734  15.833  0.37 17.79           O
ATOM   2008  O  BHOH S  56      -6.113   1.824  15.361  0.37 14.51           O
ATOM   2009  O  AHOH S  57       8.255  11.891   5.525  0.41 16.76           O
ATOM   2010  O  BHOH S  57       6.895  12.497   5.546  0.41 14.68           O
ATOM   2011  O  AHOH S  58      -8.783   0.612   6.866  0.50 19.61           O
ATOM   2012  O  BHOH S  58      -7.642  -1.064   7.393  0.50 18.98           O
ATOM   2013  O  AHOH S  59      -0.277  -0.249  14.329  0.32 12.45           O
ATOM   2014  O  BHOH S  59      -0.980   0.359  14.764  0.32 13.44           O
ATOM   2015  O  AHOH S  60       8.785   6.649   7.443  0.50 17.83           O
ATOM   2016  O  BHOH S  60       8.511   5.634   7.998  0.50 23.77           O
ATOM   2017  O  AHOH S  61      15.558   4.894   5.938  0.33 12.93           O
ATOM   2018  O  BHOH S  61      14.625   4.601   6.549  0.33 13.81           O
ATOM   2019  O  AHOH S  62      18.530   4.198   8.098  0.50 20.17           O
ATOM   2020  O  BHOH S  62      17.498   4.747   7.844  0.50 17.70           O
ATOM   2021  O  AHOH S  63       1.906  -3.682  -7.093  0.27 11.03           O
ATOM   2022  O  BHOH S  63       1.150  -3.379  -7.482  0.27 14.91           O
ATOM   2023  O  AHOH S  64      17.390   7.775   1.597  0.50 23.89           O
ATOM   2024  O  BHOH S  64      18.073   7.453   2.465  0.50 20.02           O
ATOM   2025  O  AHOH S  65       6.675  -6.244  -4.132  0.45 22.52           O
ATOM   2026  O  BHOH S  65       6.739  -5.759  -4.946  0.45 17.24           O
ATOM   2027  O  AHOH S  66      -0.058  16.274  -1.120  0.50 19.87           O
ATOM   2028  O  BHOH S  66       1.583  16.571  -1.623  0.50 23.11           O
ATOM   2029  O  AHOH S  67     -11.144   2.432   0.720  0.32 17.68           O
ATOM   2030  O  BHOH S  67     -11.744   2.638  -0.274  0.32 16.29           O
ATOM   2031  O  AHOH S  68      -9.171  -3.709  -2.415  0.50 19.99           O
ATOM   2032  O  BHOH S  68      -8.570  -4.592  -1.612  0.50 23.18           O
ATOM   2033  O  AHOH S  69      -4.751   9.681   9.804  0.24  9.62           O
ATOM   2034  O  BHOH S  69      -4.185  10.491   9.939  0.24  8.94           O
ATOM   2035  O  AHOH S  70      -4.435  -1.424  16.674  0.30 13.84           O
ATOM   2036  O  BHOH S  70      -3.306  -0.718  16.683  0.30 18.48           O
ATOM   2037  O  AHOH S  71      13.715   5.270   7.077  0.34 15.51           O
ATOM   2038  O  BHOH S  71      12.875   5.487   7.792  0.34 17.24           O
ATOM   2039  O  AHOH S  72       4.335   4.832   9.903  0.30 20.30           O
ATOM   2040  O  BHOH S  72       3.789   3.693   9.887  0.30 14.06           O
ATOM   2041  O  AHOH S  73     -10.657   2.063   5.670  0.35 13.62           O
ATOM   2042  O  BHOH S  73     -11.195   1.665   4.782  0.35 18.49           O
ATOM   2043  O  AHOH S  74      -5.127  -5.746  10.948  0.23 13.11           O
ATOM   2044  O  BHOH S  74      -4.285  -6.395  10.889  0.23 12.18           O
ATOM   2045  O  AHOH S  75      11.666  15.525   3.538  0.50 26.07           O
ATOM   2046  O  BHOH S  75      12.081  16.846   3.478  0.50 19.75           O
ATOM   2047  O  AHOH S  76       4.264   4.263 -10.449  0.50 20.58           O
ATOM   2048  O  BHOH S  76       5.157   3.561 -10.273  0.50 22.27           O
ATOM   2049  O  AHOH S  77      -9.860  -2.334   0.302  0.50 23.17           O
ATOM   2050  O  BHOH S  77      -9.858  -0.806   1.491  0.50 20.44           O"""




for line in waters.split("\n"):
    old_res_num = int(line.split()[5])
    new_res_num = str(64 + old_res_num)
    new_res_num = " "*(3-len(new_res_num))+new_res_num
    new_line = line[:23] + new_res_num + line[26:]
    assert len(new_line) == len(line), len(new_res_num)
    print(new_line)
# %%

eff_constraints = """bond {
      action = *add
      atom_selection_1 = "name O and resseq 1 and chain S and altid A"
      atom_selection_2 = "name O and resseq 1 and chain S and altid B"
      distance_ideal = 0.4436
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 2 and chain S and altid A"
      atom_selection_2 = "name O and resseq 2 and chain S and altid B"
      distance_ideal = 0.3167
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 3 and chain S and altid A"
      atom_selection_2 = "name O and resseq 3 and chain S and altid B"
      distance_ideal = 0.4769
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 4 and chain S and altid A"
      atom_selection_2 = "name O and resseq 4 and chain S and altid B"
      distance_ideal = 0.3964
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 5 and chain S and altid A"
      atom_selection_2 = "name O and resseq 5 and chain S and altid B"
      distance_ideal = 0.5163
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 6 and chain S and altid A"
      atom_selection_2 = "name O and resseq 6 and chain S and altid B"
      distance_ideal = 0.7324
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 7 and chain S and altid A"
      atom_selection_2 = "name O and resseq 7 and chain S and altid B"
      distance_ideal = 0.7803
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 8 and chain S and altid A"
      atom_selection_2 = "name O and resseq 8 and chain S and altid B"
      distance_ideal = 1.3207
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 9 and chain S and altid A"
      atom_selection_2 = "name O and resseq 9 and chain S and altid B"
      distance_ideal = 0.5708
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 10 and chain S and altid A"
      atom_selection_2 = "name O and resseq 10 and chain S and altid B"
      distance_ideal = 0.6301
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 11 and chain S and altid A"
      atom_selection_2 = "name O and resseq 11 and chain S and altid B"
      distance_ideal = 0.7798
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 12 and chain S and altid A"
      atom_selection_2 = "name O and resseq 12 and chain S and altid B"
      distance_ideal = 0.7461
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 13 and chain S and altid A"
      atom_selection_2 = "name O and resseq 13 and chain S and altid B"
      distance_ideal = 0.7094
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 14 and chain S and altid A"
      atom_selection_2 = "name O and resseq 14 and chain S and altid B"
      distance_ideal = 1.1348
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 15 and chain S and altid A"
      atom_selection_2 = "name O and resseq 15 and chain S and altid B"
      distance_ideal = 0.8058
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 16 and chain S and altid A"
      atom_selection_2 = "name O and resseq 16 and chain S and altid B"
      distance_ideal = 0.7146
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 17 and chain S and altid A"
      atom_selection_2 = "name O and resseq 17 and chain S and altid B"
      distance_ideal = 0.6214
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 18 and chain S and altid A"
      atom_selection_2 = "name O and resseq 18 and chain S and altid B"
      distance_ideal = 0.7285
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 19 and chain S and altid A"
      atom_selection_2 = "name O and resseq 19 and chain S and altid B"
      distance_ideal = 1.5855
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 20 and chain S and altid A"
      atom_selection_2 = "name O and resseq 20 and chain S and altid B"
      distance_ideal = 0.8156
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 21 and chain S and altid A"
      atom_selection_2 = "name O and resseq 21 and chain S and altid B"
      distance_ideal = 0.7867
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 22 and chain S and altid A"
      atom_selection_2 = "name O and resseq 22 and chain S and altid B"
      distance_ideal = 0.8393
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 23 and chain S and altid A"
      atom_selection_2 = "name O and resseq 23 and chain S and altid B"
      distance_ideal = 0.8325
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 24 and chain S and altid A"
      atom_selection_2 = "name O and resseq 24 and chain S and altid B"
      distance_ideal = 0.8886
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 25 and chain S and altid A"
      atom_selection_2 = "name O and resseq 25 and chain S and altid B"
      distance_ideal = 0.7608
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 26 and chain S and altid A"
      atom_selection_2 = "name O and resseq 26 and chain S and altid B"
      distance_ideal = 0.6139
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 27 and chain S and altid A"
      atom_selection_2 = "name O and resseq 27 and chain S and altid B"
      distance_ideal = 0.7522
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 28 and chain S and altid A"
      atom_selection_2 = "name O and resseq 28 and chain S and altid B"
      distance_ideal = 0.9834
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 29 and chain S and altid A"
      atom_selection_2 = "name O and resseq 29 and chain S and altid B"
      distance_ideal = 0.7141
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 30 and chain S and altid A"
      atom_selection_2 = "name O and resseq 30 and chain S and altid B"
      distance_ideal = 0.9085
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 31 and chain S and altid A"
      atom_selection_2 = "name O and resseq 31 and chain S and altid B"
      distance_ideal = 0.8265
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 32 and chain S and altid A"
      atom_selection_2 = "name O and resseq 32 and chain S and altid B"
      distance_ideal = 2.2270
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 33 and chain S and altid A"
      atom_selection_2 = "name O and resseq 33 and chain S and altid B"
      distance_ideal = 1.2347
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 34 and chain S and altid A"
      atom_selection_2 = "name O and resseq 34 and chain S and altid B"
      distance_ideal = 0.8067
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 35 and chain S and altid A"
      atom_selection_2 = "name O and resseq 35 and chain S and altid B"
      distance_ideal = 0.9120
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 36 and chain S and altid A"
      atom_selection_2 = "name O and resseq 36 and chain S and altid B"
      distance_ideal = 6.9923
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 37 and chain S and altid A"
      atom_selection_2 = "name O and resseq 37 and chain S and altid B"
      distance_ideal = 1.1617
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 38 and chain S and altid A"
      atom_selection_2 = "name O and resseq 38 and chain S and altid B"
      distance_ideal = 1.0718
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 39 and chain S and altid A"
      atom_selection_2 = "name O and resseq 39 and chain S and altid B"
      distance_ideal = 1.1046
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 40 and chain S and altid A"
      atom_selection_2 = "name O and resseq 40 and chain S and altid B"
      distance_ideal = 1.1078
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 41 and chain S and altid A"
      atom_selection_2 = "name O and resseq 41 and chain S and altid B"
      distance_ideal = 0.9540
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 42 and chain S and altid A"
      atom_selection_2 = "name O and resseq 42 and chain S and altid B"
      distance_ideal = 0.8325
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 43 and chain S and altid A"
      atom_selection_2 = "name O and resseq 43 and chain S and altid B"
      distance_ideal = 1.0343
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 44 and chain S and altid A"
      atom_selection_2 = "name O and resseq 44 and chain S and altid B"
      distance_ideal = 0.9965
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 45 and chain S and altid A"
      atom_selection_2 = "name O and resseq 45 and chain S and altid B"
      distance_ideal = 1.0373
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 46 and chain S and altid A"
      atom_selection_2 = "name O and resseq 46 and chain S and altid B"
      distance_ideal = 1.0454
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 47 and chain S and altid A"
      atom_selection_2 = "name O and resseq 47 and chain S and altid B"
      distance_ideal = 1.2136
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 48 and chain S and altid A"
      atom_selection_2 = "name O and resseq 48 and chain S and altid B"
      distance_ideal = 1.3593
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 49 and chain S and altid A"
      atom_selection_2 = "name O and resseq 49 and chain S and altid B"
      distance_ideal = 1.3604
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 50 and chain S and altid A"
      atom_selection_2 = "name O and resseq 50 and chain S and altid B"
      distance_ideal = 1.1944
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 51 and chain S and altid A"
      atom_selection_2 = "name O and resseq 51 and chain S and altid B"
      distance_ideal = 1.0083
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 52 and chain S and altid A"
      atom_selection_2 = "name O and resseq 52 and chain S and altid B"
      distance_ideal = 1.1515
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 53 and chain S and altid A"
      atom_selection_2 = "name O and resseq 53 and chain S and altid B"
      distance_ideal = 1.0003
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 54 and chain S and altid A"
      atom_selection_2 = "name O and resseq 54 and chain S and altid B"
      distance_ideal = 0.8963
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 55 and chain S and altid A"
      atom_selection_2 = "name O and resseq 55 and chain S and altid B"
      distance_ideal = 1.0567
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 56 and chain S and altid A"
      atom_selection_2 = "name O and resseq 56 and chain S and altid B"
      distance_ideal = 1.0341
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 57 and chain S and altid A"
      atom_selection_2 = "name O and resseq 57 and chain S and altid B"
      distance_ideal = 1.4880
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 58 and chain S and altid A"
      atom_selection_2 = "name O and resseq 58 and chain S and altid B"
      distance_ideal = 2.0742
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 59 and chain S and altid A"
      atom_selection_2 = "name O and resseq 59 and chain S and altid B"
      distance_ideal = 1.0331
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 60 and chain S and altid A"
      atom_selection_2 = "name O and resseq 60 and chain S and altid B"
      distance_ideal = 1.2260
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 61 and chain S and altid A"
      atom_selection_2 = "name O and resseq 61 and chain S and altid B"
      distance_ideal = 1.1611
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 62 and chain S and altid A"
      atom_selection_2 = "name O and resseq 62 and chain S and altid B"
      distance_ideal = 1.1641
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 63 and chain S and altid A"
      atom_selection_2 = "name O and resseq 63 and chain S and altid B"
      distance_ideal = 0.9210
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 64 and chain S and altid A"
      atom_selection_2 = "name O and resseq 64 and chain S and altid B"
      distance_ideal = 1.0943
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 65 and chain S and altid A"
      atom_selection_2 = "name O and resseq 65 and chain S and altid B"
      distance_ideal = 0.8738
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 66 and chain S and altid A"
      atom_selection_2 = "name O and resseq 66 and chain S and altid B"
      distance_ideal = 1.7408
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 67 and chain S and altid A"
      atom_selection_2 = "name O and resseq 67 and chain S and altid B"
      distance_ideal = 1.1932
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 68 and chain S and altid A"
      atom_selection_2 = "name O and resseq 68 and chain S and altid B"
      distance_ideal = 1.3023
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 69 and chain S and altid A"
      atom_selection_2 = "name O and resseq 69 and chain S and altid B"
      distance_ideal = 1.0004
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 70 and chain S and altid A"
      atom_selection_2 = "name O and resseq 70 and chain S and altid B"
      distance_ideal = 1.2823
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 71 and chain S and altid A"
      atom_selection_2 = "name O and resseq 71 and chain S and altid B"
      distance_ideal = 1.1208
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 72 and chain S and altid A"
      atom_selection_2 = "name O and resseq 72 and chain S and altid B"
      distance_ideal = 1.1993
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 73 and chain S and altid A"
      atom_selection_2 = "name O and resseq 73 and chain S and altid B"
      distance_ideal = 1.0859
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 74 and chain S and altid A"
      atom_selection_2 = "name O and resseq 74 and chain S and altid B"
      distance_ideal = 1.0914
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 75 and chain S and altid A"
      atom_selection_2 = "name O and resseq 75 and chain S and altid B"
      distance_ideal = 1.3890
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 76 and chain S and altid A"
      atom_selection_2 = "name O and resseq 76 and chain S and altid B"
      distance_ideal = 1.1576
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 77 and chain S and altid A"
      atom_selection_2 = "name O and resseq 77 and chain S and altid B"
      distance_ideal = 1.9518
      sigma = 0.2
    }"""

for line in eff_constraints.split("\n"):
    if not "atom_selection" in line:
        print(line)
        continue 
    old_res_num = int(line.split()[6])
    new_res_num = str(64 + old_res_num)
    new_line = " "*6+ ' '.join(line.split()[:6]) + " " + new_res_num + " " + ' '.join(line.split()[7:])
    print(new_line)
# %%

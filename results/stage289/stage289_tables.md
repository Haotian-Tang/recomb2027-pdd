# Stage 289 comparison tables

Cases analysed: 1188; cases with incomplete run records: 0.

## 1. Completion and runtime (median child CPU seconds, 2 repetitions per case)

| family              | method                | cases | complete | timeout | error | median CPU s | total CPU s |
|---------------------|-----------------------|-------|----------|---------|-------|--------------|-------------|
| exact_10k           | exact_stage43_partial | 296   | 200      | 96      | 0     | 0.023        | 11536.1182  |
| exact_10k           | dd_default            | 296   | 296      | 0       | 0     | 1.554        | 757.6926    |
| exact_10k           | dd_true_k             | 296   | 296      | 0       | 0     | 1.563        | 858.7838    |
| exact_10k           | dd_resolution         | 296   | 296      | 0       | 0     | 1.568        | 822.6662    |
| exact_biotrees      | exact_stage43_partial | 596   | 475      | 121     | 0     | 0.050        | 15993.9890  |
| exact_biotrees      | dd_default            | 596   | 596      | 0       | 0     | 3.138        | 2888.2689   |
| exact_biotrees      | dd_true_k             | 596   | 596      | 0       | 0     | 3.224        | 2951.9050   |
| exact_biotrees      | dd_resolution         | 596   | 596      | 0       | 0     | 3.230        | 2977.0703   |
| krepp_estimated_10k | exact_stage43_partial | 296   | 293      | 3       | 0     | 0.009        | 363.3036    |
| krepp_estimated_10k | dd_default            | 296   | 290      | 0       | 6     | 3.849        | 1551.3806   |
| krepp_estimated_10k | dd_true_k             | 296   | 285      | 0       | 11    | 1.634        | 809.2068    |
| krepp_estimated_10k | dd_resolution         | 296   | 289      | 0       | 7     | 2.424        | 986.0584    |

## 2. Source count returned, against the generating count

| family              | arm           | complete | k = k_true | k < k_true | k > k_true | stop reason                                                         |
|---------------------|---------------|----------|------------|------------|------------|---------------------------------------------------------------------|
| exact_10k           | dd_default    | 296      | 198        | 92         | 6          | {"low_abundance": 72, "near_zero_residue": 224}                     |
| exact_10k           | dd_true_k     | 296      | 296        | 0          | 0          | {"k_exhausted": 296}                                                |
| exact_10k           | dd_resolution | 296      | 235        | 45         | 16         | {"near_zero_residue": 294, "low_abundance": 2}                      |
| exact_biotrees      | dd_default    | 596      | 450        | 143        | 3          | {"near_zero_residue": 456, "low_abundance": 140}                    |
| exact_biotrees      | dd_true_k     | 596      | 596        | 0          | 0          | {"k_exhausted": 596}                                                |
| exact_biotrees      | dd_resolution | 596      | 570        | 12         | 14         | {"near_zero_residue": 596}                                          |
| krepp_estimated_10k | dd_default    | 290      | 17         | 122        | 151        | {"low_abundance": 276, "near_zero_residue": 14}                     |
| krepp_estimated_10k | dd_true_k     | 285      | 285        | 0          | 0          | {"k_exhausted": 285}                                                |
| krepp_estimated_10k | dd_resolution | 289      | 17         | 96         | 176        | {"low_abundance": 125, "k_exhausted": 150, "near_zero_residue": 14} |

## 3. Support recovery and exact feasibility of the returned support

| family              | method                | complete | support = truth | true edges recovered | support exactly feasible | support provably infeasible | other | median W1 | median loss |
|---------------------|-----------------------|----------|-----------------|----------------------|--------------------------|-----------------------------|-------|-----------|-------------|
| exact_10k           | dd_default            | 296      | 190             | 0.874                | 210                      | 83                          | 3     | 2.01e-05  | 3.56e-17    |
| exact_10k           | dd_true_k             | 296      | 228             | 0.943                | 256                      | 29                          | 11    | 1.66e-07  | 2.65e-17    |
| exact_10k           | dd_resolution         | 296      | 210             | 0.941                | 244                      | 42                          | 10    | 4.81e-08  | 3.41e-17    |
| exact_10k           | exact_stage43_partial | 200      | 196             | -                    | 200                      | 0                           | -     | -         | -           |
| exact_biotrees      | dd_default            | 596      | 433             | 0.917                | 446                      | 147                         | 3     | 1.54e-08  | 6.95e-17    |
| exact_biotrees      | dd_true_k             | 596      | 555             | 0.983                | 569                      | 19                          | 8     | 1.06e-08  | 5.32e-17    |
| exact_biotrees      | dd_resolution         | 596      | 547             | 0.986                | 570                      | 17                          | 9     | 1.06e-08  | 5.54e-17    |
| exact_biotrees      | exact_stage43_partial | 475      | 473             | -                    | 475                      | 0                           | -     | -         | -           |
| krepp_estimated_10k | dd_default            | 290      | 0               | 0.238                | 0                        | 284                         | 6     | 0.084     | 0.006       |
| krepp_estimated_10k | dd_true_k             | 285      | 0               | 0.200                | 0                        | 285                         | 0     | 0.076     | 0.006       |
| krepp_estimated_10k | dd_resolution         | 289      | 0               | 0.230                | 0                        | 287                         | 2     | 0.081     | 0.006       |
| krepp_estimated_10k | exact_stage43_partial | 293      | 0               | -                    | 0                        | 293                         | -     | -         | -           |

## 4. By generating source count

### exact_10k

| k_true | cases | exact complete | exact median CPU s | dd_default support=truth | dd_true_k support=truth | dd_resolution support=truth | dd_default feasible | dd_true_k feasible | dd_resolution feasible | dd_default median CPU s | dd_true_k median CPU s | dd_resolution median CPU s |
|--------|-------|----------------|--------------------|--------------------------|-------------------------|-----------------------------|---------------------|--------------------|------------------------|-------------------------|------------------------|----------------------------|
| 3      | 100   | 100            | 0.015              | 95                       | 97                      | 95                          | 98                  | 99                 | 98                     | 0.769                   | 0.764                  | 0.766                      |
| 5      | 100   | 100            | 0.047              | 86                       | 91                      | 85                          | 87                  | 91                 | 86                     | 1.458                   | 1.461                  | 1.460                      |
| 10     | 96    | 0              | 119.997            | 9                        | 40                      | 30                          | 25                  | 66                 | 60                     | 4.216                   | 5.397                  | 5.074                      |

### exact_biotrees

| k_true | cases | exact complete | exact median CPU s | dd_default support=truth | dd_true_k support=truth | dd_resolution support=truth | dd_default feasible | dd_true_k feasible | dd_resolution feasible | dd_default median CPU s | dd_true_k median CPU s | dd_resolution median CPU s |
|--------|-------|----------------|--------------------|--------------------------|-------------------------|-----------------------------|---------------------|--------------------|------------------------|-------------------------|------------------------|----------------------------|
| 2      | 120   | 120            | 0.024              | 119                      | 119                     | 119                         | 119                 | 119                | 119                    | 1.171                   | 1.174                  | 1.169                      |
| 3      | 120   | 120            | 0.029              | 117                      | 118                     | 118                         | 117                 | 118                | 118                    | 1.754                   | 1.760                  | 1.754                      |
| 5      | 120   | 120            | 0.060              | 92                       | 117                     | 114                         | 94                  | 118                | 116                    | 2.972                   | 2.981                  | 2.972                      |
| 7      | 119   | 109            | 3.102              | 75                       | 111                     | 108                         | 78                  | 113                | 113                    | 4.573                   | 4.578                  | 4.600                      |
| 10     | 117   | 6              | 119.995            | 30                       | 90                      | 88                          | 38                  | 101                | 104                    | 7.751                   | 7.969                  | 8.178                      |

### krepp_estimated_10k

| k_true | cases | exact complete | exact median CPU s | dd_default support=truth | dd_true_k support=truth | dd_resolution support=truth | dd_default feasible | dd_true_k feasible | dd_resolution feasible | dd_default median CPU s | dd_true_k median CPU s | dd_resolution median CPU s |
|--------|-------|----------------|--------------------|--------------------------|-------------------------|-----------------------------|---------------------|--------------------|------------------------|-------------------------|------------------------|----------------------------|
| 3      | 100   | 100            | 0.009              | 0                        | 0                       | 0                           | 0                   | 0                  | 0                      | 5.318                   | 0.754                  | 1.913                      |
| 5      | 100   | 100            | 0.009              | 0                        | 0                       | 0                           | 0                   | 0                  | 0                      | 3.484                   | 1.609                  | 2.376                      |
| 10     | 96    | 93             | 0.009              | 0                        | 0                       | 0                           | 0                   | 0                  | 0                      | 2.954                   | 5.467                  | 4.066                      |

## 5. By reference-set size

### exact_10k

| leaves | cases | exact complete | exact median CPU s | dd_default support=truth | dd_true_k support=truth | dd_resolution support=truth | dd_default median CPU s | dd_true_k median CPU s | dd_resolution median CPU s |
|--------|-------|----------------|--------------------|--------------------------|-------------------------|-----------------------------|-------------------------|------------------------|----------------------------|
| <=50   | 156   | 100            | 0.097              | 89                       | 108                     | 94                          | 1.190                   | 1.076                  | 1.200                      |
| 51-150 | 140   | 100            | 0.039              | 101                      | 120                     | 116                         | 1.931                   | 1.931                  | 1.930                      |

### exact_biotrees

| leaves  | cases | exact complete | exact median CPU s | dd_default support=truth | dd_true_k support=truth | dd_resolution support=truth | dd_default median CPU s | dd_true_k median CPU s | dd_resolution median CPU s |
|---------|-------|----------------|--------------------|--------------------------|-------------------------|-----------------------------|-------------------------|------------------------|----------------------------|
| <=50    | 147   | 114            | 0.058              | 94                       | 116                     | 113                         | 0.979                   | 0.979                  | 0.977                      |
| 51-150  | 249   | 199            | 0.049              | 191                      | 241                     | 241                         | 2.709                   | 2.729                  | 2.725                      |
| 151-500 | 200   | 162            | 0.081              | 148                      | 198                     | 193                         | 6.930                   | 7.037                  | 7.028                      |

### krepp_estimated_10k

| leaves | cases | exact complete | exact median CPU s | dd_default support=truth | dd_true_k support=truth | dd_resolution support=truth | dd_default median CPU s | dd_true_k median CPU s | dd_resolution median CPU s |
|--------|-------|----------------|--------------------|--------------------------|-------------------------|-----------------------------|-------------------------|------------------------|----------------------------|
| <=50   | 168   | 165            | 0.008              | 0                        | 0                       | 0                           | 4.224                   | 1.449                  | 1.935                      |
| 51-150 | 128   | 128            | 0.011              | 0                        | 0                       | 0                           | 3.223                   | 1.739                  | 3.037                      |

## 6. Minimality on exactly explainable data

An explanation is counted only when the returned support is certified exactly feasible. The reference minimum is the exact solver's answer where it completed, otherwise the generating count, which is an upper bound on the true minimum.

| family         | arm           | cases | exact and minimal | exact but non-minimal | below reference minimum | no exact explanation | incomplete |
|----------------|---------------|-------|-------------------|-----------------------|-------------------------|----------------------|------------|
| exact_10k      | dd_default    | 296   | 193               | 4                     | 13                      | 86                   | 0          |
| exact_10k      | dd_true_k     | 296   | 254               | 2                     | 0                       | 40                   | 0          |
| exact_10k      | dd_resolution | 296   | 216               | 4                     | 24                      | 52                   | 0          |
| exact_biotrees | dd_default    | 596   | 439               | 2                     | 5                       | 150                  | 0          |
| exact_biotrees | dd_true_k     | 596   | 569               | 0                     | 0                       | 27                   | 0          |
| exact_biotrees | dd_resolution | 596   | 556               | 6                     | 8                       | 26                   | 0          |

## 7. Identifiability of the returned support

| family              | arm           | complete | feasible but non-identifiable | placement at a vertex | repeated edge |
|---------------------|---------------|----------|-------------------------------|-----------------------|---------------|
| exact_10k           | dd_default    | 296      | 42                            | 17                    | 0             |
| exact_10k           | dd_true_k     | 296      | 75                            | 19                    | 0             |
| exact_10k           | dd_resolution | 296      | 62                            | 18                    | 0             |
| exact_biotrees      | dd_default    | 596      | 46                            | 20                    | 0             |
| exact_biotrees      | dd_true_k     | 596      | 69                            | 25                    | 0             |
| exact_biotrees      | dd_resolution | 596      | 75                            | 27                    | 0             |
| krepp_estimated_10k | dd_default    | 290      | 0                             | 258                   | 0             |
| krepp_estimated_10k | dd_true_k     | 285      | 0                             | 223                   | 0             |
| krepp_estimated_10k | dd_resolution | 289      | 0                             | 256                   | 0             |


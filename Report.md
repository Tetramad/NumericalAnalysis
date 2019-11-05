# LU 분해 구현 및 분석
## 구현 목록
|                 |             |         |
| :-------------: | :---------: | :-----: |
| scipy.linalg.lu | DoolittleLU | CroutLU |
|       LU        |    LUPP     |  LUCP   |
|      LULT       |  Cholesky   |         |

## 시간 분석
### 테스트 환경
|              |                       |
| :----------: | --------------------: |
|    Python    |                 3.7.4 |
|    scipy     |                 1.3.1 |
|    numpy     |                1.17.1 |
| Architecture |                 amd64 |
|     CPU      | Intel i5-4570 3.20GHz |
|     RAM      |                  12GB |
|      OS      |  Windows 10 Home 1903 |

$
\displaystyle \\
\text{평균 실행 시간} = \frac{\text{1,000회 실행 시간}}{1,000}\ (\mu s) \\
\text{PASS 조건:}\quad\forall (i,j),\quad \left|a_{i,j} - b_{i,j}\right| < 10^{-8} \\
$

### 테스트 결과
|                 | Low rank Vandermonde | Non-square(2, 3) | Non-square(3, 2) |
| :-------------: | -------------------: | ---------------: | ---------------: |
| scipy.linalg.lu |                   34 |               28 |               32 |
|   DoolittleLU   |                   36 |               20 |               26 |
|     CroutLU     |                   35 |               27 |               21 |
|       LU        |                   18 |               10 |               16 |
|      LUPP       |                   48 |               37 |               45 |
|      LUCP       |                   61 |               52 |               56 |

|                 | Single element | Symmetric | Positive definite |
| :-------------: | -------------: | --------: | ----------------: |
| scipy.linalg.lu |              8 |        32 |                33 |
|   DoolittleLU   |              4 |        37 |                38 |
|     CroutLU     |              4 |        36 |                37 |
|       LU        |              5 |        18 |                18 |
|      LUPP       |             22 |        45 |                49 |
|      LUCP       |             33 |        61 |                61 |
|      LDLT       |             11 |        25 |                25 |
|    Cholesky     |              . |         . |                19 |

|                 | Bad condition of naive LU | Bad condition of LUPP | Big(50, 50) |
| :-------------: | ------------------------: | --------------------: | ----------: |
| scipy.linalg.lu |                        33 |                    33 |        5559 |
|   DoolittleLU   |                      FAIL |                    37 |       37020 |
|     CroutLU     |                      FAIL |                    38 |       36227 |
|       LU        |                      FAIL |                    20 |       49356 |
|      LUPP       |                        47 |                  FAIL |       49310 |
|      LUCP       |                        64 |                    61 |       49550 |

|                 | High rank(16) Vandermonde | High rank(17) Vandermonde | Singular |
| :-------------: | ------------------------: | ------------------------: | -------: |
| scipy.linalg.lu |                      FAIL |                      FAIL |       37 |
|   DoolittleLU   |                      1264 |                      FAIL |     FAIL |
|     CroutLU     |                      1259 |                      FAIL |     FAIL |
|       LU        |                      FAIL |                      FAIL |     FAIL |
|      LUPP       |                      FAIL |                      FAIL |     FAIL |
|      LUCP       |                      FAIL |                      FAIL |     FAIL |

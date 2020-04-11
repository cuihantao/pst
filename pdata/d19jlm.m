% 19 Machine Benchmark
% 

% bus data format
% bus: 
% col1 number
% col2 voltage magnitude(pu)
% col3 voltage angle(degree)
% col4 p_gen(pu)
% col5 q_gen(pu),
% col6 p_load(pu)
% col7 q_load(pu)
% col8 G shunt(pu)
% col9 B shunt(pu)
% col10 bus_type
%     bus_type - 1, swing bus
%         - 2, generator bus (PV bus)
%         - 3, load bus (PQ bus)
% col11 q_gen_max(pu)
% col12 q_gen_min(pu)
% col13 v_rated (kV)
% col14 v_max  pu
% col15 v_min  pu

bus = [...
   1 1.0490    0.0000  28.3794  23.2983   0.0000   0.0000   0.00   0.00  1  25.00 -10.00    20.00  1.1  0.9;
   2 1.0490   28.1220 121.0406  62.9955   0.0000   0.0000   0.00   0.00  2  73.20 -35.00    20.00  1.1  0.9;
   3 1.0490  -57.8509 110.1829  70.1339   0.0000   0.0000   0.00   0.00  2  72.00   0.00    20.00  1.1  0.9;
   4 1.0490  -98.8765  74.7150  55.1440   0.0000   0.0000   0.00   0.00  2  93.75 -10.00    20.00  1.1  0.9;
   5 1.0490  -63.8096  68.3083  75.4348   0.0000   0.0000   0.00   0.00  2  78.00   0.00    20.00  1.1  0.9;
   6 1.0490    8.9660 127.0900  52.5457   0.0000   0.0000   0.00   0.00  2 101.80 -10.00    20.00  1.1  0.9;
   7 1.0490   16.2350  19.6897   9.1951   0.0000   0.0000   0.00   0.00  2  25.00 -10.00    20.00  1.1  0.9;
   8 1.0490  -15.9550  22.6313  10.5317   0.0000   0.0000   0.00   0.00  2  25.00 -10.00    20.00  1.1  0.9;
   9 1.0300   22.6917  12.2000  -0.5245   0.0000   0.0000   0.00   0.00  2   7.00  -5.00    20.00  1.1  0.9;
  10 1.0400   -7.8713  16.8000   0.9749   0.0000   0.0000   0.00   0.00  2  16.00   0.00    20.00  1.1  0.9;
  11 1.0490  -56.2245  27.0400  12.2507   0.0000   0.0000   0.00   0.00  2  34.00 -10.00    20.00  1.1  0.9;
  12 1.0490  -58.4070   5.6000   4.0002   0.0000   0.0000   0.00   0.00  2  16.00 -10.00    20.00  1.1  0.9;
  13 1.0490   25.3923   6.0300   4.4545   0.0000   0.0000   0.00   0.00  2  10.00 -10.00    20.00  1.1  0.9;
  14 1.0490  -55.9158  12.7300   4.9559   0.0000   0.0000   0.00   0.00  2   6.00 -10.00    20.00  1.1  0.9;
  15 1.0490    5.7456   6.6800   5.0596   0.0000   0.0000   0.00   0.00  2   8.00 -10.00    20.00  1.1  0.9;
  16 1.0490    7.4830   9.4000   4.7339   0.0000   0.0000   0.00   0.00  2  16.00 -10.00    20.00  1.1  0.9;
  17 1.0490   -2.3803  11.4100  37.4331   0.0000   0.0000   0.00   0.00  2  40.00   0.00    20.00  1.1  0.9;
  18 1.0490  -86.8303  25.4000   7.0974   0.0000   0.0000   0.00   0.00  2   7.50 -10.00    20.00  1.1  0.9;
  19 1.0490  -87.8436  21.7700   6.1294   0.0000   0.0000   0.00   0.00  2   7.25 -10.00    20.00  1.1  0.9;
  20 1.0049  -62.0468   0.0000   0.0000 105.0300  60.0000   0.00   0.00  3   0.00   0.00   500.00  1.1  0.9;
  21 0.9957   21.4671   0.0000   0.0000 107.0400  48.0000   0.00   0.00  3   0.00   0.00   500.00  1.1  0.9;
  22 0.9928  -66.8173   0.0000   0.0000  97.1000  40.0000   0.00  10.00  3   0.00   0.00   500.00  1.1  0.9;
  23 1.0430   15.5132   0.0000   0.0000  11.6800   4.0000   0.00   0.00  3   0.00   0.00   500.00  1.1  0.9;
  24 1.0421  -16.7853   0.0000   0.0000  41.5300  15.0000   0.00  10.00  3   0.00   0.00   230.00  1.1  0.9;
  25 1.0098    2.7673   0.0000   0.0000  80.8300  45.0000   0.00   0.00  3   0.00   0.00   500.00  1.1  0.9;
  26 1.0060   -3.0830   0.0000   0.0000  80.8300  35.0000   0.00  15.00  3   0.00   0.00   500.00  1.1  0.9;
  27 1.0050 -105.2786   0.0000   0.0000  97.0600  40.0000   0.00  59.00  3   0.00   0.00   500.00  1.1  0.9;
  28 1.0086 -102.1152   0.0000   0.0000  87.0600  32.0000   0.00   0.00  3   0.00   0.00   500.00  1.1  0.9;
  29 1.0243  -59.0441   0.0000   0.0000   0.0000   0.0000   0.00   0.00  3   0.00   0.00    20.00  1.1  0.9;
  30 1.0226  -60.4863   0.0000   0.0000   0.0000   0.0000   0.00   0.00  3   0.00   0.00    20.00  1.1  0.9;
  31 1.0153   23.1674   0.0000   0.0000   0.0000   0.0000   0.00   0.00  3   0.00   0.00    20.00  1.1  0.9;
  32 1.0243  -59.6333   0.0000   0.0000   0.0000   0.0000   0.00   0.00  3   0.00   0.00    20.00  1.1  0.9;
  33 1.0235    3.8774   0.0000   0.0000   0.0000   0.0000   0.00   0.00  3   0.00   0.00    20.00  1.1  0.9;
  34 1.0215    4.3326   0.0000   0.0000   0.0000   0.0000   0.00   0.00  3   0.00   0.00    20.00  1.1  0.9;
  35 1.0310   -2.6429   0.0000   0.0000   0.0000   0.0000   0.00   0.00  3   0.00   0.00    20.00  1.1  0.9;
  36 1.0171  -95.0157   0.0000   0.0000   0.0000   0.0000   0.00   0.00  3   0.00   0.00    20.00  1.1  0.9;
  37 1.0275  -92.4547   0.0000   0.0000   0.0000   0.0000   0.00   0.00  3   0.00   0.00    20.00  1.1  0.9;
  38 1.0347   18.7511   0.0000   0.0000   0.0000   0.0000   0.00   0.00  3   0.00   0.00    20.00  1.1  0.9;
  39 1.0378  -13.2266   0.0000   0.0000   0.0000   0.0000   0.00   0.00  3   0.00   0.00    20.00  1.1  0.9;
  40 0.9654  -55.4078   0.0000   0.0000   0.0000   0.0000   0.00   0.00  3   0.00   0.00   230.00  1.1  0.9;
  41 1.0312   10.2218   0.0000   0.0000   0.0000   0.0000   0.00   0.00  3   0.00   0.00   230.00  1.1  0.9;
  42 1.0146   -1.0369   0.0000   0.0000   0.0000   0.0000   0.00   0.00  3   0.00   0.00   230.00  1.1  0.9];


% line data format
% line: from bus, to bus, resistance(pu), reactance(pu),
%     line charging(pu), tap ratio, tap phase, tapmax, tapmin, tapsize

line = [...
  30  12   0.00007   0.00700   0.00000  0.0  0.0  0.0  0.0  0.0;
  11  29   0.00010   0.00200   0.00000  0.0  0.0  0.0  0.0  0.0;
  27  28   0.00010   0.00075   0.04000  0.0  0.0  0.0  0.0  0.0;
  13  31   0.00060   0.00730   0.00000  0.0  0.0  0.0  0.0  0.0;
  32  14   0.00007   0.00550   0.00000  0.0  0.0  0.0  0.0  0.0;
  22  27   0.00350   0.02460   0.29200  0.0  0.0  0.0  0.0  0.0;
  41  24   0.00850   0.09950   0.50000  0.0  0.0  0.0  0.0  0.0;
  24  40   0.00920   0.10900   0.50000  0.0  0.0  0.0  0.0  0.0;
  33  15   0.00008   0.00530   0.00000  0.0  0.0  0.0  0.0  0.0;
  34  16   0.00007   0.00630   0.00000  0.0  0.0  0.0  0.0  0.0;
  25  26   0.00025   0.00282   0.01800  0.0  0.0  0.0  0.0  0.0;
  25  26   0.00030   0.00279   0.00900  0.0  0.0  0.0  0.0  0.0;
  25  23   0.00190   0.01550   0.09520  0.0  0.0  0.0  0.0  0.0;
  42  24   0.00800   0.08400   0.33400  0.0  0.0  0.0  0.0  0.0;
  35  17   0.00002   0.00050   0.00000  0.0  0.0  0.0  0.0  0.0;
  26  21   0.00190   0.02100   0.22000  0.0  0.0  0.0  0.0  0.0;
  26  22   0.00300   0.01880   0.01800  0.0  0.0  0.0  0.0  0.0;
  36  18   0.00007   0.00600   0.00000  0.0  0.0  0.0  0.0  0.0;
  37  19   0.00007   0.00400   0.00000  0.0  0.0  0.0  0.0  0.0;
  28  20   0.00200   0.01550   0.03200  0.0  0.0  0.0  0.0  0.0;
  38   9   0.00007   0.00600   0.00000  0.0  0.0  0.0  0.0  0.0;
  39  10   0.00007   0.00600   0.00000  0.0  0.0  0.0  0.0  0.0;
  20  40   0.00000   0.02000   0.00000  1.0  0.0  0.0  0.0  0.0;
  20  29   0.00000   0.00200   0.00000  1.0  0.0  0.0  0.0  0.0;
  20  30   0.00000   0.00500   0.00000  1.0  0.0  0.0  0.0  0.0;
  20   3   0.00000   0.00070   0.00000  1.0  0.0  0.0  0.0  0.0;
  21   2   0.00000   0.00100   0.00000  1.0  0.0  0.0  0.0  0.0;
  21  31   0.00000   0.00500   0.00000  1.0  0.0  0.0  0.0  0.0;
  22  32   0.00000   0.01000   0.00000  1.0  0.0  0.0  0.0  0.0;
  22   5   0.00000   0.00080   0.00000  1.0  0.0  0.0  0.0  0.0;
  23  38   0.00000   0.00500   0.00000  1.0  0.0  0.0  0.0  0.0;
  23  41   0.00000   0.02000   0.00000  1.0  0.0  0.0  0.0  0.0;
  24  39   0.00000   0.00400   0.00000  1.0  0.0  0.0  0.0  0.0;
  25  33   0.00000   0.00300   0.00000  1.0  0.0  0.0  0.0  0.0;
  25  34   0.00000   0.00300   0.00000  1.0  0.0  0.0  0.0  0.0;
  25   6   0.00000   0.00090   0.00000  1.0  0.0  0.0  0.0  0.0;
  25  42   0.00000   0.02000   0.00000  1.0  0.0  0.0  0.0  0.0;
  26  35   0.00000   0.00070   0.00000  1.0  0.0  0.0  0.0  0.0;
  26   1   0.00000   0.00200   0.00000  1.0  0.0  0.0  0.0  0.0;
  28  36   0.00000   0.00500   0.00000  1.0  0.0  0.0  0.0  0.0;
  28  37   0.00000   0.00800   0.00000  1.0  0.0  0.0  0.0  0.0;
  28   4   0.00000   0.00080   0.00000  1.0  0.0  0.0  0.0  0.0;
  23   7   0.00000   0.00070   0.00000  1.0  0.0  0.0  0.0  0.0;
  24   8   0.00000   0.00070   0.00000  1.0  0.0  0.0  0.0  0.0];


% Machine data format
% Machine data format
%     1. machine number,
%     2. bus number,
%     3. base mva,
%     4. leakage reactance x_l(pu),
%     5. resistance r_a(pu),
%     6. d-axis sychronous reactance x_d(pu),
%     7. d-axis transient reactance x'_d(pu),
%     8. d-axis subtransient reactance x"_d(pu),
%     9. d-axis open-circuit time constant T'_do(sec),
%    10. d-axis open-circuit subtransient time constant
%          T"_do(sec),
%    11. q-axis sychronous reactance x_q(pu),
%    12. q-axis transient reactance x'_q(pu),
%    13. q-axis subtransient reactance x"_q(pu),
%    14. q-axis open-circuit time constant T'_qo(sec),
%    15. q-axis open circuit subtransient time constant
%          T"_qo(sec),
%    16. inertia constant H(sec),
%    17. damping coefficient d_o(pu),
%    18. dampling coefficient d_1(pu),
%    19. bus number
%
% note: all the following machines use sub-transient model
mac_con = [ ...
   1   1   3300   0.12  0.0   1.70   0.27   0.25   6.0   0.06   1.60   0.50   0.25   0.80   0.04   4.00   0.0   0.0   1;
   2   2  13500   0.12  0.0   0.90   0.27   0.00   9.0   0.05   0.60   0.27   0.00   0.00   0.00   4.17   0.0   0.0   2;
   3   3  13000   0.15  0.0   2.55   0.37   0.30   7.0   0.07   2.40   0.75   0.30   0.90   0.04   4.23   0.0   0.0   3;
   4   4   9000   0.12  0.0   1.70   0.27   0.25   6.0   0.06   1.60   0.50   0.25   0.80   0.04   4.44   0.0   0.0   4;
   5   5   8000   0.12  0.0   0.90   0.27   0.00   9.0   0.05   0.60   0.27   0.00   0.05   0.00   4.38   0.0   0.0   5;
   6   6  14200   0.15  0.0   1.35   0.37   0.00   9.0   0.05   0.90   0.37   0.00   0.05   0.00   4.62   0.0   0.0   6;
   7   7   2300   0.12  0.0   1.70   0.27   0.25   6.0   0.06   1.60   0.50   0.25   0.80   0.04   2.50   0.0   0.0   7;
   8   8   6000   0.12  0.0   1.70   0.27   0.25   6.0   0.06   1.60   0.50   0.25   0.80   0.04   2.67   4.0   0.0   8;
   9   9   2392   0.12  0.0   1.70   0.27   0.25   6.0   0.06   1.60   0.50   0.25   0.80   0.04   4.11   8.0   0.0   9;
  10  10   1982   0.12  0.0   0.90   0.27   0.00   9.0   0.05   0.60   0.27   0.00   0.05   0.00   2.88   5.0   0.0   10;
  11  11   3118   0.12  0.0   1.70   0.27   0.25   6.0   0.06   1.60   0.50   0.25   0.80   0.04   3.83   7.0   0.0   11;
  12  12    640   0.12  0.0   1.70   0.27   0.25   6.0   0.06   1.60   0.50   0.25   0.80   0.04   3.23   6.0   0.0   12;
  13  13    880   0.12  0.0   0.90   0.27   0.00   9.0   0.05   0.60   0.27   0.00   0.00   0.00   5.10   10.0   0.0   13;
  14  14   1500   0.12  0.0   1.70   0.27   0.25   6.0   0.06   1.60   0.50   0.25   0.80   0.04   5.36   11.0   0.0   14;
  15  15    778   0.12  0.0   0.90   0.27   0.00   9.0   0.05   0.60   0.27   0.00   0.05   0.00   5.40   11.0   0.0   15;
  16  16   1137   0.12  0.0   0.90   0.27   0.00   9.0   0.05   0.60   0.27   0.00   0.05   0.00   2.98   6.0   0.0   16;
  17  17   2208   0.12  0.0   1.70   0.27   0.25   6.0   0.06   1.60   0.50   0.25   0.80   0.04   4.82   9.0   0.0   17;
  18  18   3004   0.12  0.0   1.70   0.27   0.25   6.0   0.06   1.60   0.50   0.25   0.80   0.04   5.16   10.0   0.0   18;
  19  19   2500   0.12  0.0   1.70   0.27   0.25   6.0   0.06   1.60   0.50   0.25   0.80   0.04   3.80   8.0   0.0   19];

exc_con = [...
    0    1    0.0   25.0   0.02    0.0    0.0    6.0   -3.0;
    0    2    0.0   25.0   0.02    0.0    0.0    6.0   -3.0;
    0    3    0.0   25.0   0.02    0.0    0.0    7.0   -2.0;
    0    4    0.0   25.0   0.02    0.0    0.0    7.0   -3.0;
    0    5    0.0   25.0   0.02    0.0    0.0    7.0   -2.0;
    0    6    0.0   25.0   0.02    0.0    0.0    6.0   -3.0;
    0    7    0.0   25.0   0.02    0.0    0.0    6.0   -3.0;
    0    8    0.0   25.0   0.02    0.0    0.0    6.0   -3.0;
    0    9    0.0  100.0   0.02    0.0    0.0    6.0   -3.0;
    0   10    0.0  100.0   0.02    0.0    0.0    6.0   -3.0;
    0   11    0.0  150.0   0.02    0.0    0.0    7.0   -2.0;
    0   12    0.0  150.0   0.02    0.0    0.0    7.0   -2.0;
    0   13    0.0  100.0   0.02    0.0    0.0    6.0   -3.0;
    0   14    0.0  150.0   0.02    0.0    0.0    7.0   -2.0;
    0   15    0.0  100.0   0.02    0.0    0.0    6.0   -3.0;
    0   16    0.0  100.0   0.02    0.0    0.0    6.0   -3.0;
    0   17    0.0  100.0   0.02    0.0    0.0    6.0   -3.0;
    0   18    0.0  200.0   0.02    0.0    0.0    7.0   -3.0;
    0   19    0.0  200.0   0.02    0.0    0.0    7.0   -3.0];

% governor model

tg_con = [...
   1   9   1  20.0    1.00   0.04   0.2    0.0    1.5    5.0;
   1  10   1  20.0    1.00   0.04   6.0    0.0   -2.4    1.2;
   1  11   1  20.0    1.00   0.04   0.2    0.0    1.5    5.0;
   1  12   1  20.0    1.00   0.04   0.2    0.0    1.5    5.0;
   1  13   1  20.0    1.00   0.04   6.0    0.0   -2.4    1.2;
   1  14   1  20.0    1.00   0.04   0.2    0.0    1.5    5.0;
   1  15   1  20.0    1.00   0.04   6.0    0.0   -2.4    1.2;
   1  16   1  20.0    1.00   0.04   6.0    0.0   -2.4    1.2;
   1  17   1  20.0    1.00   0.04   6.0    0.0   -2.4    1.2;
   1  18   1  20.0    1.00   0.04   0.2    0.0    1.5    5.0;
   1  19   1  20.0    1.00   0.04   0.2    0.0    1.5    5.0];

% non-conforming load
% col 1       bus number
% col 2       fraction const active power load
% col 3       fraction const reactive power load
% col 4       fraction const active current load
% col 5       fraction const reactive current load

load_con = [...
20  0.0  0.0  0.50  0.0;
21  0.0  0.0  0.50  0.0;
22  0.0  0.0  0.50  0.0;
23  0.0  0.0  0.50  0.0;
24  0.0  0.0  0.50  0.0;
25  0.0  0.0  0.50  0.0;
26  0.0  0.0  0.50  0.0];
%27  0.0  0.0  0.00  0.0];
%28  0.0  0.0  0.50  0.0];

% load modulation control at buses 25 and 26
lmod_con = [ 1 26 1000 .1  -.1  10  0.02;
             2 25 1000 .1  -.1  10  0.05];

%Switching file defines the simulation control
% row 1 col1  simulation start time (s) (cols 2 to 6 zeros)
%     col7  initial time step (s)
% row 2 col1  fault application time (s)
%     col2  bus number at which fault is applied
%     col3  bus number defining far end of faulted line
%     col4  zero sequence impedance in pu on system base
%     col5  negative sequence impedance in pu on system base
%     col6  type of fault  - 0 three phase
%                  - 1 line to ground
%                  - 2 line-to-line to ground
%                  - 3 line-to-line
%                  - 4 loss of line with no fault
%                  - 5 loss of load at bus
%                  - 6 for no change
%     col7  time step for fault period (s)
% row 3 col1  near end fault clearing time (s) (cols 2 to 6 zeros)
%     col7  time step for second part of fault (s)
% row 4 col1  far end fault clearing time (s) (cols 2 to 6 zeros)
%     col7  time step for fault cleared simulation (s)
% row 5 col1  time to change step length (s)
%     col7  time step (s)
%
%
%
% row n col1 finishing time (s)  (n indicates that intermediate rows may be inserted)

sw_con = [...
0    0    0    0    0    0    0.02;%sets intitial time step
0.1  26   25    0.0  0.0  0    0.02;%three phase fault at bus 27
0.2  0   0    0    0    0    0.02;% clear near end
0.3   0   0    0    0    0    0.02;% clear remote end 
5.0  0   0    0    0    0    0  ];% end simulation

lmon_con = [1:length(line(:,1))]';
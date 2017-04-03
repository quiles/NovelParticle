clear all;

a = [11	28.000000	0.133504	0.555556	0.013533	0.133504	0.133504	42.544858
39	37.000000	0.153584	0.433934	0.017607	0.153584	0.153584	131.740041
48	38.000000	0.152090	0.436700	0.020007	0.152090	0.152090	360.166033
54	39.000000	0.158202	0.441296	0.018471	0.158202	0.158202	122.668348
69	12.000000	0.052080	0.590909	0.006870	0.052080	0.052080	11.736710
60	53.000000	0.190691	0.325109	0.025962	0.190691	0.190691	514.405028
33	16.000000	0.071957	0.508333	0.008669	0.071957	0.071957	24.597684
84	17.000000	0.081739	0.625000	0.008808	0.081739	0.081739	13.623010
57	29.000000	0.129911	0.539409	0.013938	0.129911	0.129911	47.407890
24	33.000000	0.138598	0.460227	0.015878	0.138598	0.138598	96.498672
26	24.000000	0.104435	0.532609	0.012013	0.104435	0.104435	40.102542
65	17.000000	0.070939	0.566176	0.008890	0.070939	0.070939	14.096719
22	19.000000	0.086284	0.549708	0.009870	0.086284	0.086284	25.986816
4	46.000000	0.185806	0.407729	0.021610	0.185806	0.185806	223.987403
38	14.000000	0.055543	0.494505	0.007846	0.055543	0.055543	21.035697
25	56.000000	0.215800	0.363636	0.026252	0.215800	0.215800	386.994423
71	25.000000	0.115248	0.606667	0.012248	0.115248	0.115248	31.961688
30	45.000000	0.177261	0.393939	0.021788	0.177261	0.177261	293.946451
83	31.000000	0.135232	0.477419	0.014967	0.135233	0.135233	86.049789
56	18.000000	0.085572	0.594771	0.009347	0.085572	0.085572	18.561708
16	24.000000	0.106955	0.594203	0.012041	0.106955	0.106955	37.877328
29	18.000000	0.077998	0.607843	0.009466	0.077998	0.077998	16.925720
66	25.000000	0.105739	0.546667	0.012514	0.105739	0.105739	42.698627
21	23.000000	0.104929	0.584980	0.011533	0.104930	0.104930	29.985457
64	26.000000	0.111869	0.489231	0.012992	0.111869	0.111869	57.081616
14	25.000000	0.105574	0.533333	0.012585	0.105574	0.105574	50.359271
77	22.000000	0.088151	0.493506	0.011537	0.088151	0.088151	58.075507
61	34.000000	0.134794	0.374332	0.016500	0.134794	0.134794	118.633857
70	13.000000	0.053116	0.487179	0.007419	0.053116	0.053116	16.893955
80	37.000000	0.153478	0.436937	0.019028	0.153478	0.153478	296.628453
85	42.000000	0.162635	0.379791	0.021640	0.162635	0.162635	411.239592
0	14.000000	0.050013	0.516484	0.008342	0.050013	0.050013	34.369647
89	10.000000	0.043038	0.800000	0.006063	0.043038	0.043038	2.294793
67	24.000000	0.103338	0.510870	0.012085	0.103338	0.103338	53.092710
35	51.000000	0.198717	0.376471	0.024091	0.198717	0.198717	338.565036
23	41.000000	0.160174	0.408537	0.019877	0.160175	0.160175	211.074711
7	21.000000	0.102664	0.661905	0.010526	0.102664	0.102664	16.022098
59	11.000000	0.053230	0.709091	0.006331	0.053230	0.053230	4.024800
62	22.000000	0.100961	0.619048	0.011064	0.100961	0.100961	22.361432
5	11.000000	0.053809	0.763636	0.006302	0.053810	0.053810	2.301075
44	34.000000	0.141600	0.452763	0.016661	0.141600	0.141600	120.857152
81	11.000000	0.043562	0.600000	0.006401	0.043562	0.043562	5.292582
88	38.000000	0.153125	0.416785	0.019438	0.153125	0.153125	287.593567
49	59.000000	0.227210	0.359439	0.027415	0.227210	0.227210	412.633107
87	40.000000	0.162955	0.425641	0.019070	0.162955	0.162955	165.741172
58	38.000000	0.158371	0.445235	0.018255	0.158371	0.158371	137.055668
15	19.000000	0.076963	0.543860	0.010266	0.076963	0.076963	37.441389
1	31.000000	0.133358	0.509677	0.015083	0.133358	0.133358	72.725009
19	16.000000	0.066852	0.583333	0.008704	0.066852	0.066852	16.608793
3	7.000000	0.027048	0.380952	0.004731	0.027048	0.027048	2.803674
53	13.000000	0.064316	0.730769	0.007143	0.064316	0.064316	3.985136
82	11.000000	0.040845	0.345455	0.006505	0.040845	0.040845	11.875615
51	11.000000	0.040109	0.436364	0.006520	0.040109	0.040109	11.095443
55	21.000000	0.089253	0.461905	0.010664	0.089253	0.089253	31.101467
43	35.000000	0.140336	0.420168	0.017467	0.140336	0.140336	190.767429
50	18.000000	0.077352	0.679739	0.009285	0.077352	0.077352	9.630143
12	16.000000	0.090317	0.766667	0.008281	0.090317	0.090317	3.896552
9	18.000000	0.085744	0.588235	0.009259	0.085744	0.085744	14.404762
46	24.000000	0.107779	0.543478	0.011907	0.107779	0.107779	38.931739
78	16.000000	0.083127	0.700000	0.008358	0.083127	0.083127	7.531771
31	8.000000	0.034825	0.464286	0.005295	0.034825	0.034825	6.447076
52	31.000000	0.131038	0.505376	0.014920	0.131038	0.131038	75.731416
32	16.000000	0.070369	0.566667	0.008570	0.070369	0.070369	19.160263
86	10.000000	0.049044	0.688889	0.005910	0.049044	0.049044	3.208040
2	18.000000	0.077408	0.549020	0.009423	0.077408	0.077408	24.104579
8	20.000000	0.092003	0.615789	0.010125	0.092003	0.092003	15.699623
17	27.000000	0.115372	0.504274	0.013199	0.115372	0.115372	53.825253
63	30.000000	0.127085	0.494253	0.014431	0.127085	0.127085	59.244024
18	22.000000	0.095002	0.528139	0.011047	0.095002	0.095002	29.706082
68	6.000000	0.031226	0.600000	0.004295	0.031226	0.031226	1.450314
37	23.000000	0.095257	0.505929	0.011622	0.095257	0.095257	39.249766
13	19.000000	0.085205	0.549708	0.010017	0.085205	0.085205	33.741609
27	10.000000	0.043093	0.711111	0.005949	0.043093	0.043093	3.700267
72	3.000000	0.014315	0.333333	0.003020	0.014315	0.014315	0.922761
40	2.000000	0.010015	1.000000	0.002557	0.010015	0.010015	0.000000
79	1.000000	0.005244	0.000000	0.002120	0.005244	0.005244	0.000000
41	1.000000	0.005256	0.000000	0.002123	0.005256	0.005256	0.000000
76	20.000000	0.085656	0.589474	0.010241	0.085656	0.085656	25.522293
36	5.000000	0.024947	0.900000	0.003815	0.024947	0.024947	0.153846
73	36.000000	0.145396	0.401587	0.017497	0.145396	0.145396	149.552020
6	1.000000	0.005570	0.000000	0.002123	0.005570	0.005570	0.000000
75	4.000000	0.018823	0.833333	0.003400	0.018823	0.018823	0.333333
47	17.000000	0.076488	0.573529	0.009194	0.076488	0.076488	22.916049
74	12.000000	0.050274	0.575758	0.006973	0.050274	0.050274	9.269491
42	9.000000	0.039164	0.583333	0.005584	0.039164	0.039164	5.493074
45	10.000000	0.048057	0.644444	0.005929	0.048057	0.048057	4.524042
10	12.000000	0.058858	0.681818	0.006774	0.058858	0.058858	5.035942
28	5.000000	0.026230	0.700000	0.003788	0.026230	0.026230	0.419472
34	1.000000	0.005209	0.000000	0.002133	0.005209	0.005209	0.000000];

[lixo indices] = sort(a(:,1),1);

x = a(indices,:);

for i=3:8
    x(:,i) = x(:,i) / max(x(:,i));
end;

plot(x(:,1),x(:,3:8));
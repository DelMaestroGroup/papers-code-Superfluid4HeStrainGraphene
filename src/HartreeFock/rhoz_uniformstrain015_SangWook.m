
rhoz_vec_SW = [...
0.000000000000000000e+00
0.000000000000000000e+00
0.000000000000000000e+00
0.000000000000000000e+00
0.000000000000000000e+00
0.000000000000000000e+00
0.000000000000000000e+00
0.000000000000000000e+00
0.000000000000000000e+00
0.000000000000000000e+00
0.000000000000000000e+00
0.000000000000000000e+00
0.000000000000000000e+00
0.000000000000000000e+00
0.000000000000000000e+00
0.000000000000000000e+00
4.244692585185185583e-12
1.501443571525926048e-10
8.757993925519998694e-09
1.627958496799404838e-07
1.846153690024632078e-06
1.317192535897784694e-05
6.438184375659552247e-05
2.274737995624642880e-04
6.077299935503631791e-04
1.271225053907265763e-03
2.154897225804811283e-03
3.049098312377923527e-03
3.699414584411924497e-03
3.938291179588886311e-03
3.752955157500743877e-03
3.255000696995933272e-03
2.606933256611258798e-03
1.950102648886002654e-03
1.376655162212069797e-03
9.246328565041717247e-04
5.948416228158295415e-04
3.687596985715033485e-04
2.213003477703510491e-04
1.292104695944461105e-04
7.358992276045553296e-05
4.105550712029582179e-05
2.248078679116319447e-05
1.211224996568479536e-05
6.419308182566637860e-06
3.369523781001308100e-06
1.755386966860803694e-06
9.040698955670894164e-07
4.646032291643021355e-07
2.428446586280858147e-07
1.226883792988240039e-07
6.576186010505476862e-08
3.382692654764891910e-08
1.979491761106371583e-08
1.208567434330222536e-08
7.311268964922960910e-09
4.772405006063704269e-09
3.839814554814815275e-09
2.808204561862222539e-09
2.594216834788148112e-09
1.698756354285925615e-09
1.508152137537777822e-09
1.498530612442962777e-09
1.016769092802963001e-09
8.856925865540740606e-10
8.337117354992592139e-10
6.117923763911110427e-10
6.619059181970371887e-10
5.747134121866665954e-10
5.160439545377778304e-10
3.310605178903703351e-10
3.101752554296296126e-10
3.112098159955555595e-10
2.351524208711111262e-10
2.458949969125926253e-10
1.731253533600000146e-10
1.221491236859259300e-10
2.174839432888889156e-10
1.048698540044444432e-10
7.564601041925926168e-11
7.700546520740740915e-11
3.570398550370370790e-11
0.000000000000000000e+00
0.000000000000000000e+00
0.000000000000000000e+00
0.000000000000000000e+00
0.000000000000000000e+00
0.000000000000000000e+00
0.000000000000000000e+00
0.000000000000000000e+00
0.000000000000000000e+00
0.000000000000000000e+00
0.000000000000000000e+00
0.000000000000000000e+00
0.000000000000000000e+00
0.000000000000000000e+00
0.000000000000000000e+00
0.000000000000000000e+00
0.000000000000000000e+00
0.000000000000000000e+00
0.000000000000000000e+00];

z_vec_SW = [0 : length(rhoz_vec_SW)-1]*9.900990e-02;

figure(1); plot(z_vec_SW,rhoz_vec_SW*0.0236/max(rhoz_vec_SW),'b')

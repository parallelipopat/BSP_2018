energy = [-5.238638341	1.68584198	-3.552796361
-5.208970593	1.639343181	-3.569627412
-5.162222297	1.600640261	-3.561582036
-5.181944444	1.625267792	-3.556676652
-5.228192149	1.687970926	-3.540221223
-5.18644769	1.638969736	-3.547477954
-5.159506128	1.611067259	-3.548438869
-5.16535578	1.619590087	-3.545765693
-5.140524905	1.59694079	-3.543584116
-5.133143296	1.580572853	-3.552570443
-5.210543722	1.688702757	-3.521840965
-5.173181327	1.633796076	-3.539385252
-5.182055203	1.65175217	-3.530303032
-5.172668669	1.624724996	-3.547943674
-5.175983304	1.645802122	-3.530181181
-5.146396447	1.597830099	-3.548566347
-5.149546018	1.611500905	-3.538045113
-5.148408409	1.616613334	-3.531795075
-5.163699012	1.641684726	-3.522014286
-5.166241561	1.645929238	-3.520312324
-5.151722222	1.631393345	-3.520328877
-5.162989932	1.645801553	-3.517188379
-5.197233035	1.690235354	-3.50699768
-5.19105368	1.689713514	-3.501340166
-5.161310613	1.649056902	-3.512253711
-5.167164469	1.656263385	-3.510901084
-5.158930125	1.654733064	-3.504197061
-5.112301587	1.602218551	-3.510083036
-5.099264411	1.57949138	-3.519773031
-5.148609012	1.645265129	-3.503343884
-5.126203026	1.613006577	-3.513196449
-5.126269857	1.61705943	-3.509210427
-5.18068473	1.695613509	-3.485071221
-5.161854069	1.665340196	-3.496513873
-5.161505924	1.67384369	-3.487662234
-5.119669094	1.625189585	-3.494479509
-5.126942275	1.64145553	-3.485486745
-5.128545028	1.640395683	-3.488149345
-5.113290673	1.629919079	-3.483371594
-5.139571261	1.653302662	-3.4862686
-5.162744466	1.691180625	-3.471563841
-5.114944219	1.635092649	-3.47985157
-5.131009246	1.666707745	-3.464301501
-5.113237335	1.647234255	-3.46600308
-5.095517258	1.616508482	-3.479008776
-5.125417482	1.659273387	-3.466144095
-5.087672229	1.605641566	-3.482030663
-5.118499942	1.650045123	-3.468454819
-5.138045773	1.675908507	-3.462137267
-5.133148149	1.670624586	-3.462523563
];

steps = linspace(1,50,50);
figure;
plot(steps, energy(:,1), '-^', steps, energy(:,2), '-v', steps, energy(:,3), '-o');
xlabel('Iterations'); ylabel('Energy (\epsilon)');
legend('Potential Energy', 'Kinetic Energy', 'Total Energy','Location','east');
title('\beta = 0.40, \gamma = 0.50');

figure;
plot(steps,energy(:,3)/mean(energy(:,3))-1,'o-');
xlabel('Iterations'); ylabel('$\frac{E-<E>}{<E>}$','Interpreter','latex','FontSize',14);
title('\beta = 0.40, \gamma = 0.50');
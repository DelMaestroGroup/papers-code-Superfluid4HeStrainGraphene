% This is functionally the same code as  rhorel_of_z_and_softening in
% folder X:\all_latexts\proposals\NASA_2018Fall_Epscor_PhysDept\paper_1\particle_density
% except that here, names are changed (AUX_ added to each name) to avoid confusion withe main code
% and plotting is suppressed. 


% LINES TO COMMENT WHEN RUNNING INSIDE   psis_... :
% 13 <- clear all
% 40 <- uniform_strain = ... (defined in psis_...) 
% 267-279, 322-326  <- plots of rho and rho_rel 
% 287 <- He_Aziz95_reduced; (loaded in psis_...) 

% clear all


% Part 1:  Compute the distribution of the relative z-positions of two neighboring He atoms
%          from single-particle distribution function rho(z).
%          Namely, we will compute: rho_rel(z) = \int rho(z') rho(z + z') dz'

% Load the single-particle density and z:

% rho(z)'s from different sources:
which_rhoz_to_load = 5;              % 1 => from our V0(z) potential as of December 2020
                                     % 2 => from Adrian's QMC, December 2020 
                                     % 3 => from V0(z) for different values of uniform strain, Feb 2021
                                     %      (costfunction sum (del_U/U)^2, "mine", eps=15.... for no strain)
                                     % 4 => from V0(z) for different values of uniform strain, Feb 2021
                                     %      (costfunction sum (del_U)^2, "Adrian's", eps=16.96 for no strain)
                                     % 5 => from Sang Wook's QMC, May 2022
                                     %      (must be the same as Andrian's 2020 data, but for several
                                     %       uniform strains and using different format of save data)
                                     
if  which_rhoz_to_load == 1
    rhoz_from1DV0z;  
    column_rho = 2;                  % the column in the file where rho is stored
                                     %  (it =2 for which_rhoz_to_load == 1,2, 
                                     %   but 2...8 for which_rhoz_to_load == 3)
elseif  which_rhoz_to_load == 2
    rhoz_fromQMC; 
    column_rho = 2;
elseif  which_rhoz_to_load == 3
    rhoz_uniformstrain_newcost;
    % uniform_strain = 0.3; % COMMENT THIS OUT WHEN RUNNING  psis_withHeHe_by_Hartree_FOCK... !!!!
                        % This value is set there.
    column_rho = round(uniform_strain/0.05) + 2;
elseif  which_rhoz_to_load == 4
    rhoz_uniformstrain_oldcost;
    % uniform_strain = 0.3; % COMMENT THIS OUT WHEN RUNNING  psis_withHeHe_by_Hartree_FOCK... !!!!
                        % This value is set there.
    column_rho = round(uniform_strain/0.05) + 2;
elseif  which_rhoz_to_load == 5
    if uniform_strain == 0 || uniform_strain == 0.05 
        runfilename = ['rhoz_uniformstrain00' int2str(uniform_strain*100) '_SangWook'];
    elseif  uniform_strain >= 0.099 && uniform_strain <= 0.3099  % i.e. for uniform_strain >= 0.10
        runfilename = ['rhoz_uniformstrain0' int2str(uniform_strain*100) '_SangWook'];
    else  % % i.e. for uniform_strain > 0.30
        runfilename = 'rhoz_uniformstrain030_SangWook';
    end
    run(runfilename);
end

% % z [A]     rho
% AUX_z_and_rho = [ ... 
% 1.95	1.6842363774762446e-11
% 1.96	3.373255519045039e-11
% 1.97	6.756302418125185e-11
% 1.98	1.353062934504598e-10
% 1.99	2.7109215574198853e-10
% 2.0	5.418382818283108e-10
% 2.01	1.0717912541863785e-9
% 2.02	2.0861475956681734e-9
% 2.03	3.98743354197179e-9
% 2.04	7.478303723050501e-9
% 2.05	1.3760510159950649e-8
% 2.06	2.484731563670723e-8
% 2.07	4.404487644877966e-8
% 2.08	7.667750819052304e-8
% 2.09	1.3115786555086714e-7
% 2.1	2.205332700076002e-7
% 2.11	3.6467444115199293e-7
% 2.12	5.933106787013195e-7
% 2.13	9.501540144304587e-7
% 2.14	1.4983903056598198e-6
% 2.15	2.327842992830389e-6
% 2.16	3.564129395631237e-6
% 2.17	5.380123266416692e-6
% 2.18	8.01000372082165e-6
% 2.19	1.1766102738266227e-5
% 2.2	1.7058655338230268e-5
% 2.21	2.44184044514654e-5
% 2.22	3.4521815526954775e-5
% 2.23	4.8218416983304634e-5
% 2.24	6.655950917179534e-5
% 2.25	9.082718887111325e-5
% 2.26	0.00012256233556535642
% 2.27	0.00016358992115176972
% 2.28	0.00021603976065873694
% 2.29	0.00028236064389176266
% 2.3	0.0003653257020128318
% 2.31	0.0004680268914719485
% 2.32	0.0005938566379414949
% 2.33	0.0007464749852207419
% 2.34	0.0009297610397691447
% 2.35	0.0011477480817536972
% 2.36	0.0014045424088948034
% 2.37	0.0017042267605388995
% 2.38	0.0020507499981093757
% 2.39	0.0024478055496197095
% 2.4	0.0028987019116886057
% 2.41	0.0034062291933137823
% 2.42	0.0039725262352276166
% 2.43	0.0045989532068537505
% 2.44	0.005285974738755662
% 2.45	0.006033058572599739
% 2.46	0.006838594396738814
% 2.47	0.007699836991011062
% 2.48	0.008612877050147712
% 2.49	0.00957264212440801
% 2.5	0.010572929051985461
% 2.51	0.011606468111083946
% 2.52	0.012665017945392829
% 2.53	0.013739489171030813
% 2.54	0.014820093509613273
% 2.55	0.015896514359344843
% 2.56	0.01695809395448953
% 2.57	0.01799403170393617
% 2.58	0.018993587961479597
% 2.59	0.019946287371883538
% 2.6	0.020842116054355007
% 2.61	0.021671707214722575
% 2.62	0.022426510296138773
% 2.63	0.023098939454869617
% 2.64	0.023682497946691573
% 2.65	0.024171875891430404
% 2.66	0.024563019808103416
% 2.67	0.024853173241868795
% 2.68	0.025040888700277042
% 2.69	0.02512601194819236
% 2.7	0.02510964045163617
% 2.71	0.024994058390240402
% 2.72	0.024782651161964576
% 2.73	0.024479802674623886
% 2.74	0.024090778955074103
% 2.75	0.023621601712574727
% 2.76	0.023078915476589513
% 2.77	0.022469851803574446
% 2.78	0.021801893827472277
% 2.79	0.021082744131860842
% 2.8	0.02032019856608982
% 2.81	0.019522028231452236
% 2.82	0.01869587144396143
% 2.83	0.017849137053840967
% 2.84	0.016988920082825248
% 2.85	0.016121930241165938
% 2.86	0.015254433516917137
% 2.87	0.014392206698353768
% 2.88	0.013540504401674706
% 2.89	0.012704037933682046
% 2.9	0.011886965124149957
% 2.91	0.011092890114597858
% 2.92	0.010324871987195608
% 2.93	0.009585441056410106
% 2.94	0.008876621622726509
% 2.95	0.008199959997708796
% 2.96	0.007556556647817715
% 2.97	0.006947101365667695
% 2.98	0.006371910456734821
% 2.99	0.0058309650221092316
% 3.0	0.0053239495192509015
% 3.01	0.0048502898888320165
% 3.02	0.004409190643104419
% 3.03	0.003999670416818792
% 3.04	0.0036205955830814814
% 3.05	0.0032707116317259102
% 3.06	0.0029486720953538397
% 3.07	0.0026530648871764184
% 3.08	0.0023824359845704864
% 3.09	0.0021353104526289704
% 3.1	0.001910210852989139
% 3.11	0.0017056731251754503
% 3.12	0.0015202600610908213
% 3.13	0.0013525725187688402
% 3.14	0.0012012585397940595
% 3.15	0.0010650205467005412
% 3.16	0.0009426208029876125
% 3.17	0.0008328853199606414
% 3.18	0.000734706392202751
% 3.19	0.0006470439378573179
% 3.2	0.0005689258117424134
% 3.21	0.0004994472492549126
% 3.22	0.000437769587613149
% 3.23	0.00038311839872152205
% 3.24	0.0003347811552371164
% 3.25	0.0002921045386282866
% 3.26	0.00025449148542567456
% 3.27	0.00022139805570579415
% 3.28	0.0001923301962914527
% 3.29	0.0001668404603298731
% 3.3	0.00014452473490523434
% 3.31	0.00012501901920895813
% 3.32	0.00010799628755039945
% 3.33	9.316346414038361e-5
% 3.34	8.02585300986841e-5
% 3.35	6.904777748754753e-5
% 3.36	5.9323220309193016e-5
% 3.37	5.0900168270624625e-5
% 3.38	4.3614965654070344e-5
% 3.39	3.7322894773460646e-5
% 3.4	3.189624118368184e-5
% 3.41	2.7222515978076652e-5
% 3.42	2.320282910129887e-5
% 3.43	1.97504065628233e-5
% 3.44	1.6789243708630335e-5
% 3.45	1.425288624648513e-5
% 3.46	1.2083330479879803e-5
% 3.47	1.0230034147628571e-5
% 3.48	8.64902935522597e-6
% 3.49	7.302129289566396e-6
% 3.5	6.156220703671295e-6
% 3.51	5.1826345196233516e-6
% 3.52	4.356587306364847e-6
% 3.53	3.6566868279195604e-6
% 3.54	3.064495313305645e-6];

% Extract  z  and  rho  from the loaded data:
if  which_rhoz_to_load < 5
    AUX_z = AUX_z_and_rho(:,1)';     % z
    AUX_rho = AUX_z_and_rho(:, column_rho)';   % rho
else  % i.e. for which_rhoz_to_load = 5
    AUX_z = z_vec_SW;     % z
    AUX_rho = rhoz_vec_SW';   % rho
end
AUX_Nz = length(AUX_z);
AUX_dz = AUX_z(2) - AUX_z(1);
AUX_rho = AUX_rho / (trapz(AUX_rho)*AUX_dz);  % makes sure that \int rho dz = 1
          % NOTE ABOUT NORMALIZATION:
          % We do *not* normalize  z  for two reasons:
          %    1) Later I compute  r^2 + z^2, where z is in also in [A]; so it is easier
          %       not to renormalize  z  for that step;
          %    2) It is not needed, since we compute the main outcome, 
          %       \int  VHeHe( sqrt(r^2+z^2) ) rho_rel(z) dz,
          %       and as long as we use the same units where  \int rho_rel(z) dz = 1

% To compute rho_rel(z) = \int rho(z') rho(z + z') dz', 
%  we need to first put  rho  on a larger  z-interval (length = 3Nz-2; see immediately below).  
AUX_rho_augm = [zeros(1,AUX_Nz-1)  AUX_rho  zeros(1,AUX_Nz-1)];
% Explanation for this step:
% - Suppose rho(kz) is originally defined (as nonzero) for kz = 1 : Nz.
% - Let rho_rel(kz_rel) be defined for kz_rel = Nmin_rel : Nmax_rel, where N{min,max}_rel are TBD. 
%   First, Nmin_rel = -(Nz-1),  since then both z' and z+z' are still defined for z' = Nz;
%   in other words,  rho_rel  will get a contribution from rho(Nz)*rho(-(Nz-1)+Nz).
%   Second, and similarly,  Nmax_rel = Nz-1, since then rho_rel gets a contribution from 
%   rho(1)*rho((Nz-1)+1).
%   *********  Thus, the size of  rho_rel = length( -(Nz-1) : Nz-1 ) = 2Nz-1.  ***************
% - Our goal is to redefine  rho  on a larger interval so that all indices within this interval
%   are defined when we compute  \int rho(z') rho(z + z') dz'.
% - For z = -(Nz-1),  min(z + z'_redefined) = 1,  =>  min(z'_redefined) = Nz, 
%   and this must correspond to z'_original = 1. Thus, I need to "pad" the original rho(z) 
%   with (Nz-1) zeros at the front.
% - For z = Nz-1,  max(z + z'_redefined) = Nz-1 + (Nz-1) + ((Nz)),
%   where:
%   Nz-1 = max(z),
%   (Nz-1) = number of the padded 0's (see the previous step), and
%   ((Nz)) = max(z') in the original rho.
%   Thus, the maximum index of the redefined  rho  must run up to 3Nz-2.

% Auxiliary vector of length (2Nz-1);  it contains values of rel distance between He atoms.
% rho_rel  will be defined over this vector.
AUX_deltavec = -(AUX_Nz-1)*AUX_dz : AUX_dz : (AUX_Nz-1)*AUX_dz;   

% Computing rho_rel(z);  the loop below goes over each value of z:
for AUX_k = 1 : 2*AUX_Nz-1
    AUX_rho_rel(AUX_k) = trapz( AUX_rho_augm(AUX_Nz : 2*AUX_Nz-1) .* ...
                            AUX_rho_augm(AUX_Nz + (AUX_k-AUX_Nz) : 2*AUX_Nz-1 + (AUX_k-AUX_Nz)) )...
                            *AUX_dz;
end

AUX_rho_rel_cut = AUX_rho_rel( AUX_rho_rel >= 1e-5 );
                  % discard values of rho_rel below 1e-5
AUX_deltavec_cut = AUX_deltavec( AUX_rho_rel >= 1e-5 );
AUX_Nz_cut = ( length(AUX_rho_rel_cut) + 1)/2;
             % rho_rel_cut(z)  in defined on z\in [ -AUX_Nz_cut, AUX_Nz_cut ];
             %  This is convenient in computing 
             %  VHeHe_softened = \int VHeHe( sqrt(r^2 + z^2) ) * rho_rel_cut(z) dz
             %   in that we can set integration limits from 0 to AUX_Nz_cut
             %   to reduce the time required for this summation by the factor of 2. 
             
rhonorm_check = trapz( AUX_rho_rel_cut )*AUX_dz;
disp([' int rhorel (z) dz = '  num2str(rhonorm_check) ])
             % Note that since AUX_dz has dimension [A], then AUX_rho_rel_cut has dimension [1/A].
             % However, we don't need to nondimensionalize either z_rel or rho_rel since below we compute
             % int Vhehe( sqrt(r^2+z_rel^2) )*rho_rel(z) dz_rel, 
             % and so such a scaling will not matter. 

% 
% figure(71);
% plot(AUX_z, AUX_rho, 'y', 'linewidth', 2)
% set(gca, 'fontsize', 20)
% xlabel('z,  [A]', 'fontsize', 20)
% ylabel('\rho (z)', 'fontsize', 20)
% xlim([min(AUX_z), max(AUX_z)])
% 
% figure(72);
% plot(AUX_deltavec_cut, AUX_rho_rel_cut, 'y', 'linewidth', 2)
% set(gca, 'fontsize', 20)
% xlabel('z,  [A]', 'fontsize', 18)
% ylabel('\rho_{relative} (z)', 'fontsize', 18)
% xlim([min(AUX_deltavec_cut), max(AUX_deltavec_cut)])
% 
%  sffa
 
 
% ------------------------------------------------------------------

% Part 2: Use rho_rel to compute the softened potential VHeHe from:
%         VHeHe_softer(r) = int VHeHe( sqrt(r^2 + del^2) ) rho_rel( del) d del

% He_Aziz95_reduced;   % <-- This is loaded by the main code. 
AUX_numr = size(r_V_95reduced, 1);      % length of VHeHe data 
AUX_dr = r_V_95reduced(2,1) - r_V_95reduced(1,1);

for  AUX_jr = 1 : AUX_numr
    
    AUX_rj_sq = r_V_95reduced(AUX_jr, 1)^2;   % r^2, same for all r2+z_rel^2 in the loop below

    % We compute 
    %   Vhehe_softened(AUX_jr) = 
    %   sum_{delz = - dz*(AUX_Nz_cut-1) : dz*(AUX_Nz_cut-1)} Vhehe_orig(r^2+dz^2)*rho_rel(delz)*d(delz) 
    %   = ( Vhehe_orig(r^2+0)*rho_rel(0) + 
    %       2*sum{delz = dz : dz*(AUX_Nz_cut-1)} Vhehe_orig(r^2+delz^2)*rho_rel(delz) ) * d(delz),
    %   where I used the fact that rho(delz) is symmetric about delz=0 and since (-delz)^2 = delz^2. 
    Vhehe_softened(AUX_jr) = r_V_95reduced(AUX_jr,2) * AUX_rho_rel_cut(AUX_Nz_cut) *AUX_dz; 
                             %  AUX_rho_rel_cut(AUX_Nz_cut) = rho_rel at z_rel = 0  - 
                             %  this is the first term in the sum above
    % Now compute the remaining (AUX_Nz_cut-1) terms in the above sum (and then multiply it by 2).
    for  AUX_k_cut = AUX_Nz_cut+1 : 2*AUX_Nz_cut-1    % since AUX_deltavec_cut(AUX_Nz_cut+1)=0+1*dz, etc
        AUX_rj_soft = sqrt( AUX_rj_sq + AUX_deltavec_cut(AUX_k_cut)^2 );
                             % the distance between He atoms, adjusted to account for delz \equiv z_rel
        AUX_jtilde = AUX_rj_soft / AUX_dr + 1;   % This is the (likely) fractional index of the
                                                 %  adjusted  r. 
                                              %  Below we will linearly interpolate VHeHe
                                              %  to this fractional index.
        AUX_jprime = floor(AUX_jtilde);       % nearest lowest integer index in VHeHe
        AUX_t = AUX_jtilde - AUX_jprime;      % used in the linear interpolation below
        if  AUX_jprime < AUX_numr
            AUX_V_k = r_V_95reduced(AUX_jprime, 2)*(1-AUX_t) + ...
                      r_V_95reduced(AUX_jprime+1, 2)*AUX_t;  % linearly interpolated value
        else     %  for rmax, don't need to interpolate, since VHeHe(rmax) is negligibly small 
            AUX_V_k = r_V_95reduced(AUX_numr, 2);
        end
        Vhehe_softened(AUX_jr) = Vhehe_softened(AUX_jr) + ...
                                 2*AUX_V_k * AUX_rho_rel_cut(AUX_k_cut) *AUX_dz;
                             % Recall that the "2*" appears because above we integrate
                             % only over z_rel > 0, and hence introduce the
                             % other "half" of rho_rel by this multiplication by 2. 
    end                     
    
end

% figure(73);
% plot(r_V_95reduced(:,1), r_V_95reduced(:,2) - Vhehe_softened')
% xlabel('r, [A]')
% ylabel('VheheORIG - VheheSOFT, [K]')
% ylim([-5 100])

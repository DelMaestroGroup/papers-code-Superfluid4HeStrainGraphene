% This code solves coupled Hartree-Fock eqs for the wave functions (WFs) psi{i}; 
% the equations are worked out in a separate pile (Derivation of Hartree-Fock for bosons,
% started 7/4/20), and the method's technical details are described in (Working out the
% AITEM for Hartree-Fock, started 7/6/20). 
% In short, equations consist of the Hartree part plus the exchange term (which is more
% expensive to compute than the Hartree part). 
% 
% (-del^2 + Vext + sum_{j~=i} int dr' VHeHe(r-r') |psi_j(r')|^2 ) psi_i   % <- Hartree
%   + sum_{j~=i} int dr' VHeHe(r-r') psi_j^*(r') psi_i^*(r') psi_j(r)     % <- Fock
%   - U <NinvU, U>^{-1} <Ninv U, L00u>  = 0                               % <- ensuring orthonormality
%
% This last term is worked out in Ref. [1] below. 
%
% The locations of the WFs in the graphene lattice are labeled along rows.
% In ..._v4 (compared to earlier versions) I automated the setup of an arbitrary
% (within reason) number of He atoms. 
%
% The organization of code follows that of spin1_BEC_ManyInitialConditions.m from my
% 2009-2010 Multi-component ITEM work [1].
%
% ********************************************************************************************
% Unlike in the "original" (= non-ME) version of this code, where I use two methods: 
% AITEM (simpler) or CGM (more complicated to code),
% here I use only the AITEM, to make the code simpler.
% ********************************************************************************************

%
% Reference:
% [1] My paper "Convergence conditions for iterative methods seeking multi-component solitary waves ...",
%     Mathematics and Computers in Simulation 81 (2011) 1572--1592.

clear all

savemyresults = 0;                    % 1 => save results in a file 
                                      % 0 => do not save results
if savemyresults == 1
    savedfilename = 'saveddata/22012501';
end

loadinitialus = 0;                    % 0 => use Gaussian-like us set up in the code;
                                      % 1 => load us from a file (this option is used when
                                      %       I had to stop a running file that had made progress
                                      %       toward the final solution, save its current variables, 
                                      %       and use those as the initial condition)
                                      
loadVext = 0;                         % 0 => use the model Vext generated inside the code
                                      % 1 => load a Vext provided by Jiang
                                      
softenVHeHe = 1;                      % 0 => use the original (Aziz) VHeHe
                                      % 1 => use VHeHe softened due to z-spread 


% ----------------------------------------------------------------------------------

% How often (after how many iterations) we display results and save data:
nplot = 200;                            % plot every nplot-th iteration 
nsave = 1000;


% - Set up the spatial domain:
% We use nondimensional units throughout; therefore, our VHeHe potential 
% (defined later) will also have to be in nondimensional units. 
uniform_strain = 0.3; 
a0 = 1.42* (1 + uniform_strain); %2.495/2.46;                            % side of graphene lattice [A]
% % % % The x- and y-periods, defined next, will be slightly adjusted later, 
% % % % to ensure that each period has an integer number of grid points:
% NO, now adjustment is needed: they are automatically defined to contain an integer number of pts...
per_x = a0*sqrt(3);                   
per_y = a0*3;                         % periods of Vext [A];
                                      % hence the dimensions of the grid must be multiple of these
space_scale = sqrt(3)*a0/pi;          % scale (in [A]) between nondimensional
                                      % and dimensional spacial variables
 
Numper_x = 6;
Numper_y = 3;                         % numbers of periods in x and y
% Total lengths along x and y, which will be adjusted slightly *soon*, to ensure that
% each period has an integer number of grid points:
xlength = Numper_x * per_x /space_scale; 
ylength = Numper_y * per_y /space_scale;  
Nx = 1/2* Numper_x/2 * 2^7;     
Ny = 1/2* Numper_y/2* 7/4* 2^7;                % number of points along x and y
                                          % (the factor 7/4 approximates per_y/per_x = sqrt(3))
dx = xlength/Nx;  dy = ylength/Ny;    % mesh sizes along x and y
% % % % Now adjust per_{x,y} to ensure that they contain an integer number of grid points,
% % % % and subsequently adjust {x,y}lengths:
% % % % per_x = round( (per_x/space_scale)/dx ) *dx*space_scale;
% % % % per_y = round( (per_y/space_scale)/dy ) *dy*space_scale;
% % % xlength = Numper_x * per_x /space_scale; 
% % % ylength = Numper_y * per_y /space_scale; 
x = [-xlength/2 :dx: xlength/2-dx];
y = [-ylength/2 :dy: ylength/2-dy];   % domains along x and y
[X,Y] = meshgrid(x,y);                % X is the Ny-by-Nx matrix that has x along each row, 
                                      %  and columns are all identical;
                                      % Y is the Ny-by-Nx matrix that has y along each column, 
                                      %  and rows are all identical.
dxdy = dx*dy;
      
xdim = x*space_scale;                 % [A] these are needed only for plotting
ydim = y*space_scale;


% - Define the spectral domain (used in convolution-type integration):
dkx = 2*pi/xlength;
dky = 2*pi/ylength;
kx = dkx*[0:Nx/2-1  -Nx/2:-1];
ky = dky*[0:Ny/2-1  -Ny/2:-1];       % k-domains along kx and ky
[KX,KY] = meshgrid(kx,ky);           % KX is the Ny-by-Nx matrix that has kx along each row, 
                                     %  and columns are identical;
                                     % KY is the Ny-by-Nx matrix that has ky along each column, 
                                     %  and rows are identical.
mfftdel2 = KX.^2 + KY.^2;            % (-del^2) in Fourier space 
                  

% - Set up the lattice vectors in direct and reciprocal space;
% notations x,y stand for the x,y-components in the rectangular coordinate system.
% Direct lattice vectors (nondimensional):
a1x = per_x/2 /space_scale;
a1y = per_y/2 /space_scale;  
a2x = - a1x;
a2y = a1y;
a1vec = [a1x; a1y];
a2vec = [a2x; a2y];

% Reciprocal lattice vectors (nondimensional):
gmat = 2*pi*inv([a1vec, a2vec].');
g1vec = gmat(:,1);
g2vec = gmat(:,2);
g3vec = g1vec + g2vec;               % auxiliary; we will be needed only for construction of Vext
                                     % Note: These notations of g{1,2,3} reciprocal lattice vectors
                                     %       have absolutely nothing to do with psi{1,2,3}.
g1x = g1vec(1);
g1y = g1vec(2);
g2x = g2vec(1);
g2y = g2vec(2);
g3x = g3vec(1);
g3y = g3vec(2);

g4vec = g1vec + g3vec;
g5vec = g2vec - g1vec;
g6vec = g2vec + g3vec;

g4x = g4vec(1);
g4y = g4vec(2);
g5x = g5vec(1);
g5y = g5vec(2);
g6x = g6vec(1);
g6y = g6vec(2);

% - Construct the external potential Vext:
xshift_Vext = 1*per_x/2 * rem(Numper_x * Numper_y, 2) /space_scale;
                                 % The minimum of Vext is shifted from 0 by this amount along x.
                                 % (The coefficient in front of per_x must be 1!
                                 %  If I set it to 0, the atoms are
                                 %  initially at the shallow saddle points
                                 %  on the ridges of graphene bonds.
                                 %  With this initial condition, the iterations still produce an error that
                                 %  decreases to ~10^-6 in the equation norm. However, then the
                                 %  iterations begin to diverge, albeit very slowly, and eventually 
                                 %  settle down to the solution with atoms at the wells of Vext.)
                                 % This shift is different than in the 1/3 and 2/3 codes,
                                 % for the reason explained before the loop where coordinates
                                 % of the centers of atoms are set up.
                                     
E_R = 9.89/(1 + 1*uniform_strain)^2;  % recoil energy [K], used as the energy scale 
                                     
if loadVext == 0                                     
%     Vext0 = 1*24/(E_R*4.5);                         % controls amplitude of Vext
%                                          % (Vext's full-swing amplitude = 4.5*Vext0)
%     Vext_nextharm = -1*Vext0*0.016;         % amplitude of the next harmonic
    
    Vext0 = 5.56/E_R;                         % controls amplitude of Vext
                                     % (Vext's full-swing amplitude = 4.5*Vext0)
    Vext_nextharm = -1.9e-2/E_R;         % amplitude of the next harmonic

    Xshifted_Vext = X - xshift_Vext;
    Vext = -Vext0*( cos( g1x*Xshifted_Vext + g1y*Y ) + ...
                    cos( g2x*Xshifted_Vext + g2y*Y ) + ...
                    cos( g3x*Xshifted_Vext + g3y*Y )) + ...
           -Vext_nextharm*( cos( g4x*Xshifted_Vext + g4y*Y ) + ...
                            cos( g5x*Xshifted_Vext + g5y*Y ) + ...
                            cos( g6x*Xshifted_Vext + g6y*Y ));
else % loadVext == 1
    V_Graphene_He;     % this loads VHeC over 2 x-periods and 1 y-period
    V_Graphene_He_x;   % this loads XforVHeC
    V_Graphene_He_y;   % this loads XforVHeC
    
    xsizeVHeC = size(VHeC,2);
    
    % Cut off the middle of the potential, to make it 1 x-period only:
    adj_nx = 0;
    VHeC_adj = VHeC(:, 1 + xsizeVHeC/4 : end - xsizeVHeC/4);
    XforVHeC_adj = XforVHeC(:, 1 + xsizeVHeC/4 : end - xsizeVHeC/4);
    YforVHeC_adj = YforVHeC(:, 1 + xsizeVHeC/4 : end - xsizeVHeC/4);
    
    % Interpolate one period of the loaded Vext onto a 1x1-period rectangle:
    xauxVext = -per_x /space_scale/2 : dx : per_x /space_scale/2 - dx;
    yauxVext = -per_y /space_scale/2 : dy : per_y /space_scale/2 - dy;
    [XauxVext, YauxVext] = meshgrid(xauxVext, yauxVext);
    VHeC_1x1period = interp2(YforVHeC_adj/2/space_scale, XforVHeC_adj/2/space_scale, ...
                             VHeC_adj, XauxVext, YauxVext);
    
    figure(200); mesh(XauxVext, YauxVext, VHeC_1x1period); view(2)
    
    % Create a tile of this VHeC:
    VHeC_tiled = repmat(VHeC_1x1period, [Numper_y, Numper_x]);
    
    figure(201); mesh(VHeC_tiled); view(2)
    
    % Now shift this VHeC by an amount defined above:
    nxshift_Vext = round(xshift_Vext/dx); 
    Vext(:, nxshift_Vext+1 : Nx) = VHeC_tiled(:, 1 : Nx - nxshift_Vext);
    Vext(:, 1 : nxshift_Vext) = VHeC_tiled(:, Nx - nxshift_Vext+1 : Nx);
    
end

    
Vext_ampl = max(max(Vext)) - min(min(Vext));
                                     % this will be used in defining preconditioning operators
% Visually examine Vext:
figure(301);
mesh(xdim,ydim,Vext); view(2)
xlabel('x'); ylabel('y');
% xlim([-3 3]); ylim([-3, 3])
axis('equal')
% 
% 
% jghjghjgjg

% - Load and prepare the He-He potential:
He_Aziz95_reduced;                   %  NOTES: 
                                     % 1:
                                     % The "code" He_Aziz95_reduced.m, which contains
                                     % r and V of the potetntial in dimensional units in the
                                     % variable called r_V_95reduced.
                                     % This is the newest Aziz potential (of 1995) provided by Adrian.
                                     % The original Adrian's data are in file He_Aziz95.dat;
                                     % there, the distance goes to 100A. In the reduced version
                                     % used here, the distance goes to 30A, which is still more
                                     % than sufficient for our purposes. 
                                     % 2:
                                     % Here we load the data of the original, i.e., NOT softened,
                                     % potential. The SOFTENED data, accounting for the z-spread of
                                     % the wavefunction, are computed a few lines below by running
                                     % the code       He_Aziz95_reduced_SOFTENED.
                                     %
% - Nondimensionalize r and V:
r_HeHeextendto = 100;                % [A] This is the distance to which I need to extend VHeHe
                                     % to make sure that it covers the entire computational domain.
rHeHe = [r_V_95reduced(:,1);  r_HeHeextendto] / space_scale;  
        % I have to use the last element in case that my numerical grid for rHeHe_arr
        % goes beyond the actual values of r stored in He_Aziz95_reduced.m.

if  softenVHeHe == 0
    VHeHe = [r_V_95reduced(:,2);  0] / E_R;
else  % if softenVHeHe = 1
    He_Aziz95_reduced_SOFTENED;
    VHeHe = [Vhehe_softened';  0] / E_R; 
end

% - Define the r and V, which originally are vectors, on a grid:
rHeHe_arr = sqrt(X.^2+Y.^2);      
rHeHe_vec = reshape(rHeHe_arr, [1 ,Ny*Nx]);
                               % r-vector corresponding to grid, *not* r from He_potential file
VHeHe_vec = interp1(rHeHe, VHeHe, rHeHe_vec);
                                            % interpolate VHeHe on the above r-vector
% % % % % VHeHe_arr = reshape(VHeHe_vec, [Ny, Nx]);
% % % % % mollify_VHeHe = 0*1; %1e-12;                  % set this to 0ish if you don't want to 
% % % % % fftVHeHe = mollify_VHeHe * fft2(VHeHe_arr);
% % % % % 
% % % % % 

scale_VHeHe = 1; %0.5e-3; %1e-12;                  % set this to 0 if there is no VHeHe
                                          % set this to < 1 to account for V_2particleWF < V_BH
VHeHe_arr = scale_VHeHe * reshape(VHeHe_vec, [Ny, Nx]);
fftVHeHe_true = fft2(VHeHe_arr);

% Cut off VHeHe at a certain height; otherwise the method will be extremely slow to converge.
Ecutoff = 2500;   % [K]            
Ecutoff_nondim = Ecutoff / E_R; 
VHeHe_arr_mollified = min(VHeHe_arr, Ecutoff_nondim);
fftVHeHe = fft2(VHeHe_arr_mollified);  



% % - Define iteration parameters:
% %
% % Constant in the preconditioning operator N (which is the same for all equations):
% cN = max(100, Vext_ampl*0.5);   % "0.5" appears to work noticeably better than "0.9"
% % The preconditioning operator (needed only in Fourier space):
% fftN_inv = 1./(cN + mfftdel2);
% fftN = 1./fftN_inv;
% 
% 
% Dt=0.5*50*0.01;   % Delta\tau (imaginary time step)
%                  % Linear case: 1000*0.0003 = 0.3 - converges well;
% dxdy = dx*dy;
% tol_D = 10^(-6+0) * 10^(-2)*Dt;                % tolerance up to which I perform computation


% - Set up initial guesses for psi's (below we will refer to psi's as u's).
%
NumHe_x = Numper_x;
NumHe_y = Numper_y*2;                   % number of He atoms per row and per column
NumHe = NumHe_x * NumHe_y;              % total number of He atoms


% This is done somewhat differently than in the 1/3 and 2/3 cases. 
% The difference is that we avoid putting the center of a cell at the corners,
% regardless of whether Numper_x is odd or even.
% Historically, this was done so as to avoid the then mysterious distortion, which I observed
% for atoms placed along the boundary. Such a distortion should be absent for periodic BC.
% Later, I found that the reason for that distortion was that I was using trapz instead of sum
% (which treats boundary points differently from those in the bulk).
% However, the way I set up He atoms is preserved. 
for nHe_y = 1 : NumHe_y
    for nHe_x = 1 : NumHe_x
        xcs( (nHe_y - 1)*NumHe_x + nHe_x ) = ( -Numper_x/2 + nHe_x - 1/2 + ...
                                               1/2*(rem(nHe_y, 2) - 1) ) * per_x/space_scale + ...
                                               0.1*randn;
        ycs( (nHe_y - 1)*NumHe_x + nHe_x ) = ( Numper_y/2 - (nHe_y - 1)/2 ) * per_y/space_scale + ...
                                             0.1*randn;
    end
end

if  loadinitialus == 0 

    %
    % Widths (make them all the same; I can't think of a situation where they are not):
    Wu = 0.4*a0/space_scale;

    pertsize = 0.0;
    for nHe = 1 : NumHe_x   % Set up initial guesses for atoms in the first (and last) half-row
                            % differently from the rest of the atoms.
        us(:,:, nHe) = exp( -( (X-xcs(nHe)).^2 + (Y-ycs(nHe)).^2 )/Wu^2 ) + ...
                       exp( -( (X-xcs(nHe)).^2 + (Y-ycs(nHe)+ylength).^2 )/Wu^2 );
                       % A Note on why I use "+ylength" above but "-xlength" later
                       % is found after I set up the last of us's. 
                       
        us(:,:, nHe) = us(:,:, nHe) + ...
                       pertsize*((2*randn-1)*(X-xcs(nHe)) + (2*randn-1)*(Y-ycs(nHe))).*...
                                exp( -( (1+pertsize*(2*randn-1))*(X-xcs(nHe)).^2 + ...
                                        (1+pertsize*(2*randn-1))*(Y-ycs(nHe)).^2 )/Wu^2 );
                                     
        us(:,:, nHe) = us(:,:, nHe) / sqrt(sum(sum( us(:,:, nHe).^2 ))*dxdy);
    end
    for nHe = NumHe_x + 1 : NumHe 
        if rem( floor((nHe-0.5)/NumHe_x), 2 ) == 0   % odd-numbered rows 
            us(:,:, nHe) = exp( -( (X-xcs(nHe)).^2 + (Y-ycs(nHe)).^2 )/Wu^2 );
        else  % i.e. if  rem( floor((nHe-0.5)/NemHe_x), 2 ) == 1   % even-numbered rows
            if  rem(nHe, NumHe_x) == 1   % first atom/cell in even-numbered row
                us(:,:, nHe) = exp( -( (X-xcs(nHe)).^2 + (Y-ycs(nHe)).^2 )/Wu^2 ) + ...
                               exp( -( (X-xcs(nHe)-xlength).^2 + (Y-ycs(nHe)).^2 )/Wu^2 );
            else
                us(:,:, nHe) = exp( -( (X-xcs(nHe)).^2 + (Y-ycs(nHe)).^2 )/Wu^2 );
            end
        end
        
        us(:,:, nHe) = us(:,:, nHe) + ...
                       pertsize*((2*randn-1)*(X-xcs(nHe)) + (2*randn-1)*(Y-ycs(nHe))).*...
                                exp( -( (1+pertsize*(2*randn-1))*(X-xcs(nHe)).^2 + ...
                                         (1+pertsize*(2*randn-1))*(Y-ycs(nHe)).^2 )/Wu^2 );
                                     
        us(:,:, nHe) = us(:,:, nHe) / sqrt(sum(sum( us(:,:, nHe).^2 ))*dxdy);
    end
    % Note on why I use "+ylength" for atoms of row 1 (top/bottom boundary)
    % but "-xlength" for atoms at the left/right boundary:
    %   - As one verifies explicitly, min(xcs) = min(x), while  max(ycs) = max(y+dy). 
    %     Thus for row 1-atoms, the 1st line in their setup creates "slightly incomplete"
    %     halves sitting at the top boundary (i.e. at max(y) rather than max(y+dy));
    %     then the 2nd line creates the remaining "complete" halves at the bottom boundary.
    %     On the contrary, for atoms sitting at the left/right boundary, the 
    %     1st line of their setup creates "complete" halves at the left boundary,
    %     while the 2nd line creates "slighty incomplete" halves at the right boundary. 
    
    
%     us(:,:,13) = us(:,:,13) + ...
%                  0.1*(X-xcs(13)).*exp( -1.1*( (X-xcs(13)).^2 + 0.8*(Y-ycs(13)).^2 )/Wu^2 );
%     us(:,:, 13) = us(:,:, 13) / sqrt(sum(sum( us(:,:, 13).^2 ))*dxdy);
    
elseif  loadinitialus == 1
    % Note that I still need xcs and ycs, defined above, because I save them at the end.
    
    load('saveddata/20072301', 'us');
    
end

% Xshifted_VextNew = X - xcs(1);
% Vext = -Vext0*( cos( g1x*Xshifted_VextNew + g1y*Y ) + ...
%                 cos( g2x*Xshifted_VextNew + g2y*Y ) + ...
%                 cos( g3x*Xshifted_VextNew + g3y*Y ) );
% 
% % Visually examine New Vext:
% figure(300);
% mesh(xdim,ydim,Vext); view(2)
% xlabel('x'); ylabel('y');
% % xlim([-3 3]); ylim([-3, 3])
% axis('equal')


                

        step4plot = 1;
        figure(3001);
        Numpan_x = NumHe_x;
        Numpan_y = NumHe_y/2;
        for npan_y = 1 : 2 : NumHe_y
            for npan_x = 1 : Numpan_x
                npan = (npan_y - 1)/2*Numpan_x + npan_x;
                n_of_u = (npan_y - 1)*Numpan_x + npan_x;
                %
                subplot(Numpan_y, Numpan_x, npan)
                mesh( xdim(1:step4plot:end), ydim(1:step4plot:end), ...
                      Vext(1:step4plot:end, 1:step4plot:end)/(5*Vext_ampl) ); view(2)
                hold on
                mesh( xdim(1:step4plot:end), ydim(1:step4plot:end), ...
                      us(1:step4plot:end, 1:step4plot:end, n_of_u) ); view(2)
                hold off
                xlabel('x'); ylabel('y');
                % xlim([-3 3]); ylim([-3, 3])
                xlim([-xlength/2, xlength/2]); 
                % ylim([-ylength/2, ylength/2]);
                axis('equal')
                title(['u' int2str(n_of_u)])
            end
        end
        
        figure(3002);
        for npan_y = 2 : 2 : NumHe_y
            for npan_x = 1 : Numpan_x
                npan = (npan_y - 2)/2*Numpan_x + npan_x;
                n_of_u = (npan_y - 2)*Numpan_x + Numpan_x + npan_x;
                %
                subplot(Numpan_y, Numpan_x, npan)
                mesh( xdim(1:step4plot:end), ydim(1:step4plot:end), ...
                      Vext(1:step4plot:end, 1:step4plot:end)/(5*Vext_ampl) ); view(2)
                hold on
                mesh( xdim(1:step4plot:end), ydim(1:step4plot:end), ...
                      us(1:step4plot:end, 1:step4plot:end, n_of_u) ); view(2)
                hold off
                xlabel('x'); ylabel('y');
                % xlim([-3 3]); ylim([-3, 3])
                xlim([-xlength/2, xlength/2]); 
                % ylim([-ylength/2, ylength/2]);
                axis('equal')
                title(['u' int2str(n_of_u)])
            end
        end
        %         
        
%         figure(3);
%          mesh( xdim(1:step4plot:end), ydim(1:step4plot:end), ...
%                       Vext(1:step4plot:end, 1:step4plot:end)/(5*Vext_ampl) ); view(2)
%                 hold on
%                 mesh( xdim(1:step4plot:end), ydim(1:step4plot:end), ...
%                       us(1:step4plot:end, 1:step4plot:end, 1) ); view(2)
%                 hold off
%                 xlabel('x'); ylabel('y');
%                 % xlim([-3 3]); ylim([-3, 3])
%                 xlim([-xlength/2, xlength/2]); 
%                 % ylim([-ylength/2, ylength/2]);
%                 axis('equal')
        

pause(1);


% ---------------------------------------------------------------------------
% Initialize auxiliary parameters for iterations:
%
% Constant in the preconditioning operator N (which is the same for all equations):
cN = max(100, Vext_ampl*0.5);   % "0.5" appears to work noticeably better than "0.9"

% The preconditioning operator (needed only in Fourier space):
fftN_inv = 1./(cN + mfftdel2);
fftN = 1./fftN_inv;

Dt=0.5*50*0.01;   % Delta\tau (imaginary time step)
                 % Linear case: 1000*0.0003 = 0.3 - converges well;
tol_D = 10^(-5+0) * 10^(-2)*Dt;                % tolerance up to which I perform computation

norm_D_eqs = 1;
norm_D = 1;                  % initialize the initial errors
counter = 0;                 % initialize iterations counter

norm_D_eqs_startME = 1e-22;  % "_ME" refers to Mode Elimination procedure, which has been
                             % found (for an unknown reason) not to work;
                             % hence I suppres its usage by setting a very low "starting error".
ME_start_counter = 50;      % counter with which ME starts (if it hasn't started earlier due to the
                             %  error becoming less than the above threshold)
gammaME_max = 10^4;          % max value of ME's parameter gamma
sME = 0.4;                   % parameter "s" of ME (how much of the slow mode we subtract)

FocFac = 1;                  % I may model just Hartree sometimes, to check if the code works properly
                             % w/o the Fock term. 
                             % 0 => Hartree only
                             % 1 => Hartree-Fock 


% Preallocate multi-dimensional arrays:
% usus = zeros(Ny, Nx, NumHe, NumHe);
fftusus  = zeros(Ny, Nx, NumHe, NumHe);
Fock = zeros(Ny, Nx, NumHe, NumHe);
Vext_VHeHe_on_us = zeros(Ny, Nx, NumHe);
L00s_Ha = zeros(Ny, Nx, NumHe);
L00s_Fo = zeros(Ny, Nx, NumHe);
L00s = zeros(Ny, Nx, NumHe);
Ninvus_cc = zeros(Ny, Nx, NumHe);
us_reshaped = zeros(Ny*Nx, NumHe);
aux_forUterm_1 = zeros(NumHe^2, NumHe^2);
aux_forUterm_2 = zeros(NumHe^2, 1);
L00minUterm = zeros(Ny, Nx, NumHe);
subtracted_mode = zeros(Ny, Nx, NumHe);


% Begin the iterations:
while norm_D > tol_D

    counter = counter+1;
%     for nHe = 1 : NumHe 
%         us_old(:,:, nHe) = us(:,:, nHe);
%     end	
    
    us_old = us;
    
    if counter > 1
        L00minUterm_old = L00minUterm;
    end
    
    
    % ***********************************************************
    %                          *** Compute new  us ***
    % ***********************************************************
            
    % -------------------------------------------------------------------------------------------
    % -------------------------------------------------------------------------------------------
    % Compute L0u = L00u - Uterm; this computation is common to the AITEM and CGM:
    
    % Step 1:   Compute the components of the vector  L^{00}u in (2.9) of [1]
    for mHe = 1 : NumHe 
        fftus(:,:, mHe) = fft2( us(:,:, mHe) );

        for  nHe = mHe : NumHe            
            % usus(:,:, mHe, nHe) = us(:,:, nHe).*conj(us(:,:, mHe));
            fftusus(:,:, mHe, nHe) = fft2( us(:,:, nHe).*conj(us(:,:, mHe)) );
        end
    end

    % This is the sum over all He atoms, sum_atoms fft(|u(atom)|^2), 
    % used to compute the Hartree contribution:
    fftusabssq_all = zeros(Ny,Nx);
    for  nHe = 1 : NumHe  
        fftusabssq_all = fftusabssq_all + fftusus(:,:, nHe, nHe);
    end

    % The following terms will be used to compute the Fock contributions:   
    for  mHe = 1 : NumHe
        for  nHe = 1 : mHe - 1
            Fock(:,:, mHe, nHe) = conj( Fock(:,:, nHe, mHe) );
        end
        for  nHe = mHe + 1 : NumHe   % this ensures that Fock(:,:, nHe,nHe) = 0
            Fock(:,:, mHe, nHe) = fftshift( ifft2( fftusus(:,:, mHe, nHe) .* fftVHeHe) )*dxdy;
        end
%         squeeze(max(max( abs(Fock(:,:, mHe, :)) )))
%         pause
    end
    % fordebug1 = zeros(NumHe, NumHe);
%     fordebug1 = squeeze(max(max( abs(Fock) )))
%     pause

    for nHe = 1 : NumHe
        Vext_VHeHe_on_us(:,:, nHe) = Vext + ...
                        fftshift( ifft2( (fftusabssq_all - fftusus(:,:, nHe, nHe) ) .* fftVHeHe) )*dxdy;
        % The Hartree terms (which include the 1-particle Hamiltonian):
        L00s_Ha(:,:, nHe) = ifft2(fftus(:,:, nHe) .* mfftdel2) + ...
                            Vext_VHeHe_on_us(:,:, nHe) .* us(:,:, nHe);
        % The Fock terms:
        L00s_Fo(:,:, nHe) = sum( Fock(:,:, :, nHe) .* us(:,:, :), 3);
        % Hartree + Fock:
        L00s(:,:, nHe) = L00s_Ha(:,:, nHe) + FocFac * L00s_Fo(:,:, nHe);
    end


    % Step2:  Compute the second term on the r.h.s. of (2.8a) and (2.9a) in [1] (see also (2.8c)).
    %         I do not compute the matrix {\mathcal U} directly, but instead use its representation
    %         using the Kronecker product at the end. Before that, I compute auxiliary matrices
    %         as described in my Notes on Hartree-Fock. 
    %
    for nHe = 1 : NumHe
        % (N^{-1} us)^* (Note that only the complex conjugate of this quantity is used below):
        Ninvus_cc(:,:, nHe) = conj( ifft2(fftus(:,:, nHe) .* fftN_inv) );   
    end

    % Matrix of <N^{-1} us,  us> and its inverse:
    for mHe = 1 : NumHe
        for nHe = 1 : mHe-1
            Ninvus_us(mHe, nHe) = conj( Ninvus_us(nHe, mHe) );
        end
        for nHe = mHe : NumHe
            Ninvus_us(mHe, nHe) = sum(sum( us(:,:, nHe).* Ninvus_cc(:,:, mHe) ))*dxdy; 
        end
    end

    Ninvus_us_inv = inv( Ninvus_us );

%     for nHe = 1 : NumHe
%         us_reshaped(:, nHe) = reshape( us(:,:, nHe), [Ny*Nx, 1] );
%     end
%     
%     aux_forUterm_1 = us_reshaped * Ninvus_us_inv;          % dimension = Ny*Nx x NumHe
%     aux_forUterm_2 = kron( aux_forUterm_1, eye(NumHe) );   % dimension = Ny*Nx*NumHe x NumHe^2

    aux_forUterm_1 = kron( Ninvus_us_inv, eye(NumHe) );   % dimension = NumHe^2 x NumHe^2 

    % <N^{-1} us, L00s>,  dimension = NumHe x NumHe
    for mHe = 1 : NumHe
        for nHe = 1 : NumHe
            Ninvus_L00s(mHe, nHe) = sum(sum( Ninvus_cc(:,:, nHe).* L00s(:,:, mHe) ))*dxdy; 
        end
    end

%     aux_forUterm_3 = aux_forUterm_2 * reshape( Ninvus_L00s, [NumHe^2, 1] );
%                      % dimension = Ny*Nx*NumHe x 1
    aux_forUterm_2 = aux_forUterm_1 * reshape( Ninvus_L00s, [NumHe^2, 1] );
                     % dimension = NumHe^2 x 1

    theUterm = zeros(Ny, Nx, NumHe);                 
    for nHe = 1 : NumHe    
        for  mHe = 1 : NumHe
           theUterm(:,:, nHe) = theUterm(:,:, nHe) + us(:,:, mHe) * aux_forUterm_2(nHe + (mHe-1)*NumHe);
        end
    end

%     theUterm( = reshape( aux_forUterm_3, [Ny, Nx, NumHe] );
%                      % U * <NinvU, U>^{-1} <NinvU, L00s> 

    if FocFac == 0
        for  nHe = 1 : NumHe
            Es(nHe) = Ninvus_L00s(nHe,nHe) / Ninvus_us(nHe,nHe);
            theUterm(:,:, nHe) = Es(nHe) * us(:,:, nHe);
        end
    end
        
        
    L00minUterm = L00s - theUterm;
	
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
    % The following is the Mode-Elimination step applied every iteration past some point
    % (implemented 3/24/22)
    if  norm_D_eqs < norm_D_eqs_startME || counter >= ME_start_counter
        
        for nHe = 1 : NumHe
            
            ffte_st_1 = fft2(e_st_1(:,:, nHe));    % Since norm_D_eqs_startME is sufficiently smaller than the
                                              % initial norm_D_eqs, and then several iterations will pass
                                              % before ME starts. Then, e_st_1  will already be defined. 
                                              % (It is defined after this loop.)
            L_e_st_1_sub = L00minUterm(:,:, nHe) - L00minUterm_old(:,:, nHe);
            e_st_1_M_e_st_1(nHe) = sum(sum( e_st_1(:,:, nHe) .*real(ifft2(fftN.*ffte_st_1)) ));
            alpha_numerator(nHe) = sum(sum( e_st_1(:,:, nHe) .*L_e_st_1_sub ));
            slowmode_numerator(nHe) = sum(sum( e_st_1(:,:, nHe) .*L00minUterm(:,:, nHe) ));
        end
        alpha_st_1 = sum(alpha_numerator)/sum(e_st_1_M_e_st_1);
        gamma_st_1_aux = 1 - sME/(alpha_st_1*Dt);
        gamma_st_1 = gamma_st_1_aux/sqrt(1+(gamma_st_1_aux/gammaME_max)^2);
        alpha_rec(counter) = alpha_st_1;
        for nHe = 1 : NumHe
            subtracted_mode(:,:, nHe) = e_st_1(:,:, nHe) *gamma_st_1 * ...
                                        sum(slowmode_numerator)/sum(e_st_1_M_e_st_1);
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         

    
   % ----------------------------------------------------------------------------------------

    %  Compute \hat{u}_{n+1} from (2.9a) in [1]:
    for nHe = 1 : NumHe
        us(:,:, nHe) = us(:,:, nHe) - Dt*real( ifft2( fftN_inv.* fft2( L00minUterm(:,:, nHe) ) ) ...
                                                - subtracted_mode(:,:, nHe) );
    % NOTE: The sign '-' instead of '+', as in the paper, in front of Dt is due to the fact
    %       that L00 (and hence the U-term) has the opposite sign as the corresponding 
    %       quantities in [1]. Recall that eigenvalues of L00 are negative in the paper
    %       but positive here. 
    end     
    
	
    % ----------------------------------------------------------------------------------------

    % Perform regular Gram-Schmidt on us.
    % Allow for the possibility to do so only every several steps (to save time):
    if norm_D_eqs > 1e-3
        do_GS = 1;
    elseif norm_D_eqs <= 1e-3 && norm_D_eqs > 1e-4
        do_GS = 1;
    else
        do_GS = 1;
    end

    if rem(counter, do_GS) == 0

        for  nHe = 1 : NumHe
            uaux = us(:,:, nHe);
            for  kHe = 1 : nHe-1
                uaux = uaux - (sum(sum( conj(us(:,:, kHe)).* uaux ))*dxdy) * us(:,:, kHe);
            end
            us(:,:, nHe) = uaux / sqrt( sum(sum( abs(uaux).^2 ))*dxdy ); 
        end

    end
    
    e_st_1 = us - us_old;
    
    % ----------------------------------------------------------------------------------------
    % ----------------------------------------------------------------------------------------
 
        
    
    % TO BE REMOVED AFTER EDITING THE ONE ABOVE
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
% % %     % The following is the Mode-Elimination acceleration step, which doesn't work
% % %     % (not due to a mistake in the code, but due to some other reason that I haven't
% % %     %  understood). 
% % %     if norm_D_eqs < norm_D_eqs_startME && do_ME == 0
% % %         do_ME = 1;
% % %         counter_start_ME = counter;
% % %     end
% % %     if  do_ME == 1 && rem(counter - counter_start_ME, timeperiod_ME) == 0
% % %         us_oldold(:,:, nHe) = us_old(:,:, nHe);
% % %     elseif  do_ME == 1 && rem(counter - counter_start_ME, timeperiod_ME) == 1
% % %         for nHe = 1 : NumHe
% % %             % Compute the decay factor of the slow mode for each u:
% % %             r_slowmode = sqrt( sum(sum( abs(us(:,:, nHe) - us_old(:,:, nHe)).^2 )) / ...
% % %                                sum(sum( abs(us_old(:,:, nHe) - us_oldold(:,:, nHe)).^2 )) );
% % %             % Subtract some part of the slow mode from the current solution:
% % %             us(:,:, nHe) = us(:,:, nHe) - ...
% % %                             0.7*r_slowmode/(r_slowmode - 1) * (us(:,:, nHe) - us_old(:,:, nHe));
% % %             % 
% % %         end
% % %     end
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         

    
    % For the purpose of stopping the iterations, compute the error,
    %  defined as the norm of the difference between the solutions at two consecutive steps.
    % First integrate (with trapz) over x & y, then sum over He atoms: 
    norm_D = 1/NumHe * sqrt( sum( sum(sum( (us - us_old).^2 ))*dxdy ) );
    rec_norm_D(counter) = norm_D;
	
	norm_D_eqs = 1/NumHe * sqrt( sum( sum(sum( abs(L00minUterm).^2 ))*dxdy ) );
	rec_norm_D_eqs(counter) = norm_D_eqs;
    
    
               
    % - Plot the solution and error:
    if rem(counter, nplot) == 0
        
%         figure(301);
%         for nHe = 1 : NumHe
%             mesh(xdim,ydim, (us(:,:, nHe) - us_old(:,:, nHe))/Dt); view(2)
%             hold on
%         end
%         hold off
%         xlabel('x'); ylabel('y');
%         xlim([-xlength/2, xlength/2]); ylim([-ylength/2, ylength/2]);
%         title(['normD = ' num2str(norm_D) ' , normDeqs = ' num2str(norm_D_eqs)])
%         pause(0.5)
         
		% step4plot = 2;
        figure(302);
        for nHe = 1 : NumHe
            mesh(xdim(1:step4plot:end),ydim(1:step4plot:end), us(1:step4plot:end,1:step4plot:end, nHe) ); view(2)
            hold on
        end
        mesh( xdim(1:step4plot:end),ydim(1:step4plot:end), Vext(1:step4plot:end,1:step4plot:end)/(3*Vext_ampl) ); view(2)
        hold off
        xlabel('x'); ylabel('y');
        xlim([-xlength/2, xlength/2]); ylim([-ylength/2, ylength/2]);
        title(['linear plot of us;  counter = ' int2str(counter) ])
         pause(0.5)
         
        figure(304);
        for nHe = 1 : NumHe
            mesh(xdim(1:step4plot:end),ydim(1:step4plot:end), log10(abs( us(1:step4plot:end,1:step4plot:end, nHe) )) ); view(2)
            hold on
        end
        hold off
        xlabel('x'); ylabel('y');
        xlim([-xlength/2, xlength/2]); ylim([-ylength/2, ylength/2]);
        title('log plot of us')

        figure(401);
        plot(log10(rec_norm_D/Dt),'b');
        hold on
        plot(log10(rec_norm_D_eqs),'r--');
        hold off
        legend('Du', 'Deqs')

        
        counter
        errs_Dt = [norm_D/Dt  norm_D_eqs  Dt]
        for nHe = 1 : NumHe
            Es_kin(nHe) = sum(sum( real( ifft2(fftus(:,:, nHe) .* mfftdel2) ).* us(:,:, nHe) ))*dxdy;
            Es_pot_relbot(nHe) = sum(sum( (Vext_VHeHe_on_us(:,:, nHe) - min(min(Vext))).* ...
                                          abs(us(:,:, nHe).^2) ))*dxdy;
            Es_pot_reltop(nHe) = sum(sum( (Vext_VHeHe_on_us(:,:, nHe) - max(max(Vext))).* ...
                                          abs(us(:,:, nHe).^2) ))*dxdy;
        end
        Es_kin_ave = mean(Es_kin);
        Es_pot_relbot_ave = mean(Es_pot_relbot);
        Es_pot_reltop_ave = mean(Es_pot_reltop);
        
%         Es_kin_arr_Kelvin = zeros(NumHe_y, NumHe_x);
%         Es_pot_arr_Kelvin = zeros(NumHe_y, NumHe_x);
        for nHe_y = 1 : NumHe_y
            for nHe_x = 1 : NumHe_x
                Es_kin_arr_Kelvin(nHe_y, nHe_x) = Es_kin( (nHe_y - 1)*NumHe_x + nHe_x ) * E_R; 
                Es_pot_relbot_arr_Kelvin(nHe_y, nHe_x) = ...
                                         Es_pot_relbot( (nHe_y - 1)*NumHe_x + nHe_x ) * E_R; 
                Es_pot_reltop_arr_Kelvin(nHe_y, nHe_x) = ...
                                         Es_pot_reltop( (nHe_y - 1)*NumHe_x + nHe_x ) * E_R;
            end
        end
        %  L00s_Fo_magnarr(nHe_y, nHe_x) = norm(L00s_Fo(:,:, (nHe_y - 1)*NumHe_x + nHe_x)); 
        d_Es_kin_arr_Kelvin = Es_kin_arr_Kelvin - Es_kin_ave*E_R;
        d_Es_pot_relbot_arr_Kelvin = Es_pot_relbot_arr_Kelvin - Es_pot_relbot_ave*E_R;
%         d_Es_pot_reltop_arr_Kelvin = Es_pot_reltop_arr_Kelvin - Es_pot_reltop_ave*E_R 
        d_Es___________kin__________pot_______arr_Kelvin = ...
            [d_Es_kin_arr_Kelvin  d_Es_pot_relbot_arr_Kelvin]
        E_kin_potrelbot_potreltop_ave = [Es_kin_ave Es_pot_relbot_ave Es_pot_reltop_ave]*E_R
        
        
        
               
        nHe_cent = NumHe_y/2 * NumHe_x + ceil(NumHe_x/2);   % pick some atom near the center of grid 
        
        Vext_VHeHe_on_uscent_true = Vext + ...
	                    fftshift( ifft2( (fftusabssq_all - fftusus(:,:, nHe_cent, nHe_cent) ) .* ...
                                         fftVHeHe_true) )*dxdy;

        Epot_cent = sum(sum( Vext_VHeHe_on_us(:,:, nHe_cent) .* abs(us(:,:, nHe_cent).^2) ))*dxdy;
        Epot_cent_true = sum(sum( Vext_VHeHe_on_uscent_true .* abs(us(:,:, nHe_cent).^2) ))*dxdy;        
        
        Epots_cent_Kelvin = [Epot_cent  Epot_cent_true] * E_R
        
        pause(0.5)
    end
    
    
    
    % - Record the solution periodically so as not to lose it due to unforeseen circumstances:
    if rem(counter, nsave) == 0
        if savemyresults == 1
            save(savedfilename, 'x', 'y', 'space_scale', 'a0', 'per_x', 'per_y', 'E_R', ...
                 'a1vec', 'a2vec', 'g1vec', 'g2vec', ... 
                 'Numper_x', 'Numper_y', 'NumHe_x', 'NumHe_y', ...
                 'xdim', 'ydim', 'xlength', 'ylength', 'Vext_ampl', ... % these are needed only for plotting
                 'Vext', 'softenVHeHe', 'VHeHe_arr', 'Ecutoff', 'scale_VHeHe', ... 
                 'cN', 'Dt', 'tol_D', 'rec_norm_D', 'rec_norm_D_eqs', ...
                 'xcs', 'ycs', 'us', ...
                 'Es_kin_arr_Kelvin', 'Es_pot_relbot_arr_Kelvin', 'Es_pot_reltop_arr_Kelvin');
        end
    end


        
end     


step4plot = 2;
figure(4001);
Numpan_x = NumHe_x;
Numpan_y = NumHe_y/2;
for npan_y = 1 : 2 : NumHe_y
    for npan_x = 1 : Numpan_x
        npan = (npan_y - 1)/2*Numpan_x + npan_x;
        n_of_u = (npan_y - 1)*Numpan_x + npan_x;
        %
        subplot(Numpan_y, Numpan_x, npan)
        mesh( xdim(1:step4plot:end), ydim(1:step4plot:end), ...
              Vext(1:step4plot:end, 1:step4plot:end)/(5*Vext_ampl) ); view(2)
        hold on
        mesh( xdim(1:step4plot:end), ydim(1:step4plot:end), ...
              us(1:step4plot:end, 1:step4plot:end, n_of_u) ); view(2)
        hold off
        xlabel('x'); ylabel('y');
        % xlim([-3 3]); ylim([-3, 3])
        xlim([-xlength/2, xlength/2]);  ylim([-ylength/2, ylength/2]);
%         axis('equal')
        title(['u' int2str(n_of_u)])
    end
end

figure(4002);
for npan_y = 2 : 2 : NumHe_y
    for npan_x = 1 : Numpan_x
        npan = (npan_y - 2)/2*Numpan_x + npan_x;
        n_of_u = (npan_y - 2)*Numpan_x + Numpan_x + npan_x;
        %
        subplot(Numpan_y, Numpan_x, npan)
        mesh( xdim(1:step4plot:end), ydim(1:step4plot:end), ...
              Vext(1:step4plot:end, 1:step4plot:end)/(5*Vext_ampl) ); view(2)
        hold on
        mesh( xdim(1:step4plot:end), ydim(1:step4plot:end), ...
              us(1:step4plot:end, 1:step4plot:end, n_of_u) ); view(2)
        hold off
        xlabel('x'); ylabel('y');
        % xlim([-3 3]); ylim([-3, 3])
        xlim([-xlength/2, xlength/2]);  ylim([-ylength/2, ylength/2]);
%         axis('equal')
        title(['u' int2str(n_of_u)])
    end
end
%         

        
figure(401);
plot(log10(rec_norm_D/Dt),'b');
hold on
plot(log10(rec_norm_D_eqs),'r--');
hold off
legend('Du', 'Deqs')


        
% Note: I started to record them in Kelvin only starting with 20061101.
% 
V_int_totperatom = (Epot_cent_true - sum(sum(Vext .* us(:,:, nHe_cent).^2 ))*dxdy ) /6 * E_R
      % This name was introduced on 2022-05-19; previously it was V_BH (incorrectly so). 
      % This is te total energy with *all* atoms due to their nonlinear interaction, divided by 
      % the number of "equivalent" interacting atoms in each "layer" of neighbors.
      % This number is 6, and hence " / 6" above. 
      % Note that *not only* the nearest neighbors, but all surrounding atoms are included.
      
% J_BC = 2 * trapz(trapz( (u32.^3) .* fftshift( ifft2( fftu22 .* fftVHeHe_true) ) )) *(dxdy)^2 * E_R
%       % The "2*" is for the same reason as above. 
      
t_BH = sum(sum( ( real( ifft2(fftus(:,:, nHe_cent).*mfftdel2) ) + Vext.*us(:,:, nHe_cent) ) .* ...
                     us(:,:, nHe_cent-NumHe_x+1) ))*dxdy * E_R 
        

figure(4102);
plot(xdim, log10(abs(us(Ny/2+1,:, nHe_cent))), 'b', 'linewidth', 2)
xlabel('x, [A]', 'fontsize', 16); ylabel('log_{10} u_{ctr}', 'fontsize', 16);
xlim([min(xdim) max(xdim)]);
set(gca, 'fontsize', 20)


% This line allows one to compute various V-like terms using saved data:
V1_BH = trapz(trapz( (us(:,:,nHe_cent).*us(:,:,nHe_cent)) .* ...
                      fftshift( ifft2( fft2(us(:,:,nHe_cent+1).*us(:,:,nHe_cent+1) ) .* ...
                      fft2(VHeHe_arr) ) ) )) ...
          *(xlength/length(x))^2 * (ylength/length(y))^2 * E_R
      % <n,n |Vint|(n+1),(n+1)>

V2_BH =  trapz(trapz( (us(:,:,nHe_cent+1).*us(:,:,nHe_cent+1)) .* ...
                       fftshift( ifft2( fft2(us(:,:,nHe_cent+1).*us(:,:,nHe_cent+2) ) .* ...
                       fft2(VHeHe_arr) ) ) )) ...
          *(xlength/length(x))^2 * (ylength/length(y))^2 * E_R
      % <(n+1),(n+1) |Vint|(n+1),(n+2)>
      
V3_BH =  trapz(trapz( (us(:,:,nHe_cent+2).*us(:,:,nHe_cent+2)) .* ...
                       fftshift( ifft2( fft2(us(:,:,nHe_cent+2).*us(:,:,nHe_cent+1) ) .* ...
                       fft2(VHeHe_arr) ) ) )) ...
          *(xlength/length(x))^2 * (ylength/length(y))^2 * E_R
      % <(n+2),(n+2) |Vint|(n+1),(n+2)>
      
V4_BH =  trapz(trapz( (us(:,:,nHe_cent-1).*us(:,:,nHe_cent-1)) .* ...
                       fftshift( ifft2( fft2(us(:,:,nHe_cent-1).*us(:,:,nHe_cent) ) .* ...
                       fft2(VHeHe_arr) ) ) )) ...
          *(xlength/length(x))^2 * (ylength/length(y))^2 * E_R
      % <(n-1),(n-1) |Vint|(n-1),n>
            
V5_BH =  trapz(trapz( (us(:,:,nHe_cent-1).*us(:,:,nHe_cent-1)) .* ...
                         fftshift( ifft2( fft2(us(:,:,nHe_cent).*us(:,:,nHe_cent) ) .* ...
                         fft2(VHeHe_arr) ) ) )) ...
          *(xlength/length(x))^2 * (ylength/length(y))^2 * E_R
      % <(n-1),(n-1) |Vint|n,n>
      
V6_BH =  trapz(trapz( (us(:,:,nHe_cent).*us(:,:,nHe_cent)) .* ...
                         fftshift( ifft2( fft2(us(:,:,nHe_cent-1).*us(:,:,nHe_cent-1) ) .* ...
                         fft2(VHeHe_arr) ) ) )) ...
          *(xlength/length(x))^2 * (ylength/length(y))^2 * E_R
      % <n,n |Vint| (n-1),(n-1)>

      Vprime = (V_int_totperatom-V1_BH)/1.63

 if savemyresults == 1
    save(savedfilename, 'x', 'y', 'Nx', 'Ny', 'space_scale', 'a0', 'per_x', 'per_y', 'E_R', ...
         'a1vec', 'a2vec', 'g1vec', 'g2vec', ... 
         'Numper_x', 'Numper_y', 'NumHe_x', 'NumHe_y', ...
         'xdim', 'ydim', 'xlength', 'ylength', 'Vext_ampl', ... % these are needed only for plotting
         'Vext', 'softenVHeHe', 'VHeHe_arr', 'Ecutoff', 'scale_VHeHe', ... 
         'cN', 'Dt', 'tol_D', 'rec_norm_D', 'rec_norm_D_eqs', ...
         'xcs', 'ycs', 'us', ...
         'Es_kin_arr_Kelvin', 'Es_pot_relbot_arr_Kelvin', 'Es_pot_reltop_arr_Kelvin', ...
         'V_int_totperatom', 't_BH', 'V1_BH', 'V2_BH', 'V3_BH', 'V4_BH', 'V5_BH', 'V6_BH', 'Vprime');
 end


 
 ufict = zeros(size(X));
 numnodes_per_x = round( (per_x/space_scale) / dx );
 for nnx = 1 : Nx - numnodes_per_x 
     ufict(:, nnx) = us(:, nnx + numnodes_per_x, 8);
 end
 figure(501);
    mesh( xdim,ydim, Vext/(5*Vext_ampl) ); view(2)
    hold on
    mesh( xdim,ydim, ufict ); view(2)
    hold off
    xlabel('x'); ylabel('y');
    % xlim([-3 3]); ylim([-3, 3])
    xlim([-xlength/2, xlength/2]); % ylim([-ylength/2, ylength/2]);
    axis('equal')
    title('u fictitious')
 
Vxxx_BH =  trapz(trapz( (us(:,:,nHe_cent-1).*us(:,:,nHe_cent-1)) .* ...
                         fftshift( ifft2( fft2(us(:,:,nHe_cent-1).*ufict ) .* ...
                         fft2(VHeHe_arr) ) ) )) ...
          *(xlength/length(x))^2 * (ylength/length(y))^2 * E_R
      
# Test-1- 

format long
%%%%%%%%%% Modify the following %%%%%%%%%%%%
E1 = 2.77e7; 
E2 = 1.44e6;
G12 = 1.13e6; 
v12 = 0.35;
v21 = v12 * E2 / E1 ;
a1 = -5e-7 ;
a2 = -1.24e-5;

deltaT = -300;%adjusted per question 
layup = [90 0 -45 45 90 90 0 -45];
Nplies = 8; 
ply_thickness =  [0.0052 0.0052 0.0052 0.0052 0.0052 0.0052 0.0052 0.0052];
tply = 0.0052;   
Nx =100; Ny=50; Nxy=-50;%%%adjusted per question 
Mx=0; My=0; Mxy=3; 

zcoords1 = [-0.0052*4 -0.0052*3 -0.0052*2 -0.0052  0 0.0052 0.0052*2 0.0052*3 0.0052*4]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% z-location ply to ply 
h = Nplies * tply ;
for i = 1:3; %%% adjust to calculate each ply ( for bottom ply 3 look at bottom matrcies in sigma/epsilon
  zcoords(i) = - (h + tply)/2 + i*tply;
end

% Q matrix (material coordinates)
Q11 = E1 / 1-v12*v21        ;
Q12 = v12 * E2 / 1-v12*v21 ;
Q22 = E2 / 1-v12*v21;
Q66 = G12;

Q = [ Q11 Q12 0; Q12 Q22 0; 0 0 Q66] ;

a = [a1 a2 0] ;
S= inv(Q) ;

% Qbar matrices (laminate coordinates) and contributions to
% ABD matrices

A = zeros(3,3);
B = zeros(3,3);
D = zeros(3,3);
NT = zeros(3,1);
MT = zeros(3,1);
alpha = zeros(3,3);
beta = zeros(3,3);
delta = zeros(3,3);
%for loop to iterate over plies
for i = 1:Nplies
   % Qbar matrices eval 
  theta  = layup(i);
  c = cosd(theta) ;
  s = sind(theta) ;
  T = [ c^2 s^2 2*c*s; s^2 c^2 -2*c*s; -c*s c*s (c^2 - s^2)];
  Qbar = inv(T) * Q * (inv(T))';
  Sbar = inv(Qbar);
  abar = T' .* a;
  % ABD matrix 
  A = A + Qbar * tply;
  B = B + Qbar * tply * zcoords(i); 
  D = D + Qbar * (tply * zcoords(i)^2  + tply^3 / 12);
  alpha = alpha + Qbar * ply_thickness(i)^2 /6;
  beta = alpha + Sbar * ply_thickness(i)^2 /6;
  delta = delta + 2 .* beta + Sbar.*(abar.^2).*tply^2 / 6;;  
  NT ;
  Qbar;
  abar;
  tply;
  deltaT;
  %thermal load added to calculated mechanical load 
  NT = NT + Qbar * abar * tply * deltaT ;
  MT = MT + Qbar * abar * tply * zcoords(i) * deltaT; 
  %midplane strains and curvatures 
  eo = A \ (NT + B * [Nx; Ny; Nxy] + D * [Mx; My; Mxy]);
  ko = alpha \ (MT + beta * [Nx; Ny; Nxy] + delta * [Mx; My; Mxy]); 
  
  %local lamina stress and strain
  e_top = eo + [0; 0; zcoords(i)];
  e_bot = eo + [0; 0; zcoords(i) - tply];
  k_top = ko + [0; 0; zcoords(i)];
  k_bot = ko + [0; 0; zcoords(i) - tply];
  e = T * [e_top e_bot];
  k = T * [k_top k_bot];

  sigma = Q* e + a* k;
  epsilon = S * sigma - a* k; 
  
    %thermal load
  NT = NT + Qbar * abar * tply * deltaT ;
  MT = MT + Qbar * abar * tply * zcoords(i) * deltaT; 

ABD    = [A B; B D];
ABDinv = inv(ABD) ;

e0k = ABDinv * [NT' MT']';
kx = e0k(1:3,1)\inv(A);
end 

Qbar
Sbar
abar
A
B
D
alpha 
beta 
delta
N_total = NT(:,1)+NT(:,2)
M_total = MT(:,1)+MT(:,2)
e0k = e0k(1:3,1) 
kx
%
sigma'
epsilon'

max_strain = [0.564 0.564  0.564  0.564  0.564  0.564  0.564  0.564 0.564 ]
ply_thickness =  [0.0052 0.0052 0.0052 0.0052 0.0052 0.0052 0.0052 0.0052]
X = 4.7e5
Xp = -2.3e5
Y = 8.9e3
Yp = -2.9e4
S = 1.09e4

x_eps = X / (X + abs(Xp));
xp_eps = Xp / (X + abs(Xp));
y_eps = Y / (Y + abs(Yp));
yp_eps = Yp / (Y + abs(Yp));
s_eps = S / (S + abs(S)); 
fail_stress = Inf;
fail_eps = Inf;
fail_ply = 0;
for i = 1:9
    if abs(max_strain(i)) > abs(fail_eps)
        fail_eps = max_strain(i);
        fail_ply = i;
    end
    if max_strain(i) > 0 && max_strain(i) < x_eps
        fail_stress = X;
        fail_eps = max_strain(i);
        fail_ply = i;
    end
    if max_strain(i) < 0 && abs(max_strain(i)) < xp_eps
        fail_stress = Xp;
        fail_eps = max_strain(i);
        fail_ply = i;
    end
    if max_strain(i) > 0 && abs(max_strain(i)) < y_eps
        fail_stress = Y;
        fail_eps = max_strain(i);
        fail_ply = i;
    end
    if max_strain(i) < 0 && abs(max_strain(i)) < yp_eps
        fail_stress = Yp;
        fail_eps = max_strain(i);
        fail_ply = i;
    end
    if abs(max_strain(i)) < s_eps
        fail_stress = S;
        fail_eps = max_strain(i);
        fail_ply = i;
    end
end 
fail_stress 
fail_eps 
fail_ply




%%%%%work for mc%%% 

% E1 = 20.6e6
% E2 = 1.5e6
% G12 = 1e6
% v12 = 0.27
% a1 = -0.5e-6 
% a2 = 15e-6 
% theta = 90
% dT = -200 
% 
% Q1=  [20.7e6, 0.407e6, 0;
% 0, 1.51e6 ,0 ;
% 0, 0 1e6]  
% T = [cosd(theta)^2, sind(theta)^2, 2 * cosd(theta) * sind(theta);
%      sind(theta)^2, cosd(theta)^2, -2 * cosd(theta) * sind(theta);
%      -cosd(theta) * sind(theta), cosd(theta) * sind(theta), cosd(theta)^2 - sind(theta)^2];
% 
% 
% Q_global = T * Q1 * T';
% 
%  ex = 0.002; ey = 0.001; exy = 0.005; 
%  strain_global = [ex ey exy]
%  strain_local = T.*strain_global
%  local_stress = Q_global.*strain_local




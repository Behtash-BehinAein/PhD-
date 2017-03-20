% Behtash Behin-Aein
clear; clc; 
%
%%%%%%%%%   Constants [MKS]   %%%%%%%%%%%
mue=9.2741e-24;            %Bohr Magneton
hbar=1.055e-34;    %Planck's Constant/2pi
gs=-2.00;                %Landau g factor
g=gs*mue/hbar;        %Gyromagnetic Ratio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%% MKS to CGS %%%%
g=g/1e4;
%%%%%%%%%%%%%%%%%%%
%
% Parameters   [CGS] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 4;                               % Total # of magnets
a =  30e-7;              % Spacing between the magents
Nt1= 20000;%#of time grid points for one bit propagation 
Nt = Nt1*(N-1);            % # of points in the time grid
t=linspace(0,50e-9,Nt);dt=t(2)- t(1);          %Time grid  
alpha = .1;                               %Damping Factor
Ku2 = 3e5;  %Second order uniaxial anisotropy constant
V0 = 10*10*10*1e-21;                              %[cm^3]
Ku2 = Ku2*V0;                                     % [erg]
Ms= 800 * V0;    %Saturation magnetization moment of NiFe
Kp = 0*pi*Ms^2;             %In plane anisotropy constant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Spatial Coordinates   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[rx(1:N,1:N) ry(1:N,1:N) rz(1:N,1:N) r(1:N,1:N) A(1:N,1:N)]...
                 = coordinate_generator(N,a);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
rtv = 2; 
rise_time = rtv*dt/1e-9;
E_disskT = zeros(1,length(rtv));
for  rr=1:length(rtv)
%
% Initialization   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MN_Initial = Initial_Magnetization1DArray(N);
% x,y and z components  
mx = MN_Initial(1:3:3*N-2)*Ms;
my = MN_Initial(2:3:3*N-1)*Ms;
mz = MN_Initial(3:3:3*N)*Ms;
m  = MN_Initial*Ms;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Applied magnetic field %%%%%%%%%%%%%%%
% Storage arrays  %%%
Hx = zeros(N,Nt);
Hy = zeros(N,Nt);
Hz = zeros(N,Nt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%rr
rt =rtv(rr);
rt1=rt-1;
cc=1;
for nn = 1 : N-2
    %%%%%%%%%%%%%
    h0=1*750;
    theta_H=pi/2;
    phi_H=0;
    %%%%%%%%%%%%%        
    Hxd = h0 * sin(theta_H) *cos(phi_H);
    Hyd = h0 * sin(theta_H) *sin(phi_H);
    Hzd = 0*h0 * cos(theta_H);                       
    ps = (cc-1)*Nt1+1000;
    pe = (cc-1)*Nt1+4000+2*rt;
    %%%
    Hx(nn,ps:ps+rt1) = linspace(0,Hxd,rt);
    Hx(nn,ps+rt:pe-rt) = Hxd;
    Hx(nn,pe-rt1:pe) = linspace(Hxd,0,rt);
    %%%
    Hy(nn,ps:ps+rt1) = linspace(0,Hyd,rt);
    Hy(nn,ps+rt:pe-rt) = Hyd;
    Hy(nn,pe-rt1:pe) = linspace(Hyd,0,rt);    
    %%%
    Hz(nn,ps:ps+rt1) = linspace(0,Hzd,rt);
    Hz(nn,ps+rt:pe-rt) = Hzd;
    Hz(nn,pe-rt1:pe) = linspace(Hzd,0,rt);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Hx(nn+1,ps:ps+rt1) = linspace(0,Hxd,rt);
    Hx(nn+1,ps+rt:pe-rt) = Hxd;
    Hx(nn+1,pe-rt1:pe) = linspace(Hxd,0,rt);
    %%%
    Hy(nn+1,ps:ps+rt1) = linspace(0,Hyd,rt);
    Hy(nn+1,ps+rt:pe-rt) = Hyd;
    Hy(nn+1,pe-rt1:pe) = linspace(Hyd,0,rt);
    %%%
    Hz(nn+1,ps:ps+rt1) = linspace(0,Hzd,rt);
    Hz(nn+1,ps+rt:pe-rt) = Hzd;
    Hz(nn+1,pe-rt1:pe) = linspace(Hzd,0,rt);
    %%%
    cc=cc+1;            
end            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%
%%%%%%%
%%%%%%%
%%%%%%%
%%%%%%%   Initialization and Storage Arrays      %%%%%%%%%
Mx = zeros(N,Nt-1); 
Mx(:,1) = mx;
My = zeros(N,Nt-1);
My(:,1) = my;
Mz = zeros(N,Nt-1);
Mz(:,1) = mz;
%%%%%%%%%%%%%%%%%%%%
Wt = zeros(3*N,3*N);
%%%%%%%%%%%%%%%%%%%%
Ht = zeros(3*N,1);
%%%%%%%%%%%%%%%%%%%%
mx_matrix = zeros(N,N);
my_matrix = zeros(N,N);
mz_matrix = zeros(N,N);
Mux_matrix = zeros(N,N);
Muy_matrix = zeros(N,N);
Muz_matrix = zeros(N,N);
%%%%%%%%%%%%%%%%%%%%
Hanix = zeros(N,Nt);
Haniy = zeros(N,Nt);
Haniz = zeros(N,Nt);
%%%%%%%%%%%%%%%%%%%%
Hanix(:,1) = zeros(N,1);
Haniy(:,1) = 2 * 1 * Ku2*my/Ms^2;  
Haniz(:,1) = - 2 * 1 * Kp*mz/Ms^2;
%%%%%%%%%%%%%%%%%%%%
                  %%%%%%%%%%%%%%%%%%%%    
                  Hdipx = zeros(N,Nt);
                  Hdipy = zeros(N,Nt);
                  Hdipz = zeros(N,Nt);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  for nn = 1 : N
                      mx_matrix(:,nn) = mx; 
                      my_matrix(:,nn) = my; 
                      mz_matrix(:,nn) = mz; 
                  end
                  mx_matrix = mx_matrix(A);
                  my_matrix = my_matrix(A);
                  mz_matrix = mz_matrix(A);   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
    %%%%%                Dipolar Field               %%%%%
    %%%           H = (3(m.r)*r - m*r^2)/r^5           %%%
    mdotr = mx_matrix.*rx + my_matrix.*ry + mz_matrix.*rz;
    Hx_Dipole_mdotr = mdotr .*rx ./ (r.^5);
    Hy_Dipole_mdotr = mdotr .*ry ./ (r.^5);
    Hz_Dipole_mdotr = mdotr .*rz ./ (r.^5);
    %%%
    Hx_Dipole_m = mx_matrix ./ (r.^3);
    Hy_Dipole_m = my_matrix ./ (r.^3);
    Hz_Dipole_m = mz_matrix ./ (r.^3);
    %%%
    Hx_Dipole_mdotr(1,:) = zeros(1,N);
    Hx_Dipole_m(1,:) = zeros(1,N);
    Hx_Dip = 3*Hx_Dipole_mdotr - Hx_Dipole_m;
%  Hx_Dip(2,2)=0;
%Hx_Dip(3,3)=0;
%Hx_Dip(3,2)=0;%
%Hx_Dip(2,3)=0;%Hx_Dip(4,2)=0;
    Hx_Dip_Collective = sum(Hx_Dip,1)';
    %
    Hy_Dipole_mdotr(1,:) = zeros(1,N);
    Hy_Dipole_m(1,:) = zeros(1,N);    
    Hy_Dip = 3*Hy_Dipole_mdotr - Hy_Dipole_m;
%  Hy_Dip(2,2)=0;
%Hy_Dip(3,3)=0;
%Hy_Dip(3,2)=0;%
%Hy_Dip(2,3)=0;%Hy_Dip(4,2)=0;
    Hy_Dip_Collective = sum(Hy_Dip,1)';
    %
    Hz_Dipole_mdotr(1,:) = zeros(1,N);
    Hz_Dipole_m(1,:) = zeros(1,N);    
    Hz_Dip = 3*Hz_Dipole_mdotr - Hz_Dipole_m;
%  Hz_Dip(2,2)=0;
%Hz_Dip(3,3)=0;
%Hz_Dip(3,2)=0;%
%Hz_Dip(2,3)=0;%Hz_Dip(4,2)=0;
    Hz_Dip_Collective = sum(Hz_Dip,1)';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             Hdipx(:,1) = Hx_Dip_Collective;
             Hdipy(:,1) = Hy_Dip_Collective;
             Hdipz(:,1) = Hz_Dip_Collective;    
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%   End of initialization    %%%%%%%%%%%%%%%%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             
% Time dependent calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for tt=2:Nt-1

    for nn = 1 : N
        mx_matrix(:,nn) = mx; 
        my_matrix(:,nn) = my; 
        mz_matrix(:,nn) = mz; 
    end
    mx_matrix = mx_matrix(A);
    my_matrix = my_matrix(A);
    mz_matrix = mz_matrix(A);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%                Dipolar Field               %%%%%
    %%%           H = (3(m.r)*r - m*r^2)/r^5           %%%
    mdotr = mx_matrix.*rx + my_matrix.*ry + mz_matrix.*rz;
    Hx_Dipole_mdotr = mdotr .*rx ./ (r.^5);
    Hy_Dipole_mdotr = mdotr .*ry ./ (r.^5);
    Hz_Dipole_mdotr = mdotr .*rz ./ (r.^5);
    %%%
    Hx_Dipole_m = mx_matrix ./ (r.^3);
    Hy_Dipole_m = my_matrix ./ (r.^3);    
    Hz_Dipole_m = mz_matrix ./ (r.^3);
    %%%
    Hx_Dipole_mdotr(1,:) = zeros(1,N);
    Hx_Dipole_m(1,:) = zeros(1,N);
    Hx_Dip = 3*Hx_Dipole_mdotr - Hx_Dipole_m;
%  Hx_Dip(2,2)=0;
%Hx_Dip(3,3)=0;
%Hx_Dip(3,2)=0;%
%Hx_Dip(2,3)=0;%Hx_Dip(4,2)=0;
    Hx_Dip_Collective = sum(Hx_Dip,1)';
    %
    Hy_Dipole_mdotr(1,:) = zeros(1,N);
    Hy_Dipole_m(1,:) = zeros(1,N);    
    Hy_Dip = 3*Hy_Dipole_mdotr - Hy_Dipole_m;
%  Hy_Dip(2,2)=0;
%Hy_Dip(3,3)=0;
%Hy_Dip(3,2)=0;%
%Hy_Dip(2,3)=0;%Hy_Dip(4,2)=0;
    Hy_Dip_Collective = sum(Hy_Dip,1)';
    %
    Hz_Dipole_mdotr(1,:) = zeros(1,N);
    Hz_Dipole_m(1,:) = zeros(1,N);    
    Hz_Dip = 3*Hz_Dipole_mdotr - Hz_Dipole_m;
%  Hz_Dip(2,2)=0;
%Hz_Dip(3,3)=0;
%Hz_Dip(3,2)=0;%
%Hz_Dip(2,3)=0;%Hz_Dip(4,2)=0;
    Hz_Dip_Collective = sum(Hz_Dip,1)';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %%%%% Uniaxial Anisotropy %%%%%%%%%%%%%%
    Huax = zeros(N,1);
    Huay = 2 * 1 * Ku2*my/Ms^2;  
    Huaz = - 2 * 1 * Kp*mz/Ms^2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Huay(1) = 0;  
    %Huax(1) = 2 * 1 * Ku2*mx(1)/Ms^2;
    %
    %%%%%   Overall field including anisotropy and the     %%%%%%%%
    %%%%%        field due to all other dipoles            %%%%%%%%
    HxNt = Hx(:,tt) + 1*Huax + 1*Hx_Dip_Collective;
    HyNt = Hy(:,tt) + 1*Huay + 1*Hy_Dip_Collective;
    HzNt = Hz(:,tt) + 1*Huaz + 1*Hz_Dip_Collective;       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %%% This is W of the precession term %%%%%%%%%%%%%%%%%%%%%%%%%%
    upper_diagonal_1 = zeros(3*N-1,1);
    upper_diagonal_2 = zeros(3*N-2,1);
    lower_diagonal_1 = zeros(3*N-1,1);
    lower_diagonal_2 = zeros(3*N-2,1);
    %%%
    %    [0     Hz  -Hy  %%%%
    %    -Hz    0    Hx  %%%%
    %     Hy   -Hx   0]  %%%%    
    upper_diagonal_1(1:3:3*N-1) = HzNt;
    upper_diagonal_1(2:3:3*N-1) = HxNt;
    %%%
    lower_diagonal_1(1:3:3*N-1) =-HzNt;
    lower_diagonal_1(2:3:3*N-1) =-HxNt;
    %%%
    upper_diagonal_2(1:3:3*N-2) =-HyNt;
    %%%
    lower_diagonal_2(1:3:3*N-2) = HyNt;
    %%%
    %%%
    %%%          
    Wt =  diag(upper_diagonal_1,1) + diag(upper_diagonal_2,2) + ...
          diag(lower_diagonal_1,-1) + diag(lower_diagonal_2,-2);
    Wt=Wt*-g/(1+alpha^2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %%% This is W of the damping term %%%%%%%%%%%%%%%%%%%%   
    Ht(1:3:3*N-2,1)= HxNt;
    Ht(2:3:3*N-1,1)= HyNt;
    Ht(3:3:3*N,1)  = HzNt;
    %%%%%%%%%%%
    % W_M(M.H) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%
    %    [M.H  0   0    %%%%
    %      0  M.H  0    %%%%
    %      0   0  M.H]  %%%%
    MMdHx = HxNt .* mx;
    Wt_MMdHx(1:3:3*N-2)= MMdHx;
    Wt_MMdHx(2:3:3*N-1)= MMdHx;
    Wt_MMdHx(3:3:3*N)  = MMdHx;
    %%%
    MMdHy = HyNt .* my;
    Wt_MMdHy(1:3:3*N-2)= MMdHy;
    Wt_MMdHy(2:3:3*N-1)= MMdHy;
    Wt_MMdHy(3:3:3*N)  = MMdHy;
    %%%
    MMdHz = HzNt .* mz;
    Wt_MMdHz(1:3:3*N-2)= MMdHz;
    Wt_MMdHz(2:3:3*N-1)= MMdHz;
    Wt_MMdHz(3:3:3*N)  = MMdHz;
    %%%%%%%%%
    Wt_MMdH = diag(Wt_MMdHx + Wt_MMdHy + Wt_MMdHz);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    %%%%%%%%%%
    % W_H(M.M) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    %%%%%%%%%%
    %%%
    above_diagonal_1 = zeros(3*N-1,1);
    above_diagonal_2 = zeros(3*N-2,1);
    below_diagonal_1 = zeros(3*N-1,1);
    below_diagonal_2 = zeros(3*N-2,1);  
    %    [HxMx  HxMy   HxMz    %%%%
    %     HyMx  HyMy   HyMz    %%%%
    %     HzMx  HzMy   HzMz    %%%%
    %%%
    Hx_my = HxNt.* my;
    Hy_mz = HyNt.* mz;
    above_diagonal_1(1:3:3*N-1)= Hx_my;
    above_diagonal_1(2:3:3*N-1)= Hy_mz;
    %%%
    Hy_mx = HyNt.* mx;
    Hz_my = HzNt.* my;
    below_diagonal_1(1:3:3*N-1)= Hy_mx;
    below_diagonal_1(2:3:3*N-1)= Hz_my;
    %%%
    Hx_mz = HxNt.* mz;
    above_diagonal_2(1:3:3*N-2)= Hx_mz;
    %%%
    Hz_mx = HzNt.* mx;
    below_diagonal_2(1:3:3*N-2)= Hz_mx;
    %%%
    %%%
    %%%
    Wt_HMdM = diag(Ht.*m) +...
              diag(above_diagonal_1,1)+...
              diag(above_diagonal_2,2)+...
              diag(below_diagonal_1,-1)+...
              diag(below_diagonal_2,-2);
    % M*M*H = M(M.H)- H(M.M)  
    Wt_McMcH = (Wt_MMdH - Wt_HMdM)* g*alpha/((1+alpha^2)*Ms);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Euler's method to get a preliminary %%%%%
          % approximation for the next value of %%%%%
          % magnetization          %%%%%%%%%%%%%%%%%%
          M_underestimate = m + dt*(Wt + Wt_McMcH)*m;              
          Mux = M_underestimate(1:3:3*N-2);
          Muy = M_underestimate(2:3:3*N-1);
          Muz = M_underestimate(3:3:3*N);
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for nn = 1 : N
        Mux_matrix(:,nn) = Mux; 
        Muy_matrix(:,nn) = Muy; 
        Muz_matrix(:,nn) = Muz; 
    end
    Mux_matrix = Mux_matrix(A);
    Muy_matrix = Muy_matrix(A);
    Muz_matrix = Muz_matrix(A);    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%                Dipolar Field               %%%%%%%%
    %%%           H = (3(m.r)*r - m*r^2)/r^5           %%%%%%
    mdotr = Mux_matrix.*rx + Muy_matrix.*ry + Muz_matrix.*rz;
    Hx_Dipole_mdotr = mdotr .*rx ./ (r.^5);
    Hy_Dipole_mdotr = mdotr .*ry ./ (r.^5);
    Hz_Dipole_mdotr = mdotr .*rz ./ (r.^5);
    %%%
    Hx_Dipole_m = Mux_matrix ./ (r.^3);
    Hy_Dipole_m = Muy_matrix ./ (r.^3);    
    Hz_Dipole_m = Muz_matrix ./ (r.^3);
    %%%
    Hx_Dipole_mdotr(1,:) = zeros(1,N);
    Hx_Dipole_m(1,:) = zeros(1,N);
    Hx_Dip = 3*Hx_Dipole_mdotr - Hx_Dipole_m;
%  Hx_Dip(2,2)=0;
%Hx_Dip(3,3)=0;
%Hx_Dip(3,2)=0;%
%Hx_Dip(2,3)=0;%Hx_Dip(4,2)=0;
    Hx_Dip_Collective = sum(Hx_Dip,1)';
    %
    Hy_Dipole_mdotr(1,:) = zeros(1,N);
    Hy_Dipole_m(1,:) = zeros(1,N);    
    Hy_Dip = 3*Hy_Dipole_mdotr - Hy_Dipole_m;
%  Hy_Dip(2,2)=0;
%Hy_Dip(3,3)=0;
%Hy_Dip(3,2)=0;%
%Hy_Dip(2,3)=0;%Hy_Dip(4,2)=0;
    Hy_Dip_Collective = sum(Hy_Dip,1)';
    %
    Hz_Dipole_mdotr(1,:) = zeros(1,N);
    Hz_Dipole_m(1,:) = zeros(1,N);    
    Hz_Dip = 3*Hz_Dipole_mdotr - Hz_Dipole_m;
%  Hz_Dip(2,2)=0;
%Hz_Dip(3,3)=0;
%Hz_Dip(3,2)=0;%
%Hz_Dip(2,3)=0;%Hz_Dip(4,2)=0;
    Hz_Dip_Collective = sum(Hz_Dip,1)';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %%%%%%%
    %%%%%%%
    %%%%%%%
    %%%%% Uniaxial Anisotropy %%%%%%%%%%%%%%%
    Huax = zeros(N,1);    
    Huay = 2 * 1 * Ku2*Muy/Ms^2;
    Huaz = - 2 * 1 * Kp*Muz/Ms^2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Huay(1) = 0;  
    %Huax(1) = 2 * 1 * Ku2*Mux(1)/Ms^2;    
    %
    %%%%%   Overall field including anisotropy and the     %%%%%%%%
    %%%%%        field due to all other dipoles            %%%%%%%%
    HxNtdt = Hx(1:N,tt+1) + 1*Huax + 1*Hx_Dip_Collective;
    HyNtdt = Hy(1:N,tt+1) + 1*Huay + 1*Hy_Dip_Collective;
    HzNtdt = Hz(1:N,tt+1) + 1*Huaz + 1*Hz_Dip_Collective;       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %%% This is W of the precession term %%%%%%%%%%%%%%%%%%%%%%%%%%
    upper_diagonal_1 = zeros(3*N-1,1);
    upper_diagonal_2 = zeros(3*N-2,1);
    lower_diagonal_1 = zeros(3*N-1,1);
    lower_diagonal_2 = zeros(3*N-2,1);
    %%%
    %    [0     Hz  -Hy  %%%%
    %    -Hz    0    Hx  %%%%
    %     Hy   -Hx   0]  %%%%    
    upper_diagonal_1(1:3:3*N-1) = HzNtdt;
    upper_diagonal_1(2:3:3*N-1) = HxNtdt;
    %%%
    lower_diagonal_1(1:3:3*N-1) =-HzNtdt;
    lower_diagonal_1(2:3:3*N-1) =-HxNtdt;
    %%%
    upper_diagonal_2(1:3:3*N-2) =-HyNtdt;
    %%%
    lower_diagonal_2(1:3:3*N-2) = HyNtdt;
    %%%
    %%%
    %%%          
    Wtdt =  diag(upper_diagonal_1,1) + diag(upper_diagonal_2,2) + ...
            diag(lower_diagonal_1,-1) + diag(lower_diagonal_2,-2);
    Wtdt=Wtdt*-g/(1+alpha^2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %%% This is W of the damping term %%%%%%%%%%%%%%%%%%%%   
    Htdt(1:3:3*N-2,1)= HxNtdt;
    Htdt(2:3:3*N-1,1)= HyNtdt;
    Htdt(3:3:3*N,1)  = HzNtdt;
    %%%%%%%%%%%
    % W_M(M.H) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%
    %    [M.H  0   0    %%%%
    %      0  M.H  0    %%%%
    %      0   0  M.H]  %%%%
    MMdHx = HxNtdt .* Mux;
    Wtdt_MMdHx(1:3:3*N-2)= MMdHx;
    Wtdt_MMdHx(2:3:3*N-1)= MMdHx;
    Wtdt_MMdHx(3:3:3*N)  = MMdHx;
    %%%
    MMdHy = HyNtdt .* Muy;
    Wtdt_MMdHy(1:3:3*N-2)= MMdHy;
    Wtdt_MMdHy(2:3:3*N-1)= MMdHy;
    Wtdt_MMdHy(3:3:3*N)  = MMdHy;
    %%%
    MMdHz = HzNtdt .* Muz;
    Wtdt_MMdHz(1:3:3*N-2)= MMdHz;
    Wtdt_MMdHz(2:3:3*N-1)= MMdHz;
    Wtdt_MMdHz(3:3:3*N)  = MMdHz;
    %%%%%%%%%
    Wtdt_MMdH = diag(Wtdt_MMdHx + Wtdt_MMdHy + Wtdt_MMdHz);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    %%%%%%%%%%
    % W_H(M.M) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    %%%%%%%%%%
    %%%
    above_diagonal_1 = zeros(3*N-1,1);
    above_diagonal_2 = zeros(3*N-2,1);
    below_diagonal_1 = zeros(3*N-1,1);
    below_diagonal_2 = zeros(3*N-2,1);  
    %    [HxMx  HxMy   HxMz    %%%%
    %     HyMx  HyMy   HyMz    %%%%
    %     HzMx  HzMy   HzMz    %%%%
    %%%
    Hx_my = HxNtdt.* Muy;
    Hy_mz = HyNtdt.* Muz;
    above_diagonal_1(1:3:3*N-1)= Hx_my;
    above_diagonal_1(2:3:3*N-1)= Hy_mz;
    %%%
    Hy_mx = HyNtdt.* Mux;
    Hz_my = HzNtdt.* Muy;
    below_diagonal_1(1:3:3*N-1)= Hy_mx;
    below_diagonal_1(2:3:3*N-1)= Hz_my;
    %%%
    Hx_mz = HxNtdt.* Muz;
    above_diagonal_2(1:3:3*N-2)= Hx_mz;
    %%%
    Hz_mx = HzNtdt.* Mux;
    below_diagonal_2(1:3:3*N-2)= Hz_mx;
    %%%
    %%%
    %%%
    Wtdt_HMdM = diag(Htdt.*M_underestimate) +...
                diag(above_diagonal_1,1) +...
                diag(above_diagonal_2,2) +...
                diag(below_diagonal_1,-1)+...
                diag(below_diagonal_2,-2);
    % M*M*H = M(M.H)- H(M.M)  
    Wtdt_McMcH = (Wtdt_MMdH - Wtdt_HMdM)* g*alpha/((1+alpha^2)*Ms);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    %%%%%
    %%%%%
    %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Heun's method
    slope_underestimate = (Wt + Wt_McMcH)*m ;     
    slope_overstimate = (Wtdt + Wtdt_McMcH )*M_underestimate;
    M = m + dt/2*(slope_underestimate + slope_overstimate);  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
    m=M;
    mx = m(1:3:3*N-2);
    my = m(2:3:3*N-1);
    mz = m(3:3:3*N);
    %%%%%%%%%%%%%%
    Mx(1:N,tt)=mx;
    My(1:N,tt)=my;
    Mz(1:N,tt)=mz;
    %%%%%%%%%%%%%%
    Hanix(1:N,tt) = Huax;
    Haniy(1:N,tt) = Huay;
    Haniz(1:N,tt) = Huaz;
    %%%%%%%%%%%%%%
    Hdipx(1:N,tt) = Hx_Dip_Collective;
    Hdipy(1:N,tt) = Hy_Dip_Collective;    
    Hdipz(1:N,tt) = Hz_Dip_Collective;    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
Heffx = Hx + Hanix + Hdipx;
Heffy = Hy + Haniy + Hdipy;
Heffz = Hz + Haniz + Hdipz;
%%%
dE_dtCGS_HeffM = zeros(N,Nt-2);
dE_dtCGS_Mdip = zeros(N,Nt-2);
dE_dtCGS_MH = zeros(N,Nt-2);
for nn=1:N
    dE_dtCGS_HeffM(nn,:) =   Heffx(nn,2:Nt-1) .* diff(Mx(nn,:)) +...
                                Heffy(nn,2:Nt-1) .* diff(My(nn,:)) +...
                                   Heffz(nn,2:Nt-1) .* diff(Mz(nn,:));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %dE_dtCGS_MH(nn,:) =      Mx(nn,1:Nt-2) .* diff(Hx(nn,1:Nt-1)) +...
    %                            My(nn,1:Nt-2) .* diff(Hy(nn,1:Nt-1)) +...
    %                               Mz(nn,1:Nt-2) .* diff(Hz(nn,1:Nt-1));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %dE_dtCGS_Mdip(nn,:) =    Mx(nn,1:Nt-2) .* diff(Hdipx(nn,1:Nt-1)) +...
    %                           My(nn,1:Nt-2) .* diff(Hdipy(nn,1:Nt-1))+...
    %                              Mz(nn,1:Nt-2) .* diff(Hdipz(nn,1:Nt-1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Power in CGS units %%%%%%%%%%%%%%%
dE_dtCGS_HeffM = - dE_dtCGS_HeffM /dt;
%dE_dtCGS_MH    = - dE_dtCGS_MH /dt;
%dE_dtCGS_Mdip  = - dE_dtCGS_Mdip /dt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%% Convert to Watts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dE_dtMKS_HeffM = dE_dtCGS_HeffM * 1e-7;
%dE_dtMKS_MH    = dE_dtCGS_MH * 1e-7;
%dE_dtMKS_Mdip  = dE_dtCGS_Mdip * 1e-7;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%% Dissipated energy in Joules  %%%%%%%%
E_diss_HeffM = sum(dE_dtMKS_HeffM,2) * dt;
%E_MH    = sum(dE_dtMKS_MH,2) * dt;
%E_Mdip  = sum(dE_dtMKS_Mdip,2) * dt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%% Dissipated/Exchange energy in eV %%%%%%%%%%%%%
E_diss_HeffM_eV = E_diss_HeffM / 1.6e-19;
%E_MH_eV    = E_MH / 1.6e-19;
%E_Mdip_eV  = E_Mdip / 1.6e-19;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%% Dissipated energy in units of kT %%%%
%E_MH_kT    = E_MH_eV / .0259;
%E_Mdip_kT  = E_Mdip_eV / .0259;
%%%%
E_diss_HeffM_kT = E_diss_HeffM_eV / .0259;
%sum(E_diss_HeffM_kT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E_disskT(rr) = sum(E_diss_HeffM_kT)
%E_diss_HeffM_kT = zeros(3,1);
end
%%%%%%%%%%%%%%
% Figures  %%%
%%%%%%%%%%%%%%
figure(1);hold on;orient('landscape')
subplot(4,1,1);
plot(t(1:Nt-1)/1e-9,My(1,1:Nt-1)/Ms,'r','linewidth',2);hold on;grid on
title('Magnet 1 ','Fontsize',18)
ylabel('M_y (y:easy axis)','Fontsize',12);xlabel('Time (ns)--->','Fontsize',12);
axis('tight')
%
subplot(4,1,2);plot(t(1:Nt-1)/1e-9,My(2,1:Nt-1)/Ms,'m','linewidth',2);grid on;
title('Magnet 2 ','Fontsize',18)
ylabel('M_y (y:easy axis)','Fontsize',12);xlabel('Time (ns)--->','Fontsize',12);
axis tight
%
subplot(4,1,3);plot(t(1:Nt-1)/1e-9,My(3,1:Nt-1)/Ms,'m','linewidth',2);grid on;
title('Magnet 3 ','Fontsize',18)
ylabel('M_y (y:easy axis)','Fontsize',12);xlabel('Time (ns)--->','Fontsize',12);
axis tight
%
subplot(4,1,4);plot(t(1:Nt-1)/1e-9,My(4,1:Nt-1)/Ms,'m','linewidth',2);grid on;
title('Magnet 3 ','Fontsize',18)
ylabel('M_y (y:easy axis)','Fontsize',12);xlabel('Time (ns)--->','Fontsize',12);
axis tight
%%% figure(11)
%%% plotyy(t(1:Nt-2)/1e-9,Hx(2,1:Nt-2)/750,t(1:Nt-2)/1e-9,My(2,1:Nt-2)/Ms)
%%%%%%%%%
figure(2)
hold on;orient('landscape')
subplot(4,1,1)
plot(t(1:Nt-2)/1e-9,-dE_dtMKS_HeffM(1,:),'r','linewidth',2);grid on;
title('Dissipated power in magnet 1 ','Fontsize',18)
ylabel('Power [w]','Fontsize',18);xlabel('Time (ns)--->','Fontsize',16);
axis tight
%
subplot(4,1,2)
plot(t(1:Nt-2)/1e-9,-dE_dtMKS_HeffM(2,:),'r','linewidth',2);grid on;
ylabel('Power [w]','Fontsize',18);xlabel('Time (ns)--->','Fontsize',16);
title('Dissipated power in magnet 2 ','Fontsize',18)
axis tight
%
subplot(4,1,3)
plot(t(1:Nt-2)/1e-9,-dE_dtMKS_HeffM(3,:),'r','linewidth',2);grid on;
ylabel('Power [w]','Fontsize',18);xlabel('Time (ns)--->','Fontsize',16);
title('Dissipated power in magnet 3 ','Fontsize',18)
axis tight
%
subplot(4,1,4)
plot(t(1:Nt-2)/1e-9,-dE_dtMKS_HeffM(4,:),'r','linewidth',2);grid on;
ylabel('Power [w]','Fontsize',18);xlabel('Time (ns)--->','Fontsize',16);
title('Dissipated power in magnet 3 ','Fontsize',18)
axis tight
%%%%%%%%%
figure(3)
hold on; orient('landscape')
subplot(4,1,1)
plot(t(1:Nt)/1e-9,Hx(1,:),'linewidth',3);hold on
title('Applied Magnetic Field on magnet 1 [Oe]','Fontsize',18);xlabel('Time (ns)--->','Fontsize',16);
axis('tight')
%
subplot(4,1,2)
plot(t(1:Nt)/1e-9,Hx(2,:),'linewidth',3);hold on
title('Applied Magnetic Field on magnet 2 [Oe]','Fontsize',18);xlabel('Time (ns)--->','Fontsize',16);
%
subplot(4,1,3)
plot(t(1:Nt)/1e-9,Hx(3,:),'linewidth',3);hold on
title('Applied Magnetic Field on magnet 3 [Oe]','Fontsize',18);xlabel('Time (ns)--->','Fontsize',16);
%
subplot(4,1,4)
plot(t(1:Nt)/1e-9,Hx(4,:),'linewidth',3);hold on
title('Applied Magnetic Field on magnet 3[Oe]','Fontsize',18);xlabel('Time (ns)--->','Fontsize',16);
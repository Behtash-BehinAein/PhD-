% Behtash Behin-Aein
% This is the code for the dimensionless LLG equation
% for a thin film single domain Ferro magnet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc;%warning off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Z                      % Z: easy axis
%                   |                      % X: in-plane hard axis
%                   |                      % Y: out-of-plane hard axis
%                   |
%                   |--------------- Y
%                 *
%               * 
%             *
%           *      
%        X
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants  -----------------------------------------------------------   
g = 1.76*1e7;                        % Gyromagnetic ratio [(rad)/(Oe.s)]   
hbar  = 1.055e-34;                   % Reduced Planck's constant [J.s]
q = 1.6e-19;                         % Charge of an electron [C]  
mu_B = g*hbar/2;                     % Bohr Magneton
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%    All equations are written in a dimensionless form          %%%%
    %%% These parameters pertain to the ordinary LLG equation       %%%%%
g = 1.76*1e7;                % Gyromagnetic ratio [(rad)/(Oe.s)]    
alpha = 0.01;                % Gilbert damping parameter
Ms = 780;                    % Saturation Magnetization [emu/cm^3]
Ku2 = 150150;%1.0*4.144e4;% Uniaxial anisotropy constant [erg/cm^3]
V = 1*(100*42*2)*1e-21;     % Volume [cm^3]
MsV = Ms*V;                   % [emu]
Ku2V = Ku2*V;                 % [erg]
KuV_kT = Ku2V*1e-7/1.6e-19/0.0259      %[kT]
Hc = 2*Ku2/Ms;                 % Switching field [Oe]
%MsV*Hc*1e-7/1.6e-19/0.0259;
%(MsV*Hc)*1e-7*2*1.6e-19/1.055e-3 4;
Hd = 1*4*pi*Ms/Hc;             % Demagnetizing field
Kd = Ms*Hd*Hc / 2;
KdV = Kd *V;
Ns = MsV/9.274e-21                                       % Number of spins
tau_c = (1+alpha^2)/(2*alpha*abs(g)*Hc);          % Scaling factor for time
% Critical Switching current
Isc = (2*q/hbar) * alpha * (2*Ku2V + KdV) *1e-7 *2;          % [A]
Icc = 10*Isc;
A1= - 5.5e-3;%-1.2*Icc%-2.7*Icc                                       % Input Bias charge current
I_H_conv=hbar/2/q/(MsV*Hc*1e-7);
% -------------------------------------------------------------------------
% Time grid   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_resolution = 40000;
time_duration = 0.2; 
t = linspace(0, time_duration, time_duration*time_resolution); 
dt = t(3) - t(2); Nt = length(t);
tspan_ns = time_duration * tau_c/1e-9;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pulse shape                                  %%%%%%%%%%%%%%%
bp = 00; bp_ns = bp*dt*tau_c/1e-9;
ap = 0*round(Nt); ap_ns = ap*dt*tau_c/1e-9;
Nramp = 2;
pw = round(Nt-ap-bp-1*Nramp); pw_ns = pw*dt*tau_c /1e-9;
pulse = [zeros(1,bp) linspace(0,1,Nramp)  ones(1,pw) linspace(1,0,Nramp) zeros(1,ap)];
% Directions   %%%%%%%%%%%%%%%
theta_p = 0;               % Pulse
theta_s = pi;      % Source Spin torque
phi_p = 0;                 % Pulse
phi_s = 0;                 % Source Spin torque
% -----------------------------------------------------------------------
Mhat = [sin(theta_s)*cos(theta_s) sin(theta_s)*sin(theta_s) cos(theta_s)]';
% -----------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%       External Fields                     %%%%%%%%%%%%%%%%%
% Pulse Magnitude   %%%%%%%%%% 
mhp = 0.0;
% The externally applied magnetic pulse  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hpx = mhp*pulse*sin(theta_p)*cos(phi_p);
hpy = mhp*pulse*sin(theta_p)*sin(phi_p);
hpz = mhp*pulse*cos(theta_p);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%    Pre-allocating memory   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Magnetization components
mx = zeros(1,Nt);
my = zeros(1,Nt);
mz = zeros(1,Nt);
% Uniaxial anisotropy field (Internal)
hux = zeros(1,Nt);
huy = zeros(1,Nt);
huz = zeros(1,Nt);
% Out of plane demagnetizing field (Internal)
hdx = zeros(1,Nt);
hdy = zeros(1,Nt);
hdz = zeros(1,Nt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the exerted torque on the magnetic moment
dmdt_Ideal = zeros(3,Nt);
% This is the precessional interaction in the z direction
Pz = zeros(1,Nt);
% This is the damping interaction in the z direction
Dz = zeros(1,Nt);
%%%%%%% Dissipated Power   %%%%%%%%%%%%%%%%%%
Pd = zeros(1,Nt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Iz = zeros(1,Nt);
Ix = zeros(1,Nt);
Iy = zeros(1,Nt);
ISL = zeros(1,Nt);
%ISL_prime = zeros(1,Nt);
ISLz = zeros(1,Nt);
ISLx = zeros(1,Nt);
ISLy = zeros(1,Nt);
IFL = zeros(1,Nt);
hsl = zeros(1,Nt);
hfl = zeros(1,Nt);
P3 = zeros(1,Nt);
ps = zeros(1,Nt);
ds = zeros(1,Nt);
%----------------------------------
PMx = zeros(1,Nt);
PMy = zeros(1,Nt);
PMz = zeros(1,Nt);
%---------------------------------
DMx = zeros(1,Nt);
DMy = zeros(1,Nt);
DMz = zeros(1,Nt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%                     Initializing                        %%%%%%%%%%%%
mz(1) =  0.999; mx(1) = sqrt(1-mz(1)^2); my(1) = 0; m = [mx(1);my(1);mz(1)];
%         -----------------------------------------------------
[h_clk_sl h_clk_fl Isl Islz Islx Isly Ifl P1 Iz2 Ix2 Iy2] =otani_ckt2(m,Mhat,A1,MsV,Ku2V);
Iz(1)=Iz2; % calculating the IsM=Iz given a M
Ix(1)=Ix2;
Iy(1)=Iy2;
%ISL_prime(1) = h_clk_sl / I_H_conv * norm(cross(m,cross(m,Mhat)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for tt = 2 : Nt
%%% Spin-Circuit LLG interface --------------------------------------------
[h_clk_sl h_clk_fl Isl Islz Islx Isly Ifl P1 Iz2 Ix2 Iy2] =otani_ckt2(m,Mhat,A1,MsV,Ku2V);
hsl(tt) = h_clk_sl * pulse(tt);
%ISL_prime(tt) = hsl(tt) / I_H_conv * norm(cross(m,cross(m,Mhat)));
hfl(tt) = h_clk_fl * pulse(tt);
ps(tt)  =  hfl(tt) - alpha*hsl(tt); 
ds(tt)  =  hfl(tt) + hsl(tt)/alpha;
%----------------------------------
PMx(tt) = ps(tt) * sin(theta_s)*cos(phi_s);
PMy(tt) = ps(tt) * sin(theta_s)*sin(phi_s);
PMz(tt) = ps(tt) * cos(theta_s);
%---------------------------------
DMx(tt) = ds(tt) * sin(theta_s)*cos(phi_s);
DMy(tt) = ds(tt) * sin(theta_s)*sin(phi_s);
DMz(tt) = ds(tt) * cos(theta_s); 
%---------------------------------
Iz(tt) = Iz2 * pulse(tt); % calculating the IsM=Iz given a M
Ix(tt) = Ix2 * pulse(tt);
Iy(tt) = Iy2 * pulse(tt);
ISL(tt) = Isl * pulse(tt);
ISLz(tt) = Islz * pulse(tt);
ISLx(tt) = Islx * pulse(tt);
ISLy(tt) = Isly * pulse(tt);
IFL(tt) = Ifl * pulse(tt); 
P3(tt)= P1*(t(tt)-t(tt-1)); 
%---------------------------------------------         
% -------------------------------------------------------------------------    
%%%%   Interaction of the magnetic moment with various fields     %%%%%%
                             M = [m';m';m']; 
%%%%   Uniaxial field %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hux(tt-1) = 0;              % This is the in-plane  x-axis anisotropy field
huy(tt-1) = 0;              % This is the OutOfPlane y-axi anisotropy field
huz(tt-1) = mz(tt-1);       % This is the in-plane  z-axis anisotropy field
hu = [hux(tt-1); huy(tt-1); huz(tt-1)];                    HU = [hu hu hu];
%%%%   Demagnetzing field  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hdx(tt-1) = 0;             % This is the in-plane x-axis demagnetzing field
hdy(tt-1) = - Hd*my(tt-1); % This is the OutOfPlane y demagnetzing field      
hdz(tt-1) = 0;             % This is the in-plane z-axis demagnetzing field
hd = [hdx(tt-1); hdy(tt-1); hdz(tt-1)];                    HD = [hd hd hd];
%%%%   External pulse      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hp = [hpx(tt-1); hpy(tt-1); hpz(tt-1)];                    HP = [hp hp hp];
%%%%   Source Spin torque                                %%%%%%%%%%%%%%%%%%
hdM = [DMx(tt-1); DMy(tt-1); DMz(tt-1)];
                                                       HDM = [hdM hdM hdM]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %%%%%%%     Precessional interaction          %%%%%%%%
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Interaction with effective uniaxial field
hum = [ 0            huz(tt-1)  -huy(tt-1)
       -huz(tt-1)    0           hux(tt-1)
        huy(tt-1)   -hux(tt-1)          0];
                                         Uniaxial = hum*m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Interaction with effective demagnetizing  field
hdm = [ 0            hdz(tt-1)   -hdy(tt-1)
       -hdz(tt-1)    0            hdx(tt-1)
        hdy(tt-1)   -hdx(tt-1)           0];
                                         Demagnetizing = hdm*m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Interaction with the externally applied pulse
hpm = [ 0            hpz(tt-1)   -hpy(tt-1)
       -hpz(tt-1)    0            hpx(tt-1)
        hpy(tt-1)   -hpx(tt-1)           0];
                                         Pulse = hpm*m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Interaction with the effective field from source spin torque
hsstm = [0               PMz(tt-1)   -PMy(tt-1)
        -PMz(tt-1)       0            PMx(tt-1)
         PMy(tt-1)    -PMx(tt-1)             0];
                                         SpinTorqueSource = hsstm*m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %     Sum of all terms        %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PI = Uniaxial + Demagnetizing + Pulse + SpinTorqueSource;
clear  Uniaxial  Demagnetizing  Pulse  SpinTorqueSource  
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %%%%%%%        Damping interaction            %%%%%%%%
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
                   %      First term           %
                   %    [m.h  0   0    [mx  %%%%
                   %      0  m.h  0     my  %%%%
                   %      0   0  m.h]   mz] %%%% 
%%%%   Interaction with effective uniaxial field
hu_mDOTh = diag([dot(m,hu) dot(m,hu) dot(m,hu)]);
                                           Uniaxial = hu_mDOTh*m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Interaction with effective demagnetizing  field
hd_mDOTh = diag([dot(m,hd) dot(m,hd) dot(m,hd)]);
                                           Demagnetizing = hd_mDOTh*m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Interaction with the externally applied pulse
hp_mDOTh = diag([dot(m,hp) dot(m,hp) dot(m,hp)]);
                                           Pulse = hp_mDOTh*m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Interaction with the source spin torque
hsst_mDOTh = diag([dot(m,hdM) dot(m,hdM) dot(m,hdM)]);
                                           SpinTorqueSource = hsst_mDOTh*m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DI1 =Uniaxial + Demagnetizing + Pulse + SpinTorqueSource;
clear  Uniaxial  Demagnetizing  Pulse  SpinTorqueSource
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               %          Second term                %
               %    [hx.mx  hx.my   hx.mz   [mx   %%%%
               %     hy.mx  hy.my   hy.mz    my   %%%%
               %     hz.mx  hz.my   hz.mz]   mz]  %%%%
%%%%   Interaction with effective uniaxial field
hu_HM = HU.*M;
                               Uniaxial = hu_HM*m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Interaction with effective demagnetizing  field
hd_HM = HD.*M;
                               Demagnetizing = hd_HM*m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Interaction with the externally applied pulse
hp_HM = HP.*M;
                               Pulse = hp_HM*m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Interaction with the effective field from linear spin torque
hsst_HM = HDM.*M;
                               SpinTorqueSource = hsst_HM*m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DI2 =Uniaxial + Demagnetizing + Pulse + SpinTorqueSource;
clear  Uniaxial  Demagnetizing  Pulse  SpinTorqueSource
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %       Damping Term         %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            DI = DI1 - DI2;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%
%%%%%%%
%%%%%%%  First estimation of m_(t+1) based on Euler's method    %%%%%%%%%%%
            dmdt_Lower = - (1/(2*alpha)) * 1 * PI - (1/2)* DI;
            m_Intermediate = m + dmdt_Lower * dt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Interaction of the magnetic moment with various fields        %%%%%%
% These are the intermediate quantities to be used by the Heun's method  %%      
%%%%   Magnetic moment                                            %%%%%%
Imx = m_Intermediate(1);
Imy = m_Intermediate(2);
Imz = m_Intermediate(3);
                                                        Im = [Imx;Imy;Imz];
                                                        IM = [Im';Im';Im'];
%%%%   Uniaxial field %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ihux = 0;
Ihuy = 0;
Ihuz = Imz;
Ihu = [Ihux; Ihuy; Ihuz];                              IHU = [Ihu Ihu Ihu];
%%%%   Demagnetzing field  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
Ihdx = 0;
Ihdy = - Hd*Imy;
Ihdz = 0;
Ihd = [Ihdx; Ihdy; Ihdz];                              IHD = [Ihd Ihd Ihd];
%%%%   External pulse      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ihp = [hpx(tt); hpy(tt); hpz(tt)];                     IHP = [Ihp Ihp Ihp];
%%%%   Source Spin torque                                %%%%%%%%%%%%%%%%%%
IhdM = [DMx(tt); DMy(tt); DMz(tt)];                IHDM = [IhdM IhdM IhdM]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %%%%%%%     Precessional interaction          %%%%%%%%
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Interaction with effective uniaxial field
Ihum = [ 0      Ihuz  -Ihuy
       -Ihuz     0     Ihux
        Ihuy   -Ihux      0];
                                                        Uniaxial = Ihum*Im;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Interaction with effective demagnetizing  field
Ihdm = [ 0      Ihdz   -Ihdy
       -Ihdz    0       Ihdx
        Ihdy   -Ihdx       0];
                                                   Demagnetizing = Ihdm*Im;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Interaction with the externally applied pulse
Ihpm = [ 0         hpz(tt)   -hpy(tt)
       -hpz(tt)    0          hpx(tt)
        hpy(tt)   -hpx(tt)         0];
                                                           Pulse = Ihpm*Im;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Source spin torque
Ihsstm = [0          PMz(tt)      -PMy(tt)
        -PMz(tt)       0           PMx(tt)
         PMy(tt)    -PMx(tt)             0];
                                              SpinTorqueSource = Ihsstm*Im;                                            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %     Sum of all terms        %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IPI =Uniaxial + Demagnetizing + Pulse + SpinTorqueSource;
clear  Uniaxial  Demagnetizing  Pulse  SpinTorqueSource
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %%%%%%%        Damping interaction            %%%%%%%%
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
                   %      First term           %
%%%%   Interaction with effective uniaxial field
Ihu_mDOTh = diag([dot(Im,Ihu) dot(Im,Ihu) dot(Im,Ihu)]);
                                                   Uniaxial = Ihu_mDOTh*Im;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Interaction with effective demagnetizing  field
Ihd_mDOTh = diag([dot(Im,Ihd) dot(Im,Ihd) dot(Im,Ihd)]);
                                              Demagnetizing = Ihd_mDOTh*Im;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Interaction with the externally applied pulse
Ihp_mDOTh = diag([dot(Im,Ihp) dot(Im,Ihp) dot(Im,Ihp)]);
                                                      Pulse = Ihp_mDOTh*Im;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Source spin torque
Ihsst_mDOTh = diag([dot(Im,IhdM) dot(Im,IhdM) dot(Im,IhdM)]);
                                         SpinTorqueSource = Ihsst_mDOTh*Im;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %     Sum of all terms        %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IDI1 =Uniaxial + Demagnetizing + Pulse + SpinTorqueSource;
clear  Uniaxial  Demagnetizing  Pulse  SpinTorqueSource
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                         
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               %          Second term                %
%%%%   Interaction with effective uniaxial field
Ihu_HM = IHU.*IM;
                                                      Uniaxial = Ihu_HM*Im;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Interaction with effective demagnetizing  field
Ihd_HM = IHD.*IM;
                                                 Demagnetizing = Ihd_HM*Im;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%
%%%%   Interaction with the externally applied pulse
Ihp_HM = IHP.*IM;
                                                         Pulse = Ihp_HM*Im;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Source spin torque
Ihsst_HM = IHDM.*IM;
                                            SpinTorqueSource = Ihsst_HM*Im;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %     Sum of all terms        %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IDI2 =Uniaxial + Demagnetizing + Pulse + SpinTorqueSource;
clear  Uniaxial  Demagnetizing  Pulse  SpinTorqueSource
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %       Damping Term          %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        IDI = IDI1 - IDI2;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                                 
%%%%%%%  Calculation of higher slope based on initial estimation   %%%%%%%%
dmdt_Higher = -(1/(2*alpha)) * 1 * IPI - (1/2) * IDI;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Heun's second order method (averaging of Lower and Higher slopes
%%%%%%            to get closer to the ideal slope)               %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = -(0.5*PI/alpha + 0.5*IPI/alpha)/2; 
D = -(0.5*DI + 0.5*IDI)/2; 
Pz(tt) = P(3);
Dz(tt) = D(3);
dmdt_Ideal(:,tt-1) = (dmdt_Higher + dmdt_Lower) /2;
Pd(tt)= (4*alpha^2/(1+alpha^2))*dot(dmdt_Ideal(:,tt-1),dmdt_Ideal(:,tt-1));
m = m + dmdt_Ideal(:,tt-1) * dt;
mx(tt) = m(1); my(tt) = m(2); mz(tt) = m(3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
sumDz = sum(Dz*dt);%/(mz(1)-mz(end))
sumPz = sum(Pz*dt);%/(mz(1)-mz(end))
%sumPz + sumDz;
% ----------------------------------------------
hdemag_y = -Hd*my;
Pz_hdemag_y = mx.*hdemag_y;
sumPz_hdemag_y = sum(Pz_hdemag_y*dt) *-0.5/alpha;
% ----------------------------------------------
% m_cross_dmdt = zeros(3,Nt);
% for tt=1:Nt
%     m_cross_dmdt(:,tt) = cross(m,dmdt_Ideal(:,tt));
% end
% sum(alpha*m_cross_dmdt(3,:)*dt)
% ----------------------------------------------
tau = t;
dtau = tau(2) - tau(1);
t = t*tau_c*1e9;
t_sec = t / 1e9;
Dissipated_Energy = trapz(tau,Pd);
Is = sqrt(Ix.^2+Iy.^2+Iz.^2);
Ns_provided = sum(Is*dt*tau_c)/q;
Nsl_provided = sum(ISL*dt*tau_c)/q;
Nslz_provided = sum(ISLz*dt*tau_c)/q;
Nslx_provided = sum(ISLx*dt*tau_c)/q;
Nsly_provided = sum(ISLy*dt*tau_c)/q;
%------------------------------------
ratio_s = Ns_provided / Ns;
ratio_sl =  Nsl_provided / Ns;
ratio_slz = Nslz_provided / Ns%(Ns*(mz(1)-mz(end)))
sl_Over_slz = ratio_sl/ratio_slz;
ratio_slx = Nslx_provided / Ns;%(Ns*(mx(1)-mx(end)))
ratio_sly = Nsly_provided / Ns;%(Ns*(my(1)-my(end)))
%--------------------------------------------------------------------------
%%%% Figures
%%%% Magnetization Dynamics
figure(10);hold on;grid on;set(gca,'Fontsize',20) 
plot(t,mx,'g','linewidth',3);
plot(t,my,'b','linewidth',3);
plot(t,mz,'r','linewidth',3);
% ylim([-1.1 1.1]);hold on;
% ylabel('m_z,m_x,m_y','Fontsize',16);xlabel('Time [ns]','Fontsize',16)
% %%%%%%%%
%%%% Spin Currents
% figure(2);hold on;grid on;set(gca,'Fontsize',20) 
% plot(t,Iz,'r','linewidth',2);
% plot(t,Ix,'g','linewidth',2);
% plot(t,Iy,'b','linewidth',2);
% ylabel('I_z,I_x,I_y,','Fontsize',16);xlabel('Time [ns]','Fontsize',16)
% %%%%%%%%%
%%%% Slonzchewski currents
figure(31);hold on;grid on;set(gca,'Fontsize',20) 
plot(t,abs(ISLz)*1e3,'r','linewidth',2);
plot(t,ISLx,'g--','linewidth',2);
plot(t,ISLy,'b--','linewidth',2);
plot(t,abs(ISL)*1e3,'k','linewidth',2);
plot(t,abs(Is)*1e3,'b','linewidth',2);
ylabel('I^z_{sl},I_{sl},I_s','Fontsize',20);xlabel('Time [ns]','Fontsize',20)
%xlim([0 .5])
%plot(t,sqrt(ISLz.^2+ISLx.^2+ISLy.^2),'k','linewidth',2);
%plot(t,ISL_prime,'k','linewidth',2);
%%%% Pulse
% % figure(4);hold on;grid on;set(gca,'Fontsize',20) 
% % plot(t,pulse,'k','linewidth',2)
% % ylabel('Pulse','Fontsize',16);xlabel('Time [ns]','Fontsize',16)
% % ylim([0 1.1])

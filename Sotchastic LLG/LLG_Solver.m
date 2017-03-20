% Behtash Behin-Aein
tic
clear;clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global alpha hd...
       max_hfx  max_hfy  max_hfz...
       b_f  f_ramp a_f pwf...
       max_hx_clk_p  max_hy_clk_p  max_hz_clk_p...
       max_hx_clk_d  max_hy_clk_d  max_hz_clk_d...
       b_clk clk_ramp a_clk pw_clk...              
       max_hx_bias_p  max_hy_bias_p  max_hz_bias_p...
       max_hx_bias_d  max_hy_bias_d  max_hz_bias_d...
       b_bias bias_ramp a_bias pw_bias clk_full f_full bias_full;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%                    Parameters                    %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These are the angles at which external sources of torque are applied
% All are defined in the standard spherical coordinate system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Clock spin-torque  %%%%%%%                  %%
theta_clk = 0+0*0.0175;                        %%
phi_clk = 0;                                     %%
%%%  Bias spin-torque   %%%%%%%                  %% 
theta_bias = 0+0*0.0175;                                  %%
phi_bias = 0;                                    %%
%%%  Magnetic field     %%%%%%%                  %%
theta_field = 0;                              %%
phi_field = 0;                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These are magnitudes of various  external  %%%%%%
% fields scaled by Hc (switching field)      %%%%%%
%h_clk_sl = 1.15;                %%
%h_clk_fl = 0.0;                   %%
%%%%%%%%%%%%%%%%%                                %% 
%h_bias_sl = 0.00;% * h_clk_sl;                   %%
%h_bias_fl = 0.0000 * h_clk_fl;                   %%
%%%%%%%%%%%%%%%%%                                %%
h_p = 0.00;                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Initialize the system %%%%%%%%%
mz0 = -.998; 
mx0 = sqrt(1-mz0^2);
my0 = 0;
m0 = [mx0 my0 mz0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants [MKS]
q = 1.6e-19;
hbar = 1.055e-34;
% Magnet parameters %%%%%%%%%%%
g = 1.76*1e7;                        % Gyromagnetic ratio [(rad)/(Oe.s)]
alpha = 0.01;                        % Gilbert damping parameter
Ms = 1000;                            % Saturatin Magnetization [emu/cm^3]
Ku2 = 0.8e6;%1.0*4.144e4; % Uni. anisotropy constant [erg/cm^3]
V = 1*(100*50*2)*1e-21;               % Volume [cm^3]
MsV = Ms*V;                          % [emu]
Ku2V = Ku2*V;                        % [erg]
Ku2V_kT = Ku2V*1e-7/1.6e-19/0.0259   % [kT]
Hc = 2*Ku2/Ms;                        % Switching field [Oe]
hd = 1*4*pi*Ms/Hc;                   % Demagnetizing field [Oe]   
hsc = alpha*(1 + 1*hd/2);
Isc = hsc * Hc * MsV * 1e-7 * 2*q/hbar
%Isc = alpha* Ku2V * 1e-7 * 2*q/hbar
%IslNorm = h_clk_sl * (2*q/hbar)*(MsV*Hc*1e-7)
%IsNorm = 0.125e-3
%CurrentDensity = IsNorm/(40*40*1e-14)/1e6
h_clk_sl =  0*Isc/((2*q/hbar)*(MsV*Hc*1e-7))          
h_clk_fl =  0*50*Isc/((2*q/hbar)*(MsV*Hc*1e-7))
%----------------------------------------------
h_bias_sl =  1.5   *Isc/((2*q/hbar)*(MsV*Hc*1e-7))
h_bias_fl =  0.2*Isc/((2*q/hbar)*(MsV*Hc*1e-7))
%0.9473
Ns = MsV/9.274e-21;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Duration of the Simulation    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NanoS = 10;                       %% Duration in units of nano-seconds %% !!!!!!!!!!!!!!!!!!!
t_span= [0 NanoS*1e-9];
tau_c = (1+alpha^2)/(2*alpha*g*Hc);
tdl = t_span /tau_c ;            %% Dimensionless time span
Ntc = 1e-9/tau_c;      %%    Number of time-constants in a nanosecond
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Clock spin-torque   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% These are the parameters that shape the pulse    [ns/tau_c]
b_clk = 0*Ntc;            %%  Time before application
clk_ramp = 0*Ntc;         %%  Ramp time of the clock
clk_full = 1;             %%  Determines between half/full pulse
a_clk = 0.025*Ntc;            %%  Time after application
pw_clk = tdl(2) - (2*clk_ramp + b_clk + a_clk); %  Pulse width of the clock
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Bias spin-torque   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% These are the parameters that shape the pulse    [ns/tau_c] 
b_bias = 0*Ntc;           %%  Time before application
bias_ramp = 0*Ntc;        %%  Ramp time of the bias
bias_full = 0;            %%  Determines between half/full pulse
a_bias = 0.0*Ntc;           %%  Time after application
pw_bias = tdl(2) - (2*bias_ramp + b_bias + a_bias); %% Pulse width of the bias
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% External Magnetic Field  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% These are the parameters that shape the pulse    [ns/tau_c]
b_f = 0*Ntc;                           %%  Time before application
f_ramp = 0*Ntc;                        %%  Ramp time
f_full = 0;                            %%Determines between half/full pulse
a_f = 0*Ntc;                           %%  Time after application
pwf = tdl(2) - (2*f_ramp + b_f + a_f); %%  Pulse width
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%              External Torques                 %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Beginning      -------------------------------------------------------
% The externally applied magnetic pulse  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_hfx = h_p* sin(theta_field)*cos(phi_field);
max_hfy = h_p* sin(theta_field)*sin(phi_field);
max_hfz = h_p* cos(theta_field);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%
% The externally applied clock spin-torque   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_hx_clk_p = (h_clk_fl - alpha*h_clk_sl) * sin(theta_clk)*cos(phi_clk);                                        
max_hy_clk_p = (h_clk_fl - alpha*h_clk_sl) * sin(theta_clk)*sin(phi_clk);                                        
max_hz_clk_p = (h_clk_fl - alpha*h_clk_sl) * cos(theta_clk);                                        
%%%%%%%%%%
max_hx_clk_d = (h_clk_fl + h_clk_sl/alpha) * sin(theta_clk)*cos(phi_clk);                                        
max_hy_clk_d = (h_clk_fl + h_clk_sl/alpha) * sin(theta_clk)*sin(phi_clk);                                        
max_hz_clk_d = (h_clk_fl + h_clk_sl/alpha) * cos(theta_clk);                                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%
% The externally applied bias spin-torque    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_hx_bias_p = (h_bias_fl - alpha*h_bias_sl) * sin(theta_bias)*cos(phi_bias);                                           
max_hy_bias_p = (h_bias_fl - alpha*h_bias_sl) * sin(theta_bias)*sin(phi_bias);                                           
max_hz_bias_p = (h_bias_fl - alpha*h_bias_sl) * cos(theta_bias);                                           
%%%%%%%%%%
max_hx_bias_d = (h_bias_fl + h_bias_sl/alpha) * sin(theta_bias)*cos(phi_bias);                                           
max_hy_bias_d = (h_bias_fl + h_bias_sl/alpha) * sin(theta_bias)*sin(phi_bias);                                           
max_hz_bias_d = (h_bias_fl + h_bias_sl/alpha) * cos(theta_bias);                                           
%    End      -------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options = odeset('RelTol',1e-13,'AbsTol',[1e-13 1e-13 1e-13]);
[T,M]= ode113('LLG_Solver_DiffEq', tdl, m0, options);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mx=M(:,1)';
my=M(:,2)';
mz=M(:,3)';
dT = diff(T)';
%%%%%%%%%%%
dmxdt = diff(mx)./dT;
dmydt = diff(my)./dT;
dmzdt = diff(mz)./dT;
%%%%%
dmdt_squared = dmxdt.*dmxdt + dmydt.*dmydt + dmzdt.*dmzdt;
%%%%%%%%%%% Intrinsic magnet dissipation %%%%%%%%%%%%%%%%%
Pd = (4*alpha^2/(1+alpha^2))* dmdt_squared;
Dissipated_Energy = trapz(T(1:end-1)',Pd)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Energy delivered by the spin-torque term  %%%%
%%%%%%%%%%% The field is m_cross_Hsl   %%%%%%%%%%%%%%%%%%%
Hx = zeros(1,length(mx));
Hy = mz.*h_clk_sl;
Hz = - my.*h_clk_sl;
Power_Delivered = 2* (Hx(1:end-1).*dmxdt + Hy(1:end-1).*dmydt +...
                                       Hz(1:end-1).*dmzdt);
Energy_Delivered = trapz(T(1:end-1)',Power_Delivered);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4)
bit = M(:,3);
bit(end);
plot(T*tau_c/1e-9,M(:,2),'b','linewidth',1);hold on
plot(T*tau_c/1e-9,M(:,1),'g','linewidth',1);
plot(T*tau_c/1e-9,M(:,3),'r','linewidth',4);hold on;
set(gca,'Fontsize',[24])
xlabel('Time [ns]','Fontsize',[24])
ylabel('Magnetization [Normalized]','Fontsize',[24])
%xlim([0 30])
%plot(T(1:end-1)*tau_c/1e-9,Pd,'r','linewidth',3);hold on;grid on
ylim([-1.1 1.1])











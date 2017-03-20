% Behtash Behin-Aein
clear all;clc
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
       b_clk clk_ramp a_clk pw_clk...              
       max_hx_bias_p  max_hy_bias_p  max_hz_bias_p...
       max_hx_bias_d  max_hy_bias_d  max_hz_bias_d...
       b_bias bias_ramp a_bias pw_bias clk_full f_full bias_full...
       theta_clk phi_clk h_clk_fl h_clk_sl ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%                    Parameters                    %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These are the angles at which external sources of torque are applied
% All are defined in the standard spherical coordinate system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Clock spin-torque  %%%%%%%                  %%
theta_clk = 0*pi-0*pi/2+0*0.0175;                        %%
phi_clk = 0;                                    %%
%%%  Bias spin-torque   %%%%%%%                  %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change pi to 0 for the opposite direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta_bias = 0+0*0.0175;                                  %%
phi_bias = 0;                                    %%
%%%  Magnetic field     %%%%%%%                  %%
load mz_IN
theta_field_a =  acos(mzz);
for aa =  1:length(theta_field_a)
    theta_field = theta_field_a(aa);               %%
phi_field = 0;                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Initialize the system %%%%%%%%%
mz0 = -0.999; 
mx0 = sqrt(1-mz0^2);
my0 = 0;
m0 = [mx0 my0 mz0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These are magnitudes of various  external  %%%%%%
% fields scaled by Hc (switching field)      %%%%%%
h_clk_sl_a = [0.5:0.1:8];
h_clk_sl = 0;
%mzz = zeros(1,length(h_clk_sl_a));
%for ss = 1: length(h_clk_sl_a)   %BBBBBBBBBBBBBB
%    h_clk_sl = h_clk_sl_a(ss);   %BBBBBBBBBBBBBB
%--------------------------------------------------
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Changed
h_clk_fl = 0.0;         % Substitute 0 with H_s
                                % obtained from V_s
                     %(Must be scaled by Hc=2Ku/Ms)
%--------------------------------------------------
%%%%%%%%%%%%%%%%%                                %% 
h_bias_sl = 0.0;% * h_clk_sl;                   %%
h_bias_fl = 0.0 * h_clk_fl;                   %%
%%%%%%%%%%%%%%%%%                                %%
h_p = 1.1;                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Magnet parameters %%%%%%%%%%%
g = 1.76*1e7;                        % Gyromagnetic ratio [(rad)/(Oe.s)]
alpha = 0.01;                        % Gilbert damping parameter
Ms =800;                            % Saturatin Magnetization [emu/cm^3]
Ku2 = 3.0*1e5;%1.0*4.144e4; % Uni. anisotropy constant [erg/cm^3]
V = 1*(50*50*2)*1.5e-21;               % Volume [cm^3]
MsV = Ms*V;                          % [emu]
Ku2V = Ku2*V;                        % [erg]
Ku2V_kT = Ku2V*1e-7/1.6e-19/0.0259;   % [kT]
Hc = 2*Ku2/Ms;                        % Switching field [Oe]
hd = 1*4*pi*Ms/Hc;                   % Demagnetizing field [Oe]        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Duration of the Simulation    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NanoS = 5;                       %% Duration in units of nano-seconds
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
a_clk = 0.0*Ntc;            %%  Time after application
pw_clk = tdl(2) - (2*clk_ramp + b_clk + a_clk); %  Pulse width of the clock
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Bias spin-torque   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% These are the parameters that shape the pulse    [ns/tau_c]
b_bias = 0*Ntc;           %%  Time before application
bias_ramp = 0*Ntc;        %%  Ramp time of the bias
bias_full = 0;            %%  Determines between half/full pulse
a_bias = 0*Ntc;           %%  Time after application
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
[T,M]= ode113('LLG_Solver_DiffEq2', tdl, m0, options);
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
Dissipated_Energy = trapz(T(1:end-1)',Pd);
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
bit = M(:,3);
%mzz(ss) = bit(end);
mzz(aa) = bit(end);
%end
%save mz_IN mzz 


end



figure(4)
%%%%plot(T*tau_c/1e-9,M(:,2),'b','linewidth',3);grid on;hold on
%%%%plot(T*tau_c/1e-9,M(:,1),'g','linewidth',3);
plot(T*tau_c/1e-9,(180/pi)*acos(M(:,3)),'r','linewidth',3);hold on;grid on
%plot(T(1:end-1)*tau_c/1e-9,Pd,'r','linewidth',3);hold on;grid on
ylim([0 180])

figure(5)
bit = M(:,3);
bit(end);
%%%%plot(T*tau_c/1e-9,M(:,2),'b','linewidth',3);grid on;hold on
%%%%plot(T*tau_c/1e-9,M(:,1),'g','linewidth',3);
plot(T*tau_c/1e-9,M(:,3),'g','linewidth',3);hold on;grid on
%plot(T(1:end-1)*tau_c/1e-9,Pd,'r','linewidth',3);hold on;grid on
ylim([-1 1])

figure(6)
plot(h_clk_sl_a,mzz,'b','linewidth',4);hold on;grid on
%plot(h_clk_sl_a,tanh(mzz),'k--','linewidth',4)
%plot(-h_clk_sl_a,-tanh(mzz),'k--','linewidth',4)
ylim([-1 1])


figure(7)
plot(cos(theta_field_a),mzz,'b','linewidth',4);hold on;grid on
%plot(h_clk_sl_a,tanh(mzz),'k--','linewidth',4)
%plot(-h_clk_sl_a,-tanh(mzz),'k--','linewidth',4)
ylim([-1 1])




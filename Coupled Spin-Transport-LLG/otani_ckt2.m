function [h_clk_sl,h_clk_fl,Isl,Islz,Islx,Isly,Ifl,P1,Iz2,Ix2,Iy2]=otani_ckt2(m,Mhat,Ic,MsV,Ku2V)
% function for Self consistent modelling of Otani Circuit with:
% (1) 4 Component spin conductances
% (2) LLG block attached to output FM
% Srikant Srinivasan Sept. 28, 2010

%%% First few lines required for a self contained code rather than a function
%%% call
% clear all; clc;
ii=1;
zdir=[1 0 0];
% var=logspace(-3,0,301);
% % for ii=1:301
% mz0 = 1-var(ii);%1e-4; 
% mx0 = sqrt(1-mz0^2);
% my0 = 0;
% m = [mx0 my0 mz0];
% Ic=6.1e-3;


m=m([ end 1:end-1 ]);
IN(1,1)= Mhat(3,1); 
IN(1,2)= Mhat(1,1);
IN(1,3)= Mhat(2,1); 

% Constants (all MKS, except energy which is in eV)
Z=zeros(4);
q=1.6e-19; h=6.626e-34; hbar=h/2/pi;
% mub=9.274e-24; mu0=4*pi*1e-7;
% g=1.76e11;  % (rad/T/s)      % Gyromagnetic ratio
% alpha=0.007;                 % Gilbert damping parameter
% Ms=780e3;         %A/m       % Saturatin Magnetization
% Ku2=3.14e3;   %J/m^3         % Uni. anisotropy constant
% Vol=2.5*(170*80*4)*1e-27;          % Volume [m^3]
% MsV=Ms*Vol;                  %
% Ku2V=Ku2*Vol;                % J
% Ku2V_kT=Ku2V/1.6e-19/0.0259; % [kT]
Hk=2*Ku2V/MsV;             % Coercive field
% Hd = 4*pi*Ms/Hc;                   % Demagnetizing field [Oe]
%%
%%%%%%% Expt. Ckt. Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Magnet
% Cross section
PF1=0.5;PF2=0.5;
AF1=50*110e-18; lF1=2e-9;   % Area, length of Magnet 1
AF2=50*110e-18; lF2=2e-9;    % Area, length of Magnet 2
lambdaF=5e-9; rhoF=17.1e-8;
RF1=lambdaF*rhoF/AF1; RF2=lambdaF*rhoF/AF2; %Parameters of magnets
LF1=lF1/lambdaF; LF2=lF2/lambdaF; %lF2=4e-9 or 20e-9 % Normalized magnet length
kf=1.36e10; Modes=kf*kf/2/pi; RqF=h/q/q;

% Channel
t=50e-9; AN=100e-9*t;   % thickness, cross sectional area of Channel
lambdaN=.5e-6; rhoN=0.69e-8; RN=lambdaN*rhoN;
RN1=RN/AN; RN2=RN1; RN3=RN1; %RN=RN1;%RN/AN; RN3=RN/AN;
lN2=50e-9; LN2=lN2/lambdaN; LN1=4; LN3=1e-3;

% Gold lead
lambdaG=1e-8;rhoG=7e-8;RGol=lambdaG*rhoG;Rau=RGol/AF1;Lau=5;
%%


%%
%%%%% Spin Ckt Description
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Conductances
% Ferromagnetic contacts
[GBF1,G0BF1]=GBr2(RqF,0,PF1,Modes*AF1,0);
[GBF2,G0BF2]=GBr2(RqF,0,PF2,Modes*AF2,0);

[GF1,G0F1]=GBr2(RF1,LF1,PF1,0,0);
[GF2,G0F2]=GBr2(RF2,LF2,PF2,0,0);

if max(max(GBF1))<max(max(GF1))
    GF1=GBF1; GF2=GBF2;
    G0F1=G0BF1; G0F2=G0BF2;
else
    G0F1=G0F1+G0BF1; G0F2=G0F2+G0BF2;
end

U1=rotmat(zdir,IN);
GF1=U1*GF1*U1';
G0F1=U1*G0F1*U1';


U=rotmat(zdir,m); 
GF2=U*GF2*U';
G0F2=U*G0F2*U';

% Normal channel
[GN1,G0N1]=GBr2(RN1,LN1,0,0,0);
[GN2,G0N2]=GBr2(RN2,LN2,0,0,0);
[GN3,G0N3]=GBr2(RN3,LN3,0,0,0);

% Top Gold Contacts
[GA1,G0A1]=GBr2(Rau,Lau,0,0,0);
[GA2,G0A2]=GBr2(Rau,Lau,0,0,0);


% Non-local computation
% Conductance matrix
G=[G0A1+GA1 -GA1 Z Z Z Z Z;
    -GA1 G0A1+GA1+G0F1+GF1 -GF1 Z Z Z Z;
    Z -GF1 G0F1+GF1+GN2+G0N2+G0N1+GN1 -GN2 Z Z Z ;
    Z Z -GN2 G0N2+GN2+G0F2+GF2+GN3+G0N3 -GF2 Z -GN3;
    Z Z Z -GF2 GF2+G0F2+GA2+G0A2 -GA2 Z;
    Z Z Z Z -GA2 GA2+G0A2 Z;
    Z Z Z -GN3 Z Z GN3+G0N3];

C = [Ic;zeros(27,1)];%  Terminal currents
V=G\C;V=reshape(V,4,7);% Terminal voltages
% delV=V(1,6)-V(1,7); %Non-Local voltage measured
IF2=    -      GF2*(V(:,4)-V(:,5))    -     (G0F2)*V(:,4);% current entering FM2
%%

%%%% Extracting spin currents, Field and Slonczewski components
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Iz2(ii)=IF2(2);
Ix2(ii)=IF2(3);
Iy2(ii)=IF2(4);
% -------------
Is=IF2(2:4);
%if norm(cross(m,Is))== 0
%    m = [-0.999999 sqrt(1-0.999999^2) 0];
%end
Ist = -cross(m,cross(m,Is));
FLhat = - cross(m,IN)/norm(cross(m,IN));
SLhat = - cross(m,cross(m,IN))/norm(cross(m,cross(m,IN)));
Isl = dot(Ist,SLhat);
Islv = dot(Ist,SLhat) * SLhat;
Islz = dot(Islv,[1 0 0]);
Islx = dot(Islv,[0 1 0]);
Isly = dot(Islv,[0 0 1]);
% ----------------------
Ifl=dot(Ist,FLhat); 
I_H_conv = hbar/2/q/(MsV*Hk*1e-7);
h_clk_sl(ii)= Isl * I_H_conv / norm(cross(m,cross(m,IN)));
h_clk_fl(ii)= Ifl * I_H_conv / norm(cross(m,IN));
P1=Ic*V(1,1);

end


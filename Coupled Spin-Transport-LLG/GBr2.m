function [Gmat, G0mat] = GBr(R,L,P,f,ang)
% Srikant Srinivasan
%This function written on Aug 31,2010 for Multiterminal_ASL_4x4.m
%R=rho*lambda, L=L/lambda, P=polarization (fraction),
% theta, phi magnetization vector: theta wrt z and phi wrt x
% f=mixing conductance fraction, ang=mixing angle.
a=f*cos(ang);b=f*sin(ang);
MS=diag([0 1 1 1]);
Gmat=[1 P 0 0;
    P P^2 0 0;
    0 0 0 0;
    0 0 0 0];
if f==0
    % FM, NM etc
    if L==0
        %tunnel barrier
        Gmat=(1/R)*(Gmat+(1-P^2)*MS);
        G0mat=[];
    else
        if P==0
            % NM
            Gmat=(1/R/L)*(Gmat+(1-P^2)*L*csch(L)*MS);
            G0mat=1/R*tanh(L/2)*MS;
        else
            % FM
            % Gmat=(1/R/L)*(Gmat+(1-P^2)*L*csch(L)*MS);
            Gmat=(1/R/L)*(Gmat+(1-P^2)*L*csch(L)*diag([0 1 0 0]));
            G0mat=(1-P^2)/R*tanh(L/2)*diag([0 1 0 0]);
        end
    end
else
    %Bauer Interface shunt
    Gmat=f/R*[1 P 0 0; P 1 0 0; 0 0 0 0; 0 0 0 0];
%     Gmat=f/R*[1 P 0 0; P 1 0 0; 0 0 1 0; 0 0 0 1];
    G0mat=1/R*[0 0 0 0; 0 0 0 0; 0 0 a b; 0 0 -b a];
end
end
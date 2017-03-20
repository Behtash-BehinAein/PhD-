% Suppose Is is the spin current entering the magnet
% calculated using NEGF. These lines will then calculate 
% the quantities that will go into the LLG code.
% Note that you need to import MsV and Hc based on 
% magnet parameters and also M based on the direction of the 
% injecting magnet. 

Ist = -cross(m,cross(m,Is));
FLhat = - cross(m,M)/norm(cross(m,M));
SLhat = - cross(m,cross(m,M))/norm(cross(m,cross(m,M)));
% ----------------------------
Isl = dot(Ist,SLhat);
Ifl=dot(Ist,FLhat); 
% ----------------------------
I_H_conv=hbar/2/q/(MsV*Hc*1e-7);
h_clk_sl= Isl * I_H_conv / norm(cross(m,cross(m,M)));
h_clk_fl= Ifl * I_H_conv / norm(cross(m,M));
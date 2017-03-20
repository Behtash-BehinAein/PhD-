function [MN_initial] = Initial_Magnetization1DArray(N)
%  Up %%%%%%%%
mu = [sqrt(1-0.9^2);-0.9;0];
%%%%%%%%%%%%%%
%  Down  %%
md = [sqrt(1-0.5^2);0.5;0];
%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Column vector to hold the x,y,z components of
% magnetic moments of individual magnets
MN_initial =zeros(3*N,1);
cc=1;
for nn = 1 : N
    a=mod(nn,2);
    if a~=0
        MN_initial(cc:cc+2) = mu;    
    else
        MN_initial(cc:cc+2) = md;
    end        
    cc = cc+3;
end
%MN_initial(1:3) = [sqrt(1-1^2);-1;0];
% MN_initial(1:3) = [-1;0;0];
% MN_initial(4:6) = [1;0;0];
% MN_initial(7:9) = [0;-1;0];
% MN_initial(10:12) = [0;1;0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
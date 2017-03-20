function  [rx ry rz r A]= coordinate_generator(N,a)
%  Spatial Coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%
%Nx = 3;           % # of magnets along the x-axis
%Ny = 2;            % # of magnets along the y-axis
%Nz = 0;            % # of magnets along the z-axis 
%N=Nx+Ny+Nz;
%N=2;
% This is for a cubic structure %%%%%%%%%%
% Spacing between nano-magnets [CGS] %%%%%
%a =  24 *5e-7;
%a =  100e-7;              %Along the x-axis  
%b =  34 *5e-7;            %Along the y-axis 
%e =  70e-7;
%c =  23*1e-7;             %Along the z-axis  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Spatial structure
%       Note that the magnets along x-axis
%                   will occupy the origin
%               y-axis
%                 o  
%                 o 
%                 * * *  x-axis
%               ^
%             ^
%         z-axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%XCx = -a: a : (Nx-2)*a;  
%XCy=  zeros(1,length(XCx));
%XCz = zeros(1,length(XCx));
%%%%%
%YCy = -b : 2*b :(Ny-1)*b;
%YCx = zeros(1,length(YCy));
%YCz = zeros(1,length(YCy));
%%%%%
%ZCz = c : c : Nz*c;
%ZCx = zeros(1,length(ZCz));
%ZCy = zeros(1,length(ZCz));
%%%%%
% X coordinate of all magnets %%%%%%%%%%%%
%CNx = [XCx YCx ZCx];
%CNx = [-a 0 a    0     0     -e     -(e+a)    -e];
CNx = (0:a:(N-1)*a);
% Y coordinate of all magnets %%%%%%%%%%%%
%CNy = [XCy YCy ZCy];
%CNy = [ 0 0 0   -b     b    e+b       e    -(e+b)];
CNy = zeros(1,N);
% Z coordinate of all magnets %%%%%%%%%%%%
CNz = zeros(1,N);%[XCz YCz ZCz];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Row vector holding magnet numbers
RN = 1:N;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Matrices holding all circular shifts  %%%
A  = zeros(N,N);
A(1,:) = RN;
cx = zeros(N,N);
cx(1,:) = RN;
cy = zeros(N,N);
cy(1,:) = RN;
cz = zeros(N,N);
cz(1,:) = RN;
Cx = zeros(N,N);
Cx(1,:) = RN;
Cy = zeros(N,N);
Cy(1,:) = RN;
Cz = zeros(N,N);
Cz(1,:) = RN;
%%%%
%%%%
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii=1:N-1    
    A (ii+1,:)= circshift(RN,[0 ii]);
    cx(ii+1,:)= circshift(CNx,[0 ii]);
    cy(ii+1,:)= circshift(CNy,[0 ii]);
    cz(ii+1,:)= circshift(CNz,[0 ii]);   
    Cx(ii+1,:)= CNx;
    Cy(ii+1,:)= CNy;
    Cz(ii+1,:)= CNz;    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Matrices holding coordinates of r (r is the %% 
% vector connecting two magnets) for pairwise
% interaction of magnets           %%%%%%%%%%
rx = Cx - cx;
ry = Cy - cy;
rz = Cz - cz;
r = sqrt(rx.^2 + ry.^2 + rz.^2);
r(1,:)=5e-7*ones(1,N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R= rotmat(a,b)
%Srikant Srinivasan
% Defining a rotation matrix to rotate output magnet from a to b
% clear all; clc
% a=[1;0;0]; b=[0;1;0];
a=a/norm(a); b=b/norm(b);
%theta=acos(dot(a,b)); c=cos(theta); s=sin(theta);
c=dot(a,b); s=sqrt(1-c^2);
if s==0
    u=[0 0 0];
else
    u=cross(a,b)/norm(cross(a,b));
end
%%% X,Y,Z coordinate system
% ux=u(1); uy=u(2); uz=u(3);
% R=[ux^2+(1-ux^2)*c      ux*uy*(1-c)-uz*s    ux*uz*(1-c)+uy*s;
%    ux*uy*(1-c)+uz*s     uy^2+(1-uy^2)*c     uy*uz*(1-c)-ux*s;
%    ux*uz*(1-c)-uy*s     uy*uz*(1-c)+ux*s    uz^2+(1-uz^2)*c]

%%% Z,X,Y coordinate system
ux=u(2); uy=u(3); uz=u(1);
R=[1        0                   0                   0;
    0 uz^2+(1-uz^2)*c      uz*ux*(1-c)-uy*s    uz*uy*(1-c)+ux*s;
    0 uz*ux*(1-c)+uy*s     ux^2+(1-ux^2)*c     ux*uy*(1-c)-uz*s;
    0 uz*uy*(1-c)-ux*s     ux*uy*(1-c)+uz*s    uy^2+(1-uy^2)*c];
end

% R=[uz^2+(1-uz^2)*c      uz*ux*(1-c)-uy*s    uz*uy*(1-c)+ux*s;
%    uz*ux*(1-c)+uy*s     ux^2+(1-ux^2)*c     ux*uy*(1-c)-uz*s;
%    uz*uy*(1-c)-ux*s     ux*uy*(1-c)+uz*s    uy^2+(1-uy^2)*c];
% R*b
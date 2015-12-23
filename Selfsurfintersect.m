function [ D tt] =  Selfsurfintersect( n1 )

tt=-2+4*rand(n1,1);     % Using  spherical coordinates

x=tt.^3-3*tt;
y=-tt.^2;
z=3*(-.5+rand(n1,1));
D=[x y z];  




% n1=2500;
% teta=(2*pi)*randn(n1,1);     % Using  spherical coordinates
% x=cos(3*teta).*sin(teta);
% y=cos(3*teta).*cos(teta);
% D2=[x y];

%  D=[D1;D2;D3];

end

% [group,path]=Path_Based_Cluster_LandMarks( D ,30 , 70 , 15 , 2 );
% scatter3(D(:,1),D(:,2),D(:,3),'.','b');
% axis equal
% view(-145,45)
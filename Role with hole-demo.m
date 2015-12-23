clc
clear all
rng(0)
n=1000;            % number of data points
knn=60;           % number of nearest neighbors
ep=.3;             % parameter for epsilon graph
ang=90;             % angle or cervature constrained
kk=2;              % desired dimension of intrinsic space
%%%% m is the number of points
[ D th] = Roll(n);
n=size(D,1);
ll=n;              % number of landmarks
% % load X;
% load five_affine_subspaces
% D=X';
noiseLevel=0.001;
noise = noiseLevel * rand(size(D));
D = D + noise;
%%%%%%%%%%%%%%%Plot the original figure%%%%%%%%%%%
% [a b]=sort(D(:,1),1);
% DD=zeros(m,2);
% for i=1:m
%     DD(i,:)=D(b(i),:);  % sorting data points based on x-axis
% end
% D=DD;

colormap jet;
color=th;
figure(1);
scatter3(D(:,1), D(:,2),D(:,3), 10,color, 'filled');
%%%%%%%%%%%%%%%%%%% build in KNN graph %%%%%%%
 [edge_matrix, weights]=Build_KNN(D,knn);        %  KNN grapoh
% [edge_matrix, weights]=Build_Epsilon(D,ep);    %  epsilon graph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sn=(1:n)';
%sn=ceil(n*rand(ll,1));   % landmark points
%%%%%%%%%%%%%%%%%%%%  Apply constraint shortest path %%%%%%%%%%%%%
 % [Binary_matrix,path] = Constraint_Dijkstra(D,edge_matrix,weights,sn,ang);
 [Binary_matrix,path] = Constraint_Dijkstra(D,edge_matrix,weights,sn,ang);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dist=zeros(ll,n);
for i=1:ll
    for j=1:n
        h=path{i,j};
        v=size(h,2);
        if v~=1
                    for t=1:v-1
            dd=norm(D(h(t),:)-D(h(t+1),:));  %%% distance matrix
            dist(i,j)=dist(i,j)+dd;
                    end
        else
            dist(i,j)=inf;
        end
          dist(i,sn(i))=0;
    end
end

 dd=(dist+dist')/2;
 dd(dd>40)=nan;
    
    
Y = mdscale(dd,2);
figure (3);
scatter(Y(:,1),Y(:,2),10,color, 'filled');
title('Output of MDS')
axis equal
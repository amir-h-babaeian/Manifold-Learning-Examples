clc
clear all
n=8000;           % number of data points
ll=4;              % number of landmarks
knn=60;           % number of nearest neighbors
ep=.3;             % parameter for epsilon graph
ang=17;            % angle or cervature constrained
kk=2;              % desired dimension of intrinsic space
rng(0);
%%%% m is the number of points
[ D th] = Selfsurfintersect(n);
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
axis equal
%%%%%%%%%%%%%%%%%%% build in KNN graph %%%%%%%
 [edge_matrix, weights]=Build_KNN(D,knn);        %  KNN grapoh
% [edge_matrix, weights]=Build_Epsilon(D,ep);    %  epsilon graph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sn=ceil(n*rand(ll,1));   % landmark points
%%%%%%%%%%%%%%%%%%%%  Apply constraint shortest path %%%%%%%%%%%%%
 % [Binary_matrix,path] = Constraint_Dijkstra(D,edge_matrix,weights,sn,ang);
 [Binary_matrix,path] = Constraint_Dijkstra(D,edge_matrix,weights,sn,ang);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2);
scatter3(D(:,1),D(:,2),D(:,3),'.','b');
hold on            
scatter3(D(sn,1),D(sn,2),D(sn,3),'*','R');  %%% indicates Landmarks

hold on;
hh=path{3,3451};
plot3(D(hh,1),D(hh,2),D(hh,3),'--rs','LineWidth',2,...
                'MarkerEdgeColor','g',...
                'MarkerFaceColor','g',...
                'MarkerSize',5);
hold on;
hh=path{2,7};
plot3(D(hh,1),D(hh,2),D(hh,3),'--rs','LineWidth',2,...
                'MarkerEdgeColor','g',...
                'MarkerFaceColor','g',...
                'MarkerSize',5);
hold on;
hh=path{3,50};
plot3(D(hh,1),D(hh,2),D(hh,3),'--rs','LineWidth',2,...
                'MarkerEdgeColor','g',...
                'MarkerFaceColor','g',...
                'MarkerSize',5);
hold on;
hh=path{4,30};
plot3(D(hh,1),D(hh,2),D(hh,3),'--rs','LineWidth',2,...
                'MarkerEdgeColor','g',...
                'MarkerFaceColor','g',...
                'MarkerSize',5);            
                        
            
hold off           


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Landmark MDS %%%%%%%%%%%%%%%%%%%%%%%
DL=zeros(ll,ll);
for i=1:ll
    for j=1:ll
        DL(i,j)=dist(i,sn(j)); %% Distance matrix of landmark points
    end
end
DL=(DL+DL')/2;

%%%% applying classical MDS
DLTA=DL.^2;
Mi=mean(DLTA,2);
M=mean(Mi);
B=zeros(ll,ll);
for i=1:ll
    for j=1:ll
        B(i,j)=-(DLTA(i,j)-Mi(i)-Mi(j)+M)/2;
    end
end
[V,lambda] = eig(B);             %% eigen value decomposition
d=diag(lambda);
[a b]=sort(d,'descend');         %% sort the eigenvalues
num=size(d(d>0),1);              %% number of positive eigen values;
k=min(kk,num);                   %% take the minumum of our dsired out put  
   % dimension and and number of positive eigen values as the dimension of
   % the intrinsic space
                
%%%%%%        CMDS MATLAB VERSION         %%%%%%%%%%%%%%%%%%%%            
% [V,a] = cmdscale(DLTA);
% num=size(a(a>0),1);            %% number of positive eigen values;
% k=min(kk,num);                 %% take the minumum of our dsired out put  
% b=1:1:k;



% [V,a] = mdscale(D,2);
% % num=size(a(a>0),1);            %% number of positive eigen values;
% % k=min(kk,num);                 %% take the minumum of our dsired out put  
% % b=1:1:k; 
% 


L=zeros(k,ll);
for i=1:k
    L(i,:)=sqrt(a(i))*V(:,b(i))';
end
%%%% finding the rest of points using triangulation
Linv=zeros(k,ll);
for i=1:k
    Linv(i,:)=(1/sqrt(a(i)))*V(:,b(i))';
end

Mn=mean(DLTA,2);
X=zeros(k,n);
DLTdist=dist.^2;
for j=1:n
    X(:,j)=-Linv*(DLTdist(:,j)-Mn)/2;
end
%%%%%%%%%%%%%% Display the embeded points
% Y=zeros(m,1);
X=X';
figure(3);
scatter(X(:,1),X(:,2), 10 ,color, 'filled');
axis equal

%%%%%%%%%%%%% Take care of outliers(points on the intersection) %%%%%%%%%
kt=10;
[n1,d1] = knnsearch(D,D,'K',kt+1,'Distance','euclidean');
edge=n1(:,2:end);             % matrix of conected edges to each node
weights=d1(:,2:end);          % matrix of distance between conected nodes
epsilon=.3;
%  [a, b]=sort(X(:,2));
%  dif=a-[0;a(1:end-1)];
Y=zeros(n,2);
X(isnan(X))=inf;
tt=0;
for i=1:n
    if  -1.5-epsilon<X(i,2)&& X(i,2)<1.5+epsilon
        Y(i,:)=X(i,:);
    else
        td=[0 0];
        nn=0;
        tt=tt+1;
        outD(tt)=i;                 %%% points embedded uncorrectly
        for j=1:kt
            if -.5-epsilon<X(edge(i,j),2) && X(edge(i,j),2)<.5+epsilon
                td=td+X(edge(i,j),:);
                nn=nn+1;
            end
        end
        Y(i,:)=td/nn;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4);
scatter(Y(:,1),Y(:,2), 10 ,color, 'filled');
axis equal

%%%%%%%% Show points that are not embedded very well %%%%%
figure(6);
scatter3(D(:,1),D(:,2),D(:,3),'.','b');
hold on
s=40;
scatter3(D(outD,1),D(outD,2),D(outD,3),s,'r','fill');  %%% indicates Landmarks
hold off

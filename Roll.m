
function [ D th] =  Roll( N )

fact=1; % Length of roll


tt = sort((fact*pi/2)*(1+2*rand(1,N)));  height = fact*7*rand(1,N);
X = [tt.*cos(tt); height; tt.*sin(tt)];
D=X-repmat(mean(X,2),1,N);
D=D';


h1=find(1.5<D(:,2));
h2=find(-1.5>D(:,2));
h3=find(D(:,1)>-1);
h4=[h1 ; h2 ; h3];
h5= unique(h4);
th=tt(h5);
D=D(h5,:);


end

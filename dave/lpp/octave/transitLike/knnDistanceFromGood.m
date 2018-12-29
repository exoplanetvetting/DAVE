function [ dymean] = knnDistanceFromGood( x,y,knn )
%x is a n x dim  matrix of those objects known to be of the same class
%   y are  n2 x dim matrix of those objects wanting to be classified
% knn is the integer number of nearest neighbors.
% returns the mean distance to the classified objects (x) for each n2 in y
% Also returns the 1 sigma standard deviation of performing the same knn
% test on the x set.

%Perform knn on the x set.
%[n,dx]=knnsearch(x,x,'k',knn,'distance','minkowski','p',3);
%dxmean=mean(dx');
%dxstd=std(dxmean);
%dxmax=max(dxmean);

[n,dy]=knnsearch(x,y,'k',knn,'distance','minkowski','p',3);
dymean=mean(dy');

% figure;
% hist(dxmean,60);
% 
% 
% figure;
% [ny,by]=hist(dymean,floor(length(y)/5));
% ny=ny/max(ny);
% [nx,bx]=hist(dxmean,60);
% nx=nx/max(nx);
% bar(by,ny,'k');
% hold on
% bar(bx,nx,'r');
% hold off
% xlim([0 dxstd*15]);
end


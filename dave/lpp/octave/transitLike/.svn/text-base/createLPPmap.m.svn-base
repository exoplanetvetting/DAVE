function [ Ymap,Zmap ] = createLPPmap(sample,want,nDim,nneighbors)
%Input a sample Matrix from createMatrixByType with the .X field
%Map using LPP based on those in the "want"ed sample.
%Do this for both 2 dimensions and nDim
%Based on using nneighbors.
%Return an array for the n Dimensions and 2 dimensions.

X=sample.X;
%want=sample.mes>mesLimit & sample.d==dtype;

[mappedX,mappingX]=compute_mapping(X(want,:), 'LPP',nDim,nneighbors);
[mappedY,mappingY]=compute_mapping(X(want,:),'LPP',2,nneighbors);

%Compute mapping for all in the sample
[Yall]=maplle_oos(sample.X,mappingX,nDim);
[Zall]=maplle_oos(sample.X,mappingY,2);

Ymap={};
Ymap.nDim=nDim;
Ymap.mapping=mappingX;
Ymap.mapped=Yall;
Ymap.want=want;


Zmap={};
Zmap.nDim=2;
Zmap.mapping=mappingY;
Zmap.mapped=Zall;
Zmap.want=want;

end




% %Compute the standard deviations nad 68% contours for each dimension as a function of mes.
% %Find bins based on the lowest 90% of mes in log space.
% %Limits shoudl be based only on the types of PC and EB.
% 
% Yuse=Yall(good,:);
% Zuse=Zall(good,:);
% goodmes=sample.mes(good);
% 
% sortmes=sort(goodmes);
% swant=sample.mes<sortmes(end-floor(0.15*length(goodmes)));
% [h,bins]=hist(log10(sample.mes(swant)),4);
% bins=10.^bins;
% 
% [stdY, meanY, meanStatsY,minY,maxY,upperY,lowerY] = lpp_similartest(Yuse,goodmes, bins, nDim);
% [stdZ, meanZ, meanStatsZ,minZ,maxZ,upperZ,lowerZ] = lpp_similartest(Zuse, goodmes, bins, 2);
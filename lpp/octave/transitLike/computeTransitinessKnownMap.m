function [ outputInfo ] = computeTransitinessKnownMap( info,map,outfile )
%Using a previous map and population (map) of known transits
%Compute the lpp-knn transit metric on a new set of folded bin light curves
%map is the result of the code compute_Transitiness.m
%info is a result of the code createMatrixByType.m
%
%It then prints the reslts to a file one line per TCE.
%This gets both the number of dimensions and the number near neighbors from
%the map struct.
%The map struct also contains a vector of which elements to match
%to.(map.knnGood)

%Apply the LPP map to the old sample
%[Ygood]=maplle_oos(map.X,map.Ymap,map.nDim);
%[Zgood]=maplle_oos(map.X,map.Zmap,map.nDim);
Yorig=map.Ymap.mapped;
Zorig=map.Zmap.mapped;

%Apply the LPP map for the out of sample
[Yall]=maplle_oos(info.X,map.Ymap.mapping,map.nDim);
[Zall]=maplle_oos(info.X,map.Zmap.mapping,2);

%x are known transits
%y are those that need to be classified
x=Yorig(map.knnGood,:); 
y=Yall; 

[ dymean, dxmean, dxstd, dxmax ] = knnDistance_inclass( x,y,map.knn );

fid=fopen(outfile,'w');

fprintf(fid,'#Date = %s\n#1SigmaDistance = %f\n#NDim = %i\n#knn = %i\n#type =%s \n',date,dxstd,map.nDim,map.knn,info.dettype);
fprintf(fid,'#TCE    MeanDistance   sampleType\n');
for i=1:length(y)
    
    fprintf(fid,'%s  %f  %i\n', info.tce{i}, dymean(i),info.d(i));
    
end

fclose(fid);

outputInfo=struct([]);
outputInfo(1).Ymap=Yall;
outputInfo.Zmap=Zall;
outputInfo.knn=map.knn;
outputInfo.outfile=outfile;
outputInfo.dymean=dymean;
outputInfo.dxmean=dxmean;
outputInfo.dxstd=dxstd;
outputInfo.dxmax=dxmax;
outputInfo.transitMetric=dymean;
outputInfo.transit1sigmacut=dxstd;
outputInfo.nDim=map.nDim;
outputInfo.svnVersion='$Id: computeTransitinessKnownMap.m 58926 2015-04-09 18:26:39Z sethomps $';

oi=catstruct(outputInfo,info);

outputInfo=oi;

end


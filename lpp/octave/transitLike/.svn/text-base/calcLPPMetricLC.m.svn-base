function [ Tlpp, Y, binnedFlux ] = calcLPPMetricLC(dataStruct, mapStructFile )
%This function is meant to calaculate the LPP transit metric on one
%individual light curve. This has hard coded into it, the values used to
%bin the light curve for the Q1-Q17 DR24 calculation of the LPP transit
%metric.  Be sure to give the mapping that corresponds to your detrending.
%
%dataStruct Contains
% time - vector, days 
% flux - vector, Detrended counts
% duration - float, hours 
% period - float, days
% phase - float, BKJD 
% mapStructFile - string containing a struct containing
%               map with:
%                 - nDim number of reduced dimensions (integer)
%                 - mapping   - M (Nxn), mean(1xN), name='LPP' (matlab
%                             dimensionality reduction toolbox)
%                 - knnGood -- vector of logicals for knn step (good
%                 transits) (not the 'want' vector)
%                 - mapped -- original vectors in the mapping
%                 - knn
%%%

ntrfr=5;
npts=100;
knn=15;

[ binnedFlux ] = foldBinLightCurve( dataStruct,ntrfr,npts);

load(mapStructFile);

[ Tlpp, Y ] = computeOneTransMetricKnownMap(binnedFlux',map );


end


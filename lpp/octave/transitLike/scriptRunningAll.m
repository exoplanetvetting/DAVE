%Steps to create the transitiness metric.
%

other(1).dir='/soc/nfs/so-nfs/Q17/TrapezoidFitPDC/DATS/';
other(2).dir='/soc/nfs/so-nfs/Q17/DVSummaries/DATS/';
other(1).name='TrapezoidDet';
other(2).name='TpsMedDet';
other(1).tag='-fulltime.dat';
other(2).tag='-fulltime.dat';

storeDir='/home/sethomps/Kepler/RoboVetter/Q1Q17knnlpp/DataStore/';
tceList='tces.txt'; %'tces.txt';  %Input File
dvFileList='getdvpaths.txt';
%tceList='smalldetrtest.txt'
order=2;

[ output ] = createMedPDCLightCurves( dvFileList, tceList, storeDir, order, other,1 );


%%
%Create the sample based on certain inputs
%Then set aside 10% for testing.
%type is an integer with 1= unknown  2=pc 3=s-flag 4=c-flag 5=e-flag
%6=n-flag 7=no flag : then code to add 100 to the type to create the test set from a
%random population of randPer percentage of the sample.

federationFile='koimatch_Q17ops_FINAL_112614.txt'; %Only call federated best matches =1
cumulativeTable='cumulative.psv';  %psv created from makeCumulativeTable on svn
outfile='labeledTceListDR24.dat';
tceFile='tces.txt';
triagefails='triagefailsQ12Q16.txt';
randPer=10; %Use randPer to set aside some percentage for testing.

types=createLabeledTCEList(tceFile,federationFile,cumulativeTable,triagefails,outfile,randPer);


%%

tceFile='/home/smullall/Kepler/RoboVetter/Q1Q17knnlpp/labeledTceListDR24.dat';
storedir='/home/sethomps/Kepler/RoboVetter/Q1Q17knnlpp/DataStore/';
npts=100;
ntrfr=2; %was 5 for production run
%over=10;

type='TrapezoidDet';
[infoTrap] = createMatrixByType(tceFile,storedir,npts,ntrfr,type);

type='PdcMedDet';
[infoPdc] = createMatrixByType(tceFile,storedir,npts,ntrfr,type);

type='TpsMedDet';
[infoTps] = createMatrixByType(tceFile,storedir,npts,ntrfr,type);


outfile=['inputQ1Q17DR24_Matrix',num2str(datenum(date)-730000),'.mat'];
save(outfile,'infoTrap','infoPdc','infoTps')

%%
%In case you needed to run createLabeledTceList more than once.
%infoPdc.d=types;
%infoTps.d=types;
%infoTrap.d=types;

%%
%Do the mapping and create a transitiness metric using Pdc light curves
info=infoPdc;
nDim=20;
knn=15;
types=info.d;
%want=(types==2 | types==3 | types==4 | types==5 | types==7) & info.mes>7.8;
want=info.mes>8;
length(want(want))
savefile=['mappingQ1Q17DR24pdc_',num2str(datenum(date)-730000),'.mat'];

koilike=(types==2 | types==3 | types==4 | types==5 | types==7);

outfile=['lppknnDistancePdc_',info.dettype,num2str(datenum(date)-730000),'.dat'];

[ outputInfo ] = compute_Transitiness( info,want,koilike,nDim,knn,outfile );

save(savefile,'outputInfo')
%
%Set limit and run test
nsig=4;
decision=outputInfo.transitMetric > nsig* outputInfo.transit1sigmacut;
d=info.d;
print_tmetricResults(decision,d','Test Results for 4 sigma cut PDC')

%Create a plot of mapping and transitMetric
y=log10(outputInfo.transitMetric);
x=log10(outputInfo.mes);
interactive_plotting(x,y,outputInfo)
xlabel('log(MES)')
%figure
y=log10(outputInfo.transitMetric);
x=log10(outputInfo.periods);
interactive_plotting(x,y,outputInfo)
xlabel('log(Period)')

%%
%Do the mapping and create a transitiness metric using Tps light curves
info=infoTps;
%info.X=info.Xa;
nDim=20;
knn=15;
types=info.d;
%want=(types==2 | types==3 | types==4 | types==5 | types==7) & info.mes>7.8;
want=info.mes>8;% & info.mes<30;
length(want(want))
savefile=['mappingQ1Q17DR24dv_',num2str(datenum(date)-730000),'.mat'];

koilike=(types==2 | types==3 | types==4 | types==5 | types==7);

outfile=['lppknnDistanceDv_',info.dettype,num2str(datenum(date)-730000),'.dat'];

[ outputInfo ] = compute_Transitiness( info,want,koilike,nDim,knn,outfile );

save(savefile,'outputInfo')
%%
%Set limit and run test
nsig=3;
meanpc=median(outputInfo.transitMetric(outputInfo.knnGood));
decision=outputInfo.transitMetric > nsig* outputInfo.transit1sigmacut+meanpc;
d=outputInfo.d';
print_tmetricResults(decision,d,'Test Results for 3 sigma cut TPS')

%Create a plot of mapping and transitMetric
y=log10(outputInfo.transitMetric);
x=log10(outputInfo.mes);
interactive_plotting(x,y,outputInfo)
xlabel('log(MES)')
%figure
y=log10(outputInfo.transitMetric);
x=log10(outputInfo.periods);
interactive_plotting(x,y,outputInfo)
xlabel('log(Period)')
%%
%Plot a normalized histogram
figure
outputInfo.dall=outputInfo.d;
outputInfo.d=outputInfo.d(~isinf(log10(outputInfo.transitMetric)));
outputInfo.logTM=log10(outputInfo.transitMetric(~isinf(log10(outputInfo.transitMetric))));
histogramParamLabels(outputInfo,'logTM',100)
title('DV LPP Metric')
xlabel('DV LPP Metric')
ylabel('Normalized Histogram')
%print('-dpng','DVMetricHist.png')
%%
%Do the mapping and create a transitiness metric using Trap light curves
info=infoTrap;
nDim=20;
knn=15;
types=info.d;
%want=(types==2 | types==3 | types==4 | types==5 | types==7) & info.mes>7.8;
want=info.mes>8;
length(want(want))
savefile=['mappingQ1Q17DR24trap_',num2str(datenum(date)-730000),'.mat'];

koilike=(types==2 | types==3 | types==4 | types==5 | types==7);

outfile=['lppknnDistanceTrap_',info.dettype,num2str(datenum(date)-730000),'.dat'];

[ outputInfo ] = compute_Transitiness( info,want,koilike,nDim,knn,outfile );

save(savefile,'outputInfo')
%
%Set limit and run test
nsig=4;
decision=outputInfo.transitMetric > nsig* outputInfo.transit1sigmacut;
d=info.d';
print_tmetricResults(decision,d,'Test Results for 4 sigma cut Trap')

%Create a plot of mapping and transitMetric
y=log10(outputInfo.transitMetric);
x=log10(outputInfo.mes);
interactive_plotting(x,y,outputInfo)
xlabel('log(MES)')
%figure
y=log10(outputInfo.transitMetric);
x=log10(outputInfo.periods);
interactive_plotting(x,y,outputInfo)
xlabel('log(Period)')
%%
%Let's try to compare values between the types of detrending.
load('mappingQ1Q17DR24trap_6084.mat');
outputInfoTrap=outputInfo;
load('mappingQ1Q17DR24dv_6083.mat');
outputInfoDV=outputInfo;
wantTrap=zeros(length(outputInfoTrap.tce),1);
wantDV=zeros(length(outputInfoDV.tce),1);

for i=1:length(outputInfoTrap.tce)
    atce=outputInfoTrap.tce{i};
    [trapi]=find(strcmp(atce,outputInfoDV.tce));
    if length(trapi)>=1
        wantTrap(i)=1;
        wantDV(trapi)=1;
    end
end
commontces=outputInfoTrap.tce(wantTrap==1);

common={};
common.tce=commontces;
common.DVmetric=outputInfoDV.transitMetric(wantDV==1);
common.Trapmetric=outputInfoTrap.transitMetric(wantTrap==1);
common.periods=outputInfoDV.periods(wantDV==1);
common.mes=outputInfoDV.mes(wantDV==1);
common.d=outputInfoDV.d(wantDV==1);
common.dur=outputInfoDV.durs(wantDV==1);
common.XTrap=outputInfoTrap.X(wantTrap==1,:);
common.XDV=outputInfoDV.X(wantDV==1,:);
common.XaDV=outputInfoDV.Xa(wantDV==1,:);
common.XbDV=outputInfoDV.Xb(wantDV==1,:);
common.wsmes=outputInfoTrap.wsmes(wantTrap==1);
common.ses=outputInfoTrap.ses(wantTrap==1);
common.dettype='common';
common.svnVersion=outputInfoDV.svnVersion;
common.X=[common.DVmetric;common.Trapmetric;common.wsmes/(max(common.wsmes));common.mes/max(common.mes);common.dur./(common.periods*24*60);common.ses/max(common.ses)]';
%%
%Do the mapping and create a transitiness metric using common light curves
info=common;
nDim=3;
knn=15;
types=info.d;
%want=(types==2 | types==3 | types==4 | types==5 | types==7) & info.mes>7.8;
want=info.mes>8;
length(want(want))
savefile=['mappingQ1Q17DR24common_',num2str(datenum(date)-730000),'.mat'];

koilike=(types==2 | types==3 | types==4 | types==5 | types==7);

outfile=['lppknnDistanceCommon_',info.dettype,num2str(datenum(date)-730000),'.dat'];

[ outputInfo ] = compute_Transitiness( info,want,koilike,nDim,knn,outfile );

save(savefile,'outputInfo')
%%
%Set limit and run test
nsig=1.1;
decision=outputInfo.transitMetric > nsig* outputInfo.transit1sigmacut;
%outputInfo.d=infoTrap.d(wantTrap==1)
d=outputInfo.d';
print_tmetricResults(decision,d,'Test Results for 4 sigma cut Common Data')

%Create a plot of mapping and transitMetric
y=log10(outputInfo.transitMetric);
x=log10(outputInfo.mes);
interactive_plotting(x,y,outputInfo)
xlabel('log(MES)')
%figure
y=log10(outputInfo.transitMetric);
x=log10(outputInfo.periods);
interactive_plotting(x,y,outputInfo)
xlabel('log(Period)')

%Plot a normalized histogram
figure
outputInfo.dall=outputInfo.d;
outputInfo.d=outputInfo.d(~isinf(log10(outputInfo.transitMetric)));
outputInfo.logTM=log10(outputInfo.transitMetric(~isinf(log10(outputInfo.transitMetric))));
histogramParamLabels(outputInfo,'logTM',100)
title('Common LPP Metric')
xlabel('LPP Metric')
ylabel('Normalized Histogram')
print('-dpng','CommonMetricHist.png')
%Script to calculate the transitness for injected transits.
%Using an old mapping and population of known KOIs.
cd /home/smullall/Kepler/RoboVetter/transitInversion
%Step one. Create the median detrended light curves for each detrender.
other(1).dir='/soc/nfs/so-nfs/Q17-Inversion/TrapezoidFitPDC/DATS/';
other(2).dir='/soc/nfs/so-nfs/Q17-Inversion/DVSummaries/DATS/';
other(1).name='TrapDet'; %Do not change
other(2).name='DVMedDet'; %Do not change
other(1).tag='-fulltime.dat';
other(2).tag='-fulltime.dat';

storeDir='/home/sethomps/Kepler/RoboVetter/transitInversion/DataStore/';
tceList='TCEList.txt';  %Input File
%tceList='difftcelist.txt';
dvFileList='getDvPaths.txt';  %Run Fergal's code to get this.
%tceList='smalldetrtest.txt'
order=2;

[ output ] = createMedPDCLightCurves( dvFileList, tceList, storeDir, order, other, 1 );

%%
%Step two: Create a data matrix for each detrender.
tceFile='TCEList.txt';
storeDir='/home/sethomps/Kepler/RoboVetter/transitInversion/DataStore/';
npts=100;
ntrfr=5;
%over=10;

type='TrapDet';
[infoTrap] = createMatrixByType(tceFile,storeDir,npts,ntrfr,type);

%type='PdcMedDet';
%[infoPdc] = createMatrixByType(tceFile,storedir,npts,ntrfr,type);

type='DVMedDet';
[infoDV] = createMatrixByType(tceFile,storeDir,npts,ntrfr,type);


outfile=['inputInvQ1Q17DR24_Matrix',num2str(datenum(date)-730000),'.mat'];
save(outfile,'infoTrap','infoDV')

%%
%Step three: Gather the right information to apply the same mapping and
%comparision
%
%
origMapFile='../Q1Q17knnlpp/ProductionRun/mappingQ1Q17DR24dv_6083.mat';
load(origMapFile);
origMap=outputInfo;
%These lines can be removed for future runs.
%types=origMap.d;
%origMap.knnGood=(types==2 | types==3 | types==4 | types==5 | types==7);

%load(outfile)
info=infoDV;

savefile=['InvQ1Q17_LPP',info.dettype,num2str(datenum(date)-730000),'.dat'];
InvInfoDV=computeTransitinessKnownMap(info,origMap, savefile);

outfile=['InvQ1Q17_LPP',info.dettype,num2str(datenum(date)-730000),'.mat'];
save(outfile,'InvInfoDV')
%%
%Step three: Gather the right information to apply the same mapping and
%comparision
%
%
origMapFile='../Q1Q17knnlpp/ProductionRun/mappingQ1Q17DR24trap_6084.mat';
load(origMapFile);
origMap=outputInfo;
%These lines can be removed for future runs.
%types=origMap.d;
%origMap.knnGood=(types==2 | types==3 | types==4 | types==5 | types==7);

%load(outfile)
info=infoTrap;

savefile=['InvQ1Q17_LPP',info.dettype,num2str(datenum(date)-730000),'.dat'];
InvInfoTrap=computeTransitinessKnownMap(info,origMap, savefile);

outfile=['InvQ1Q17_LPP',info.dettype,num2str(datenum(date)-730000),'.mat'];
save(outfile,'InvInfoTrap')
%%
figure
plot(log10(InvInfoDV.periods),log10(InvInfoDV.transitMetric),'.k')
hold on;plot([-1,3],[-2.5,-2.5],'r-')
text(-0.25,0,'Fail LPP')
text(-0.25,-3.3,'Pass LPP')
xlabel('log(Period)')
title('Transit Inversion, LPP DV Detrender')
print('-dpng','dvInvQ1Q17period.png')

figure
plot(log10(InvInfoTrap.periods),log10(InvInfoTrap.transitMetric),'.k')
hold on;plot([-.5,3],[-2.5,-2.5],'r-')
text(-0.25,0,'Fail LPP')
text(-0.25,-3.3,'Pass LPP')
xlabel('log(Period)')
title('Transit Inversion, LPP Alt Detrender')
print('-dpng','altInvQ1Q17period.png')

%%
% %Apply the code to create the common lpp metric, integrating both results.
% %
% load('TrInQ1Q17DR24_LPPTrapDet6084.mat');
% outputInfoTrap=trInInfoTrap;
% load('TrInQ1Q17DR24_LPPDVMedDet6084.mat');
% outputInfoDV=trInInfoDV;
% %---
% wantTrap=zeros(length(outputInfoTrap.tce),1);
% wantDV=zeros(length(outputInfoDV.tce),1);
% 
% for i=1:length(outputInfoTrap.tce)
%     atce=outputInfoTrap.tce{i};
%     [trapi]=find(strcmp(atce,outputInfoDV.tce));
%     if length(trapi)>=1
%         wantTrap(i)=1;
%         wantDV(trapi)=1;
%     end
% end
% commontces=outputInfoTrap.tce(wantTrap==1);
% 
% common={};
% common.tce=commontces;
% common.DVmetric=outputInfoDV.transitMetric(wantDV==1);
% common.Trapmetric=outputInfoTrap.transitMetric(wantTrap==1);
% common.periods=outputInfoDV.periods(wantDV==1);
% common.mes=outputInfoDV.mes(wantDV==1);
% common.d=outputInfoDV.d(wantDV==1);
% common.dur=outputInfoDV.durs(wantDV==1);
% common.XTrap=outputInfoTrap.X(wantTrap==1,:);
% common.XDV=outputInfoDV.X(wantDV==1,:);
% common.XaDV=outputInfoDV.Xa(wantDV==1,:);
% common.XbDV=outputInfoDV.Xb(wantDV==1,:);
% common.wsmes=outputInfoTrap.wsmes(wantTrap==1);
% common.ses=outputInfoTrap.ses(wantTrap==1);
% common.dettype='common';
% common.svnVersion=outputInfoDV.svnVersion;
% common.X=[common.DVmetric;common.Trapmetric;common.wsmes/(max(common.wsmes));common.mes/max(common.mes);common.dur./(common.periods*24*60);common.ses/max(common.ses)]';
% %%
% %Do the mapping and create a transitiness metric using common light curves
% info=common;
% nDim=3;
% knn=15;
% types=info.d;
% %load known map
% previnfo=load('mappingQ1Q17DR24common_6084.mat');
% origMap=previnfo.outputInfo;
% %want=(types==2 | types==3 | types==4 | types==5 | types==7) & info.mes>7.8;
% want=info.mes>8;
% length(want(want))
% savefile=['TrInmappingQ1Q17DR24common_',num2str(datenum(date)-730000),'.mat'];
% 
% koilike=(types==2 | types==3 | types==4 | types==5 | types==7);
% 
% outfile=['TrInQ1Q17DR24_LPP',info.dettype,num2str(datenum(date)-730000),'.dat'];
% 
% trInInfoCommon=computeTransitinessKnownMap(info,origMap, savefile);
% 
% save(savefile,'trInInfoCommon')
% %%
% %This code is for updating the types after the fact. good for plotting.
% %
% tcefile='TrInTcesFed.txt';
% [tce,d]=textread(tcefile,'%s %f');
% 
% for i=1:length(trInInfoDV.d)
%     t=trInInfoDV.tce{i};
%     a=find(strcmp(t,tce));
%     trInInfoDV.d(i)=d(a);
% end
% %%
% for i=1:length(trInInfoTrap.d)
%     t=trInInfoTrap.tce{i};
%     a=find(strcmp(t,tce));
%     trInInfoTrap.d(i)=d(a);
% end
% 
% for i=1:length(trInInfoCommon.d)
%     t=trInInfoCommon.tce{i};
%     a=find(strcmp(t,tce));
%     trInInfoCommon.d(i)=d(a);
% end

%%
%Plot up some results
info=InvInfoDV;
x=log10(info.periods);
y=log10(info.transitMetric);
interactive_plotting(x,y,info)
title('DV Detrending')
%%
info=InvInfoTrap;
x=log10(info.periods);
y=log10(info.transitMetric);
interactive_plotting(x,y,info)
title('Trap Detrending')
% info=trInInfoCommon;
% x=log10(info.periods);
% y=log10(info.transitMetric);
% interactive_plotting(x,y,info)
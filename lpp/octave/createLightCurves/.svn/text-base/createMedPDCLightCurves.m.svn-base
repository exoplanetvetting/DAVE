function [ TCE ] = createMedPDCLightCurves( dvFileList, tceList, storeDir, order, otherdir, start )
%createLightCurves creates one .mat file per objects in tceList containing
%all the light curves you need for the LPP/KNN analysis on TCEs.
%Input 3 strings:
%dvFileList is created by running Fergal's getDvPaths
%tceList is a file of tces. It will also read in the second column
%giving an integer type (2=planet or what you might like)
%order is the order to remove the flux offsets in the quarterly pdc data
%first.
%other is a struct array contining the directories of other people's light
%curves. Usually the TPS median detrended and Chris Detrended.
%other needs a dir, name and tag for each entry. All strings.
%start is an integer to give the starting point from your tceList, when
%things crash early. usually this should be set to 1

[KIC,dvdir]=textread(dvFileList, '%f %s');
[TCE,type]=textread(tceList, '%s %f');

%
for i=start:length(TCE)
    clear output;
    clear added;
    
    t=TCE(i)
    kic=str2double(t{1}(1:9));
    pn=str2double(t{1}(11:12));
    
    j=find(KIC == kic);
    
  try
    [ dvOutput ] = getDvStructInfo(dvdir{j},kic,pn);
    [ dvInput ] = getDvInputLightCurves( dvdir{j});
    
    %Determine the Filter Length
    filterLen1=10*floor(dvOutput.pulse*3600/(1765.0));
    filterLen2=floor(dvOutput.period*24*3600/(1765.0));
    if filterLen1 > filterLen2
        filterLen=filterLen2;
    else
        filterLen=filterLen1;
    end
    
    %Do median detrending of the pdc light curve
    gaps=dvInput.pdcGaps;
    filled=dvInput.pdcFilled;
    if ~isempty(filled)
        gaps(filled+1)=1;  %Appears to be zero based indicies in filled
    end
    pdcFlux=dvInput.pdcTimeSeries(~gaps);
    quarters=dvInput.quarters(~gaps);
        
    [ medDet, polyDet ] = medDetPdcData( pdcFlux,filterLen,quarters, order);

    %Notice that there is no outlier rejection. only poly and median
    %detrending
    time=dvOutput.bjdTime(~gaps);
    added=struct([]);
    added(1).PdcMedDetTime=time;
    added.PdcMedDet=medDet;
    added.PdcPolyDet=polyDet;
    added.PdcAllGaps=gaps;
    added.medfilterLen=filterLen;
    

    %Get the other information from otherdir
    for k=1:length(otherdir)
        
        [time,flux] = getlcData(otherdir(k).dir,t,otherdir(k).tag);

        fname=otherdir(k).name;
        tname=[otherdir(k).name,'Time'];
        added.(fname)=flux;
        added.(tname)=time;
        
    end

    output=catstruct(dvOutput,dvInput,added);
    
    %A tsting plot
    %figure
    %subplot(1,1,1)
    %plot(added.TrapezoidDetTime,added.TrapezoidDet,'bo')
    %hold on
    %plot(added.TpsMedDetTime,added.TpsMedDet,'r.')
    %plot(added.PdcMedDetTime,added.PdcMedDet,'gs')
    
    
    %Write to the appropriate file
    file=[storeDir,t{1},'_lcinfo.mat']
    save(file,'output')

  catch exception
     'Failed to create file for'
     TCE{i}
     exception.identifier
     file=[storeDir,t{1},'_lcinfo.mat']
     save(file,'t')
  end
    
    
end


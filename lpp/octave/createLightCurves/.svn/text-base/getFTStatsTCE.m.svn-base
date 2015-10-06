function [ ] =getFTStatsTCE( tceFile,storedir,N )
%Add FTStatsPeriod list of harmonics to the output that is read in from the
%lc_info file in DataStore.
%N is number of period harmonics to fit.

[TCE,type]=textread(tceFile, '%s %f');

info=struct([]);

for j=1:length(TCE)

    t=TCE(j)
    kic=str2double(t{1}(1:9));
    pn=str2double(t{1}(11:12));

    fullpath=[storedir,t{1},'_lcinfo.mat'];
    
    if exist(fullpath,'file') == 2
        
        %Get information from the matfile
        load(fullpath);
        
         output.harmonicN=N;
         t=output.PdcMedDetTime;
         y=output.PdcMedDet;
         
         if length(t)>100
            [a,noiseLevel]=getFTStatsPeriod(t,y,output.period,N);
            output.harmonicAmpsPdc=a;
            output.harmonicNoisePdc=noiseLevel;
            output.harmonicSigPdc=a/noiseLevel;
         else
            output.harmonicAmps=[0,0,0,0];
            output.harmonicNoise=0;
            output.harmonicSig=[0,0,0,0];
            output.harmonicSigChange=[0,0,0,0];            
         end
         
         t=output.TpsMedDetTime;
         y=output.TpsMedDet;
         if length(t)>100
            [aTps,noiseLevelTps]=getFTStatsPeriod(t,y,output.period,N);
            output.harmonicAmpsTps=aTps;
            output.harmonicNoiseTps=noiseLevelTps;
            output.harmonicSigTps=a/noiseLevelTps;
            %output.harmonicSigChange=output.harmonicSigPdc-output.harmonicSigTps;
         else
            output.harmonicAmpsTps=[0,0,0,0];
            output.harmonicNoiseTps=0;
            output.harmonicSigTps=[0,0,0,0];
            %output.harmonicSigChange2=[0,0,0,0];
         end
         
         t=output.TrapezoidDetTime;
         y=output.TrapezoidDet;
         if length(t)>100
           [aTrap,noiseLevelTrap]=getFTStatsPeriod(t,y,output.period,N);
           output.harmonicAmpsTrap=aTrap;
           output.harmonicNoiseTrap=noiseLevelTrap;
           output.harmonicSigTrap=a/noiseLevelTrap;
           %output.harmonicSigChange2=output.harmonicSigTrap-output.harmonicSigTps;
         else
            output.harmonicAmpsTrap=[0,0,0,0];
            output.harmonicNoiseTrap=0;
            output.harmonicSigTrap=[0,0,0,0];
            %output.harmonicSigChange=[0,0,0,0];
            %output.harmonicSigChange2=[0,0,0,0];
         end
        
         
         save(fullpath,'output');
        
    end
end

end


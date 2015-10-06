function [info] = createMatrixByType(tceFile,storedir,npts,ntrfr,dettype)
%Function to create a sample matrix to feed to LLE or LPP
%Return a detrended and whitened output of a folded binned light curve.
%get inputs from the ddir location. 
%Second column of tceFile is an integer type.
%dettype is the name that is used in the lcinfo file
%PdcMedDet or TrapezoidDet for instance.
%ntrfr is the number of transit fraction

%over=10;
clear X
clear Xwhole
clear storeKIC
clear storePul
clear storeDur


%[datafiles]=textread([ddir,'lcinfo'], '%s');
[TCE,type]=textread(tceFile, '%s %f');

info=struct([]);
info(1).dettype=dettype;
n=1;
for j=1:length(TCE)

    t=TCE(j)
    kic=str2double(t{1}(1:9));
    pn=str2double(t{1}(11:12));

    fullpath=[storedir,t{1},'_lcinfo.mat'];
    
    if exist(fullpath,'file') == 2
        
        %Get information from the matfile
        load(fullpath);
        
        period=output.period;
        phase=output.phase;
        pulse=output.pulse;
        dur=output.dur;
        
        %Time and flux from the median detrended TPS light curve
        fname=dettype;
        tname=[dettype,'Time'];
        
        time=output.(tname);
        flux=output.(fname);
        

        if (mean(flux) ~=0 && ~isnan(mean(flux)) && length(time)>50)

          %Create phased light curve
          phaselc=mod((time-(phase-0.5*period))/period,1);
          
          %establish the fraction of the phased light curve that has the
          %transit, can't be more than half the period.
          if ~isnan(dur) && dur > 0
              trdur=dur; %Use measured duration if we have it.
          else
              trdur=pulse;
          end
          
          trfr=(trdur/24)/period;
          if ntrfr*trfr > 0.5;
              trfr=0.5/ntrfr;
          end
 
           %Run median average on two sets of bins. 
           a=[.03:1/npts:0.47,0.53:1/npts:0.97];
           b=( (0.5-ntrfr*trfr):(4*ntrfr*trfr)/(npts):(0.5+ntrfr*trfr) );
        
           [runta,runya]=running_median(phaselc,flux,1.25/(npts),a);  %1.25 was 1
           [runtb,runyb]=running_median(phaselc,flux,(3*ntrfr*trfr)/(npts),b); %3 was 2

           %Combine the two sets iof bins
           runymess=[runya;runyb];
           runtmess=[runta,runtb];
           [bins,ibins]=sort(runtmess);
           runy=runymess(ibins);
           runt=bins;

            %Scale by the minimum value of each light curve
            %mflux=median(runyb);
            scale=-1*min(runyb);
            if scale ~= 0
                scaledFlux=(runy)/scale;
                scaledA=runya/scale;
                scaledB=runyb/scale;
            else
                scaledFlux=runy;
                scaledA=runya;
                scaledB=runyb;
            end

           
%Comment out the plot for speed
%                 subplot(2,1,1)
%                  plot(phaselc,flux,'.')
%                  hold on
%                  plot(runt,runy,'gx--')
%                  hold off
%                 subplot(2,1,2)
%                   plot(scaledFlux,'o-')
%                 drawnow


            %shortscaledFlux=shortScaledFlux-mean(shortScaledFlux);
            %Add the trimmed flux to the data structre for the lle function
            info(1).X(n,:)=scaledFlux;
            info(1).Xa(n,:)=scaledA;
            info(1).Xb(n,:)=scaledB;
            info(1).kics(n)=kic;
            info(1).periods(n)=output.period;
            info(1).snrs(n)=output.snr;
            info(1).pulse(n)=output.pulse;
            info(1).durs(n)=output.dur;
            info(1).d(n)=type(j);
            info(1).ses(n)=output.maxses;
            info(1).mes(n)=output.maxmes;
            info(1).wsmes(n)=output.wsmaxmes;
            info(1).wsmin(n)=output.wsminmes;
            info(1).wsmaxphase(n)=output.wsmaxphase;
            info(1).wsminphase(n)=output.wsminphase;
            info(1).logg(n)=output.logg.value;
            info(1).effTemp(n)=output.effTemp.value;
            info(1).tce(n)=t;
            info(1).svnVersion='$Id: createMatrixByType.m 59152 2015-04-27 17:29:22Z smullall $';
            n=n+1;
  
        end
    end
end


%X(isnan(X))=0;
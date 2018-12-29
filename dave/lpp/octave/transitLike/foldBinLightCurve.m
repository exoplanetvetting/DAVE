function [ binnedFlux ] = foldBinLightCurve( dataStruct,ntrfr,npts)
%Fold and bin light curve for input to LPP metric calculation.
%   Detailed explanation goes here
%   duration in hours, period in days, time in days phase in bkjd
%

    time=dataStruct.time;
    phase=dataStruct.phase;
    dur=dataStruct.duration;
    period=dataStruct.period;
    flux=dataStruct.flux - 1;

          %Create phased light curve
          phaselc=mod((time-(phase-0.5*period))/period,1);
          
          %establish the fraction of the phased light curve that has the
          %transit, can't be more than half the period.
          if ~isnan(dur) && dur > 0
              trdur=dur; %Use measured duration if we have it.
          else
              trdur=.02*period/(24);
          end
          
          trfr=(trdur/24)/period;
          if ntrfr*trfr > 0.5;
              trfr=0.5/ntrfr;
          end
 
           %Run median average on two sets of bins. 
           a=[.03:1/npts:0.47,0.53:1/npts:0.97];
           b=( (0.5-ntrfr*trfr):(4*ntrfr*trfr)/(npts):(0.5+ntrfr*trfr) );
        
           [runta,runya]=running_median(phaselc,flux,1.5/(npts),a);  %1.25 was 1
           [runtb,runyb]=running_median(phaselc,flux,(5*ntrfr*trfr)/(npts),b); %3 was 2

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
                %scaledA=runya/scale;
                %scaledB=runyb/scale;
            else
                scaledFlux=runy;
                %scaledA=runya;
                %scaledB=runyb;
            end

        binnedFlux=scaledFlux;

end


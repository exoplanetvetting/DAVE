function [ info ] = foldBinLPPMatrix( time,flux,outputSt,ntrfr,npts )
%Does the same thing as the guts of createMatrixByType
%This function was created to do some plotting when not doing a production
%run

period=outputSt.period;
phase=outputSt.phase;
pulse=outputSt.pulse;
dur=outputSt.dur;

info=struct([]);

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
    
info(1).scaledFlux=scaledFlux;
info(1).scaledA=scaledA;
info(1).scaledB=scaledB;
info(1).phases=phaselc;
info(1).bins=bins;
info(1).abins=runta;
info(1).bbins=runtb;
info(1).aFlux=runya;
info(1).bFlux=runyb;
info(1).binFlux=runy;
info(1).scaleFactor=scale;
info(1).kic=outputSt.kic;
info(1).period=outputSt.period;
info(1).snr=outputSt.snr;
info(1).pulse=outputSt.pulse;
info(1).durs=outputSt.dur;


end


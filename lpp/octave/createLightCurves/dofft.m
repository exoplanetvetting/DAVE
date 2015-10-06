function [ freq,amp] = dofft( t,y,oversample )
%dofft takes the fft by respacing the unequally spaced data.
%It uses the first two points to determine the spacing of the data.

%resample data.
dtime=t(2)-t(1);
newt=t(1):dtime:t(end);
newd=pchip(t,y,newt);

NFFT=2^nextpow2(length(y)*oversample);

ampfft=fft(newd,NFFT);

amp=2*abs(ampfft(1:NFFT/2+1))/length(y);

freq =1/(2*dtime)*linspace(0,1,NFFT/2+1);

% Plot single-sided amplitude spectrum.
%plot(freq,amp) 
%drawnow
end


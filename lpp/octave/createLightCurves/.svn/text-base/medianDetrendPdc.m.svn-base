function [time, fluxlc, gaps]=medianDetrend(flux,times,gaps,filterLen)
%median detrend by filterLen cadences for a light curve.

%Check for outlier sigma limits.
upper=5.5;
lower=-30;

%Check that you have a reasonable filter length
if filterLen<0
    filterLen=96;
end

%smooth and normalize by requested filter length
want=~gaps & ~isnan(flux);
useflux=flux(want);
usetimes=times(want);

smoothfunc=medfilt1(useflux,filterLen);  %96-2days
smoothdata=useflux./smoothfunc - 1;  


outlier=smoothdata>upper*std(smoothdata) | smoothdata<lower*std(smoothdata);

fluxlc=useflux(~outlier);


        %concatenate these results onto the full vector;
        fluxlc=[fluxlc; smoothPdcData];
        time=[time; qtime{1}(~isnan(pdcdata{1}))];
        qinfo=[qinfo; qv];
        pdcraw=[pdcraw;pdcdata];
       
        gapflags=[gapflags;bitor(outlier,qual)];

        
function [time, detlc, gaps]=medianDetrend(flux,times,gaps,filterLen)
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
useflux=flux;
time=times;

smoothfunc=medfilt1(useflux,filterLen);  %96=2days
smoothdata=useflux./smoothfunc - 1;  

outlier=smoothdata>upper*std(smoothdata(want)) | smoothdata<lower*std(smoothdata(want));

gaps=~outlier & want;

detlc=smoothdata;

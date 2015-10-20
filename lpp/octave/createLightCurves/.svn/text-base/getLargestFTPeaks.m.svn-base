function [ f,a,noiseLevel ] = getLargestFTPeaks(t,y,npeaks )
%Take a light curve and return the n largest peak frequency and amplitudes.
%Need to prewhiten between. Would be quicker without prewhitening 
%

f=zeros(npeaks,1);
a=zeros(npeaks,1);

yresid=y;

for i=1:npeaks
    [ftf,fta]=dofft(t,yresid,8);

    [a(i),ind]=max(fta);
    f(i)=ftf(ind);
    v=[a(i),0];
    P=1/f(i);
    N=1;
    [ c_end,resid,c_amp,c_tzero,c_amp_err,c_tzero_err,c_per,yfit ] = fit_nsine_Pfixed( t,y,v,P,N );
    
    a(i)=c_amp(1);
    
    yresid=yresid-yfit;
    
    plot(ftf,fta,'r')
    drawnow
end

noiseLevel=(mean(fta.^2))^(1/2);


end


function [ a,noiseLevel ] = getFTStatsPeriod(t,y,period,N )
%Take a light curve and return the amplitdue at a particular period plus n
%harmonics.
%Also return the noise level overall.
%

a=0;
[ftf,fta]=dofft(t,y,8);
noiseLevel=(mean(fta.^2))^(1/2);

 v=[noiseLevel,0];
 P=period;
 
 [ c_end,resid,c_amp,c_tzero,c_amp_err,c_tzero_err,c_per,yfit ] = fit_nsine_Pfixed( t,y,v,P,N );
    
 a=c_amp;
  
 %figure
 %plot(ftf,fta,'r')
 %drawnow

end


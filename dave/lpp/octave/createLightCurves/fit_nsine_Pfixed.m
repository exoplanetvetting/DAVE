function [ c_end,resid,c_amp,c_tzero,c_amp_err,c_tzero_err,c_per,fity ] = fit_nsine_Pfixed( t,y,v,P,N )
%Fit N sine waves equally spaced in frequency.  for input frequency of 1/P
% Period is fixed to the input period.
%v contains the initial conditions of 1 frequency, 1 amplitude and 1 phase.
% frequency will be the equally spacing.
%This will return an array of N amplitudes and N phases for each harmonic

%size(t);
%size(y);

niter=100;
options=statset('MaxIter',niter);

%,'Robust','off','TolX',1e-9);

%Can fit one at a time since we are fixing the period
n=1;
c_amp=zeros(N,1);
c_amp_err=zeros(N,1);
c_tzero=zeros(N,1);
c_tzero_err=zeros(N,1);
c_per=zeros(N,1);
coef=zeros(2*N,1);
for i=1:2:N*2
    Pn= P/n;
    [cend,resid,J]=nlinfit(t,y,@sine_fun,[v(1),v(2)],options);
    ci = nlparci(cend,resid,'jacobian',J);
    coef(i)=cend(1); 
    coef(i+1)=cend(2);
    if cend(1) < 0
      c_tzero(n)=mod(cend(2)-Pn/2,Pn);
    else
      c_tzero(n)=mod(cend(2),Pn);
    end
    c_amp(n)=abs(cend(1));
    c_amp_err(n)=abs(ci(1,1)-ci(1,2))/4; 
    c_tzero_err(n)=abs(ci(2,1)-ci(2,2))/4;
    c_per(n)=Pn;
    
    fity=sine_fun([c_amp(n),c_tzero(n)],t);
    %y=y-fity;
    %plot(y)
    n=n+1;    
end
c_end=coef;


function [y]=sine_fun(v,t)

    y=v(1)*sin(2*pi*(1/Pn)*(t-v(2)));
end

end
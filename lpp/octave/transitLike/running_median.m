function [runt,runy]=running_median(t,y,dt,runt)
%Take a running meidan over chunks of dt.
%runt is an array of the points to return.
% dt is the size that will be combined and runt determines the points to
% return.

%
newy=zeros(length(y),1);
newt=zeros(length(y),1);

% set up output time.  Over resolve by 10.
%over=npts*dt/(max(t)-min(t));
%runt=min(t):dt/over:max(t);

%sort and get out indicies.
[st,stin]=sort(t);
%Create sorted arrays.
for i=1:length(y)
    newy(i)=y(stin(i));
    newt(i)=t(stin(i));
end

runy=zeros(length(runt),1);
for i=1:length(runt)
    tmp=(0);
    sum=0;
    c=1;
    for j=1:length(newt)
        if (newt(j) >=(runt(i)-dt) && newt(j) <= (runt(i)+dt))
            tmp(c)=newy(j);
            sum=sum+newy(j);
            c=c+1;
        end
    end
    %runy(i)=sum/(c-1);
    runy(i)=median(tmp);
end



%end


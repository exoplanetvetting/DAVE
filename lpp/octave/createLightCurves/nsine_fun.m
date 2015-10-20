function [y]=nsine_fun(v,t)
    N=v(1);
    f=v(2);
    y=zeros(size(t));
    for i=3:2:N*2+2
        y=y+v(i)*sin(2*pi*(f*(i-1)/2)*(t-v(i+1)));
    end
   

end
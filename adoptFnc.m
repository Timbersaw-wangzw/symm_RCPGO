function [f,w]=adoptFnc(x,a,c)
% hold on
% grid on
% axis equal
% fplot(@(x)adoptFnc(x,-1e5,1),[-5,5]);
% fplot(@(x)adoptFnc(x,-2,1),[-5,5]);
% fplot(@(x)adoptFnc(x,0,1),[-5,5]);
% fplot(@(x)adoptFnc(x,0.5,1),[-5,5]);
% fplot(@(x)adoptFnc(x,1,1),[-5,5]);
% fplot(@(x)adoptFnc(x,1.5,1),[-5,5]);
% fplot(@(x)adoptFnc(x,2,1),[-5,5]);
    r=x/c;
    r2=r.^2;
    N=length(x);
    w=zeros(N,1);
    f=w;
    for i=1:N
        if a==2
            f(i)=0.5*r2(i);
            w(i)=1/c^2;
        elseif a==0
            f(i)=log(0.5*r2(i)+1);
            w(i)=2/(x(i)^2+2*c^2);
        elseif a<-1e8
            f(i)=1-exp(-0.5*r2(i));
            w(i)=exp(-0.5*r2(i))/c^2;
        else
            a1=abs(a-2)/a;
            f(i)=a1*((r2(i)/abs(a-2)+1)^(a/2)-1);
            w(i)=(r2(i)/abs(a-2)+1)^(a/2-1);
        end
    end
end
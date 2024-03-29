function b= traceDiff(X,A)

b=zeros(6,1);

h=1e-8;



Xd=X.double();
for i=1:6
    I=zeros(6,1);
    I(i)=1;
    skew_I=skewa(I);
    b(i)=trace(Xd'*(skew_I'+skew_I)*Xd)-2*trace(A'*skew_I*Xd);

%     xi=h*I;
%     XX=SE3.exp(xi)*X;
%     j1=trace((XX-A)'*(XX-A));
%     j2=trace((X-A)'*(X-A));
%     z(i)=(j1-j2)/h;

end
end


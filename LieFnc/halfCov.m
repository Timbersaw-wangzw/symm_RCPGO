function half_cov=halfCov(pts)
N=length(pts(1,:));
M=zeros(6,6);
for i=1:N
    pi=pts(1:3,i);
    mi=dotVec(pi);
    M=M+mi'*mi;
end
half_cov=chol(M);
end
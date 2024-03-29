function [T_group,result,con_d]=multiViewICP(G,coeff,max_icp,disType,robType)

anchor_idx=G.numnodes;
edge_idx=G.Edges(:,1).EndNodes;
T_group=cell(G.numnodes,1);
%% initalization 
for i=1:G.numnodes
    T_group{i}=SE3;
end


x=zeros(6*(G.numnodes-1),1);    
res=zeros(G.numnodes-1,1);
for e=1:G.numnodes-1
    truth_T=G.Nodes{e,"T"}{1}.inv();
    res(e,1)=norm(truth_T-T_group{e},'fro');
end

con_d=[];
% res=[];
u12=1;

ww1=0.5;
ww2=1-ww1;
w1=ones(G.numedges,1);
w2=ones(G.numedges,1);
if contains(disType,'symm')
    [A1,b1,F01,L01]=getWJacobianAndResRevi(G,T_group,anchor_idx,w1,x,'ptp');
    [A2,b2,F02,L02]=getWJacobianAndResRevi(G,T_group,anchor_idx,w2,x,'qtq');
    A=ww1*A1+ww2*A2;
    b=ww1*b1+ww2*b2;
    F0=(F01+F02)/2;
    L0=(L01+L02)/2;
else
    [A,b,F0,L0]=getWJacobianAndResRevi(G,T_group,anchor_idx,w1,x,disType);
end
mu=coeff*max(diag(A));
% mu=100;
v=2;
A=A+mu*eye(6*(G.numnodes-1));
update_group=cell(G.numnodes-1,1,1);
% L0=F0;
% u1=2;
% u2=1e-2;
a_vec=[1/2,1/4, 0, -1/4, -1/2, -1, -2, -4, -8, -16, -32,-64,-128,-256];
a=a_vec(1);
for icp=1:max_icp
    opts.SYM = true;
    opts.POSDEF = true;
    x=linsolve(A,-b,opts);
    max_x=max(x);
    str=sprintf('iteration:%d-%d,max variation:%.3e',icp,max_icp,max_x);
    str=[disType,'-',robType,'-',str];
    disp(str);
    if max_x<1e-5
        break;
    end
    %% Update
    res_vec=zeros(G.numnodes-1,1);
    for i=1:G.numnodes-1
        xi=x(6*(i-1)+1:6*i);
        dTi=SE3.exp(xi);
        update_group{i}=dTi;
        T_group{i}=dTi*T_group{i};
        truth_T=G.Nodes{i,"T"}{1}.inv();
        res_vec(i)=norm(truth_T-T_group{i},'fro');
    end
    res=[res,res_vec];

    r1=zeros(G.numedges,1);
    r2=zeros(G.numedges,1);
    Fh=0;
    for e=1:G.numedges
        idx_i=edge_idx(e,1);
        idx_j=edge_idx(e,2);
        Ti=T_group{idx_i};
        Tj=T_group{idx_j};
        Tij=G.Edges{e,'Tij'};
        if strcmp(disType,'ptp')
            half_cov=G.Edges{e,"half_cov_p"}{1};
            Te=Tij.inv()*Tj.inv()*Ti;
            be=half_cov*vexa(Te.log());
            r1(e)=be'*be;
        end
        if strcmp(disType,'qtq')
            Te=Tj.inv()*Ti*Tij.inv();
            half_cov=G.Edges{e,"half_cov_q"}{1};
            be=half_cov*vexa(Te.log());
            r1(e)=be'*be;
            Fh=Fh+w1(e)*norm(r1(e))^2;
        end
        Fh=Fh+w1(e)*norm(r1(e))^2;
        if contains(disType,'symm')
            Fh=Fh+w2(e)*norm(r2(e))^2;
            half_cov=G.Edges{e,"half_cov_p"}{1};
            Te=Tij.inv()*Tj.inv()*Ti;
            be=half_cov*vexa(Te.log());
            r1(e)=be'*be;

            Te=Tj.inv()*Ti*Tij.inv();
            half_cov=G.Edges{e,"half_cov_q"}{1};
            be=half_cov*vexa(Te.log());
            r2(e)=be'*be;
 
        end
    end
    if icp==1
        u21 = 1e-2;
        u11 = 4*median(r1);
        if contains(disType,'symm')
            u22 = 1e-2;
            u12 = 4*median(r2);
        end
        mu=1e4;
    end
    if contains(disType,'symm')
        [~,~,Fh1,Lh1]=getWJacobianAndResRevi(G,T_group,anchor_idx,w1,x,'ptp');
        [~,~,Fh2,Lh2]=getWJacobianAndResRevi(G,T_group,anchor_idx,w2,x,'qtq');
        Fh=Fh1+Fh2;
        Lh=Lh1+Lh2;
    else
        [~,~,Fh,Lh]=getWJacobianAndResRevi(G,T_group,anchor_idx,w1,x,disType);
    end
    %% update weights
    if contains(disType,'symm')
        if contains(robType,'welsch')
            w1=weight(r1,u11);
            w2=weight(r2,u12);
        else
            [~,w1]=adoptFnc(r1,a,u11);
            [~,w2]=adoptFnc(r2,a,u12);
        end
    else
        if contains(robType,'welsch')
            w1=weight(r1,u11);
        else
            [~,w1]=adoptFnc(r1,a,u11);
        end
    end
    rho=(F0-Fh)/(L0-Lh);
    %% gets new jacobians and residuals
    if contains(disType,'symm')
        [A1,b1,F01,~]=getWJacobianAndResRevi(G,T_group,anchor_idx,w1,zeros(6*(G.numnodes-1),1),'ptp');
        [A2,b2,F02,~]=getWJacobianAndResRevi(G,T_group,anchor_idx,w2,zeros(6*(G.numnodes-1),1),'qtq');
        A=ww1*A1+ww2*A2;
        b=ww1*b1+ww2*b2;
        F0=(F01+F02)/2;
    else
        [A,b,F0,~]=getWJacobianAndResRevi(G,T_group,anchor_idx,w1,zeros(6*(G.numnodes-1),1),'ptp');
    end
    L0=F0;
    if mu<1e7
        if rho>0
            mu=mu*max(1/3,1-(2*rho-1)^3);
            v=2;
        else
            mu=mu*2;
            v=2*v;
        end
    end
    idx=4;
    if(mod(icp,4)==0)
        u11 = u11/2;
        u12 = u12/2;
    end
    
    if(mod(icp,idx)==0)

        idx=floor(icp/idx);
        if idx>length(a_vec)
            a = a_vec(end)/idx;
        else
            a = a_vec(idx);
        end
    end
    if u11<u21
        u11 = u21;
    end
    if contains(disType,'symm')
        if u12<u22
            u12 = u22;
        end
    end
    A=A+1e4*eye(6*(G.numnodes-1));
%     A=A+mu*eye(6*(G.numnodes-1));
end
method=[robType,'_',disType];
result.method=method;
result.res=res;
end
function [f,g]=quadric_f(x,A,b)
    f=0.5*x'*A*x+b'*x;
    g=A*x+b;
end
function [c,ceq,gradc,gradceq]=quadric_con(x,A_group,b_group,bd_vec)
len=length(bd_vec(:,1));
c=zeros(2*len,1);
gradceq=[];
ceq=[];
% gradc=zeros(1,len);
for e=1:len
    Ae=A_group{e};
    be=b_group{e};
    fe=Ae*x+be;
    lb=bd_vec(e,1);
    ub=bd_vec(e,2);
%     lb=0;
%     ub=300;
    c(2*(e-1)+1)=fe'*fe-(ub)^2-1e-3;
    c(2*e)=(lb)^2-fe'*fe-1e-3;
    ge=2*Ae'*fe;
    gradc(:,2*(e-1)+1)=ge;
    gradc(:,2*e)=-ge;
end
end
function Hout = hessianfcn(x,lambda,A0,A_group)
    H=A0;
    len=length(A_group);
    Hout=H;
    for e=1:len
        Ae=A_group{e};
        Hg1=2*(Ae'*Ae);
        Hg2=-2*(Ae'*Ae);
        Hout = Hout + lambda.ineqnonlin(2*(e-1)+1)*Hg1+lambda.ineqnonlin(2*e)*Hg2;
    end
end
%% Get constriant graph 
function w = weight(r,u)
w=zeros(1,length(r));
for i=1:length(r)
    w(i) = exp(-r(i)/(2*u^2));
end
end
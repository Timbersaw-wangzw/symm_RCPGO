function [T_group,result,con_d]=conwMultiViewICP(G,CG,max_icp,tol,disType,robType,info,lambda)

% last pose is anchor point
anchor_idx=G.numnodes;
edge_idx=G.Edges(:,1).EndNodes;
T_group=cell(G.numnodes,1);
%% initalization 
for i=1:G.numnodes
    T_group{i}=SE3;
end



x=zeros(6*(G.numnodes-1),1);    

res=zeros(G.numnodes-1,1);
con_d=zeros(CG.numedges,1);
initial_d=zeros(CG.numedges,1);
ratio=zeros(CG.numedges,1);
for ce=1:CG.numedges
    O1=CG.Edges(ce,"Oi").Oi;
    O2=CG.Edges(ce,"Oj").Oj;
    d=CG.Edges(ce,4).d;
    initial_d(ce)=norm(O1-O2);
    max_con_iter=floor(max_icp/4);
    % ratio of dynamic slack of constraint
    ratio(ce)=abs(d-1e-7-initial_d(ce))/max_con_iter;
    con_d(ce,1)=abs(initial_d(ce)-d);
end


for e=1:G.numnodes-1
    truth_T=G.Nodes{e,"T"}{1}.inv();
    res(e,1)=norm(truth_T-T_group{e},'fro');
end

con_d=[];
% res=[];

u12=1;
w1=ones(G.numedges,1);
w2=ones(G.numedges,1);
if contains(disType,'symm')
    [A1,b1,F01,L01]=getWJacobianAndResRevi(G,T_group,anchor_idx,w1,x,'ptp');
    [A2,b2,F02,L02]=getWJacobianAndResRevi(G,T_group,anchor_idx,w2,x,'qtq');
    A=A1+A2;
    b=b1+b2;
    F0=F01+F02;
    L0=L01+L02;
else
    [A,b,F0,L0]=getWJacobianAndResRevi(G,T_group,anchor_idx,w1,x,disType);
end

[A_group,b_group,bd_vec,C]=getConstraints(CG,anchor_idx,initial_d,ratio,1,tol);
mu=1e-6*max(diag(A));
% mu=100;
v=2;
A=A+mu*eye(6*(G.numnodes-1));
update_group=cell(G.numnodes-1,1,1);
F0=F0+100*C;
L0=L0+100*C;
% u1=2;
% u2=1e-2;
a_vec=[1/2,1/4, 0, -1/4, -1/2, -1, -2, -4, -8, -16, -32,-64,-128,-256];
a=a_vec(1);
for icp=1:max_icp
    opts.SYM = true;
    opts.POSDEF = true;
%     x0=linsolve(A,-b,opts);
    x0=zeros(6*(G.numnodes-1),1);
%     fprintf("iteration:%d-%d\n",icp,max_icp) 
%     fun=@(x)quadric_f(x,A,b,A_group,b_group);
    fun=@(x)quadric_f(x,A,b);
    nonlcon=@(x)quadric_con(x,A_group,b_group,bd_vec);
    options = optimoptions('fmincon','Algorithm','interior-point',...
    'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true, ...
    'HessianFcn',@(x,lambda) hessianfcn(x,lambda,A,A_group), ...
    'CheckGradients',false,...
    'Display','none');
    [x,Lh] = fmincon(fun,x0,[],[],[],[],[],[],nonlcon,options);
    max_x=max(x);
    fprintf("con-welsch-iteration:%d-%d,max variation:%.3e\n",icp,max_icp,max_x);
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
    %% Update constraints
    con_d_vec=zeros(CG.numedges,1);
    Ch=0;
    Dlh=0;
    for ce=1:CG.numedges
        idx_i=CG.Edges.EndNodes(ce,1);
        idx_j=CG.Edges.EndNodes(ce,2);
        if idx_i~=anchor_idx
            dTi=update_group{idx_i};
        else
            dTi=SE3;
        end
        if idx_j~=anchor_idx
            dTj=update_group{idx_j};
            
        else
            dTj=SE3;
        end
        Oi=dTi*CG.Edges(ce,"Oi").Oi';
        Oj=dTj*CG.Edges(ce,"Oj").Oj';
        CG.Edges(ce,"Oi").Oi=Oi';
        CG.Edges(ce,"Oj").Oj=Oj';
        d=CG.Edges(ce,4).d;
        con_d_vec(ce)=abs(norm(Oi-Oj)-d);
        Ch=Ch+con_d_vec(ce);

        Ae=A_group{ce};
        be=b_group{ce};
        fe=Ae*x+be;
        Dlh=Dlh+abs(norm(fe)-d);
    end
    con_d=[con_d,con_d_vec];
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
        u11 = 5*median(r1);
        if contains(disType,'symm')
            u22 = 1e-2;
            u12 = 5*median(r2);
        end
        mu=1e4;
    end
    if contains(disType,'symm')
        [~,~,Fh1,Lh1]=getWJacobianAndResRevi(G,T_group,anchor_idx,w1,x,'ptp');
        [~,~,Fh2,Lh2]=getWJacobianAndResRevi(G,T_group,anchor_idx,w2,x,'qtq');
        Fh=Fh1+Fh2;
        Lh=Lh1+Lh2;
    else
        [~,~,Fh,Lh]=getWJacobianAndResRevi(G,T_group,anchor_idx,w1,x,'ptp');
    end
    Fh=Fh+100*Ch;
    Lh=Lh+100*Dlh;
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

    [A_group,b_group,bd_vec,C]=getConstraints(CG,anchor_idx,initial_d,ratio,icp,tol);
    if contains(disType,'symm')
        [A1,b1,F11,~]=getWJacobianAndResRevi(G,T_group,anchor_idx,w1,x,'ptp');
        [A2,b2,F12,~]=getWJacobianAndResRevi(G,T_group,anchor_idx,w2,x,'qtq');
        A=A1+A2;
        b=b1+b2;
        F1=F11+F12;
    else
        [A,b,F1,Lh]=getWJacobianAndResRevi(G,T_group,anchor_idx,w1,x,'ptp');
    end
    rho=(F0-Fh)/(L0-Lh);
    F0=F1+100*C;
    L0=100*C;
    if mu<1e7
        if rho>0&&Fh<F0
            mu=mu*max(1/3,1-(2*rho-1)^3);
            v=2;
        else
            mu=mu*v;
            v=2*v;
        end
        if(mod(icp,4)==0)
            u11 = u11/2;
            u12 = u12/2;
        end
        idx=2;
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
    end
%     A=A+mu*eye(6*(G.numnodes-1));
    A=A+lambda*eye(6*(G.numnodes-1));
end
% info= sprintf('_RCPGO_Num%d',CG.numedges);

method=[disType,info];
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
    c(2*(e-1)+1)=fe'*fe-(ub)^2-1e-6;
    c(2*e)=(lb)^2-fe'*fe-1e-6;
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
function [A_group,b_group,bd_vec,C]=getConstraints(CG,anchor_idx,initial_d,ratio,icp,tol)
Edges=CG.Edges;
edge_idx=Edges(:,1).EndNodes;
A_group=cell(CG.numedges,1);
b_group=cell(CG.numedges,1);
bd_vec=zeros(CG.numedges,2);
threshold=0.001;
C=0;
for e=1:CG.numedges
    Ae=zeros(4,6*(CG.numnodes-1));
    idx_i=edge_idx(e,1);
    idx_j=edge_idx(e,2);
    O1=CG.Edges(e,"Oi").Oi';
    O2=CG.Edges(e,"Oj").Oj';
    Mi=dotVec(O1);
    Mj=dotVec(O2);
    d=CG.Edges(e,4).d;
    if idx_i~=anchor_idx&&idx_j~=anchor_idx
        Ae(1:4,6*(idx_i-1)+1:6*idx_i)=Mi;
        Ae(1:4,6*(idx_j-1)+1:6*idx_j)=-Mj;
    else
        Ae(1:4,6*(idx_i-1)+1:6*idx_i)=Mi;
    end
    new_d=norm(O1-O2);
    C=C+abs(new_d-d);
    if new_d<d-tol% get closer to the lower bound and constant upper bound
        bd_vec(e,2)=d;
        if abs(norm(O1-O2)-d)<threshold
            lb=d-tol;
        else
            lb_d=initial_d(e)+icp*ratio(e);
            if(lb_d>d||new_d>d)
                lb=d-tol;
            else
                lb=max([new_d,lb_d]);
            end
        end
        bd_vec(e,1)=lb;
    else% get closer to the upper bound and lower bound is constant
        bd_vec(e,1)=d;
        if abs(norm(O1-O2)-d)<threshold
            ub=d+tol;
        else
            ub_d=initial_d(e)-icp*ratio(e);
            if(ub_d<d||new_d<d)
                ub=d+tol;
            else
                ub=min([new_d,ub_d]);
            end
        end
        bd_vec(e,2)=ub;
    end
    A_group{e}=Ae;
    b_group{e}=[O1-O2;0];
end

end

function w = weight(r,u)
w=zeros(1,length(r));
for i=1:length(r)
    w(i) = exp(-r(i)/(2*u^2));
end
end
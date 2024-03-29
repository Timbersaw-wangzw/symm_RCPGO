function [A,b,F0,L0]=getWJacobianAndResRevi(G,T_group,anchor_idx,w,x,type)
A=zeros(6*(G.numnodes-1),6*(G.numnodes-1));
b=zeros(6*(G.numnodes-1),1);
Edges=G.Edges;
edge_idx=Edges(:,1).EndNodes;
L0=0;
F0=0;
w = w/sum(w); 
w2 = sqrt(w);
if strcmp(type,'ptp')||strcmp(type,'qtq')
    for e=1:G.numedges
        Ae=zeros(6,6*(G.numnodes-1));
        idx_i=edge_idx(e,1);
        idx_j=edge_idx(e,2);
        Ti=T_group{idx_i};
        Tj=T_group{idx_j};
        Tij=Edges{e,'Tij'};
        [Mi,Mj,be]=disFnc(G,e,anchor_idx,Tij,Ti,Tj,x,type);
        if idx_i~=anchor_idx&&idx_j~=anchor_idx
            Ae(1:6,6*(idx_i-1)+1:6*idx_i)=Mi(1:6,1:6);
            Ae(1:6,6*(idx_j-1)+1:6*idx_j)=Mj(1:6,1:6);
        else
            Ae(1:6,6*(idx_i-1)+1:6*idx_i)=Mi(1:6,1:6);
        end
        if strcmp(type,'ptp')
            half_cov=G.Edges{e,"half_cov_p"}{1};
        else
            half_cov=G.Edges{e,"half_cov_q"}{1};
        end
        cov_Ae=w2(e)*half_cov*Ae;
        cov_be=w2(e)*half_cov*be;
        
        Le=(cov_Ae*x+cov_be);
        A=A+cov_Ae'*cov_Ae;
        b=b+cov_Ae'*cov_be;
        F0=F0+cov_be'*cov_be;
        L0=L0+Le'*Le;
    end
elseif strcmp(type,'pln')
    for e=1:G.numedges
        Ae=zeros(6,6*(G.numnodes-1));
        idx_i=edge_idx(e,1);
        idx_j=edge_idx(e,2);
        Ti=T_group{idx_i};
        Tj=T_group{idx_j};
        Tij=Edges{e,'Tij'};
        [Mi,Mj,be]=disFnc(G,e,anchor_idx,Tij,Ti,Tj,x,'ptp');
        if idx_i~=anchor_idx&&idx_j~=anchor_idx
            Ae(1:6,6*(idx_i-1)+1:6*idx_i)=Mi(1:6,1:6);
            Ae(1:6,6*(idx_j-1)+1:6*idx_j)=Mj(1:6,1:6);
        else
            Ae(1:6,6*(idx_i-1)+1:6*idx_i)=Mi(1:6,1:6);
        end
        ptsN_p=G.Edges{e,"Pi"}(1);
        ptsN_q=G.Edges{e,"Pj"}(1);
        
        ptsP=ptsN_p.pts(1:3,:);
        np=Ti.SO3*ptsN_p.pts(4:6,:);
        Tq=Tj*Tij;
%         Rij=Tj.SO3.double()*Tij.SO3.double();
        nq=Tq.SO3*np;
        nn=np+nq;
        Pi=Ti*ptsP;
        Pj=Tj*Tij*ptsP;
        N=length(ptsP(1,:));
        nn=[nn;ones(1,N)];
        le=Ae*x+be;
        for i=1:N
            Mi=nn(:,i)'*(Tq.double())*dotVec(Pi(:,i));
            Ai=w2(e)*Mi*Ae;
            A=A+Ai'*Ai;
            bi=w2(e)*Mi*be;
            b=b+Ai'*bi;
            F0=F0+bi'*bi;
            lei=w2(e)*Mi*le;
            L0=L0+lei'*lei;
        end
        F0=L0;
    end
end
end
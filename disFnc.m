function [Mi,Mj,be] = disFnc(G,e,anchor_idx,Tij,Ti,Tj,x,type)
switch type
    case 'ptp'
        Te=Tij.inv()*Tj.inv()*Ti;
        invTi=Ti.inv();
        AdinvTi=invTi.Ad();
        be=vexa(Te.log());
        invJr=InvJr(be);
        Mi=invJr*AdinvTi;
        Mj=-1*invJr*AdinvTi;
    case 'qtq'
        Te=Tj.inv()*Ti*Tij.inv();
        Tq=Ti*Tij.inv();
        invTq=Tq.inv();
        AdinvTq=invTq.Ad();
        be=vexa(Te.log());
        invJr=InvJr(be);
        Mi=invJr*AdinvTq;
        Mj=-1*invJr*AdinvTq;
    case 'pln'
        Edges=G.Edges;
        edge_idx=Edges(:,1).EndNodes;
        
%         Mi=zeros(6*(G.numnodes-1),6*(G.numnodes-1));
        idx_i=edge_idx(e,1);
        idx_j=edge_idx(e,2);
        be=zeros(6*(G.numnodes-1),1);

        
        
end
end


function con_d= updateCon(CG,T_group)
%% Update constraints
C=0;
con_d=zeros(CG.numedges,2);
for ce=1:CG.numedges
    idx_i=CG.Edges.EndNodes(ce,1);
    idx_j=CG.Edges.EndNodes(ce,2);
    Ti=T_group{idx_i};
    Tj=T_group{idx_j};
    d=CG.Edges(ce,4).d;
    Oi=CG.Edges(ce,"Oi").Oi;
    Oj=CG.Edges(ce,"Oj").Oj;
    con_d(ce,1)=abs(norm(Oi-Oj)-d);
    Oii=Ti*CG.Edges(ce,"Oi").Oi';
    Ojj=Tj*CG.Edges(ce,"Oj").Oj';
    con_d(ce,2)=abs(norm(Oii-Ojj)-d);
end
end


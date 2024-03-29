function plotGRMSE(G,CG,x,y)
%PLOTGRMSE 此处显示有关此函数的摘要
%   此处显示详细说明
for e=1:G.numedges
    ii=G.Edges(e,1).EndNodes(1);
    jj=G.Edges(e,1).EndNodes(2);
%     Ti=G.Nodes(ii,1).T{1};
%     Tj=G.Nodes(jj,1).T{1};
    Ti=G.Nodes(ii,1).T{1};
    Tj=G.Nodes(jj,1).T{1};
    TTij=Tj*inv(Ti);
    Tij=G.Edges{e,"Tij"}(1);
    rmse(e)=norm(TTij.double()-Tij.double(),'fro');
end
% figure
% plot(log10(rmse),'*-');
% xlim([1,G.numedges]);
% xticks([1,10,20,30,40,50,60,72]);

figure
axis off
axis equal

% p=plot(G,'XData',x,'YData',y);

p=plot(G,'XData',x,'YData',y);
p.EdgeCData=log10(rmse);
p.LineWidth=2;
p.EdgeFontSize=12;
colormap jet
caxis([-4,-2])
% colorbar southoutside
ax=gca;
ax.XAxis.Visible='off';
ax.YAxis.Visible='off';

for ce=1:CG.numedges
%     highlight(p,CG_idx(i,:),'EdgeColor','black','LineWidth',5,'LineStyle','-.','NodeColor','red')
    O1=CG.Edges(ce,2).Oi{1,1};
    O2=CG.Edges(ce,3).Oj{1,1};
    d=CG.Edges(ce,4).d;
    edge_rmse(ce)=abs(norm(O1-O2)-d);
end
figure
p=plot(CG,'XData',x,'YData',y);
p.LineWidth=2;
p.EdgeFontSize=12;
p.EdgeCData=edge_rmse;
ax=gca;
ax.XAxis.Visible='off';
ax.YAxis.Visible='off';
% CG_idx=CG.Edges.EndNodes;
colormap jet

% colorbar southoutside
end


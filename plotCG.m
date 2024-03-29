function plotCG(CG,x,y,edge_rmse)
% Ni=sqrt(G.numnodes);
% Nj=Ni;
% xx = 0:1:Ni-1;
% yy= Nj-1:-1:0;
% x=repmat(xx,1,Nj);
% y=repmat(yy,Ni,1);
% y=reshape(y,25,1);

if nargin==3
    for ce=1:CG.numedges
    %     highlight(p,CG_idx(i,:),'EdgeColor','black','LineWidth',5,'LineStyle','-.','NodeColor','red')
        O1=CG.Edges(ce,2).Oi{1,1};
        O2=CG.Edges(ce,3).Oj{1,1};
        d=CG.Edges(ce,4).d;
        edge_rmse(ce)=abs(norm(O1-O2)-d);
    end
end

figure
p=plot(CG,'XData',x,'YData',y);
p.LineWidth=2;
p.EdgeFontSize=12;
p.EdgeCData=log10(edge_rmse);
clim([-4,1]);
xlim([0,8.5])
ylim([0,8.5])
ax=gca;
ax.XAxis.Visible='off';
ax.YAxis.Visible='off';
% CG_idx=CG.Edges.EndNodes;
colormap jet
colorbar southoutside
fontname(gcf,"Times New Roman")
fontsize(gcf,8, "points")
set(gcf,'unit','centimeters','Position',[0 0,5,5])
% set(gca,'looseInset',[0 0 0.03 0.03]);
% set(gca,'LooseInset', get(gca,'TightInset'))
% set(gca,'position', [0.05 0.05 0.8 0.8]);
% set(gcf, 'PaperPosition', [0.5 0.5 5.5 5.5])
% set(gcf, 'PaperSize', [5.5 5.5])
% print('123.emf','-demf','-r600')
exportgraphics(gcf,'CG.emf','Resolution',600)
end


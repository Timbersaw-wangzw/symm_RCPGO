function [CG,info]= diffEdges(G,CG,CGType)
%DIFFEDGES 此处显示有关此函数的摘要
%   此处显示详细说明
s=[1 6 12 2 4 5 7];
t=[20 25 23 10 6 11 15];

switch CGType
    case 1
        x=1:72;
        CG=CG.rmedge(x);
        info='_RCPGO_CG1';
    case 2
        x=1:72;
        x(60)=[];
        x(10)=[];
        x(45)=[];
        x(20)=[];
        x(25)=[];
        x(22)=[];
        x(23)=[];
        x(24)=[];
        CG=CG.rmedge(x);
        info='_RCPGO_CG2';
        % edge_idx=floor(linspace(1,72,9));
        
        % x=1:72;
        % for i=1:length(edge_idx)
        %     x(edge_idx(i))=nan;
        % end
        % xx=x(~isnan(x));
        % CG=CG.rmedge(xx);
    case 3
        edge_idx=floor(linspace(1,72,20));
        info='_RCPGO_CG3';
        x=1:72;
        for i=1:length(edge_idx)
            x(edge_idx(i))=nan;
        end
        xx=x(~isnan(x));
        CG=CG.rmedge(xx);
    case 4
        edge_idx=floor(linspace(1,72,31));
        info='_RCPGO_CG4';
        x=1:72;
        for i=1:length(edge_idx)
            x(edge_idx(i))=nan;
        end
        xx=x(~isnan(x));
        CG=CG.rmedge(xx);
    case 5
        edge_idx=floor(linspace(1,72,42));
        info='_RCPGO_CG5';
        x=1:72;
        for i=1:length(edge_idx)
            x(edge_idx(i))=nan;
        end
        xx=x(~isnan(x));
        CG=CG.rmedge(xx);
    case 6
        x=1:4:72;
        CG=CG.rmedge(x);
        info='_RCPGO_CG6';
end

CG = addCG(G,CG,s,t);
% plotCG(CG,x,y)
% plot(G)
end

function CG= addCG(G,CG,s,t)
all_O=zeros(3,G.numnodes);
all_TO=zeros(3,G.numnodes);
for i=1:G.numnodes
    Ti=G.Nodes{i,"T"}{1};
    pts=G.Nodes{i,"pts"}{1};

    truth_pts=inv(Ti)*pts(1:3,:);
    mean_Oi=mean(pts(1:3,:),2);
    all_O(:,i)=mean_Oi;
    mean_TOi=mean(truth_pts,2);
    all_TO(:,i)=mean_TOi;
end
for i=1:length(s)
    ii=s(i);
    jj=t(i);
    O1=all_O(:,ii);
    O2=all_O(:,jj);
    TO1=all_TO(:,ii);
    TO2=all_TO(:,jj);
    d=norm(TO1-TO2);
    newEdge=table([ii jj], O1',O2',d,'VariableNames',{'EndNodes','Oi','Oj','d'});
    CG=addedge(CG,newEdge); 
end

end

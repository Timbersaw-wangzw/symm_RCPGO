% CG.Edges{1,"Oi"}
clc
clear all
close all
addpath('github_repo');
addpath('cylinderData\');
addpath('LieFnc\')
addpath('realData\')
addpath('welsch\')
experiments="real";
% experiments="simulation";
if contains(experiments,"simulation")
    %% simulated Graph
    load("sim_G.mat")
    load("sim_CG.mat")
    [CG,info]= diffEdges(G,CG,4);% change CG tolopogy with different edges
    lambda=1e4;% regularization parameter in simulation
    xx = 0:1:4;
    yy= 4:-1:0;
    x=repmat(xx,1,5);% coordinates used in plot CG 
    y=repmat(yy,5,1);% 
    y=reshape(y,25,1);% coordinates used in plot CG
    max_icp=100;
else
    %% real graph
    load("real_CG.mat");
    load("real_G.mat");
    % load("result_t_rcpgo.mat")
    % load("result_t_pgo.mat")
    %% lambda 
    lambda=5e1;% regularization parameter in real experiments
    y=8*[1,1,1,0,0,0,0,0,0,0,0.5,0,1];% coordinates used in plot CG
    x=[0,1,2,0,1,2,3,4,5,6,7,8,8];% coordinates used in plot CG
%     max_icp=80;
    max_icp=300;
    info="real_CG";
end
%% distance type
% disType='ptp';
% disType='qtq';
disType='symm';
%% robust function type
robType = 'welsch';
% robType = 'adopt';
%% PGO with or without CG
method= sprintf('con_CG_%d',CG.numedges);
% method='non_con';
% method='all'
coeff=1e-3;
%% symm_RCPGO
if contains(method,'non')
    [T_group,result]=multiViewICP(G,coeff,max_icp,disType,robType);
    con_d= updateCon(CG,T_group);
else
    [T_group,result,con_d]=conwMultiViewICP(G, ...
        CG,max_icp,coeff, ...
        disType,robType,info,5e1);
    
end
%% plot error in CG
plotCG(CG,x,y,con_d(:,end));
figure
hold on
axis off
for i=1:G.numnodes
    % truth_T=G.Nodes{e,"T"}{1}.inv();
    TT=T_group{i};
    pts=G.Nodes{i,"pts"}{1}(1:3,:);
    mv_pts=TT*pts;
%     mv_pts=pts;
    plot3(mv_pts(1,:),mv_pts(2,:),mv_pts(3,:),'.');
end

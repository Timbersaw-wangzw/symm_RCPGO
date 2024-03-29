# High Accuracy and Robust Robotic Inspection by Constrainted Pose Graph Optimization
This repository includes the source code of the paper High Accuracy and Robust Robotic Inspection by Constrainted Pose Graph Optimization.
This code was written in MATLAB. 
# Introduction and Usage
The directory `github_repo` is Lie algebra library;
The directory `cylinderData` contains pose graph 'sim_G.m' and constrainted pose graph 'sim_CG.m' in simulation;
The directory `LieFnc` is the right and left jacobian and other operators in Lie algebra which is need in optimization.
The directory `realData` contains pose graph 'real_G.m' and constrainted pose graph 'real_CG.m' in real experiments;
The directory `welsch` containts robust traditional PGO 'multiViewICP.m' and ours symm_RCGPO 'conwMultiViewICP.m'.
You can run directly `main.m` about simulation and real experiments.
There are some results shown below 
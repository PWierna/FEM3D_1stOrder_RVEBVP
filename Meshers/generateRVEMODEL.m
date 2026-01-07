% This script takes a raw RVE-mesh model (usually obtained after reading an
% ANSYS-generated mesh) and gives it the proper structure for the posterior
% RVE-Problem solver. It also generates the data associated to the node
% pairs necessary for the application of Periodic BCs. A prismatic RVE
% model is assumed.
clear,clc,close all

% Load Input file
load('RVE3L_28el_raw.mat');

% Output filename
filename = 'RVE3L_28elv2';%'RVE0_coarsex1'%;'rveansys_ref';%'rveansys_materials';%'microPorto01_fine';

% Add useful paths
addpath("D:\PWierna\MyCodes\FEMRVE3D_1stOrder_BVP\PreProcessing");

% Just make input consistent with solver input's structure
MODEL.Coordinates = [ (1:size(MODEL.Coordinates,1))' MODEL.Coordinates ];
MODEL.Conectivity = [ (1:size(MODEL.Conectivity,1))' MODEL.Conectivity ];

% Identify node pairs and add to model structure
RVEMODEL = identifyRVENodePairs( MODEL );

% Save final model file
save([filename,'.mat'],'RVEMODEL')

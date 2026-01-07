clear,clc
project_path    = "Example/";
model_fname     = "micro_model.mat";
materials_fname = "micro_materials.mat";

%% LOAD INPUT FILES ============================================================== %%
load(project_path+model_fname,'RVEMODEL'); % RVE-MODEL 
load(project_path+materials_fname,'MESO_MATERIALS'); % Materials data

%% INPUTS/SETTINGS =============================================================== %%

%Set general options:
Options.kinadm_conditions_type = "Periodic"; %Options: "Periodic" , "Linear" , "Taylor" , "Minimal"
Options.quadrature = "Gauss 3x3x3";
Options.Meso.Verbosity = 1;
Options.vtk_filename = 'vtk_output';
Options.insertion_type = "Strains";
%RVEMODEL.Conectivity(:,2)=1;%Force homogeneous

%Set macro-scale strains and displacements values to be inserted at the RVE:
StrainM = 1e-1*[ 0.71, 0.03, 0.28, 0.05, 0.10, 0.82 ]; %[xx,yy,zz,xy,xz,yz]
DisplM  = zeros(3, 1); %[Ux,Uy,Uz]


%% ADD PATHS ===================================================================== %%        
addpath(genpath('Routines/'));         %Add main routines path
addpath(genpath('ConstitutiveLaws/')); %Add ConstitutiveLaws' path
addpath(genpath('Elements/'));         %Add Elements Libraries path
addpath(genpath('Utilities/'));        %Add Utilities path
addpath(genpath('PostProcessing/'));   %Add Post-Processing scripts path
addpath(genpath('PreProcessing/'));    %Add Pre-Processing scripts path


%% INITIALIZATIONS =============================================================== %%

%Initialize fluctuating displacements vector and internal variables matrix
[d_0, IntVarsPrev] = InitializeRVEProblem(RVEMODEL, Options);

%Compute (linearized) Macro-Displacement Field:
d_M = MicroBVP3D_InsertMacroDisplacements_at_nodes(DisplM, StrainM, RVEMODEL);


%% SOLVE BVP ===================================================================== %%
[ Strains, FluctStrains, Stresses, IntVarsNew, d_New, MesoState ] = MicroRVE3D_HFLaw(d_M, StrainM, IntVarsPrev, d_0, RVEMODEL, MESO_MATERIALS, Options);


%% POST-PROCESSING =============================================================== %%
RVEProblem3D_VTKOutput(d_M ,StrainM, d_New, Strains, FluctStrains, Stresses, IntVarsNew.Damage, RVEMODEL, project_path, Options);


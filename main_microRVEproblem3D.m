project_path    = "Example/";
model_fname     = "micro_model.mat";
materials_fname = "micro_materials.mat";

%% LOAD FILES -------------------------------------------------------------------- %%
load(project_path+model_fname,'RVEMODEL');
load(project_path+materials_fname,'MESO_MATERIALS');

%% INPUTS ------------------------------------------------------------------------ %%
%Set options:
Options.kinadm_conditions_type = "Periodic"; %Options: "Periodic" , "Linear" , "Taylor" , "Minimal"
Options.quadrature = "Gauss 3x3x3";
Options.Meso.Verbosity = 1;
Options.vtk_filename = 'vtk_output';
Options.insertion_type = "Strains";

%RVEMODEL.Conectivity(:,2)=1;%Force homogeneous

%Set macro- strain and displacement values:
StrainM = 1e-1*[ 0.71 ;  %Exx
                 0.03 ;  %Eyy
                 0.28 ;  %Ezz
                 0.05 ;  %Exy
                 0.10 ;  %Exz
                 0.82 ]; %Eyz
DisplM = zeros( 3 , 1 );


%% ADD PATHS --------------------------------------------------------------------- %%        
addpath(genpath('Routines/'));         %Add main routines path
addpath(genpath('ConstitutiveLaws/')); %Add ConstitutiveLaws' path
addpath(genpath('Elements/'));         %Add Elements Libraries path
addpath(genpath('Utilities/'));        %Add Utilities path
addpath(genpath('PostProcessing/'));   %Add Post-Processing scripts path
addpath(genpath('PreProcessing/'));    %Add Pre-Processing scripts path


%% INITIALIZATIONS --------------------------------------------------------------- %%
dofpn = 3;	%dofs per node( Hard setting )
npe   = size(RVEMODEL.Conectivity(:,2:end),2); %Nodes per element
ndofs = dofpn * size(RVEMODEL.Coordinates, 1); %Total N°dofs

%-Get element type (we assume uniform mesh):
switch npe
    case 8 %8-noded bilinear element
        elem_type = "Hexa8N";
    case 27
        elem_type = "Hexa27N";
end

%-Get quadrature data:
[ n_lips , ~ , ~ ] = QuadratureData3D( Options.quadrature );

%-Get total number of global ips: 
n_gips  = n_lips * size(RVEMODEL.Conectivity,1); %(I assume all elements use the same quadrature)

%Initialize displacements and internal variables:
d_0 = zeros( ndofs , 1 );
IntVarsPrev = zeros( n_gips , 2 );

%Compute (linearized) Macro-Displacement Field:
d_M = MicroBVP3D_InsertMacroDisplacements_at_nodes( DisplM , StrainM , RVEMODEL );

%% SOLVE BVP --------------------------------------------------------------------- %%
[ Strains , FluctStrains , Stresses , IntVarsNew , d_New , MesoState ] = MicroRVE3D_HFLaw( d_M , StrainM , IntVarsPrev , d_0 , RVEMODEL , MESO_MATERIALS , Options );


%% POST-PROCESSING --------------------------------------------------------------- %%
RVEProblem3D_VTKOutput( d_M , StrainM , d_New , Strains , FluctStrains , Stresses , IntVarsNew.Damage , RVEMODEL , project_path , Options );


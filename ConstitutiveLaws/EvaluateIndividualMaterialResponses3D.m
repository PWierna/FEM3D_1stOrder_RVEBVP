function [ Stresses , TgModulis , IntVarsNew ] = EvaluateIndividualMaterialResponses3D( Strains , IntVarsHist , GIPsCon , MODEL , MATERIALS ) 


%% INITIALIZATIONS =============================================================== %%
% IP Arrays:
n_gips    = size(Strains,1);
Stresses  = zeros(n_gips,6);  %Vector with InPlane Stresses at Gauss Points 
TgModulis = zeros(n_gips,36); %Matrix with In-plane Ctg coefficients at IPs 
                              %Rows = GlobalIP , Cols = [ Axialx , Axialy , Shearxy ]
IntVarsNew.StrainLike = zeros(n_gips,2);   %Vector with r Internal Variable at IPs [r_fiber , r_matrix]
IntVarsNew.Damage     = zeros(n_gips,2);   %Vector with Damage Variable at IPs [d_fiber , d_matrix]

%% GET "TRUE" MATERIALS ========================================================== %%
true_materials_ids = unique( MATERIALS.LIST(MODEL.Conectivity(:,2)) ); %Unique materials

%% LOOP OVER MATERIALS =========================================================== %%
for mat_id = 1:length(true_materials_ids)    
    
    %Get indexes of the ips corresponding to this material:
    idxmat_gips = find( MATERIALS.LIST(MODEL.Conectivity(GIPsCon(:,1),2)) == true_materials_ids(mat_id) );
    
    %Get material properties:
    matprops    = MATERIALS.(true_materials_ids(mat_id));
    
    %Compute individual material response for every ip with this material:
    switch matprops.ConstitutiveLaw
        case "IsotContDamage3D"
            %TODO
        case "OrthotropicElastic3D"
            [ Stresses(idxmat_gips,:),...
              TgModulis(idxmat_gips,:) ...
              ] = OrthotropicElastic3D_V( Strains(idxmat_gips,:), matprops );
            %IntVarsNew.Damage(idxmat_gips,:)     = 0.0;
            %IntVarsNew.StrainLike(idxmat_gips,:) = 0.0;
        case "VoigtFibersMatrix_Damage"
            %TODO
        case "VoigtFibersMatrix_Damage_Uncoupled"
            %TODO
        case "IsotLinearElastic3D"
            %TODO
        case "InterfacesUniaxialDamage_Reg"
            %TODO
    end
    
end

return
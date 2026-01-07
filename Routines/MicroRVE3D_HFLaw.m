function [ Strains, FluctStrains, Stresses, IntVarsNew , d_New , MesoState ] = MicroRVE3D_HFLaw( DisplM , StrainM , IntVarsPrev , d_Old , RVEMODEL , MESO_MATERIALS , Options )
%% =============================================================================== %%
% - High-Fidelity Mesoscale Equilibrium Problem: 3DRVE for plates (4DoF elements) - % 
% ------------------------ Newton-Raphson Iterative Scheme ------------------------ %
%                                                                                   %
%                                                                                   %
%                                                                  PWierna  20-X-22 %
%% =============================================================================== %%

% HARD SETTINGS ------------------------------------------------------------------- %
res_tol    = 1e-03;
max_iter   = 15;
elem_macro = 1;
ip_macro   = 1;
dofpn      = 3;	%dofs per node( Hard setting )
% --------------------------------------------------------------------------------- %

%% 1st COMPUTATIONS ============================================================== %% 
n_dofs  = dofpn * size(RVEMODEL.Coordinates, 1); %Total N°dofs
StrainM = reshape(StrainM,6,1); %Assert macro-strain vector is in column fashion

%-Compute Boundary conditions Matrix:
[ L , dofd ] = build_RVE3D_BCmatrix( RVEMODEL , Options );   
dofs = 1:n_dofs;              %DoFs vector
doff = dofs; doff(dofd) = []; %Free DoFs vector

%% N-R NON-LINEAR ITERATIVE LOOP ================================================= %%
    
%Initializizations before iterative loop 
d_New     = d_Old; %Fluctuant displacements vector (solution)
iter      = 0;     %Initialize iteration counter
norm_res  = 1;     %Initialize residual norm value
norm_res0 = 1;     %Reference Residual Norm to evaluate convergence criteria

%Report RVE Equilibrium Problem ID:
if Options.Meso.Verbosity == 1
    %fprintf(ProcessData.OutputID,'\n\t\t MESO PROBLEM: Macro Element N°%d-IP %d\n',elem_macro,ip_macro);
    fprintf(1, '\n\t\t MESO PROBLEM: Macro Element N°%d-IP %d\n',elem_macro,ip_macro);
end


%Start iterative loop until convergence-------------------------------------------- %
while ( norm_res >= res_tol && iter <= max_iter ) 

    if Options.Meso.Verbosity == 1
        %fprintf(ProcessData.OutputID,'\t\t\t Meso [NEWTON-RAPHSON] Iteration %d:',iter);
        fprintf(1, '\t\t\t Meso [NEWTON-RAPHSON] Iteration %d:',iter);
    end
    
    %Evaluate tg stiffness and internal force vector:
    [ Ktg , Fint , Stresses , Strains , FluctStrains, IntVarsNew ] = Hexa_StiffnessFint( d_New , StrainM , DisplM , IntVarsPrev , RVEMODEL , MESO_MATERIALS , Options );

    %Compute Meso-Residual forces vector:
    Res = Fint; %In absence of external forces
    if iter == 0
        norm_res0 = norm( L'*Res(dofd) + Res(doff) );
        if norm_res0 <= 1e-6; norm_res0 = 1.0; end %Prevent 0/0 or NaN in case 1st Res = 0 (trivial problem)
    end
    
    %Compute Convergence Criteria:
    norm_res = norm( L'*Res(dofd) + Res(doff) ) / norm_res0; 
    
    %Report
    if Options.Meso.Verbosity == 1
        %fprintf(ProcessData.OutputID,'\t| res | = %4.6e\n', norm_res);
        fprintf(1, '\t| res | = %4.6e\n', norm_res);
    end
    
    %If converged, break while loop
    if (norm_res < res_tol)
        break
    end
    
    %Update iteration counter
    iter = iter + 1; 

    %Solve linearized equations system
    Ddf = - ( L'*Ktg(dofd,dofd)*L + L'*Ktg(dofd,doff) + Ktg(doff,dofd)*L + Ktg(doff,doff) ) \ ( L'*Res(dofd) + Res(doff) );
    Ddd = L * Ddf;

    %Update displacement fluctuations
    d_New(doff) = d_New(doff) + Ddf;
    d_New(dofd) = d_New(dofd) + Ddd;
        
end %<----------------------------------------------- <--- [Ends Iterative Loop] -- %
    

%Report convergence and update Meso-State:
if iter > max_iter
    MesoState = 0;
    if Options.Meso.Verbosity == 1
        %fprintf(ProcessData.OutputID, '\n\t [NEWTON-RAPHSON] WARNING - Max N° of iterations exceeded\n');
        %fprintf(ProcessData.OutputID, '\t                  Solver did not converge. ElemMacro%d IP%d \n',elem_macro,ip_macro);
        %fprintf(ProcessData.OutputID, '\t                  |res| = %1.4e > %1.4e = tol\n', norm_res,res_tol);
        fprintf(2, '\n\t [NEWTON-RAPHSON] WARNING - Max N° of iterations exceeded\n');
        fprintf(2, '\t                  Solver did not converge. ElemMacro%d IP%d \n',elem_macro,ip_macro);
        fprintf(2, '\t                  |res| = %1.4e > %1.4e = tol\n', norm_res,res_tol);
    end
else
    MesoState = 1;
    if Options.Meso.Verbosity == 1
        %fprintf(ProcessData.OutputID,'\t\t\t Meso Problem EL%d - IP%d Converged in [ %d / %d ] iterations .-\n', elem_macro, ip_macro, iter, max_iter);
        fprintf(1, '\t\t\t Meso Problem EL%d - IP%d Converged in [ %d / %d ] iterations .-\n', elem_macro, ip_macro, iter, max_iter);
    end
end


%% HOMOGENIZATION ================================================================ %%

%Get ips global coordinates and weights:
[ Wgips , ~ ] = Hexa_GlobalIPWeightsandCoords( RVEMODEL , Options );

%Compute integral of the fluctuating displacements over the rve:
int_strains_fluct = Wgips' * FluctStrains;

check_kinadmconds = boolean(zeros(1,6));
check_kinadmconds( abs(int_strains_fluct)<=1e-8 ) = true;

if all(check_kinadmconds)
    fprintf('\nKinematic admisibility conditions: OK\n')
else
    fprintf('\nKinematic adm conditions not satisfied. Check:\n')
    disp(string(check_kinadmconds))
    disp(int_strains_fluct)
end

% Reminder: Arrangement of Stresses and Constitutive properties:
    %           Sigmax = Stresses(:,1);
    %           Sigmay = Stresses(:,2);
    %           Sigmay = Stresses(:,3);
    %           Taoxy  = Stresses(:,4);
    %           Taoxz  = Stresses(:,5);
    %           Taoyz  = Stresses(:,6);
    
% %-Get quadrature data:
% [ n_lips , w_lips , c_lips ] = QuadratureData3D( Options.quadrature );
% 
% %-Global IP Conectivity : (WE ASSUME ALL ELEMENTS USE THE SAME QUADRATURE)
% %n_gips  = n_lips * n_elem; %Total number of global ips 
% GIPsCon = [repelem((1:n_elem)',n_lips),repmat((1:n_lips)',n_elem,1)]; %[ Element N° , Local IP N° ]
% 
% %-Compute Local Interpolation data:
% switch elem_type
%     case "Hexa8N"
%         LocalInterpData = LagrangeInterp3D_C08N_Local( c_lips );
%     case "Hexa27N"
%         [ LocalInterpData , ~ ] = LagrangeInterp3D_C127N_Local( c_lips );
% end

% %-Compute Global Interpolating Data:
% [ N_gips , detJ_gips , ~ , ~ , dNdX_gips ] = LagrangeInterp3D_Global_Vect( GIPsCon , LocalInterpData , RVEMODEL ); %Think this should be inside WeightsandCoords
% 
% %-Get Global IP Weights and Coordinates:
% [ GIPWeights , GIPCoords ] = RVE1D_GlobalWeightsAndCoords( gip_con , LocalInterpData(:,1,:) , detJ_Gips , w_lips , RVEMODEL ); 

return
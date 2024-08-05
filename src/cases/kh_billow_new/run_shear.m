
%% driver_mode1_shoal.m
% Driver script for mode1 shoal case.
% Parameter file 'mode1_shoal.txt' contains list of cases to run
% with different parameters. This script makes the corresponding
% spins.conf files and the initial u,v,w,rho fields for the model.

clearvars

test = true;                       % set to false to write data to disk

% Case independent parameters
% Spatial parameters
Lx = 10.54075393;
Ly = 1.0;
Lz = 0.3;
Nx = 1024;
Ny = 1;
Nz = 256;
min_x = 0;
min_y = 0;
min_z = 0;

% Expansion types
type_x = 'FOURIER';
type_y = 'FREE_SLIP';
type_z = 'FREE_SLIP';
mapped_grid = 'false';
% Physical parameters
g = 9.81;
rot_f = 0.0;
rho_0 = 1000.0;
visco = 1.e-6;
kappa_rho = 1.e-7;

pert = 1.e-3;

dz_u = Lz/10;
dz_rho = dz_u;
a = 0;
R = 5;
rho_perturb = 0.0e-8;
perturb_k = 2.38434;

% Temporal Parameters
final_time = 10;
plot_interval = .5;
dt_max = 0.0;

% Restart Options
restart = 'false';
restart_time = 0.0;
restart_sequence = 0;
restart_from_dump = 'false';
compute_time = -1;

% Perturbation Parameter
perturb = 0.5e-3;

% Filter Parameters
f_cutoff = 0.6;
f_order = 2.0;
f_strength = 20.0;

%% These are specific to this experiment

% Recreate profiles from Carpenter et al.,
% Constants
delRho = 0.13;
Nj = 7; % Number of J points
Nlambda = 5; % Number of Wavelengths to test

% Variables
J = linspace(-.251 , 1, Nj);
J = J(~J==0);

% Assume del rho fixed:
delRho = delRho.*sign(J);
g_prime = delRho*g/rho_0;
delU = sqrt(g_prime*dz_rho./J)';
%delta_u = 0.3616628264
%delta_rho = 0.01;

% Secondary diagnostics
compute_stresses_bottom = true;

cdir = pwd; % Remember directory script is run from
%% Loop through different cases
for numcase = 1:Nj
    casename{numcase} = ['KH_Case_', num2str(numcase)];
    if numcase == 1
        mkdir(casename{1});
        cd(casename{1});
    end    
    % write data to disk
    mkdir(['../' casename{numcase}]), cd(['../' casename{numcase}])
    
    %% Write params to spins.conf. Change what is written according to your experiment.
    fid = fopen('spins.conf','wt');
    fprintf(fid, '## KH Configuration File \n');
    fprintf(fid, 'name = %s', casename{numcase});
    fprintf(fid, '\n # Spatial Parameters \n');
    fprintf(fid,'Lx = %6.2f \n', Lx);
    fprintf(fid,'Ly = %6.2f \n', Ly);
    fprintf(fid,'Lz = %6.2f \n', Lz);
    fprintf(fid,'Nx = %d \n', Nx);
    fprintf(fid,'Ny = %d \n', Ny);
    fprintf(fid,'Nz = %d \n', Nz);
    
    fprintf(fid,'min_x = %6.2f \n', min_x);
    fprintf(fid,'min_y = %6.2f \n', min_y);
    fprintf(fid,'min_z = %6.2f \n', min_z);
    
    fprintf(fid, '\n # Expansion types \n');
    fprintf(fid,'type_x = %s \n', type_x);
    fprintf(fid,'type_y = %s \n', type_y);
    fprintf(fid,'type_z = %s \n', type_z);
    fprintf(fid,'mapped_grid = %s \n', mapped_grid);
    
    fprintf(fid, '\n # Physical Parameters \n');
    fprintf(fid,'g = %12.3f \n', g);
    fprintf(fid,'rho_0 = %12.1f \n', rho_0);
    fprintf(fid,'visco = %.2e \n', visco);
    fprintf(fid,'kappa_rho = %.2e \n', kappa_rho);
    
    fprintf(fid, '\n # Problem Parameters \n');
    fprintf(fid,'delta_rho = %4.3f \n', delRho(numcase));
    fprintf(fid, 'delta_u = %4.3f \n', delU(numcase));
    fprintf(fid, 'dz_u = %4.3f \n', dz_u);
    fprintf(fid, 'dz_rho = %4.3f \n', dz_rho);
    fprintf(fid, 'a = %4.3f \n', a);
    fprintf(fid, 'R = %4.3f \n', R);
    fprintf(fid, 'rho_perturb = %4.3f \n', rho_perturb);
    fprintf(fid, 'perturb_k = %4.3f \n', perturb_k);

    fprintf(fid, '\n # Temporal Parameters \n');
    fprintf(fid,'final_time = %12.8f\n',final_time);
    fprintf(fid,'plot_interval = %12.8f \n',plot_interval);
    
    fprintf(fid, '\n # Restart Parameters \n');
    fprintf(fid,'restart = false \n');
    fprintf(fid,'restart_time = 0.0 \n');
    fprintf(fid,'restart_sequence=0 \n');
    fprintf(fid,'restart_from_dump = false \n');
    fprintf(fid,'compute_time = -1 \n');
    
    fprintf(fid, '\n # Filter Parameters \n');
    fprintf(fid,'f_cutoff = %12.8f \n', f_cutoff);
    fprintf(fid,'f_order = %12.8f \n', f_order);
    fprintf(fid,'f_strength = %12.8f \n', f_strength);
    
    fprintf(fid, '\n # Diagnostics \n');
    if compute_stresses_bottom
        fprintf(fid, 'compute_stresses_bottom = false');
    end
    
    fclose(fid);
    
    % Copy executable files across to directory
    if isunix
        copyfile('../kh_billow_new/kh_billow_new.x', '.');
    end 
end

fid = fopen('../case_list','wt');
for i = 1:Nj
    fprintf(fid, [casename{i}, '\n']);
end




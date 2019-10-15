function load_HYDROTHERM_output(folder)

if nargin==0
%     folder = 'F:\Projects\HT2D\Benchmarks\Hydrotherm HTI Data\Sim_1';
    folder = '../examples/Ex.08';
    suffix='.test';
end

% ------------------------ fluid properties ------------------------

% pressure
[P,istep,time,varnames,x,y,z,nx,ny,nz] = read_HYDROTHERM_file([folder '/Out_pressure' suffix]);
VAR.P = P{1}';

% enthalpy
[H,istep,time,varnames] = read_HYDROTHERM_file([folder '/Out_enthalpy' suffix]);
VAR.H = H{1}';

% temperature
[T,istep,time,varnames] = read_HYDROTHERM_file([folder '/Out_temperature' suffix]);
VAR.T = T{1}';

% density
[Rho,istep,time,varnames] = read_HYDROTHERM_file([folder '/Out_density' suffix]);
VAR.Rho_l = 1e3*Rho{1}'; % conversion from g/cm^3 to SI units
VAR.Rho_v = 1e3*Rho{2}'; % conversion from g/cm^3 to SI units

% viscosity
[Mu,istep,time,varnames] = read_HYDROTHERM_file([folder '/Out_viscosity' suffix]);
VAR.Mu_l = 1e-1*Mu{1}'; % conversion from g/cm-s ("Poise" ) to SI units
VAR.Mu_v = 1e-1*Mu{2}'; % conversion from g/cm-s ("Poise" ) to SI units

% velocity
[V,istep,time,varnames] = read_HYDROTHERM_file([folder '/Out_velocity' suffix]);
VAR.Vx_l = 0.01*V{1}'; % conversion from cm/s to SI units
VAR.Vz_l = 0.01*V{2}'; % conversion from cm/s to SI units
VAR.Vx_v = 0.01*V{3}'; % conversion from cm/s to SI units
VAR.Vz_v = 0.01*V{4}'; % conversion from cm/s to SI units

% ------------------------ matrix properties ------------------------
% porosity
[Props,varname] = read_HYDROTHERM_Probdefine_file([folder '/Out_Probdefine'  suffix],size(VAR.T,2));
VAR.Phi   = Props{1};
VAR.Perm  = Props{2};
VAR.K_m   = Props{3};
VAR.Cp_m  = Props{4};
VAR.Rho_m = Props{5};

% ------------------------ save as matlab data ------------------------
save([folder '/HTI_data.mat'],'istep','time','x','y','z','nx','ny','nz','VAR');

end
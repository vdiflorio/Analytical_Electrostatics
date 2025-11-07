  
% Copyright (C) 2024-2025 Sergii V. Siryk
% Copyright (C) 2024-2025 Vincenzo Di Florio
% 
% This program is free software: you can redistribute it and/or modify  
% it under the terms of the GNU General Public License as published by  
% the Free Software Foundation, either version 3 of the License, or  
% (at your option) any later version.  
% 
% This program is distributed in the hope that it will be useful,  
% but WITHOUT ANY WARRANTY; without even the implied warranty of  
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the  
% GNU General Public License for more details.  
% 
% You should have received a copy of the GNU General Public License  
% along with this program. If not, see <https://www.gnu.org/licenses/>.  
 
 
function analytical_pb(pqr_path, mode)
% analytical_pb.m
%
% This script computes electrostatic energy and, optionally, the
% electrostatic potential for a system of charged spheres defined in a PQR file.
%
% -------------------------------------------------------------------------
% USAGE EXAMPLES:
%
%   % 1. Default usage (reads '../data/30spheres.pqr', computes energy only)
%   matlab -nodisplay -r "compute_potential"
%
%   % 2. Specify PQR path, compute energy only
%   matlab -nodisplay -r "compute_potential('../data/30spheres.pqr', 'energy')"
%
%   % 3. Compute potential on given points from a file (x y z columns)
%   matlab -nodisplay -r "compute_potential('../data/30spheres.pqr', 'points.txt')"
%
%   % 4. Generate a mesh box around the spheres and compute potential map
%   matlab -nodisplay -r "compute_potential('../data/30spheres.pqr', 'mesh')"
%
% -------------------------------------------------------------------------
% INPUT ARGUMENTS:
%   pqr_path : (optional) Path to PQR file (default: '../data/30spheres.pqr')
%   mode     : (optional) One of:
%       'energy'   - only compute energies (default)
%       'points.txt' - compute potential at points in the file
%       'mesh'     - build bounding box around spheres and compute potential map
%
% -------------------------------------------------------------------------
% OUTPUT:
%   - Prints real parts of total, reaction, ionic, and Coulombic energies.
%   - Optionally computes potential if requested.
% -------------------------------------------------------------------------


%% ==== Handle input arguments and defaults ====
if nargin < 1 || isempty(pqr_path)
    pqr_path = '../data/30spheres.pqr';
end
if nargin < 2 || isempty(mode)
    mode = 'energy'; % default behavior
end

fprintf('>>> Using PQR file: %s\n', pqr_path);
fprintf('>>> Mode: %s\n', mode);

%% ==== Physical and medium parameters ====
medium_params.epsilon = 80;

e_0  = 8.85418781762e-12; 
kb   = 1.380649e-23;    
T    = 273.15 + 25; 
e    = 1.602176634e-19;
N_av = 6.022e23;    
Angs = 1e-10;       
e_m  = 2.0;
e_w  = 80.;
ionic_strength = 0.145;
C_0 = 1.0e3 * N_av * ionic_strength;

eps_in  = 4.0*pi * e_0 * e_m * kb*T*Angs/(e^2); 
eps_out = 4.0*pi * e_0 * e_w * kb*T*Angs/(e^2); 
k = sqrt(2.0 * C_0 * Angs^2 * e^2 / (e_0 * e_w * kb * T));
t_eps0_inv = (e^2)/(4.0*pi * e_0 * kb*T*Angs);
medium_params.kappa = k;

%% ==== Read spheres from PQR file ====
fprintf('>>> Reading PQR file...\n');
fid = fopen(pqr_path);
if fid < 0
    error('Cannot open PQR file: %s', pqr_path);
end
C = textscan(fid, '%*s %*s %*s %*s %*s %f %f %f %f %f', ...
    'Delimiter', {' '}, 'MultipleDelimsAsOne', true, 'CollectOutput', true);
fclose(fid);
pqr = C{1};

num_particles = size(pqr, 1);
fprintf('>>> Found %d particles.\n', num_particles);

for i = 1:num_particles
    particle.center = [pqr(i,1), pqr(i,2), pqr(i,3)];
    particle.charge = pqr(i,4);
    particle.radius = pqr(i,5);
    particle.dielectric_constant = 2;
    particles_params(i) = particle;
end

%% ==== Compute multipole coefficients and total energy ====
fprintf('>>> Loading Clebsch-Gordan coefficients...\n');
load array_cg15.mat array_clebschgordan cg_n_max;

n_max = 15;
tol_gmres = 1e-12;

fprintf('>>> Computing multipole coefficients...\n');
tic;
[particles_coefficients, energy] = multi_spheres_cg5( ...
    particles_params, medium_params, n_max, ...
    array_clebschgordan, cg_n_max, tol_gmres, t_eps0_inv);
toc;

fprintf('>>> Computing energy components...\n');
[energy_reaction, energy_ionic, energy_Coulombic] = ...
    calc_energy_components0(particles_params, medium_params, n_max, ...
    particles_coefficients, energy, t_eps0_inv);

%% ==== Print energy components ====
fprintf('\n=== Energy summary ===\n');
fprintf('Total energy:      %e\n', real(energy));
fprintf('Reaction energy:   %e\n', real(energy_reaction));
fprintf('Ionic energy:      %e\n', real(energy_ionic));
fprintf('Coulombic energy:  %e\n', real(energy_Coulombic));
fprintf('=====================================\n\n');

%% ==== Optionally compute potential ====
if strcmpi(mode, 'energy')
    fprintf('>>> Energy-only mode: skipping potential calculation.\n');
else
    if strcmpi(mode, 'mesh')
        fprintf('>>> Generating uniform grid around spheres...\n');

        % Determine bounding box for all spheres
        centers = reshape([particles_params.center], 3, []).';
        radii = [particles_params.radius]';
        min_corner = min(centers - radii, [], 1);
        max_corner = max(centers + radii, [], 1);
        padding = 5.0;
        min_corner = min_corner - padding;
        max_corner = max_corner + padding;

        % Create uniform grid
        step = 2.0; % Å
        [x, y, z] = meshgrid(min_corner(1):step:max_corner(1), ...
                             min_corner(2):step:max_corner(2), ...
                             min_corner(3):step:max_corner(3));
        nodes = [x(:), y(:), z(:)];
    else
        % Assume mode is a file containing [x y z] points
        fprintf('>>> Loading potential points from file: %s\n', mode);
        nodes = load(mode);
        if size(nodes,2) ~= 3
            error('Points file must have 3 columns [x y z].');
        end
    end

    fprintf('>>> Computing potential for %d points...\n', size(nodes,1));
    tic;
    [potential, is_inside_solute, ~] = calculate_potentials1( ...
        nodes, particles_params, ...
        particles_coefficients, medium_params, n_max, t_eps0_inv);
    toc;

    fprintf('>>> Potential computation complete.\n');

    %% ==== Save results ====
    if strcmpi(mode, 'mesh')
        % ---- Save as Gaussian .cube format ----
        fprintf('>>> Saving potential map to file: potential_map.cube\n');

        % Reshape potential back to 3D grid
        nx = size(x,1);
        ny = size(x,2);
        nz = size(x,3);
        pot_grid = reshape(potential, [nx, ny, nz]);

        % Conversion constants
        coeff = 0.5291772108;  % 1 Bohr = 0.5291772108 Å
        step_bohr = step / coeff;  % convert step to Bohr
        origin_bohr = min_corner / coeff;  % convert origin to Bohr
        centers_bohr = reshape([particles_params.center], 3, []).' / coeff;

        % ---- Write the cube header ----
        fid = fopen('potential_map.cube','w');
        fprintf(fid, 'Cube file generated by compute_potential.m\n');
        fprintf(fid, 'Electrostatic potential map (values in internal units)\n');
        fprintf(fid, '%5d %12.6f %12.6f %12.6f\n', 0, origin_bohr(1), origin_bohr(2), origin_bohr(3));
        fprintf(fid, '%5d %12.6f %12.6f %12.6f\n', nx, step_bohr, 0, 0);
        fprintf(fid, '%5d %12.6f %12.6f %12.6f\n', ny, 0, step_bohr, 0);
        fprintf(fid, '%5d %12.6f %12.6f %12.6f\n', nz, 0, 0, step_bohr);

        % ---- Atom list (dummy atoms for visualization tools) ----
        for i = 1:numel(particles_params)
            c = centers_bohr(i,:);
            fprintf(fid, '%5d %12.6f %12.6f %12.6f %12.6f\n', 1, 0.0, c(1), c(2), c(3));
        end

        % ---- Write potential values (6 per line) ----
        count = 0;
        for ix = 1:nx
            for iy = 1:ny
                for iz = 1:nz
                    fprintf(fid, '%13.5e', real(pot_grid(ix,iy,iz)));
                    count = count + 1;
                    if mod(count,6) == 0
                        fprintf(fid, '\n');
                    end
                end
            end
        end
        fclose(fid);
        fprintf('>>> Saved potential map to "potential_map.cube" (in Bohr units)\n');
    else
        % ---- Save potential at given points as .dat ----
        fprintf('>>> Saving potential data to file: potential_points.dat\n');
        data = [nodes, real(potential)];
        fid = fopen('potential_points.dat','w');
        fprintf(fid, '# x(Å) y(Å) z(Å) potential\n');
        fprintf(fid, '%12.6f %12.6f %12.6f %13.6e\n', data.');
        fclose(fid);
        fprintf('>>> Saved potential data to "potential_points.dat"\n');
    end
end

%% ==== Visualization ====
fprintf('>>> Displaying sphere configuration...\n');
figure;
[Xid, Yid, Zid] = sphere;
for i = 1:numel(particles_params)
    X = Xid * particles_params(i).radius;
    Y = Yid * particles_params(i).radius;
    Z = Zid * particles_params(i).radius;
    surf(X + particles_params(i).center(1), ...
         Y + particles_params(i).center(2), ...
         Z + particles_params(i).center(3));
    hold on;
end
axis equal;
xlabel('X (Å)'); ylabel('Y (Å)'); zlabel('Z (Å)');
title('Sphere configuration');
view(3);
grid on;

fprintf('>>> Computation complete.\n');
end
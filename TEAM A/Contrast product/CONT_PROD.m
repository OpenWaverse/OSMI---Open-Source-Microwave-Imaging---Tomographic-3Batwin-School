close all;clear;clc
%

%%SETUP
here = fileparts(mfilename('fullpath'));                       % folder of this example
if ~isfolder(fullfile(here, '..', '..', 'Forward_solvers'))      % Live Editor: use the open file
    try, here = fileparts(matlab.desktop.editor.getActiveFilename); catch, here = pwd; end
end
addpath(fullfile(here, '..', '..'));
setup_osmi;
resultsDir = fullfile(here, 'results');
%% Add measurements
s_inc = sparameters('meas_1_background.s12p');
s_tot = sparameters('meas_12_pla_circle_pos_1.s12p');
%
%% FREQUENCIES
frequencies = [0.7e9 0.8e9 0.9e9 1e9 1.1e9]; 
%% MATERIAL
f = s_inc.Frequencies; %[0.5 - 1.25] GHz
ep_r_bg = linspace(45,60,51);  %[45-60]
sigma_bg = linspace(0.3,0.8,51); %[0.3 - 0.8] S/m
ep_r_t = linspace(1,10,51);  %[45-60]
sigma_t = linspace(0.01,0.06,51); %[0.3 - 0.8] S/m
idx = [11;16;21;26;31]; 
%
%% GEOMETRY
geom.doiRadius   = 0.10;
geom.antennaRing = 0.12;
geom.nAntennas   = 12;
geom.targetRadius = 0.01;
geom.targetCenter = [0.03, 0.03];
%
opts.name      = "domain";
opts.outputDir = resultsDir;
opts.current   = 1e-3;                         % A
%
load("Scenario4DS.mat")  %COMMON MESH
DELTA_S = zeros(12,12,5);
for i = 1:length(idx)
    sim.f = frequencies(i);
    materials.background.epsr  = ep_r_bg(idx(i));
    materials.background.sigma = sigma_bg(idx(i));             % S/m
    %
    materials.doi.epsr  = ep_r_bg(idx(i));
    materials.doi.sigma = sigma_bg(idx(i));
    %
    materials.target.epsr  = ep_r_bg(idx(i));
    materials.target.sigma = sigma_bg(idx(i));
    % INCIDENT FIELDS
    domain_inc = generate_inc_fields(sim, materials, geom, meshControl, opts, scenario);
    % TOTAL FIELDS
    opts.name      = "target";
    opts.outputDir = resultsDir;
    opts.current   = 1e-3;   % A
    %
    materials.target.epsr  = ep_r_t(idx(i));
    materials.target.sigma = sigma_t(idx(i));
    domain_tot = generate_tot_fields(sim, materials, geom, meshControl, opts, scenario);
    %
    DELTA_S(:,:,i) = domain_tot.Eprobes - domain_inc.Eprobes;
end
%% OPERATOR
REC = zeros(18134,1,5);
load("scenario4OP.mat")
for i = 1:length(idx)
    sim.f = frequencies(i);
    materials.background.epsr  = ep_r_bg(idx(i));
    materials.background.sigma = sigma_bg(idx(i));             % S/m
    %
    materials.doi.epsr  = ep_r_bg(idx(i));
    materials.doi.sigma = sigma_bg(idx(i));
    %
    materials.target.epsr  = ep_r_t(idx(i));
    materials.target.sigma = sigma_t(idx(i));
    %
    domain = generate_inc_fields(sim, materials, geom, meshControl, opts,scenario);
    Ecells = domain.Ecells;
    %% Build the Born Linear Operator K (matching OSMI's exact convention)
    doiIdx    = domain.mesh.cellsDOI;
    E_inc_doi = domain.Ecells(doiIdx, :);      % [Ncell x Nant]
    dA        = domain.mesh.area(doiIdx);      % [Ncell x 1]
    %
    N_ant  = geom.nAntennas;
    [n, m] = find(triu(true(N_ant)));          % n = Tx (col), m = Rx (row), includes n == m
    %    
    eps0  = 8.854187817e-12;
    omega = 2 * pi * sim.f;
    %  
    E_tx = E_inc_doi(:, n);                    % [Ncell x Npairs]
    E_rx = E_inc_doi(:, m);                    % [Ncell x Npairs]
    %
    opts.threshold = 0.999;
    opts.plotSVD   = true;
    %
    rec = born_tsvd(DELTA_S(:,:,i), domain.Ecells(doiIdx, :), dA, ...
                    sim.f, opts.current, opts);
    REC(:,1,i) = rec.delta_eps;
end
%% PRODUCT
DX_prod = ones(size(REC(:,1,1)));
for i = 1:length(idx)
    DX_prod = DX_prod.*REC(:,1,i);
end
DX_prod_max = max(DX_prod,[],'all');
DX_CP = abs(DX_prod)/DX_prod_max;
%
delta_eps = DX_CP;
%
doiIdx    = domain.mesh.cellsDOI;
verts = domain.mesh.verts;
cells = domain.mesh.cells(doiIdx, :); 
plotContrast(verts, cells, real(rec.delta_eps), "Recovered $\mathrm{Re}\{\Delta\varepsilon_r\}$");
plotContrast(verts, cells, imag(rec.delta_eps), "Recovered $\mathrm{Im}\{\Delta\varepsilon_r\}$");
plotContrast(verts, cells, abs(rec.delta_eps), "Recovered $\mathrm{Im}\{\Delta\varepsilon_r\}$");
    




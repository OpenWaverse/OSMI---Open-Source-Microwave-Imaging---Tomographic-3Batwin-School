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
ep_r_t = linspace(1,10,51);  %[1-10]
sigma_t = linspace(0.01,0.06,51); %[0.01 - 0.06] S/m
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
    %%%%%%%%%%%%%%%%%% SIMULATION %%%%%%%%%%%%%%%%%%
    % sim.f = frequencies(i);
    % materials.background.epsr  = ep_r_bg(idx(i));
    % materials.background.sigma = sigma_t(idx(i));             % S/m
    % %
    % materials.doi.epsr  = ep_r_bg(idx(i));
    % materials.doi.sigma = sigma_bg(idx(i));
    % %
    % materials.target.epsr  = ep_r_bg(idx(i));
    % materials.target.sigma = sigma_bg(idx(i));
    % % INCIDENT FIELDS
    % domain_inc = generate_inc_fields(sim, materials, geom, meshControl, opts, scenario);
    % % TOTAL FIELDS
    % opts.name      = "target";
    % opts.outputDir = resultsDir;
    % opts.current   = 1e-3;   % A
    % %
    % materials.target.epsr  = ep_r_t(idx(i));
    % materials.target.sigma = sigma_t(idx(i));
    % domain_tot = generate_tot_fields(sim, materials, geom, meshControl, opts, scenario);
    % %
    % DELTA_S(:,:,i) = domain_tot.Eprobes - domain_inc.Eprobes;
    %%%%%%%%%%%%%%%%%%%%%%% REAL MEAUSUREMENTS %%%%%%%%%%%%%%%%%%%
    
    DELTA_S(:,:,i) = s_tot.Parameters(:,:,idx(i)) - s_inc.Parameters(:,:,idx(i));
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
plotContrast(verts, cells, abs(rec.delta_eps), "Recovered $\mathrm{Abs}\{\Delta\varepsilon_r\}$");
%
%%METRICS
% 
% here = pwd;
% % library = fullfile(here,"Forward_solvers\lib");
% % addpath(library)
% results = fullfile(here,"\results");
% addpath(results)
% operatorName  = "operator";   % E_inc on the reconstruction mesh
% stateName     = "total";      % measured state (with target)
% referenceName = "incident";   % reference state (no target)
% 
% operator  = load_fields(operatorName,  results);
% state     = load_fields(stateName,     results);
% reference = load_fields(referenceName, results);
% %
% verts   = operator.mesh.verts;
% cells   = operator.mesh.cells;
% doi     = operator.mesh.cellsDOI;
% dA      = operator.mesh.area(doi);
% freq    = operator.freq;
% current = operator.current;
% omega   = 2 * pi * freq;
% eps0    = 8.854187817e-12;
% %
% % Baseline background permittivity (matches Main2.m: epsr = 20, sigma = 0.02)
% epsc_bg       = 20 - 1j * 0.02 / (omega * eps0);
% perm_curr     = repmat(epsc_bg, size(cells, 1), 1);
% %
% [E_sc, truth] = differential_data(state, reference, operator);
% % Extract ground truth complex contrast vector B on the operator DOI cells
% if isstruct(truth) && isfield(truth, 'delta_eps')
%     truth_vec = truth.delta_eps(:);
% elseif isnumeric(truth) && numel(truth) == numel(doi)
%     truth_vec = truth(:);
% else
%     % Fallback: exact ground truth on operator centroids (Main2.m: epsr=40, sigma=0.06)
%     centroids = (verts(cells(doi,1),:) + verts(cells(doi,2),:) + verts(cells(doi,3),:)) / 3;
%     in_tgt    = vecnorm(centroids - [0.03, 0.03], 2, 2) <= 0.01;
%     truth_vec = zeros(numel(doi), 1);
%     truth_vec(in_tgt) = (40 - 1j * 0.06 / (omega * eps0)) - epsc_bg;
% end
% %
% maxIter   = 10;
% stepAlpha = 0.65;
% 
% opts.threshold = 0.999;
% opts.plotSVD   = false;
% 
% nrmseReHist  = zeros(maxIter, 1);
% nrmseImHist  = zeros(maxIter, 1);
% nrmseAvgHist = zeros(maxIter, 1);
% [nrmseReHist(iter), nrmseImHist(iter), nrmseAvgHist(iter)] = compute_nrmse(delta_eps, truth_vec);
% 
% 
% 

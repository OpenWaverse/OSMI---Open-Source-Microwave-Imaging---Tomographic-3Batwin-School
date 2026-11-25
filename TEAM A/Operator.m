clc; close all; clear;
here = fileparts(mfilename('fullpath'));                       % folder of this example
if ~isfolder(fullfile(here, '..', '..', 'Forward_solvers'))      % Live Editor: use the open file
    try, here = fileparts(matlab.desktop.editor.getActiveFilename); catch, here = pwd; end
end
addpath(fullfile(here, '..', '..'));
setup_osmi;
resultsDir = fullfile(here, 'results');
sim.f = 1e9;                                   % Hz
materials.background.epsr  = 20;
materials.background.sigma = 0.02;             % S/m

materials.doi.epsr  = 20;
materials.doi.sigma = 0.02;
geom.doiRadius   = 0.10;
geom.antennaRing = 0.12;
geom.nAntennas   = 12;
c0         = 299792458;
lambda_doi = c0/(sim.f*sqrt(materials.doi.epsr));
lambda_bg  = c0/(sim.f*sqrt(materials.background.epsr));

meshControl.hDOI = lambda_doi/12;
meshControl.hBG  = lambda_bg/12;
meshControl.hAnt = lambda_bg/15;
meshControl.rAnt = meshControl.hAnt;
opts.name      = "domain";
opts.outputDir = resultsDir;
opts.current   = 1e-3;                         % A
[domain, scenario] = generate_fields(sim, materials, geom, meshControl, opts);

Ecells = domain.Ecells;

%% Build the Born Linear Operator K (matching OSMI's exact convention)
doiIdx    = domain.mesh.cellsDOI;
E_inc_doi = domain.Ecells(doiIdx, :);      % [Ncell x Nant]
dA        = domain.mesh.area(doiIdx);      % [Ncell x 1]

N_ant  = geom.nAntennas;
[n, m] = find(triu(true(N_ant)));          % n = Tx (col), m = Rx (row), includes n == m

eps0  = 8.854187817e-12;
omega = 2 * pi * sim.f;

E_tx = E_inc_doi(:, n);                    % [Ncell x Npairs]
E_rx = E_inc_doi(:, m);                    % [Ncell x Npairs]

% Operator K: [Npairs x Ncell]
operator = (1j * omega * eps0 / opts.current) * (E_tx.' .* E_rx.') .* dA(:).';

%doi = operator.mesh.cellsDOI;

opts.threshold = 0.999;
opts.plotSVD   = true;

load("DELTA_S.mat")

rec = born_tsvd(DELTA_S, domain.Ecells(doiIdx, :), dA, ...
      sim.f, opts.current, opts);
%
verts = domain.mesh.verts;
cells = domain.mesh.cells(doiIdx, :); 
plotContrast(verts, cells, real(rec.delta_eps), "Recovered $\mathrm{Re}\{\Delta\varepsilon_r\}$");
plotContrast(verts, cells, imag(rec.delta_eps), "Recovered $\mathrm{Im}\{\Delta\varepsilon_r\}$");

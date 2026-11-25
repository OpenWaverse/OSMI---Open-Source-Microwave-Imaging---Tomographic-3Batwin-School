clear all; close; clc
%%SETUP
here = fileparts(mfilename('fullpath'));                       % folder of this example
if ~isfolder(fullfile(here, '..', '..', 'Forward_solvers'))      % Live Editor: use the open file
    try, here = fileparts(matlab.desktop.editor.getActiveFilename); catch, here = pwd; end
end
addpath(fullfile(here, '..', '..'));
setup_osmi;
resultsDir = fullfile(here, 'results');
%% FREQUENCY
sim.f = 1e9;                                   % Hz
%% INCIDENT FIELDS
% MATERIAL
materials.background.epsr  = 20;
materials.background.sigma = 0.02;             % S/m

materials.doi.epsr  = 20;
materials.doi.sigma = 0.02;
% GEOMETRY
geom.doiRadius   = 0.10;
geom.antennaRing = 0.12;
geom.nAntennas   = 12;
% MESH CONTROL
c0         = 299792458;
lambda_doi = c0/(sim.f*sqrt(materials.doi.epsr));
lambda_bg  = c0/(sim.f*sqrt(materials.background.epsr));

meshControl.hDOI = lambda_doi/20;
meshControl.hBG  = lambda_bg/15;
meshControl.hAnt = lambda_bg/25;
meshControl.rAnt = meshControl.hAnt;
%
opts.name      = "domain";
opts.outputDir = resultsDir;
opts.current   = 1e-3;                         % A
% incident fields
[domain_inc, scenario] = generate_inc_fields(sim, materials, geom, meshControl, opts);
%% TOTAL FIELDS
% MATERIAL
materials.background.epsr  = 20;
materials.background.sigma = 0.01;             % S/m

materials.doi.epsr  = 20;
materials.doi.sigma = 0.02;

materials.target.epsr  = 60;
materials.target.sigma = 0.5;
%GEOMETRY
geom.doiRadius    = 0.10;
geom.antennaRing  = 0.12;
geom.nAntennas    = 12;
geom.targetRadius = 0.01;
geom.targetCenter = [0.03, 0.03];
%
opts.name      = "target";
opts.outputDir = resultsDir;
opts.current   = 1e-3;                         % A
domain_tot = generate_tot_fields(sim, materials, geom, meshControl, opts, scenario);

DELTA_S = domain_tot.Eprobes - domain_inc.Eprobes;

save("DELTA_S.mat","DELTA_S")
%%
% figure
% imagesc(abs(DELTA_S)) 
% axis square
% axis tight
% cb = colorbar;
% ylabel(cb,'Magnitude (dB)')
% cb.Label.Interpreter = 'latex';
% cb.Label.FontSize = 12;
% 
% xlabel('Port j')
% ylabel('Port i')
% title('Scattering Parameters (dB)')

% %% MODIFICATION PATTERN
% DELTA_S = DELTA_S - diag(diag(DELTA_S));
% figure
% imagesc(abs(DELTA_S)) 
% axis square
% cb = colorbar;
% ylabel(cb,'Magnitude (dB)')
% cb.Label.Interpreter = 'latex';
% cb.Label.FontSize = 12;
% 
% xlabel('Port j')
% ylabel('Port i')
% title('Scattering Parameters (dB)')

%




%INVERSION
doi = operator.mesh.cellsDOI;

opts.threshold = 0.999;
opts.plotSVD   = true;
rec = born_tsvd(E_sc, operator.Ecells(doi, :), operator.mesh.area(doi), ...
                operator.freq, operator.current, opts);
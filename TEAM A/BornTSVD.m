% To be ran after the Main2

clc; close all; clear;
here = fileparts(mfilename('fullpath'));                       % folder of this example
if ~isfolder(fullfile(here, '..', '..', 'Forward_solvers'))      % Live Editor: use the open file
    try, here = fileparts(matlab.desktop.editor.getActiveFilename); catch, here = pwd; end
end
addpath(fullfile(here, '..', '..'));
setup_osmi;
resultsDir = fullfile(here, 'results');

operatorName  = "operator";   % E_inc on the reconstruction mesh
stateName     = "total";      % measured state (with target)
referenceName = "incident";   % reference state (no target)

operator  = load_fields(operatorName,  resultsDir);
state     = load_fields(stateName,     resultsDir);
reference = load_fields(referenceName, resultsDir);

[E_sc, truth] = differential_data(state, reference, operator);

doi = operator.mesh.cellsDOI;

opts.threshold = 0.999;
opts.plotSVD   = true;
rec = born_tsvd(E_sc, operator.Ecells(doi, :), operator.mesh.area(doi), ...
    operator.freq, operator.current, opts);

verts = operator.mesh.verts;
cells = operator.mesh.cells(doi, :);

plotContrast(verts, cells, real(rec.delta_eps), "Recovered $\mathrm{Re}\{\Delta\varepsilon_r\}$");
plotContrast(verts, cells, imag(rec.delta_eps), "Recovered $\mathrm{Im}\{\Delta\varepsilon_r\}$");

img_norm = abs(rec.delta_eps) ./ max(abs(rec.delta_eps));
plotContrast(verts, cells, img_norm, "Normalised $|\Delta\varepsilon_r|$");

binThreshold = 10^(-5/20);
plotContrast(verts, cells, double(img_norm > binThreshold), "Binarised $|\Delta\varepsilon_r|$ ($-5$ dB)");
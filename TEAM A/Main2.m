clc; close all; clear;
here = "C:\Users\ibeb\OneDrive - Ciências ULisboa\Doutoramento\Coisas uteis\School_Turim\Matlab\TEAM A";                       % folder of this example
addpath(fullfile(here, '..', '..'));
setup_osmi;
resultsDir = fullfile(here, 'results');

sim.f = 1e9; 

materials.background.epsr  = 50;
materials.background.sigma = 0.01;             % S/m

materials.doi.epsr  = 50;
materials.doi.sigma = 0.01;

materials.target.epsr  = 10;
materials.target.sigma = 0.5;

geom.doiRadius    = 0.078;
geom.antennaRing  = 0.078;
geom.nAntennas    = 12;
geom.targetRadius = 0.01;
geom.targetCenter = [0.03, 0.03];

c0         = 299792458;
lambda_doi = c0/(sim.f*sqrt(materials.doi.epsr));
lambda_bg  = c0/(sim.f*sqrt(materials.background.epsr));
lambda_tgt = c0/(sim.f*sqrt(materials.target.epsr));

meshData.hDOI    = lambda_doi/30;
meshData.hBG     = lambda_bg/24;
meshData.hAnt    = lambda_bg/30;
meshData.hTarget = lambda_tgt/30;
meshData.rAnt    = meshData.hAnt;

meshRec         = meshData;
meshRec.hDOI    = lambda_doi/24;
meshRec.hTarget = lambda_tgt/24;

opts.outputDir     = resultsDir;
opts.current       = 1e-3;                     % A
opts.name          = "incident";
opts.includeTarget = false;
incident = generate_fields(sim, materials, geom, meshData, opts);

opts.name          = "total";
opts.includeTarget = true;
[total, scenario] = generate_fields(sim, materials, geom, meshData, opts);

opts.name          = "operator";
opts.includeTarget = false;
operator = generate_fields(sim, materials, geom, meshRec, opts);

tx = 1;

plotField(incident, tx, 'abs', 'Title', "$|E_{inc}|$, Tx " + tx);
plotField(total,    tx, 'abs', 'Title', "$|E_{tot}|$, Tx " + tx);
plotField(total,    tx, 'abs', 'Reference', incident, 'Title', "$|E_{sc}|$, Tx " + tx);
plotField(total,    tx, 'phase', 'Reference', incident, 'Title', "$\angle E_{sc}$, Tx " + tx);

E_sc = total.Eprobes - incident.Eprobes;
figure;
imagesc(20*log10(abs(E_sc)/max(abs(E_sc(:)))));
axis image; colorbar;
xlabel('Transmitter $n$'); ylabel('Receiver $m$');
title('$|E_{sc}(m,n)|$ (dB)');
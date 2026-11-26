clc; close all; clear;
here = "C:\Users\ibeb\OneDrive - Ciências ULisboa\Doutoramento\Coisas uteis\School_Turim\Matlab\TEAM A";                       % folder of this example
addpath(fullfile(here, '..', '..'));
setup_osmi;
resultsDir = fullfile(here, 'results');

sim.f = 800e6; 

materials.background.epsr  = 50;
materials.background.sigma = 0.3;             % S/m

materials.doi.epsr  = 50;
materials.doi.sigma = 0.3;

materials.target.epsr  = 10;
materials.target.sigma = 0.06;

geom.doiRadius    = 0.078;
geom.antennaRing  = 0.078;
geom.nAntennas    = 12;
geom.targetRadius = 0.01;
geom.targetCenter = [0, 0];

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
opts.name          = "incident_coeff";
opts.includeTarget = false;
incident = generate_fields(sim, materials, geom, meshData, opts);

opts.name          = "total_coeff";
opts.includeTarget = true;
[total, scenario] = generate_fields(sim, materials, geom, meshData, opts);

opts.name          = "operator_coeff";
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


%% Computing the experimental scattering matrix

[S_bg, freqs] = read_s12p('C:\Users\ibeb\OneDrive - Ciências ULisboa\Doutoramento\Coisas uteis\School_Turim\Matlab\Examples\Measurements_lab\meas_1_background.s12p');
[S_target, freqs] = read_s12p('C:\Users\ibeb\OneDrive - Ciências ULisboa\Doutoramento\Coisas uteis\School_Turim\Matlab\Examples\Measurements_lab\meas_2_pla_circle_center.s12p');

S_bg = S_bg(:, :, 1:numel(freqs));

% Find the index matching sim.f (make sure sim.f is between 500e6 and 840e6)
f_idx = find(abs(freqs - sim.f) < 1e3, 1);
assert(~isempty(f_idx), 'sim.f (%.1f MHz) is not in the measured range (500-840 MHz).', sim.f/1e6);

S_bg_f     = S_bg(:, :, f_idx);       % 12x12 measured background
S_target_f = S_target(:, :, f_idx);   % 12x12 measured PLA target
S_sc_f     = S_target_f - S_bg_f;     % 12x12 measured differential scattering

%% 2. Mask diagonal self-reflections and unmeasured 0+0j switch pairs
% Odd-to-odd and even-to-even off-diagonal entries are 0+0j due to the RF switch matrix,
% and diagonal entries (m == n) represent antenna mismatch rather than transmission.
valid_mask = (abs(S_bg_f) > 0) & (abs(S_target_f) > 0) & ~eye(geom.nAntennas);

% Compute calibration matrix using high-SNR total fields (instead of noisy differential S_sc_f)
C_tot = zeros(geom.nAntennas);
C_tot(valid_mask) = total.Eprobes(valid_mask) ./ S_target_f(valid_mask);

% (Alternative: background-only calibration if your simulated PLA parameters are uncertain)
% C_inc = zeros(geom.nAntennas);
% C_inc(valid_mask) = incident.Eprobes(valid_mask) ./ S_bg_f(valid_mask);

%% Save calibration matrix for use in other scripts
cal_matrix = C_tot;
f_cal      = sim.f;

if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

saveFile = fullfile(resultsDir, 'calibration_matrix.mat');
save(saveFile, 'cal_matrix', 'valid_mask', 'f_cal', 'S_bg_f');
fprintf('Saved calibration matrix to: %s\n', saveFile);





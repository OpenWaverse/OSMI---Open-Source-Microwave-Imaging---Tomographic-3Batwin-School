% To be ran after Main2.m (hybrid DBIM + BornTSVD reconstruction)

clc; close all; clear;
here = fileparts(mfilename('fullpath'));
if ~isfolder(fullfile(here, '..', '..', 'Forward_solvers'))
    try, here = fileparts(matlab.desktop.editor.getActiveFilename); catch, here = pwd; end
end
addpath(fullfile(here, '..', '..'));
setup_osmi;
resultsDir = fullfile(here, 'results');

% Load calibration matrix and background slice already saved at f_cal
load(fullfile(resultsDir, 'calibration_matrix.mat'), 'cal_matrix', 'valid_mask', 'f_cal', 'S_bg_f');

% Load only the unknown target measurement and slice at f_cal
[S_unk, freqs] = read_s12p('C:\Users\ibeb\OneDrive - Ciências ULisboa\Doutoramento\Coisas uteis\School_Turim\Matlab\Examples\Measurements_lab\meas_2_pla_circle_center.s12p');
f_idx   = find(abs(freqs - f_cal) < 1e3, 1);
assert(~isempty(f_idx), 'f_cal (%.1f MHz) not found in target file.', f_cal/1e6);
S_unk_f = S_unk(:, :, f_idx);

% Apply element-wise calibration and zero-fill unmeasured entries (no NaNs)
E_sc = cal_matrix .* (S_unk_f - S_bg_f);
E_sc(~valid_mask) = 0;

%% 1. Load precomputed fields from Main2.m (exact same setup as BornTSVD.m)
operatorName = "operator_coeff";  % E_inc on the reconstruction mesh
stateName     = "total_coeff";      % measured state (with target)
referenceName = "incident_coeff";   % reference state (no target)

operator  = load_fields(operatorName,  resultsDir);
state     = load_fields(stateName,     resultsDir);
reference = load_fields(referenceName, resultsDir);

% Symmetrise scattered matrix via Lorentz reciprocity (+3 dB SNR)
E_sc = 0.5 * (E_sc + E_sc.');

%% 2. Extract mesh info & locate antenna node indices inside operator.mesh
verts   = operator.mesh.verts;
cells   = operator.mesh.cells;
doi     = operator.mesh.cellsDOI;
dA      = operator.mesh.area(doi);
freq    = operator.freq;
current = operator.current;
omega   = 2 * pi * freq;
eps0    = 8.854187817e-12;

% Mask to consider the available antenna pairs (odd in first switching
% matrics, even in the second switching matrix)
N_ant         = size(E_sc, 1);
switchID      = mod(0:N_ant-1, 2) + 1;
mask_bistatic = (switchID(:) ~= switchID(:).');   % 72 active entries
mask_unique   = triu(mask_bistatic, 1);

% Locate probe node indices stored inside operator.mesh by generate_fields
meshFields = fieldnames(operator.mesh);
probeField = meshFields(contains(lower(meshFields), 'probe'));
assert(~isempty(probeField), 'No probe field found in operator.mesh. Available fields: %s', strjoin(meshFields, ', '));
probeNodes = operator.mesh.(probeField{1});

% Baseline background permittivity (matches Main2.m: epsr = 20, sigma = 0.02)
epsc_bg       = 50 - 1j * 0.3 / (omega * eps0);
perm_curr     = repmat(epsc_bg, size(cells, 1), 1);
delta_eps_acc = zeros(numel(doi), 1);

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

%% 3. Hybrid DBIM + BornTSVD Loop
maxIter   = 10;
stepAlpha = 0.85;

opts.threshold = 0.98;
opts.plotSVD   = false;

nrmseReHist  = zeros(maxIter, 1);
nrmseImHist  = zeros(maxIter, 1);
nrmseAvgHist = zeros(maxIter, 1);

currFields = operator;
E_inc_oper = operator.Eprobes;
tStart     = tic;

for iter = 1:maxIter
    % Residual data: Iteration 1 = full E_sc (exact linear BornTSVD)
    if iter == 1
        E_res = E_sc;
    else
        E_sim_sc = 0.5 * ((currFields.Eprobes - E_inc_oper) + (currFields.Eprobes - E_inc_oper).');
        E_res    = E_sc - E_sim_sc;
    end

    opts.plotSVD = (iter == 1);
    %rec_step = born_tsvd(E_res, currFields.Ecells(doi, :), dA, freq, current, opts);
    %rec_step = born_tsvd_masked(E_res, currFields.Ecells(doi, :), dA, freq, current, mask_unique, opts);
    rec_step = born_tsvd_masked(E_res, currFields.Ecells(doi, :), dA, freq, current, mask_unique, opts);

    % Accumulate update & clamp physical bounds (Re >= 1, Im <= 0)
    delta_eps_acc = delta_eps_acc + stepAlpha * rec_step.delta_eps;
    eps_tot       = epsc_bg + delta_eps_acc;
    delta_eps_acc = (max(1.0, real(eps_tot)) + 1j * min(0.0, imag(eps_tot))) - epsc_bg;
    

    if iter == 1
        delta_eps_born = delta_eps_acc;
    end

    % Evaluate Jury NRMSE(q) for real and imaginary components
    %[nrmseReHist(iter), nrmseImHist(iter), nrmseAvgHist(iter)] = compute_nrmse(delta_eps_acc, truth_vec);

    fprintf('Iter %d/%d | Res: %.4f | NRMSE(Re): %.4f | NRMSE(Im): %.4f | Mean: %.4f\n', ...
            iter, maxIter, norm(E_res,'fro')/norm(E_sc,'fro'), ...
            nrmseReHist(iter), nrmseImHist(iter), nrmseAvgHist(iter));

    % Re-solve forward problem on operator mesh for iterations 2..maxIter
    if iter < maxIter
        perm_curr(doi) = epsc_bg + delta_eps_acc;
        currFields = solve_fields(verts, cells, perm_curr, freq, probeNodes, ...
                                  epsc_bg, current, doi, "dbim_step");
    end
end

elapsedTime = toc(tStart);

fprintf('\n=== Jury Metric Summary (BornTSVD vs Hybrid DBIM) ===\n');
fprintf('Execution time    : %.2f s\n', elapsedTime);
fprintf('Iter 1 (BornTSVD) -> NRMSE(Re): %.4f | NRMSE(Im): %.4f | Mean: %.4f\n', ...
        nrmseReHist(1), nrmseImHist(1), nrmseAvgHist(1));
fprintf('Iter %d (DBIM)     -> NRMSE(Re): %.4f | NRMSE(Im): %.4f | Mean: %.4f\n', ...
        maxIter, nrmseReHist(end), nrmseImHist(end), nrmseAvgHist(end));

%% 4. Visualisation (matching BornTSVD.m)
cellsDOI = cells(doi, :);

plotContrast(verts, cellsDOI, real(delta_eps_born), ...
    sprintf("Iter 1 (BornTSVD) $\\mathrm{Re}\\{\\Delta\\varepsilon_r\\}$ (NRMSE=%.3f)", nrmseReHist(1)));
plotContrast(verts, cellsDOI, real(delta_eps_acc), ...
    sprintf("Iter %d (DBIM) $\\mathrm{Re}\\{\\Delta\\varepsilon_r\\}$ (NRMSE=%.3f)", maxIter, nrmseReHist(end)));
plotContrast(verts, cellsDOI, imag(delta_eps_acc), ...
    sprintf("Iter %d (DBIM) $\\mathrm{Im}\\{\\Delta\\varepsilon_r\\}$ (NRMSE=%.3f)", maxIter, nrmseImHist(end)));

img_norm = abs(delta_eps_acc) ./ max(abs(delta_eps_acc));
plotContrast(verts, cellsDOI, img_norm, "DBIM Normalised $|\Delta\varepsilon_r|$");

binThreshold = 10^(-5/20);
plotContrast(verts, cellsDOI, double(img_norm > binThreshold), "DBIM Binarised $|\Delta\varepsilon_r|$ ($-5$ dB)");


function [S, freqs] = read_s12p(filename)
% READ_S12P Parses a 12-port Touchstone file into a 12x12xN_freq complex array.
% Outputs:
%   S     - 12x12xN_freq complex S-parameter matrix ordered as S(rx, tx, freq)
%   freqs - 1xN_freq frequency vector in Hz

lines = strtrim(readlines(filename));

% Default Touchstone header settings (# GHz S MA R 50)
freq_scale = 1e9;
fmt = 'MA';

% Parse option line (#)
opt_idx = find(startsWith(lines, '#'), 1);
if ~isempty(opt_idx)
    tokens = upper(split(lines(opt_idx)));
    if any(tokens == "HZ"),  freq_scale = 1;   end
    if any(tokens == "KHZ"), freq_scale = 1e3; end
    if any(tokens == "MHZ"), freq_scale = 1e6; end
    if any(tokens == "GHZ"), freq_scale = 1e9; end

    if any(tokens == "RI"), fmt = 'RI'; end
    if any(tokens == "MA"), fmt = 'MA'; end
    if any(tokens == "DB"), fmt = 'DB'; end
end

% Keep only data lines
data_lines = lines(~startsWith(lines, {'!', '#'}) & strlength(lines) > 0);

freqs = [];
raw_vals = [];
for k = 1:numel(data_lines)
    % Strip inline comments if present
    clean_line = extractBefore(data_lines(k) + "!", "!");
    nums = str2double(split(strtrim(clean_line)));
    nums = nums(~isnan(nums));
    if isempty(nums), continue; end

    % A frequency header line always has an odd number of tokens (1 freq + 2*N pairs)
    if mod(numel(nums), 2) == 1
        freqs(end+1) = nums(1) * freq_scale; %#ok<AGROW>
        raw_vals = [raw_vals; nums(2:end)]; %#ok<AGROW>
    else
        raw_vals = [raw_vals; nums]; %#ok<AGROW>
    end
end

% Convert value pairs into complex numbers
v1 = raw_vals(1:2:end);
v2 = raw_vals(2:2:end);
switch fmt
    case 'RI'
        c_vals = v1 + 1j * v2;
    case 'MA'
        c_vals = v1 .* exp(1j * deg2rad(v2));
    case 'DB'
        c_vals = (10 .^ (v1 / 20)) .* exp(1j * deg2rad(v2));
end

% Keep only complete 12x12 blocks (144 complex entries per frequency)
n_freq = floor(numel(c_vals) / 144);
freqs = freqs(1:n_freq);

% Reshape row-major Touchstone data into S(rx, tx, freq)
S = permute(reshape(c_vals(1:n_freq*144), [12, 12, n_freq]), [2, 1, 3]);
end
% Compare two binary dumps of complex symbols (float32 interleaved I,Q per sample).
% Default: bypass dump vs Gardner on-time dump (same format as gardner_ontime.bin).
%
% Usage:
%   compare_viterbi_gardner_dump
%   Or set file names below, optional maxLagSymbols and maxSymbols.

clear; close all;

opts.bypassFile   = '../data/viterbi_input_bypass.bin';
opts.gardnerFile  = '../data/gardner_ontime.bin';
opts.maxLagSymbols = 512;   % search lag in [-maxLagSymbols, +maxLagSymbols]
opts.maxSymbols    = 1e7;   % max complex samples read from each file (from start)
dumpFirst = 0;

%% Read (float32 interleaved I,Q per row in file -> 2 x N matrix)
fid = fopen(opts.bypassFile, 'rb');
if fid < 0
    zb = []; nb = 0;
else
    raw = fread(fid, [2, dumpFirst], 'float32');
    raw = fread(fid, [2, inf], 'float32');
    fclose(fid);
    if isempty(raw)
        zb = []; nb = 0;
    else
        ntake = min(size(raw, 2), opts.maxSymbols);
        raw = raw(:, 1:ntake);
        zb = raw(1,:) + 1i * raw(2,:);
        zb = zb(:);
        nb = numel(zb);
    end
end

fid = fopen(opts.gardnerFile, 'rb');
if fid < 0
    zg = []; ng = 0;
else
    raw = fread(fid, [2, dumpFirst], 'float32');
    raw = fread(fid, [2, inf], 'float32');
    fclose(fid);
    if isempty(raw)
        zg = []; ng = 0;
    else
        ntake = min(size(raw, 2), opts.maxSymbols);
        raw = raw(:, 1:ntake);
        zg = raw(1,:) + 1i * raw(2,:);
        zg = zg(:);
        ng = numel(zg);
    end
end

fprintf('Bypass: %d complex samples, Gardner: %d complex samples\n', nb, ng);

if nb < 16 || ng < 16
    error('Files too short or missing. Check paths and BYPASS_GARDNER_DUMP_VITERBI_INPUT / DEBUG_GARDNER_OUTPUTS.');
end

%% Best integer-symbol lag: maximize |sum zb .* conj(zg_shifted)| / (||zb|| ||zg||)
% Convention: positive lag means Gardner is delayed vs bypass (compare zb(n) to zg(n+lag)).
L = min(opts.maxLagSymbols, floor(min(nb, ng) / 4));
bestLag = 0;
bestScore = -inf;
for lag = -L:L
    if lag >= 0
        n = min(nb, ng - lag);
        n=min(n,1000);
        if n < 8, continue; end
        a = zb(1:n);
        b = zg((1:n) + lag);
    else
        lg = -lag;
        n = min(nb - lg, ng);
        n=min(n,1000);
        if n < 8, continue; end
        a = zb((1:n) + lg);
        b = zg(1:n);
    end
    score = abs(sum(a .* conj(b))) / (sqrt(sum(abs(a).^2)) * sqrt(sum(abs(b).^2)) + eps);
    if score > bestScore
        bestScore = score;
        bestLag = lag;
    end
end

fprintf('Best lag (symbol indices, Gardner relative to bypass): %d\n', bestLag);
fprintf('Best normalized correlation magnitude: %.6f\n', bestScore);

%% Overlay after alignment
if bestLag >= 0
    nShow = min(10000000, min(nb, ng - bestLag));
    a = zb(1:nShow);
    b = zg((1:nShow) + bestLag);
else
    lg = -bestLag;
    nShow = min(10000000, min(nb - lg, ng));
    a = zb((1:nShow) + lg);
    b = zg(1:nShow);
end

figure;
subplot(2,1,1);
plot(real(a), 'b'); hold on; plot(real(b), 'r--'); grid on;
title('Real (blue=bypass, red=Gardner aligned)');
legend('bypass', 'gardner');
subplot(2,1,2);
plot(imag(a), 'b'); hold on; plot(imag(b), 'r--'); grid on;
title('Imag'); grid on;

figure;
plot(abs(a-b), 'b'); hold on; grid on;
title('abs err aligned)');

plotConst=0;
if plotConst==1
    
    figure;
    scatter(real(a), imag(a), 0.5, 'b'); hold on;
    scatter(real(b), imag(b), 0.5, 'r');
    axis equal; grid on;
    title('Constellations (aligned window)');
    legend('bypass', 'gardner');

end

%% NMSE in aligned window
nmse = mean(abs(a - b).^2) / (mean(abs(a).^2) + eps);
fprintf('NMSE (complex) on aligned window: %.3e\n', nmse);

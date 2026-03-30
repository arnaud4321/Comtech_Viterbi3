function plot_symbol_rate_mags()
% Plot SymbolRateEstimator mags dump (peak detection spectrum around DC).
%
% Reads ../data/symbol_rate_mags.bin written when DEBUG_SYMBOL_RATE_ESTIMATOR_DUMP is enabled.
% Format:
%   int32 len
%   int32 kMax
%   int32 bestK
%   double FsHzUsed
%   int32 NfftUsed
%   double FsHzBase
%   int32 NfftBase
%   double OsFactor
%   int32 nTime
%   float32 I_time[nTime]
%   float32 Q_time[nTime]
%   float32 mags_zoom[len]   where mags_zoom corresponds to k = -kMax:kMax
%   float32 mags_full[NfftUsed]  where mags_full corresponds to k = 0:NfftUsed-1 (raw FFT bins)

fname = fullfile('..','data','symbol_rate_mags.bin');
fid = fopen(fname,'rb');
if fid < 0
    error('Cannot open %s (run with DEBUG_SYMBOL_RATE_ESTIMATOR_DUMP).', fname);
end

len   = fread(fid, 1, 'int32');
kMax  = fread(fid, 1, 'int32');
bestK = fread(fid, 1, 'int32');
FsHzUsed  = fread(fid, 1, 'double');
NfftUsed  = fread(fid, 1, 'int32');
FsHzBase  = fread(fid, 1, 'double');
NfftBase  = fread(fid, 1, 'int32');
OsFactor  = fread(fid, 1, 'double');
nTime     = fread(fid, 1, 'int32');
I_time    = fread(fid, double(nTime), 'single');
Q_time    = fread(fid, double(nTime), 'single');
mags_zoom  = fread(fid, double(len), 'single');
mags_full  = fread(fid, double(NfftUsed), 'single');
fclose(fid);

if length(mags_zoom) ~= len
    error('Short read: expected %d mags_zoom, got %d.', len, length(mags_zoom));
end
if length(mags_full) ~= NfftUsed
    error('Short read: expected %d mags_full, got %d.', NfftUsed, length(mags_full));
end
if len ~= (2*kMax + 1)
    warning('len != 2*kMax+1 (len=%d, kMax=%d).', len, kMax);
end
if length(I_time) ~= nTime || length(Q_time) ~= nTime
    error('Short read: expected nTime=%d, got I=%d Q=%d.', nTime, length(I_time), length(Q_time));
end

% Full-band frequency axis (unshifted, [0, Fs)).
k_full = (0:(NfftUsed-1)).';
f_full = double(k_full) * double(FsHzUsed) / double(NfftUsed); % Hz
mag_full_db = 10*log10(double(mags_full) + 1e-30);

% Zoom around FsBase/2 (physical), expressed on the (possibly oversampled) grid.
k0 = round((0.5 * FsHzBase) * double(NfftUsed) / double(FsHzUsed)); % 0-based bin index
idx_zoom = (k0 - kMax):(k0 + kMax);   % 0-based indices
idx_zoom1 = idx_zoom + 1;            % MATLAB 1-based
f_zoom = f_full(idx_zoom1);
mag_zoom_db = mag_full_db(idx_zoom1);

% Optional sanity check: compare dumped mags_zoom with the extracted zoom from full-band.
mag_zoom_dump_db = 10*log10(double(mags_zoom) + 1e-30);
if length(mag_zoom_dump_db) == length(mag_zoom_db)
    diff_db = max(abs(mag_zoom_dump_db - mag_zoom_db));
    if diff_db > 1e-6
        warning('Zoom mismatch: dumped mags_zoom vs extracted full-band differs by up to %.3g dB.', diff_db);
    end
end

figure;
plot(f_full/1e6, mag_full_db, 'Color', [0.7 0.7 0.7], 'LineWidth', 1.0); hold on;
plot(f_zoom/1e6, mag_zoom_db, 'b', 'LineWidth', 1.6);
grid on;
xlabel('Frequency (MHz)');
ylabel('Magnitude (dB, arbitrary)');
title(sprintf('SymbolRateEstimator mags (FsUsed=%.3f MHz, NfftUsed=%d, bestK=%d, Os=%.2f, FsBase=%.3f MHz, NfftBase=%d)', ...
    FsHzUsed/1e6, NfftUsed, bestK, OsFactor, FsHzBase/1e6, NfftBase));
% Mark bestK
fBest = double(bestK) * double(FsHzUsed) / double(NfftUsed);
yl = ylim;
plot([fBest fBest]/1e6, yl, '--r', 'LineWidth', 1.0);
fHalfBase = 0.5 * FsHzBase;
plot([fHalfBase fHalfBase]/1e6, yl, ':k', 'LineWidth', 1.0);
legend('full-band mags','zoom (FsBase/2 \pm kMax)','bestK','FsBase/2','Location','best');

figure;
subplot(2,1,1);
plot(double(I_time)); grid on; xlabel('n'); ylabel('I');
title(sprintf('I_{time} used for FFT (nTime=%d)', nTime));
subplot(2,1,2);
plot(double(Q_time)); grid on; xlabel('n'); ylabel('Q');
title('Q_{time} used for FFT');

sigIn=complex(I_time,Q_time);
amp = abs(sigIn);
figure;
subplot(2,1,1);
plot(double(amp)); grid on; xlabel('n'); ylabel('|I+jQ|');
title('Amplitude |I+jQ| used for non-linearity');
subplot(2,1,2);
N = length(amp);
f_amp = double(0:(N-1)) * double(FsHzUsed) / double(N); % Hz
plot(f_amp/1e6, abs(fft(double(amp)))); grid on;
xlabel('Frequency (MHz)'); ylabel('|FFT(amplitude)|');
title('FFT of amplitude (no DC removal)');


% NB: le code C++ applique déjà le suréchantillonnage (si OsFactor>1).
% Ici on se contente d'afficher quelques indicateurs utiles.
fprintf('Header: FsUsed=%.6f Hz, NfftUsed=%d, FsBase=%.6f Hz, NfftBase=%d, OsFactor=%.3f\n', ...
    FsHzUsed, NfftUsed, FsHzBase, NfftBase, OsFactor);
end


% Plot soft symbols constellation, phase error and phase estimate (DEBUG_PHASE_CORR dumps).
% softSymbols.bin: 3 streams per batch, each stream BatchSize1=1024 I,Q floats; batch = 3072 symbols = 6144 floats.
% phaseError.bin: one float (rad) per phase update (one per filter batch of SPB samples).
% phaseEstRad.bin: one float (rad) per phase update - applied phase correction.
% Symbol rate 10.71e6 sym/s; phase updates at Fs/SPB per second.

SymbolRate = 10.71e6;
Fs = 21.42e6;
SPB = 8192;
BatchSize1 = 1024;
BatchSize3 = 3 * BatchSize1;   % 3072 symbols per batch
NumSymbols = 1000000;
NumPhaseSamples = 2000;        % number of phase error samples to plot

time_sec = input('Start time (seconds, relative to run start): ');
if isempty(time_sec)
    time_sec = 0;
end

%% Soft symbols: 1000 symbols starting at time_sec (wrt symbol rate)
start_symbol = round(time_sec * SymbolRate);
fid = fopen('softSymbols.bin', 'rb');
if fid == -1
    fid = fopen('../data/softSymbols.bin', 'rb');
end
if fid == -1
    error('Cannot open softSymbols.bin');
end
% Layout: symbol index s -> float index 2*s (I,Q); byte offset = s*8
offset_bytes = start_symbol * 2 * 4;  % 2 floats per symbol, 4 bytes per float
fseek(fid, offset_bytes, 'bof');
d = fread(fid, NumSymbols * 2, 'float32');
fclose(fid);
if numel(d) < NumSymbols * 2
    warning('Requested %d symbols, got %d floats. File may be short.', NumSymbols, numel(d)/2);
    NumSymbols = floor(numel(d) / 2);
    d = d(1:NumSymbols*2);
end
I = d(1:2:end);
Q = d(2:2:end);

%% Phase error and phase estimate: evolution for same time range and a bit more
% Phase update rate = Fs/SPB per second; index i -> time i * (SPB/Fs)
phase_dt = SPB / Fs;
start_phase_idx = round(time_sec / phase_dt);
offset_phase_bytes = start_phase_idx * 4;  % one float per sample

fid2 = fopen('phaseError.bin', 'rb');
if fid2 == -1
    fid2 = fopen('../data/phaseError.bin', 'rb');
end
if fid2 == -1
    warning('Cannot open phaseError.bin; skipping phase error plot.');
    phase_err = [];
    phase_time = [];
else
    fseek(fid2, offset_phase_bytes, 'bof');
    phase_err = fread(fid2, NumPhaseSamples, 'float32');
    fclose(fid2);
    phase_time = (start_phase_idx + (0:numel(phase_err)-1))' * phase_dt;
end

fid3 = fopen('phaseEstRad.bin', 'rb');
if fid3 == -1
    fid3 = fopen('../data/phaseEstRad.bin', 'rb');
end
if fid3 == -1
    warning('Cannot open phaseEstRad.bin; skipping phase estimate plot.');
    phase_est_rad = [];
else
    fseek(fid3, offset_phase_bytes, 'bof');
    phase_est_rad = fread(fid3, NumPhaseSamples, 'float32');
    fclose(fid3);
    if isempty(phase_time)
        phase_time = (start_phase_idx + (0:numel(phase_est_rad)-1))' * phase_dt;
    end
end

%% Plot
figure('Name', 'Soft symbols & phase (DEBUG_PHASE_CORR)');

subplot(3,1,1);
plot(I, Q, '.');
xlabel('I');
ylabel('Q');
title(sprintf('Constellation (soft symbols at Viterbi input), 1000 syms from t=%.3f s', time_sec));
grid on;
axis equal;

subplot(3,1,2);
if ~isempty(phase_err)
    plot(phase_time, abs(phase_err), 'b.-');
    xlabel('Time (s)');
    ylabel('|Phase error| (rad)');
    title(sprintf('Phase error magnitude from t=%.3f s (%d samples)', time_sec, numel(phase_err)));
    grid on;
else
    text(0.5, 0.5, 'phaseError.bin not available', 'HorizontalAlignment', 'center');
end

subplot(3,1,3);
if ~isempty(phase_est_rad)
    plot(phase_time, phase_est_rad, 'r.-');
    xlabel('Time (s)');
    ylabel('Phase est. (rad)');
    title(sprintf('Phase estimate (applied correction) from t=%.3f s', time_sec));
    grid on;
else
    text(0.5, 0.5, 'phaseEstRad.bin not available', 'HorizontalAlignment', 'center');
end

% Plot FFT-based center freq estimation on frequency-corrected signal.
% Reads 1e6 samples from outFreqCorr.bin starting at user-specified time (sec).
% Same display as plot_outfft: abscissa = offset frequency (Hz), magnitude.

Fs = 21.42e6;   % Sampling frequency (Hz)
NumSamples = 1e6;  % complex samples to read

time_sec = input('Start time (seconds, relative to run start): ');
if isempty(time_sec)
    time_sec = 0;
end

fid = fopen('/home/agueguen/Documents/Projets/Comtech_Viterbi3/data/outFreqCorr.bin', 'rb');
if fid == -1
    fid = fopen('../data/outFreqCorr.bin', 'rb');
end
if fid == -1
    error('Cannot open outFreqCorr.bin');
end

% File layout: interleaved I,Q (I0,Q0,I1,Q1,...), 4 bytes per float
% Skip to start: time_sec * Fs complex samples = time_sec*Fs * 2 floats
offset_complex = round(time_sec * Fs);
offset_floats = offset_complex * 2;
offset_bytes = offset_floats * 4;
fseek(fid, offset_bytes, 'bof');

% Read 1e6 complex = 2e6 floats
d = fread(fid, NumSamples * 2, 'float32');
fclose(fid);

if numel(d) < NumSamples * 2
    warning('Requested %g samples, got %g. File may be short.', NumSamples, numel(d)/2);
    NumSamples = floor(numel(d) / 2);
    d = d(1:NumSamples*2);
end

I = d(1:2:end);
Q = d(2:2:end);
z = I + 1i*Q;

% z^4 (QPSK 4th-power), then FFT
z4 = z.^4;
fft_out = fft(z4);
N = numel(fft_out);

% FFT bin -> frequency in z^4 spectrum (Hz), then offset = f_fft/4
f_fft_Hz = (0:N-1)' * (Fs / N);
f_fft_Hz(f_fft_Hz > Fs/2) = f_fft_Hz(f_fft_Hz > Fs/2) - Fs;
f_offset_Hz = f_fft_Hz / 4;

% Peak (skip DC if desired; here include 0)
[~, peak_idx] = max(abs(fft_out));
peak_offset_Hz = f_offset_Hz(peak_idx);
fprintf('Peak at offset frequency: %g Hz\n', peak_offset_Hz);

figure('Name', 'Freq-corrected signal FFT');
subplot(2,1,1);
plot(f_offset_Hz, abs(fft_out));
xlabel('Offset frequency (Hz)');
ylabel('|FFT(z^4)|');
title(sprintf('Magnitude vs offset (t=%.2f s, %g samples)', time_sec, NumSamples));
grid on;

subplot(2,1,2);
plot(f_offset_Hz, 20*log10(max(abs(fft_out), 1e-12)));
xlabel('Offset frequency (Hz)');
ylabel('|FFT(z^4)| (dB)');
title('Magnitude (dB) vs offset frequency');
grid on;
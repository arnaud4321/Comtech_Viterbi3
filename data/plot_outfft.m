% Plot FFT from outfft.bin with abscissa = offset frequency (Hz).
% FFT is on z^4 = (I+j*Q)^4, so tone is at 4*f_offset -> f_offset = f_fft/4.
% Sampling frequency and power-4 are taken into account.

Fs = 21.42e6;   % Sampling frequency (Hz)

fid = fopen('/home/agueguen/Documents/Projets/Comtech_Viterbi3/data/outfft.bin', 'rb');
if fid == -1
    fid = fopen('../data/outfft.bin', 'rb');
end
if fid == -1
    error('Cannot open outfft.bin');
end
d = fread(fid, Inf, 'float32');
fclose(fid);

I = d(1:2:end);
Q = d(2:2:end);
fft_out = I + 1i*Q;
N = numel(fft_out);

% FFT bin k -> frequency in z^4 spectrum (Hz): [-Fs/2, Fs/2)
% bin 0..N/2-1 -> 0 to Fs/2; bin N/2..N-1 -> -Fs/2 to 0
f_fft_Hz = (0:N-1)' * (Fs / N);
f_fft_Hz(f_fft_Hz > Fs/2) = f_fft_Hz(f_fft_Hz > Fs/2) - Fs;

% Offset frequency (carrier): z^4 has tone at 4*f_offset => f_offset = f_fft/4
f_offset_Hz = f_fft_Hz / 4;

figure('Name', 'FFT estimation');
subplot(2,1,1);
plot(f_offset_Hz, abs(fft_out));
xlabel('Offset frequency (Hz)');
ylabel('|FFT(z^4)|');
title('Magnitude vs offset frequency (power 4)');
grid on;

subplot(2,1,2);
plot(f_offset_Hz, 20*log10(max(abs(fft_out), 1e-12)));
xlabel('Offset frequency (Hz)');
ylabel('|FFT(z^4)| (dB)');
title('Magnitude (dB) vs offset frequency');
grid on;

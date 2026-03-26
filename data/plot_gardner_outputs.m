% plot_gardner_outputs.m
% Offline viewer for Gardner debug dumps:
%   - gardner_early.bin
%   - gardner_ontime.bin
%   - gardner_late.bin
%   - gardner_err.bin
%   - gardner_omega.bin
%
% Files contain interleaved float32 I,Q complex samples.
% User provides start time (seconds) and number of monitored samples.

clear; close all;

Fs_in = 21.42e6;            % input sampling frequency (2 sps)
Fs_out = Fs_in / 2;         % Gardner outputs one sample per symbol

start_sec=0.0;
num_samples=10000;

if (num_samples==0)
    start_sec = input('Start time [s] (wrt Fs_in): ');
    num_samples = input('Number of monitored samples: ');
end

if isempty(start_sec) || isempty(num_samples) || num_samples <= 0
    error('Invalid user inputs.');
end

start_idx = max(0, round(start_sec * Fs_out));

fn_early = 'gardner_early.bin';
fn_ontime = 'gardner_ontime.bin';
fn_late = 'gardner_late.bin';
fn_err = 'gardner_err.bin';
fn_omega = 'gardner_omega.bin';

% Read early
fidE = fopen(fn_early, 'rb');
if fidE == -1
    error('Cannot open file: %s', fn_early);
end
fseek(fidE, start_idx * 2 * 4, 'bof');
rawE = fread(fidE, [2, num_samples], 'float32=>double');
fclose(fidE);
zE = complex(rawE(1,:), rawE(2,:));

% Read on-time
fidO = fopen(fn_ontime, 'rb');
if fidO == -1
    error('Cannot open file: %s', fn_ontime);
end
fseek(fidO, start_idx * 2 * 4, 'bof');
rawO = fread(fidO, [2, num_samples], 'float32=>double');
fclose(fidO);
zO = complex(rawO(1,:), rawO(2,:));

% Read late
fidL = fopen(fn_late, 'rb');
if fidL == -1
    error('Cannot open file: %s', fn_late);
end
fseek(fidL, start_idx * 2 * 4, 'bof');
rawL = fread(fidL, [2, num_samples], 'float32=>double');
fclose(fidL);
zL = complex(rawL(1,:), rawL(2,:));

% Read err (float32, 1 per symbol)
fidErr = fopen(fn_err, 'rb');
if fidErr == -1
    warning('Cannot open file: %s (skipping err plot)', fn_err);
    err = [];
else
    fseek(fidErr, start_idx * 4, 'bof');
    err = fread(fidErr, [1, num_samples], 'float32=>double');
    fclose(fidErr);
end

% Read omega (float32, 1 per symbol)
fidOmega = fopen(fn_omega, 'rb');
if fidOmega == -1
    warning('Cannot open file: %s (skipping omega plot)', fn_omega);
    omega = [];
else
    fseek(fidOmega, start_idx * 4, 'bof');
    omega = fread(fidOmega, [1, num_samples], 'float32=>double');
    fclose(fidOmega);
end

n = min([numel(zE), numel(zO), numel(zL)]);
if n == 0
    error('No data available in requested interval.');
end
zE = zE(1:n);
zO = zO(1:n);
zL = zL(1:n);
t = (0:n-1) / Fs_out;
if ~isempty(err)
    err = err(1:min(n, numel(err)));
end
if ~isempty(omega)
    omega = omega(1:min(n, numel(omega)));
end

fprintf('Loaded %d samples starting at %.6f s (Fs_out = %.3f MHz)\n', n, start_sec, Fs_out/1e6);

figure('Name', 'Gardner constellation');
subplot(1,3,1);
plot(real(zE), imag(zE), '.'); axis equal; grid on;
title('Early constellation'); xlabel('I'); ylabel('Q');
subplot(1,3,2);
plot(real(zO), imag(zO), '.'); axis equal; grid on;
title('On-time constellation'); xlabel('I'); ylabel('Q');
subplot(1,3,3);
plot(real(zL), imag(zL), '.'); axis equal; grid on;
title('Late constellation'); xlabel('I'); ylabel('Q');

figure('Name', 'Gardner amplitudes');
plot(t, abs(zE), 'DisplayName', 'Early'); hold on;
plot(t, abs(zO), 'DisplayName', 'On-time');
plot(t, abs(zL), 'DisplayName', 'Late');
grid on; xlabel('Time [s]'); ylabel('|z|');
title('Amplitude evolution');
legend('Location', 'best');

err2 = real(zO).*real(zL-zE) + imag(zO).*imag(zL-zE);


if ~isempty(err)
    figure('Name', 'Gardner timing error');
    plot(t(1:numel(err)), err, 'DisplayName', 'err'); grid on;
    xlabel('Time [s]'); ylabel('Timing error');
    title('Gardner error evolution');
    hold on;
    plot(t(1:numel(err2)), err2, 'DisplayName', 'err2'); grid on;
    xlabel('Time [s]'); ylabel('Timing error2');
end

if ~isempty(omega)
    figure('Name', 'Gardner omega');
    plot(t(1:numel(omega)), omega, 'DisplayName', 'omega'); grid on;
    xlabel('Time [s]'); ylabel('omega [samples/symbol]');
    title('Gardner omega evolution');
end


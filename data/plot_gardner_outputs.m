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

start_sec=0.5;
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
fn_phi = 'gardner_phi.bin';
fn_vit_in = 'viterbi_input_from_gardner.bin';
fn_gin = 'gardner_input_iq.bin';

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
    fseek(fidErr, start_idx * 8, 'bof');
    err = fread(fidErr, [1, num_samples], 'double');
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

% Read phi (float32, 1 per symbol): phase within symbol, in [0,1)
fidPhi = fopen(fn_phi, 'rb');
if fidPhi == -1
    warning('Cannot open file: %s (skipping phi plot)', fn_phi);
    phi = [];
else
    fseek(fidPhi, start_idx * 4, 'bof');
    phi = fread(fidPhi, [1, num_samples], 'float32=>double');
    fclose(fidPhi);
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
if ~isempty(phi)
    phi = phi(1:min(n, numel(phi)));
end

% Read Gardner input IQ stream at Fs_in (float32 IQ interleaved)
start_idx_in = max(0, round(start_sec * Fs_in));
fidGin = fopen(fn_gin, 'rb');
if fidGin == -1
    warning('Cannot open file: %s (skipping Gardner input eye plot)', fn_gin);
    zGin = [];
else
    fseek(fidGin, start_idx_in * 2 * 4, 'bof');
    rawGin = fread(fidGin, [2, 2*num_samples], 'float32=>double');
    fclose(fidGin);
    zGin = complex(rawGin(1,:), rawGin(2,:));
end

% Read symbols actually sent to Viterbi (float32 IQ interleaved)
fidVit = fopen(fn_vit_in, 'rb');
if fidVit == -1
    warning('Cannot open file: %s (skipping Viterbi input plot)', fn_vit_in);
    zVit = [];
else
    fseek(fidVit, start_idx * 2 * 4, 'bof');
    rawVit = fread(fidVit, [2, num_samples], 'float32=>double');
    fclose(fidVit);
    zVit = complex(rawVit(1,:), rawVit(2,:));
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

fprintf('mean err=%3.3f, mean err2=%3.3f\n',mean(err), mean(err2));

w = 2000;                 % fenêtre en symboles
b = ones(1,w)/w;
err_ma = filter(b, 1, err);
figure; plot(err_ma); grid on; title('moving mean(err) via filter');

if ~isempty(omega)
    figure('Name', 'Gardner omega');
    plot(t(1:numel(omega)), omega, 'DisplayName', 'omega'); grid on;
    xlabel('Time [s]'); ylabel('omega [samples/symbol]');
    title('Gardner omega evolution');
    
    h = zoom;
    % On définit une fonction qui s'exécute APRÈS chaque zoom
    set(h, 'ActionPostCallback', @(obj, evd) ...
        set(evd.Axes, 'YTickLabel', num2str(get(evd.Axes, 'YTick')', '%.6f')));
    
    % On le fait aussi pour le Pan (la main pour se déplacer)
    hPan = pan;
    set(hPan, 'ActionPostCallback', @(obj, evd) ...
        set(evd.Axes, 'YTickLabel', num2str(get(evd.Axes, 'YTick')', '%.6f')));
    
    % On force le premier affichage
    set(gca, 'YTickLabel', num2str(get(gca, 'YTick')', '%.6f'));
end

if ~isempty(phi)
    figure('Name', 'Gardner phase within symbol (phi)');
    plot(t(1:numel(phi)), phi, 'DisplayName', 'phi'); grid on;
    xlabel('Time [s]'); ylabel('phi in [0,1)');
    title('Gardner phase within symbol: phi = frac(t/2)');

    if ~isempty(err)
        figure('Name', 'err vs phi');
        scatter(phi(1:numel(err)), err, 2, '.'); grid on;
        xlabel('phi in [0,1)'); ylabel('err');
        title('Gardner timing error vs phase within symbol');
    end
end

if ~isempty(zVit)
    nVit = numel(zVit);
    figure('Name', 'Viterbi input from Gardner');
    subplot(1,2,1);
    plot(real(zVit), imag(zVit), '.'); axis equal; grid on;
    title('Symbols sent to Viterbi'); xlabel('I'); ylabel('Q');
    subplot(1,2,2);
    nCmp = min(n, nVit);
    plot(real(zO(1:nCmp)), imag(zO(1:nCmp)), '.'); axis equal; grid on;
    title('Gardner on-time (reference)'); xlabel('I'); ylabel('Q');
end

if ~isempty(zGin)
    % Eye diagram at Gardner input (2 sps): overlap 2 symbols windows.
    sps = 2;
    win = 2 * sps; % 2 symbols
    nTrace = min(2000, floor((numel(zGin) - win) / sps));
    if nTrace > 0
        tEye = (0:win-1) / sps;
        figure('Name', 'Gardner input eye diagram');
        subplot(1,2,1); hold on;
        for k = 0:nTrace-1
            idx = 1 + k * sps;
            plot(tEye, real(zGin(idx:idx+win-1)), 'b-');
        end
        grid on; xlabel('Symbol time'); ylabel('I');
        title('Eye diagram at Gardner input (I)');

        subplot(1,2,2); hold on;
        for k = 0:nTrace-1
            idx = 1 + k * sps;
            plot(tEye, imag(zGin(idx:idx+win-1)), 'r-');
        end
        grid on; xlabel('Symbol time'); ylabel('Q');
        title('Eye diagram at Gardner input (Q)');
    end
end


% Plot DD phase tracking debug dumps (ReceiverPhaseTrackingDD, DEBUG_PHASE_DD).
%
% Files (float32):
%   - phase_dd_in.bin    : interleaved I,Q per symbol (input symbols to DD loop)
%   - phase_dd_out.bin   : interleaved I,Q per symbol (output after phase correction)
%   - phase_dd_phase.bin : phase estimate (rad) per symbol
%   - phase_dd_err.bin   : phase error (rad) per symbol
%   - phase_dd_err_filt.bin : filtered error magnitude used for lock: EMA(|err|) (rad), per symbol
%   - phase_dd_thresh.bin   : lock threshold (rad), per symbol
%
% This script plots the input/output constellations and the phase/phase-error evolution
% starting at a user-specified time (sec) for N symbols.
%
% NOTE: time is referenced to symbol rate (1 sample per symbol at this point).

SymbolRate = 10.71e6; % symbols per second (QPSK)

start_sec = input('Start time (seconds, wrt symbol clock): ');
if isempty(start_sec), start_sec = 0; end

Nsym = input('Number of symbols to plot: ');
if isempty(Nsym), Nsym = 2000; end

start_sym = round(start_sec * SymbolRate);
offset_iq_bytes = start_sym * 2 * 4;   % 2 floats per symbol, 4 bytes per float
offset_f_bytes  = start_sym * 4;       % 1 float per symbol

% Open in current dir, fallback to ../data (no local function to keep compatibility)
fidIn = fopen('phase_dd_in.bin', 'rb');
if fidIn == -1, fidIn = fopen(fullfile('..','data','phase_dd_in.bin'), 'rb'); end

fidOut = fopen('phase_dd_out.bin', 'rb');
if fidOut == -1, fidOut = fopen(fullfile('..','data','phase_dd_out.bin'), 'rb'); end

fidPh = fopen('phase_dd_phase.bin', 'rb');
if fidPh == -1, fidPh = fopen(fullfile('..','data','phase_dd_phase.bin'), 'rb'); end

fidEr = fopen('phase_dd_err.bin', 'rb');
if fidEr == -1, fidEr = fopen(fullfile('..','data','phase_dd_err.bin'), 'rb'); end

fidErF = fopen('phase_dd_err_filt.bin', 'rb');
if fidErF == -1, fidErF = fopen(fullfile('..','data','phase_dd_err_filt.bin'), 'rb'); end

fidTh = fopen('phase_dd_thresh.bin', 'rb');
if fidTh == -1, fidTh = fopen(fullfile('..','data','phase_dd_thresh.bin'), 'rb'); end

if fidIn == -1 || fidOut == -1
    error('Cannot open phase_dd_in.bin or phase_dd_out.bin (enable DEBUG_PHASE_DD and run).');
end

fseek(fidIn, offset_iq_bytes, 'bof');
din = fread(fidIn, 2*Nsym, 'float32');
fclose(fidIn);
Iin = din(1:2:end); Qin = din(2:2:end);

fseek(fidOut, offset_iq_bytes, 'bof');
dout = fread(fidOut, 2*Nsym, 'float32');
fclose(fidOut);
Iout = dout(1:2:end); Qout = dout(2:2:end);

phase = [];
perr = [];
err_filt = [];
thresh = [];
if fidPh ~= -1
    fseek(fidPh, offset_f_bytes, 'bof');
    phase = fread(fidPh, Nsym, 'float32');
    fclose(fidPh);
end
if fidEr ~= -1
    fseek(fidEr, offset_f_bytes, 'bof');
    perr = fread(fidEr, Nsym, 'float32');
    fclose(fidEr);
end
if fidErF ~= -1
    fseek(fidErF, offset_f_bytes, 'bof');
    err_filt = fread(fidErF, Nsym, 'float32');
    fclose(fidErF);
end
if fidTh ~= -1
    fseek(fidTh, offset_f_bytes, 'bof');
    thresh = fread(fidTh, Nsym, 'float32');
    fclose(fidTh);
end

t = (start_sym + (0:Nsym-1))' / SymbolRate;

figure('Name','DD Phase Tracking Debug');

subplot(2,2,1);
plot(Iin, Qin, '.', 'MarkerSize', 6);
axis equal; grid on;
xlabel('I'); ylabel('Q');
title(sprintf('Input constellation (t=%.4f s, N=%d)', start_sec, Nsym));

subplot(2,2,2);
plot(Iout, Qout, '.', 'MarkerSize', 6);
axis equal; grid on;
xlabel('I'); ylabel('Q');
title('Output constellation (phase-corrected)');

subplot(2,2,3);
if ~isempty(phase)
    plot(t, phase, 'r-'); grid on;
    xlabel('Time (s)'); ylabel('Phase est (rad)');
    title('Phase estimate');
else
    text(0.5,0.5,'phase_dd_phase.bin not available','HorizontalAlignment','center');
end

subplot(2,2,4);
if ~isempty(perr)
    plot(t, perr, 'b-'); hold on; grid on;
    if ~isempty(err_filt)
        plot(t, err_filt, 'm-');
    end
    if ~isempty(thresh)
        plot(t, thresh, 'k--');
    end
    hold off;
    xlabel('Time (s)'); ylabel('Phase error (rad)');
    title('DD phase error + filtered error + threshold');
    legend_entries = {'err (rad)'};
    if ~isempty(err_filt), legend_entries{end+1} = 'err\_filt (EMA|err|)'; end
    if ~isempty(thresh), legend_entries{end+1} = 'threshold'; end
    legend(legend_entries, 'Location', 'best');
else
    text(0.5,0.5,'phase_dd_err.bin not available','HorizontalAlignment','center');
end


import sys
import os
import time
import math
import queue
import threading
import traceback
from pathlib import Path


def write_troubleshoot(msg: str):
    try:
        repo = Path(__file__).resolve().parent.parent
        data = repo / "data"
        data.mkdir(parents=True, exist_ok=True)
        p = data / "constellation_display.troubleshoot.log"
        with open(p, "wt") as f:
            f.write(msg)
            if not msg.endswith("\n"):
                f.write("\n")
    except Exception:
        pass


try:
    import matplotlib
    matplotlib.use("TkAgg")
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    import numpy as np
except Exception:
    write_troubleshoot(
        "Constellation display failed to start (matplotlib/TkAgg import).\n"
        f"Python: {sys.executable}\n"
        f"DISPLAY={os.environ.get('DISPLAY','')}\n\n"
        + traceback.format_exc()
        + "\n\n"
        "A troubleshooting log was created. Please send this file to the integrator:\n"
        "  data/constellation_display.troubleshoot.log\n\n"
        "Client steps:\n"
        "- Ensure an X11 display is available (DISPLAY set, e.g. ':0').\n"
        "- Install python3-tk and matplotlib:\n"
        "    sudo apt install -y python3-tk python3-matplotlib\n"
        "- If running over SSH: use X11 forwarding (ssh -X) or run locally.\n"
    )
    raise


def parse_kv(line: str):
    out = {}
    parts = line.strip().split()
    for p in parts[1:]:
        if "=" in p:
            k, v = p.split("=", 1)
            out[k] = v
    return out


def fmt_lock_status(raw: str) -> str:
    """Map 1/0 from C++ to readable English for lock flags."""
    s = str(raw).strip()
    if s == "1":
        return "locked"
    if s == "0":
        return "unlocked"
    return s


def fmt_ready_status(raw: str) -> str:
    """centralReady: 1/0 -> ready / not_ready."""
    s = str(raw).strip()
    if s == "1":
        return "ready"
    if s == "0":
        return "not_ready"
    return s


def format_sync_event_line(kv: dict) -> str:
    """Une ligne lisible pour le panneau (kind + horodatage + détails utiles)."""
    kind = kv.get("kind", "?")
    tw = kv.get("t_wall_iso", "")
    ts = kv.get("t_sim", "")
    parts = [tw, kind, f"t_sim={ts}s"]
    if kind == "viterbi_desync":
        parts.append(f"ema={kv.get('emaGrowth', '?')}")
        parts.append(f"thr={kv.get('thr', '?')}")
        parts.append(f"gard={kv.get('gardLocked', '?')} phase={kv.get('phaseLocked', '?')}")
        parts.append(f"evm={kv.get('evmRms', '?')}")
        parts.append(f"gPpm={kv.get('gardnerPpm', '?')}")
        if "hints" in kv:
            parts.append(f"hints={kv['hints']}")
    elif kind == "viterbi_lock":
        parts.append(f"minM=({kv.get('minMargin0', '?')},{kv.get('minMargin1', '?')},{kv.get('minMargin2', '?')})")
    elif kind == "viterbi_acq_timeout":
        parts.append(f"gard={kv.get('gardLocked', '?')} phase={kv.get('phaseLocked', '?')}")
    elif kind == "prbs_desync":
        parts.append(f"berBatch={kv.get('berBatch', '?')}")
        if "hints" in kv:
            parts.append(f"hints={kv['hints']}")
    elif kind == "prbs_lock":
        parts.append(f"seed={kv.get('seed', '?')}")
    line = "  ".join(str(p) for p in parts)
    # Largeur affichage (une ligne par événement) — évite les débordements horizontaux.
    if len(line) > 132:
        return line[:129] + "..."
    return line


max_abs = 2.0
draw_period = 0.20  # seconds (max GUI refresh rate), may be overridden by CONFIG

plt.ion()
fig = plt.figure(num="Constellation (Viterbi input)", constrained_layout=False, facecolor="#f5f5f5")
try:
    # Hauteur généreuse : les 2 rangées de courbes ne doivent pas être mangées par la barre d’outils.
    fig.set_size_inches(15.0, 11.5, forward=True)
except Exception:
    pass

# bottom ~11% : réserve pour la toolbar matplotlib (sinon les courbes du bas sont tronquées).
fig.subplots_adjust(left=0.048, right=0.985, top=0.93, bottom=0.11, wspace=0.20, hspace=0.42)

# Ratios : moins de poids sur la rangée 0 → PPM / SNR / BER restent lisibles.
gs = fig.add_gridspec(
    nrows=3,
    ncols=3,
    height_ratios=[2.85, 1.22, 1.22],
    width_ratios=[1.85, 2.25, 3.4],
)
ax = fig.add_subplot(gs[0, 0])
ax_spec = fig.add_subplot(gs[0, 1])
# Colonne infos : Channel un peu plus court, sync plus compact (moins de lignes coupées par clip).
gs_info = gs[0, 2].subgridspec(3, 1, height_ratios=[0.88, 1.05, 0.58], hspace=0.36)
ax_info_left = fig.add_subplot(gs_info[0, 0])
ax_info_right = fig.add_subplot(gs_info[1, 0])
ax_info_sync = fig.add_subplot(gs_info[2, 0])
for _a in (ax_info_left, ax_info_right, ax_info_sync):
    _a.set_axis_off()
    _a.set_facecolor("#fafafa")

hist_spec1 = gs[1, :].subgridspec(1, 3, wspace=0.20)
ax_ppm = fig.add_subplot(hist_spec1[0, 0])
ax_gain = fig.add_subplot(hist_spec1[0, 1], sharex=ax_ppm)
ax_freq = fig.add_subplot(hist_spec1[0, 2], sharex=ax_ppm)

hist_spec2 = gs[2, :].subgridspec(1, 3, wspace=0.20)
ax_snr = fig.add_subplot(hist_spec2[0, 0], sharex=ax_ppm)
ax_ber = fig.add_subplot(hist_spec2[0, 1], sharex=ax_ppm)
ax_lock = fig.add_subplot(hist_spec2[0, 2], sharex=ax_ppm)

sc = ax.scatter([], [], s=4)
ax.set_title("Constellation (Viterbi input)", pad=6)
ax.set_xlabel("I")
ax.set_ylabel("Q")
ax.grid(True, alpha=0.3)
ax.set_aspect("equal", adjustable="box")

ax_spec.set_title("RX Spectrum (Baseband)", pad=6)
ax_spec.set_xlabel("Normalized Frequency", labelpad=6)
ax_spec.set_ylabel("Magnitude (dB)")
ax_spec.grid(True, alpha=0.3)
ax_spec.set_xlim(-0.5, 0.5)
ax_spec.set_ylim(-60, 20)
line_spec, = ax_spec.plot([], [], "-", lw=1.0, color='blue')

# clip_on=True : obligatoire — sinon les bbox dépassent et recouvrent le panneau du dessous (sync).
TEXT_KW = dict(
    family="monospace",
    va="top",
    clip_on=True,
)
info_text_left = ax_info_left.text(
    0.02,
    0.99,
    "",
    transform=ax_info_left.transAxes,
    ha="left",
    fontsize=6.75,
    linespacing=1.05,
    bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="#333333", linewidth=0.7, alpha=0.98),
    **TEXT_KW,
)
info_text_rx_l = ax_info_right.text(
    0.02,
    0.99,
    "",
    transform=ax_info_right.transAxes,
    ha="left",
    fontsize=6.5,
    linespacing=1.05,
    bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="#333333", linewidth=0.7, alpha=0.98),
    **TEXT_KW,
)
info_text_rx_r = ax_info_right.text(
    0.51,
    0.99,
    "",
    transform=ax_info_right.transAxes,
    ha="left",
    fontsize=6.5,
    linespacing=1.05,
    bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="#333333", linewidth=0.7, alpha=0.98),
    **TEXT_KW,
)
info_text_sync = ax_info_sync.text(
    0.02,
    0.99,
    "",
    transform=ax_info_sync.transAxes,
    ha="left",
    fontsize=6.0,
    linespacing=1.08,
    bbox=dict(boxstyle="round,pad=0.18", fc="#fff8f0", ec="#aa6633", linewidth=0.7, alpha=0.98),
    **TEXT_KW,
)


def apply_limits():
    ax.set_xlim(-max_abs, max_abs)
    ax.set_ylim(-max_abs, max_abs)


def raise_constellation_window():
    """TkAgg: la fenêtre peut rester derrière le terminal ou ne pas se mapper sans un lift explicite."""
    try:
        w = getattr(fig.canvas.manager, "window", None)
        if w is None:
            return
        try:
            # Taille mini pour que la grille (3 rangées + toolbar) ne soit pas écrasée verticalement.
            w.minsize(1280, 920)
        except Exception:
            pass
        try:
            w.lift()
        except Exception:
            pass
        try:
            w.focus_force()
        except Exception:
            pass
    except Exception:
        pass


apply_limits()
plt.show(block=False)
plt.pause(0.001)
raise_constellation_window()

# Quit only when the user closes the window — not when stdin closes (Ctrl+C on C++, clean shutdown, EOF).
_close_state = [False]


def _on_figure_close(_evt):
    _close_state[0] = True


fig.canvas.mpl_connect("close_event", _on_figure_close)
_feed_end_notified = [False]

ax_ppm.set_title("Total PPM: applied vs estimated", pad=8)
ax_ppm.set_ylabel("ppm")
ax_ppm.grid(True, alpha=0.3)
ax_ppm.tick_params(axis="both", labelsize=8)

ax_gain.set_title("Gain: applied vs corrected", pad=8)
ax_gain.set_ylabel("dB", labelpad=6)
ax_gain.grid(True, alpha=0.3)
ax_gain.tick_params(axis="both", labelsize=8)

ax_freq.set_title("Frequency: applied vs corrected", pad=8)
ax_freq.set_ylabel("Hz", labelpad=14)
ax_freq.grid(True, alpha=0.3)
ax_freq.tick_params(axis="y", pad=8)
ax_freq.yaxis.set_label_coords(-0.07, 0.5)

ax_snr.set_title("SNR (dB)", pad=6)
ax_snr.set_ylabel("dB", labelpad=6)
ax_snr.grid(True, alpha=0.3)
ax_snr.set_xlabel("t_sim (s)")
ax_snr.tick_params(axis="both", labelsize=8)

ax_ber.set_title("BER", pad=6)
ax_ber.set_ylabel("BER", labelpad=6)
ax_ber.set_yscale("log")
ax_ber.grid(True, alpha=0.3)
ax_ber.set_xlabel("t_sim (s)")
ax_ber.tick_params(axis="both", labelsize=8)

ax_lock.set_title("Lock Status", pad=6)
ax_lock.set_yticks([0, 1])
ax_lock.set_yticklabels(['Unlock', 'Lock'])
ax_lock.grid(True, alpha=0.3)
ax_lock.set_xlabel("t_sim (s)")
ax_lock.tick_params(axis="both", labelsize=8)

LEG_FS = 7
LEG_KW = dict(fontsize=LEG_FS, framealpha=0.92, borderpad=0.4)

line_ppm_applied, = ax_ppm.plot([], [], "-", lw=1.2, label="Channel totalPpm")
line_ppm_sum, = ax_ppm.plot([], [], "--", lw=1.2, label="Rx totalPpm corrected")
ax_ppm.legend(loc="lower left", **LEG_KW)

line_gain_applied, = ax_gain.plot([], [], "-", lw=1.2, label="Channel gainDb")
line_gain_corr, = ax_gain.plot([], [], "--", lw=1.2, label="Rx gainDb")
ax_gain.legend(loc="lower left", **LEG_KW)

line_f_applied, = ax_freq.plot([], [], "-", lw=1.2, label="Channel freqHz")
line_f_corr, = ax_freq.plot([], [], "--", lw=1.2, label="Rx (NCO+PhaseDD) Hz")
ax_freq.legend(loc="upper left", **LEG_KW)

line_snr_applied, = ax_snr.plot([], [], "-", lw=1.2, label="Channel SNR")
line_snr_est, = ax_snr.plot([], [], "--", lw=1.2, label="Est SNR (Viterbi)")
ax_snr.legend(loc="lower left", **LEG_KW)

line_ber, = ax_ber.plot([], [], "-", lw=1.2, color="red", label="BER (PRBS)")
line_ber_vit, = ax_ber.plot([], [], "--", lw=1.2, color="darkorange", label="raw BER (Vit)")
ax_ber.legend(loc="lower left", **LEG_KW)

line_lock_vit, = ax_lock.step([], [], "-", lw=1.2, label="Viterbi", where="post")
line_lock_prbs, = ax_lock.step([], [], "--", lw=1.2, label="PRBS", where="post")
ax_lock.legend(loc="lower left", **LEG_KW)

hist_max = 3000
# Spectrum FFT: use at most this many samples (subsample the frame) — full N≈2048 FFT each draw is slow in Tk.
spec_fft_max = 1024
t_hist = []
ppm_applied = []
ppm_sum = []
ppm_sum_ema = None
spec_pwr_ema = None
gain_applied = []
gain_corr = []
f_applied = []
f_corr = []
snr_applied = []
snr_est = []
ber_hist = []
ber_vit_hist = []
lock_vit = []
lock_prbs = []
sync_event_log = []
sync_event_max = 6

stdin_closed = False
latest_pts = None
latest_kv = None
dirty = False
last_draw = 0.0
last_ylim_update = 0.0
x_window_sec = 180.0

# Événements stdin → file (découplage Tk/matplotlib). Le thread lit en continu pour ne pas
# mélanger select/readline TextIOWrapper avec la boucle GUI (souvent: aucune donnée affichée).
evt_q = queue.Queue()


def _decode_line(raw: bytes) -> str:
    return raw.decode("ascii", errors="replace").strip()


def _stdin_reader():
    """Lit le protocole C++ sur stdin (binaire) sans partager le GIL avec les redraws."""
    bio = sys.stdin.buffer
    while True:
        try:
            line = bio.readline()
        except Exception:
            evt_q.put(("EOF",))
            return
        if not line:
            evt_q.put(("EOF",))
            return
        s = _decode_line(line)
        if not s:
            continue
        if s.startswith("QUIT"):
            evt_q.put(("QUIT",))
            return
        if s.startswith("CONFIG"):
            evt_q.put(("CONFIG", s))
            continue
        if s.startswith("EVENT"):
            kv = parse_kv(s)
            evt_q.put(("EVENT", kv))
            continue
        if s.startswith("FRAME"):
            kv = parse_kv(s)
            pts_x = []
            pts_y = []
            while True:
                try:
                    l2b = bio.readline()
                except Exception:
                    evt_q.put(("EOF",))
                    return
                if not l2b:
                    evt_q.put(("EOF",))
                    return
                l2 = _decode_line(l2b)
                if l2 == "END":
                    break
                try:
                    x, y = l2.split()
                    pts_x.append(float(x))
                    pts_y.append(float(y))
                except Exception:
                    continue
            evt_q.put(("FRAME", kv, pts_x, pts_y))
            continue


threading.Thread(target=_stdin_reader, name="constellation_stdin", daemon=True).start()


def _apply_frame_history(kv):
    global ppm_sum_ema
    try:
        if "t_sim" not in kv:
            return
        t = float(kv["t_sim"])
        t_hist.append(t)
        ppm_applied.append(float(kv.get("chanTotalPpm", "nan")))
        vppm = float(kv.get("rxTotalPpmCorr", "nan"))
        if vppm == vppm:
            if ppm_sum_ema is None:
                ppm_sum_ema = vppm
            else:
                ppm_sum_ema = 0.15 * vppm + 0.85 * ppm_sum_ema
        ppm_sum.append(ppm_sum_ema if ppm_sum_ema is not None else vppm)
        gain_applied.append(float(kv.get("chanGainDb", kv.get("chanRangeDb", "nan"))))
        gain_corr.append(float(kv.get("gainDb", "nan")))
        f_a = float(kv.get("chanFreqHz", "nan"))
        f_c = float(kv.get("centralNcoHz", "nan")) + float(kv.get("phaseFreqHz", "0.0"))
        f_applied.append(f_a)
        f_corr.append(f_c)

        snr_applied.append(float(kv.get("chanSnrDb", "nan")))
        snr_est.append(float(kv.get("viterbiSurvivorRawCodedSnrDb", "nan")))

        ber_v = float(kv.get("ber", "nan"))
        if ber_v == 0.0:
            ber_v = 1e-9
        ber_hist.append(ber_v)

        ber_vit_v = float(kv.get("viterbiSurvivorRawCodedBer", "nan"))
        if ber_vit_v == ber_vit_v and ber_vit_v == 0.0:
            ber_vit_v = 1e-9
        ber_vit_hist.append(ber_vit_v)

        vit = str(kv.get("vitLocked", "0")).strip()
        prbs = str(kv.get("prbsLocked", "0")).strip()
        lock_vit.append(0.95 if vit == "1" else 0.0)
        lock_prbs.append(1.0 if prbs == "1" else 0.05)

        if len(t_hist) > hist_max:
            t_hist[:] = t_hist[-hist_max:]
            ppm_applied[:] = ppm_applied[-hist_max:]
            ppm_sum[:] = ppm_sum[-hist_max:]
            gain_applied[:] = gain_applied[-hist_max:]
            gain_corr[:] = gain_corr[-hist_max:]
            f_applied[:] = f_applied[-hist_max:]
            f_corr[:] = f_corr[-hist_max:]
            snr_applied[:] = snr_applied[-hist_max:]
            snr_est[:] = snr_est[-hist_max:]
            ber_hist[:] = ber_hist[-hist_max:]
            ber_vit_hist[:] = ber_vit_hist[-hist_max:]
            lock_vit[:] = lock_vit[-hist_max:]
            lock_prbs[:] = lock_prbs[-hist_max:]
    except Exception:
        pass


while not _close_state[0]:
    # Laisser Tk traiter les événements même sans nouvelles données
    plt.pause(0.02)

    if stdin_closed and not _feed_end_notified[0]:
        _feed_end_notified[0] = True
        try:
            print(
                "[constellation_display] Data feed ended (sim stopped or pipe closed). "
                "Close the plot window to exit.",
                file=sys.stderr,
            )
        except Exception:
            pass

    try:
        while True:
            evt = evt_q.get_nowait()
            if evt[0] == "EOF":
                stdin_closed = True
                break
            if evt[0] == "QUIT":
                stdin_closed = True
                break
            if evt[0] == "EVENT":
                _, kv = evt
                sync_event_log.append(format_sync_event_line(kv))
                if len(sync_event_log) > sync_event_max:
                    sync_event_log[:] = sync_event_log[-sync_event_max:]
                dirty = True
            elif evt[0] == "CONFIG":
                _, line = evt
                kv = parse_kv(line)
                if "maxAbs" in kv:
                    try:
                        max_abs = float(kv["maxAbs"])
                    except Exception:
                        pass
                if "drawPeriodSec" in kv:
                    try:
                        v = float(kv["drawPeriodSec"])
                        if v > 0.0:
                            draw_period = v
                    except Exception:
                        pass
                apply_limits()
                dirty = True
            elif evt[0] == "FRAME":
                _, kv, pts_x, pts_y = evt
                latest_pts = (pts_x, pts_y)
                latest_kv = kv
                dirty = True
                _apply_frame_history(kv)
    except queue.Empty:
        pass

    # Throttled rendering (keeps GUI responsive without blocking the simulation).
    now = time.time()
    if dirty and (now - last_draw) >= draw_period:
        dirty = False
        last_draw = now

        # Constellation points
        if latest_pts is not None:
            px, py = latest_pts
            sc.set_offsets(list(zip(px, py)))

            # Spectrum: drawn on ax_spec (middle top). Uses only the current FRAME symbols (not the whole run).
            # Subsample before FFT so redraw stays cheap; constellation above still uses every point.
            if len(px) > 0:
                c_pts = np.asarray(px, dtype=np.float64) + 1j * np.asarray(py, dtype=np.float64)
                n = len(c_pts)
                if n > spec_fft_max:
                    step = max(1, n // spec_fft_max)
                    c_spec = c_pts[::step][:spec_fft_max]
                else:
                    c_spec = c_pts
                win = np.hanning(len(c_spec))
                spec = np.fft.fftshift(np.fft.fft(c_spec * win))
                pwr = np.abs(spec)**2
                if spec_pwr_ema is None or len(spec_pwr_ema) != len(pwr):
                    spec_pwr_ema = pwr
                else:
                    alpha_spec = 0.4  # Increased from 0.1 to adapt faster (e.g. if 1 update/sec)
                    spec_pwr_ema = alpha_spec * pwr + (1.0 - alpha_spec) * spec_pwr_ema
                
                mag = 10 * np.log10(spec_pwr_ema + 1e-12)
                if len(mag) > 0:
                    mag -= np.max(mag)  # Normalize peak to 0 dB
                freqs = np.fft.fftshift(np.fft.fftfreq(len(c_spec), d=1.0))
                line_spec.set_data(freqs, mag)
                
                # Autoscale Y if needed, but usually -60 to 5 is fine.
                ax_spec.set_ylim(-60, 5)

        # Right panel text (English)
        if latest_kv is not None:
            kv = latest_kv
            seg_map = {"0": "stable1", "1": "acc1", "2": "stable2", "3": "acc2"}
            chan_lines = ["[Channel]"]
            if "t_sim" in kv and "t_rate" in kv:
                chan_lines.append(
                    f"t_sim={float(kv['t_sim']):.3f}s  t_rate={float(kv['t_rate']):.3f}s"
                )
            if "chanSeg" in kv:
                chan_lines.append(f"segment = {seg_map.get(kv['chanSeg'], kv['chanSeg'])}")
            if "chanFreqHz" in kv and "chanTotalPpm" in kv:
                chan_lines.append(f"freqHz={kv['chanFreqHz']}  totalPpm={kv['chanTotalPpm']}")
            elif "chanFreqHz" in kv:
                chan_lines.append(f"freqHz = {kv['chanFreqHz']}")
            elif "chanTotalPpm" in kv:
                chan_lines.append(f"totalPpm= {kv['chanTotalPpm']}")
            if "chanGainDb" in kv or "chanRangeDb" in kv:
                v = kv.get("chanGainDb", kv.get("chanRangeDb", ""))
                snr = ""
                if "chanSnrDb" in kv:
                    try:
                        sv = float(kv["chanSnrDb"])
                        if math.isfinite(sv):
                            snr = f"  SNR={sv:.2f}dB"
                    except Exception:
                        pass
                chan_lines.append(f"gainDb={v}{snr}")
            elif "chanSnrDb" in kv:
                try:
                    v = float(kv["chanSnrDb"])
                    if math.isfinite(v):
                        chan_lines.append(f"SNR     = {v:.2f} dB (applied)")
                except Exception:
                    pass
            if "samplerMsps" in kv:
                chan_lines.append(f"samplerMsps = {kv['samplerMsps']}")

            rx_lines = ["[Receiver]"]
            for key, label in [
                ("symRateMsps", "symRateMsps"),
                ("symRateCorrPpm", "symRateCorrPpm"),
                ("gardnerPpm", "gardnerPpm"),
                ("rxTotalPpmCorr", "rxTotalPpmCorr"),
                ("gardLocked", "gardLocked"),
                ("phaseLocked", "phaseLocked"),
                ("vitLocked", "vitLocked"),
                ("prbsLocked", "prbsLocked"),
                ("numBits", "numBits"),
                ("numErrors", "numErrors"),
                ("ber", "BER"),
                ("centralReady", "centralReady"),
                ("centralTargetHz", "centralTargetHz"),
                ("centralNcoHz", "centralNcoHz"),
                ("phaseFreqHz", "phaseFreqHz"),
                ("gainDb", "gainDb"),
                ("rxFilterMsps", "rxFilterMsps"),
            ]:
                if key in kv:
                    val = kv[key]
                    if key in ("gardLocked", "phaseLocked", "vitLocked", "prbsLocked"):
                        val = fmt_lock_status(val)
                    elif key == "centralReady":
                        val = fmt_ready_status(val)
                    rx_lines.append(f"{label:14s}= {val}")
            if "viterbiSurvivorRawCodedEvmRms" in kv and "viterbiSurvivorRawCodedSnrDb" in kv:
                try:
                    er = float(kv["viterbiSurvivorRawCodedEvmRms"])
                    snr_e = float(kv["viterbiSurvivorRawCodedSnrDb"])
                    if er == er and snr_e == snr_e:
                        rx_lines.append(f"{'EVM_rms (Vit)':14s}= {er:.6f} (frame)")
                        rx_lines.append(f"{'SNR_est (Vit)':14s}= {snr_e:.2f} dB")
                except Exception:
                    pass
            if "viterbiSurvivorRawCodedBer" in kv:
                try:
                    brv = float(kv["viterbiSurvivorRawCodedBer"])
                    if brv == brv:
                        rx_lines.append(f"{'rawBER (Vit)':14s}= {brv:.6e}")
                except Exception:
                    pass
            info_text_left.set_text("\n".join(chan_lines))
            rx_body = rx_lines[1:]
            if not rx_body:
                info_text_rx_l.set_text(rx_lines[0])
                info_text_rx_r.set_text("")
            else:
                k = (len(rx_body) + 1) // 2
                info_text_rx_l.set_text("\n".join([rx_lines[0]] + rx_body[:k]))
                info_text_rx_r.set_text("\n".join(rx_body[k:]))

        if sync_event_log:
            info_text_sync.set_text("[Sync / désync — horodaté]\n" + "\n".join(sync_event_log))
        else:
            info_text_sync.set_text("")

        # History lines (no relim/autoscale every frame; update limits at low rate)
        if t_hist:
            line_ppm_applied.set_data(t_hist, ppm_applied)
            line_ppm_sum.set_data(t_hist, ppm_sum)
            line_gain_applied.set_data(t_hist, gain_applied)
            line_gain_corr.set_data(t_hist, gain_corr)
            line_f_applied.set_data(t_hist, f_applied)
            line_f_corr.set_data(t_hist, f_corr)
            
            line_snr_applied.set_data(t_hist, snr_applied)
            line_snr_est.set_data(t_hist, snr_est)
            line_ber.set_data(t_hist, ber_hist)
            line_ber_vit.set_data(t_hist, ber_vit_hist)
            line_lock_vit.set_data(t_hist, lock_vit)
            line_lock_prbs.set_data(t_hist, lock_prbs)

            # X window
            tmax = t_hist[-1]
            xmin = max(0.0, tmax - x_window_sec)
            ax_ppm.set_xlim(xmin, tmax + 1e-6)
            ax_snr.set_xlim(xmin, tmax + 1e-6)

            # Y limits update only once per second
            if (now - last_ylim_update) >= 1.0:
                last_ylim_update = now
                def set_ylim(ax_, ys):
                    vals = [v for v in ys if v == v]  # drop NaN
                    if not vals:
                        return
                    lo = min(vals)
                    hi = max(vals)
                    if hi - lo < 1e-9:
                        hi = lo + 1.0
                    m = 0.10 * (hi - lo)
                    ax_.set_ylim(lo - m, hi + m)

                set_ylim(ax_ppm, ppm_applied + ppm_sum)
                set_ylim(ax_gain, gain_applied + gain_corr)
                set_ylim(ax_freq, f_applied + f_corr)
                set_ylim(ax_snr, snr_applied + snr_est)
                
                # BER
                vals_ber = [v for v in (ber_hist + ber_vit_hist) if v == v and v > 0]
                if vals_ber:
                    min_b = max(1e-8, min(vals_ber) * 0.5)
                    max_b = min(1.0, max(vals_ber) * 2.0)
                    if max_b <= min_b:
                        max_b = 1.0
                        min_b = 1e-6
                    ax_ber.set_ylim(min_b, max_b)
                
                ax_lock.set_ylim(-0.1, 1.1)

        apply_limits()
        fig.canvas.draw_idle()
        fig.canvas.flush_events()
    elif stdin_closed:
        # Pipe fermé : plus de nouvelles trames, mais on garde la dernière image et la boucle Tk vivante.
        fig.canvas.draw_idle()
        fig.canvas.flush_events()

import sys
import select
import os
import time
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


max_abs = 2.0
draw_period = 0.20  # seconds (max GUI refresh rate), may be overridden by CONFIG

plt.ion()
fig = plt.figure(num="Constellation (Viterbi input)", constrained_layout=False)
try:
    fig.set_size_inches(11.0, 7.0, forward=True)
except Exception:
    pass

fig.subplots_adjust(left=0.06, right=0.99, top=0.93, bottom=0.08, wspace=0.06, hspace=0.22)

# Layout:
# - Top: constellation (left) + state/status (right)
# - Bottom: 3 history plots side-by-side
gs = fig.add_gridspec(
    nrows=2,
    ncols=2,
    height_ratios=[3.4, 1.25],
    width_ratios=[4.0, 2.4],
)
ax = fig.add_subplot(gs[0, 0])
ax_info = fig.add_subplot(gs[0, 1])

hist_spec = gs[1, :].subgridspec(1, 3, wspace=0.20)
ax_ppm = fig.add_subplot(hist_spec[0, 0])
ax_gain = fig.add_subplot(hist_spec[0, 1], sharex=ax_ppm)
ax_freq = fig.add_subplot(hist_spec[0, 2], sharex=ax_ppm)

sc = ax.scatter([], [], s=4)
ax.set_title("Constellation (Viterbi input)")
ax.set_xlabel("I")
ax.set_ylabel("Q")
ax.grid(True, alpha=0.3)
ax.set_aspect("equal", adjustable="box")

ax_info.set_axis_off()
info_text = ax_info.text(
    0.02,
    0.98,
    "",
    transform=ax_info.transAxes,
    va="top",
    ha="left",
    fontsize=9,
    family="monospace",
    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", alpha=0.9),
)


def apply_limits():
    ax.set_xlim(-max_abs, max_abs)
    ax.set_ylim(-max_abs, max_abs)


apply_limits()
plt.show(block=False)
plt.pause(0.001)

ax_ppm.set_title("Total PPM: applied vs estimated")
ax_ppm.set_ylabel("ppm")
ax_ppm.grid(True, alpha=0.3)
ax_ppm.tick_params(axis="both", labelsize=8)

ax_gain.set_title("Gain: applied vs corrected")
ax_gain.set_ylabel("dB", labelpad=6)
ax_gain.grid(True, alpha=0.3)
ax_gain.tick_params(axis="both", labelsize=8)

ax_freq.set_title("Frequency: applied vs corrected")
ax_freq.set_ylabel("Hz", labelpad=4)
ax_freq.grid(True, alpha=0.3)
ax_freq.set_xlabel("t_sim (s)")
# Place the y-axis label more inside its own subplot to avoid overlapping ax_gain.
ax_freq.yaxis.set_label_coords(-0.05, 0.5)

line_ppm_applied, = ax_ppm.plot([], [], "-", lw=1.2, label="Channel totalPpm")
line_ppm_sum, = ax_ppm.plot([], [], "--", lw=1.2, label="Rx totalPpm corrected")
ax_ppm.legend(loc="upper right", fontsize=8)

line_gain_applied, = ax_gain.plot([], [], "-", lw=1.2, label="Channel gainDb")
line_gain_corr, = ax_gain.plot([], [], "--", lw=1.2, label="Rx gainDb")
ax_gain.legend(loc="upper right", fontsize=8)

line_f_applied, = ax_freq.plot([], [], "-", lw=1.2, label="Channel freqHz")
line_f_corr, = ax_freq.plot([], [], "--", lw=1.2, label="Rx (NCO+PhaseDD) Hz")
ax_freq.legend(loc="upper right", fontsize=8)
ax_freq.tick_params(axis="both", labelsize=8)

hist_max = 600
t_hist = []
ppm_applied = []
ppm_sum = []
ppm_sum_ema = None
gain_applied = []
gain_corr = []
f_applied = []
f_corr = []

stdin_closed = False
latest_pts = None
latest_kv = None
dirty = False
last_draw = 0.0
last_ylim_update = 0.0
x_window_sec = 30.0

while True:
    # Let GUI process events even when no data arrives
    plt.pause(0.001)
    if stdin_closed:
        # End of simulation: keep windows open until user closes them.
        if not plt.fignum_exists(fig.number):
            break
        continue

    # Read and parse as fast as possible to drain the pipe (prevents slowing down the simulation).
    r, _, _ = select.select([sys.stdin], [], [], 0.05)
    if r:
        line = sys.stdin.readline()
        if not line:
            stdin_closed = True
            continue
        line = line.strip()
        if not line:
            continue
        if line.startswith("QUIT"):
            break
        if line.startswith("CONFIG"):
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
            continue
        if line.startswith("FRAME"):
            kv = parse_kv(line)
            pts_x = []
            pts_y = []
            while True:
                l2 = sys.stdin.readline()
                if not l2:
                    stdin_closed = True
                    break
                l2 = l2.strip()
                if l2 == "END":
                    break
                try:
                    x, y = l2.split()
                    pts_x.append(float(x))
                    pts_y.append(float(y))
                except Exception:
                    continue

            latest_pts = (pts_x, pts_y)
            latest_kv = kv
            dirty = True

            # Update history arrays (cheap) but do not redraw yet.
            try:
                if "t_sim" in kv:
                    t = float(kv["t_sim"])
                    t_hist.append(t)
                    ppm_applied.append(float(kv.get("chanTotalPpm", "nan")))
                    vppm = float(kv.get("rxTotalPpmCorr", "nan"))
                    # Light smoothing for display only (EMA).
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
                    if len(t_hist) > hist_max:
                        t_hist[:] = t_hist[-hist_max:]
                        ppm_applied[:] = ppm_applied[-hist_max:]
                        ppm_sum[:] = ppm_sum[-hist_max:]
                        gain_applied[:] = gain_applied[-hist_max:]
                        gain_corr[:] = gain_corr[-hist_max:]
                        f_applied[:] = f_applied[-hist_max:]
                        f_corr[:] = f_corr[-hist_max:]
            except Exception:
                pass

    # Throttled rendering (keeps GUI responsive without blocking the simulation).
    now = time.time()
    if dirty and (now - last_draw) >= draw_period and plt.fignum_exists(fig.number):
        dirty = False
        last_draw = now

        # Constellation points
        if latest_pts is not None:
            px, py = latest_pts
            sc.set_offsets(list(zip(px, py)))

        # Right panel text (English)
        if latest_kv is not None:
            kv = latest_kv
            seg_map = {"0": "stable1", "1": "acc1", "2": "stable2", "3": "acc2"}
            chan_lines = ["[Channel]"]
            if "t_sim" in kv and "t_rate" in kv:
                chan_lines.append(f"t_sim   = {float(kv['t_sim']):.3f} s")
                chan_lines.append(f"t_rate  = {float(kv['t_rate']):.3f} s")
            if "chanSeg" in kv:
                chan_lines.append(f"segment = {seg_map.get(kv['chanSeg'], kv['chanSeg'])}")
            if "chanFreqHz" in kv:
                chan_lines.append(f"freqHz  = {kv['chanFreqHz']}")
            if "chanTotalPpm" in kv:
                chan_lines.append(f"totalPpm= {kv['chanTotalPpm']}")
            if "chanGainDb" in kv or "chanRangeDb" in kv:
                v = kv.get("chanGainDb", kv.get("chanRangeDb", ""))
                chan_lines.append(f"gainDb  = {v}")
            if "samplerMsps" in kv:
                chan_lines.append(f"samplerMsps = {kv['samplerMsps']}")

            rx_lines = ["", "[Receiver]"]
            for key, label in [
                ("symRateMsps", "symRateMsps"),
                ("symRatePpm", "symRatePpm"),
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
            info_text.set_text("\n".join(chan_lines + rx_lines))

        # History lines (no relim/autoscale every frame; update limits at low rate)
        if t_hist:
            line_ppm_applied.set_data(t_hist, ppm_applied)
            line_ppm_sum.set_data(t_hist, ppm_sum)
            line_gain_applied.set_data(t_hist, gain_applied)
            line_gain_corr.set_data(t_hist, gain_corr)
            line_f_applied.set_data(t_hist, f_applied)
            line_f_corr.set_data(t_hist, f_corr)

            # X window
            tmax = t_hist[-1]
            xmin = max(0.0, tmax - x_window_sec)
            ax_ppm.set_xlim(xmin, tmax + 1e-6)

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

        apply_limits()
        fig.canvas.draw_idle()
        fig.canvas.flush_events()

plt.close(fig)

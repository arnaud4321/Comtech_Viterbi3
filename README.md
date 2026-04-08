# Comtech Viterbi3 (`viterbi3`)

Multithreaded C++ simulation of a QPSK-oriented transmit/receive chain with Viterbi decoding.

## What it is

The simulation runs an end-to-end pipeline:

- **TX**: PRBS/file source → scrambling → 3 parallel branches (alignment hypotheses) → convolutional coding → pulse shaping
- **Channel**: AWGN + time-varying impairments (frequency profile, gain profile) + sampling-clock offset model
- **RX**: matched filter → symbol-rate estimation + resampling → coarse frequency estimator (VIVA) + NCO + AGC →
  Gardner timing recovery → decision-directed carrier phase tracking → 3× Viterbi + alignment lock →
  PRBS lock + descrambling

Algorithm documentation lives in the source headers (`src/*.h`) and file blocks (`@file`) and is rendered via Doxygen.

## Build

### Dependencies

The project is built with CMake and links against:

- **UHD** (USRP Hardware Driver, ≥ 4.7) — CMake uses `find_package(UHD 4.7.0 REQUIRED)`; see [Building and Installing UHD from source](https://files.ettus.com/manual/page_build_guide.html) for upstream dependency lists and platform notes.
- **Intel IPP** (`ipps`, `ippcore`; expected under `/opt/intel/oneapi/ipp/latest` unless you override `IPP_ROOT`)
- **Boost** (system, program_options)
- **pthread**

**Using a packaged UHD (typical on Ubuntu):**

```bash
sudo apt-get install libuhd-dev
```

If CMake does not find UHD, point it at your install (e.g. custom prefix): set environment variable `UHD_DIR` or pass `-DCMAKE_PREFIX_PATH=/path/to/uhd/prefix`.

**Build dependencies for compiling UHD itself** (from the Ettus manual — *Setting up the dependencies on Ubuntu*):

```bash
sudo apt-get install autoconf automake build-essential ccache cmake cpufrequtils doxygen ethtool \
  g++ git inetutils-tools libboost-all-dev libncurses5 libncurses5-dev libusb-1.0-0 libusb-1.0-0-dev \
  libusb-dev python3-dev python3-mako python3-numpy python3-requests python3-scipy python3-setuptools \
  python3-ruamel.yaml
```

Fedora/RHEL-style systems: see the same guide for `yum` / `dnf` package lists. Other requirements (compiler, CMake, Boost, LibUSB, Python, Mako, etc.) are summarized under *Build Dependencies* on that page.

On Ubuntu you can install that set in one step with `scripts/install_uhd_build_deps_ubuntu.sh` (runs `sudo apt-get install …`).

**Other packages for this repository:**

- `libboost-system-dev` `libboost-program-options-dev`
- oneAPI IPP installed under `/opt/intel/oneapi/ipp/latest` (or set `IPP_ROOT` in CMake)

### Configure & build

From the repository root:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

The executable is produced under:

- `bin/release/viterbi3` (Release)
- `bin/debug/viterbi3` (Debug)

## Run

`viterbi3` expects a single argument: the path to a JSON configuration file.

```bash
./bin/release/viterbi3 src/config.json
```

If you run in `FILE_TX` mode (see `config.json`), the main thread can also drain receiver output.

## Configuration (`src/config.json`)

The configuration file is read by `Params::ReadParams` and defines:

- **Transmitter**
  - `Method`: 0 = PRBS, 1 = file (`FileName`)
  - `RollOff`: RRC roll-off
- **Channel**
  - `Esn0`: Base Es/N0 in dB
  - `ApplyGainBeforeNoise`: if true, channel gain is applied to the signal before adding noise (effectively varying the SNR). If false, gain is applied to signal+noise.
  - `InitialGainDb`: initial channel gain (dB)
  - `DynamicRangeDb`: gain span applied over the ramp segments (dB)
  - `InitialFrequencyShift`, `FrequencyShift`: carrier offset profile (Hz)
  - `ClockMismatchPpm`: fixed sampling-clock mismatch (ppm)
  - `CarrierToSymbolRateRatio`: maps carrier offset (Hz) to ppm “Doppler” term
  - `AccelerationPeriod`, `StablePeriod`: piecewise schedule (stable → ramp → stable → ramp)
- **CentralFrequency** (coarse frequency estimator + NCO + AGC)
  - `Enable`: 1/0 enables the coarse estimator stage.
  - `EstimationBlockSamples`: number of complex samples (2 sps domain) accumulated per estimate.
  - `EstimatePeriodSec`: re-estimation cadence (seconds).
  - `MaxOffsetHz`: search half-width for the carrier offset (Hz).
  - `PeakToMedianThreshold`: detector threshold (peak/median) on the VIVA FFT magnitude.
  - `FftSizeMultiplier`: zero-padding ratio (\(N_\mathrm{fft} = \text{mult} \cdot N_\mathrm{in}\)).
  - `NcoHzEmaAlpha`: exponential smoothing of the NCO target (0..1).
  - `MaxHzSlewRate`: limits NCO changes smoothly (Hz/s); set to 0 or negative to disable.
  - `PowerEmaAlpha`: AGC power EMA smoothing coefficient.
  - `TargetAvgPower`: if \(>0\), absolute AGC target; otherwise seeded from first chunk and slowly tracks VIVA `AvgPower`.
- **SymbolRateEstimator** (symbol-rate acquisition / tracking)
  - `FftSize`: FFT size used in the estimator window.
  - `MaxOffsetHz`: search half-width around the reference line (Hz).
  - `EstimatePeriodSec`: how often the RX filter thread attempts a new estimate (seconds).
  - `PeakToMedianThreshold`: detection threshold (peak/median).
  - `MaxRelativeJump`: max allowed relative jump vs previous estimate (gates re-estimates).
- **ConstellationDisplay** (optional UI / live plots)
  - `PeriodSec`: how often the receiver sends a FRAME to the UI; set to 0 or negative to disable.
  - `DrawPeriodSec`: GUI redraw throttle (seconds). Independent from `PeriodSec`.
  - `NumSymbols`: how many points to display per frame (capped by the frame size in the sender).
  - `MaxAbs`: plot bounds for I/Q in \([-MaxAbs, +MaxAbs]\).
  - `Backend`: currently `"matplotlib"` is supported.
  - `PythonExe`, `PythonScript`: python interpreter and script path (default `src/constellation_display.py`).
  - `XDisplay`: optional DISPLAY override (e.g. `":0"`). Empty inherits the environment.
  - `Width`, `Height`, `ClearScreen`: legacy ASCII backend options (ignored by matplotlib).
- **Simulation**
  - `Method`: simulation mode selector (see `SimModes` in `src/definitions.h`).
  - `NumErrors`: error budget used as a stop criterion in some modes.
  - `Seed`: RNG seed (0 triggers a hardware-random seed when available).
  - `Debug`: debug pacing / slow-down flags (component-specific).
  - `DisplayPeriodSec`: status print period (seconds). Set to 0 or negative to disable periodic logs (events only).

## Doxygen documentation

Generate HTML docs:

```bash
cmake --build build --target doc
```

Then open:

- `docs/doxygen/html/index.html`

The “main blocks” index is also in `docs/doxygen_mainpage.md`.

## Notes

- The channel publishes “applied” telemetry (segment, frequency, total ppm, gain dB) that can be displayed in the
  constellation GUI.
- The repository contains additional scripts/logs under `data/` used for debugging/analysis.


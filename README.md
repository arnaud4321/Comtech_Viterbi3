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

- **FFTW3f**
- **Boost** (system, program_options)
- **pthread**
- **Intel IPP** (optional: auto-detected if installed under `/opt/intel/oneapi/ipp/latest`)

On Ubuntu-like distributions you typically need packages similar to:

- `fftw3-dev` (or `libfftw3-dev` + float variant if split)
- `libboost-system-dev` `libboost-program-options-dev`
- (optional) oneAPI IPP installed under `/opt/intel/oneapi/ipp/latest`

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
  - `Esn0`: Es/N0 in dB
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


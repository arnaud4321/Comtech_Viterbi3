@mainpage Comtech Viterbi3

## Overview

**Comtech Viterbi3** is a multithreaded C++ simulation of a QPSK-oriented receive chain
with Viterbi decoding, intended for algorithm and link-budget style experiments.

### Where the algorithms are documented

**Signal processing and algorithms are described in the source headers** (`src/*.h`), in the
`@file` / class blocks, not on this page. Use the **Classes**, **Files**, or **Namespaces** tabs
in the generated HTML, or follow the references below.

### Chain (index of main blocks)

- @ref Transmitter
- @ref AWGNChannel
- @ref ChannelSamplingClockOffset (optional sampling-clock offset in the channel)
- @ref Sampler
- @ref Receiver (orchestrates the blocks below)
- @ref SymbolRateEstimator
- @ref ReceiverResampler and @ref Resampler
- @ref ReceiverFreqCorrector and @ref CentralFreqEstimatorViva
- @ref ReceiverTimingTracking and @ref GardnerTiming
- @ref ReceiverPhaseTrackingDD
- @ref Viterbi
- @ref ConstellationDisplay (optional)

### DSP helpers and shared utilities

- @ref Lagrange4Simd — Lagrange-4 interpolation (AVX), used in Gardner timing, SCO resampling, and symbol-rate oversampling.
- @ref IppComplexDft1d — forward complex DFT (float) via Intel IPP; backs @ref SymbolRateEstimator and @ref CentralFreqEstimatorViva (no FFTW dependency).
- @ref ConsoleAlert.h — optional red `[ALERT]` lines for ring backpressure; compile with `ENABLE_CONSOLE_ALERTS` to enable (default: off).

### Building

Configure with **CMake** (see `CMakeLists.txt` in the repository root):

| Dependency | Notes |
|------------|--------|
| **UHD** | `find_package(UHD 4.7.0 REQUIRED)`. Install e.g. `libuhd-dev` (Ubuntu). [Ettus UHD build guide](https://files.ettus.com/manual/page_build_guide.html) lists packages to build UHD from source; optional helper: `scripts/install_uhd_build_deps_ubuntu.sh`. |
| **Intel IPP** | `ipps`, `ippcore`; default `IPP_ROOT=/opt/intel/oneapi/ipp/latest` (CMake cache). |
| **Boost** | `system`, `program_options`. |
| **pthread** | POSIX threads. |

Optional **compile definitions** (documented in source where used):

- `ENABLE_CONSOLE_ALERTS` — emit backpressure / saturation console messages.
- `PHASE_DD_USE_LIBM_ATAN2` — use `std::atan2` instead of the fast polynomial in decision-directed phase error (@ref ReceiverPhaseTrackingDD).
- `ENABLE_RAW_PREVITERBI_METRICS` — enables the continuous raw pre-Viterbi intercorrelation, SER, and BER computation blocks in the Receiver. **WARNING: Must be DISABLED when running with USRP (Mode 1 or 3) as the O(N^2) brute-force correlation causes severe CPU choking and catastrophic t_rate drop.**
- `DEBUG_STATISTICS` — tracks Viterbi decoding metric growth to compute an EMA. It is used to dynamically unlock and resync the Viterbi decoder under very low SNR conditions (e.g., when the channel drops below the "cliff" at ~1.9dB, using a threshold of `kViterbiDesyncGrowthThr = 40.0`).

Generate this HTML reference from the repository root:

```bash
doxygen Doxyfile
```

or, if the CMake target is configured:

```bash
cmake --build build --target doc
```

Output directory: `docs/doxygen/html/` — open `index.html` in a browser.

### API documentation

Class and file reference is generated from Doxygen comments in `src/*.h` and `src/*.cpp`
(third-party `json.hpp` is excluded via `Doxyfile`).

### Configuration (`config.json`)

The simulation is controlled via a JSON configuration file (e.g. `src/config.json`) read by `Params::ReadParams`. The configuration encompasses several major sections:

- **Transmitter**: defines the data source (`Method`: PRBS or File) and RRC `RollOff`.
- **Channel**: controls the AWGN and impairments:
  - `Esn0`: Base Es/N0 in dB.
  - `ApplyGainBeforeNoise`: whether to apply the channel gain directly to the signal (varying the SNR) or to signal+noise.
  - `InitialGainDb`, `DynamicRangeDb`: for piecewise dynamic gain ramps.
  - `InitialFrequencyShift`, `FrequencyShift`: for piecewise dynamic frequency offset ramps.
  - `ClockMismatchPpm`, `CarrierToSymbolRateRatio`: static sampling clock mismatches and Doppler scaling.
  - `AccelerationPeriod`, `StablePeriod`: defines the duration of the piecewise ramps and stable states.
- **CentralFrequency**: parameters for the coarse frequency estimator (VIVA), NCO smoothing, and AGC targets.
- **SymbolRateEstimator**: defines the FFT properties and search window to find the exact symbol rate.
- **TimingTracking**: configures the Gardner timing recovery loop parameters (Kp, Ki).
- **PhaseTracking**: configures the decision-directed phase PLL parameters (Kp, Ki).
- **ConstellationDisplay**: UI and Python plotting arguments.
- **Simulation**: defines the simulation stop modes (e.g., number of errors) and `DisplayPeriodSec` for console telemetry.

For a full parameter listing, refer to the project `README.md`.

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

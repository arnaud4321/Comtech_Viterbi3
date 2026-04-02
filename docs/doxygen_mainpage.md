@mainpage Comtech Viterbi3

## Overview

**Comtech Viterbi3** is a multithreaded C++ simulation of a QPSK-oriented receive chain
with Viterbi decoding, intended for algorithm and link-budget style experiments.

### Where the algorithms are documented

**Signal processing and algorithms are described in the source headers** (`src/*.h`), in the
`@file` / class blocks, not on this page. Use the **Classes** or **Files** tabs in the generated HTML,
or follow the references below.

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

### Building

Configure with CMake (see project `CMakeLists.txt`): Intel IPP, Boost, FFTW3f, pthread.

### API documentation

Class and file reference is generated from Doxygen comments in `src/*.h` (third-party
`json.hpp` is excluded). Run:

    cmake --build build --target doc

Then open `docs/doxygen/html/index.html`.

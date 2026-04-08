/**
 * @file ConsoleAlert.h
 * @brief ANSI red styling for saturation / backpressure @c [ALERT] lines on stderr/stdout.
 *
 * @details By default, **no** alert text is emitted: define @c ENABLE_CONSOLE_ALERTS (e.g.
 * @c -DENABLE_CONSOLE_ALERTS) to print backpressure / ring-saturation messages. Fatal buffer
 * misconfiguration messages before @c std::abort() follow the same switch (abort still occurs).
 */
#pragma once

namespace ConsoleAlert
{
/// Opening: red `[ALERT]` then red message body (reset with @ref kReset after the line).
inline constexpr const char* kRedOpen = "\033[31m[ALERT]\033[0m \033[31m";
/// Green styling for non-alert status markers (reset with @ref kReset).
inline constexpr const char* kGreenOpen = "\033[32m";
inline constexpr const char* kReset = "\033[0m";
} // namespace ConsoleAlert

#if defined(ENABLE_CONSOLE_ALERTS)
/// Expands to @c do { … } while (0) so saturation / backpressure logging runs.
#define CONSOLE_ALERT_STMT(...)                                                                                        \
    do                                                                                                                 \
    {                                                                                                                  \
        __VA_ARGS__                                                                                                    \
    } while (0)
#else
/// No-op: suppress @c [ALERT] lines while keeping real-time behavior unchanged.
#define CONSOLE_ALERT_STMT(...)                                                                                        \
    do                                                                                                                 \
    {                                                                                                                  \
    } while (0)
#endif

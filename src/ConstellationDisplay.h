/**
 * @file ConstellationDisplay.h
 * @brief Optional child-process bridge: sends FRAME lines + IQ samples to Python/matplotlib or ASCII UI.
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <cerrno>
#include <cstring>
#include <fcntl.h>
#include <signal.h>
#include <string>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

namespace
{
inline std::string getSelfExePath()
{
    char buf[4096];
    const ssize_t n = ::readlink("/proc/self/exe", buf, sizeof(buf) - 1);
    if (n <= 0)
        return {};
    buf[n] = '\0';
    return std::string(buf);
}

inline std::string dirnameOf(const std::string& path)
{
    const std::string::size_type pos = path.find_last_of('/');
    if (pos == std::string::npos)
        return {};
    return path.substr(0, pos);
}

inline std::string joinPath(const std::string& a, const std::string& b)
{
    if (a.empty())
        return b;
    if (!a.empty() && a.back() == '/')
        return a + b;
    return a + "/" + b;
}

inline std::string normalizeRepoRootFromExeDir(const std::string& exeDir)
{
    // Expected layout: <repo>/bin/<cfg>/viterbi3 → repo root is two levels up.
    // exeDir = <repo>/bin/<cfg>
    std::string d = exeDir;
    for (int i = 0; i < 2; ++i)
    {
        const auto pos = d.find_last_of('/');
        if (pos == std::string::npos)
            return {};
        d = d.substr(0, pos);
    }
    return d;
}

inline std::string resolvePathFromRepo(const std::string& repoRoot, const std::string& relOrAbs)
{
    if (relOrAbs.empty())
        return relOrAbs;
    if (!relOrAbs.empty() && relOrAbs[0] == '/')
        return relOrAbs;
    if (repoRoot.empty())
        return relOrAbs;
    return joinPath(repoRoot, relOrAbs);
}
} // namespace

/** @brief Backend, cadence, and plot bounds for @ref ConstellationDisplay. */
struct ConstellationDisplayConfig
{
    double PeriodSec = 0.0;   // <=0 disables
    // Max GUI redraw rate (seconds). Independent from PeriodSec (data send rate).
    // <=0: use script default.
    double DrawPeriodSec = 0.0;
    int NumSymbols = 2048;    // number of points to plot from latest frame
    double MaxAbs = 2.0;      // plot range: I,Q in [-MaxAbs, +MaxAbs]
    // Backend selection:
    // - "matplotlib": open a dedicated window via python3 + matplotlib (recommended)
    // - "ascii": fallback in-terminal (legacy)
    std::string Backend = "matplotlib";
    std::string PythonExe = "python3";
    std::string PythonScript = "src/constellation_display.py";
    // Optional X11 display override (e.g. ":0"). Empty = inherit environment.
    std::string XDisplay = "";
    // ASCII-only options (ignored for matplotlib)
    int Width = 61;
    int Height = 31;
    bool ClearScreen = false;
};

/**
 * @brief Optional live constellation / monitoring: forked Python (matplotlib) fed by a simple text protocol.
 *
 * @details **Protocol (stdin of child):** a @c FRAME line with key=value fields (tags, rates, lock flags,
 * channel snapshot, etc.), followed by interleaved ASCII or raw float IQ samples depending on mode, then an
 * @c END line. The child script (see @ref ConstellationDisplayConfig::PythonScript) parses frames and
 * updates scatter plots. No DSP here—only serialization of receiver/channel state for visualization.
 */
class ConstellationDisplay
{
public:
    ConstellationDisplay() = default;
    ~ConstellationDisplay() { Stop(); }

    /** @brief Child PID after successful @ref Start (for logs). */
    pid_t GetChildPid() const { return childPid_; }

    /** @brief Fork @c python matplotlib reader; opens pipe for @ref UpdateEx. */
    bool Start(const ConstellationDisplayConfig& cfg)
    {
        Stop();
        cfg_ = cfg;
        if (cfg_.Backend != "matplotlib")
            return false;

        const std::string exePath = getSelfExePath();
        const std::string exeDir = dirnameOf(exePath);
        const std::string repoRoot = normalizeRepoRootFromExeDir(exeDir);
        const std::string scriptPath = resolvePathFromRepo(repoRoot, cfg_.PythonScript);
        const std::string dataDir = resolvePathFromRepo(repoRoot, "data");
        const std::string stderrPath = joinPath(dataDir, "constellation_display.stderr.log");
        const std::string launchPath = joinPath(dataDir, "constellation_display.launch.log");

        std::string mkdirCmd = "mkdir -p ";
        mkdirCmd += dataDir.empty() ? "data" : dataDir;
        mkdirCmd += " >/dev/null 2>&1 || true";
        std::system(mkdirCmd.c_str());

        // Write a small launch log for client troubleshooting.
        {
            FILE* f = std::fopen(launchPath.c_str(), "wt");
            if (f)
            {
                const char* disp = ::getenv("DISPLAY");
                std::fprintf(f, "pythonExe=%s\n", cfg_.PythonExe.c_str());
                std::fprintf(f, "pythonScript=%s\n", scriptPath.c_str());
                std::fprintf(f, "cfgXDisplay=%s\n", cfg_.XDisplay.c_str());
                std::fprintf(f, "envDISPLAY=%s\n", disp ? disp : "");
                std::fprintf(f, "stderrLog=%s\n", stderrPath.c_str());
                std::fclose(f);
            }
        }

        int fds[2] = {-1, -1};
        if (pipe(fds) != 0)
            return false;

        const pid_t pid = fork();
        if (pid < 0)
        {
            close(fds[0]);
            close(fds[1]);
            return false;
        }

        if (pid == 0)
        {
            // Child
            // stdin from pipe
            dup2(fds[0], STDIN_FILENO);
            close(fds[0]);
            close(fds[1]);

            // Redirect stderr to a file (keeps the main terminal clean).
            int fdErr = open(stderrPath.c_str(), O_CREAT | O_TRUNC | O_WRONLY, 0644);
            if (fdErr >= 0)
            {
                dup2(fdErr, STDERR_FILENO);
                close(fdErr); 
            }

            if (!cfg_.XDisplay.empty())
                setenv("DISPLAY", cfg_.XDisplay.c_str(), 1);

            // exec: python -u script
            const char* argv0 = cfg_.PythonExe.c_str();
            const char* argv1 = "-u";
            const char* argv2 = scriptPath.c_str();
            char* const argv[] = {const_cast<char*>(argv0),
                                  const_cast<char*>(argv1),
                                  const_cast<char*>(argv2),
                                  nullptr};
            execvp(argv0, argv);
            std::fprintf(stderr, "execvp failed: %s\n", std::strerror(errno));
            _exit(127);
        }

        // Parent
        close(fds[0]);
        fdWrite_ = fds[1];
        childPid_ = pid;

        // Send initial config line
        if (dprintf(fdWrite_, "CONFIG maxAbs=%f drawPeriodSec=%f\n", cfg_.MaxAbs, cfg_.DrawPeriodSec) < 0)
        {
            Stop();
            return false;
        }
        return true;
    }

    void Stop()
    {
        if (fdWrite_ >= 0)
        {
            dprintf(fdWrite_, "QUIT\n");
            close(fdWrite_);
            fdWrite_ = -1;
        }
        if (childPid_ > 0)
        {
            int status = 0;
            // Give the child a moment to exit cleanly.
            for (int i = 0; i < 20; ++i)
            {
                pid_t r = waitpid(childPid_, &status, WNOHANG);
                if (r == childPid_)
                {
                    childPid_ = -1;
                    return;
                }
                usleep(10000);
            }
            kill(childPid_, SIGTERM);
            waitpid(childPid_, &status, 0);
            childPid_ = -1;
        }
    }

    // Close the pipe without sending QUIT and without waiting/killing the child.
    // This lets the GUI keep running after normal end-of-simulation (EOF on stdin).
    void ClosePipeKeepAlive()
    {
        if (fdWrite_ >= 0)
        {
            close(fdWrite_);
            fdWrite_ = -1;
        }
        // Detach: do not reap/kill in destructor.
        childPid_ = -1;
    }

    /** @brief Send one constellation frame without extra key=value overlay text. */
    void Update(const float* i, const float* q, int n, double t_sim, double t_rate)
    {
        UpdateEx(i, q, n, t_sim, t_rate, std::string{});
    }

    /**
     * @brief Write @c FRAME line (with optional @a extraKvs), up to @a n IQ pairs, then @c END.
     * @param i In-phase samples (at least @a n).
     * @param q Quadrature samples (at least @a n).
     * @param n Number of complex samples offered (capped by config @c NumSymbols).
     * @param t_sim Wall time since start (seconds), printed on the FRAME line.
     * @param t_rate Time-axis in samples at nominal rate (seconds), printed on the FRAME line.
     * @param extraKvs Space-separated @c key=value tokens appended to the FRAME header (may be empty).
     */
    void UpdateEx(const float* i, const float* q, int n, double t_sim, double t_rate, const std::string& extraKvs)
    {
        if (fdWrite_ < 0 || !i || !q || n <= 0)
            return;
        const int count = std::min(n, std::max(1, cfg_.NumSymbols));
        if (!extraKvs.empty())
        {
            if (dprintf(fdWrite_, "FRAME t_sim=%f t_rate=%f n=%d %s\n", t_sim, t_rate, count, extraKvs.c_str()) < 0)
            {
                Stop();
                return;
            }
        }
        else
        {
            if (dprintf(fdWrite_, "FRAME t_sim=%f t_rate=%f n=%d\n", t_sim, t_rate, count) < 0)
            {
                Stop();
                return;
            }
        }
        for (int k = 0; k < count; ++k)
        {
            if (dprintf(fdWrite_, "%f %f\n", i[k], q[k]) < 0)
            {
                Stop();
                return;
            }
        }
        if (dprintf(fdWrite_, "END\n") < 0)
        {
            Stop();
            return;
        }
    }

private:
    ConstellationDisplayConfig cfg_{};
    int fdWrite_ = -1;
    pid_t childPid_ = -1;
};

inline void RenderConstellationAscii(const float* i, const float* q, int n,
                                     const ConstellationDisplayConfig& cfg,
                                     double t_sim, double t_rate)
{
    if (!i || !q || n <= 0)
        return;
    if (cfg.Width < 21 || cfg.Height < 11 || cfg.MaxAbs <= 0.0)
        return;

    const int W = cfg.Width;
    const int H = cfg.Height;
    std::vector<std::string> grid(static_cast<size_t>(H), std::string(static_cast<size_t>(W), ' '));

    const int cx = W / 2;
    const int cy = H / 2;
    for (int x = 0; x < W; ++x)
        grid[static_cast<size_t>(cy)][static_cast<size_t>(x)] = '-';
    for (int y = 0; y < H; ++y)
        grid[static_cast<size_t>(y)][static_cast<size_t>(cx)] = '|';
    grid[static_cast<size_t>(cy)][static_cast<size_t>(cx)] = '+';

    const double inv = 1.0 / cfg.MaxAbs;
    const int count = std::min(n, std::max(1, cfg.NumSymbols));
    for (int k = 0; k < count; ++k)
    {
        const double ii = static_cast<double>(i[k]);
        const double qq = static_cast<double>(q[k]);
        const int x = static_cast<int>(std::llround((ii * inv) * static_cast<double>(cx) + static_cast<double>(cx)));
        const int y = static_cast<int>(std::llround(((-qq) * inv) * static_cast<double>(cy) + static_cast<double>(cy)));
        if (x >= 0 && x < W && y >= 0 && y < H)
            grid[static_cast<size_t>(y)][static_cast<size_t>(x)] = '*';
    }

    if (cfg.ClearScreen)
        std::cout << "\033[H\033[2J";

    std::cout << "[Constellation] t_sim=" << t_sim << " s"
              << " t_rate=" << t_rate << " s"
              << " N=" << count
              << " maxAbs=" << cfg.MaxAbs
              << " (" << W << "x" << H << ")"
              << std::endl;
    for (const auto& row : grid)
        std::cout << row << std::endl;
}


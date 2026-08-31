// pi_calculator.cpp
//
// Computes Pi to any desired number of decimal digits using the
// Chudnovsky brothers' algorithm, evaluated via BINARY SPLITTING and
// parallelized across CPU cores with std::thread. Arbitrary-precision
// arithmetic is provided by GMP.
//
// Why binary splitting (instead of summing term-by-term)?
//   Summing the Chudnovsky series one term at a time is an inherently
//   sequential chain (each term depends on the running total) and
//   doesn't parallelize. Binary splitting instead represents a range
//   of terms [a, b) as an exact integer triple (P, Q, T) and combines
//   two half-ranges with a handful of big-integer multiplications:
//
//     combine(P1,Q1,T1, P2,Q2,T2) = (P1*P2, Q1*Q2, T1*Q2 + P1*T2)
//
//   The two halves are computed independently, so the recursion tree
//   is embarrassingly parallel -- and as a bonus, binary splitting is
//   also asymptotically FASTER than naive summation even on a single
//   core, because it keeps the big-integer multiplications balanced
//   in size (close to GMP's fast multiplication sweet spot) instead
//   of repeatedly multiplying a huge accumulator by small numbers.
//
// While it runs, a text progress bar is printed to stderr (never to
// stdout, so piping/redirecting the digits themselves is unaffected).
// It only reports progress in units of "leaf terms and combine steps
// completed" -- a reasonable but approximate stand-in for wall-clock
// progress, since the handful of combines nearest the root are far
// more expensive than everything below them. Very fast runs (well
// under a second) print nothing at all, since the first progress poll
// never gets a chance to fire.
//
// Build (Linux/macOS):
//   You need the GMP development headers/library installed:
//     Debian/Ubuntu:  sudo apt-get install libgmp-dev
//     Fedora:         sudo dnf install gmp-devel
//     macOS (brew):   brew install gmp
//
//   Compile with (note -pthread for multithreading):
//     g++ -O2 -std=c++17 -pthread pi_calculator.cpp -o pi_calculator -lgmpxx -lgmp
//
// Build (Windows, native .exe, via MSYS2/MinGW-w64 -- recommended):
//   1. Install MSYS2 from https://www.msys2.org/ and run it once so it
//      can update itself.
//   2. Open the "MSYS2 UCRT64" shell from the Start menu (MSYS2's own
//      recommended default environment for new builds).
//   3. Install the compiler and GMP:
//        pacman -Syu                                  (then reopen the shell if prompted)
//        pacman -S mingw-w64-ucrt-x86_64-gcc mingw-w64-ucrt-x86_64-gmp
//   4. From that same shell, cd to the folder with this file (Windows
//      C:\Users\you\... is reachable as /c/Users/you/...) and compile:
//        g++ -O2 -std=c++17 -static -static-libgcc -static-libstdc++ \
//            pi_calculator.cpp -o pi_calculator.exe -lgmpxx -lgmp -lwinpthread
//      The -static flags bundle the runtime into the .exe so it also
//      runs outside the MSYS2 shell (plain cmd.exe/PowerShell, or on a
//      different machine) with no extra DLLs required.
//   5. Run it from cmd.exe or PowerShell: pi_calculator.exe 1000
//      (double-clicking just flashes a console window and closes it,
//      since this is a command-line program).
//
// Build (Windows, native .exe, via Visual Studio + vcpkg -- alternative):
//   vcpkg install gmp:x64-windows, then vcpkg integrate install, then
//   build pi_calculator.cpp in a normal Visual Studio C++ console
//   project. std::thread needs no special linker flag on MSVC (drop
//   -pthread, which is POSIX-only). If gmpxx.h isn't found, the vcpkg
//   gmp port's C++ bindings support has been known to lag -- ask for a
//   version of this file that uses GMP's plain C API instead.
//
// A note on very large digit counts (several hundred million and up) on
// WINDOWS specifically: GMP's own internal size accounting has a real bug
// on Windows builds (32-bit "long" there vs. 64-bit on Linux/Mac) that can
// surface as a "Cannot allocate memory" crash for an impossible, non-real
// size, late into an otherwise-successful run, regardless of how much RAM
// is actually free. It's not specific to this program or to how the
// prebuilt .exe was built. See the accompanying README's appendix for the
// full evidence. Practical effect: below a couple hundred million digits,
// Windows is solid; past that, prefer the WSL/Linux build (a 1-billion-
// digit run has been verified correct there against a published
// reference), since Linux's 64-bit "long" doesn't have room for this bug
// to occur.
//
// Usage:
//   ./pi_calculator <digits> [output_file] [threads]
//   ./pi_calculator <digits> [-o|--output FILE] [-t|--threads N] [--debug-mem]
//
//   <digits>       Number of decimal digits of Pi to compute (after "3.").
//   [output_file]  Optional. Write digits to this file instead of stdout
//                  (recommended for large digit counts). Pass "-" to
//                  keep stdout while still supplying a thread count.
//   [threads]      Optional. Number of worker threads to use. Defaults
//                  to the number of hardware threads detected.
//   --debug-mem    Optional. Logs every allocation of 20 MB or more to
//                  stderr as it happens. Off by default; only useful for
//                  chasing down an allocation failure (see memtrack below).
//
// Examples:
//   ./pi_calculator 1000
//   ./pi_calculator 1000000 pi_1M.txt
//   ./pi_calculator 1000000 -o pi_1M.txt -t 8
//   ./pi_calculator 1000000000 -o pi_1B.txt -t 8 --debug-mem

#include <gmpxx.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>

namespace chudnovsky {

// Chudnovsky series constants:
//   1/Pi = 12 * sum_{k=0}^inf (-1)^k (6k)! (A + B*k) / ((3k)! (k!)^3 C^(3k+3/2))
constexpr long A = 13591409;
constexpr long B = 545140134;
constexpr long C = 640320;

// A leaf term (for a single index a >= 1) is:
//   P(a) = (6a-5)(2a-1)(6a-1)
//   Q(a) = a^3 * C^3 / 24        (always an exact integer)
//   T(a) = (-1)^a * P(a) * (A + B*a)
// These come from simplifying the ratio between consecutive Chudnovsky
// terms; see the accompanying derivation notes in the project history.
// Two ranges combine as:
//   P(a,c) = P(a,b) * P(b,c)
//   Q(a,c) = Q(a,b) * Q(b,c)
//   T(a,c) = T(a,b) * Q(b,c) + P(a,b) * T(b,c)

struct PQT {
    mpz_class P, Q, T;
};

// Maximum number of worker threads allowed to be "in flight" at once
// (including the ones already spawned by earlier recursive calls).
// Recursion falls back to plain sequential calls once this budget is
// exhausted or a range gets too small to be worth the thread overhead.
std::atomic<int> g_threadBudget{0};
constexpr long MIN_TERMS_FOR_THREAD = 4000;

// Progress tracking: each leaf and each internal combine reports its
// *range size* (b-a) as its weight, rather than counting every node
// equally. A big-integer multiplication's cost grows with operand size,
// and a combine near the root spans (and thus multiplies numbers sized
// close to) the whole term range, so weighting by range size is a much
// closer proxy for actual work than a flat "1 unit per node" count --
// which would otherwise hit ~100% while the handful of most expensive
// combines nearest the root are still running, making the bar look
// stuck. It's still an approximation (multiplication cost grows faster
// than linearly with size), so the reporter additionally clamps the
// displayed percentage below 100% until the computation has genuinely
// finished.
std::atomic<long long> g_workDone{0};

PQT binarySplit(long a, long b) {
    if (b - a == 1) {
        PQT r;
        if (a == 0) {
            r.P = 1;
            r.Q = 1;
        } else {
            mpz_class za(a);
            r.P = (6 * za - 5) * (2 * za - 1) * (6 * za - 1);
            mpz_class a3 = za * za * za;
            static const mpz_class C3_24 = mpz_class(C) * C * C / 24;
            r.Q = a3 * C3_24;
        }
        r.T = r.P * (mpz_class(A) + mpz_class(B) * mpz_class(a));
        if (a % 2 != 0) r.T = -r.T;
        g_workDone.fetch_add(1, std::memory_order_relaxed);  // range size = 1
        return r;
    }

    long m = a + (b - a) / 2;

    bool spawnThread = (b - a) >= MIN_TERMS_FOR_THREAD &&
                        g_threadBudget.fetch_sub(1, std::memory_order_relaxed) > 0;
    if (!spawnThread) {
        // Undo the speculative decrement if we didn't actually use it.
        if ((b - a) >= MIN_TERMS_FOR_THREAD) g_threadBudget.fetch_add(1, std::memory_order_relaxed);

        PQT left = binarySplit(a, m);
        PQT right = binarySplit(m, b);
        PQT out;
        out.P = left.P * right.P;
        out.Q = left.Q * right.Q;
        out.T = left.T * right.Q + left.P * right.T;
        g_workDone.fetch_add(b - a, std::memory_order_relaxed);
        return out;
    }

    PQT left, right;
    std::thread worker([&]() { left = binarySplit(a, m); });
    right = binarySplit(m, b);
    worker.join();
    g_threadBudget.fetch_add(1, std::memory_order_relaxed);

    PQT out;
    out.P = left.P * right.P;
    out.Q = left.Q * right.Q;
    out.T = left.T * right.Q + left.P * right.T;
    g_workDone.fetch_add(b - a, std::memory_order_relaxed);
    return out;
}

// Approximates the total weight (sum of range sizes over every node --
// leaves and internal combines -- in the binary-splitting recursion
// tree for N leaves) as N * (depth + 2), a slight overestimate of the
// true sum. Cheap to compute regardless of how large N gets, unlike
// exactly simulating the recursion shape.
long long approxTotalWork(long numTerms) {
    long n = std::max(2L, numTerms);
    double depth = std::ceil(std::log2(static_cast<double>(n)));
    return static_cast<long long>(numTerms) * static_cast<long long>(depth + 2);
}

// Prints a simple text progress bar to stderr (so it never contaminates
// digits written to stdout) until `done` is set. Polls rather than
// blocks, so a computation that finishes before the first poll interval
// elapses never prints anything at all -- fast runs stay silent. The
// displayed percentage is clamped below 100% until `done` actually
// becomes true, since the work-weighting above is only approximate and
// should never claim completion early.
void progressReporter(long long totalWork, std::atomic<bool>& done,
                       std::chrono::steady_clock::time_point start) {
    using namespace std::chrono_literals;
    const int barWidth = 30;
    bool printedAnything = false;

    while (!done.load(std::memory_order_relaxed)) {
        std::this_thread::sleep_for(150ms);
        if (done.load(std::memory_order_relaxed)) break;

        long long workDone = g_workDone.load(std::memory_order_relaxed);
        double pct = totalWork > 0
                         ? (100.0 * std::min(workDone, totalWork) / totalWork)
                         : 0.0;
        pct = std::min(pct, 99.9);  // never claim done until we actually are
        double elapsed =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
                .count();

        int filled = static_cast<int>(pct / 100.0 * barWidth);
        std::ostringstream line;
        line << "\r[" << std::string(filled, '#')
             << std::string(barWidth - filled, '-') << "] " << std::fixed
             << std::setprecision(1) << std::setw(5) << pct << "%  ("
             << elapsed << "s elapsed)";
        std::cerr << line.str() << std::flush;
        printedAnything = true;
    }

    if (printedAnything) {
        // Overwrite the bar with a completed one, then move to a fresh line.
        std::ostringstream line;
        line << "\r[" << std::string(barWidth, '#') << "] 100.0%  (done)";
        std::cerr << line.str() << std::string(10, ' ') << "\n";
    }
}

// Runs `work` while printing a live "<label>... N.Ns" status to stderr
// every ~150ms, so a slow single-shot GMP call (a huge sqrt, multiply,
// divide, or decimal conversion -- none of which have any natural
// notion of partial progress) doesn't look like a silent hang. Like the
// tree progress bar, a stage that finishes before the first poll prints
// nothing at all.
void runStage(const std::string& label, const std::function<void()>& work) {
    using namespace std::chrono_literals;
    std::atomic<bool> stageDone{false};
    auto stageStart = std::chrono::steady_clock::now();
    bool printedAnything = false;

    std::thread spinner([&]() {
        while (!stageDone.load(std::memory_order_relaxed)) {
            std::this_thread::sleep_for(150ms);
            if (stageDone.load(std::memory_order_relaxed)) break;
            double elapsed = std::chrono::duration<double>(
                                  std::chrono::steady_clock::now() - stageStart)
                                  .count();
            std::ostringstream line;
            line << "\r" << label << "... " << std::fixed
                 << std::setprecision(1) << elapsed << "s";
            std::cerr << line.str() << std::string(10, ' ') << std::flush;
            printedAnything = true;
        }
    });

    work();

    stageDone.store(true, std::memory_order_relaxed);
    spinner.join();  // join() synchronizes-with the loop above, so reading
                      // printedAnything below is race-free.

    if (printedAnything) {
        double elapsed = std::chrono::duration<double>(
                              std::chrono::steady_clock::now() - stageStart)
                              .count();
        std::ostringstream line;
        line << "\r" << label << "... done (" << std::fixed
             << std::setprecision(2) << elapsed << "s)";
        std::cerr << line.str() << std::string(10, ' ') << "\n";
    }
}

// Prints a rough heads-up about memory requirements for large digit
// counts. The exact peak usage depends on GMP internals and thread
// scheduling and isn't something we can predict precisely, but the
// qualitative shape is predictable: the intermediate P/Q/T integers
// binary splitting works with are noticeably bigger than the final
// output (their bit length grows roughly with an extra log2(numTerms)
// factor), several of them are alive at once per combine step, and
// each concurrent thread multiplies that further. This is printed
// instead of finding out the hard way, several minutes into a run,
// via a GMP "Cannot allocate memory" crash.
void printMemoryNotice(unsigned long digits, unsigned long precisionBits,
                        unsigned int numThreads) {
    double outputMB = precisionBits / 8.0 / (1024.0 * 1024.0);
    std::cerr << "Target precision: ~" << std::fixed << std::setprecision(0)
               << outputMB << " MB. ";
    if (digits >= 10'000'000) {
        std::cerr << "\nNote: at " << digits
                   << " digits, the intermediate numbers used internally "
                      "run noticeably larger than that, several are held "
                      "in memory at once, and this multiplies with thread "
                      "count (currently "
                   << numThreads
                   << "). If you hit a 'Cannot allocate memory' error on a "
                      "run with genuinely limited RAM, re-running with a "
                      "lower thread count (-t N) trades speed for less peak "
                      "memory use. On Windows, well past a couple hundred "
                      "million digits, that same error can also come from a "
                      "known GMP bug that thread count does NOT fix -- see "
                      "the README's appendix; the WSL/Linux build avoids "
                      "it.\n";
    }
}

// Computes Pi to `digits` decimal digits and returns it as a decimal
// string (first character is the integer part "3").
std::string computePi(unsigned long digits, unsigned int numThreads) {
    const double DIGITS_PER_TERM = 14.1816474627254776555;
    const unsigned long guardDigits = 20;
    const unsigned long totalDigits = digits + guardDigits;
    const unsigned long precisionBits =
        static_cast<unsigned long>(totalDigits * 3.3219280948873623) + 64;

    mpf_set_default_prec(precisionBits);

    printMemoryNotice(digits, precisionBits, numThreads);

    const long numTerms =
        static_cast<long>(totalDigits / DIGITS_PER_TERM) + 2;

    // Each spawned thread "uses up" one unit of budget; the calling
    // thread itself is free. Using (numThreads - 1) as the budget caps
    // total concurrency at roughly numThreads threads overall.
    g_threadBudget.store(numThreads > 1 ? static_cast<int>(numThreads - 1) : 0,
                          std::memory_order_relaxed);

    g_workDone.store(0, std::memory_order_relaxed);
    const long long totalWork = approxTotalWork(numTerms);
    std::atomic<bool> progressDone{false};
    auto progressStart = std::chrono::steady_clock::now();
    std::thread reporter(progressReporter, totalWork, std::ref(progressDone),
                          progressStart);

    PQT r = binarySplit(0, numTerms);

    progressDone.store(true, std::memory_order_relaxed);
    reporter.join();

    // pi = (Q * 426880 * sqrt(10005)) / T
    // (426880*sqrt(10005) = C^(3/2)/12, using C = 640320 = 64*10005)
    //
    // Each of these is a single big GMP call with no internal notion of
    // partial progress, so they're wrapped individually with runStage()
    // rather than left to run silently -- for large digit counts, the
    // division and the final decimal conversion in particular can each
    // take a while on their own.
    // mpf_set_default_prec(precisionBits) was already called above, so
    // these default-constructed mpf_class values all pick it up.
    mpf_class Qf, Tf, sqrt10005;
    runStage("Converting integers to floating point", [&]() {
        Qf = r.Q;
        Tf = r.T;
    });
    runStage("Computing sqrt(10005)", [&]() {
        mpf_sqrt(sqrt10005.get_mpf_t(), mpf_class(10005).get_mpf_t());
    });

    mpf_class numerator;
    runStage("Multiplying",
             [&]() { numerator = Qf * mpf_class(426880) * sqrt10005; });

    mpf_class pi;
    runStage("Dividing (often the slowest single step)",
             [&]() { pi = numerator / Tf; });

    // Extract `digits` decimal digits after the leading "3" by
    // truncating (chopping) rather than rounding -- see the note in
    // the original single-threaded version for why: the standard
    // convention for "Pi to N digits" is truncation, matching
    // published digit tables, not round-to-nearest at the boundary.
    mp_exp_t exponent;
    std::string raw;
    runStage("Converting to decimal digits",
             [&]() { raw = pi.get_str(exponent, 10, totalDigits + 1); });

    std::string result;
    if (exponent <= 0) {
        result = std::string(-exponent, '0') + raw;
    } else {
        result = raw;
    }

    if (result.size() > digits + 1) result.resize(digits + 1);
    while (result.size() < digits + 1) result += '0';

    return result;
}

}  // namespace chudnovsky

// Diagnostic allocator: intercepts every malloc/realloc/free GMP makes so
// that when (if) an allocation fails, we know exactly how big it was and
// how much was already allocated at the time -- instead of just knowing
// which stage of the computation we were in. Logs any single allocation
// or growth of at least LOG_THRESHOLD_MB to stderr as it happens, so the
// last few lines before a crash are the direct lead-up to it. This is
// opt-in (pass --debug-mem) rather than on by default: it was what
// originally pinned down a real GMP-on-Windows bug (see the README), and
// is kept around in case anything similar ever needs chasing down again,
// but it adds console noise and a small overhead that most runs don't
// need.
namespace memtrack {

constexpr long long LOG_THRESHOLD_BYTES = 20LL * 1024 * 1024;  // 20 MB
std::atomic<long long> g_currentBytes{0};
std::atomic<long long> g_peakBytes{0};

void updatePeak(long long current) {
    long long peak = g_peakBytes.load(std::memory_order_relaxed);
    while (current > peak &&
           !g_peakBytes.compare_exchange_weak(peak, current,
                                               std::memory_order_relaxed)) {
    }
}

void logLine(const char* tag, double sizeMB, double currentMB, double peakMB) {
    char buf[256];
    std::snprintf(buf, sizeof(buf), "[mem] %-7s %10.1f MB   (in use: %8.1f MB, peak: %8.1f MB)\n",
                  tag, sizeMB, currentMB, peakMB);
    std::fputs(buf, stderr);
    std::fflush(stderr);
}

void* trackedMalloc(size_t size) {
    void* p = std::malloc(size);
    long long cur = g_currentBytes.fetch_add(static_cast<long long>(size),
                                              std::memory_order_relaxed) +
                     static_cast<long long>(size);
    updatePeak(cur);
    if (!p) {
        logLine("MALLOC FAILED", size / (1024.0 * 1024.0),
                cur / (1024.0 * 1024.0), g_peakBytes.load() / (1024.0 * 1024.0));
    } else if (static_cast<long long>(size) >= LOG_THRESHOLD_BYTES) {
        logLine("malloc", size / (1024.0 * 1024.0), cur / (1024.0 * 1024.0),
                g_peakBytes.load() / (1024.0 * 1024.0));
    }
    return p;
}

void* trackedRealloc(void* ptr, size_t oldSize, size_t newSize) {
    void* p = std::realloc(ptr, newSize);
    long long delta =
        static_cast<long long>(newSize) - static_cast<long long>(oldSize);
    long long cur =
        g_currentBytes.fetch_add(delta, std::memory_order_relaxed) + delta;
    updatePeak(cur);
    if (!p && newSize > 0) {
        logLine("REALLOC FAILED (old->new)", newSize / (1024.0 * 1024.0),
                cur / (1024.0 * 1024.0), g_peakBytes.load() / (1024.0 * 1024.0));
    } else if (delta >= LOG_THRESHOLD_BYTES) {
        logLine("realloc", newSize / (1024.0 * 1024.0), cur / (1024.0 * 1024.0),
                g_peakBytes.load() / (1024.0 * 1024.0));
    }
    return p;
}

void trackedFree(void* ptr, size_t size) {
    std::free(ptr);
    g_currentBytes.fetch_sub(static_cast<long long>(size), std::memory_order_relaxed);
}

}  // namespace memtrack

void printUsage(const char* prog) {
    std::cerr << "Usage: " << prog << " <digits> [output_file] [threads]\n";
    std::cerr << "   or: " << prog
               << " <digits> [-o|--output FILE] [-t|--threads N] [--debug-mem]\n";
    std::cerr << "  digits:  number of decimal digits of Pi to compute\n";
    std::cerr << "  output_file / -o FILE : write digits to this file "
                 "(default: stdout; use \"-\" to keep stdout explicitly)\n";
    std::cerr << "  threads / -t N        : worker thread count (default: "
                 "hardware_concurrency, auto-lowered for very large digit "
                 "counts -- see below)\n";
    std::cerr << "  --debug-mem           : log every allocation of 20 MB "
                 "or more to stderr as it happens (off by default; only "
                 "useful for chasing down an allocation failure)\n";
    std::cerr << "The -o/-t flags and the older positional form can't be "
                 "mixed for the same slot; use one style or the other.\n";
}

int main(int argc, char** argv) {
    // Scanned separately, and before any GMP object is constructed, since
    // mp_set_memory_functions() must be called before GMP allocates
    // anything and the full option parser below runs too late for that.
    bool debugMem = false;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--debug-mem") {
            debugMem = true;
            break;
        }
    }
    if (debugMem) {
        mp_set_memory_functions(memtrack::trackedMalloc,
                                 memtrack::trackedRealloc,
                                 memtrack::trackedFree);
    }

    if (argc < 2) {
        printUsage(argv[0]);
        return 1;
    }

    if (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help") {
        printUsage(argv[0]);
        return 0;
    }

    long digitsArg = std::strtol(argv[1], nullptr, 10);
    if (digitsArg <= 0) {
        std::cerr << "Error: digits must be a positive integer.\n";
        return 1;
    }
    unsigned long digits = static_cast<unsigned long>(digitsArg);

    unsigned int numThreads = std::thread::hardware_concurrency();
    if (numThreads == 0) numThreads = 1;
    bool threadsExplicit = false;
    std::string outputFile;
    bool haveOutputFile = false;

    // Two argument styles are accepted at once, on purpose: the original
    // positional form (<digits> [output_file] [threads]) stays working
    // for backward compatibility, but a bare word in that 3rd slot is an
    // easy way to silently lose a thread count -- e.g. running just
    // "pi_calculator 500000000 8" without an output file passes "8" as
    // the OUTPUT FILENAME, not the thread count, and the thread count
    // quietly falls back to the default/auto-cap with no error at all.
    // -o/--output and -t/--threads are explicit and never ambiguous;
    // prefer them. Any argument that doesn't parse into either slot (via
    // either style) is now a hard error instead of being silently
    // ignored or misinterpreted.
    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-t" || arg == "--threads") {
            if (i + 1 >= argc) {
                std::cerr << "Error: " << arg << " requires a value.\n";
                return 1;
            }
            long t = std::strtol(argv[++i], nullptr, 10);
            if (t <= 0) {
                std::cerr << "Error: thread count must be a positive integer, got '"
                           << argv[i] << "'.\n";
                return 1;
            }
            numThreads = static_cast<unsigned int>(t);
            threadsExplicit = true;
        } else if (arg == "-o" || arg == "--output") {
            if (i + 1 >= argc) {
                std::cerr << "Error: " << arg << " requires a value.\n";
                return 1;
            }
            outputFile = argv[++i];
            haveOutputFile = true;
        } else if (arg == "--debug-mem") {
            // Already handled above (before mp_set_memory_functions), so
            // there's nothing left to do here except recognize it and not
            // fall through to the "unrecognized option" error below.
        } else if (!arg.empty() && arg[0] == '-' && arg != "-") {
            std::cerr << "Error: unrecognized option '" << arg << "'.\n";
            printUsage(argv[0]);
            return 1;
        } else if (!haveOutputFile) {
            // Legacy positional slot 1: output file (or "-" for stdout).
            outputFile = arg;
            haveOutputFile = true;
        } else if (!threadsExplicit) {
            // Legacy positional slot 2: thread count. Unlike the old
            // behavior, a value that doesn't parse is now a hard error
            // instead of being silently dropped.
            long t = std::strtol(arg.c_str(), nullptr, 10);
            if (t <= 0) {
                std::cerr << "Error: expected a thread count here, got '"
                           << arg << "'. If you meant an option, use -t N.\n";
                return 1;
            }
            numThreads = static_cast<unsigned int>(t);
            threadsExplicit = true;
        } else {
            std::cerr << "Error: unexpected extra argument '" << arg << "'.\n";
            printUsage(argv[0]);
            return 1;
        }
    }

    // Safety default: an unspecified thread count on a very large digit
    // count is exactly the combination that risks running out of memory
    // (each concurrent thread can hold its own large intermediate
    // numbers at once). Auto-cap in that case rather than defaulting to
    // every hardware thread; an explicit thread count (-t/--threads or
    // the legacy positional slot) always overrides this and is never
    // capped, since that's a deliberate choice.
    constexpr unsigned long AUTO_CAP_DIGITS_THRESHOLD = 50'000'000;
    constexpr unsigned int AUTO_CAP_THREADS = 4;
    if (!threadsExplicit && digits >= AUTO_CAP_DIGITS_THRESHOLD &&
        numThreads > AUTO_CAP_THREADS) {
        std::cerr << "Note: " << digits << " digits with no thread count given "
                   << "-- defaulting to " << AUTO_CAP_THREADS
                   << " threads instead of all " << numThreads
                   << " hardware threads, to limit peak memory use. Pass an "
                      "explicit thread count with -t N (e.g. a higher count "
                      "if you have plenty of RAM, or lower if you still run "
                      "out of memory).\n";
        numThreads = AUTO_CAP_THREADS;
    }

    auto start = std::chrono::steady_clock::now();
    std::string piDigits = chudnovsky::computePi(digits, numThreads);
    auto end = std::chrono::steady_clock::now();
    double seconds = std::chrono::duration<double>(end - start).count();

    std::string formatted = piDigits.substr(0, 1) + "." + piDigits.substr(1);

    if (haveOutputFile && outputFile != "-") {
        std::ofstream out(outputFile);
        if (!out) {
            std::cerr << "Error: could not open output file " << outputFile << "\n";
            return 1;
        }
        out << formatted << "\n";
        out.close();
        std::cerr << "Wrote " << digits << " digits of Pi to " << outputFile
                   << " using " << numThreads << " thread(s) in " << seconds
                   << " s\n";
    } else {
        std::cout << formatted << "\n";
        std::cerr << "(" << digits << " digits computed using " << numThreads
                   << " thread(s) in " << seconds << " s)\n";
    }

    return 0;
}

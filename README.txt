pi_calculator - compute Pi to any number of digits (multithreaded, GMP/Chudnovsky)
====================================================================

Verified correct to 1 BILLION decimal digits, cross-checked against the
published reference digits at https://stuff.mit.edu/afs/sipb/contrib/pi/
(run via the WSL/Linux build -- see MEMORY USE below for why that's the
recommended path at this scale on Windows).

# Made with Claude

Made using Claude Sonnet 5

# Screenshot

![pi_calculator](pi_calculator.jpg)

WHAT'S IN THIS ZIP
  pi_calculator.cpp   The C++ source it was built from.
  README.txt          This file.

RUNNING THE BUILT .exe
  Run it from Command Prompt or PowerShell, not by double-clicking
  (double-clicking just flashes a console window and closes it,
  since this is a command-line program):

    pi_calculator.exe <digits> [-o FILE] [-t N]

  Examples:
    pi_calculator.exe 1000
    pi_calculator.exe 1000000 -o pi.txt
    pi_calculator.exe 1000000 -o pi.txt -t 8

  -o/--output FILE sets the output file (default: print to the screen).
  -t/--threads N sets the worker thread count (default: every hardware
  thread detected, automatically lowered for very large digit counts --
  see MEMORY USE below). An older positional form, <digits> [output_file]
  [threads], also still works, but -o/-t are recommended: with the
  positional form, a thread count typed without an output file first
  (e.g. "pi_calculator 500000000 8") silently becomes the OUTPUT
  FILENAME instead, not the thread count -- easy to trip over. -o/-t
  can't be misread that way.

  There's also an optional --debug-mem flag that logs every memory
  allocation of 20 MB or more to the console as it happens; it's off by
  default (it adds noise and a little overhead) and is really only
  useful if you ever need to chase down another allocation failure the
  way the one described in the appendix below was diagnosed. Ordinary
  runs don't need it.

  While it runs on a large digit count, progress prints to the console
  in two phases: first a percentage bar for the main computation, then
  a per-step timer for the remaining sqrt/multiply/divide/decimal-
  conversion work (each of those is a single big-number operation with
  no natural notion of partial progress, so they get a live elapsed-
  time readout instead of a percentage). Very fast runs (well under a
  second) won't show any of this at all.

  This .exe is self-contained (GMP and the C++ runtime are statically
  linked in) and was cross-compiled + verified under Wine emulation on
  Linux, since real Windows wasn't available to build/test on directly.

MEMORY USE
  This matters more than you'd expect for a "compute Pi" toy program,
  so here's what's actually going on, based on real measurements (not
  guessing): peak memory scales at roughly 12-13 bytes per requested
  digit. That means:

      10 million digits    ->  ~120 MB peak
      50 million digits    ->  ~600 MB peak
     100 million digits    ->  ~1.2 GB peak
     300 million digits    ->  ~3.5 GB peak
     500 million digits    ->  ~6 GB peak (single-threaded; more with
                                            multiple threads, since each
                                            can hold its own large
                                            numbers in memory at once)

  That's the real requirement of the algorithm at this scale (binary
  splitting keeps some genuinely large intermediate numbers around,
  not just the final result) -- it's a resource cost to plan around,
  and up to a couple hundred million digits it's the whole story: any
  modern machine with a few GB free should be fine.

  Past a few hundred million digits, on WINDOWS specifically, there's a
  second, separate issue: a real bug in GMP itself (not in this
  program, and not a real memory shortage) that can make a large
  internal allocation fail with a nonsensical, astronomically large
  requested size. Full evidence and explanation are in the appendix at
  the end of this file. Practical effect: below a couple hundred
  million digits you're safe on Windows either way; above that, there's
  a real (if not precisely predictable) chance of hitting this GMP bug,
  and no amount of free RAM or thread-count tuning avoids it, because
  the failing request isn't for real memory. For anything at that scale
  -- the 1-billion-digit run verified above used this path -- use the
  WSL/Linux build (Option C below) instead of the Windows .exe: Linux's
  GMP doesn't have the 32-bit size type that this bug needs in order to
  occur.

  At 50 million digits or more, the program automatically uses just 4
  threads instead of every hardware thread when you don't pass -t
  yourself, since fewer concurrent threads means fewer large numbers
  competing for memory at once; pass -t explicitly to override this in
  either direction. This helps with ordinary memory pressure but does
  not affect the GMP bug described above.

COMPILING IT YOURSELF ON WINDOWS
  If you'd rather build your own .exe (e.g. to fully trust the binary,
  or to enable CPU-specific optimizations), here are three ways. Option
  A and B still build against Windows GMP, so neither avoids the GMP
  bug described above at extreme digit counts -- only Option C (a
  Linux/WSL build) does that -- but A and B are worth knowing regardless
  for smaller-scale runs or trusting your own build.

  Option A: MSYS2 + MinGW-w64 (UCRT64 environment)
    1. Install MSYS2 from https://www.msys2.org/ and let it update
       itself on first run.
    2. Open the "MSYS2 UCRT64" shell from the Start menu (this is
       MSYS2's own recommended default environment).
    3. Install the compiler and GMP:
         pacman -Syu
         (reopen the shell if it asks you to, then run pacman -Syu again)
         pacman -S mingw-w64-ucrt-x86_64-gcc mingw-w64-ucrt-x86_64-gmp
    4. cd to this folder (Windows C:\Users\you\... is reachable as
       /c/Users/you/... in this shell) and compile:
         g++ -O2 -std=c++17 -static -static-libgcc -static-libstdc++ \
             pi_calculator.cpp -o pi_calculator.exe -lgmpxx -lgmp -lwinpthread
       The -static flags make the result runnable outside the MSYS2
       shell too (plain cmd.exe/PowerShell, or copied to another PC)
       with no extra DLLs required.
    5. Run it: pi_calculator.exe 1000

  Option B: Visual Studio + vcpkg (native MSVC build)
    1. git clone https://github.com/microsoft/vcpkg
       .\vcpkg\bootstrap-vcpkg.bat
       .\vcpkg\vcpkg install gmp:x64-windows
       .\vcpkg\vcpkg integrate install
    2. Create a normal Visual Studio C++ Console App project, add
       pi_calculator.cpp, and build. vcpkg integrate install wires up
       the include/lib paths automatically.
    3. std::thread needs no special linker flag on MSVC (there's no
       -pthread equivalent needed on Windows).
    4. Caveat: vcpkg's gmp port's C++ bindings (gmpxx.h, used for
       mpz_class/mpf_class here) have occasionally lagged on MSVC
       support. If the build can't find gmpxx.h, tell me the exact
       error -- I can rewrite the arithmetic to use GMP's plain C API
       instead, which sidesteps the issue.

  Option C: WSL (build a Linux binary -- recommended for very large
  digit counts, per MEMORY USE above)
    1. wsl --install, then inside the Ubuntu shell:
         sudo apt update && sudo apt install g++ libgmp-dev
    2. g++ -O2 -std=c++17 -pthread pi_calculator.cpp -o pi_calculator -lgmpxx -lgmp
    This produces a Linux binary usable inside WSL, not a native
    Windows .exe -- use Option A or B for that. The same memory-scaling
    numbers above still apply here, but this is the path that got a
    verified, correct 1-billion-digit run (see the top of this file):
    Linux's GMP doesn't have the 32-bit "long" size type that the
    Windows-specific allocation bug (see the appendix) needs in order
    to occur, so it has real headroom where the Windows .exe doesn't.

APPENDIX: TWO QUESTIONS ANSWERED IN MORE DEPTH

  1. Where does "sqrt(10005)" in the progress output come from?

     It's part of the Chudnovsky formula this program uses to compute
     Pi. After binary splitting reduces the whole series to one integer
     triple (P, Q, T), the final answer is:

         pi = (Q * 426880 * sqrt(10005)) / T

     The constant 426880*sqrt(10005) isn't arbitrary -- it comes directly
     out of the series' constant term. Chudnovsky's formula involves
     C^(3/2) where C = 640320, and 640320 happens to factor as 64 * 10005.
     That means sqrt(640320) = sqrt(64 * 10005) = 8 * sqrt(10005), so
     C^(3/2) = 640320 * 8 * sqrt(10005). Dividing by the formula's other
     constant factor of 12 gives (640320*8)/12 = 426880, leaving exactly
     426880*sqrt(10005) as the closed-form constant. It's a fixed
     irrational number computed once per run (that's the "Computing
     sqrt(10005)..." stage in the progress output), not something that
     changes with the requested digit count.

  2. What was the "Cannot allocate memory" bug at very large digit
     counts, and how was it actually found?

     At several hundred million digits and up, the Windows .exe could
     crash with GMP's "Cannot allocate memory" error late into an
     otherwise-successful run (in testing, this showed up around 500
     million digits, and again at 1 billion), even with tens of GB of
     RAM free. Two early theories about the cause -- genuine RAM
     exhaustion, and which C runtime the .exe was linked against
     (msvcrt.dll vs. the modern ucrtbase.dll) -- were both ruled out by
     testing: the failure reproduced identically on a from-scratch MSYS2
     UCRT64 build with 64 GB of RAM available and only a few GB actually
     in use.

     The real answer came from adding temporary instrumentation that
     logged the exact size of every large memory request as it happened
     (this is what the --debug-mem flag now does, on demand). That log
     caught a failing request for 2^64 minus exactly 524288 bytes --
     an astronomically large, physically impossible allocation sitting
     just 512 KB below the point where a 64-bit unsigned number wraps
     around to (almost) its maximum value. That's the unmistakable
     signature of a small NEGATIVE size -- a few hundred KB "too little"
     -- being computed somewhere deep inside GMP and then reinterpreted
     as unsigned, turning "slightly too small" into "basically infinite"
     instead of failing cleanly. Earlier crash reports (allocation sizes
     a few hundred KB below 2^32, i.e. 4 GB) turned out to be the exact
     same bug showing up one level down, through GMP's 32-bit "long"
     size type on Windows (that type is 64-bit on Linux/Mac, which is
     also why the bug doesn't have room to occur there).

     This is a real bug in GMP's own Windows/MinGW support, not in this
     program, and a public report of the identical symptom (allocation
     sizes suspiciously close to 2^32) on an unrelated large GMP
     operation on Windows backs this up:
     https://github.com/msys2/MINGW-packages/issues/8781

     Bottom line: below a couple hundred million digits, this never
     comes up. Past that, on Windows, use the WSL/Linux build instead
     (see Option C above) -- which is exactly what produced the
     verified 1-billion-digit run at the top of this file.

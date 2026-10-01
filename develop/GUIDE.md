# *DOCS*
## DEVELOPER GUIDE

#### 10/08/2026
mikael.mieskolainen@cern.ch

----

Run repository commands from the repository root after loading the project and Conda environments:

```bash
conda activate graniitti
source install/setenv.sh
```

### Short style guide

- Program internal constants with CAPITALS.
- Model parameters typically lowercase.
- Namespaces with small letters.
- Head namespace is `gra`. Example: `gra::math::pow2`.
- All library classes start with M, as do .h/.cc files, except program files.
- Classes, methods and functions with `PascalStyle` except mathematical such as `pow2`. Short class names
- Hungarian_ for class member variables usually only if confusion is possible.
- STRUCT names with capitals, if contain mostly constant paremeters.
- Short comments, with `[REFERENCE: ]` for the original models / theory.
- Clang-tidy is used for C++ code, max 100 columns per line.
- Ruff is used for Python code, max 100 columns per line.

### Sampling failure policy

- Normal phase space rejection, cuts and zero amplitude support must return a suitable value and must not throw.
- Set `kinematics_ok = false` only for ordinary phase space failure.
- Reserve `PhaseSpaceFailure` for an exceptional internal phase space invariant failure. `MProcess` must catch it, reject the point and mark a technical failure.
- Reserve `KinematicsFailure` for an exceptional failure inside amplitude kinematics, such as nonfinite input, a degenerate projection axis, failed momentum closure or a failed on shell projection.
- Reserve `AmplitudeFailure` for exceptional amplitude failures, such as invalid physical momenta, nonfinite spin or boost normalization, inconsistent topology, helicity or color dimensions, missing required model state, failed color flow construction or inconsistent screening state.
- `MProcess` must convert both amplitude status failures to the generic amplitude exception boundary, reject the event, increment the amplitude and technical failure counters and reset all amplitude and screening state.
- Invalid function arguments and invalid configuration choices may throw `std::invalid_argument` and must not be treated as sampled event failures.
- A physical zero amplitude is not a failure. A recurring failure rate indicates a code, process setup or numerical problem and must be investigated.
- Every amplitude implementation must use the common evaluation status and cleanup interfaces.

### C++ DEBUG and MEMORY diagnostics

```bash
valgrind --tool=callgrind ./(Your binary)
```

Generates `callgrind.out.<pid>`, which can be opened with `kcachegrind`. Compile with `-g` for source annotations, as described in the [Callgrind manual](https://valgrind.org/docs/manual/cl-manual.html).

### Assembler output

#### Create assembler code

```bash
c++ -S -fverbose-asm -g -O2 test.cc -o test.s
```
#### Delete old files

```bash
find . \( -type f -name '*._old*' -o -type d -name '*._old' \) -prune -exec rm -rf -- {} +
```

#### Delete EOS path

```bash
eos root://eoshome-m.cern.ch rm -r /eos/user/m/<path>
```

#### Create asm interlaced with source lines

```bash
as -alhnd test.s > test.lst
```

### Python profiler

```bash
python -m cProfile -o out.profile "myscript.py" --arguments

# Visualize
snakeviz out.profile
```

### Bash commands

#### Search for a variable

```bash
rg -n 'variablename' src
```

Test commands and measurement result tables are documented in [tests/README.md](../tests/README.md).

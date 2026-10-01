#!/usr/bin/env bash
# Install Rivet and its analysis dependencies using the GRANIITTI HepMC3 library
# Run with: bash install/installRIVET.sh
set -eo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$REPO_ROOT/tests/environment.sh"
activate_conda_environment graniitti
source "$REPO_ROOT/install/setenv.sh"
# Isolate Rivet Python modules from external ROOT and other Python installations
unset PYTHONHOME PYTHONPATH

JOBS="${JOBS:-4}"
[[ "$JOBS" =~ ^[1-4]$ ]] || { echo "JOBS must be between 1 and 4" >&2; exit 2; }
[[ -f "$HEPMC3SYS/include/HepMC3/GenEvent.h" ]] || {
  echo "Install HepMC3 first with bash install/installHEPMC3.sh" >&2
  exit 1
}
export INSTALL_PREFIX="${INSTALL_PREFIX:-$REPO_ROOT/.local/rivet}"
export BUILD_PREFIX="${BUILD_PREFIX:-$REPO_ROOT/tmp/rivet-install}"
mkdir -p "$INSTALL_PREFIX/bin" "$BUILD_PREFIX"
INSTALL_PREFIX="$(cd "$INSTALL_PREFIX" && pwd)"
BUILD_PREFIX="$(cd "$BUILD_PREFIX" && pwd)"
# Supply the small LaTeX font sizing package required by Matplotlib when needed
TEX_DIR="$INSTALL_PREFIX/share/texmf/tex/latex/cm-super"
if command -v kpsewhich >/dev/null && ! kpsewhich type1ec.sty >/dev/null && [[ ! -s "$TEX_DIR/type1ec.sty" ]]; then
  mkdir -p "$TEX_DIR"
  curl --fail --location --silent --show-error https://mirrors.ctan.org/fonts/ps-type1/cm-super/type1ec.sty -o "$BUILD_PREFIX/type1ec.sty"
  cp "$BUILD_PREFIX/type1ec.sty" "$TEX_DIR/"
fi
if command -v latex >/dev/null && ! command -v dvipng >/dev/null && [[ ! -x "$INSTALL_PREFIX/bin/dvipng" ]]; then
  command -v apt-get >/dev/null && command -v dpkg-deb >/dev/null || {
    echo "Install dvipng to render Rivet plots with LaTeX" >&2
    exit 1
  }
  (cd "$BUILD_PREFIX" && apt-get download dvipng && dpkg-deb -x dvipng_*.deb "$INSTALL_PREFIX")
  cp "$INSTALL_PREFIX/usr/bin/dvipng" "$INSTALL_PREFIX/bin/"
fi
export MAKE=make MAKEFLAGS="-j$JOBS" HEPMCPATH="$HEPMC3SYS" INSTALL_HEPMC=0
export CFLAGS="${CFLAGS:--O2} -g0" CXXFLAGS="${CXXFLAGS:--O2} -g0"
# YODA files and Rivet plots do not require the optional HDF5 output library
export INSTALL_HDF5="${INSTALL_HDF5:-0}"
if [[ "$INSTALL_HDF5" == 0 ]]; then
  export YODA_CONFFLAGS="${YODA_CONFFLAGS:-} --disable-h5"
fi
# Use Conda zlib headers with the Conda compiler sysroot
export YODA_CONFFLAGS="${YODA_CONFFLAGS:-} --with-zlib=$CONDA_PREFIX"
export RIVET_CONFFLAGS="${RIVET_CONFFLAGS:---disable-static} --with-zlib=$CONDA_PREFIX"
# Select individual analyses when only the repository comparisons are needed
RIVET_ANALYSES="${RIVET_ANALYSES:-all}"
if [[ "$RIVET_ANALYSES" != all ]]; then
  RIVET_CONFFLAGS+=" --disable-analyses"
fi
export PIP_NO_CACHE_DIR=1
cd "$BUILD_PREFIX"
bash "$REPO_ROOT/install/rivet-bootstrap" 2>&1 | tee -a install.log
source "$INSTALL_PREFIX/rivetenv.sh"
if [[ -f "$INSTALL_PREFIX/lib/Rivet/RivetTests.so" ]]; then
  mv "$INSTALL_PREFIX/lib/Rivet/RivetTests.so" "$INSTALL_PREFIX/lib/Rivet/RivetTests.so.$(date +%s)._old"
fi
if [[ "$RIVET_ANALYSES" != all ]]; then
  mkdir -p "$INSTALL_PREFIX/lib/Rivet" "$INSTALL_PREFIX/share/Rivet"
  SOURCES=()
  for name in $RIVET_ANALYSES; do
    files=("$BUILD_PREFIX/Rivet-$(rivet-config --version)"/analyses/plugin*/"$name.cc")
    [[ ${#files[@]} == 1 && -f "${files[0]}" ]] || { echo "Unknown Rivet analysis: $name" >&2; exit 2; }
    SOURCES+=("${files[0]}")
    for file in "${files[0]%.cc}".{cc,info,plot,yoda,yoda.gz}; do
      [[ ! -f "$file" ]] || cp "$file" "$INSTALL_PREFIX/share/Rivet/"
    done
  done
  rivet-build -j "$JOBS" "$INSTALL_PREFIX/lib/Rivet/RivetTests.so" "${SOURCES[@]}"
fi
rivet --version
python -c 'import rivet, yoda'
printf '\nActivate Rivet with: source %q/rivetenv.sh\n' "$INSTALL_PREFIX"

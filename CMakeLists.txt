# GRANIITTI CMake driven compilation
#
# In the root folder, execute:
#
# mkdir build
# cd build
#
# cmake .. (optionally -DWITH_ROOT=ON -DWITH_SYSTEMPATH=ON -DWITH_TEST=ON [see more options below])
# make -j4
# 
# (c) 2017-2026 Mikael Mieskolainen.
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.
# ------------------------------------------------------------- 

cmake_minimum_required(VERSION 3.18 FATAL_ERROR)
project(GRANIITTI LANGUAGES CXX)

# Remove Conda runtime injection when the external C++ ABI is selected
if(DEFINED ENV{CONDA_CXX_ABI} AND "$ENV{CONDA_CXX_ABI}" STREQUAL "0")
  get_filename_component(GRANIITTI_TOOLCHAIN_BIN "${CMAKE_CXX_COMPILER}" DIRECTORY)
  get_filename_component(GRANIITTI_TOOLCHAIN_PREFIX "${GRANIITTI_TOOLCHAIN_BIN}" DIRECTORY)
  set(GRANIITTI_TOOLCHAIN_LIB "${GRANIITTI_TOOLCHAIN_PREFIX}/lib")
  foreach(LINKER_FLAG_VARIABLE
          CMAKE_EXE_LINKER_FLAGS CMAKE_SHARED_LINKER_FLAGS CMAKE_MODULE_LINKER_FLAGS)
    foreach(INJECTED_FLAG
            "-Wl,-rpath,${GRANIITTI_TOOLCHAIN_LIB}"
            "-Wl,-rpath-link,${GRANIITTI_TOOLCHAIN_LIB}"
            "-L${GRANIITTI_TOOLCHAIN_LIB}")
      string(REPLACE "${INJECTED_FLAG}" "" ${LINKER_FLAG_VARIABLE}
                     "${${LINKER_FLAG_VARIABLE}}")
    endforeach()
  endforeach()
endif()

# Use an optimized single-configuration build unless the user selects another type
if(NOT CMAKE_CONFIGURATION_TYPES)
  if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
  endif()
  set_property(
    CACHE CMAKE_BUILD_TYPE
    PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel
  )
endif()

# ------------------------------------------------------------------------------
# Set the installation prefix to the source directory's bin folder
# Running "make -j4" will automatically install executables to graniitti/bin
# ------------------------------------------------------------------------------
set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/bin" CACHE PATH "Installation Directory" FORCE)

# ------------------------------------------------------------------------------
# Global output directories for final targets:
#   - Executables, shared libraries, and static libraries will be placed in "bin"
#     (within the build directory)
# ------------------------------------------------------------------------------
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
file(MAKE_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}" "${CMAKE_INSTALL_PREFIX}")

# ------------------------------------------------------------------------------
# Options
# ------------------------------------------------------------------------------
# Use tri-state strings for optional external frameworks
set(WITH_ROOT "AUTO" CACHE STRING "Enable ROOT-dependent targets (AUTO, ON, or OFF)")
set_property(CACHE WITH_ROOT PROPERTY STRINGS AUTO ON OFF)
set(WITH_LIBTORCH "OFF" CACHE STRING "Enable LIBTORCH-dependent targets (AUTO, ON, or OFF)")
set_property(CACHE WITH_LIBTORCH PROPERTY STRINGS AUTO ON OFF)
set(MARCH "native" CACHE STRING "Target CPU ISA passed to -march")
set_property(CACHE MARCH PROPERTY STRINGS x86-64 x86-64-v2 x86-64-v3 native)

option(WITH_SYSTEMPATH "Include system paths (/usr/include,/usr/local/include)" ON)
option(WITH_TEST "Enable building TEST-dependent targets" OFF)
option(WITH_VALGRIND "Enable VALGRIND flags" OFF)
option(WITH_MADLOOP "Build generated MadLoop process runtimes" ON)

# C++ tests include clang-tidy checks backed by the generated compilation database
if(WITH_TEST)
  find_program(CLANG_TIDY_EXECUTABLE NAMES clang-tidy)
  if(NOT CLANG_TIDY_EXECUTABLE)
    message(FATAL_ERROR
      "WITH_TEST=ON requires clang-tidy; install the dependencies from environment.yml"
    )
  endif()
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
endif()

# Enable safe build accelerators by default
option(WITH_CCACHE "Enable ccache for faster recompiles" ON)
option(WITH_PCH "Enable precompiled headers for common includes" ON)

# Keep unity builds opt-in because merging translation units can violate anonymous namespace isolation
option(WITH_UNITY_BUILD "Enable unity builds for faster compilation" OFF)

if(NOT WITH_ROOT MATCHES "^(AUTO|ON|OFF)$")
  message(FATAL_ERROR "WITH_ROOT must be AUTO, ON, or OFF")
endif()
if(NOT WITH_LIBTORCH MATCHES "^(AUTO|ON|OFF)$")
  message(FATAL_ERROR "WITH_LIBTORCH must be AUTO, ON, or OFF")
endif()
if(MARCH STREQUAL "")
  message(FATAL_ERROR "MARCH must name a compiler-supported -march target")
endif()

# ------------------------------------------------------------------------------
# Print option statuses
# ------------------------------------------------------------------------------
message(STATUS "Configuration options:")
message(STATUS "  CMAKE_BUILD_TYPE:   ${CMAKE_BUILD_TYPE}")
message(STATUS "  WITH_ROOT:         ${WITH_ROOT} <AUTO|ON|OFF>")
message(STATUS "  WITH_TEST:         ${WITH_TEST} <ON|OFF>")
message(STATUS "  WITH_LIBTORCH:     ${WITH_LIBTORCH} <AUTO|ON|OFF>")
message(STATUS "  WITH_SYSTEMPATH:   ${WITH_SYSTEMPATH} <ON|OFF>")
message(STATUS "  WITH_VALGRIND:     ${WITH_VALGRIND} <ON|OFF>")
message(STATUS "  WITH_MADLOOP:      ${WITH_MADLOOP} <ON|OFF>")
message(STATUS "  WITH_CCACHE:       ${WITH_CCACHE} <ON|OFF>")
message(STATUS "  WITH_UNITY_BUILD:  ${WITH_UNITY_BUILD} <ON|OFF>")
message(STATUS "  WITH_PCH:          ${WITH_PCH} <ON|OFF>")
message(STATUS "  MARCH:             ${MARCH}")

# ------------------------------------------------------------------------------
# Graniitti requires C++20 for all core, ROOT and test targets
# ------------------------------------------------------------------------------
if(NOT DEFINED CXX_STANDARD)
  set(CXX_STANDARD 20)
endif()
if(CXX_STANDARD LESS 20)
  message(FATAL_ERROR "Graniitti requires CXX_STANDARD >= 20")
endif()
set(CMAKE_CXX_STANDARD ${CXX_STANDARD})
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ------------------------------------------------------------------------------
# Check for required external library environment variables
# ------------------------------------------------------------------------------
if(NOT DEFINED ENV{HEPMC3SYS})
  message(FATAL_ERROR "Environment variable HEPMC3SYS is not set; please set it (e.g. source ./install/setenv.sh)")
endif()
if(NOT DEFINED ENV{LHAPDFSYS})
  message(FATAL_ERROR "Environment variable LHAPDFSYS is not set; please set it (e.g. source ./install/setenv.sh)")
endif()

set(HEPMC3SYS $ENV{HEPMC3SYS})
set(LHAPDFSYS $ENV{LHAPDFSYS})

# Resolve ROOT from ROOTSYS first and then PATH
set(ROOT_ENABLED OFF)
if(NOT WITH_ROOT STREQUAL "OFF")
  if(DEFINED ENV{ROOTSYS} AND EXISTS "$ENV{ROOTSYS}/bin/root-config")
    set(ROOTSYS "$ENV{ROOTSYS}")
    set(ROOT_CONFIG_EXECUTABLE "${ROOTSYS}/bin/root-config")
  else()
    find_program(ROOT_CONFIG_EXECUTABLE root-config)
  endif()

  if(ROOT_CONFIG_EXECUTABLE)
    execute_process(
      COMMAND ${ROOT_CONFIG_EXECUTABLE} --prefix
      RESULT_VARIABLE ROOT_PREFIX_RESULT
      OUTPUT_VARIABLE ROOTSYS
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(ROOT_PREFIX_RESULT EQUAL 0 AND NOT ROOTSYS STREQUAL "")
      set(ROOT_ENABLED ON)
      message(STATUS "ROOT found at ${ROOTSYS}, enabling ROOT support.")
    endif()
  endif()

  if(NOT ROOT_ENABLED AND WITH_ROOT STREQUAL "ON")
    message(FATAL_ERROR "WITH_ROOT=ON requires a working root-config via ROOTSYS or PATH")
  elseif(NOT ROOT_ENABLED)
    message(STATUS "ROOT not found, disabling optional ROOT support.")
  endif()
endif()

# Resolve AUTO without making PyTorch mandatory
set(LIBTORCH_ENABLED OFF)
if(NOT WITH_LIBTORCH STREQUAL "OFF")
  set(LIBTORCHSYS "")
  set(TORCH_PYTHON_ERROR "")

  if(DEFINED ENV{LIBTORCHSYS} AND NOT "$ENV{LIBTORCHSYS}" STREQUAL "")
    set(LIBTORCHSYS $ENV{LIBTORCHSYS})
  else()
    find_package(Python3 COMPONENTS Interpreter QUIET)
    if(Python3_Interpreter_FOUND)
      execute_process(
        COMMAND ${Python3_EXECUTABLE} -c "import pathlib, torch; print(pathlib.Path(torch.__file__).parent)"
        RESULT_VARIABLE TORCH_PYTHON_RESULT
        OUTPUT_VARIABLE LIBTORCHSYS
        ERROR_VARIABLE TORCH_PYTHON_ERROR
        OUTPUT_STRIP_TRAILING_WHITESPACE
      )
    endif()
  endif()

  if(NOT LIBTORCHSYS STREQUAL "")
    list(APPEND CMAKE_PREFIX_PATH "${LIBTORCHSYS}")
  endif()

  # Conda CUDA toolkits keep headers and nvcc under the target-specific prefix
  if(DEFINED ENV{CONDA_PREFIX} AND
     EXISTS "$ENV{CONDA_PREFIX}/targets/${CMAKE_SYSTEM_PROCESSOR}-linux/include/cuda_runtime.h")
    set(CUDA_TOOLKIT_ROOT_DIR "$ENV{CONDA_PREFIX}/targets/${CMAKE_SYSTEM_PROCESSOR}-linux"
        CACHE PATH "CUDA toolkit root used by PyTorch" FORCE)
    set(CUDAToolkit_ROOT "$ENV{CONDA_PREFIX}/targets/${CMAKE_SYSTEM_PROCESSOR}-linux"
        CACHE PATH "CUDA toolkit root used by PyTorch" FORCE)
  endif()

  # Also search LCG views and other standard CMake prefixes
  find_package(Torch QUIET CONFIG)

  if(Torch_FOUND)
    set(LIBTORCH_ENABLED ON)
    message(STATUS "PyTorch found, enabling LIBTORCH support.")
  elseif(WITH_LIBTORCH STREQUAL "ON")
    message(FATAL_ERROR
      "WITH_LIBTORCH=ON could not find TorchConfig.cmake via LIBTORCHSYS, "
      "Python ${Python3_EXECUTABLE}, or CMAKE_PREFIX_PATH: ${TORCH_PYTHON_ERROR}")
  else()
    message(STATUS
      "PyTorch not found via LIBTORCHSYS, Python ${Python3_EXECUTABLE}, "
      "or CMAKE_PREFIX_PATH; disabling optional LIBTORCH support.")
  endif()
endif()

if(LIBTORCH_ENABLED)
  separate_arguments(TORCH_CXX_OPTIONS NATIVE_COMMAND "${TORCH_CXX_FLAGS}")
endif()

# ------------------------------------------------------------------------------
# Common compiler flags
# MARCH defaults to the build host ISA
# ------------------------------------------------------------------------------
include(CheckCXXCompilerFlag)
unset(GRANIITTI_MARCH_SUPPORTED CACHE)
check_cxx_compiler_flag("-march=${MARCH}" GRANIITTI_MARCH_SUPPORTED)
if(NOT GRANIITTI_MARCH_SUPPORTED)
  message(FATAL_ERROR "The C++ compiler does not support -march=${MARCH}")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fPIC -pipe -march=${MARCH}")
# Keep Debug unoptimized even when an external toolchain injects a baseline -O flag
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
if(WITH_VALGRIND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIE -pg")
endif()

# ------------------------------------------------------------------------------
# Enable ccache if option is set and ccache is found
# ------------------------------------------------------------------------------
if(WITH_CCACHE)
  find_program(CCACHE_FOUND ccache)
  if(CCACHE_FOUND)
    set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_FOUND}")
    message(STATUS "ccache found and enabled for faster recompiles.")
  else()
    message(WARNING "ccache not found; install it for faster incremental builds.")
  endif()
endif()

# ------------------------------------------------------------------------------
# Enable unity builds if option is set (batches sources for faster compilation)
# ------------------------------------------------------------------------------
if(WITH_UNITY_BUILD)
  set(CMAKE_UNITY_BUILD ON)
  set(CMAKE_UNITY_BUILD_BATCH_SIZE 8)  # Adjust batch size as needed; smaller for less memory usage
  message(STATUS "Unity builds enabled for faster compilation.")
endif()

# ------------------------------------------------------------------------------
# Include directories
# ------------------------------------------------------------------------------
include_directories(
  ${HEPMC3SYS}/include
  ${LHAPDFSYS}/include
  ${CMAKE_CURRENT_SOURCE_DIR}                 # Project root
  ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Treat vendored headers as third-party so compiler diagnostics stay focused on GRANIITTI code
include_directories(SYSTEM
  ${CMAKE_CURRENT_SOURCE_DIR}/libs
  ${CMAKE_CURRENT_SOURCE_DIR}/libs/FTensor
  ${CMAKE_CURRENT_SOURCE_DIR}/libs/Eigen
  ${CMAKE_CURRENT_SOURCE_DIR}/libs/Eigen/unsupported
)
if(WITH_SYSTEMPATH)
  include_directories(/usr/include /usr/local/include)
endif()

# ------------------------------------------------------------------------------
# Link directories
# ------------------------------------------------------------------------------
link_directories(
  ${HEPMC3SYS}/lib ${HEPMC3SYS}/lib64
  ${LHAPDFSYS}/lib ${LHAPDFSYS}/lib64
)
if(ROOT_ENABLED)
  link_directories(${ROOTSYS}/lib ${ROOTSYS}/lib/root)
endif()

# ------------------------------------------------------------------------------
# External libraries (specified as lists)
# ------------------------------------------------------------------------------
set(HEPMC3_LIBS -lHepMC3 -lHepMC3search)
set(LHAPDF6_LIBS -lLHAPDF)
set(STANDARD_LIBS -lstdc++ -lm -pthread -lrt)
find_package(ZLIB REQUIRED)

if(LIBTORCH_ENABLED)
  set(LIBTORCH_LIBS ${TORCH_LIBRARIES})
else()
  set(LIBTORCH_LIBS "")
endif()
if(ROOT_ENABLED)
  execute_process(
    COMMAND ${ROOT_CONFIG_EXECUTABLE} --glibs
    OUTPUT_VARIABLE ROOT_GLIBS_RAW
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  separate_arguments(ROOT_LIBS UNIX_COMMAND "${ROOT_GLIBS_RAW}")
  # Harmonic fits use the reentrant Minuit2 interface
  list(APPEND ROOT_LIBS -lMinuit -lMinuit2)
  message(STATUS "ROOT link flags from root-config: ${ROOT_GLIBS_RAW}")
endif()
set(COMMON_LIBS ${HEPMC3_LIBS} ${LHAPDF6_LIBS} ${STANDARD_LIBS}
                ${CMAKE_DL_LIBS} ZLIB::ZLIB)
set(GENERATOR_PROGRAMS gr minbias xscan ampfit)
set(MODEL_TESTS
  test_amplitudes
  test_core
  test_durham
  test_eikonal
  test_form
  test_lbyl
  test_photoproduction
  test_process
  test_regge
  test_multi
  test_soft_model
  test_spin
  test_sudakov
  test_tensor_pomeron
  test_tensor_vector
  test_tensor_spin3
)
set(COMMON_PCH_HEADERS
  <algorithm> <cmath> <complex> <fstream> <iomanip> <iostream> <map>
  <memory> <random> <set> <sstream> <string> <unordered_map> <vector>
)

# Register installation and update the source bin only after relinking
function(register_graniitti_executable target_name)
  install(TARGETS ${target_name} RUNTIME DESTINATION .)
  add_custom_command(TARGET ${target_name} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_INSTALL_PREFIX}
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:${target_name}>
            ${CMAKE_INSTALL_PREFIX}/$<TARGET_FILE_NAME:${target_name}>
    COMMENT "Updating installed ${target_name}"
  )
endfunction()

# Keep ROOT executables on the same ABI selected during configuration
function(configure_root_runtime target_name)
  set_target_properties(${target_name} PROPERTIES
    BUILD_WITH_INSTALL_RPATH TRUE
    INSTALL_RPATH
      "${ROOTSYS}/lib;${ROOTSYS}/lib/root;${HEPMC3SYS}/lib;${HEPMC3SYS}/lib64;${LHAPDFSYS}/lib;${LHAPDFSYS}/lib64"
    INSTALL_RPATH_USE_LINK_PATH FALSE
  )
endfunction()

# ------------------------------------------------------------------------------
# Build GRANIITTI-common static library
# ------------------------------------------------------------------------------
# Keep hand-written sources explicit by prospective C++ module domain
set(SRC_MAIN
  "${CMAKE_CURRENT_SOURCE_DIR}/src/MGlobals.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/MModelTune.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/MUserCuts.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/MUserHistograms.cc"
)
set(SRC_MATH
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Math/MDistribution.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Math/MInterpolation.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Math/MIntegration.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Math/MSpecialFunctions.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Math/MPolarFourier.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Math/MPolarQuadrature.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Math/MPolarCut.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Math/MStatistics.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Math/MTransport.cc"
)
set(SRC_KINEMATICS
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Kinematics/MDecaySampling.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Kinematics/MCentral.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Kinematics/MCollinear.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Kinematics/MEventKinematics.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Kinematics/MFactorized.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Kinematics/MHardDiffraction.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Kinematics/MQuasiElastic.cc"
)
set(SRC_NUCLEAR
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Nuclear/MBeam.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Nuclear/MBreakup.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Nuclear/MFinal.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Nuclear/MEvaporation.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Nuclear/MMass.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Nuclear/MEMD.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Nuclear/MIncoherent.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Nuclear/MConfig.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Nuclear/MDensity.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Nuclear/MFluctuation.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Nuclear/MGlauber.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Nuclear/MGGCF.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Nuclear/MHotSpot.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Nuclear/MNucleus.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Nuclear/MPhoto.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Nuclear/MPhoton.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Nuclear/MRecord.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Nuclear/MUPCScreen.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Nuclear/MSetup.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Nuclear/MShadow.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Nuclear/MSteering.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Nuclear/MUPC.cc"
)
set(SRC_SPIN
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Spin/MHelicityDecay.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Spin/MHelicityNorm.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Spin/MDirac.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Spin/MHelicity.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Spin/MHelicityScatter.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Spin/MPoleLS.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Spin/MSpin.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Spin/MSpinDensity.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Spin/MQMetrics.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Spin/MWigner.cc"
)
set(SRC_ANALYSIS
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Analysis/MSpherical.cc"
)
set(SRC_PARTICLE
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Particle/MForm.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Particle/MPDG.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Particle/MResonance.cc"
)
set(SRC_EIKONAL
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Eikonal/MEikonal.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Eikonal/MProtonScreen.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Eikonal/MEikonalMatrix.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Eikonal/MElasticCNI.cc"
)
set(SRC_PDF
  "${CMAKE_CURRENT_SOURCE_DIR}/src/PDF/MHardPomeronPDF.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/PDF/MLHAPDF.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/PDF/MSudakov.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/PDF/MSudakovIR.cc"
)
set(SRC_REGGE
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Regge/MFragment.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Regge/MReggePhoto.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Regge/MReggeInit.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Regge/MReggeInclusive.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Regge/MReggeMPInit.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Regge/MReggeXPInit.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Regge/MReggeGPInit.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Regge/MRegge.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Regge/MReggeCEP.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Regge/MReggeForm.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Regge/MReggeLoop.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Regge/MReggeMP.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Regge/MReggeMulti.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Regge/MReggeNumerics.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Regge/MReggeParam.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Regge/MReggeSig.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Regge/MReggeSource.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Regge/MReggeMPXP.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Regge/MReggeXP.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Regge/MReggeGP.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Regge/MSoftModel.cc"
)
set(SRC_TENSOR
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Tensor/MTensorExchange.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Tensor/MTensorForward.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Tensor/MTensorInit.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Tensor/MTensorPhoto.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Tensor/MTensorPomeron.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Tensor/MTensorParam.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Tensor/MTensorVector.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Tensor/MTensorSpin3.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Tensor/MTensorProcess.cc"
)
set(SRC_PHOTON
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Photon/MPhotoDiss.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Photon/MFlux.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Photon/MGamma.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Photon/MPhotoQCD.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Photon/MPhotoVM.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Photon/MPhotoZ.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Photon/MQED.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Photon/MRadiative.cc"
)
set(SRC_QCD
  "${CMAKE_CURRENT_SOURCE_DIR}/src/QCD/MDurham.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/QCD/MPartonFlavour.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/QCD/MPartonProposal.cc"
)
set(SRC_PROCESS
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Process/MAmpMatch.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Process/MDecayTree.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Process/MEventRecord.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Process/MForwardExcitation.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Process/MHelicityConfig.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Process/MProcess.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Process/MProcessConfig.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Process/MProcessScreening.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Process/MProcessState.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Process/MProcessTable.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Process/MSubProc.cc"
)
set(SRC_SAMPLING
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Sampling/MRandom.cc"
)
set(SRC_TECH
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Tech/MAux.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Tech/MJson.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Tech/MH1.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Tech/MH2.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Tech/MLHE.cc"
)
set(SRC_CORE
  ${SRC_MAIN}
  ${SRC_ANALYSIS}
  ${SRC_MATH}
  ${SRC_KINEMATICS}
  ${SRC_NUCLEAR}
  ${SRC_SPIN}
  ${SRC_PARTICLE}
  ${SRC_EIKONAL}
  ${SRC_PDF}
  ${SRC_REGGE}
  ${SRC_TENSOR}
  ${SRC_PHOTON}
  ${SRC_QCD}
  ${SRC_PROCESS}
  ${SRC_SAMPLING}
  ${SRC_TECH}
)

# Discover only converter-owned generated MG5 translation units recursively
file(GLOB_RECURSE SRC_AMPLITUDE CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Amplitude/MG5/*.cc"
)
# Discover native analytic amplitude translation units recursively
file(GLOB_RECURSE SRC_ANALYTIC_AMPLITUDE CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Amplitude/Photon/*.cc"
)
# Retained converter outputs are never active translation units
list(FILTER SRC_AMPLITUDE EXCLUDE REGEX "/[^/]*\\._old[^/]*/")
set(SRC_GENERATOR
  "${CMAKE_CURRENT_SOURCE_DIR}/src/MGraniitti.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Sampling/MNeuroJac.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Sampling/MNeuroJacLoss.cc"
)
set(SRC_COMMON ${SRC_CORE} ${SRC_AMPLITUDE} ${SRC_ANALYTIC_AMPLITUDE})

# Build each converter-owned MadLoop source bundle once per build tree
set(MADLOOP_TARGETS)
set(MADLOOP_WRAPPERS)
if(WITH_MADLOOP)
  find_program(MADLOOP_FORTRAN_COMPILER NAMES gfortran)
  if(MADLOOP_FORTRAN_COMPILER)
    find_package(Python3 COMPONENTS Interpreter REQUIRED)
    file(GLOB_RECURSE MADLOOP_SOURCE_DATA CONFIGURE_DEPENDS
      "${CMAKE_CURRENT_SOURCE_DIR}/MG5cards/*/*/source/source.json"
    )
    foreach(MADLOOP_DATA ${MADLOOP_SOURCE_DATA})
      get_filename_component(MADLOOP_SOURCE ${MADLOOP_DATA} DIRECTORY)
      get_filename_component(MADLOOP_OWNER ${MADLOOP_SOURCE} DIRECTORY)
      get_filename_component(MADLOOP_NAME ${MADLOOP_OWNER} NAME)
      get_filename_component(MADLOOP_PROJECTION_ROOT ${MADLOOP_OWNER} DIRECTORY)
      get_filename_component(MADLOOP_PROJECTION ${MADLOOP_PROJECTION_ROOT} NAME)
      set(MADLOOP_WORK "${CMAKE_CURRENT_BINARY_DIR}/MadLoop/${MADLOOP_NAME}/work")
      set(MADLOOP_RUNTIME "${CMAKE_CURRENT_BINARY_DIR}/MadLoop/${MADLOOP_NAME}/runtime")
      set(MADLOOP_OUTPUTS
        "${MADLOOP_RUNTIME}/runtime.json"
        "${MADLOOP_RUNTIME}/libMadLoop.so"
        "${MADLOOP_RUNTIME}/libMadLoopAdapter.so"
      )
      add_custom_command(
        OUTPUT ${MADLOOP_OUTPUTS}
        COMMAND ${Python3_EXECUTABLE} -B
          ${CMAKE_CURRENT_SOURCE_DIR}/develop/MG2GRA/modules/madloop.py bundle
          --source ${MADLOOP_SOURCE}
          --work ${MADLOOP_WORK}
          --output ${MADLOOP_RUNTIME}
          --compiler ${MADLOOP_FORTRAN_COMPILER}
          --cxx ${CMAKE_CXX_COMPILER}
        DEPENDS
          ${CMAKE_CURRENT_SOURCE_DIR}/develop/MG2GRA/modules/madloop.py
          ${MADLOOP_DATA}
          ${MADLOOP_SOURCE}/madloop.tar.gz
          ${MADLOOP_SOURCE}/cuttools.tar.gz
          ${MADLOOP_SOURCE}/ninja.tar.gz
          ${MADLOOP_SOURCE}/oneloop.tar.gz
          ${MADLOOP_SOURCE}/MADGRAPH_LICENSE
        COMMENT "Building MadLoop runtime ${MADLOOP_NAME}"
        VERBATIM
      )
      set(MADLOOP_TARGET "madloop_${MADLOOP_NAME}")
      add_custom_target(${MADLOOP_TARGET} DEPENDS ${MADLOOP_OUTPUTS})
      set(MADLOOP_CARD_TARGET "madloop_${MADLOOP_NAME}_card")
      add_custom_target(${MADLOOP_CARD_TARGET}
        COMMAND ${CMAKE_COMMAND} -E make_directory
          "${MADLOOP_RUNTIME}/MadLoop5_resources"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
          "${MADLOOP_OWNER}/param_card.dat"
          "${MADLOOP_RUNTIME}/param_card.dat"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
          "${MADLOOP_OWNER}/param_card.dat"
          "${MADLOOP_RUNTIME}/MadLoop5_resources/param_card.dat"
        DEPENDS "${MADLOOP_OWNER}/param_card.dat"
        COMMENT "Syncing MadLoop parameter card ${MADLOOP_NAME}"
        VERBATIM
      )
      add_dependencies(${MADLOOP_CARD_TARGET} ${MADLOOP_TARGET})
      list(APPEND MADLOOP_TARGETS ${MADLOOP_CARD_TARGET})
      list(APPEND MADLOOP_WRAPPERS
        "${CMAKE_CURRENT_SOURCE_DIR}/src/Amplitude/MG5/${MADLOOP_PROJECTION}/AMP_MG5_${MADLOOP_NAME}.cc"
      )
    endforeach()
  else()
    message(WARNING "No gfortran compiler was found, disabling MadLoop runtimes.")
  endif()
endif()

# Suppress harmless warnings produced by generated MadGraph amplitude code only
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set_source_files_properties(${SRC_AMPLITUDE} PROPERTIES
    COMPILE_FLAGS "-Wno-unused-variable -Wno-int-in-bool-context"
  )
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set_source_files_properties(${SRC_AMPLITUDE} PROPERTIES
    COMPILE_FLAGS "-Wno-unused-variable"
  )
endif()

# Compile common sources directly into the static library
add_library(graniitti_common STATIC ${SRC_COMMON})
target_include_directories(graniitti_common PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
if(MADLOOP_TARGETS)
  add_dependencies(graniitti_common ${MADLOOP_TARGETS})
  set_property(SOURCE ${MADLOOP_WRAPPERS} APPEND PROPERTY COMPILE_DEFINITIONS
    GRANIITTI_MADLOOP_ROOT="${CMAKE_CURRENT_BINARY_DIR}/MadLoop"
  )
endif()

# Build one common precompiled header and explicitly reuse it downstream
if(WITH_PCH)
  target_precompile_headers(graniitti_common PRIVATE
    ${COMMON_PCH_HEADERS}
  )
  message(STATUS "Precompiled headers enabled for graniitti_common.")
endif()

# Compile generator and neural-adaptation code separately from ROOT analysis consumers
add_library(graniitti_generator STATIC ${SRC_GENERATOR})
target_include_directories(graniitti_generator PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
if(LIBTORCH_ENABLED)
  target_include_directories(graniitti_generator SYSTEM PRIVATE ${TORCH_INCLUDE_DIRS})
  target_compile_options(graniitti_generator PRIVATE ${TORCH_CXX_OPTIONS})
  target_compile_definitions(graniitti_generator PRIVATE GRANIITTI_USE_LIBTORCH=1)
endif()
if(WITH_PCH)
  target_precompile_headers(graniitti_generator PRIVATE ${COMMON_PCH_HEADERS})
  message(STATUS "Precompiled headers enabled for graniitti_generator.")
endif()

# ------------------------------------------------------------------------------
# Build ROOT-dependent static library (if enabled)
# ------------------------------------------------------------------------------
if(ROOT_ENABLED)
  file(GLOB SRC_ROOT
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Analysis/*.cc"
  )
  list(FILTER SRC_ROOT EXCLUDE REGEX "/MSpherical\\.cc$")
  add_library(graniitti_root STATIC ${SRC_ROOT})
  target_include_directories(graniitti_root PUBLIC ${ROOTSYS}/include ${ROOTSYS}/include/root)

  # ROOT consumers reuse one project-standard precompiled header
  if(WITH_PCH)
    target_precompile_headers(graniitti_root PRIVATE ${COMMON_PCH_HEADERS})
    message(STATUS "Precompiled headers enabled for graniitti_root.")
  endif()
endif()

# ------------------------------------------------------------------------------
# Build CATCH3-TEST-dependent static library (if enabled)
# ------------------------------------------------------------------------------
if(WITH_TEST)
  enable_testing()
  set(SRC_TEST "${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp/support/catchmain.cc")
  add_library(graniitti_test STATIC ${SRC_TEST})
  target_include_directories(graniitti_test PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

  # Reuse the already-built common PCH for the Catch entry point
  if(WITH_PCH)
    target_precompile_headers(graniitti_test REUSE_FROM graniitti_common)
    message(STATUS "Precompiled headers enabled for graniitti_test (reusing from common).")
  endif()
endif()

# ------------------------------------------------------------------------------
# Build executables from src/Program/
# ------------------------------------------------------------------------------
file(GLOB SRC_PROGRAM "${CMAKE_CURRENT_SOURCE_DIR}/src/Program/*.cc")
foreach(srcfile ${SRC_PROGRAM})
  get_filename_component(exename ${srcfile} NAME_WE)
  add_executable(${exename} ${srcfile})
  if(exename IN_LIST GENERATOR_PROGRAMS)
    target_link_libraries(${exename}
      graniitti_generator graniitti_common graniitti_common
      ${COMMON_LIBS} ${LIBTORCH_LIBS}
    )
  else()
    target_link_libraries(${exename}
      graniitti_common graniitti_common ${COMMON_LIBS}
    )
  endif()

  # Reuse precompiled headers if enabled
  if(WITH_PCH)
    target_precompile_headers(${exename} REUSE_FROM graniitti_common)
  endif()

  register_graniitti_executable(${exename})
endforeach()

# ------------------------------------------------------------------------------
# Build ROOT-dependent executables from src/Program/Analysis/ (if enabled)
# ------------------------------------------------------------------------------
if(ROOT_ENABLED)
  file(GLOB SRC_PROGRAM_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/src/Program/Analysis/*.cc")
  foreach(srcfile ${SRC_PROGRAM_ROOT})
    get_filename_component(exename ${srcfile} NAME_WE)
    add_executable(${exename} ${srcfile})
    target_include_directories(${exename} PRIVATE ${ROOTSYS}/include ${ROOTSYS}/include/root)
    target_link_libraries(${exename}
      graniitti_root graniitti_common graniitti_common ${COMMON_LIBS} ${ROOT_LIBS}
    )
    configure_root_runtime(${exename})

    # Reuse the ROOT-compatible PCH
    if(WITH_PCH)
      target_precompile_headers(${exename} REUSE_FROM graniitti_root)
    endif()

    register_graniitti_executable(${exename})
  endforeach()
endif()

# ------------------------------------------------------------------------------
# Build test executables (if enabled)
# ------------------------------------------------------------------------------
if(WITH_TEST)
  file(GLOB CPP_TEST_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp/test_*.cc"
  )
  set(CPP_ROOT_TESTS test_analyzer test_harmonic_analysis)
  set(CPP_INTEGRATION_TESTS ${MODEL_TESTS} test_hard_pomeron test_nuclear_processes test_sampling)
  set(CPP_BENCHMARKS test_neurojac)


  foreach(srcfile ${CPP_TEST_SOURCES})
    get_filename_component(exename ${srcfile} NAME_WE)
    if(exename IN_LIST CPP_ROOT_TESTS AND NOT ROOT_ENABLED)
      continue()
    endif()

    add_executable(${exename} ${srcfile})
    if(exename IN_LIST CPP_ROOT_TESTS)
      target_include_directories(${exename} PRIVATE ${ROOTSYS}/include ${ROOTSYS}/include/root)
      target_link_libraries(${exename}
        graniitti_root graniitti_common graniitti_test graniitti_common
        ${COMMON_LIBS} ${ROOT_LIBS}
      )
      configure_root_runtime(${exename})
    elseif(exename IN_LIST CPP_BENCHMARKS)
      target_link_libraries(${exename}
        graniitti_generator graniitti_common ${COMMON_LIBS} ${LIBTORCH_LIBS}
      )
    else()
      # All non-ROOT tests can use the generator and shared model helpers
      target_link_libraries(${exename}
        graniitti_generator graniitti_common graniitti_test graniitti_common
        ${COMMON_LIBS} ${LIBTORCH_LIBS}
      )
    endif()

    # Reuse precompiled headers if enabled
    if(WITH_PCH)
      if(exename IN_LIST CPP_ROOT_TESTS)
        target_precompile_headers(${exename} REUSE_FROM graniitti_root)
      else()
        target_precompile_headers(${exename} REUSE_FROM graniitti_common)
      endif()
    endif()

    register_graniitti_executable(${exename})

    if(exename IN_LIST CPP_INTEGRATION_TESTS)
      set(test_label "integration")
    elseif(exename IN_LIST CPP_BENCHMARKS)
      set(test_label "benchmark")
    else()
      set(test_label "unit")
    endif()
    add_test(NAME "cpp.${exename}" COMMAND $<TARGET_FILE:${exename}>)
    set_tests_properties("cpp.${exename}" PROPERTIES
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      LABELS "${test_label}"
    )
  endforeach()

  # Check hand-written generator and command-line translation units in parallel through CTest
  set(CPP_CLANG_TIDY_SOURCES ${SRC_CORE} ${SRC_GENERATOR} ${SRC_ANALYTIC_AMPLITUDE} ${SRC_PROGRAM})
  foreach(srcfile IN LISTS CPP_CLANG_TIDY_SOURCES)
    file(RELATIVE_PATH lint_source "${CMAKE_CURRENT_SOURCE_DIR}" "${srcfile}")
    string(REPLACE "/" "." lint_name "${lint_source}")
    string(REGEX REPLACE "\\.cc$" "" lint_name "${lint_name}")
    add_test(
      NAME "cpp.clang_tidy.${lint_name}"
      COMMAND
        ${CLANG_TIDY_EXECUTABLE}
        --config-file=${CMAKE_CURRENT_SOURCE_DIR}/.clang-tidy
        --quiet
        -p=${CMAKE_BINARY_DIR}
        ${srcfile}
    )
    set_tests_properties("cpp.clang_tidy.${lint_name}" PROPERTIES
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      LABELS "unit;lint"
    )
  endforeach()
endif()

# ------------------------------------------------------------------------------
# Add custom target: superclean
# This target cleans the contents of the bin folder (without removing the folder itself)
# ------------------------------------------------------------------------------
add_custom_target(superclean
  COMMAND find ${CMAKE_BINARY_DIR}/bin -mindepth 1 -delete
  COMMENT "Cleaning contents of bin folder"
)

# ------------------------------------------------------------------------------
# Final message
# ------------------------------------------------------------------------------
message(STATUS "
  Makefile creation via CMake complete.
  
  For faster builds, consider using Ninja: cmake -G Ninja ..
  Then build with: ninja -j4
  
  Compile: make -j4 (or ninja -j4).
  Install: updated automatically for each relinked executable
  
  Final binaries are placed in: ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
  Installation prefix is set to: ${CMAKE_INSTALL_PREFIX}")

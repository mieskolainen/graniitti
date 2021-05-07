# See: https://anaconda.org/conda-forge/gcc_linux-64
#
# Show all: ls $CONDA_PREFIX/bin

conda install -c conda-forge gcc_linux-64 
conda install -c conda-forge gxx_linux-64 

# Set compilers before compiling anything
export CC=$CONDA_PREFIX/bin/x86_64-conda_cos6-linux-gnu-gcc
export CXX=$CONDA_PREFIX/bin/x86_64-conda_cos6-linux-gnu-g++

# Compile HepMC3 and LHAPDF
# cd install && source autoinstall.sh && cd ..

# Compile
# make -j8 CXX=$CXX CXX_STANDARD=c++17 TEST=TRUE

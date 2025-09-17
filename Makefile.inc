# PYTHIA configuration file.
# Generated on Sun Sep 14 07:30:10 UTC 2025 with the following command:
# ./configure --prefix=/usr/local --with-hepmc3=/usr/local --with-lhapdf6=/usr/local --with-rivet=/usr/local --with-gzip

# Install directory prefixes.
PREFIX_BIN=/usr/local/bin
PREFIX_INCLUDE=/usr/local/include
PREFIX_LIB=/usr/local/lib
PREFIX_SHARE=/usr/local/share/Pythia8

# Compilation flags (see ./configure --help for further documentation).
CXX=g++
CXX_DTAGS=-Wl,--disable-new-dtags
CXX_COMMON=-O2 -std=c++11 -pedantic -W -Wall -Wshadow -fPIC -pthread -DGZIP
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname,
LIB_SUFFIX=.so
OBJ_COMMON=

EVTGEN_USE=false
EVTGEN_CONFIG=
EVTGEN_BIN=
EVTGEN_INCLUDE=
EVTGEN_LIB=

FASTJET3_USE=false
FASTJET3_CONFIG=
FASTJET3_BIN=
FASTJET3_INCLUDE=
FASTJET3_LIB=

HEPMC2_USE=false
HEPMC2_CONFIG=
HEPMC2_BIN=
HEPMC2_INCLUDE=
HEPMC2_LIB=

HEPMC3_USE=true
HEPMC3_CONFIG=HepMC3-config
HEPMC3_BIN=/usr/local/bin/
HEPMC3_INCLUDE=-I/usr/local/include
HEPMC3_LIB=-L/usr/local/lib -Wl,-rpath,/usr/local/lib -lHepMC3

LHAPDF5_USE=false
LHAPDF5_CONFIG=
LHAPDF5_BIN=
LHAPDF5_INCLUDE=
LHAPDF5_LIB=

LHAPDF6_USE=true
LHAPDF6_CONFIG=lhapdf-config
LHAPDF6_BIN=/usr/local/bin/
LHAPDF6_INCLUDE=-I/usr/local/include
LHAPDF6_LIB=-L/usr/local/lib -Wl,-rpath,/usr/local/lib -lLHAPDF

POWHEG_USE=false
POWHEG_CONFIG=
POWHEG_BIN=
POWHEG_INCLUDE=
POWHEG_LIB=

RIVET_USE=true
RIVET_CONFIG=rivet-config
RIVET_BIN=/usr/local/bin/
RIVET_INCLUDE=-I/usr/local/include
RIVET_LIB=-L/usr/local/lib -Wl,-rpath,/usr/local/lib -lRivet

YODA_USE=false
YODA_CONFIG=
YODA_BIN=
YODA_INCLUDE=
YODA_LIB=

ROOT_USE=false
ROOT_CONFIG=
ROOT_BIN=
ROOT_INCLUDE=
ROOT_LIB=

GZIP_USE=true
GZIP_CONFIG=
GZIP_BIN=
GZIP_INCLUDE=
GZIP_LIB=-lz

PYTHON_USE=false
PYTHON_CONFIG=
PYTHON_BIN=
PYTHON_INCLUDE=
PYTHON_LIB=

MG5MES_USE=false
MG5MES_CONFIG=
MG5MES_BIN=
MG5MES_INCLUDE=
MG5MES_LIB=

OPENMP_USE=false
OPENMP_CONFIG=
OPENMP_BIN=
OPENMP_INCLUDE=
OPENMP_LIB=

MPICH_USE=false
MPICH_CONFIG=
MPICH_BIN=
MPICH_INCLUDE=
MPICH_LIB=

HDF5_USE=false
HDF5_CONFIG=
HDF5_BIN=
HDF5_INCLUDE=
HDF5_LIB=

HIGHFIVE_USE=false
HIGHFIVE_CONFIG=
HIGHFIVE_BIN=
HIGHFIVE_INCLUDE=
HIGHFIVE_LIB=

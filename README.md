# Jets + Pythia

## Environment
Easiest way to get the environment is via Apptainer (Singularity) container:

```bash
apptainer pull docker://hepstore/rivet-pythia:main
```

Command to run the container:

```bash
apptainer exec -B /gpfs01 rivet-pythia_main.sif bash
```

## Compile + Local Run
```bash
make
./makeTree
```

Parameters can be tuned in `makeTree.cc`
```cpp
   const double R = 0.4;
   const double jetEtaMax = 1.0 - R;
   const double dPhiMin = 0.75 * M_PI; // back-to-back requirement
   const double jetPtMin = 3.0;
   // particle parameters
   const double partPtMin = 0.15;
   const double partEtaMax = 1.0;
```

## Using Batchfarm (HTCondor)
```bash
./submit/submit_all.sh
```
- It will submit jobs using definied ptHat bins in `submit/ptHatBins.list`
- The output will be stored in `submit/output/`
- After jobs are done, the script will merge trees (Histogram `stats` contains `cross section` and `nEvents`, which are additive)
- And execute analysis macro `anaTrees/anaTrees.cpp+`


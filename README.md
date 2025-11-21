# Jet–Medium Interaction Simulation Using Color-Tracking LBT

This repository provides a full simulation pipeline for studying **jet–medium interactions in the quark–gluon plasma (QGP)** using a **color-tracking Linear Boltzmann Transport (LBT)** model.

## Introduction

In this framework, partons from medium are categorized into **positive partons and negative partons**:

### **Positive Partons**
Positive partons include:
- **Jet shower partons** produced in the initial hard scattering  
- **Recoil partons** generated when QGP constituents receive momentum from jet–medium elastic or inelastic interactions  

These partons represent the *actual physical energy and momentum flow* and are propagated through the medium during the LBT simulation.

### **Negative Partons**
Negative partons represent the **energy–momentum holes** left in the QGP when medium partons are kicked out as recoil partons.  
They carry:
- The same magnitude of momentum as the recoil partons  
- The *opposite sign*, representing loss of local medium energy  

Negative partons **must be subtracted during final analysis** to ensure:
- **Energy–momentum conservation**
- Proper treatment of medium backreaction
- Correct reconstruction of final hadron spectra and jet observables  

This positive–negative scheme is essential for physically consistent modeling of jet-induced medium response.

---

## Repository Structure

```
.
├── PYTHIA8/             # PYTHIA8 event generator
├── LBT/                 # Linear Boltzmann Transport model with color tracking
├── lpposi_Frag/         # Hadronization module (positive parton)
│   └── ALICE-cjet/      # Analysis Module for jets spectrum
├── lpnega_Frag/         # Hadronization module (negative parton)
├── pp_lpposi_Frag/      # Hadronization for pp baseline
│   └── ALICE-cjet/      # Analysis Module for jets spectrum
├── RAA-D0-PbPb5020/     # Analysis Module for hadrons spectrum
├── hydrofile/           # Prepared hydrofiles for different collisions
└── README.md
```

---

## Features
- Full jet → QGP → hadron simulation chain  
- Color-tracking during jet transport  
- Configurable medium and jet parameters  
- Independent compilation of all modules  
- Suitable for heavy-ion jet-quenching phenomenology  
---

## Important Note on Paths

Most folders (**PYTHIA**, **lpposi_Frag**, **lpnega_Frag**, **pp_lpposi_Frag**) contain their own **compile** file.  
These "compile" files may include **hard-coded paths** to:

- PYTHIA / fastjet installation  
- Local include / lib directories  

⚠️ **Users MUST edit these paths to match their own environment before compiling.**

---

## Requirements
- `g++`
- `cmake`  
- **PYTHIA8** 
- **fastjet**
- **fjcontrib**

---

## Installation & Build

### **1. Event Generation (PYTHIA)**

```
cd PYTHIA8/
./compile.sh
./driver.sh
```

---

### **2. Jet–Medium Transport (LBT)**

```
cd LBT/
rm -rf build
./compile
./driver.sh
```

This reads PYTHIA output and produces medium-modified partons with color tracking.

---

### **3. Hadronization Stage**

Each hadronization module must be compiled separately:

```
cd lpposi_Frag/
./compile.sh
./driver.sh
```
This reads LBT output and produces hadrons.

Repeat similarly for:

- `lpnega_Frag/`
- `pp_lpposi_Frag/`

---

## Usage Workflow

```
1. Run PYTHIA
2. Run LBT
3. Run hadronization modules
```

---

## Outputs

- **PYTHIA** → initial parton list: `jet_shower_parton.dat`
 cross section information: `jet_sigmainfo.dat`
- **LBT** → medium-modified partons with color information: `lp-posi.dat` `lp-nega.dat`
- **Hadronization** → final hadron list: `hadrons-posi.dat` `hadrons-nega.dat`

---

## Hydrodynamic Profile (QGP Background)
To simulate QGP produced in different collision systems or centralities, you must change the hydrodynamic input files used by LBT.

The required hydrofiles typically include:

- `bulk.dat`
- `mc_glauber.dat`

You can find prepared hydrofiles for different collisions in the folder `hydrofile/`.

To apply a specific QGP background:

1. Choose the corresponding hydrofile set from `hydrofile/`.
2. Copy the files(e.g. `bulk.dat`, `mc_glauber.dat`) into: `LBT/hydroProfile/`

---

## Improving Efficiency with PYTHIA Event Selection
To improve overall program efficiency, users are encouraged to enable event selection in PYTHIA. By selecting only relevant events (e.g. based on kinematic cuts or production of specific jets), the simulation can avoid generating unnecessary events, significantly reducing runtime and computational cost without affecting the physics goals. 

**Current example in this repository:** The `PYTHIA` folder is configured to select only events that contain inclusive jets within a given jet momentum (pT) and rapidity range, serving as a practical example of event pre-selection to improve efficiency

---

## pTHat Range Splitting for Efficiency

If necessary, users can further improve efficiency by generating events in different ranges of transverse momentum exchange using:
```
PhaseSpace:pTHatMin
PhaseSpace:pTHatMax
```
Each pTHat interval can be simulated separately. The final physical results should then be obtained by summing contributions from all pTHat ranges, with each range properly scaled by its corresponding cross section. This approach preserves physical correctness while significantly reducing computational waste in regions of low relevance.

The pTHat range is configured in the PYTHIA command file:
```
PYTHIA8/m01.cmnd
```

---

## Known Limitations
- Modules build independently  
- Input/output formats must match predefined structures  
- Medium profile must be manually configured  

---

## Citation
If you use this code for published work, please cite the related LBT papers and associated publications.

---

## Contact
For issues, questions, or contributions, please contact the maintainer or open a GitHub issue.

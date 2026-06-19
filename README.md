# BACHELOR-S-THESIS-FVG-nsPipeFlow-Unsteady-Driving

This repository contains the MATLAB scripts developed for the post-processing and generation of asymmetric Reynolds-number cycles used in the Bachelor's Thesis:

**Numerical study on the variation of heat transfer in a turbulent pipe flow with asymmetric flow-rate cycles featuring slow acceleration and rapid deceleration**

The work focuses on turbulent pipe-flow simulations subjected to asymmetric flow-rate cycles, with particular attention to the analysis of heat-transfer response through the Nusselt number and the skin-friction coefficient.

## Repository scope

This repository only includes the MATLAB scripts used to:

- generate asymmetric Reynolds-number profiles;
- normalize the corresponding bulk-velocity forcing signal;
- post-process Nusselt-number time series;
- compare DNS-based Nusselt-number data with the Gnielinski correlation;
- post-process skin-friction time series.

The full DNS solver is not included in this repository. Only the MATLAB scripts developed for the thesis workflow are provided.

## Repository structure

```text
matlab/
├── cycle_generation/
│   ├── demo_generate_Re_profile_asym_raisedcos_double_hold.m
│   └── generate_Re_profile_ndim_multicycle_asym_raisedcos_doble_hold.m
│
├── nusselt_postprocessing/
│   ├── process_nusselt_timeseries.m
│   └── compare_nusselt_gnielinski_steady.m
│
└── friction_postprocessing/
    └── process_friction_timeseries.m
````

## MATLAB scripts

### Cycle generation

The folder `matlab/cycle_generation/` contains the scripts used to generate smooth asymmetric Reynolds-number profiles based on raised-cosine transitions.

The generated cycle is composed of four stages:

1. lower holding stage;
2. smooth acceleration;
3. upper holding stage;
4. smooth deceleration.

The main function writes the normalized forcing signal required by the simulation workflow.

Generated files may include:

```text
Profile_Data.txt
Reynolds_Average.txt
Resumen_Caso.txt
```

These files are considered output files and are not intended to be version-controlled.

### Nusselt-number post-processing

The folder `matlab/nusselt_postprocessing/` contains scripts for processing Nusselt-number timeseries obtained from DNS simulations.

The scripts allow:

* cleaning and reconstruction of time series;
* removal of initial transient cycles;
* cycle-based averaging;
* comparison against steady DNS reference values;
* comparison against the Gnielinski correlation.

### -friction post-processing

The folder `matlab/friction_postprocessing/` contains the script used to process the friction timeseries.

The script performs:

* time-series reconstruction;
* interpolation over a uniform time grid;
* detection of the steady cyclic regime;
* cycle-based averaging of the skin-friction coefficient;
* optional sign correction during rapid deceleration stages.

## Basic usage

Clone the repository:

```bash
git clone https://github.com/vanogonzalezfrancisco-commits/BACHELOR-S-THESIS-FVG-nsPipeFlow-Unsteady-Driving-.git
```

Open MATLAB and add the repository folders to the MATLAB path:

```matlab
addpath(genpath('matlab'))
```

To generate an asymmetric Reynolds-number profile, run:

```matlab
demo_generate_Re_profile_asym_raisedcos_double_hold
```

To process Nusselt-number data, use:

```matlab
out = process_nusselt_timeseries(filename);
```

To process skin-friction data, use:

```matlab
out = process_friction_timeseries(filename);
```

where `filename` is the corresponding time-series file exported from the simulation workflow.

## Output files

The scripts may generate text files, figures, or intermediate data files. These files are intentionally excluded from version control through `.gitignore`.

Typical generated files include:

```text
Profile_Data.txt
Reynolds_Average.txt
Resumen_Caso.txt
```

Large simulation outputs and local post-processing results should not be committed to the repository.

## Requirements

The scripts are written in MATLAB and are intended for post-processing numerical simulation data from turbulent pipe-flow cases.

The code has been organized for clarity and reproducibility within the context of the thesis workflow.

## Author

**Francisco Vañó González**

Bachelor's Thesis
Grado en Ingeniería en Tecnologías Industriales
Universidad de Málaga

## License

This repository is intended for academic and research purposes.

If a `LICENSE` file is included in the repository, the use of this code is governed by the terms specified in that file.

```
```

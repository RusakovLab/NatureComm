# NatureComm
Computational Model of Glutamate Diffusion and NMDA Receptor Activation
# Glutamate Diffusion and NMDA Receptor Modeling

This repository contains the MATLAB code and associated files used for the computational modeling of glutamate diffusion in the extracellular space and its interaction with NMDA receptors. The results of this study have been submitted to *Nature Communications*. The model simulates the diffusion of glutamate molecules released from astrocytes, their interaction with surrounding structures, and the subsequent activation of NMDA receptors.

## Overview of the Model

The model consists of several key components:

* **Glutamate Diffusion Simulation:** The `FixBallsAstrogliaSurfaceRelease.m` script simulates the diffusion of glutamate molecules in a 3D space filled with spherical obstacles representing astrocytes and other cellular structures. The script uses Brownian motion to model the random movement of glutamate molecules and accounts for adhesion to astrocyte surfaces.
* **Input Parameters:** The `InputParametersSR.m` file defines the parameters for the simulation, including the number of particles, the size of the computational domain, and the properties of the astrocytes and adhesion zones.
* **NMDA Receptor Dynamics:** The `NMDA.m` and `NMDA_SpaceSR.m` scripts model the dynamics of NMDA receptor activation in response to glutamate binding. The model includes binding, unbinding, desensitization, and resensitization processes, based on kinetic parameters from the literature.
* **Data Visualization:** The `PlotDataControl.m` script is used to visualize the spatial distribution of glutamate molecules and the state of NMDA receptors over time. It generates 3D scatter plots and contour plots to illustrate the results of the simulation.
* **Statistical Analysis:** The `statisticSR.txt` file contains the parameters used for the statistical analysis of the simulation results, including the number of trials and the probability of adhesion.

## Key Features

* **3D Diffusion Simulation:** The model simulates glutamate diffusion in a realistic 3D environment with obstacles representing astrocytes and other cellular structures.
* **Adhesion and Uptake:** The model includes mechanisms for glutamate adhesion to astrocyte surfaces and uptake by transporters, with stochastic processes governing the binding and unbinding events.
* **NMDA Receptor Activation:** The model captures the complex dynamics of NMDA receptor activation, including the effects of glutamate concentration and receptor desensitization.
* **Visualization Tools:** The repository includes scripts for visualizing the spatial distribution of glutamate and the state of NMDA receptors, providing insights into the spatiotemporal dynamics of the system.

## Usage

To run the simulation, follow these steps:

1.  **Set Parameters:** Modify the parameters in `InputParametersSR.m` and `statisticSR.txt` to match your experimental setup.
2.  **Run Simulation:** Execute `FixBallsAstrogliaSurfaceRelease.m` to simulate glutamate diffusion and adhesion.
3.  **Analyze NMDA Dynamics:** Use `NMDA_SpaceSR.m` to model NMDA receptor activation based on the glutamate distribution generated in the previous step.
4.  **Visualize Results:** Use `PlotDataControl.m` to generate 3D plots and contour maps of the simulation results.

## Files Included

* `FixBallsAstrogliaSurfaceRelease.m`: Main script for simulating glutamate diffusion.
* `InputParametersSR.m`: Defines input parameters for the simulation.
* `NMDA.m`: Models NMDA receptor dynamics.
* `NMDA_SpaceSR.m`: Simulates NMDA receptor activation based on glutamate distribution.
* `PlotDataControl.m`: Script for visualizing simulation results.
* `statisticSR.txt`: Contains statistical parameters for the simulation.

## Citation

If you use this code in your research, please cite the associated paper submitted to *Nature Communications*.

# Extended KKS Code using Phase-Field Method : A Phase-field model for growth and coarsening of Si precipitate in AlSi10Mg SLM in a super-saturated matrix

## Overview
This framework extends the  KKS model using the Phase-Field method to reproduce DSC tests. It serves as a software realization of the research described in [the asoociated paper], focusing on [growth and coarseing of precipitates in PBF AlSi10Mg alloy]. This implementation includes various computational tools and data to replicate and expand upon the paper's findings.

## related Article (under consideration)
- Extension of a KKS model to predict the microstructure evolution in LPBF AlSi10Mg alloy submitted to non isothermal processes
- more details will be available 

## Code Structure
The project is organized into two main parts, focusing on different aspects of the KKS code implementation, and consists of three primary Jupyter notebooks:

### Jupyter Notebooks
1. **Main Simulation Code**
   - `pyFi_origin.ipynb`: This notebook contains the core KKS code used for running the simulations. It is designed to handle the main computational processes including the evolution of microstructures.

2. **Post-Processing Tools**
   - `Postprocess.ipynb`: This notebook is utilized for the post-processing of simulation outputs. It deals with:
     - Mapping of microstructures.
     - Analysis of Differential Scanning Calorimetry (DSC) curves.
     - Qualitative analysis including microstructure evolution.

3. **DSC Analysis Tool**
   - `compute_DSC.ipynb`: Dedicated to computing quantitative data from simulations, specifically focusing on DSC heat curves.

### Additional Details
Each notebook utilizes output data which includes various simulation details like:
- Concentrations of phases.
- Order parameters.
- Energy residuals.
These information are captured at regular intervals throughout the simulation process, enabling a thorough analysis of both the physical properties and the computational performance.

Note: The maps allow also the continuation/restart/bakcup of the simulation.

  

### Data and Scripts
- **Data Files**: 
  - `micro_final.npz`: Contains the Calphad data used within the simulations.
  - Initialization data are stored in a text file for easy access and manipulation.
- **Script Files**: 
  - Python script files (`*.py`) contain the various routines required to run the code. These scripts are automatically loaded by the main Jupyter notebook.

### Saved Outputs
The output maps are saved in associated folders, organized by:
- `big_grids/`: Contains maps used in the article for content representation.
- `small_grids/`: Contains maps for sensitivity analysis.
- `nucleation/`: Stores outputs related to the injection of nucleons within the initialization phase.

## Installation Requirements
- **Python Version**: Python 3.7 or more recent.
- **Dependencies**: The code uses modern Python libraries and acceleration and parallelization packages such as:
  - `numba`
  - `pyfftw`
  - Modern coding techniques like vectorization and broadcasting are extensively employed.

## Usage
To use this framework:
1. Ensure that you have Python 3.7 or higher installed along with the required libraries.
2. Clone this repository to your local machine.
3. Open the `pyFi_origin.ipynb` to run the KKS code or `Postprocess.ipynb` for post-processing tasks.

## Licensing
You are free to use, edit, and distribute this code under [Apache License v2.0]. 

## Contributors
- **Code Developers**: [Jocelyn Delahaye and Seifallah El Fetni]

---

## Supplementary video:
- Simulation of the microsctructure evolution of AlSi10Mg SLM material during anisothermal heat load : matrix desaturation and growth and coalescence of d-Si precipitates. Further post-processing allows the obtention of the DSC curves
 as detailled in the related article and achieved by the post-process tool.
- ![](https://github.com/SFETNI/KKS-Phase-Field-code/blob/main/save_fig/Supp_Mat_computed_microstructure_evolution.gif)

For more information on how to setup and run the code, refer to the detailed comments within each script and notebook. If you encounter any issues or have questions, feel free to open an issue in this repository.

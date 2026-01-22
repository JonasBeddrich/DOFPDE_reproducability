# A kernel compression method for distributed-order fractional partial differential equations 

This repository is the code base of the manuscript: 
A kernel compression method for distributed-order fractional partial differential equations by Beddrich and Wohlmuth (2025) available [here, reference will be updated](https://arxiv.org/abs/2508.13631). 

## Project dependencies

Core dependencies are:

- Julia (1.11.2 or above)
- Python (3.8.10 or above)
- [CMake](https://cmake.org/) (3.19 or above)
- [MFEM](https://mfem.org/) (4.8 or above)
- [Hypre](https://computing.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods) (necessary for MFEM)
- [Metis](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) (necessary for MFEM)

See also [BUILD-MFEM.md](./BUILD-MFEM.md) on how to install the dependencies.

## Build instructions
1. Clone the repository
   ```bash
   git clone https://github.com/JonasBeddrich/DOFPDE_reproducability.git
   ```
2. Create a build directory
   ```bash
   mkdir build && cd build
   ```
3. Run CMake and point it to the source directory
   ```bash
   cmake ..
   ```
4. Build the project
   ```bash
   make -j
   ```

## Code structure 
- Exponential sum approximations (Julia) 
    - **ExpSum:** scripts and visualization
- Test cases for distributed-order fractional differential equations (Python)
    - **DOFDE:** scripts and visualization    
- Implementations for distributed-order fractional partial differential equations (MFEM)
    - **applications:** scripts to run simulations
    - **data:** preprocessing (mesh creation, parameter fields) 
    - **Images:** postprocessing (visualization) 
    - **src:** code base of the implementations of the differential operators

## Usage

### 3. Exponential sum approximation (**ExpSum**) 

- Run computeRADOKernels.jl to
    - Compute rational approximations
    - Transform to exponential sum approximations
    - Compute L1 errors
- Visualization using plots_RA_DO.ipynb

### 5.1 Influence of the kernel approximation (**DOFDE**)

- Run DOFDE_Example1.ipynb

### 5.2 Non-zero initial conditions (**DOFDE**)

- Run DOFDE_Example2.ipynb

### 5.3 Distributed-order diffusion-wave equation 

#### 5.3.1 Convergence Study 

- Run executable convergence_study_dodiffwaveRIIA2
- Temporal discretization has to be chosen in dodiffusionwaveoperatorRIIA2
    - Implicit Euler
    - An L-stable DIRK method
    - Radau IIA scheme  

#### 5.3.2 Lake Starnberg 
*Attention: Due to the spatial resolution, the simulations require significant computational resources.*

- Create mesh with data/Starnbergersee/create_mesh.ipynb
- Run executable DOwave_Starnbergersee
    * Scenario 0: $\phi_{2,0.1}$
    * Scenario 1: $\phi_{2,0.5}$
    * Scenario 2: $\phi(\alpha) = 1$ (not in the manuscript)
- Move .pvd files into Images/DOwave_Starnbergersee/data
- Visualization using Graphics_Starnbergersee.py 

#### 5.3.3 Space-dependent operator 
*Attention: Due to the spatial resolution, the simulations require significant computational resources.*

- Create spatial dependencies with data/Marmousi/marmousi.ipynb
- Run executable DOwave_Marmousi 
    * Scenario 0: Marmousi
    * Scenario 1: Random field 
    * Scenario 2: Geometric
- Move .pvd files into Images/DOwave_Marmousi/data
- Visualization using Graphics_Marmousi.ipynb 

### 5.4 Nonlinear elasticity with DO damping 

- Run executable compression
    * Scenario 0: $\tau = 0.1$ 
    * Scenario 1: $\tau = 0.5$ (not in the manuscript)
    * Scenario 2: $\tau = 2$
    * Scenario 3: $\phi_{1,0.1}$ (not in the manuscript)
- Visualization using Paraview 

## Developers

- [Jonas Beddrich](mailto:jonas.beddrich@tum.de)

## License

Copyright (c) 2026 Jonas Beddrich

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
documentation files (the “Software”), to deal in the Software without restriction, including without limitation the
rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit
persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

See also [online (MIT License)](https://opensource.org/license/mit/)

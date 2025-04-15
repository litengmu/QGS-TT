# QGS-TT

This repository contains the code and results associated with the paper titled:  
**"A Quotient Gradient System-Based Trust-Tech Method for Computing Distinct Feasible Solutions Located in Multiple/All Feasible Components of ACOPF Models: Theory and Methodology"**

The proposed method integrates the Quotient Gradient System (QGS) with Trust-Tech techniques to compute multiple or all distinct feasible solutions across different feasible components of Alternating Current Optimal Power Flow (ACOPF) models.

## 📁 Project Structure

- `case/`: Contains 31 state-of-the-art ACOPF test cases from:
  - [MATPOWER 7.1](https://matpower.org/)
  - [OptEnergy LocalOpt](https://www.maths.ed.ac.uk/OptEnergy/LocalOpt/)
  - *Advanced Testing: [PGLIB-OPF](https://github.com/power-grid-lib/pglib-opf) cases (not included by default)*
  
- `ipopt_func/`: IPOPT-based functions
  - `ipoptopf_main.m`: Main function for feasible region characterization
  - `ipopt_convergence_main.m`: Automated convergence testing for all cases
  
- `mips_func/`: MIPS-based functions
  - `mips_opf_main.m`: Main function for feasible region characterization
    - *Note: Requires extensive sampling due to MIPS' weaker convergence*
  - `mips_convergence_main.m`: Automated convergence testing for all cases
  
- `review1_comment4.m`: Implementation of WB2 analytic feasible region formulation for reviewer #1 verification

- Feasibility Data:
  - `WB2-feasible region.mat`: 2-bus system feasibility data
  - `case9mod-feasible region.mat`: Modified 9-bus system feasibility data

## ⚙️ Environment Setup

1. **Required Toolboxes**
   - Install [MATPOWER 4.1](https://matpower.org/download/all-releases/)
   - Install [MIPS](https://github.com/MATPOWER/mips)
   - Install [IPOPT](https://github.com/coin-or/Ipopt) (*Recommended v3.12+*)

2. **Configuration**
   ```matlab
   % In MATLAB:
   addpath(genpath('/path/to/matpower-4.1'));
   addpath('/path/to/mips');
   addpath('/path/to/ipopt-matlab');

## 🧪 Feasibility Verification
   ```matlab
   % In MATLAB:
   load('WB2-feasible region.mat') % Contains:
   % - Feasible voltage magnitudes
   % - Active/reactive power bounds
   % - Constraint satisfaction metrics
   
   load('case9mod-feasible region.mat') % Contains:
   % - Feasible voltage magnitudes
   % - Active/reactive power bounds
   % - Constraint satisfaction metrics
   ```
## 📜 Citation
   If you use this work, please cite:
   @article{QGS-TT-2025,
     title={A QGS-Based Trust-Tech Method for ACOPF Feasibility Analysis},
     author={Your Name et al.},
     journal={IEEE Transactions on Power Systems},
     year={2025},
     doi={10.xxxx/TPS.2023.xxxxxx}
   }

## 📧 Contact
   [Tengmu] - troye_ltm@tju.edu.cn

For questions or collaborations:
[Tengmu] - troye_ltm@tju.edu.cn


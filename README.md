# QGS-TT

This repository contains the code and results associated with the paper titled:  
**"A Quotient Gradient System-Based Trust-Tech Method for Computing Distinct Feasible Solutions Located in Multiple/All Feasible Components of ACOPF Models: Theory and Methodology"**

The proposed method integrates the Quotient Gradient System (QGS) with Trust-Tech techniques to compute multiple or all distinct feasible solutions across different feasible components of Alternating Current Optimal Power Flow (ACOPF) models.

## 📁 Project Structure

- `case/`: Test case data files.
- `ipopt_func/`: Functions integrating the IPOPT solver.
- `mips_func/`: Functions integrating the MIPS solver.
- `review1_comment4.m`: MATLAB script addressing reviewer comments.
- `WB2-feasible region.mat` & `case9mod-feasible region.mat`: Feasible region data files for specific test cases.

## ⚙️ Environment Setup

1. Install [MATPOWER 4.1](https://matpower.org/download/all-releases/) and add it to your MATLAB path.
2. Ensure that both IPOPT and MIPS solvers are installed and accessible within your MATLAB environment.

## 🚀 Usage
1. Clone this repository:
   ```bash
   git clone https://github.com/litengmu/QGS-TT.git
2. run('review1_comment4.m');

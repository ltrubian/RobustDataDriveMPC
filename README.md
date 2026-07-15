# Robust Data Driven Model Predictive Control

> **Course:** Learning Dynamical Systems – Constrol System Engineering
> **University:** University of Padua
> **Academic Year:** 2025/2026
> **Authors:** Lorenzo Favaron, Luca Galeazzo, Lorenzo Trubian
> **Professor:** Mattia Zorzi

## Abstract [TODO]
*Briefly summarize the project in 3-5 sentences. State the problem addressed, the methodology used (e.g., simulation, theoretical analysis, experimental setup), and the key results obtained. Mention if this project replicates a specific paper or solves a specific exercise.*

## Objectives
Implement the one-optimizer scheme for the optimal selection of the tolerance $c$ and propose:
- Implementation a MPC with the LFM as nominal model using MATLAB.
- Comparison the performance of the new controller vs MPC with just the predicion of the filters as starting point.

## Theoretical Background [TODO]
*Briefly outline the core theoretical concepts required to understand the project. This section demonstrates your grasp of the course material.*
- **Key Equations:**  State-space representation of a system
    $$x(k+1) = A x(k) + Bv(k)+ K u(k) \\ y(k) = C x(k) + Dv(k)$$
    Idea to implement
    $$\min_{c, u, \hat{y}, \hat{x}} \sum_{t=-L}^{0} \beta^t\| y_t - \hat{y}_{t|t-1}(c) \|^2 + J(x(k), \mathbf{u}_k, k) $$
    QP problem solving the MPC
    $$J(x(k), \mathbf{u}_k, k) = \sum_{i=0}^{N-1} \left[ \|x(k+i)\|^2_{C^\top Q C} + \|u(k+i)\|^2_R\right] + \|x(k+N)\|^2_{C^\top P_f C}$$
    Mention the governing equations (e.g., State-space representation, Transfer functions).
- **Assumptions:** List any simplifications made (e.g., linearization, neglecting friction).
- **References:** B. Levy and R. Nikoukhah, "Robust state-space filtering under incremental
    model perturbations subject to a relative entropy tolerance," *IEEE Trans.
    Autom. Control*, vol. 58, pp. 682–695, Mar. 2013

## Methodology & Implementation
### System Architecture [TODO]
*Describe the workflow. If applicable, include a block diagram (you can embed your TikZ image here or a link to it).*
![General block scheme](./docs/general_block_scheme.png)

### Tools & Dependencies
- **Language:** MATLAB R2025b
- **Libraries:** Control Systems Toolbox, Optimization Toolbox, OSQP (optional)

### Algorithm Description [TODO]
*Step-by-step logic of your implementation:*
1. Data preprocessing / Model definition.
2. Controller design / Simulation setup.
3. Optimization / Solving process.
4. Data analysis.

## Results
### Simulation / Experimental Setup [TODO]
*Describe the parameters used (e.g., sampling time, initial conditions, input signals).*

### Key Findings
*Present the main outcomes. Use bullet points for clarity.*
- **Stability:** The system is [stable/unstable] with a margin of...
- **Performance:** Rise time = [X]s, Overshoot = [Y]%.
- **Comparison:** Method A showed a [Z]% improvement over Method B in terms of...

### Visualizations [TODO]
*Embed your most important plots here.*
*Figure 1: Step response of the closed-loop system.*

## How to Reproduce [TODO]
*Instructions for running your code or compiling your report. This is crucial for academic reproducibility.*

### Prerequisites (optional)
Execute the **./install_osqp.m** script to install the OSQP solver, which is faster than the default MATLAB one *quadprog*

## License
This project is submitted for academic purposes only.
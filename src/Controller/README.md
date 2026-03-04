# Different Controllers

COMMON INTERFACE: NO

## MPC with outer uncertainty estimation
The controller estimate the optimal uncertainty *c* and THEN a MPC compute the controller

- status:       IMPLEMENTED
- simulated:    YES (first_simulation.m)
- file:         mpc_optimizer.m
- helper files:
  - estimate_uncertainty.m  (minimizing estimation error)


## MPC with inner uncertainty estimation
The controller estimate the optimal uncertainty using the MPC optimizer itself

- status:       NOT IMPLEMENTED
- simulated:    NO ()
- file:         Controller.m
- helper files:

# Different Robust Kalman Filter implementations

COMMON INTERFACE: NO

## rkinteration
- origin:       prof's package
- helper file:
- tested:       TRUSTED, other implementations are tested against it


## robust_kalman_step
- origin:       oracle
- helper file:  calculate_gamma.m
- tested:       NO


## RobustKalmanFilter
- origin:       handmade
- helper file:  LagrangeMultiplier.m
- tested:       YES


## Other files

### maxtol.m
This file computes the maximum tolerance a system can handle.
It comes from the prof's MATLAB package

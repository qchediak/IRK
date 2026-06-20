# Starting Approximations for Implicit Runge-Kutta Methods Applied to Ordinary Differential Equations and Differential-Algebraic Equations

__Author:__ Q. Chediak

The MATLAB files in this repository illustrate and verify the analytical results of the PhD dissertation listed below and an associated publication currently in preparation.

The code includes implementations of implicit Runge-Kutta (IRK) methods and partitioned Runge-Kutta (PRK) methods used to solve autonomous ordinary differential equations (ODEs) and index 3 differential-algebraic equations (DAEs); specifically, autonomous semi-explicit index 3 DAEs in Hessenberg form.  It also contains numerical experiments investigating the order of starting approximations for the nonlinear systems arising in these methods.

This code includes implementations of Newton's method and the simplified Newton method for solving nonlinear systems. As a third option, it can interface with `nsold`, a nonlinear solver developed by C. T. Kelley, through the `kelley` option.  Users wishing to use `nsold` should obtain it separately.

## Repository Structure

- `solve_ODE`: Solves ODEs by an IRK method.

- `solve_pODE`: Solves partitioned ODEs by a PRK method.

- `solve_index3`: Solves index 3 DAEs by (implicit) PRK methods.  This works for several preset problems and could easily be modified to work for an arbitrary problem of the correct form.

- `starting_approx_order_ODE`: Tests the order of starting approximations in the nonlinear system of the IRK method applied to an ODE.

- `starting_approx_order_pODE`: Tests the order of starting approximations in the nonlinear system of the PRK method applied to a partitioned ODE.

- `starting_approx_order_DAE`: Tests the order of starting approximations in the nonlinear system of the PRK method applied to a DAE.

- `shared_functions`: Functions used by multiple directories.

## Running the Code

To run any of the experiments, enter a directory and run the command `main`.  To change the parameters, edit the `init` file.

## Associated Work

This repository accompanies the PhD dissertation _Starting Approximations for Implicit Runge-Kutta Methods Applied to Ordinary Differential Equations and Differential-Algebraic Equations_ by Q. Chediak (2026).

## License

See the `LICENSE.md` file for licensing information.


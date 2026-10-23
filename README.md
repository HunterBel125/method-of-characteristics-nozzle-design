# Method of Characteristics Nozzle Design

**Deliverable:** A MATLAB program that uses the method of characteristics to generate a two-dimensional supersonic nozzle contour, visualize the flow field, evaluate operating conditions, and study numerical refinement.

**Completion:** This project was completed as a two-person gas dynamics coursework project at the University of Arizona.

**Objective:** Develop a numerical nozzle design tool for a target exit Mach number of 6. The project combines characteristic calculations, nozzle wall construction, and flow visualization with an analysis of pressure, temperature, mass flow, and sensitivity to the number of characteristics.

## Part 1: Nozzle Geometry and Flow Visualization

[`part1_nozzle_geometry.m`](part1_nozzle_geometry.m) calculates the nozzle contour and characteristic network using the method of characteristics.

The code uses Prandtl–Meyer relations and characteristic compatibility equations to determine local flow angles and Mach numbers. Centerline symmetry and wall tangency conditions determine the characteristic intersections and nozzle wall coordinates. The resulting Mach values are interpolated to produce a contour map of the flow field.

Outputs include:

- Nozzle wall contour and initial throat arc.
- Both characteristic families and their intersections.
- Mach number contour map.
- Final wall angle, exit Mach values, nozzle length, and exit-to-throat height ratio.

This file contains the shared numerical solver used by Parts 2 and 3.

## Part 2: Operating Conditions and Performance

[`part2_operating_conditions.m`](part2_operating_conditions.m) uses the nozzle solution to evaluate flow properties and operating conditions.

Isentropic relations are used to calculate pressure and temperature distributions from the local Mach numbers. The geometry is scaled using a specified exit height and nozzle width. The code then evaluates choked mass flow over a range of stagnation pressures and calculates the original project’s momentum-flow and wall-pressure-force estimates.

Outputs include:

- Centerline Mach number, pressure ratio, and temperature ratio.
- Pressure and temperature contour maps.
- Dimensional throat height and nozzle length.
- Stagnation pressure required to match the specified exit back pressure.
- Stagnation temperature required to satisfy the assumed exit-temperature criterion.
- Mass flow as a function of stagnation pressure.
- Momentum-flow and wall-pressure-force estimates.

## Part 3: Numerical Refinement

[`part3_grid_convergence.m`](part3_grid_convergence.m) repeats the nozzle calculation using different numbers of initial characteristics.

The computed exit-to-throat height ratio is compared with the theoretical isentropic area–Mach relation. The code records the absolute difference, exit height, and nozzle length for each discretization, then plots the area-ratio error against the number of characteristics.

This study examines the sensitivity of the calculated geometry to characteristic resolution and provides a numerical consistency check. It does not constitute experimental validation.

## Teamwork and Contributions

We divided the project into functional components, reviewed and debugged each other’s work, and integrated our contributions into a completed MATLAB program.

My primary contributions focused on generating the nozzle contour and characteristic network and creating the Mach number contour map.

The three files in this repository were reconstructed from the completed program to present its components clearly. They illustrate the project’s organization rather than preserve the original intermediate versions.

## Running the Code

Save all three `.m` files in the same folder and set that folder as MATLAB’s Current Folder. Run any component by entering:

```matlab
part1_nozzle_geometry
part2_operating_conditions
part3_grid_convergence
```

Parts 2 and 3 automatically call the solver in Part 1 with its plotting disabled. You do not need to run Part 1 beforehand.

To retrieve the numerical solution without generating figures:

```matlab
result = part1_nozzle_geometry(50, false);
```

The first argument specifies the number of initial characteristics. The second controls whether Part 1 produces its plots.

Edit the design inputs in Part 1 to change the target Mach number or gas specific heat ratio. Operating inputs are defined in Part 2, and the refinement counts are defined in Part 3.

## Assumptions and Interpretation

The model assumes steady, two-dimensional, inviscid flow of a calorically perfect gas, with a default specific heat ratio of 1.4. It does not account for boundary layers, viscous losses, real-gas effects, or three-dimensional sidewall effects.

The following details are relevant when interpreting the retained project calculations:

- The pressure calculated by matching exit pressure to back pressure is an ideal design condition, not a separately calculated nozzle-starting limit.
- The quantity labeled `thrust` in the code is `mdot × V_exit`. It represents exit momentum flow and does not include a complete thrust balance.
- The wall-force estimate integrates internal absolute pressure; it does not subtract external pressure.
- The geometry covers the centerline-to-wall region. The intended full-height or half-height convention should be established before interpreting dimensional areas and mass flow.
- The minimum exit temperature of 55 K is an assumed project criterion; the code does not model condensation.

## Skills Applied

MATLAB programming, compressible-flow analysis, method of characteristics, numerical root finding, flow-field interpolation, numerical integration, refinement studies, collaborative debugging, and technical visualization.

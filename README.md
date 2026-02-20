Stellar Remnant Orbital Evolution in a Milky Way Potential
Overview

This project simulates the dynamical evolution of massive stars in a Milky Way–like gravitational potential, including stellar evolution, supernova kicks, and remnant orbital redistribution.

The simulation integrates:

Stellar evolution using SEVN

Galactic dynamics using a composite Milky Way potential

Supernova natal kicks

Symplectic leapfrog orbital integration

Final spatial distribution analysis of neutron stars and black holes

The final output is a CSV file containing initial and final spatial coordinates of compact remnants under different kick prescriptions and metallicities.

Scientific Goal

The code aims to study:

How supernova kicks affect the spatial distribution of compact remnants

The influence of metallicity on remnant populations

The long-term dynamical evolution of neutron stars and black holes in a realistic galactic potential

This is particularly useful for:

Compact object population synthesis studies

Galactic structure modeling

Gravitational wave progenitor environment studies

Remnant retention analysis

Physical Model
Milky Way Potential

The galaxy model consists of:

NFW Dark Matter Halo

Miyamoto–Nagai Disk

Truncated Power-Law Bulge

These components are combined using MultiPotential.

All quantities are converted to N-body units using:

Mass scale: 
10
10
𝑀
⊙
10
10
M
⊙
	​


Length scale: 1000 pc

Stellar Population

Stars are initialized with:

Exponential radial disk distribution (Rd = 2.6 kpc)

Exponential vertical distribution (zd = 0.3 kpc)

Circular velocities derived from local gravitational acceleration

Kroupa Initial Mass Function (10–150 M☉)

Each star is evolved individually using SEVN.

Supernova Kicks

Two kick configurations are supported:

hobbs

hobbs_pure

Kick velocity dispersion:

265 km/s

Natal kicks are applied when a remnant forms:

RemnantType 5 → Neutron Star

RemnantType 6 → Black Hole

Numerical Integration

Orbital integration is performed using:

Symplectic Leapfrog Integrator

Half-step velocity update
Full-step position update
Final half-step velocity update

This ensures:

Good energy conservation

Stability over long integration times

Suitable for collisionless galactic dynamics

Time step:

dt = 0.001 (N-body units)
Code Structure
1. MilkyWayPotential

Defines the composite galactic potential.

2. generate_stellar_population()

Generates initial positions and velocities

Samples stellar masses from Kroupa IMF

Initializes SEVN Star objects

3. evolve_stars()

Evolves stars physically using SEVN

Updates masses dynamically

Applies SN kicks

Integrates motion during stellar evolution

4. integrate_orbits()

Continues orbital evolution after stellar evolution

Integrates positions up to final simulation time

5. analyze_results()

Identifies neutron stars and black holes

Extracts final positions

6. save_results_to_csv()

Saves:

Initial coordinates

Final coordinates

Remnant type

Metallicity

Kick prescription

Kick dispersion

Output

The simulation produces:

stellar_simulation_results.csv

Columns:

Column	Description
Type	Neutron Star / Black Hole
Initial X,Y,Z	Initial galactocentric position
Final X,Y,Z	Final position after evolution
Metallicity	Stellar metallicity
SN Kicks	Kick prescription
Kick Std	Kick velocity dispersion
Requirements

You must install:

numpy

pandas

sevnpy

fireworks (nbodylib)

scipy (optional but recommended)

Example:

pip install numpy pandas scipy

SEVN and fireworks must be installed from their respective repositories.

How to Run

Simply execute:

python your_script_name.py

The simulation will:

Loop over kick configurations

Loop over metallicities

Generate 100 stars per configuration

Evolve stars and orbits

Save results to CSV

Adjustable Parameters

You may modify:

N_stars = 100
metallicities = [0.02, 0.002]
kick dispersion
integration time
disk scale lengths
IMF limits

For large-scale studies, increase N_stars.

Limitations

No self-gravity between stars (test particle approximation)

Fixed galactic potential (static MW model)

No dynamical friction

No binary evolution

Single-star evolution only

Possible Extensions

Include binary population synthesis

Add time-varying galactic potential

Include dynamical friction for massive BHs

Add spiral arms perturbations

Perform statistical analysis of radial distributions

Compare scale heights between NS and BH populations

Add metallicity-dependent IMF

Scientific Interpretation

This framework allows you to explore:

Radial migration of compact remnants

Vertical heating of neutron stars

Retention fraction in galactic disk

Metallicity effects on compact object formation

It can be extended to connect with:

Gravitational wave source distributions

Pulsar population synthesis

Galactic compact object demographics

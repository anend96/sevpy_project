# Stellar Remnant Orbital Evolution in a Milky Way Potential

## Overview

This project simulates the dynamical evolution of massive stars in a Milky Way–like gravitational potential, including stellar evolution, supernova natal kicks, and long-term orbital redistribution of compact remnants.

The simulation combines:

- Stellar evolution (via SEVN)
- Composite Milky Way gravitational potential
- Supernova kick prescriptions
- Symplectic leapfrog orbital integration
- Final remnant spatial distribution analysis

The final output is a CSV file containing the initial and final galactocentric positions of neutron stars and black holes under different metallicities and kick models.

---

## Scientific Goal

This project aims to investigate:

- How supernova natal kicks reshape the spatial distribution of compact remnants
- The impact of metallicity on neutron star and black hole populations
- Long-term orbital evolution in a realistic galactic potential

Applications include:

- Compact object population synthesis  
- Galactic structure modeling  
- Gravitational wave progenitor environment studies  
- Remnant retention analysis  

---

## Physical Model

### Milky Way Potential

The galaxy model includes:

1. NFW Dark Matter Halo  
2. Miyamoto–Nagai Disk  
3. Truncated Power-Law Bulge  

These components are combined using a `MultiPotential` framework.

All quantities are converted to N-body units:

- Mass scale: \( 10^{10} M_\odot \)  
- Length scale: 1000 pc  

---

### Stellar Population

Stars are initialized with:

- Exponential radial disk profile (Rd = 2.6 kpc)  
- Exponential vertical scale height (zd = 0.3 kpc)  
- Circular velocities from local gravitational acceleration  
- Kroupa Initial Mass Function (10–150 \(M_\odot\))  

Each star is evolved individually using SEVN.

---

## Supernova Kicks

Supported kick prescriptions:

- `hobbs`
- `hobbs_pure`

Kick velocity dispersion:

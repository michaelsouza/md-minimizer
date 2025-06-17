<img src="https://r2cdn.perplexity.ai/pplx-full-logo-primary-dark%402x.png" class="logo" width="120"/>

# Using Breakable Harmonic Bonds in LAMMPS

## Overview

In LAMMPS, to use a breakable harmonic bond style, you actually combine the standard harmonic bond style with additional commands or packages that allow bond breaking, since the plain `bond_style harmonic` does not support bond breaking by itself.

## Using `bond_style harmonic`

### Definition

The `harmonic` bond style defines a bond potential with:

* Equilibrium distance $r_0$
* Force constant $K$

### Syntax

```lammps
bond_style harmonic
bond_coeff 1 K r0
```

This style does not inherently support bond breaking; it models bonds as permanent harmonic springs[^6][^8].

## Using Breakable Bonds in LAMMPS

### Available Options

For breakable bonds, LAMMPS offers specialized bond styles and fixes:

* `bond_style bpm/rotational` or `bpm/spring` from the BPM package
* `bond_style quartic` which breaks bonds when stretched beyond a cutoff
* `fix bond/break` command that can dynamically break bonds during simulation

## Combining Harmonic Bonds with Breakability

### Method 1: Using `fix bond/break`

There is no direct `harmonic/breakable` bond style in LAMMPS. To simulate harmonic bonds that can break:

1. Use `bond_style harmonic` for the bond potential
2. Use `fix bond/break` to specify breaking criteria

### Method 2: Using Breakable Bond Styles

Alternatively, use breakable bond styles:

```lammps
bond_style bpm/rotational break yes
bond_coeff ...
```

## Summary

* Use `bond_style harmonic` for standard harmonic bonds (no breakage)
* For breakable harmonic-like bonds:
	+ Use `bpm/rotational` or `quartic` bond styles
	+ Or combine `harmonic` with `fix bond/break`

<div style="text-align: center">⁂</div>

## References

[^1]: [LAMMPS Bond Styles](https://docs.lammps.org/bond_style.html)
[^2]: [fix bond/break](https://docs.lammps.org/fix_bond_break.html)
[^3]: [LAMMPS Documentation](https://www.afs.enea.it/software/lammps/doc17/html/bond_style.html)
[^4]: [Broken Bonds Howto](https://docs.lammps.org/Howto_broken_bonds.html)
[^5]: [BPM Rotational Bonds](https://docs.lammps.org/bond_bpm_rotational.html)
[^6]: [Harmonic Bonds](https://docs.lammps.org/bond_harmonic.html)
[^7]: [Quartic Bonds](https://gensoft.pasteur.fr/docs/lammps/2020.03.03/bond_quartic.html)
[^8]: [Harmonic Bonds (legacy)](https://www.afs.enea.it/software/lammps/doc17/html/bond_harmonic.html)
[^9]: [fix bond/break (legacy)](https://www.afs.enea.it/software/lammps/doc17/html/fix_bond_break.html)

---

## Package Requirements

The `fix bond/break` command requires the **MC package**, which must be explicitly installed. Similarly, BPM package styles require separate installation.

## compute count/type command

### Syntax

```lammps
compute ID group-ID count/type mode
```

### Modes

* atom
* bond
* angle
* dihedral
* improper

### Examples

```lammps
compute 1 all count/type atom
compute 1 flowmols count/type bond
```

Added in version 15Jun2023. Counts:

* Atoms per atom type
* Bonds per bond type
* Angles per angle type
* Dihedrals per dihedral type
* Impropers per improper type

For mode = bond, broken bonds (type 0) are also counted.

### Output

Calculates:

* Global vector of counts
* For bond mode, also a global scalar count of broken bonds

Values are "intensive".
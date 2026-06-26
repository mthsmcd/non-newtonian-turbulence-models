# non-newtonian-turbulence-models
OpenFOAM implementation of Reynolds-averaged turbulence models for non-Newtonian fluids.
Developed and tested in [OpenFOAM v2306]().

The four models currently available were originally proposed for power-law (PL), Bingham (Bn) 
and Herschel-Bulkley (HB) fluids in the [reference papers](#references), and investigated in our paper:

- **_"Turbulence modeling for viscoplastic fluids: state-of-the-art review and open-source investigation comparing with new DNS"_
 Matheus S. S. Macedo, Gustavo E. O. Celis, Hamidreza Anbarlooei, Roney L. Thompson,
Daiane M. Iceri, Letícia Bizarre, Vanessa Richon and Marcelo S. Castro,
accepted for publishing at the _Journal of Non-Newtonian Fluid Mechanics_ in June 2026.**

Below you can find a list of the models, some details on them, their usage and the references
where they were originally proposed, as well as references to the original Newtonian models upon which they were built.

<!-- TOC -->
* [non-newtonian-turbulence-models](#non-newtonian-turbulence-models)
  * [List of Models](#list-of-models)
    * [1. MalinKE](#1-malinke)
    * [2. BartosikKE](#2-bartosikke)
    * [3. GEkEpsilonZetaF](#3-gekepsilonzetaf)
    * [4. LKTSkOmegaSST](#4-lktskomegasst)
  * [Usage](#usage)
  * [Example cases](#example-cases)
  * [References](#references)
    * [Non-Newtonian Models](#non-newtonian-models)
    * [Baseline Newtonian Models](#baseline-newtonian-models)
<!-- TOC -->

## List of Models

### 1. MalinKE

Proposed in [**Malin (1998)**](), it is a low-Reynolds (Re) $k-\varepsilon$ model for Herschel-Bulkley fluids
based on the Newtonian model of [**Lam & Bremhorst (1981)**]().

### 2. BartosikKE

Proposed in [**Bartosik (2010)**](), it is a low-Re $k-\varepsilon$ model for Herschel-Bulkley fluids
based on the Newtonian model of [**Launder & Sharma (1974)**]().

### 3. GEkEpsilonZetaF

Proposed in [**Gavrilov & Rudyak (2016)**](), it is a four equation low-Re $k-\varepsilon-\zeta-f$ model
originally proposed for PL fluids, but which can also be extended to HB.
It is based on the Newtonian model of [**Hanjalic, Popovac & Hadziabdic (2004)**]().

### 4. LKTSkOmegaSST

Proposed in [**Lovato *et al.* (2022)**](), a $k-\omega$ SST model for Herschel-Bulkley fluids
based on the Newtonian model of [**Menter (1994)**]().

## Usage

## Example cases

## References

### Non-Newtonian Models

### Baseline Newtonian Models
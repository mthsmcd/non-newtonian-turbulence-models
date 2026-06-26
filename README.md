# non-newtonian-turbulence-models
OpenFOAM implementation of Reynolds-averaged (RANS) turbulence models for non-Newtonian fluids.
Developed and tested in OpenFOAM v2306.

Currently there are four models available, they were originally proposed for turbulent flows of
power-law (PL), Bingham (Bn) and Herschel-Bulkley (HB) fluids 
in the [reference papers](#references), and investigated in our review paper:

- **_"Turbulence modeling for viscoplastic fluids: state-of-the-art review and open-source investigation comparing with new DNS"_
 Matheus S. S. Macedo, Gustavo E. O. Celis, Hamidreza Anbarlooei, Roney L. Thompson,
Daiane M. Iceri, Letícia Bizarre, Vanessa Richon and Marcelo S. Castro,
accepted for publishing at the _Journal of Non-Newtonian Fluid Mechanics_ in June 2026.**

Below you can find a [list of the models](#turbulence-models) with brief descriptions, 
details on their [usage](#usage) and the [references](#non-newtonian-models)
where they were originally proposed. 
[References](#baseline-newtonian-models) to the baseline Newtonian models upon which they were built are also presented.
For full details please refer to the [review paper](#review-paper).

If you find this repository useful please cite our work, as well as the [original works](#references) of the implemented models.

## Repository structure

All turbulence and transport models are compiled together as part of the library
`libNonNewtonianRAS`.
The folders [`library/turbulence-models`](library/turbulence-models) and [`library/transport-models`](library/transport-models) contain the files
for all models.
To compile just run the `Allwmake` script.

The repository also contains some OpenFOAM example cases for the channel and pipe flows,
for more details please refer to the [usage section](#usage) and to our review paper [Macedo _et al._ (2026)](#review-paper).

## Library models

### Transport models

#### 1. RASHerschelBulkley

In the RANS models of this repository, the average molecular viscosity $\nu$ depends on turbulence, for this reason we have
created this separate transport model for turbulent HB fluids.
However, because each turbulence model has a particular way of modeling $\nu$,
the fluid parameters will be read by the **turbulence model** and $\nu$ will also be calculated
within the **turbulence model**.

The fluid parameters are alike the standard HB, with three additional ones `m`, `nu0` and `xi`. 
They should be defined in the `transportProperties` dictionary in the following manner:

```
RASHerschelBulkley
{
  // consistency index
  k     <value>;
  
  // flow behaviour index
  n     <value>; 
  
  // yield stress
  tau0  <value>; 
  
  // papanastasiou regularization parameter
  // only required in GRkEpsilonZetaF and LKTSkOmegaSST models
  m     <value>;
  
  // maximum viscosity of bi-viscosity regularization
  // only required in MalinKE
  nu0   <value>;
  
  // optional, yield stress ratio tau0/tauw
  // only used in BartosikKE
  xi    <value>;  // 0 < xi < 1.0
}
```

Note that to reduce the `RASHerschelBulkley` to Newtonian, PL and Bn fluids you just need to set
the parameters accordingly:

- Newtonian: $n = 1$ and $\tau_0 = 0$
- PL: $\tau_0 = 0$
- Bn: $n = 1$

All turbulence models require the definition of `k`, `n` and `tau0`. 
Other additional parameters are individually required 
by each RANS, please check them out below in the listing of the [turbulence models](#turbulence-models).

### Turbulence models

All models may be used in simulations with Newtonian, PL, Bn and HB fluids, but each has a fluid it is best suited for.

#### 1. `MalinKE`

Proposed in **Malin (1998)**, it is a low-Reynolds (Re) $k$-$\varepsilon$ model for Herschel-Bulkley fluids
based on the Newtonian model of **Lam & Bremhorst (1981)**.

Best suited for PL, transformed the damping function $f_\mu$ into a function of `n`.
Requires the definition of a maximum regularizing viscosity `nu0` in `transportProperties`.

#### 2. `BartosikKE`

Proposed in **Bartosik (2010)**, it is a low-Re $k$-$\varepsilon$ model for Herschel-Bulkley fluids
based on the Newtonian model of **Launder & Sharma (1974)**.

Best suited for Bn, transformed the damping function $f_\mu$ into a function of `xi`, the ratio of the yield stress to the wall shear stress
(_i.e._ $\xi = \tau_0 / \tau_w$).
`xi` may be calculated iteratively or a fixed value may be imposed at `transportProperties`.
It imposes a constant $\nu$ equal to the value at the wall.


#### 3. `GRkEpsilonZetaF`

Proposed in **Gavrilov & Rudyak (2016)**, it is a four equation low-Re $k$-$\varepsilon$-$\zeta$-$f$ model
originally proposed for PL fluids, but which can also be extended to HB.
It is based on the Newtonian model of **Hanjalić, Popovac & Hadžiabdić (2004)**.

Best suited for PL.
Requires the definition of a regularizing parameter `m` in `transportProperties`.

#### 4. `LKTSkOmegaSST`

Proposed in **Lovato *et al.* (2022)**, a $k$-$\omega$ SST model for Herschel-Bulkley fluids
based on the Newtonian model of **Menter (1994)**.

Best suited for Bn and HB.
Requires the definition of a regularizing parameter `m` in `transportProperties`.

## Usage

### Example cases

## References

### Review paper

1. Macedo, M. S. S., Celis, G. E. O., Anbarlooei, H. R., Thompson, R. L.,
   Iceri, D. M., Bizarre, L., Richon, V. and Castro M. S. 
"Turbulence modeling for viscoplastic fluids: state-of-the-art review and open-source investigation comparing with new DNS"
_Journal of Non-Newtonian Fluid Mechanics_ _**accepted for publishing**_ (2026).

### Non-Newtonian models

2. Malin, M. R.
  "Turbulent pipe flow of Herschel-Bulkley fluids"
  _International Communications in Heat and Mass Transfer_, 25(3), 321-330 (1998).
  DOI: [10.1016/S0735-1933(98)00019-0](https://doi.org/10.1016/S0735-1933(98)00019-0)

3. Bartosik, A.
   "Application of rheological models in prediction of turbulent slurry flow"
   _Flow, Turbulence and Combustion_, 84, 277-293 (2010).
   DOI: [10.1007/s10494-009-9234-y](https://doi.org/10.1007/s10494-009-9234-y)

4. Gavrilov, A. A. and Rudyak, V. Y.
   "Reynolds-averaged modeling of turbulent flows of power-law fluids"
   _Journal of Non-Newtonian Fluid Mechanics_, 227, 45-55 (2016).
   DOI: [10.1016/j.jnnfm.2015.11.006](https://doi.org/10.1016/j.jnnfm.2015.11.006)

5. Lovato, S., Keetels, G. H., Toxopeus, S. L. and Settels, J. W.
   "An eddy-viscosity model for turbulent flows of Herschel-Bulkley fluids"
   _Journal of Non-Newtonian Fluid Mechanics_, 301, 104729 (2022).
   DOI: [10.1016/j.jnnfm.2021.104729](https://doi.org/10.1016/j.jnnfm.2021.104729)

### Baseline Newtonian models

6. Lam, C. K. G. and Bremhorst, K. (1981)
   "A modified form of the $k$-$\varepsilon$ model for predicting wall turbulence"
   _Journal of Fluids Engineering_, 103(3), 456-460.
   DOI: [10.1115/1.3240815](https://doi.org/10.1115/1.3240815)

7. Launder, B. E. and Sharma, B. I. (1974)
   "Application of the energy-dissipation model of turbulence to the calculation of the flow near a spinning disc"
   _Letters in Heat and Mass Transfer_, 1(2), 131-137.
   DOI: [10.1016/0094-4548(74)90150-7](https://doi.org/10.1016/0094-4548(74)90150-7)

8. Hanjalić, K., Popovac, M. and Hadžiabdić, M. (2004)
   "A robust near-wall elliptic-relaxation eddy-viscosity turbulence model for CFD"
   _International Journal of Heat and Fluid Flow_, 25, 1047-1051.
   DOI: [10.1016/j.ijheatfluidflow.2004.07.005](https://doi.org/10.1016/j.ijheatfluidflow.2004.07.005)

9. Menter, F. (1994)
   "Two-equation eddy-viscosity turbulence models for engineering applications"
   _AIAA Journal_, 32(8), 1598-1605.
   DOI: [10.2514/3.12149](https://doi.org/10.2514/3.12149)
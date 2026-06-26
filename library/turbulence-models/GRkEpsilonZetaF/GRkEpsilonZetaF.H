/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.

    2026 M. S. S. Macedo, Coppe/UFRJ
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::RASModels::GRkEpsilonZetaF

Group
    grpRASTurbulence

Description
    Gavrilov & Rudyak k-epsilon-zeta-f turbulence model for power-law fluids.

    This turbulence model is described in:
    \verbatim
        Gavrilov, A. A., Rudyak, V. Y. (2016).
        Reynolds-averaged modeling of turbulent flows of power-law fluids.
        Journal of non-Newtonian fluid mechanics, 227, 45-55.
    \endverbatim

    The model is built upon the Newtonian Hanjalic', Popovac & Hadziabdic'
    k-epsilon-zeta-f model:

    \verbatim
        Hanjalic', K., Popovac, M., M. Hadziabdic' (2004).
        A robust near-wall elliptic-relaxation eddy-viscosity turbulence model
        for CFD.
        International journal of heat and fluid flow, 25, 1047-1051.
    \endverbatim

    It was investigated with other viscoplastic models and validated against
    DNS in:

    \verbatim
        Macedo, M. S. S., Celis, G. E. O., Anbarlooei, H. R., Thompson, R. L.,
        Iceri, D. M., Bizarre, L., Richon, V., Castro, M. S. (2026).
        Turbulence modeling for viscoplastic fluids: state-of-the-art review
        and open-source investigation comparing with new DNS.
        Journal of non-Newtonian fluid mechanics
    \endverbatim

SourceFiles
    GRkEpsilonZetaF.C

SeeAlso
    kEpsilon.C

\*---------------------------------------------------------------------------*/

#ifndef GRkEpsilonZetaF_H
#define GRkEpsilonZetaF_H

#include "turbulentTransportModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                        Class GRkEpsilonZetaF Declaration
\*---------------------------------------------------------------------------*/

class GRkEpsilonZetaF
:
    public eddyViscosity<incompressible::RASModel>
{
    // Private Member Functions

        // Viscoplastic model
        void correctApparentViscosity(const volScalarField& etta);
        dimensionedScalar wallFriction() const;
        tmp<volScalarField> dNudG(const volScalarField& sr) const;
        tmp<volScalarField> strainRate() const;
        tmp<volScalarField> calcEtta(const volScalarField& sr) const;
        tmp<volScalarField> calcNuNN(const volScalarField& sr) const;

protected:

    // Protected Data

        // Model coefficients

            dimensionedScalar Cmu_;
            dimensionedScalar Ceps1a_;
            dimensionedScalar Ceps1b_;
            dimensionedScalar Ceps2_;
            dimensionedScalar Cf1_;
            dimensionedScalar Cf2_;
            dimensionedScalar CL_;
            dimensionedScalar Ceta_;
            dimensionedScalar Ctau_;
            dimensionedScalar sigmaK_;
            dimensionedScalar sigmaEps_;
            dimensionedScalar sigmaZeta_;

        // non-Newtonian coefficients
            dimensionedScalar Cn_;
            dimensionedScalar Cepsn_;

        // Fields

            //- Turbulent kinetic energy [m2/s2]
            volScalarField k_;

            //- Turbulent kinetic energy dissipation rate [m2/s3]
            volScalarField epsilon_;

            //- Normalised wall-normal fluctuating velocity scale [-]
            volScalarField zeta_;

            //- Elliptic relaxation factor [1/s]
            volScalarField f_;

            //- Turbulent time scale [s]
            volScalarField T_;

            //- non-Newtonian viscosity
            volScalarField nuNN_;

        // Bounding values

            dimensionedScalar zetaMin_;
            dimensionedScalar fMin_;
            dimensionedScalar TMin_;
            dimensionedScalar L2Min_;
            dimensionedScalar ettaMin_;

            // flag to control if epsilon is solved with zeroed BC
            // standard is true.
            Switch epsilonZeroAtWall_;

            // Herschel Bulkley variables
            IOdictionary transportProperties_;
            dictionary HerschelBulkleyCoeffs_;
            dimensionedScalar Kind_;
            dimensionedScalar n_;
            dimensionedScalar tau0_;
            dimensionedScalar m_;
            Switch NNFlag_;
            Switch printXi_;

            // number of internal nu calculations per iteration
            label nCorrNu_;

            //- Apparent viscosity
            volScalarField  etta_;

            //- Near wall epsilon
            volScalarField epsw_;

    // Protected Member Functions

        //- Update nut with the latest available k, zeta, and T
        virtual void correctNut();

        //- Return the turbulent time scale, T
        tmp<volScalarField> Ts() const;

        //- Return the turbulent length scale, L
        tmp<volScalarField> Ls() const;

        virtual void correctApparentViscosity();

        //- Update non-Newtonian viscosity
        virtual void correctNuNN();

public:

    //- Runtime type information
    TypeName("GRkEpsilonZetaF");


    // Constructors

        //- Construct from components
        GRkEpsilonZetaF
        (
            const geometricOneField& alpha,
            const geometricOneField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~GRkEpsilonZetaF() = default;


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the turbulent kinetic energy field
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulent kinetic energy dissipation rate field
        virtual tmp<volScalarField> epsilon() const
        {
            if (epsilonZeroAtWall_)
            {
                return epsilon_ + epsw_;
            }
            return epsilon_;
        }

        //- Return the normalised wall-normal fluctuating velocity scale field
        virtual tmp<volScalarField> zeta() const
        {
            return zeta_;
        }

        //- Return the elliptic relaxation factor field
        virtual tmp<volScalarField> f() const
        {
            return f_;
        }

        virtual tmp<volScalarField> nu() const
        {
            return etta_;
        }

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DkEff",
                    nut_/sigmaK_ + etta_
                )
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DepsilonEff",
                    nut_/sigmaEps_ + etta_
                )
            );
        }

        //- Return the effective diffusivity for zeta (LUU:Eq. 17)
        tmp<volScalarField> DZetaEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DzetaEff",
                    nut_/sigmaZeta_ + etta_
                )
            );
        }

        tmp<volScalarField> nuEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "nuEff",
                    nut_ + nuNN_ + etta_
                )
            );
        }

        //- Return the source term into the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff(volVectorField& U) const;

        //- Solve the transport equations and correct the turbulent viscosity
        virtual void correct();

    };


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

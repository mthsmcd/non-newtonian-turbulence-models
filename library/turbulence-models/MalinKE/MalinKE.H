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
    Foam::incompressible::RASModels::MalinKE

Group
    grpIcoRASTurbulence

Description
    Malin low-Reynolds number k-epsilon turbulence model for
    viscoplastic fluids of the Herschel-Bulkley model.

    This turbulence model is described in:
    \verbatim
        Malin, M. R. (1998).
        Turbulent pipe flow of Herschel-Bulkley fluids.
        International communications in heat and mass transfer, 25(3), 321-330.
    \endverbatim

    The model is built upon the Newtonian Lam & Bremhorst k-epsilon model:

    \verbatim
        Lam, C. K. G., & Bremhorst, K. (1981).
        A modified form of the k-ε model for predicting wall turbulence.
        Journal of Fluids Engineering, 103(3), 456-460.
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
    MalinKE.C

\*---------------------------------------------------------------------------*/

#ifndef MalinKE_H
#define MalinKE_H

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
                      Class MalinKE Declaration
\*---------------------------------------------------------------------------*/

class MalinKE
:
    public eddyViscosity<incompressible::RASModel>
{
    // Private Member Functions

        tmp<volScalarField> Rt() const;
        tmp<volScalarField> fMu(const volScalarField& Rt) const;
        tmp<volScalarField> f1(const volScalarField& fMu) const;
        tmp<volScalarField> f2(const volScalarField& Rt) const;

        void correctNut(const volScalarField& fMu);
        void correctApparentViscosity(const volScalarField& etta);
        dimensionedScalar wallFriction() const;

        // Viscoplastic model
        tmp<volScalarField> strainRate() const;
        tmp<volScalarField> calcEtta(const volScalarField& sr) const;


protected:

    // Protected data

        // model constants
        dimensionedScalar Cmu_;
        dimensionedScalar Ceps1_;
        dimensionedScalar Ceps2_;
        dimensionedScalar Ceps3_;
        dimensionedScalar sigmaEps_;

        // unyielded regions on/off flag
        Switch unyielded_;

        // turbulent quantities
        volScalarField k_;
        volScalarField epsilon_;

        //- Wall distance
        //  Note: different to wall distance in parent RASModel
        //  which is for near-wall cells only
        const volScalarField& y_;

        // Herschel Bulkley variables
        IOdictionary transportProperties_;
        dictionary HerschelBulkleyCoeffs_;
        dimensionedScalar Kind_;
        dimensionedScalar n_;
        dimensionedScalar tau0_;
        dimensionedScalar nu0_;
        volScalarField  etta_;

    // Protected Member Functions

        virtual void correctNut();
        virtual void correctApparentViscosity();

public:

    //- Runtime type information
    TypeName("MalinKE");


    // Constructors

        //- Construct from components
        MalinKE
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
    virtual ~MalinKE() = default;


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        // return effective viscosity
        tmp<volScalarField> nuEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "nuEff",
                    nut_ + nu()
                )
            );
        }

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", nut_ + nu())
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DepsilonEff", nut_/sigmaEps_ + nu())
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Return the source term into the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff(volVectorField& U) const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        virtual tmp<volScalarField> nu() const
        {
            return etta_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

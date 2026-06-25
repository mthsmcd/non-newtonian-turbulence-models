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
    Foam::incompressible::RASModels::BartosikKE

Group
    grpIcoRASTurbulence

Description
    Bartosik low-Reynolds number k-epsilon turbulence model for
    viscoplastic fluids of the Herschel-Bulkley model.

    This turbulence model is described in:
    \verbatim
        Bartosik, A. (2010).
        Application of Rheological Models in Prediction of Turbulent Slurry Flow.
        Flow, turbulence and combustion, 84, 277-293.
    \endverbatim

    The model is built upon the Newtonian Launder & Sharma k-epsilon model:

    \verbatim
        Launder, B. E., & Sharma, B. I. (1974).
        Application of the energy-dissipation model of turbulence to the
        calculation of flow near a spinning disc.
        Letters in heat and mass transfer, 1(2), 131-137.
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
    BartosikKE.C

\*---------------------------------------------------------------------------*/

#ifndef BartosikKE_H
#define BartosikKE_H

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
                      Class BartosikKE Declaration
\*---------------------------------------------------------------------------*/

class BartosikKE
:
    public eddyViscosity<incompressible::RASModel>
{
    // Private Member Functions

        tmp<volScalarField> fMu(const dimensionedScalar& tauw) const;
        tmp<volScalarField> f2() const;

        void correctNut(const volScalarField& fMu);
        void correctApparentViscosity(const dimensionedScalar& etta);

        // Viscoplastic model
        dimensionedScalar wallFriction() const;
        dimensionedScalar calcEtta(const dimensionedScalar& tauw) const;

protected:

    // Protected data

        // model constants
        dimensionedScalar Cmu_;
        dimensionedScalar C1_;
        dimensionedScalar C2_;
        dimensionedScalar C3_;
        dimensionedScalar sigmak_;
        dimensionedScalar sigmaEps_;

        // turbulent quantities
        volScalarField k_;
        volScalarField epsilon_;

        // Herschel Bulkley variables
        IOdictionary transportProperties_;
        dictionary HerschelBulkleyCoeffs_;
        dimensionedScalar Kind_;
        dimensionedScalar n_;
        dimensionedScalar tau0_;
        dimensionedScalar xi_;

        // Apparent viscosity
        volScalarField etta_;

    // Protected Member Functions

        virtual void correctNut();
        virtual void correctApparentViscosity();


public:

    //- Runtime type information
    TypeName("BartosikKE");


    // Constructors

        //- Construct from components
        BartosikKE
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
    virtual ~BartosikKE() = default;


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
                    nut_ + etta_
                )
            );
        }

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", nut_ + etta_)
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DepsilonEff", nut_/sigmaEps_ + etta_)
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

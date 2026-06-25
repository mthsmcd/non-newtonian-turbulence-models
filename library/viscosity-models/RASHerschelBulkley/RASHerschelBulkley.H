/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2019 OpenCFD Ltd.
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
    Foam::viscosityModels::RASHerschelBulkley

Description
     Herschel-Bulkley non-Newtonian viscosity model.

SourceFiles
    RASHerschelBulkley.C

\*---------------------------------------------------------------------------*/

#ifndef HerschelBulkley_H
#define HerschelBulkley_H

#include "viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "turbulentTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class RASHerschelBulkley Declaration
\*---------------------------------------------------------------------------*/

class RASHerschelBulkley
:
    public viscosityModel
{
    // Private data

        dictionary RASHerschelBulkleyCoeffs_;

        dimensionedScalar k_;
        dimensionedScalar n_;
        dimensionedScalar tau0_;
        dimensionedScalar nu0_;


protected:

    // Protected data

//        volScalarField nu_;


    // Protected Member Functions

        //- Calculate and return the laminar viscosity
//        tmp<volScalarField> calcNu() const;


public:

    //- Runtime type information
    TypeName("RASHerschelBulkley");


    // Constructors

        //- Construct from components
        RASHerschelBulkley
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~RASHerschelBulkley() = default;


    // Member Functions

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const
        {
            // Incompressible
            {
                const auto* turbPtr =
                        U_.db().findObject<incompressible::turbulenceModel>
                                (
                                        turbulenceModel::propertiesName
                                );

                if (turbPtr)
                {
                    return tmp<volScalarField>(turbPtr->nu());
                }
            }

            FatalErrorInFunction
                    << " The turbulence model is not found in the database."
                    << exit(FatalError);

            return nullptr;
        }
//
        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const
        {
//            // Incompressible
//            {
//                const auto* turbPtr =
//                        U_.db().findObject<incompressible::turbulenceModel>
//                                (
//                                        turbulenceModel::propertiesName
//                                );
//
//                if (turbPtr)
//                {
//                    return (turbPtr->nu()->boundaryField()[patchi]);
//                }
//            }
//
//            FatalErrorInFunction
//                    << " The turbulence model is not found in the database."
//                    << exit(FatalError);
//
//            return nullptr;
//            return nu_.boundaryField()[patchi];
            volScalarField nu_(nu());
            return nu_.boundaryField()[patchi];
        }
//
        //- Correct the laminar viscosity
        virtual void correct()
        {

        }

        // Return power law index
        virtual dimensionedScalar n() const
        {
            return n_;
        }

        // Return consistency index
        virtual dimensionedScalar K() const
        {
            return k_;
        }

        // Return yield stress
        virtual dimensionedScalar tau0() const
        {
            return tau0_;
        }

        // Return unyielded viscosity
        virtual dimensionedScalar nu0() const
        {
            return nu0_;
        }

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& viscosityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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

\*---------------------------------------------------------------------------*/

#include "RASHerschelBulkley.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(RASHerschelBulkley, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        RASHerschelBulkley,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::RASHerschelBulkley::RASHerschelBulkley
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    RASHerschelBulkleyCoeffs_
    (
        viscosityProperties.optionalSubDict(typeName + "Coeffs")
    ),
        k_("k", dimViscosity, RASHerschelBulkleyCoeffs_),
        n_("n", dimless, RASHerschelBulkleyCoeffs_),
        tau0_("tau0", dimViscosity/dimTime, RASHerschelBulkleyCoeffs_),
        nu0_("nu0", dimViscosity, RASHerschelBulkleyCoeffs_)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::RASHerschelBulkley::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    RASHerschelBulkleyCoeffs_ =
        viscosityProperties.optionalSubDict(typeName + "Coeffs");

    RASHerschelBulkleyCoeffs_.readEntry("k", k_);
    RASHerschelBulkleyCoeffs_.readEntry("n", n_);
    RASHerschelBulkleyCoeffs_.readEntry("tau0", tau0_);
    RASHerschelBulkleyCoeffs_.readEntry("nu0", nu0_);

    return true;
}


// ************************************************************************* //

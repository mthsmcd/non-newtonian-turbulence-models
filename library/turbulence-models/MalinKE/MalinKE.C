/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2020,2023 OpenCFD Ltd.

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

#include "MalinKE.H"
#include "wallDist.H"
#include "bound.H"
#include "addToRunTimeSelectionTable.H"
#include "fvOptions.H"
#include "wallPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(MalinKE, 0);
addToRunTimeSelectionTable(RASModel, MalinKE, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> MalinKE::Rt() const
{
    return sqr(k_)/(nu()*epsilon_ + dimensionedScalar ("SMALL", pow4(dimVelocity), SMALL));
}

tmp<volScalarField> MalinKE::fMu(const volScalarField& Rt) const
{
    tmp<volScalarField> Ry(sqrt(k_)*y_/(nu() + dimensionedScalar ("SMALL", dimViscosity, SMALL)));
    tmp<volScalarField> fmu = sqr(scalar(1) - exp(-0.0165*Ry/pow(n_, 0.25)))*(scalar(1) + 20.5/(Rt + SMALL));

    // boundings imposed in
    // R. Willems - Numerical Analysis of Turbulent Heat Transfer within an axial-flux permanent magnet machine
    // - Master Thesis - T.U. Eindhoven (2020)
    return max(fmu, scalar(0.01116225));
}

tmp<volScalarField> MalinKE::f1(const volScalarField& fMu) const
{
    return scalar(1) + pow3(0.05/(fMu + SMALL));
}

tmp<volScalarField> MalinKE::f2(const volScalarField& Rt) const
{
    return scalar(1) - exp(-sqr(Rt));
}

void MalinKE::correctNut(const volScalarField& fMu)
{
    nut_ = Cmu_*fMu*sqr(k_)/epsilon_;
    nut_.correctBoundaryConditions();
}

void MalinKE::correctApparentViscosity(const volScalarField& etta)
{
    etta_ = etta;
    etta_.correctBoundaryConditions();
}

tmp<volScalarField> MalinKE::strainRate() const
{
    return sqrt(2.0)*mag(symm(fvc::grad(U_)));
}

tmp<volScalarField> MalinKE::calcEtta(const volScalarField& sr) const
{

    dimensionedScalar tone("tone", dimless/sr.dimensions(), 1.0);
    dimensionedScalar rtone("rtone", dimless/dimTime, 1.0);

    if (!unyielded_)
    {
        return (tau0_ + Kind_*rtone*pow(tone*sr, n_))
               /(max(sr, dimensionedScalar ("SMALL", sr.dimensions(), SMALL)));
    }

    return min
    (
        nu0_,
        (tau0_ + Kind_*rtone*pow(tone*sr, n_))
        /(max(sr, dimensionedScalar ("SMALL", sr.dimensions(), SMALL)))
    );

}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void MalinKE::correctNut()
{
    correctNut(fMu(Rt()));
}

void MalinKE::correctApparentViscosity()
{
    correctApparentViscosity(calcEtta(strainRate()));

    if (tau0_.value() > 0.0)
    {
        dimensionedScalar tw = wallFriction();
        Info << "Yield stress ratio tau0/tauw = " << tau0_.value()/tw.value() << endl;
    }

    Info << "Average apparent viscosity etta = " << gAverage(etta_) << endl;
}

dimensionedScalar MalinKE::wallFriction() const
{

    volScalarField shearMag("shearMag",
                            sqrt(2.0)*mag(symm(fvc::grad(U_)))
                            );
    shearMag.correctBoundaryConditions();

    // shear magnitude boundary field
    volScalarField::Boundary& sMBf = shearMag.boundaryFieldRef();

    // etta boundary field
    const volScalarField::Boundary& nuBf = etta_.boundaryField();

    // get mesh boundaries
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    scalar totalArea = 0.0;
    scalar tauSum = 0.0;
    forAll (pbm, patchi)
    {
        if (isA<wallPolyPatch>(pbm[patchi]))
        {
            // shear rate magnitude at wall
            fvPatchField<double>& gDotw = sMBf[patchi];
            const fvPatchField<double>& nuw = nuBf[patchi];

            // faces areas at patchi
            fvsPatchField<double> area = mesh_.magSf().boundaryField()[patchi];

            // integrate the friction factor over all wall faces
            forAll(pbm[patchi], facei)
            {
                scalar a = area[facei];
                totalArea += a;
                tauSum += a*gDotw[facei]*nuw[facei];
            } //faces
        } // if wall patch
    } // patches

    dimensionedScalar tauw ("tauw", tau0_.dimensions(), tauSum/totalArea);

    // writing tauw
    if (mesh_.time().writeTime())
    {
        IOdictionary propsDict
        (
            IOobject
            (
                typeName + "Properties",
                mesh_.time().timeName(),
                "uniform",
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            )
        );
        propsDict.add("tauw", tauw);
        propsDict.regIOobject::write();
    }

    return tauw;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

MalinKE::MalinKE
(
    const geometricOneField& alpha,
    const geometricOneField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<incompressible::RASModel>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    Cmu_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cmu",
            coeffDict_,
            0.09
        )
    ),

    Ceps1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Ceps1",
            coeffDict_,
            1.44
        )
    ),

    Ceps2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Ceps2",
            coeffDict_,
            1.92
        )
    ),

    Ceps3_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Ceps3",
            coeffDict_,
            0
        )
    ),

    sigmaEps_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaEps",
            coeffDict_,
            1.3
        )
    ),

    unyielded_(coeffDict_.getOrDefault<Switch>("unyielded", true)),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    epsilon_
    (
        IOobject
        (
            IOobject::groupName("epsilon", alphaRhoPhi.group()),
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    y_(wallDist::New(mesh_).y()),

    transportProperties_
    (
        IOobject
        (
            "transportProperties",
            runTime_.constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    HerschelBulkleyCoeffs_
    (
        transportProperties_.subDict("RASHerschelBulkleyCoeffs")
    ),

    // consistency index
    Kind_("k", dimViscosity, HerschelBulkleyCoeffs_),

    // power law index
    n_("n", dimless, HerschelBulkleyCoeffs_),

    // yield stress
    tau0_("tau0", dimViscosity/dimTime, HerschelBulkleyCoeffs_),

    // lowest apparent viscosity value
    nu0_("nu0", dimViscosity, HerschelBulkleyCoeffs_),

    etta_
    (
        IOobject
        (
            IOobject::groupName("nu", nut_.group()),
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        Kind_
    )

{
    bound(k_, kMin_);
    bound(epsilon_, epsilonMin_);

    if (type == typeName)
    {
        printCoeffs(type);
        Info << "Herschel-Bulkley parameters\n" << "K = " << Kind_ << endl;
        Info << "n = " << n_ << endl;
        Info << "tau0 = " << tau0_ << endl;
        Info << "nu0 = " << nu0_ << endl;

        if (unyielded_)
        {
            Info << "Unyielded flow active." << endl;
        }
        else
        {
            Info << "No unyielded flow." << endl;
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool MalinKE::read()
{
    if (eddyViscosity<incompressible::RASModel>::read())
    {
        Cmu_.readIfPresent(coeffDict());
        Ceps1_.readIfPresent(coeffDict());
        Ceps2_.readIfPresent(coeffDict());
        Ceps3_.readIfPresent(coeffDict());
        sigmaEps_.readIfPresent(coeffDict());

        unyielded_.readIfPresent("unyielded", coeffDict());

        HerschelBulkleyCoeffs_ = coeffDict().topDict().subDict("HerschelBulkleyCoeffs");
        HerschelBulkleyCoeffs_.readEntry("k", Kind_);
        HerschelBulkleyCoeffs_.readEntry("n", n_);
        HerschelBulkleyCoeffs_.readEntry("tau0", tau0_);
        HerschelBulkleyCoeffs_.readEntry("nu0", nu0_);

        return true;
    }

    return false;
}


void MalinKE::correct()
{
    if (!turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<incompressible::RASModel>::correct();

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField G(this->GName(), nut*(tgradU() && devTwoSymm(tgradU())));
    tgradU.clear();

    // Update epsilon and G at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();
    // Push any changed cell values to coupled neighbours
    epsilon_.boundaryFieldRef().evaluateCoupled<coupledFvPatch>();

    const volScalarField Rt(this->Rt());
    const volScalarField fMu(this->fMu(Rt));

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon_)
        + fvm::div(alphaRhoPhi, epsilon_)
        - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
        ==
        Ceps1_*f1(fMu)*alpha*rho*G*epsilon_/k_
        - fvm::Sp(Ceps2_*f2(Rt)*alpha*rho*epsilon_/k_, epsilon_)
        + fvOptions(alpha, rho, epsilon_)
    );

    epsEqn.ref().relax();
    fvOptions.constrain(epsEqn.ref());
    epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
    solve(epsEqn);
    fvOptions.correct(epsilon_);
    bound(epsilon_, this->epsilonMin_);

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
        + fvm::div(alphaRhoPhi, k_)
        - fvm::laplacian(alpha*rho*DkEff(), k_)
        ==
        alpha*rho*G
        - fvm::Sp(alpha*rho*epsilon_/max(k_, dimensionedScalar ("SMALL", k_.dimensions(),SMALL)), k_)
        + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

// Update nut with latest available k,epsilon
    correctNut();

// update etta_
    correctApparentViscosity();
}

Foam::tmp<Foam::fvVectorMatrix>
MalinKE::divDevRhoReff
(
    volVectorField& U
) const
{

    return
    (
        - fvc::div(this->alpha_*this->rho_*this->nuEff()*dev2(T(fvc::grad(U))))
        - fvm::laplacian(this->alpha_*this->rho_*this->nuEff(), U)
    );

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //

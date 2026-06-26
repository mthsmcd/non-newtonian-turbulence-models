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

#include "BartosikKE.H"
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

defineTypeNameAndDebug(BartosikKE, 0);
addToRunTimeSelectionTable(RASModel, BartosikKE, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> BartosikKE::fMu(const dimensionedScalar& tauw) const
{
    dimensionedScalar xi(xi_);
    if (xi.value() == 0.0)
    {
        xi = tau0_/tauw;
    }
    return exp(-3.4*(1 + xi)/sqr(scalar(1) + sqr(k_)/(etta_*epsilon_)/50.0));;
}

tmp<volScalarField> BartosikKE::f2() const
{
    return
            scalar(1)
          - 0.3*exp(-min(sqr(sqr(k_)/(etta_*epsilon_)), scalar(50)));
}

void BartosikKE::correctNut(const volScalarField& fMu)
{
    nut_ = Cmu_*fMu*sqr(k_)/epsilon_;
    nut_.correctBoundaryConditions();
}

// Viscoplastic model functions
void BartosikKE::correctApparentViscosity(const dimensionedScalar& etta)
{
    etta_ = etta;
    etta_.correctBoundaryConditions();
}

dimensionedScalar BartosikKE::wallFriction() const
{
    volScalarField shearMag("shearMag",
                            sqrt(2.0)*mag(symm(fvc::grad(U_)))
                            );
    shearMag.correctBoundaryConditions();

    // shear magnitude boundary field
    volScalarField::Boundary& sMBf = shearMag.boundaryFieldRef();

    // get mesh boundaries
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    scalar totalArea = 0.0;
    scalar tauSum = 0.0;
    dimensionedScalar ettaAvg = gAverage(etta_);
    forAll (pbm, patchi)
    {
        if (isA<wallPolyPatch>(pbm[patchi]))
        {
            // shear rate magnitude at wall
            fvPatchField<double>& gDotw = sMBf[patchi];

            // faces areas at patchi
            fvsPatchField<double> area = mesh_.magSf().boundaryField()[patchi];

            // integrate the friction factor over all wall faces
            forAll(pbm[patchi], facei)
            {
                scalar a = area[facei];
                totalArea += a;
                tauSum += a*gDotw[facei]*ettaAvg.value();
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


dimensionedScalar BartosikKE::calcEtta(const dimensionedScalar& tauw) const
{
    dimensionedScalar tone("tone", dimViscosity/tauw.dimensions(), 1.0);
    dimensionedScalar rtone("rtone", tauw.dimensions()/Kind_.dimensions(), 1.0);

    dimensionedScalar tw(tauw);
    if (xi_.value() > 0.0)
    {
        tw = tau0_/xi_;
    }
    return tw * tone * pow(rtone*Kind_/(tw - tau0_), 1.0/n_);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void BartosikKE::correctNut()
{
    correctNut(fMu(wallFriction()));
}

void BartosikKE::correctApparentViscosity()
{

    dimensionedScalar tw = wallFriction();
    correctApparentViscosity(calcEtta(tw));

    if (xi_.value() != 0.0)
    {
        Info << "Imposed yield stress ratio xi = " << xi_.value() << endl;
    }
    Info << "Yield stress ratio tau0/tauw = " << tau0_.value()/tw.value() << endl;

    Info << "Apparent viscosity etta = " << gAverage(etta_) << endl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

BartosikKE::BartosikKE
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

    C1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C1",
            coeffDict_,
            1.44
        )
    ),

    C2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C2",
            coeffDict_,
            1.92
        )
    ),

    C3_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C3",
            coeffDict_,
            0
        )
    ),

    sigmak_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "sigmak",
            coeffDict_,
            1.0
        )
    ),

    sigmaEps_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "sigmaEps",
            coeffDict_,
            1.3
        )
    ),

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

    xi_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "xi",
            HerschelBulkleyCoeffs_,
            dimless,
            0.0
        )
    ),

    etta_
    (
        IOobject
        (
            "nu",
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

        if (xi_.value() == 0.0)
        {
            Info << "Yield stress ratio xi not specified, will be calculated iteratively!" << endl;
        }
        else
        {
            Info << "Specified xi = " << xi_.value() << endl;
        }
        if (n_.value() == 1.0 && tau0_.value() == 0.0)
        {
            Info << "\nNewtonian fluid being used in the simulation!" << endl;
        }
    }

    if
    (
        Kind_.value() < 0.0
        || tau0_.value() < 0.0
        || n_.value() < 0.0
        || xi_.value() < 0.0
    )
    {
        FatalErrorInFunction
            << "Negative value used for K, tau0, n or xi:" << nl
            << "K = " << Kind_.value() << nl
            << "tau0 = " << tau0_.value() << nl
            << "n = " << n_.value() << nl
            << "xi = " << xi_.value() << nl
            << exit(FatalError);
    }

    if
    (
        mag(sigmaEps_.value()) < VSMALL
    )
    {
        FatalErrorInFunction
            << "Non-zero values are required for the model constants:" << nl
            << "sigmaEps = " << sigmaEps_ << nl
            << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool BartosikKE::read()
{
    if (eddyViscosity<incompressible::RASModel>::read())
    {
        Cmu_.readIfPresent(coeffDict());
        C1_.readIfPresent(coeffDict());
        C2_.readIfPresent(coeffDict());
        C3_.readIfPresent(coeffDict());
        sigmak_.readIfPresent(coeffDict());
        sigmaEps_.readIfPresent(coeffDict());

        HerschelBulkleyCoeffs_ = transportProperties_.subDict("RASHerschelBulkleyCoeffs");
        HerschelBulkleyCoeffs_.readEntry("k", Kind_);
        HerschelBulkleyCoeffs_.readEntry("n", n_);
        HerschelBulkleyCoeffs_.readEntry("tau0", tau0_);
        HerschelBulkleyCoeffs_.readEntry("xi", xi_);

        return true;
    }

    return false;
}


void BartosikKE::correct()
{
    if (!this->turbulence_)
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

    // Calculate parameters and coefficients for Launder-Sharma low-Reynolds
    // number model
    volScalarField E(2.0*this->etta_*nut*fvc::magSqrGradGrad(U));
    volScalarField D(2.0*this->etta_*magSqr(fvc::grad(sqrt(k_))));

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField G(this->GName(), nut*(tgradU() && devTwoSymm(tgradU())));
    tgradU.clear();

    // Update epsilon and G at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();
    // Push any changed cell values to coupled neighbours
    epsilon_.boundaryFieldRef().evaluateCoupled<coupledFvPatch>();

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon_)
      + fvm::div(alphaRhoPhi, epsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
     ==
        C1_*alpha*rho*G*epsilon_/k_
      - fvm::Sp(C2_*f2()*alpha*rho*epsilon_/k_, epsilon_)
      + alpha*rho*E
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
      - fvm::Sp(alpha*rho*(epsilon_ + D)/k_, k_)
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    correctNut();
    correctApparentViscosity();

}

Foam::tmp<Foam::fvVectorMatrix>
BartosikKE::divDevRhoReff
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

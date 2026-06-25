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

#include "LKTSkOmegaSST.H"
#include "fvOptions.H"
#include "bound.H"
#include "addToRunTimeSelectionTable.H"
#include "wallPolyPatch.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(LKTSkOmegaSST, 0);
addToRunTimeSelectionTable(RASModel, LKTSkOmegaSST, dictionary);

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> LKTSkOmegaSST::F1
(
    const volScalarField& CDkOmega
) const
{
    tmp<volScalarField> CDkOmegaPlus = max
    (
        CDkOmega,
        dimensionedScalar(dimless/sqr(dimTime), 1.0e-21)
    );

        tmp<volScalarField> arg1 = min
        (
            max
            (
                (scalar(1)/betaStar_)*sqrt(k_)/(omega_*y_),
                scalar(500)*etta_/(sqr(y_)*omega_)
            ),
            (4*alphaOmega2_)*k_/(CDkOmegaPlus*sqr(y_))
        );

    return tanh(pow4(arg1));
}

tmp<volScalarField> LKTSkOmegaSST::F2() const
{
    tmp<volScalarField> arg2 = max
    (
        (scalar(2)/betaStar_)*sqrt(k_)/(omega_*y_),
        scalar(500)*etta_/(sqr(y_)*omega_)
    );

    return tanh(sqr(arg2));
}


void LKTSkOmegaSST::correctNut()
{
    volScalarField magS(sqrt(2.0*magSqr(symm(fvc::grad(U_)))));
    nut_ = a1_*k_/max(a1_*omega_, b1_*F2()*magS);
    nut_.correctBoundaryConditions();
}

tmp<volScalarField> LKTSkOmegaSST::Pk
(
    const volScalarField& G
) const
{
    // the  2.0 factor comes from the TMR website
    return min(G, 2.0*(c1_*betaStar_)*k_*omega_);
}

tmp<volScalarField> LKTSkOmegaSST::epsilonByk
(

) const
{
    return betaStar_*omega_;
}

tmp<volScalarField> LKTSkOmegaSST::GbyNu0
(
    const volTensorField& gradU
) const
{
    return tmp<volScalarField>::New
    (
        IOobject::scopedName(this->type(), "GbyNu"),
        gradU && devTwoSymm(gradU)
    );
}

tmp<volScalarField> LKTSkOmegaSST::GbyNu
(
    const volScalarField& S2,
    const volScalarField& F2
) const
{
    return min
    (
        S2,
        (c1_/a1_)*betaStar_*omega_*max(a1_*omega_, b1_*F2*sqrt(S2))
    );
}

void LKTSkOmegaSST::correctApparentViscosity()
{
    Info << "Updating apparent viscosity etta" << endl;
    volScalarField etta0 = etta_;
    for (int i = 0; i < nCorrNu_; ++i)
    {
        correctApparentViscosity(calcEtta(strainRate()));
    }

    Info << "Etta residual after " << nCorrNu_ << " iterations = " << sum(mag(etta_ - etta0)).value() << endl;

    if (tau0_.value() > 0.0)
    {
        dimensionedScalar tw = wallFriction();
        Info << "Yield stress ratio tau0/tauw = " << tau0_.value()/tw.value() << endl;
    }
    Info << "Average apparent viscosity etta = " << gAverage(etta_) << endl;
}

void LKTSkOmegaSST::correctNuNN()
{
    nuNN_ = calcNuNN(strainRate());
    nuNN_.correctBoundaryConditions();
}


dimensionedScalar LKTSkOmegaSST::wallFriction() const
{

    volScalarField shearMag("shearMag",
                            sqrt(2.0)*mag(symm(fvc::grad(U_)))
                            );

    shearMag.correctBoundaryConditions();

    // shear magnitude boundary field
    volScalarField::Boundary& sMBf = shearMag.boundaryFieldRef();

    volScalarField nueff("nuEff", nuEff());
    const volScalarField::Boundary& nuBf = nueff.boundaryField();

    // get mesh boundaries
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    scalar totalArea = 0.0;
    scalar tauSum = 0.0;
    forAll (pbm, patchi)
    {
        if (isA<wallPolyPatch>(pbm[patchi]))
        {
            // shear rate magnitude at wall = gamma_dot_w
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

tmp<volScalarField> LKTSkOmegaSST::dNudG(const volScalarField& sr) const
{
    dimensionedScalar tone("tone", dimless/sr.dimensions(), 1.0);
    dimensionedScalar rtone("rtone", dimless/dimTime, 1.0);
    return Cflag_*(-tau0_*(1.0 - exp(-m_*tone*sr)) + m_*sr*tone*tau0_*exp(-m_*tone*sr) + (n_ - 1.0)*Kind_*rtone*pow(tone*sr, n_))/max(sqr(sr), dimensionedScalar("SMALL", sqr(sr.dimensions()), SMALL));
}


void LKTSkOmegaSST::correctApparentViscosity(const volScalarField& etta)
{
    etta_ = etta;
    bound(etta_, ettaMin_);
    etta_.correctBoundaryConditions();
}

tmp<volScalarField> LKTSkOmegaSST::strainRate() const
{
    return sqrt(2.0*magSqr(symm(fvc::grad(U_))) + Cflag_*Cbeta_*epsilon()/etta_);
}

tmp<volScalarField> LKTSkOmegaSST::calcEtta(const volScalarField& sr) const
{

    dimensionedScalar tone("tone", dimless/sr.dimensions(), 1.0);
    dimensionedScalar rtone("rtone", dimless/dimTime, 1.0);

    return min
    (
        nu0_,
        (tau0_ + Kind_*rtone*pow(tone*sr, n_))
        /max(sr, dimensionedScalar ("SMALL", sr.dimensions(), SMALL))
    );
}

tmp<volScalarField> LKTSkOmegaSST::calcNuNN(const volScalarField& sr) const
{
    return Cflag_*Cbeta_*dNudG(sr)*epsilon()/max(sr*etta_, dimensionedScalar("SMALL", sr.dimensions()*etta_.dimensions(), SMALL));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

LKTSkOmegaSST::LKTSkOmegaSST
(
    const alphaField& alpha,
    const rhoField& rho,
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

    Ce1_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Ce1",
            coeffDict_,
            2.5
        )
    ),

    Ce2_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Ce2",
            coeffDict_,
            1.85
        )
    ),

    Cbeta_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Cbeta",
            coeffDict_,
            0.667
        )
    ),

    Ctau_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Ctau",
            coeffDict_,
            0.6
        )
    ),

    Cd_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Cd",
            coeffDict_,
            0.4
        )
    ),

    Cg_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Cg",
            coeffDict_,
            0.6
        )
    ),

    NNFlag_
    (
        Switch::getOrAddToDict
        (
            "NNFlag",
            coeffDict_,
            true
        )
    ),


    Cflag_
    (
        NNFlag_.operator bool()
    ),

    alphaK1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaK1",
            coeffDict_,
            0.85
        )
    ),

    alphaK2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaK2",
            coeffDict_,
            1.0
        )
    ),

    alphaOmega1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaOmega1",
            coeffDict_,
            0.5
        )
    ),

    alphaOmega2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaOmega2",
            coeffDict_,
            0.856
        )
    ),

    gamma1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "gamma1",
            coeffDict_,
            5.0/9.0
        )
    ),

    gamma2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "gamma2",
            coeffDict_,
            0.44
        )
    ),

    beta1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "beta1",
            coeffDict_,
            0.075
        )
    ),

    beta2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "beta2",
            coeffDict_,
            0.0828
        )
    ),

    betaStar_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "betaStar",
            coeffDict_,
            0.09
        )
    ),

    a1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "a1",
            coeffDict_,
            0.31
        )
    ),

    b1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "b1",
            coeffDict_,
            1.0
        )
    ),

    c1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "c1",
            coeffDict_,
            10.0
        )
    ),

    F3_
    (
        Switch::getOrAddToDict
        (
            "F3",
            coeffDict_,
            false
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
    ),

    ettaMin_(dimensionedScalar("ettaMin", dimViscosity, SMALL)),

    nuNN_
    (
        IOobject
        (
            "nun",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar (dimViscosity, Zero)
    ),

    y_(wallDist::New(mesh_).y()),

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

    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", alphaRhoPhi.group()),
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    decayControl_
    (
        Switch::getOrAddToDict
        (
            "decayControl",
            coeffDict_,
            false
        )
    ),

    kInf_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "kInf",
            coeffDict_,
            k_.dimensions(),
            0
        )
    ),

    omegaInf_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "omegaInf",
            coeffDict_,
            omega_.dimensions(),
            0
        )
    ),

    nCorrNu_(coeffDict_.getOrDefault<label>("nInternalCorrectors", 1)),

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

    m_("m", dimless, HerschelBulkleyCoeffs_),

    Fe_("Fe", 0.5*tanh(8.0*(n_ - 0.75)) + 0.5),

    Ce_("Ce", Ce1_*Fe_ + Ce2_*(1.0 - Fe_))

{
    bound(k_, kMin_);
    bound(omega_, omegaMin_);
    bound(etta_, ettaMin_);

    setDecayControl(coeffDict_);

    if (type == typeName)
    {
        printCoeffs(type);

        Info << "Herschel-Bulkley parameters\n" << "K = " << Kind_ << endl;
        Info << "n = " << n_ << endl;
        Info << "tau0 = " << tau0_ << endl;
        Info << "m = " << m_ << endl;

        if (!NNFlag_)
        {
            Info << "\nNon-Newtonian modeling turned off in the strain-rate and all PDEs!" << endl;
            Info << "NNFlag  = " << NNFlag_ << endl;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void LKTSkOmegaSST::setDecayControl
(
    const dictionary& dict
)
{
    decayControl_.readIfPresent("decayControl", dict);

    if (decayControl_)
    {
        kInf_.read(dict);
        omegaInf_.read(dict);

        Info<< "    Employing decay control with kInf:" << kInf_
            << " and omegaInf:" << omegaInf_ << endl;
    }
    else
    {
        kInf_.value() = 0;
        omegaInf_.value() = 0;
    }
}

bool LKTSkOmegaSST::read()
{
    if (eddyViscosity<incompressible::RASModel>::read())
    {
        alphaK1_.readIfPresent(coeffDict());
        alphaK2_.readIfPresent(coeffDict());
        alphaOmega1_.readIfPresent(coeffDict());
        alphaOmega2_.readIfPresent(coeffDict());
        gamma1_.readIfPresent(coeffDict());
        gamma2_.readIfPresent(coeffDict());
        beta1_.readIfPresent(coeffDict());
        beta2_.readIfPresent(coeffDict());
        betaStar_.readIfPresent(coeffDict());
        a1_.readIfPresent(coeffDict());
        b1_.readIfPresent(coeffDict());
        c1_.readIfPresent(coeffDict());
        F3_.readIfPresent("F3", coeffDict());

        setDecayControl(coeffDict());

        Ce1_.readIfPresent(coeffDict());
        Ce2_.readIfPresent(coeffDict());
        Cbeta_.readIfPresent(coeffDict());
        Ctau_.readIfPresent(coeffDict());
        Cd_.readIfPresent(coeffDict());
        Cg_.readIfPresent(coeffDict());

        NNFlag_.readIfPresent("NNFlag", coeffDict());

        return true;
    }

    return false;
}


void LKTSkOmegaSST::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Construct local convenience references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    volScalarField& nuNN = this->nuNN_;

    tmp<volTensorField> tgradU(fvc::grad(U));
    tmp<volSymmTensorField> tS(dev(symm(tgradU)));

    volScalarField S2("S2", 2.0*(dev(tS()) && tS()));
    volScalarField G(this->GName(), nut*S2);

// - boundary condition changes a cell value
    // - normally this would be triggered through correctBoundaryConditions
    // - which would do
    //      - fvPatchField::evaluate() which calls
    //      - fvPatchField::updateCoeffs()
    // - however any processor boundary conditions already start sending
    //   at initEvaluate so would send over the old value.
    // - avoid this by explicitly calling updateCoeffs early and then
    //   only doing the boundary conditions that rely on initEvaluate
    //   (currently only coupled ones)

    //- 1. Explicitly swap values on coupled boundary conditions
    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();
    // omegaWallFunctions change the cell value! Make sure to push these to
    // coupled neighbours. Note that we want to avoid the re-updateCoeffs
    // of the wallFunctions so make sure to bypass the evaluate on
    // those patches and only do the coupled ones.
    omega_.boundaryFieldRef().template evaluateCoupled<coupledFvPatch>();

    //- 2. Make sure the boundary condition calls updateCoeffs from
    //     initEvaluate
    //     (so before any swap is done - requires all coupled bcs to be
    //      after wall bcs. Unfortunately this conflicts with cyclicACMI)
    omega_.correctBoundaryConditions();

    const volScalarField CDkOmega
    (
        (2*alphaOmega2_)*(fvc::grad(k_) & fvc::grad(omega_))/omega_
    );

    const volScalarField F1(this->F1(CDkOmega));

    const volScalarField gamma(this->gamma(F1));
    const volScalarField beta(this->beta(F1));


    volScalarField sr("sr", strainRate());

    // non-Newtonian additional viscosity for the k equation
    volScalarField DNNu("DNNu",Cd_*dNudG(sr)*S2/max(sr, dimensionedScalar("SMALL", sr.dimensions(), SMALL)));
    // non-Newtonian DN term (explicit)
    volScalarField DN("DN",fvc::laplacian(DNNu, k_));
    // non-Newtonian gamma term
    volScalarField GN("GN", -Cg_*nuNN*S2);

    //non-Newtonian E term
    volScalarField EN("EN", Ce_*(DN + GN)/nut);

    // Turbulent frequency equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(alpha, rho, omega_)
      + fvm::div(alphaRhoPhi, omega_)
      - fvm::laplacian(alpha*rho*DomegaEff(F1), omega_)
     ==
        alpha*rho*gamma*(G/nut + EN)
      - fvm::Sp(alpha*rho*beta*omega_, omega_)
      + fvm::Sp
        (
            alpha*rho*(scalar(1) - F1)*CDkOmega/omega_,
            omega_
        )
      + fvOptions(alpha, rho, omega_)
    );
    omegaEqn.ref().relax();
    fvOptions.constrain(omegaEqn.ref());
    omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
    solve(omegaEqn);
    fvOptions.correct(omega_);
    bound(omega_, omegaMin_);

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(F1), k_)
     ==
        alpha*rho*(Pk(G) + GN + DN)
      - fvm::Sp(alpha*rho*epsilonByk(), k_)
      + fvOptions(alpha, rho, k_)
    );

    tgradU.clear();
    tS.clear();

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, kMin_);

    correctApparentViscosity();
    correctNut();
    correctNuNN();
}

Foam::tmp<Foam::fvVectorMatrix>
LKTSkOmegaSST::divDevRhoReff
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

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

#include "GRkEpsilonZetaF.H"
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

defineTypeNameAndDebug(GRkEpsilonZetaF, 0);
addToRunTimeSelectionTable(RASModel, GRkEpsilonZetaF, dictionary);

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //


void GRkEpsilonZetaF::correctNut()
{
    nut_ = Cmu_*zeta_*k_*T_;
    nut_.correctBoundaryConditions();
}

void GRkEpsilonZetaF::correctApparentViscosity()
{
    volScalarField etta0 = etta_;
    for (int i = 0; i < nCorrNu_; ++i)
    {
        correctApparentViscosity(calcEtta(strainRate()));
    }

    Info << "Nu residual after " << nCorrNu_ << " iterations = " << sum(mag(etta_ - etta0)).value() << endl;

    if (tau0_.value() > 0.0 && printXi_)
    {
        dimensionedScalar tw = wallFriction();
        Info << "Yield stress ratio tau0/tauw = " << tau0_.value()/tw.value() << endl;
    }
    Info << "Average apparent viscosity nu = " << gAverage(etta_) << endl;
}

void GRkEpsilonZetaF::correctNuNN()
{
    nuNN_ = calcNuNN(strainRate());
    nuNN_.correctBoundaryConditions();
}

tmp<volScalarField> GRkEpsilonZetaF::Ts() const
{
    return
        max
        (
            min(
                    k_/epsilon(),
                    1.0/(Cmu_*sqrt(6*magSqr(dev(symm(fvc::grad(U_)))))*zeta_ + dimensionedScalar("SMALL", dimless/dimTime, SMALL))
                ),
            Ctau_*sqrt(etta_/epsilon())
        );
}

tmp<volScalarField> GRkEpsilonZetaF::Ls() const
{
    return
        CL_*max
        (
            min(
                    pow(k_, 1.5)/epsilon(),
                    sqrt(k_)/(Cmu_*sqrt(6*magSqr(dev(symm(fvc::grad(U_)))))*zeta_ + dimensionedScalar("SMALL", dimless/dimTime, SMALL))
                ),
            Ceta_*pow025(pow3(etta_)/epsilon())
        );
}

dimensionedScalar GRkEpsilonZetaF::wallFriction() const
{

    volScalarField shearMag("shearMag",
                            sqrt(2.0)*mag(symm(fvc::grad(U_)))
                            );
    shearMag.correctBoundaryConditions();

    // shear magnitude boundary field
    volScalarField::Boundary& sMBf = shearMag.boundaryFieldRef();

    // etta boundary field
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

tmp<volScalarField> GRkEpsilonZetaF::dNudG(const volScalarField& sr) const
{
    dimensionedScalar tone("tone", dimless/sr.dimensions(), 1.0);
    dimensionedScalar rtone("rtone", dimless/dimTime, 1.0);

    if (!NNFlag_)
    {
        dimensionedScalar tzero ("tzero", dimViscosity/sqr(sr.dimensions()), 0.0);
        return tzero*sr;
    }

    if (m_.value() > 0.0)
    {
        return (-tau0_*(1.0 - exp(-m_*tone*sr)) + m_*sr*tone*tau0_*exp(-m_*tone*sr) + (n_ - 1.0)*Kind_*rtone*pow(tone*sr, n_))/max(sqr(sr), dimensionedScalar("SMALL", sqr(sr.dimensions()), SMALL));
    }

    return (-tau0_  + (n_ - 1.0)*Kind_*rtone*pow(tone*sr, n_))/max(sqr(sr), dimensionedScalar("SMALL", sqr(sr.dimensions()), SMALL));
}


void GRkEpsilonZetaF::correctApparentViscosity(const volScalarField& etta)
{
    etta_ = etta;
    bound(etta_, ettaMin_);
    etta_.correctBoundaryConditions();
}

tmp<volScalarField> GRkEpsilonZetaF::strainRate() const
{
    if (NNFlag_)
    {
        return sqrt(2.0*magSqr(symm(fvc::grad(U_))) + (k_/T_)/etta_);
    }
    return sqrt(2.0*magSqr(symm(fvc::grad(U_))));
}

tmp<volScalarField> GRkEpsilonZetaF::calcEtta(const volScalarField& sr) const
{

    dimensionedScalar tone("tone", dimless/sr.dimensions(), 1.0);
    dimensionedScalar rtone("rtone", dimless/dimTime, 1.0);

    if (m_.value() > 0.0)
    {
        (1.0 - exp(-m_*tone*sr))*(tau0_ + Kind_*rtone*pow(tone*sr, n_))
        /max(sr, dimensionedScalar ("SMALL", sr.dimensions(), SMALL));
    }

    return (tau0_ + Kind_*rtone*pow(tone*sr, n_))
           /max(sr, dimensionedScalar ("SMALL", sr.dimensions(), SMALL));

}

tmp<volScalarField> GRkEpsilonZetaF::calcNuNN(const volScalarField& sr) const
{
    return Cn_*dNudG(sr)*(k_/T_)/max(sr*etta_, dimensionedScalar("SMALL", sr.dimensions()*dimViscosity, SMALL));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

GRkEpsilonZetaF::GRkEpsilonZetaF
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

    Cmu_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Cmu",
            coeffDict_,
            0.22
        )
    ),

    Ceps1a_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Ceps1a",
            coeffDict_,
            1.4
        )
    ),

    Ceps1b_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Ceps1b",
            coeffDict_,
            0.012
        )
    ),

    Ceps2_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Ceps2",
            coeffDict_,
            1.9
        )
    ),

    Cf1_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Cf1",
            coeffDict_,
            1.4
        )
    ),

    Cf2_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Cf2",
            coeffDict_,
            0.65
        )
    ),

    CL_
    (
        dimensionedScalar::getOrAddToDict
        (
            "CL",
            coeffDict_,
            0.36
        )
    ),

    Ceta_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Ceta",
            coeffDict_,
            85.0
        )
    ),

    Ctau_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Ctau",
            coeffDict_,
            6.0
        )
    ),

    sigmaK_
    (
        dimensionedScalar::getOrAddToDict
        (
            "sigmaK",
            coeffDict_,
            1.0
        )
    ),

    sigmaEps_
    (
        dimensionedScalar::getOrAddToDict
        (
            "sigmaEps",
            coeffDict_,
            1.3
        )
    ),

    sigmaZeta_
    (
        dimensionedScalar::getOrAddToDict
        (
            "sigmaZeta",
            coeffDict_,
            1.2
        )
    ),

    Cn_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Cn",
            coeffDict_,
            0.5
        )
    ),

    Cepsn_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Cepsn",
            coeffDict_,
            0.65
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

    zeta_
    (
        IOobject
        (
            IOobject::groupName("zeta", alphaRhoPhi.group()),
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    f_
    (
        IOobject
        (
            IOobject::groupName("f", alphaRhoPhi.group()),
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    T_
    (
        IOobject
        (
            "T",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh_,
        dimensionedScalar(dimTime, Zero)
    ),

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

    zetaMin_(dimensionedScalar("zetaMin", zeta_.dimensions(), SMALL)),
    fMin_(dimensionedScalar("fMin", f_.dimensions(), SMALL)),
    TMin_(dimensionedScalar("TMin", dimTime, SMALL)),
    L2Min_(dimensionedScalar("L2Min", sqr(dimLength), SMALL)),
    ettaMin_(dimensionedScalar("ettaMin", dimViscosity, SMALL)),

    epsilonZeroAtWall_(coeffDict_.getOrDefault<Switch>("epsilonZeroAtWall", true)),

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

    m_(HerschelBulkleyCoeffs_.getOrDefault<label>("m", 0.0)),

    NNFlag_(true),

    printXi_(coeffDict_.getOrDefault<label>("printYieldStressRatio", false)),

    nCorrNu_(coeffDict_.getOrDefault<label>("nInternalCorrectors", 3)),

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

    epsw_
    (
        IOobject
        (
            IOobject::groupName("epsw", alphaRhoPhi.group()),
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        2.0*etta_*magSqr(fvc::grad(sqrt(k_)))
    )

{
    bound(k_, kMin_);
    bound(epsilon_, epsilonMin_);
    bound(zeta_, zetaMin_);
    bound(f_, fMin_);

    if (type == typeName)
    {
        printCoeffs(type);

        if (!epsilonZeroAtWall_)
        {
            Info << "\n\nEpsilon not being solved with zeroed BC!" << endl;
        }

        Info << "Herschel-Bulkley parameters\n" << "K = " << Kind_ << endl;
        Info << "n = " << n_ << endl;
        Info << "tau0 = " << tau0_ << endl;

        if(m_.value() == 0.0)
        {
            Info << "Papanastasiou regularization turned off." << endl;
        }
        else
        {
            Info << "Papanastasiou regularization turned on with m = " << m_.value() << endl;
        }

        if (n_.value() == 1.0 && tau0_.value() == 0.0)
        {
            NNFlag_ = false;
            Info << "\nNewtonian fluid being used in the simulation" << endl;
            Info << "Non-Newtonian modeling turned off!" << endl;
        }
    }

    if
    (
        Kind_.value() < 0.0
        || tau0_.value() < 0.0
        || n_.value() < 0.0
        || m_.value() < 0.0
    )
    {
        FatalErrorInFunction
            << "Negative value used for K, tau0, n or m:" << nl
            << "K = " << Kind_.value() << nl
            << "tau0 = " << tau0_.value() << nl
            << "n = " << n_.value() << nl
            << "m = " << m_.value() << nl
            << exit(FatalError);
    }

    if
    (
        mag(sigmaK_.value()) < VSMALL
     || mag(sigmaEps_.value()) < VSMALL
     || mag(sigmaZeta_.value()) < VSMALL
    )
    {
        FatalErrorInFunction
            << "Non-zero values are required for the model constants:" << nl
            << "sigmaK = " << sigmaK_ << nl
            << "sigmaEps = " << sigmaEps_ << nl
            << "sigmaZeta = " << sigmaZeta_ << nl
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool GRkEpsilonZetaF::read()
{
    if (eddyViscosity<incompressible::RASModel>::read())
    {
        Cmu_.readIfPresent(coeffDict());
        Ceps1a_.readIfPresent(coeffDict());
        Ceps1b_.readIfPresent(coeffDict());
        Ceps2_.readIfPresent(coeffDict());
        Cf1_.readIfPresent(coeffDict());
        Cf2_.readIfPresent(coeffDict());
        CL_.readIfPresent(coeffDict());
        Ceta_.readIfPresent(coeffDict());
        Ctau_.readIfPresent(coeffDict());
        sigmaK_.readIfPresent(coeffDict());
        sigmaEps_.readIfPresent(coeffDict());
        sigmaZeta_.readIfPresent(coeffDict());
        Cn_.readIfPresent(coeffDict());
        Cepsn_.readIfPresent(coeffDict());

        epsilonZeroAtWall_.readIfPresent("epsilonZeroAtWall", coeffDict());
        printXi_.readIfPresent("printYieldStressRatio", coeffDict());

        return true;
    }

    return false;
}


void GRkEpsilonZetaF::correct()
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

    // epsilon value at the wall vicinity
    if (epsilonZeroAtWall_)
    {
        epsw_ = 2.0*this->etta_*magSqr(fvc::grad(sqrt(k_)));
    }

    // f at wall vicinity
    volScalarField fw("fw", 2.0*this->etta_*magSqr(fvc::grad(sqrt(zeta_))));

    T_ = Ts();
    bound(T_, TMin_);
    const volScalarField L2(type() + "L2", sqr(Ls()) + L2Min_);

    volScalarField Ceps1("Ceps1", Ceps1a_*(scalar(1.0) + Ceps1b_/zeta_));

    tmp<volSymmTensorField> tS(dev(symm(fvc::grad(U))));
    volScalarField G(this->GName(), nut*(2.0*(dev(tS()) && tS())));

    // Update epsilon and G at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();
    // Push any changed cell values to coupled neighbours
    epsilon_.boundaryFieldRef().template evaluateCoupled<coupledFvPatch>();

    // non-Newtonian additional terms
    // initialized with zeros and updated if NNFlag
    dimensionedScalar tzero ("tzero", dimless, Zero);
    // non-Newtonian additional viscosity for the k equation
    volScalarField DNNu("DNNu", tzero*etta_);
    // non-Newtonian D term (explicit)
    volScalarField DN("DN", tzero*G);
    // non-Newtonian gamma term
    volScalarField GN("GN", tzero*G);

    if (NNFlag_)
    {
        volScalarField sr("sr", strainRate());
        DNNu = Cn_*dNudG(sr)*2.0*magSqr(tS())/max(sr, dimensionedScalar("SMALL", sr.dimensions(), SMALL));
        DN = fvc::laplacian(DNNu, k_);
        GN = -2.0*nuNN_*magSqr(tS());
    }

    //- epsilon equation corrector to solve with zero BC
    volScalarField E
    (
        IOobject
        (
            "E",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(epsilon_.dimensions()/dimTime, Zero)
    );
    if (epsilonZeroAtWall_)
    {
        // correcting term that drives epsilon to zero at wall
        E = Ceps2_*epsw_/T_;
    }

    tS.clear();

    // Turbulent kinetic energy dissipation rate equation
    // k/T ~ epsilon
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon_)
      + fvm::div(alphaRhoPhi, epsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
      ==
        alpha*rho*Ceps1*(G + Cepsn_*(GN + DN))/T_
      - fvm::Sp(alpha*rho*Ceps2_/T_, epsilon_)
      + alpha*rho*E
      + fvOptions(alpha, rho, epsilon_)
    );

    epsEqn.ref().relax();
    fvOptions.constrain(epsEqn.ref());
    epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
    solve(epsEqn);
    fvOptions.correct(epsilon_);
    bound(epsilon_, epsilonMin_);

    // Turbulent kinetic energy equation
    // epsilon/k ~ 1/Ts
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*(DkEff() + DNNu), k_)
      ==
        alpha*rho*(G + GN)
      - fvm::Sp(alpha*rho*(1.0/T_), k_)
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, kMin_);

    // Elliptic relaxation function equation
    // All source terms are non-negative functions
    tmp<fvScalarMatrix> fEqn
    (
      - fvm::laplacian(f_)
      ==
      - fvm::Sp(1.0/L2, f_)
      - (
            (Cf1_ - 1.0 + Cf2_*G/epsilon_)*(zeta_ - 2.0/3.0)/T_
            - fw
        )/L2
    );
    fEqn.ref().relax();
    solve(fEqn);
    bound(f_, fMin_);

    // wall normal fluctuation ratio equation
    tmp<fvScalarMatrix> zetaEqn
    (
        fvm::ddt(alpha, rho, zeta_)
      + fvm::div(alphaRhoPhi, zeta_)
      - fvm::laplacian(alpha*rho*DZetaEff(), zeta_)
      ==
        alpha*rho*(f_ - fw)
      - fvm::SuSp
        (
            alpha*rho*
            (
                (G + GN + DN)/k_
            )
          , zeta_
        )
      + fvOptions(alpha, rho, zeta_)
    );

    zetaEqn.ref().relax();
    fvOptions.constrain(zetaEqn.ref());
    solve(zetaEqn);
    fvOptions.correct(zeta_);
    bound(zeta_, zetaMin_);

    correctNut();
    if (NNFlag_)
    {
        correctApparentViscosity();
        correctNuNN();
    }
}

Foam::tmp<Foam::fvVectorMatrix>
GRkEpsilonZetaF::divDevRhoReff
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

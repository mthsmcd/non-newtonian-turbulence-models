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
    Foam::RASModels::LKTSkOmegaSST

Group
    grpRASTurbulence

Description
    Lovato et al. low-Reynolds number k-omega SST turbulence model for
    viscoplastic fluids of the Herschel-Bulkley model.

    This turbulence model is described in:
    \verbatim
        Lovato, S., Keetels, G. H., Toxopeus, S. L., Settels, J. W. (2022).
        An eddy-viscosity model for turbulent flows of Herschel-Bulkley fluids.
        Journal of non-Newtonian fluid mechanics, 301, 104729.
    \endverbatim

    The model is built upon the Newtonian k-omega SST model:

    \verbatim
        Menter, F. (1994).
        Two-equation eddy-viscosity turbulence models for engineering applications.
        AIAA Journal, 32(8), 1598-1605.
    \endverbatim

    Note that the original implementation of the k-omega SST in OpenFOAM (v2306) had
    contained minor discrepancies described at the turbulence modeling resource webpage
    https://tmbwg.github.io/turbmodels/openfoam_issues.html

    These issues were corrected in this implementation.

    The viscoplastic k-omega SST of Lovato et al. was investigated with other viscoplastic models
    and validated against DNS in:

    \verbatim
        Macedo, M. S. S., Celis, G. E. O., Anbarlooei, H. R., Thompson, R. L.,
        Iceri, D. M., Bizarre, L., Richon, V., Castro, M. S. (2026).
        Turbulence modeling for viscoplastic fluids: state-of-the-art review
        and open-source investigation comparing with new DNS.
        Journal of non-Newtonian fluid mechanics
    \endverbatim

SourceFiles
    LKTSkOmegaSST.C

\*---------------------------------------------------------------------------*/

#ifndef LKTSkOmegaSST_H
#define LKTSkOmegaSST_H

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
                        Class LKTSkOmegaSST Declaration
\*---------------------------------------------------------------------------*/

class LKTSkOmegaSST
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

        // non-Newtonian coefficients
            dimensionedScalar Ce1_;
            dimensionedScalar Ce2_;
            dimensionedScalar Cbeta_;
            dimensionedScalar Ctau_;
            dimensionedScalar Cd_;
            dimensionedScalar Cg_;

            // flag and associated constant to turn off non-Newtonian modeling
            Switch NNFlag_;
            dimensionedScalar Cflag_;

        // non-Newtonian blending function

        // Model coefficients

            dimensionedScalar alphaK1_;
            dimensionedScalar alphaK2_;

            dimensionedScalar alphaOmega1_;
            dimensionedScalar alphaOmega2_;

            dimensionedScalar gamma1_;
            dimensionedScalar gamma2_;

            dimensionedScalar beta1_;
            dimensionedScalar beta2_;

            dimensionedScalar betaStar_;

            dimensionedScalar a1_;
            dimensionedScalar b1_;
            dimensionedScalar c1_;

            //- Flag to include the F3 term
            Switch F3_;

        // Fields

            //- non-Newtonian viscosity
            volScalarField nuNN_;

            //- Wall distance
            //  Note: different to wall distance in parent RASModel
            //  which is for near-wall cells only
            const volScalarField& y_;

            //- Turbulent kinetic energy [m2/s2]
            volScalarField k_;

            //- Specific dissipation rate field [1/s]
            volScalarField omega_;

            // Decay control
            //- Flag to include the decay control
            Switch decayControl_;
            dimensionedScalar kInf_;
            dimensionedScalar omegaInf_;

        // viscoplastic model parameters

            // numberf of internal nu calculations per iteration
            label nCorrNu_;

            // Herschel Bulkley variables
            IOdictionary transportProperties_;
            dictionary HerschelBulkleyCoeffs_;
            dimensionedScalar Kind_;
            dimensionedScalar n_;
            dimensionedScalar tau0_;
            dimensionedScalar nu0_;
            dimensionedScalar m_;

            //- Apparent viscosity
            volScalarField  etta_;

            //- Bounding value
            dimensionedScalar ettaMin_;

            // non-Newtonian blending function used to calculate Ce_
            dimensionedScalar Fe_;
            dimensionedScalar Ce_;


    // Protected Member Functions

        //- Set decay control with kInf and omegaInf
        void setDecayControl(const dictionary& dict);

        virtual tmp<volScalarField> F1(const volScalarField& CDkOmega) const;
        virtual tmp<volScalarField> F2() const;

        //- Return the blended field
        tmp<volScalarField> blend
        (
            const volScalarField& F1,
            const dimensionedScalar& psi1,
            const dimensionedScalar& psi2
        ) const
        {
            return F1*(psi1 - psi2) + psi2;
        }

        tmp<volScalarField> alphaK(const volScalarField& F1) const
        {
            return blend(F1, alphaK1_, alphaK2_);
        }

        tmp<volScalarField> alphaOmega(const volScalarField& F1) const
        {
            return blend(F1, alphaOmega1_, alphaOmega2_);
        }

        tmp<volScalarField> beta
        (
            const volScalarField& F1
        ) const
        {
            return tmp<volScalarField>::New
            (
                IOobject::scopedName(this->type(), "beta"),
                blend(F1, beta1_, beta2_)
            );
        }

        tmp<volScalarField> gamma
        (
            const volScalarField& F1
        ) const
        {
            return tmp<volScalarField>::New
            (
                IOobject::scopedName(this->type(), "gamma"),
                blend(F1, gamma1_, gamma2_)
            );
        }

        //- Return k production rate
        virtual tmp<volScalarField> Pk
        (
            const volScalarField& G
        ) const;

        //- Return epsilon/k which for standard RAS is betaStar*omega
        virtual tmp<volScalarField> epsilonByk
        (

        ) const;

        //- Return (G/nu)_0
        virtual tmp<volScalarField> GbyNu0
        (
            const volTensorField& gradU
        ) const;

        //- Return G/nu
        virtual tmp<volScalarField> GbyNu
        (
            const volScalarField& GbyNu0,
            const volScalarField& F2
        ) const;

        //- Update nut with the latest available k and omega
        virtual void correctNut();

        //- Update the apparent viscosity
        virtual void correctApparentViscosity();

        //- Update non-Newtonian additional viscosity
        virtual void correctNuNN();

public:

    //- Runtime type information
    TypeName("LKTSkOmegaSST");

    // Constructors

        //- Construct from components
        LKTSkOmegaSST
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
    virtual ~LKTSkOmegaSST() = default;


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the turbulent kinetic energy field
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        virtual tmp<volScalarField> nu() const
        {
            return etta_;
        }

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff(const volScalarField& F1) const
        {
            return tmp<volScalarField>::New
            (
                "DkEff",
                alphaK(F1)*this->nut_ + this->etta_
            );
        }

        //- Return the effective diffusivity for omega
        tmp<volScalarField> DomegaEff(const volScalarField& F1) const
        {
            return tmp<volScalarField>::New
            (
                "DomegaEff",
                alphaOmega(F1)*this->nut_ + this->etta_
            );
        }

        tmp<volScalarField> nuEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "nuEff",
                    nut_ + Ctau_*nuNN_ + this->etta_
                )
            );
        }

        virtual tmp<volScalarField> epsilon() const
        {
            return betaStar_*omega_*k_;
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

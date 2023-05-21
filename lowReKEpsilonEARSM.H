/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::RASModels::lowReKEpsilonEARSM

Description
    Standard k-epsilon turbulence model for incompressible and compressible
    flows including rapid distortion theory (RDT) based compression term.

SourceFiles
    lowReKEpsilonEARSM.C

\*---------------------------------------------------------------------------*/

#ifndef lowReKEpsilonEARSM_H
#define lowReKEpsilonEARSM_H

#include "RASModel.H"
#include "nonlinearEddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                          Class lowReKEpsilonEARSM Declaration
\*---------------------------------------------------------------------------*/
template<class BasicMomentumTransportModel>
class lowReKEpsilonEARSM
:
    public nonlinearEddyViscosity<RASModel<BasicMomentumTransportModel>>
{
protected:

    // Protected data

        // Model coefficients

    dimensionedScalar Ctau_;
    dimensionedScalar CY1prime_;
    dimensionedScalar CY2prime_;
    dimensionedScalar c1_;
    dimensionedScalar C1prime_;
    dimensionedScalar B2_;
    dimensionedScalar C1_;
    dimensionedScalar C2_;
    dimensionedScalar C3_;
    dimensionedScalar sigmak_;
    dimensionedScalar sigmaEps_;

// Fields
    volScalarField k_;
    volScalarField epsilon_;
    volSymmTensorField aij_;
    
    //- Wall distance
    const volScalarField& y_;
    //turbulent timescale
    volScalarField tau_;
    volScalarField Cmu_;


// Protected Member Functions

    virtual void correctNut();
    virtual void correctNonlinearStress(const volTensorField& gradU);
    virtual tmp<fvScalarMatrix> kSource() const;
    virtual tmp<fvScalarMatrix> epsilonSource() const;
    //tmp<volScalarField> f1() const;

public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;
    typedef typename BasicMomentumTransportModel::transportModel transportModel;


    //- Runtime type information
    TypeName("lowReKEpsilonEARSM");


    // Constructors

        //- Construct from components
        lowReKEpsilonEARSM
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& type = typeName
        );

        //- Disallow default bitwise copy construction
        lowReKEpsilonEARSM(const lowReKEpsilonEARSM&) = delete;


    //- Destructor
    virtual ~lowReKEpsilonEARSM()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return volScalarField::New
            (
                "DkEff",
                (this->nut_/sigmak_ + this->nu())
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return volScalarField::New
            (
                "DepsilonEff",
                (this->nut_/sigmaEps_ + this->nu())
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
             
        
        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
        
        
    	// Member Operators

        //- Disallow default bitwise assignment
        void operator=(const lowReKEpsilonEARSM&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "lowReKEpsilonEARSM.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

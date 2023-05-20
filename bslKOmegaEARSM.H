/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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
    Foam::RASModels::bslKOmegaEARSM

Group
    grpRASTurbulence

See also
    Foam::kOmegaSST

SourceFiles
    bslKOmegaEARSM.C

\*---------------------------------------------------------------------------*/

#ifndef bslKOmegaEARSM_H
#define bslKOmegaEARSM_H

#include "kOmegaSSTBase.H"
#include "RASModel.H"
#include "nonlinearEddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                          Class bslKOmegaEARSM Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class bslKOmegaEARSM
:
    public Foam::kOmegaSST
    <
        nonlinearEddyViscosity<RASModel<BasicMomentumTransportModel>>,
        BasicMomentumTransportModel
    >
{

protected:

	    volScalarField Cmu_;	
            dimensionedScalar alphaK1_;
            dimensionedScalar alphaOmega2_;
            dimensionedScalar gamma1_;
            dimensionedScalar gamma2_;
            dimensionedScalar Prt_;
            dimensionedScalar kappa_;
            dimensionedScalar C1_;
            dimensionedScalar A1_;

            //Fields
     	    volSymmTensorField aij_;
            //turbulent timescale
            volScalarField tau_;
            
            virtual void correctNut();
            virtual void correctNonlinearStress(const volTensorField& gradU);
 
public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;
    typedef typename BasicMomentumTransportModel::transportModel transportModel;


    //- Runtime type information
    TypeName("bslKOmegaEARSM");


    // Constructors

        //- Construct from components
        bslKOmegaEARSM
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& type = typeName
        );


    //- Destructor
    virtual ~bslKOmegaEARSM()
    {}
    
    //- Solve the turbulence equations and correct the turbulence viscosity
    virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "bslKOmegaEARSM.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //

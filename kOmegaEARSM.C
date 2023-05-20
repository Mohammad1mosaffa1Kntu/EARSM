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

\*---------------------------------------------------------------------------*/

#include "kOmegaEARSM.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "bound.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{
// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
template<class BasicMomentumTransportModel>
void kOmegaEARSM<BasicMomentumTransportModel>::correctNut()
{
    correctNonlinearStress(fvc::grad(this->U_));
    
}

template<class BasicMomentumTransportModel>
void kOmegaEARSM<BasicMomentumTransportModel>::correctNonlinearStress(const volTensorField& gradU)
{
    this->nut_ =(Cmu_/betaStar_)*k_/omega_;//TODO:check
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
    this->nonlinearStress_ = k_*aij_;   
}
template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix> kOmegaEARSM<BasicMomentumTransportModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            k_,
            dimVolume*this->rho_.dimensions()*k_.dimensions()
            /dimTime
        )
    );
}


template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix> kOmegaEARSM<BasicMomentumTransportModel>::omegaSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omega_,
            dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
        )
    );
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template<class BasicMomentumTransportModel> 
kOmegaEARSM<BasicMomentumTransportModel>::kOmegaEARSM
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& type
)
:
    nonlinearEddyViscosity<RASModel<BasicMomentumTransportModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport
    ),
    Cmu_
    (
	IOobject
	(
		"Cmu",
		this->runTime_.timeName(),
		this->mesh_,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
        this->mesh_,
        dimensionedScalar("Cmu",dimless,0.09)
    ),
    betaStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),
    beta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta",
            this->coeffDict_,
            0.072
        )
    ),
    gamma_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma",
            this->coeffDict_,
            0.52
        )
    ),
    alphaK_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK",
            this->coeffDict_,
            0.5
        )
    ),
    alphaOmega_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega",
            this->coeffDict_,
            0.5
        )
    ),
    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    A1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "A1",
            this->coeffDict_,
            1.245
        )
    ),
    aij_
    (
    IOobject
    (
        "aij",
        this->runTime_.timeName(),
        this->mesh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
        this->mesh_,
        dimensionedSymmTensor ("aij",dimless, symmTensor::zero)
    ),
    tau_
    (
        max((1./(this->betaStar_*this->omega_)),6.*sqrt( this->nu() / (this->betaStar_*this->k_*this->omega_)) ) 
    )
{
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
        Info << "kOmegaEARSM Turbulence Model." << endl;
    }
}
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class BasicMomentumTransportModel>
bool kOmegaEARSM<BasicMomentumTransportModel>::read()
{
    if (eddyViscosity<RASModel<BasicMomentumTransportModel>>::read())
    {
        betaStar_.readIfPresent(this->coeffDict());
        beta_.readIfPresent(this->coeffDict());
        gamma_.readIfPresent(this->coeffDict());
        alphaK_.readIfPresent(this->coeffDict());
        alphaOmega_.readIfPresent(this->coeffDict());
        A1_.readIfPresent(this->coeffDict());


        return true;
    }
    else
    {
        return false;
    }
}

template<class BasicMomentumTransportModel>
void kOmegaEARSM<BasicMomentumTransportModel>::correct()
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
    const volScalarField& nu = this->nu();
    const Foam::fvModels& fvModels(Foam::fvModels::New(this->mesh_)); 
    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(this->mesh_)
    );

    BasicMomentumTransportModel::correct();

    volScalarField::Internal divU(fvc::div(fvc::absolute(this->phi(), U)));

    volTensorField gradU = fvc::grad(U);


	tau_ = max((1./(this->betaStar_*this->omega_)),6.*sqrt( nu /(this->betaStar_*this->k_*this->omega_)));
	 

	#include "../EARSMCalc.H"


	volScalarField Q = 5.0/6*(sqr(N)-2*IIw)*(2*sqr(N)-IIw);
	volScalarField beta1=-(N*(2*sqr(N)-7*IIw)/Q); 
	volScalarField beta3=-12*(tr(SWW)/(Q*N));
	volScalarField beta4=-2*((sqr(N)-2*IIw)/Q);  
	volScalarField beta6=-6*(N/Q);
	volScalarField beta9= 6.0/Q;
	
	aij_ = symm(
          beta3*(WW-(1./3.)*IIw*I)
        + beta4*(SW-WS)
        + beta6*(SWW + WWS - (2./3.)*tr(SWW)*I - IIw*S)
        //+ beta9*(WSWW - WWSW + 0.5*IIw*(SW - WS)) //TODO:check
        + beta9*(WSWW-WWSW)
        );
    

    Cmu_= -0.5*(beta1 + IIw * beta6);

    volScalarField::Internal G 
    (
        this->GName(),
        (this->nut_*twoSymm(gradU) - this->nonlinearStress_) && gradU
    );

    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();

     // Turbulence specific dissipation rate equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(alpha, rho, omega_)
      + fvm::div(alphaRhoPhi, omega_)
      - fvm::laplacian(alpha*rho*DomegaEff(), omega_)
     ==
        gamma_*alpha()*rho()*G*omega_()/k_()
      - fvm::SuSp(((2.0/3.0)*gamma_)*alpha()*rho()*divU, omega_)
      - fvm::Sp(beta_*alpha()*rho()*omega_(), omega_)
      + omegaSource()
      + fvModels.source(alpha, rho, omega_)
    );

    omegaEqn.ref().relax();
    fvConstraints.constrain(omegaEqn.ref());
    omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
    solve(omegaEqn);
    fvConstraints.constrain(omega_);
    bound(omega_, this->omegaMin_);
        
    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha()*rho()*G
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
      - fvm::Sp(betaStar_*alpha()*rho()*omega_(), k_)
      + kSource()
      + fvModels.source(alpha, rho, k_)
    );


    kEqn.ref().relax();
    fvConstraints.constrain(kEqn.ref());
    solve(kEqn);
    fvConstraints.constrain(k_);
    bound(k_, this->kMin_);


    // Re-calculate viscosity and non-linear stress
    correctNonlinearStress(gradU);
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

//

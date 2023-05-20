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

#include "bslKOmegaEARSM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{
// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
template<class BasicMomentumTransportModel>
void bslKOmegaEARSM<BasicMomentumTransportModel>::correctNut()
{
    correctNonlinearStress(fvc::grad(this->U_));
    
}

template<class BasicMomentumTransportModel>
void bslKOmegaEARSM<BasicMomentumTransportModel>::correctNonlinearStress(const volTensorField& gradU)
{
    this->nut_ = (Cmu_/this->betaStar_)*this->k_/this->omega_;
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
    this->nonlinearStress_ = this->k_*aij_;

}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template<class BasicMomentumTransportModel> 
bslKOmegaEARSM<BasicMomentumTransportModel>::bslKOmegaEARSM
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
    Foam::kOmegaSST 
    <
        nonlinearEddyViscosity<RASModel<BasicMomentumTransportModel>>,
        BasicMomentumTransportModel
    >
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
    alphaK1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK1",
            this->coeffDict_,
            0.5
        )
    ),
    alphaOmega2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega2",
            this->coeffDict_,
            0.85616
        )
    ),
    gamma1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma1",
            this->coeffDict_,
            0.5532
        )
    ),
    gamma2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma2",
            this->coeffDict_,
            0.4403
        )
    ),
    Prt_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Prt",
            this->coeffDict_,
            1.0
        )
    ),
    kappa_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kappa",
            this->coeffDict_,
            0.41
        )
    ),
    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            this->coeffDict_,
            1.8
        )
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
    if (type == typeName)
    {
        this->printCoeffs(type);
        Info << "BSL kOmega EARSM Turbulence Model." << endl;
    }
}
template<class BasicMomentumTransportModel>
void bslKOmegaEARSM<BasicMomentumTransportModel>::correct()
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

    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(this->mesh_)
    );

    BasicMomentumTransportModel::correct();

    volScalarField::Internal divU(fvc::div(fvc::absolute(this->phi(), U)));

    volTensorField gradU = fvc::grad(U);//dev();
	tau_ = max((1./(this->betaStar_*this->omega_)),6.*sqrt( nu /(this->betaStar_*this->k_*this->omega_)));
	

	#include "../EARSMCalc.H"

        volScalarField Q = 5.0/6*(sqr(N)-2*IIw)*(2*sqr(N)-IIw);
	volScalarField beta1=-(N*(2*sqr(N)-7*IIw)/Q); 
	volScalarField beta3=-12*(tr(SWW)/(Q*N));
	volScalarField beta4=-2*((sqr(N)-2*IIw)/Q);  
	volScalarField beta6=-6*(N/Q);
	volScalarField beta9= 6.0/Q;
	
	aij_ = symm(
	  //beta1*S
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
    this->omega_.boundaryFieldRef().updateCoeffs();

    
    volScalarField CDkOmega
    (
        (2*this->alphaOmega2_)*(fvc::grad(this->k_) & fvc::grad(this->omega_))/this->omega_
    );

    volScalarField F1(this->F1(CDkOmega));
    {
        volScalarField::Internal gamma(this->gamma(F1));
        volScalarField::Internal beta(this->beta(F1));

        tmp<fvScalarMatrix> omegaEqn
        (
            fvm::ddt(alpha, rho, this->omega_)
          + fvm::div(alphaRhoPhi, this->omega_)
          - fvm::laplacian(alpha*rho*this->DomegaEff(F1), this->omega_) 
          ==
		    alpha()*rho()*gamma
           *min
            (
                G,
                this->c1_*this->betaStar_*this->k_()*this->omega_() 
            )*this->omega_()/this->k_()
		  - fvm::SuSp((2.0/3.0)*alpha()*rho()*gamma*divU, this->omega_)
          - fvm::Sp(alpha()*rho()*beta*this->omega_(), this->omega_)
          - fvm::SuSp 
            (
                alpha()*rho()*(F1() - scalar(1))*CDkOmega()/this->omega_(),
                this->omega_
            )
        );

        omegaEqn.ref().relax();
        fvConstraints.constrain(omegaEqn.ref());
        omegaEqn.ref().boundaryManipulate(this->omega_.boundaryFieldRef());
        solve(omegaEqn);
        fvConstraints.constrain(this->omega_);
        bound(this->omega_, this->omegaMin_);
    }

    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, this->k_)
      + fvm::div(alphaRhoPhi, this->k_)
      - fvm::laplacian(alpha*rho*this->DkEff(F1), this->k_) 
     ==
        min
	    (
		alpha()*rho()*G, 
		(this->c1_*this->betaStar_)*alpha()*rho()*this->k_()*this->omega_() //TODO:check
	    )
	  - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, this->k_)
      - fvm::Sp(alpha()*rho()*this->betaStar_*this->omega_(),this->k_)
    );

    kEqn.ref().relax();
    fvConstraints.constrain(kEqn.ref());
    solve(kEqn);
    fvConstraints.constrain(this->k_);
    bound(this->k_, this->kMin_);

    // Re-calculate viscosity and non-linear stress
    correctNonlinearStress(gradU);
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

//

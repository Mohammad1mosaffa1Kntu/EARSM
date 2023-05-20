/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
#include "ChienKEpsilonEARSM.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "bound.H"
#include "wallDist.H"
#include "wallFvPatch.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Static Data Members * * * * * * * * * * //


// * * * * * * * * * * Private Member Functions * * * * * * * * * //
template<class BasicMomentumTransportModel>
void ChienKEpsilonEARSM<BasicMomentumTransportModel>::correctNut()
{
    correctNonlinearStress(fvc::grad(this->U_));
    
}

template<class BasicMomentumTransportModel>
void ChienKEpsilonEARSM<BasicMomentumTransportModel>::correctNonlinearStress(const volTensorField& gradU)
{
    //this->nut_ = Cmu_*sqr(k_)/epsilon_;
    this->nut_ = Cmu_*k_*tau_;
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
    this->nonlinearStress_ = k_*aij_;

}

/*template<class BasicMomentumTransportModel>
tmp<volScalarField> ChienKEpsilonEARSM<BasicMomentumTransportModel>::fMu() const
{
	volScalarField y=wallDist(this->mesh_).y();
	volScalarField Rey = pow(k_,0.5)*(y/this->nu());
	volScalarField yStar=pow(Rey,2)*0.003 + pow(Rey, 0.5)*2.4;
	return	scalar(1)- exp(-0.0115*yStar);
}*/
template<class BasicMomentumTransportModel>
tmp<volScalarField> ChienKEpsilonEARSM<BasicMomentumTransportModel>::f2() const
{
	return	scalar(1)- 0.22*exp(-sqr(sqr(k_)/(this->nu()*epsilonTilda_*6))); //TODO:check
}
// * * * * * * * * * * * * * * Constructors * * * * * * * * * * * //
template<class BasicMomentumTransportModel>
ChienKEpsilonEARSM<BasicMomentumTransportModel>::ChienKEpsilonEARSM
(
	const alphaField& alpha,
	const rhoField& rho,
	const volVectorField& U,
	const surfaceScalarField& alphaRhoPhi,
	const surfaceScalarField& phi,
	const transportModel& transport,
	const word& type 
):
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
Cmu_(
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
C1_(
	dimensioned<scalar>::lookupOrAddToDict
	(
		"C1",
		this->coeffDict_,
		1.35
	)
	),
C2_(
	dimensioned<scalar>::lookupOrAddToDict
	(
		"C2",
		this->coeffDict_,
		1.8
	)
	),
sigmaEps_(
	dimensioned<scalar>::lookupOrAddToDict
	(
		"sigmaEps",
		this->coeffDict_,
		1.3
	)
	),
sigmak_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmak",
            this->coeffDict_,
            1.0
        )
    ),
k_(
	IOobject
	(
		"k",
		this->runTime_.timeName(),
		this->mesh_,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
	this->mesh_
	),
epsilonTilda_(
	IOobject
	(
		"epsilon",
		this->runTime_.timeName(),
		this->mesh_,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
	this->mesh_
	),
aij_(
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
tau_(
	IOobject
	(
		"tau",
		this->runTime_.timeName(),
		this->mesh_,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	k_/epsilonTilda_
)
{
	bound(k_, this->kMin_);
	bound(epsilonTilda_, this->epsilonMin_);
	
	if (type == typeName)
    	{
        	this->printCoeffs(type);
        }
}
// * * * * * * * * * * * * * Member Functions * * * * * * * * * * //
template<class BasicMomentumTransportModel>
bool ChienKEpsilonEARSM<BasicMomentumTransportModel>::read()
{
if (nonlinearEddyViscosity<RASModel<BasicMomentumTransportModel>>::read())
{
	C1_.readIfPresent(this->coeffDict());
	C2_.readIfPresent(this->coeffDict());
	sigmaEps_.readIfPresent(this->coeffDict());
	sigmak_.readIfPresent(this->coeffDict());
	return true;
}
else
{
	return false;
}
}
template<class BasicMomentumTransportModel>
void ChienKEpsilonEARSM<BasicMomentumTransportModel>::correct()
{
	if (!this->turbulence_)
	{
		return;
	}
	const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    const volScalarField& nu = this->nu();
    volTensorField gradU = fvc::grad(U); 

	tau_ = max(k_/epsilonTilda_, 6.0*sqrt(nu/epsilonTilda_));
	
	#include "../EARSMCalc.H"
	
	volScalarField Q = 5.0/6*(sqr(N)-2*IIw)*(2*sqr(N)-IIw);
	volScalarField beta_1=-(N*(2*sqr(N)-7*IIw)/Q); 
	volScalarField beta_3=-12*(tr(SWW)/(Q*N));
	volScalarField beta_4=-2*((sqr(N)-2*IIw)/Q);  
	volScalarField beta_6=-6*(N/Q);
	volScalarField beta_9= 6.0/Q;

		
	volScalarField y=wallDist(this->mesh_).y();
	volScalarField Re_t = pow(k_,0.5)*(y/nu);
	volScalarField yStar = 2.4*pow(Re_t,0.5) + 0.003*pow(Re_t,2.0);
	volScalarField f1= 1.0 - exp(-yStar/26.0);
	
	volScalarField beta_2_lowre=(3.0*1.8-4.0)/(max(IIs,5.74))*(1-pow(f1,2));    
	volScalarField beta_4_lowre = pow(f1,2)*beta_4 - 1.8/(2*max(IIs,5.74))*(1-pow(f1,2)); 
	

	aij_= symm( 
		//f1*beta_1*S 
		//+ 
		  beta_2_lowre*(SS - (1.0/3.0)*IIs * I)
		+ pow(f1,2)*beta_3*(WW-(1.0/3.0)*IIw * I)
		+ beta_4_lowre*(SW-WS)
		+ f1*beta_6*(SWW+WWS-2.0/3*tr(SWW)*I)
		+ pow(f1,2)*beta_9*(WSWW-WWSW)
		);
	Info << "test new" <<endl;
	
	Cmu_=-f1*(beta_1 + IIw * beta_6)/2.0;
	
	volScalarField::Internal G //TODO: check gradU or T(gradU)
    (
        this->GName(),
        (this->nut_*twoSymm(gradU) - this->nonlinearStress_) && gradU
    );                    
	//volScalarField G(this->GName(), -R_ && T(gradU));
	//volSymmTensorField P(-twoSymm(R_ & gradU));
    //volScalarField G(this->GName(), 0.5*mag(tr(P)));

	
	//calculate yPlus
	volScalarField yPlus 
	(
		IOobject
		(	
			"yPlus",
			this->runTime_.timeName(),
			this->mesh_,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
        	this->mesh_,
        	dimensionedScalar ("yPlus",dimless, 0)
	);
	

	const fvPatchList& patches = this->mesh_.boundary();
	forAll(patches, patchi)
    {
		const fvPatch& patch = patches[patchi];
		
		if (isA<wallFvPatch>(patch))
		{
			yPlus.boundaryFieldRef()[patchi] = y.boundaryField()[patchi]*sqrt
			(
				nu.boundaryField()[patchi]
				*mag(U.boundaryField()[patchi].snGrad())
			)/nu.boundaryField()[patchi];
		}
    }

	const volScalarField E(-2.0*nu*epsilonTilda_/pow(y,2)*exp(-0.5*yStar)*exp(-0.04*yPlus)); //TODO:CHECK
	const volScalarField D(2.0*nu*k_/pow(y,2)*exp(-0.04*yPlus)); //TODO:CHECK

// Dissipation rate equation
 	
 	// Update epsilon and G at the wall
    epsilonTilda_.boundaryFieldRef().updateCoeffs();
    
	tmp<fvScalarMatrix> epsEqn
	(
	    fvm::ddt(alpha, rho, epsilonTilda_)
      + fvm::div(alphaRhoPhi, epsilonTilda_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilonTilda_)
      ==
        C1_*alpha()*rho()*G*epsilonTilda_/k_
      - fvm::Sp(C2_*f2()*alpha()*rho()*epsilonTilda_/k_, epsilonTilda_) 
      + alpha()*rho()*E //TODO:check
    );
	epsEqn.ref().relax(); 
	solve(epsEqn);
	bound(epsilonTilda_, this->epsilonMin_);

// Turbulent kinetic energy equation
	tmp<fvScalarMatrix> kEqn
	(
		fvm::ddt(alpha,rho,k_)
		+ fvm::div(alphaRhoPhi, k_)
		- fvm::laplacian(alpha*rho*DkEff(), k_)
		==
		alpha()*rho()*G 
		- fvm::Sp(alpha()*rho()*(epsilonTilda_ + D)/k_, k_) //TODO:check
    );
	kEqn.ref().relax(); 
	solve(kEqn);
	bound(k_, this->kMin_);
	
	// Re-calculate viscosity and non-linear stress
    correctNonlinearStress(gradU); 
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace RASModels_
//} // End namespace incompressible
} // End namespace Foam
// *************************************************************** //

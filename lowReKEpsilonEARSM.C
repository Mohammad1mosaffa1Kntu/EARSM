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

#include "lowReKEpsilonEARSM.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "bound.H"
#include "wallDist.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
namespace RASModels
{


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
template<class BasicMomentumTransportModel>
void lowReKEpsilonEARSM<BasicMomentumTransportModel>::correctNut()
{
    correctNonlinearStress(fvc::grad(this->U_));
    
}

template<class BasicMomentumTransportModel>
void lowReKEpsilonEARSM<BasicMomentumTransportModel>::correctNonlinearStress(const volTensorField& gradU)
{
    //this->nut_ = Cmu_*sqr(k_)/epsilon_;
    this->nut_ = Cmu_*k_*tau_;
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
    this->nonlinearStress_ = k_*aij_;

}
template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix> lowReKEpsilonEARSM<BasicMomentumTransportModel>::kSource() const
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
tmp<fvScalarMatrix> lowReKEpsilonEARSM<BasicMomentumTransportModel>::epsilonSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            epsilon_,
            dimVolume*this->rho_.dimensions()*epsilon_.dimensions()
            /dimTime
        )
    );
}
/*template<class BasicMomentumTransportModel>
tmp<volScalarField> lowReKEpsilonEARSM<BasicMomentumTransportModel>::f1() const
{
    tmp<volScalarField> Rey = sqr(k_)*y_/(this->nu());

    return scalar(1.0) - exp(-1.0/26*(CY1prime_*sqrt(Rey)+ CY2prime_*sqr(Rey)));	
    
}*/


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
lowReKEpsilonEARSM<BasicMomentumTransportModel>::lowReKEpsilonEARSM
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

    Ctau_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ctau",
            this->coeffDict_,
            6.0
        )
    ),
    CY1prime_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CY1prime",
            this->coeffDict_,
            2.4 
        )
    ),
    CY2prime_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CY2prime",
            this->coeffDict_,
            0.003
        )
    ),
    c1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "c1",
            this->coeffDict_,
            1.8
        )
    ),
    C1prime_ //TODO:check
    (
    	dimensioned<scalar>::lookupOrAddToDict
        (
            "C1prime_",
            this->coeffDict_,
            (9.0/4*(c1_ - 1.0)).value()
        )
    ),
    B2_ 
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "B2",
            this->coeffDict_,
	    1.8
        )
    ),
    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            this->coeffDict_,
            1.44
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            this->coeffDict_,
            1.92
        )
    ),
    C3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C3",
            this->coeffDict_,
            0
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
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            this->coeffDict_,
            1.3
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
    epsilon_
    (
        IOobject
        (
            IOobject::groupName("epsilon", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    aij_
    (
        IOobject
        (
            IOobject::groupName("aij", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedSymmTensor("aij",dimless,symmTensor::zero)
    ),
    y_(wallDist::New(this->mesh_).y()),
    tau_
    (
        IOobject
        (
            IOobject::groupName("tau", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        k_/epsilon_
    ),
    Cmu_
    (
        IOobject
        (
            IOobject::groupName("Cmu", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("Cmu_",dimless,0.09)
    )
{
    bound(k_, this->kMin_);
    bound(epsilon_, this->epsilonMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class BasicMomentumTransportModel>
bool lowReKEpsilonEARSM<BasicMomentumTransportModel>::read()
{
    if (nonlinearEddyViscosity<RASModel<BasicMomentumTransportModel>>::read())
    {
        C1_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());
        C3_.readIfPresent(this->coeffDict());
        c1_.readIfPresent(this->coeffDict());
        B2_.readIfPresent(this->coeffDict());
        CY1prime_.readIfPresent(this->coeffDict());
        CY2prime_.readIfPresent(this->coeffDict());
        sigmak_.readIfPresent(this->coeffDict());
        sigmaEps_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
void lowReKEpsilonEARSM<BasicMomentumTransportModel>::correct()
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
    volTensorField gradU = fvc::grad(U); 
    
 
    const Foam::fvModels& fvModels(Foam::fvModels::New(this->mesh_));
    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(this->mesh_)
    );

    nonlinearEddyViscosity<RASModel<BasicMomentumTransportModel>>::correct();

    volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))()
    );

    //tau_ = max(k_/epsilon_,Ctau_*sqrt(this->nu()/epsilon_));
    tau_ = max(k_/epsilon_, 6.0*sqrt(nu/epsilon_));
	
	#include "../EARSMCalc.H"
	
	volScalarField Q = 5.0/6*(sqr(N)-2*IIw)*(2*sqr(N)-IIw);
	volScalarField beta1=-(N*(2*sqr(N)-7*IIw)/Q); 
	volScalarField beta3=-12*(tr(SWW)/(Q*N));
	volScalarField beta4=-2*((sqr(N)-2*IIw)/Q);  
	volScalarField beta6=-6*(N/Q);
	volScalarField beta9= 6.0/Q;

        
        volScalarField Rey = sqrt(k_)*y_/(this->nu());
   	volScalarField f1(scalar(1.0) - exp(-1.0/26*(CY1prime_*sqrt(Rey)+ CY2prime_*sqr(Rey))));
	//volScalarField Re_t = pow(k_,0.5)*(y_/nu);
	//volScalarField yStar = 2.4*pow(Re_t,0.5) + 0.003*pow(Re_t,2.0);
	//volScalarField f1= 1.0 - exp(-yStar/26.0);
	
	volScalarField beta2_lowre=(3.0*1.8-4.0)/(max(IIs,5.74))*(1-pow(f1,2));    
	volScalarField beta4_lowre = pow(f1,2)*beta4 - 1.8/(2*max(IIs,5.74))*(1-pow(f1,2)); 
	

	aij_= symm( 
		//f1*beta1*S 
		//+ 
		  beta2_lowre*(SS - (1.0/3.0)*IIs * I)
		+ pow(f1,2)*beta3*(WW-(1.0/3.0)*IIw * I)
		+ beta4_lowre*(SW-WS)
		+ f1*beta6*(SWW+WWS-2.0/3*tr(SWW)*I)
		+ pow(f1,2)*beta9*(WSWW-WWSW)
		);
	
	
	Cmu_=-0.5*f1*(beta1 + IIw * beta6);
	
	//volScalarField::Internal G(this->GName(), -R_ && T(gradU));
    volScalarField::Internal G //TODO: check gradU or T(gradU)
    (
        this->GName(),
        (this->nut_*twoSymm(gradU) - this->nonlinearStress_) && gradU
    );

    // Update epsilon and G at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon_)
      + fvm::div(alphaRhoPhi, epsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
     ==
        C1_*alpha()*rho()*G*epsilon_()/k_()
      - fvm::SuSp(((2.0/3.0)*C1_ - C3_)*alpha()*rho()*divU, epsilon_)
      - fvm::Sp(C2_*alpha()*rho()*epsilon_()/k_(), epsilon_)
      + epsilonSource()
      + fvModels.source(alpha, rho, epsilon_)
    );

    epsEqn.ref().relax();
    fvConstraints.constrain(epsEqn.ref());
    epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
    solve(epsEqn);
    fvConstraints.constrain(epsilon_);
    bound(epsilon_, this->epsilonMin_);

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha()*rho()*G
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
      - fvm::Sp(alpha()*rho()*epsilon_()/k_(), k_)
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

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "EDC.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::combustionModels::EDC<Type>::EDC
(
    const word& modelType,
    const fvMesh& mesh,
    const word& phaseName
)
:
    laminar<Type>(modelType, mesh, phaseName),
    Cepsilon_(readScalar(this->coeffs().lookup("Cepsilon"))),
    Ctau_(readScalar(this->coeffs().lookup("Ctau"))),
    turbulentReaction_(this->coeffs().lookup("turbulentReaction")),
    tau_
    (
        IOobject
        (
            IOobject::groupName("EDC:tau", phaseName),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("tau", dimless, 0.0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::combustionModels::EDC<Type>::~EDC()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
void Foam::combustionModels::EDC<Type>::correct()
{
    if (this->active())
    {
        laminar<Type>::correct();

        if (turbulentReaction_)
        {
            tmp<volScalarField> tepsilon(this->turbulence().epsilon());
            const volScalarField& epsilon = tepsilon();
            tmp<volScalarField> tk(this->turbulence().k());
            const volScalarField& k = tk();
            tmp<volScalarField> tmuEff(this->turbulence().muEff());
            const volScalarField& muEff = tmuEff();
            tmp<volScalarField> trho(this->rho());
            const volScalarField& rho = trho();

            forAll(epsilon, i)
            {
                scalar epsilonStar =
                    Cepsilon_*sqrt(sqrt(muEff[i]/rho[i]*epsilon[i]/k[i]));

                scalar tauStar = 
                    Ctau_**sqrt(max(muEff[i]/rho[i]/(epsilon[i] + SMALL), 0));
            }
        }
    }
}


template<class Type>
Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::EDC<Type>::R(volScalarField& Y) const
{
    // R = rho*epsilonStar*epsilonStar / tauStar / (1-pow(epsilonStar,3) * (Ystar - Y)
    return laminar<Type>::R(Y);
}


template<class Type>
Foam::tmp<Foam::volScalarField>
Foam::combustionModels::EDC<Type>::dQ() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject::groupName("EDC:dQ", this->phaseName_),
            laminar<Type>::dQ()
        )
    );
}


template<class Type>
Foam::tmp<Foam::volScalarField>
Foam::combustionModels::EDC<Type>::Sh() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject::groupName("EDC:Sh", this->phaseName_),
            laminar<Type>::Sh()
        )
    );
}


template<class Type>
bool Foam::combustionModels::EDC<Type>::read()
{
    if (laminar<Type>::read())
    {
        this->coeffs().lookup("Ctau") >> Ctau_;
        this->coeffs().lookup("Cepsilon") >> Cepsilon_;
        this->coeffs().lookup("turbulentReaction") >> turbulentReaction_;
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //

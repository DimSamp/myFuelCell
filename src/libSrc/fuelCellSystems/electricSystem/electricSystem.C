/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by the original author
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

#include "electricSystem.H"
#include "regionTypeList.H"
#include "dissolvedModel.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(electricSystem, 0);
}

const Foam::word Foam::electricSystem::dictName("electricProperties");


const Foam::Enum<Foam::electricSystem::algorithmType>
Foam::electricSystem::algorithmTypeNames
({
    { algorithmType::segregated, "segregated" },
    { algorithmType::coupled, "coupled" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::regionTypes::electric& Foam::electricSystem::lookupRegion
(
    regionTypeList& regions,
    const word& name
)
{
    forAll(regions, i)
    {
        if (regions[i].name() == name)
        {
            if (!isA<regionTypes::electric>(regions[i]))
            {
                FatalErrorInFunction
                    << "Region " << name << " is not an electric region"
                    << exit(FatalError);
            }

            return refCast<regionTypes::electric>(regions[i]);
        }
    }

    FatalErrorInFunction
        << "Cannot find region " << name << " (" << dictName << ")"
        << exit(FatalError);

    return refCast<regionTypes::electric>(regions[0]);
}


void Foam::electricSystem::solveRegion
(
    regionTypes::electric& r,
    const scalar refRelax
)
{
    Info << "\nSolve for region " << r.name() << ":\n" << endl;

    tmp<fvScalarMatrix> phiEqn
    (
      - fvm::laplacian(r.sigmaField_, r.phi_, "laplacian(sigma,phi)")
      - r.j_
    );

    //- Set reference values
    if (r.phi_.needReference())
    {
        if (&r != &iR_)
        {
            FatalErrorInFunction
                << "The potential of the electron region " << r.name()
                << " has no fixed-value patch"
                << exit(FatalError);
        }

        //- Update the potential field.
        //- Electro-Neutral, the total electron+proton flux should be zero
        const scalarField& source = r.j_;
        const scalarField& volume = r.V();
        scalarField sum = source*volume;

        scalar iDot(Foam::gSum(sum)/Foam::gSum(volume));

        for (label id = 0; id < r.nCells(); id++)
        {
            //- The option has to be true (forceReference)
            //- otherwise, errors happen in parallel simulations
            phiEqn->setReference(id, r.phi_[id] - iDot*refRelax, true);
        }
    }

    phiEqn->solve();

    r.i_ = -r.sigmaField_ * fvc::grad(r.phi_);
    r.i_.correctBoundaryConditions();

    if (r.dissolved_.valid())
    {
        r.dissolved_->solve();
    }
}


void Foam::electricSystem::solveSegregated()
{
    solveRegion(eR_, 0);
    solveRegion(iR_, ionReferenceRelax_);
}


void Foam::electricSystem::solveCoupled()
{
    //- Unknowns [phiE; phiI] in one matrix.
    //- J of both regions holds the base-point current J* (= realised
    //- current of the previous iteration), ButlerVolmer has written
    //- dJdEta = G'(eta*) and etaE = eta* + E in the CL cells.
    Info<< "\nCoupled solve for regions " << eR_.name() << " + " << iR_.name()
        << ":\n" << endl;

    const scalarField& etaE = etaE_();

    matrix_->assemble(dJdEta_(), etaE);

    scalarField psi(matrix_->psi());

    solverPerformance perf = matrix_->solve
    (
        psi,
        matrix_->source(),
        iR_.solverDict("phiCoupled")
    );

    Info<< perf.solverName() << ":  Solving for phiCoupled"
        << ", Initial residual = " << perf.initialResidual()
        << ", Final residual = " << perf.finalResidual()
        << ", No Iterations " << perf.nIterations() << endl;

    matrix_->scatter(psi);

    control_->afterSolve(matrix_());

    //- Realised current, base point of the next iteration
    matrix_->realise(etaE);

    //- Current densities and dissolved model, as in the segregated solve
    eR_.i_ = -eR_.sigmaField_*fvc::grad(eR_.phi_);
    eR_.i_.correctBoundaryConditions();
    iR_.i_ = -iR_.sigmaField_*fvc::grad(iR_.phi_);
    iR_.i_.correctBoundaryConditions();

    if (eR_.dissolved_.valid())
    {
        eR_.dissolved_->solve();
    }
    if (iR_.dissolved_.valid())
    {
        iR_.dissolved_->solve();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::electricSystem::electricSystem
(
    const Time& runTime,
    regionTypeList& regions
)
:
    IOdictionary
    (
        IOobject
        (
            dictName,
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    eR_(lookupRegion(regions, get<word>("electronRegion"))),
    iR_(lookupRegion(regions, get<word>("ionRegion"))),
    algorithm_(algorithmTypeNames.get("algorithm", *this)),
    ionReferenceRelax_(0),
    omega_(1)
{
    Info<< nl << "Electric system: " << eR_.name() << " + " << iR_.name()
        << ", algorithm " << algorithmTypeNames[algorithm_] << nl;

    if (coupled())
    {
        omega_ =
            subDict("coupledCoeffs").getOrDefault<scalar>("omega", 1.0);

        if (omega_ <= 0 || omega_ > 1)
        {
            FatalIOErrorInFunction(subDict("coupledCoeffs"))
                << "omega = " << omega_ << " is not in (0, 1]"
                << exit(FatalIOError);
        }

        Info<< "    omega = " << omega_ << nl;
    }
    else
    {
        ionReferenceRelax_ =
            subDict("segregatedCoeffs").get<scalar>("ionReferenceRelax");

        Info<< "    ionReferenceRelax = " << ionReferenceRelax_ << nl;
    }

    control_.reset(new electricControl(eR_, subDict("control")));

    if (control_->needsCoupled() && !coupled())
    {
        FatalIOErrorInFunction(*this)
            << "The galvanostatic mode needs 'algorithm coupled'"
            << exit(FatalIOError);
    }

    if (coupled())
    {
        //- Their values are set by ButlerVolmer
        dJdEta_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "dJdEta",
                    iR_.time().timeName(),
                    iR_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                iR_,
                dimensionedScalar(iR_.j_.dimensions()/iR_.phi_.dimensions(), Zero)
            )
        );

        etaE_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "etaE",
                    iR_.time().timeName(),
                    iR_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                iR_,
                dimensionedScalar(iR_.phi_.dimensions(), Zero)
            )
        );

        matrix_.reset(new coupledMatrix(eR_, iR_, omega_));
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::electricSystem::correct()
{
    control_->beforeSolve();
}


void Foam::electricSystem::solve()
{
    if (coupled())
    {
        solveCoupled();
    }
    else
    {
        solveSegregated();
    }
}


// ************************************************************************* //

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

#include "ButlerVolmer.H"
#include "butlerVolmerInverse.H"
#include "electricSystem.H"
#include "phaseModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::activationOverpotentialModels::ButlerVolmer<Thermo>::ButlerVolmer
(
    const phaseModel& phase,
    const dictionary& dict
)
:
    ActivationOverpotentialModel<Thermo>(phase, dict),
    relaxNoted_(false),
    inverseTol_(dict.template getOrDefault<scalar>("inverseTolerance", 1e-12)),
    inverseMaxIter_(dict.template getOrDefault<label>("inverseMaxIter", 100)),
    kinOld_(phase.mesh().nCells(), Zero)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Thermo>
Foam::activationOverpotentialModels::ButlerVolmer<Thermo>::~ButlerVolmer()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class Thermo>
void Foam::activationOverpotentialModels::ButlerVolmer<Thermo>::correct()
{
    //- Create Nernst model
    this->createNernst();

    //- Correct
    this->nernst_->correct();

    //- Get sub regions
    //- Refer to regionType
    //- Including: fluid, electron (BPP + GDL + CL), ion (CLs + membrane)
    const regionType& fluidPhase = this->region
    (
        word(this->regions_.subDict("fluid").lookup("name"))
    );
    const regionType& electronPhase = this->region
    (
        word(this->regions_.subDict("electron").lookup("name"))
    );
    const regionType& ionPhase = this->region
    (
        word(this->regions_.subDict("ion").lookup("name"))
    );

    //- Source/sink terms for electron and ion fields
    //- See names in regions/electronIon/electronIon.C
    scalarField& SE = const_cast<volScalarField&>
    (
        electronPhase.template lookupObject<volScalarField>("J")
    );
    scalarField& SI = const_cast<volScalarField&>
    (
        ionPhase.template lookupObject<volScalarField>("J")
    );

    //- Potential fields
    const scalarField& phiE = electronPhase.template
        lookupObject<volScalarField>
        (
            word(this->phiNames_["electron"])
        );
    const scalarField& phiI = ionPhase.template
        lookupObject<volScalarField>
        (
            word(this->phiNames_["ion"])
        );

    //- Reference
    scalarField& eta = this->eta_;
    //- Nernst field
    scalarField& nernst = this->nernst_()();
    //- Current density
    scalarField& j = this->j_;

    //- Access temperature and reactant
    const scalarField& T = this->thermo_.T();
    const scalarField& s = this->phase_;

    //- Mixture mole fraction
    const scalarField W(this->thermo_.W()/1000.0);
    //- Mixture density
    const scalarField& rho= this->phase_.thermo().rho();
    //- Store values for all species

    scalarField coeff(this->phase_.mesh().nCells(), 1.0);
    //- Loop for all relative species
    forAllConstIter(dictionary, this->species_, iter)
    {
        if (iter().isDict())
        {
            const word& name = iter().keyword();
            const dictionary& dict0 = iter().dict();

            scalar ksi = readScalar(dict0.lookup("ksi"));
            scalar cRef = readScalar(dict0.lookup("cRef"));

            const scalarField& X = this->phase_.X(name);

            coeff *= Foam::pow(X*rho/W/cRef, ksi);
        }
    }

    //- Only consider catalyst zone
    label znId = fluidPhase.cellZones().findZoneID(this->zoneName_);
    const labelList& cells = fluidPhase.cellZones()[znId];

    //- electron transfer
    //- cathode side, < 0
    //- anode side, >0
    scalar n = this->nernst_->rxnList()["e"];
    scalar sign = n/mag(n);

    //- The total current: volume averaged
    scalar Rj(0.0);

    //- Coupled phiE-phiI solve (see electricSystem, coupledMatrix): the
    //- electric system then holds dJdEta and etaE. Instead of evaluating j
    //- from the lagged potentials, the base point is rebuilt from the
    //- realised current of the last coupled solve, SI, rescaled with the
    //- change of the rate constant at the same overpotential:
    //- Jr = SI kin/kinOld, written back to SI and SE (-Jr), the base current
    //- of the assembly. Then (inverse update)
    //- eta* = G^-1(Jr), dJdEta = G'(eta*), etaE = eta* + E, j = sign(n) Jr
    //- with G(eta) = sign(n) kin (exp(alpha y) - exp(-(1 - alpha) y)),
    //- y = n F eta/(R T). G has no closed-form inverse for a general
    //- transfer coefficient alpha; it is inverted cell by cell with
    //- Newton-Raphson (butlerVolmerInverse)
    const electricSystem& electric = this->electricSys();

    if (electric.coupled())
    {
        if (ionPhase.name() != electric.ionRegionName())
        {
            FatalErrorInFunction
                << "Ion region " << ionPhase.name() << " of zone "
                << this->zoneName_ << " is not the ion region "
                << electric.ionRegionName() << " of "
                << electricSystem::dictName
                << exit(FatalError);
        }

        if (this->alpha_ <= 0 || this->alpha_ >= 1)
        {
            FatalErrorInFunction
                << "Zone " << this->zoneName_ << ": alpha = " << this->alpha_
                << ". The coupled phiE-phiI solve needs 0 < alpha < 1"
                << exit(FatalError);
        }

        if (!relaxNoted_ && this->dict_.found("relax"))
        {
            Info<< "ButlerVolmer " << this->zoneName_ << ": 'relax' is not"
                << " used by the coupled electric algorithm" << endl;
        }
        relaxNoted_ = true;

        scalarField& dJdEta =
            const_cast<volScalarField&>(electric.dJdEta());
        scalarField& etaE =
            const_cast<volScalarField&>(electric.etaE());

        //- Cells where the Newton-Raphson inversion did not converge
        label nFailed = 0;

        forAll(cells, cellI)
        {
            const label fluidId = cells[cellI];
            const label ionId =
                ionPhase.cellMap()[fluidPhase.cellMapIO()[fluidId]];

            const scalar kin = Foam::max
            (
                this->j0_.value()*coeff[fluidId]
               *Foam::pow(s[fluidId], this->gamma_),
                SMALL*this->j0_.value()
            );
            //- nf = n F/(R T): y = nf eta
            const scalar nf = n*F/Rgas/T[fluidId];

            const label electronId =
                electronPhase.cellMap()[fluidPhase.cellMapIO()[fluidId]];

            //- Rescale the carried current with the change of the rate
            //- constant, at the same overpotential: Jr was realised with
            //- kinOld. kinOld = 0 (first iteration of a run, also after a
            //- restart): ratio 1
            const scalar kinRatio =
                kinOld_[fluidId] > 0 ? kin/kinOld_[fluidId] : scalar(1);
            kinOld_[fluidId] = kin;

            const scalar Jr = SI[ionId]*kinRatio;
            SI[ionId] = Jr;
            SE[electronId] = -Jr;

            //- Solve exp(alpha y) - exp(-(1 - alpha) y) = sign(n) Jr/kin
            scalar dhdy = 1;
            label nIter = 0;
            const scalar y = butlerVolmerInverse
            (
                sign*Jr/kin,
                this->alpha_,
                dhdy,
                nIter,
                inverseTol_,
                inverseMaxIter_
            );

            if (nIter > inverseMaxIter_)
            {
                ++nFailed;
            }

            eta[fluidId] = y/nf;
            dJdEta[ionId] = kin*mag(nf)*dhdy;
            etaE[ionId] = eta[fluidId] + nernst[fluidId];

            //- Species/heat see no reverse current (as the segregated max)
            j[fluidId] = Foam::max(sign*Jr, scalar(0));

            Rj += fluidPhase.V()[fluidId] * j[fluidId];
        }

        reduce(Rj, sumOp<scalar>());
        reduce(nFailed, sumOp<label>());

        if (nFailed)
        {
            WarningInFunction
                << "Zone " << this->zoneName_ << ": the Newton-Raphson"
                << " inversion of the Butler-Volmer relation did not"
                << " converge in " << nFailed << " cells (maxIter "
                << inverseMaxIter_ << ", tolerance " << inverseTol_ << ")"
                << endl;
        }

        Info << "Total current (A) at " << this->zoneName_
             << ": " << Rj << " (coupled)" << endl;

        return;
    }

    forAll(cells, cellI)
    {
        //- get cell IDs
        label fluidId = cells[cellI];
        label electronId = electronPhase.cellMap()[fluidPhase.cellMapIO()[fluidId]];
        label ionId = ionPhase.cellMap()[fluidPhase.cellMapIO()[fluidId]];

        //- activation overpotential
        eta[fluidId] = 
        (
            eta[fluidId]*(scalar(1) - this->relax_)
          + (phiE[electronId] - phiI[ionId] - nernst[fluidId])
          * this->relax_
        );

        //- Buttler-volmer relation
        j[fluidId] = Foam::max
        (
            this->j0_.value()*
            coeff[fluidId]*
            Foam::pow(s[fluidId], this->gamma_)*
            (
                Foam::exp(n*this->alpha_*F*eta[fluidId]/Rgas/T[fluidId])
              - Foam::exp(-n*(scalar(1) - this->alpha_)*F*eta[fluidId]/Rgas/T[fluidId])
            )
            ,
            scalar(0)
        );

        //- anode side: SE = Rj, SI = -Rj
        //- cathode side: SE = -Rj, SI = Rj
        SE[electronId] = -sign*j[fluidId];
        SI[ionId] = -SE[electronId];

        Rj += fluidPhase.V()[fluidId] * j[fluidId];
    }

    reduce(Rj, sumOp<scalar>());

    Info << "Total current (A) at " << this->zoneName_
         << ": " << Rj << endl;
}

// ************************************************************************* //

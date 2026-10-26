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

#include "electricControl.H"
#include "coupledMatrix.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum<Foam::electricControl::controlType>
Foam::electricControl::controlTypeNames
({
    { controlType::galvanostatic, "galvanostatic" },
    { controlType::potentiostatic, "potentiostatic" },
});


const Foam::Enum<Foam::electricControl::galvanostaticMode>
Foam::electricControl::galvanostaticModeNames
({
    { galvanostaticMode::explicitUpdate, "explicit" },
    { galvanostaticMode::implicitUpdate, "implicit" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::electricControl::zoneMean(const scalarField& J) const
{
    const scalarField& source = J;
    const scalarField& volume = eR_.V();

    scalar iDot(0.0);

    if (cellZoneIDs_.size())
    {
        forAll(cellZoneIDs_, zoneI)
        {
            const labelList& cells = eR_.cellZones()[cellZoneIDs_[zoneI]];

            forAll(cells, i)
            {
                const label cellI = cells[i];
                iDot += source[cellI] * volume[cellI];
            }
        }

        reduce(iDot, sumOp<scalar>());
    }
    else
    {
        iDot = Foam::gSum(source*volume);
    }

    return iDot/Foam::gSum(eR_.magSf().boundaryField()[patchID_]);
}


Foam::scalar Foam::electricControl::ibar() const
{
    return zoneMean(eR_.j_);
}


void Foam::electricControl::setPatchValue(const scalar V)
{
    eR_.phi_.boundaryFieldRef()[patchID_] == V;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::electricControl::electricControl
(
    regionTypes::electric& electron,
    const dictionary& dict
)
:
    eR_(electron),
    type_(controlTypeNames.get("type", dict)),
    patchName_(dict.get<word>("patch")),
    patchID_(electron.boundaryMesh().findPatchID(patchName_)),
    zoneName_(dict.getOrDefault<word>("zone", word::null)),
    cellZoneIDs_(),
    ibar_(nullptr),
    voltage_(nullptr),
    mode_(galvanostaticMode::explicitUpdate),
    rE_(0),
    relax_(1),
    maxStep_(0.05),
    Vmin_(-GREAT),
    Vmax_(GREAT),
    gMin_(1e-8),
    psiC_()
{
    if (patchID_ == -1)
    {
        FatalIOErrorInFunction(dict)
            << "Cannot find patch " << patchName_
            << " in region " << eR_.name() << nl
            << "    Valid patches: " << eR_.boundaryMesh().names()
            << exit(FatalIOError);
    }

    if (zoneName_.size())
    {
        cellZoneIDs_ = eR_.cellZones().indices(zoneName_);

        if (cellZoneIDs_.empty())
        {
            FatalIOErrorInFunction(dict)
                << "Cannot find cell zone " << zoneName_
                << " in region " << eR_.name()
                << exit(FatalIOError);
        }
    }

    Info<< "    control: " << controlTypeNames[type_]
        << " on patch " << patchName_ << " of " << eR_.name()
        << ", current collected over ";
    if (zoneName_.size())
    {
        Info<< "zone " << zoneName_ << nl;
    }
    else
    {
        Info<< "the whole region" << nl;
    }

    if (galvanostatic())
    {
        ibar_.reset(Function1<scalar>::New("ibar", dict, &eR_));

        const dictionary& gDict = dict.subDict("galvanostatic");

        mode_ = galvanostaticModeNames.get("mode", gDict);
        rE_ = gDict.get<scalar>("rE");

        Info<< "    galvanostatic: mode " << galvanostaticModeNames[mode_]
            << ", rE = " << rE_;

        if (mode_ == galvanostaticMode::implicitUpdate)
        {
            relax_ = gDict.getOrDefault<scalar>("relax", 1.0);
            maxStep_ = gDict.getOrDefault<scalar>("maxStep", 0.05);
            Vmin_ = gDict.get<scalar>("Vmin");
            Vmax_ = gDict.get<scalar>("Vmax");
            gMin_ = gDict.getOrDefault<scalar>("gMin", 1e-8);

            if (relax_ <= 0 || relax_ > 1 || maxStep_ <= 0 || Vmin_ >= Vmax_)
            {
                FatalIOErrorInFunction(gDict)
                    << "Need 0 < relax <= 1, maxStep > 0 and Vmin < Vmax"
                    << exit(FatalIOError);
            }

            Info<< " (fallback), relax = " << relax_
                << ", maxStep = " << maxStep_
                << " V, V in [" << Vmin_ << ", " << Vmax_ << "] V"
                << ", gMin = " << gMin_;
        }

        Info<< endl;
    }
    else
    {
        voltage_.reset(Function1<scalar>::New("voltage", dict, &eR_));
        Info<< endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::electricControl::beforeSolve()
{
    volScalarField& phi = eR_.phi_;

    scalarField phiBoundary = phi.boundaryFieldRef()[patchID_];

    const scalar ibar0 = ibar();
    const scalar t = eR_.time().value();

    if (galvanostatic())
    {
        //- Implicit: the patch value is set after the coupled solve
        if (mode_ == galvanostaticMode::explicitUpdate)
        {
            phi.boundaryFieldRef()[patchID_] ==
                phiBoundary + rE_*(ibar0 - ibar_->value(t));
        }
    }
    else
    {
        phi.boundaryFieldRef()[patchID_] == voltage_->value(t);
    }

    Info << "ibar: " << ibar0 << "\t"<< "voltage: " << Foam::gAverage(phiBoundary) << endl;
}


void Foam::electricControl::afterSolve
(
    const coupledMatrix& matrix,
    scalarField& psi,
    const scalarField& etaE,
    const dictionary& solverControls
)
{
    //- Explicit update and potentiostatic: the patch value was set before
    //- the solve
    if (!needsCoupled())
    {
        return;
    }

    //- Patch value of the solve (psi = phi1)
    const scalar Vold = Foam::gAverage(eR_.phi_.boundaryField()[patchID_]);

    //- Unit response M phic = c, warm start from the last iteration
    const scalarField c(matrix.unitSource(patchID_));

    if (psiC_.size() != psi.size())
    {
        psiC_.setSize(psi.size(), Zero);
    }

    solverPerformance perf =
        matrix.solve(psiC_, c, solverControls, "phiUnit");

    Info<< perf.solverName() << ":  Solving for phiUnit"
        << ", Initial residual = " << perf.initialResidual()
        << ", Final residual = " << perf.finalResidual()
        << ", No Iterations " << perf.nIterations() << endl;

    //- Current of phi1 and its change per volt, in the electron region:
    //- J_E = -(J* + s (phiE - phiI - eta* - E)), s = r/V
    const regionTypes::electric& iR = matrix.ion();
    const scalarField& Vi = iR.V();
    const scalarField& Jstar = iR.j_;
    const labelList& pairE = matrix.pairE();
    const labelList& pairI = matrix.pairI();
    const scalarField& r = matrix.r();
    const label nE = matrix.nElectron();

    scalarField J1(eR_.j_);
    scalarField dJ(eR_.nCells(), Zero);

    forAll(pairI, p)
    {
        const label e = pairE[p];
        const label i = pairI[p];
        const scalar sp = r[p]/Vi[i];

        J1[e] = -(Jstar[i] + sp*(psi[e] - psi[nE + i] - etaE[i]));
        dJ[e] = -sp*(psiC_[e] - psiC_[nE + i]);
    }

    const scalar ibar1 = zoneMean(J1);
    const scalar g = zoneMean(dJ);
    const scalar target = ibar_->value(eR_.time().value());

    //- No usable sensitivity: keep phi1, explicit update of the patch
    if (g > -gMin_)
    {
        const scalar V = Vold + rE_*(ibar1 - target);

        setPatchValue(V);

        WarningInFunction
            << "d ibar/dV = " << g << " >= -gMin = " << -gMin_
            << ", explicit update with rE = " << rE_ << endl;

        Info<< "Implicit control: ibar1 = " << ibar1
            << ", g = " << g
            << ", fallback V = " << V << endl;

        return;
    }

    //- Newton step on the patch value, limited
    scalar dV = relax_*(target - ibar1)/g;
    bool limited = false;

    if (mag(dV) > maxStep_)
    {
        dV = sign(dV)*maxStep_;
        limited = true;
    }

    scalar V = Vold + dV;

    if (V < Vmin_ || V > Vmax_)
    {
        V = min(max(V, Vmin_), Vmax_);
        limited = true;
    }

    dV = V - Vold;

    //- Superposition phi = phi1 + dV phic, new patch value
    psi += dV*psiC_;

    setPatchValue(V);

    Info<< "Implicit control: ibar1 = " << ibar1
        << ", g = " << g << " A/m2/V"
        << ", -1/g = " << -1/g << " Ohm m2"
        << ", dV = " << dV
        << ", V = " << V
        << ", ibar = " << ibar1 + g*dV
        << (limited ? " (limited)" : "") << endl;
}


// ************************************************************************* //

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
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::electricControl::ibar() const
{
    const scalarField& source = eR_.j_;
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
    rE_(0)
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
            << ", rE = " << rE_ << endl;
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
        phi.boundaryFieldRef()[patchID_] ==
            phiBoundary + rE_*(ibar0 - ibar_->value(t));
    }
    else
    {
        phi.boundaryFieldRef()[patchID_] == voltage_->value(t);
    }

    Info << "ibar: " << ibar0 << "\t"<< "voltage: " << Foam::gAverage(phiBoundary) << endl;
}


void Foam::electricControl::afterSolve(coupledMatrix& matrix)
{
    //- Explicit update: the patch value was set before the solve
}


// ************************************************************************* //

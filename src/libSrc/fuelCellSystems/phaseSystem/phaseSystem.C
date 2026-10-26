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

#include "fvCFD.H"
#include "phaseSystem.H"
#include "surfaceInterpolate.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

#include "pimpleControl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseSystem, 0);
    defineRunTimeSelectionTable(phaseSystem, dictionary);
}

const Foam::word Foam::phaseSystem::propertiesName("regionProperties");


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::IOobject Foam::phaseSystem::createIOobject
(
    const fvMesh& mesh
) const
{
    IOobject io
    (
        propertiesName,
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (io.typeHeaderOk<IOdictionary>(true))
    {
        Info<< "Creating phaseSystem from " << io.name() << nl << endl;

        io.readOpt(IOobject::MUST_READ_IF_MODIFIED);
    }
    else
    {
        Info<< "No phaseSystem present" << nl << endl;
        Info<< "Single phase solver assumed" << nl << endl;

        io.readOpt(IOobject::NO_READ);
    }

    return io;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseSystem::phaseSystem
(
    const fvMesh& mesh
)
:
    IOdictionary(createIOobject(mesh)),

    mesh_(mesh),

    porousZone_
    (
        new porousZoneList
        (
            mesh
        )
    ),

    phaseModels_(lookup("phases"), phaseModel::iNew(*this)),

    MRF_(mesh_),

    continuous_(lookupOrDefault<word>("continuous", phaseModels_[0].name())),

    p_rgh_
    (
        IOobject
        (
            "p_rgh",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    pimple_(const_cast<fvMesh&>(mesh))
{
    mesh.setFluxRequired(p_rgh_.name());

    // Create LTS field, if LTS is enabled
    bool LTS = fv::localEulerDdt::enabled(mesh);

    if (LTS)
    {
        Info<< "Using LTS" << endl;

        trDeltaT_ = autoPtr<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    fv::localEulerDdt::rDeltaTName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("one", dimless/dimTime, 1),
                extrapolatedCalculatedFvPatchScalarField::typeName
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseSystem::~phaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::phaseSystem::correct()
{
    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].correct();
    }
}


void Foam::phaseSystem::correctThermo()
{
    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].correctThermo();
    }
}


void Foam::phaseSystem::correctTurbulence()
{
    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].correctTurbulence();
    }
}


void Foam::phaseSystem::correctEnergyTransport()
{
    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].correctEnergyTransport();
    }
}


void Foam::phaseSystem::setRDeltaT()
{
    if (!trDeltaT_.valid())
    {
        return;
    }

    volScalarField& rDeltaT = trDeltaT_();

    const dictionary& LTSdict = this->subDict("LTS");

    scalar maxCo
    (
        LTSdict.lookupOrDefault<scalar>("maxCo", 0.2)
    );

    scalar maxDeltaT
    (
        LTSdict.lookupOrDefault<scalar>("maxDeltaT", GREAT)
    );

    scalar rDeltaTSmoothingCoeff
    (
        LTSdict.lookupOrDefault<scalar>("rDeltaTSmoothingCoeff", 0.02)
    );

    // Set the reciprocal time-step from the local Courant number
    rDeltaT.ref() = max
    (
        1/dimensionedScalar("maxDeltaT", dimTime, maxDeltaT),
        fvc::surfaceSum(phiMagMax())()()
       /((2*maxCo)*mesh_.V())
    );

    // Update tho boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();

    fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);

    Info<< "Flow time scale min/max = "
        << gMin(1/rDeltaT.primitiveField())
        << ", " << gMax(1/rDeltaT.primitiveField()) << endl;
}


// ************************************************************************* //

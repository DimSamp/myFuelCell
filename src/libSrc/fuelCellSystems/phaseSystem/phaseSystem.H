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

Class
    Foam::phaseSystem

Description
    Base class of the (single-phase) flow system of a fluid region.

SourceFiles
    phaseSystem.C
    newPhaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef phaseSystem_H
#define phaseSystem_H

#include "IOdictionary.H"

#include "phaseModel.H"
#include "PtrListDictionary.H"

#include "IOMRFZoneList.H"
#include "fvOptions.H"
#include "pimpleControl.H"

#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatricesFwd.H"

#include "porousZoneList.H"
#include "localEulerDdtScheme.H"

#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class pimpleControl;

/*---------------------------------------------------------------------------*\
                         Class phaseSystem Declaration
\*---------------------------------------------------------------------------*/

class phaseSystem
:
    public IOdictionary
{
public:

    // Public typedefs

        typedef PtrListDictionary<phaseModel> phaseModelList;


protected:

    // Protected data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Porous zone list
        autoPtr<porousZoneList> porousZone_;

        //- Phase models
        phaseModelList phaseModels_;

        //- Optional MRF zones
        IOMRFZoneList MRF_;

        //- Added by Shidong, 27.07
        //- continuous phase
        word continuous_;

        //- Pressure p_rgh
        volScalarField p_rgh_;

        //- local time step
        autoPtr<volScalarField> trDeltaT_;

        //- Pimple control
        pimpleControl pimple_;

        //- Ended add


    // Protected member functions

        //- Create IO dictionary
        IOobject createIOobject(const fvMesh& mesh) const;


public:

    //- Runtime type information
    TypeName("phaseSystem");

    //- Default name of the phase properties dictionary
    static const word propertiesName;

    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            phaseSystem,
            dictionary,
            (
                const fvMesh& mesh
            ),
            (mesh)
        );

    // Constructors

        //- Construct from fvMesh
        phaseSystem(const fvMesh& mesh);


    //- Destructor
    virtual ~phaseSystem();



    //- selector
        static autoPtr<phaseSystem> New
        (
            const fvMesh& mesh
        );


    // Member Functions

        // Access

            //- Return the mesh
            inline const fvMesh& mesh() const;

            //- Return the phase models
            inline const phaseModelList& phases() const;

            //- Access the phase models
            inline phaseModelList& phases();

            //- Return MRF zones
            inline const IOMRFZoneList& MRF() const;

            //- Access the fvOptions
            inline fv::options& fvOptions() const;

            //- Return porousZones
            inline const porousZoneList& porousZone() const;

            inline porousZoneList& porousZone();

            //- Name of continuous phase
            inline const word& continuous() const;

            inline word& continuous();

            //- Pressure
            inline const volScalarField& p_rgh() const;

            inline volScalarField& p_rgh();


        // Properties

            //- Return relative flux magnitude
            virtual tmp<surfaceScalarField> phirMag() const = 0;

            //- Return max magnitude flux
            virtual tmp<surfaceScalarField> phiMagMax() const = 0;


        // Evolution

            //- Solve for the phase fractions
            virtual void solve() = 0;

            //- Correct the fluid properties other than those listed below
            virtual void correct();

            //- Correct the thermodynamics
            virtual void correctThermo();

            //- Correct the turbulence
            virtual void correctTurbulence();

            //- Correct the energy transport e.g. alphat
            virtual void correctEnergyTransport();

            //- Set delta T
            void setRDeltaT();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "phaseSystemI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

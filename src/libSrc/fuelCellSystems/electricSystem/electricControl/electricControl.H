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
    Foam::electricControl

Description
    Operating control of the cell: the potential of a Dirichlet patch of the
    electron region (e.g. interconnect0) is either prescribed
    (potentiostatic) or adjusted such that the mean current density through
    that patch reaches the target (galvanostatic).

    The current is collected from J of the electron region, over a cell zone
    (e.g. ccl) or over the whole region.

    Galvanostatic modes:
        explicit:  V += rE (ibar - ibar*) before the solve of each iteration

    \verbatim
    control
    {
        type    galvanostatic;      // galvanostatic | potentiostatic
        patch   interconnect0;
        zone    ccl;                // optional, default: whole region
        ibar    constant 8000;      // potentiostatic: voltage <Function1>

        galvanostatic
        {
            mode    explicit;
            rE      1e-5;
        }
    }
    \endverbatim

SourceFiles
    electricControl.C

\*---------------------------------------------------------------------------*/

#ifndef electricControl_H
#define electricControl_H

#include "electric.H"
#include "Function1.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class coupledMatrix;

/*---------------------------------------------------------------------------*\
                       Class electricControl Declaration
\*---------------------------------------------------------------------------*/

class electricControl
{
public:

    // Public data types

        //- Operating mode
        enum class controlType
        {
            galvanostatic,
            potentiostatic
        };

        static const Enum<controlType> controlTypeNames;

        //- Galvanostatic update
        enum class galvanostaticMode
        {
            explicitUpdate
        };

        static const Enum<galvanostaticMode> galvanostaticModeNames;


private:

    // Private data

        //- Electron region, holds the controlled patch
        regionTypes::electric& eR_;

        //- Operating mode
        controlType type_;

        //- Controlled patch
        word patchName_;
        label patchID_;

        //- Zone over which the current is collected (optional)
        word zoneName_;
        labelList cellZoneIDs_;

        //- Target mean current density (galvanostatic)
        autoPtr<Function1<scalar>> ibar_;

        //- Prescribed voltage (potentiostatic)
        autoPtr<Function1<scalar>> voltage_;

        //- Galvanostatic update
        galvanostaticMode mode_;

        //- Gain of the explicit update [V m2/A]
        scalar rE_;


    // Private Member Functions

        //- Mean current density through the controlled patch
        scalar ibar() const;

        //- No copy construct
        electricControl(const electricControl&) = delete;

        //- No copy assignment
        void operator=(const electricControl&) = delete;


public:

    // Constructors

        //- Construct from the control dictionary
        electricControl
        (
            regionTypes::electric& electron,
            const dictionary& dict
        );


    //- Destructor
    ~electricControl() = default;


    // Member Functions

        //- Galvanostatic?
        bool galvanostatic() const
        {
            return type_ == controlType::galvanostatic;
        }

        //- Does the mode need the coupled algorithm?
        bool needsCoupled() const
        {
            return false;
        }

        //- Update the patch value before the potential solve
        void beforeSolve();

        //- Correct the potentials after the coupled solve
        //  (nothing to do for the explicit update)
        void afterSolve(coupledMatrix& matrix);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

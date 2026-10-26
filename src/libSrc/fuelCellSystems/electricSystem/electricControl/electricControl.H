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
        implicit:  (coupled algorithm only) the coupled system is linear in
                   the patch value V, M phi = b0 + V c. After the normal
                   solve (phi1, with V_old) a second solve M phic = c gives
                   the response to +1 V, so phi(V) = phi1 + (V - V_old) phic
                   and the mean current ibar(V) = ibar1 + g (V - V_old),
                   g = d ibar/dV. The update
                       dV = relax (ibar* - ibar1)/g,
                   limited to +-maxStep and to [Vmin, Vmax], is applied to
                   phi and to the patch value within the same iteration.
                   If g >= -gMin the update is skipped and the explicit
                   update with rE is used instead.

    \verbatim
    control
    {
        type    galvanostatic;      // galvanostatic | potentiostatic
        patch   interconnect0;
        zone    ccl;                // optional, default: whole region
        ibar    constant 8000;      // potentiostatic: voltage <Function1>

        galvanostatic
        {
            mode    explicit;       // explicit | implicit
            rE      1e-5;           // explicit; implicit: fallback

            // implicit only
            relax   1;              // (0, 1]
            maxStep 0.05;           // [V]
            Vmin    0;              // [V]
            Vmax    1.2;            // [V]
            gMin    1e-8;           // optional [A/m2/V]
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
            explicitUpdate,
            implicitUpdate
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

        //- Implicit: relaxation of the voltage step, (0, 1]
        scalar relax_;

        //- Implicit: largest voltage step per iteration [V]
        scalar maxStep_;

        //- Implicit: bounds of the patch value [V]
        scalar Vmin_;
        scalar Vmax_;

        //- Implicit: smallest accepted -d ibar/dV [A/m2/V]
        scalar gMin_;

        //- Implicit: unit response phic of the last iteration
        //  (initial guess of the next unit solve)
        scalarField psiC_;


    // Private Member Functions

        //- Mean current density through the controlled patch for the
        //  current J of the electron region
        scalar zoneMean(const scalarField& J) const;

        //- Mean current density through the controlled patch
        scalar ibar() const;

        //- Set a uniform patch value
        void setPatchValue(const scalar V);

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
            return
                galvanostatic()
             && mode_ == galvanostaticMode::implicitUpdate;
        }

        //- Update the patch value before the potential solve
        void beforeSolve();

        //- Correct the solution psi = [phiE; phiI] of the coupled solve
        //  and the patch value (implicit mode), before psi is scattered
        //  and the current is realised
        void afterSolve
        (
            const coupledMatrix& matrix,
            scalarField& psi,
            const scalarField& etaE,
            const dictionary& solverControls
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

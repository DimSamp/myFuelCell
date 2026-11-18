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
    Foam::ButlerVolmer.C

Description
    The Butler-Volmer relation is an abstract concept in electrochemistry
    that describes the mathematical relationship between the rate of
    an electrochemical reaction and the overpotential
    (the deviation from the equilibrium potential) experienced
    by an electrode during that reaction.

    \f[
        j = j_0 \prod_i\left(\frac{p_i}{p_{ref}}\right)^\ksi (1-s)^\gamma
        \left(
            \exp{\frac{n\alpha F}{RT}\eta}
          - \exp{-\frac{n(1 - \alpha)F}{RT}\eta}
        \right)
    \f]

    The transfer coefficient alpha (default 0.5) multiplies the dominant
    branch on both electrodes, because n < 0 and eta < 0 on the cathode.

    With the coupled electric algorithm (electricSystem) the relation is
    inverted in every catalyst cell, eta = G^-1(J), for 0 < alpha < 1 by
    Newton-Raphson (butlerVolmerInverse.H). Optional entries:

    \verbatim
        inverseTolerance    1e-12;  // relative step of y = n F eta/(R T)
        inverseMaxIter      100;
    \endverbatim

SourceFiles
    ButlerVolmer.C
    butlerVolmerInverse.H

\*---------------------------------------------------------------------------*/

#ifndef ButlerVolmer_H
#define ButlerVolmer_H

#include "ActivationOverpotentialModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phaseModel;

namespace activationOverpotentialModels
{

/*---------------------------------------------------------------------------*\
                           Class ButlerVolmer Declaration
\*---------------------------------------------------------------------------*/
template<class Thermo>
class ButlerVolmer
:
    public ActivationOverpotentialModel<Thermo>
{
    // Private data

        //- The note that 'relax' is unused has been printed (coupled)
        bool relaxNoted_;

        //- Newton-Raphson inversion (coupled): relative step tolerance
        scalar inverseTol_;

        //- Newton-Raphson inversion (coupled): maximum iterations per cell
        label inverseMaxIter_;

        //- Rate constant kin of the last coupled iteration, per fluid cell
        //- (coupled; catalyst zone cells). Zero where it was never set:
        //- the rate change of such a cell is taken as 1. Read on restart
        volScalarField kinOld_;


public:
    //- Runtime type information
    TypeName("ButlerVolmer");

    // Constructors

    ButlerVolmer
    (
        const phaseModel& phase,
        const dictionary& dict
    );

    //- Destructor
        virtual ~ButlerVolmer();


    // Member Functions

        virtual void correct();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace ActivationOverpotentialModels
} // End namespace Foam
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ButlerVolmer.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

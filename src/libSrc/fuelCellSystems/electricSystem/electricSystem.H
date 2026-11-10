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
    Foam::electricSystem

Description
    Solves the potentials of one electron and one ion region and applies the
    operating control (constant/electricProperties).

    Algorithms:
        segregated: phiE and phiI are solved one after the other with the
                    current J of the activation overpotential model. phiI
                    has no Dirichlet patch; its level is set by a reference
                    that moves with the net current (ionReferenceRelax).
        coupled:    one matrix for [phiE; phiI] with the Butler-Volmer
                    kinetics linearised in the catalyst cells (coupledMatrix).
                    Requires the ButlerVolmer model, 0 < alpha < 1.

    \verbatim
    electronRegion  phiE;
    ionRegion       phiI;

    algorithm       coupled;    // segregated | coupled

    segregatedCoeffs
    {
        ionReferenceRelax  -4e-11;
    }

    coupledCoeffs
    {
        omega   1;              // Newton damping, (0, 1]
    }

    control
    {
        ...                     // see electricControl
    }
    \endverbatim

    The coupled linear solver is read from system/<ionRegion>/fvSolution,
    entry phiCoupled.

SourceFiles
    electricSystem.C

\*---------------------------------------------------------------------------*/

#ifndef electricSystem_H
#define electricSystem_H

#include "IOdictionary.H"
#include "electric.H"
#include "electricControl.H"
#include "coupledMatrix.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class regionTypeList;

/*---------------------------------------------------------------------------*\
                       Class electricSystem Declaration
\*---------------------------------------------------------------------------*/

class electricSystem
:
    public IOdictionary
{
public:

    // Public data types

        //- Solution algorithm
        enum class algorithmType
        {
            segregated,
            coupled
        };

        static const Enum<algorithmType> algorithmTypeNames;


private:

    // Private data

        //- Electron and ion regions
        regionTypes::electric& eR_;
        regionTypes::electric& iR_;

        //- Solution algorithm
        algorithmType algorithm_;

        //- Segregated: relaxation of the ion potential reference
        scalar ionReferenceRelax_;

        //- Coupled: Newton damping of the coupling coefficient
        scalar omega_;

        //- Operating control
        autoPtr<electricControl> control_;

        //- Coupled: combined matrix
        autoPtr<coupledMatrix> matrix_;

        //- Coupled: written by ButlerVolmer in the ion catalyst cells,
        //  dJ/deta at the base point [A/V/m3] and eta* + E [V]
        autoPtr<volScalarField> dJdEta_;
        autoPtr<volScalarField> etaE_;


    // Private Member Functions

        //- Find an electric region by name
        static regionTypes::electric& lookupRegion
        (
            regionTypeList& regions,
            const word& name
        );

        //- Segregated solve of one region
        void solveRegion(regionTypes::electric& r, const scalar refRelax);

        //- Segregated algorithm
        void solveSegregated();

        //- Coupled algorithm
        void solveCoupled();

        //- No copy construct
        electricSystem(const electricSystem&) = delete;

        //- No copy assignment
        void operator=(const electricSystem&) = delete;


public:

    //- Runtime type information
    TypeName("electricSystem");

    //- Name of the dictionary (constant/electricProperties)
    static const word dictName;


    // Constructors

        //- Construct from the (already created) regions
        electricSystem(const Time& runTime, regionTypeList& regions);


    //- Destructor
    virtual ~electricSystem() = default;


    // Member Functions

        //- Coupled algorithm?
        bool coupled() const
        {
            return algorithm_ == algorithmType::coupled;
        }

        //- Name of the ion region
        const word& ionRegionName() const
        {
            return iR_.name();
        }

        //- Coupled: dJ/deta at the base point (ion region)
        const volScalarField& dJdEta() const
        {
            return dJdEta_();
        }

        //- Coupled: eta* + E (ion region)
        const volScalarField& etaE() const
        {
            return etaE_();
        }

        //- Operating control, before the potential solve
        void correct();

        //- Solve the potentials
        void solve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

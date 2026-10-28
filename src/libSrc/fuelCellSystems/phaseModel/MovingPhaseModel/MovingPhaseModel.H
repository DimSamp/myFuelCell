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
    Foam::MovingPhaseModel

Description
    Class which represents a moving fluid phase. Holds the velocity, fluxes and
    turbulence model and can generate the momentum equation. The interface is
    quite restrictive as it also has to support an equivalent stationary model,
    which does not store motion fields or a turbulence model.

    Possible future extensions include separating the turbulent fuctionality
    into another layer.

See also
    StationaryPhaseModel

SourceFiles
    MovingPhaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef MovingPhaseModel_H
#define MovingPhaseModel_H

#include "phaseModel.H"
#include "phaseCompressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class MovingPhaseModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseModel>
class MovingPhaseModel
:
    public BasePhaseModel
{
protected:

    // Protected data

        //- Velocity field
        volVectorField U_;

        //- Flux
        surfaceScalarField phi_;

        //- Mass flux
        surfaceScalarField alphaRhoPhi_;

        //- Turbulence model
        autoPtr<phaseCompressibleTurbulenceModel> turbulence_;

        //- Continuity error due to the flow
        volScalarField continuityErrorFlow_;

        //- Continuity error due to any sources
        volScalarField continuityErrorSources_;

private:

    // Private static member functions

        //- Calculate and return the flux field
        tmp<surfaceScalarField> phi(const volVectorField& U) const;


public:

    // Constructors

        MovingPhaseModel
        (
            const phaseSystem& fluid,
            const word& phaseName,
            const label index
        );


    //- Destructor
    virtual ~MovingPhaseModel();


    // Member Functions

        //- Correct the phase properties other than the thermo and turbulence
        virtual void correct();

        //- Correct the turbulence
        virtual void correctTurbulence();

        //- Correct the energy transport e.g. alphat
        virtual void correctEnergyTransport();


        // Momentum

            //- Return the momentum equation
            virtual tmp<fvVectorMatrix> UEqn();

            //- Return the velocity
            virtual tmp<volVectorField> U() const;

            //- Access the velocity
            virtual volVectorField& URef();

            //- Return the volumetric flux
            virtual tmp<surfaceScalarField> phi() const;

            //- Access the volumetric flux
            virtual surfaceScalarField& phiRef();

            //- Return the mass flux of the phase
            virtual tmp<surfaceScalarField> alphaRhoPhi() const;

            //- Access the mass flux of the phase
            virtual surfaceScalarField& alphaRhoPhiRef();

            //- Return the continuity error
            virtual tmp<volScalarField> continuityError() const;


        // Turbulence

            //- Return the effective dynamic viscosity
            virtual tmp<volScalarField> muEff() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MovingPhaseModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

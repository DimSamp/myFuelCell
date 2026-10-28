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
    Foam::phaseModel

Description
    Base class for diversed phases

SourceFiles
    phaseModel.C
    newPhaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef phaseModel_H
#define phaseModel_H

#include "dictionary.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatricesFwd.H"
#include "rhoThermo.H"
#include "phaseCompressibleTurbulenceModelFwd.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phaseSystem;

/*---------------------------------------------------------------------------*\
                           Class phaseModel Declaration
\*---------------------------------------------------------------------------*/

class phaseModel
:
    public volScalarField
{
    // Private data

        //- Reference to the phaseSystem to which this phase belongs
        const phaseSystem& fluid_;

        //- Name of phase
        word name_;

        //- Index of phase
        label index_;

public:

    //- Runtime type information
    ClassName("phaseModel");

    //- Default name of the water name
    static const word water;

    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            phaseModel,
            dictionary,
            (
                const phaseSystem& fluid,
                const word& phaseName,
                const label index
            ),
            (fluid, phaseName, index)
        );


    // Constructors

        phaseModel
        (
            const phaseSystem& fluid,
            const word& phaseName,
            const label index
        );

        //- Return clone
        autoPtr<phaseModel> clone() const;


    // Selectors

        static autoPtr<phaseModel> New
        (
            const phaseSystem& fluid,
            const word& phaseName,
            const label index
        );

        //- Return a pointer to a new phase created on freestore
        //  from Istream
        class iNew
        {
            const phaseSystem& fluid_;
            mutable label indexCounter_;

        public:

            iNew
            (
                const phaseSystem& fluid
            )
            :
                fluid_(fluid),
                indexCounter_(-1)
            {}

            autoPtr<phaseModel> operator()(Istream& is) const
            {
                indexCounter_++;
                return autoPtr<phaseModel>
                (
                    phaseModel::New(fluid_, word(is), indexCounter_)
                );
            }
        };


    //- Destructor
    virtual ~phaseModel();


    // Member Functions

        //- Return the name of this phase
        const word& name() const;

        //- Return the name of the phase for use as the keyword in PtrDictionary
        const word& keyword() const;

        //- Return the index of the phase
        label index() const;

        //- Return the system to which this phase belongs
        const phaseSystem& fluid() const;

        //- Correct the phase properties
        virtual void correct();

        //- Correct the thermodynamics
        virtual void correctThermo();

        //- Correct the turbulence
        virtual void correctTurbulence();

        //- Correct the energy transport
        virtual void correctEnergyTransport();

        //- Return the total mass transfer rate
        virtual tmp<volScalarField> dmdt() const = 0;


        // Thermo

            //- Return whether the phase is isothermal
            virtual bool isothermal() const = 0;

            //- Return the heat source/sink (as explicit form)
            virtual tmp<volScalarField> Qdot() const = 0;

            //- Return the thermophysical model
            virtual const rhoThermo& thermo() const = 0;

            //- Access the thermophysical model
            virtual rhoThermo& thermoRef() = 0;

            //- Return the density field
            virtual tmp<volScalarField> rho() const = 0;


        // Species

            //- Return whether the phase is pure (i.e., not multi-component)
            virtual bool pure() const = 0;

            //- Return the species fraction equation
            virtual tmp<fvScalarMatrix> YiEqn(volScalarField& Yi) = 0;

            //- Return the source/sink term for species by name
            virtual const volScalarField& iDmdt(const word& name) const = 0;

            //- Return the species mass fractions
            virtual const PtrList<volScalarField>& Y() const = 0;

            //- Return a species mass fraction by name
            virtual const volScalarField& Y(const word& name) const = 0;

            //- Access the species mass fractions
            virtual PtrList<volScalarField>& YRef() = 0;

            //- Return the active species mass fractions
            virtual const UPtrList<volScalarField>& YActive() const = 0;

            //- Access the active species mass fractions
            virtual UPtrList<volScalarField>& YActiveRef() = 0;

            //- Return a species mass fraction by name
            virtual const volScalarField& X(const word& name) const = 0;

        // Momentum

            //- Return the momentum equation
            virtual tmp<fvVectorMatrix> UEqn() = 0;

            //- Return the velocity
            virtual tmp<volVectorField> U() const = 0;

            //- Access the velocity
            virtual volVectorField& URef() = 0;

            //- Return the volumetric flux
            virtual tmp<surfaceScalarField> phi() const = 0;

            //- Access the volumetric flux
            virtual surfaceScalarField& phiRef() = 0;

            //- Return the mass flux of the phase
            virtual tmp<surfaceScalarField> alphaRhoPhi() const = 0;

            //- Access the mass flux of the phase
            virtual surfaceScalarField& alphaRhoPhiRef() = 0;

            //- Return the continuity error
            virtual tmp<volScalarField> continuityError() const = 0;

        // Transport

            //- Return the laminar dynamic viscosity
            virtual tmp<volScalarField> mu() const = 0;

            //- Return the laminar dynamic viscosity on a patch
            virtual tmp<scalarField> mu(const label patchi) const = 0;

            //- Return the laminar kinematic viscosity
            virtual tmp<volScalarField> nu() const = 0;

            //- Return the laminar kinematic viscosity on a patch
            virtual tmp<scalarField> nu(const label patchi) const = 0;

            //- Return the laminar thermal conductivity
            virtual tmp<volScalarField> kappa() const = 0;

            //- Return the laminar thermal conductivity on a patch
            virtual tmp<scalarField> kappa(const label patchi) const = 0;

            //- Thermal diffusivity for energy of mixture [kg/m/s]
            virtual tmp<volScalarField> alphahe() const = 0;

            //- Thermal diffusivity for energy of mixture for patch [kg/m/s]
            virtual tmp<scalarField> alphahe(const label patchi) const = 0;

            //- Return the effective thermal conductivity, given the turbulent
            //  thermal diffusivity
            virtual tmp<volScalarField> kappaEff
            (
                const volScalarField& alphat
            ) const = 0;

            //- Return the effective thermal conductivity on a patch, given the
            //  turbulent thermal diffusivity
            virtual tmp<scalarField> kappaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const = 0;

            //- Return the laminar thermal diffusivity for enthalpy
            virtual tmp<volScalarField> alpha() const = 0;

            //- Return the laminar thermal diffusivity for enthalpy on a patch
            virtual tmp<scalarField> alpha(const label patchi) const = 0;

            //- Return the effective thermal diffusivity for enthalpy, given the
            //  turbulent thermal diffusivity
            virtual tmp<volScalarField> alphaEff
            (
                const volScalarField& alphat
            ) const = 0;

            //- Return the effective thermal diffusivity for enthalpy on a
            //  patch, given the turbulent thermal diffusivity
            virtual tmp<scalarField> alphaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const = 0;


        // Turbulence

            //- Return the effective dynamic viscosity
            virtual tmp<volScalarField> muEff() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

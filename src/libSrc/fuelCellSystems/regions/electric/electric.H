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
    Foam::electric.C

Description
    The 'electric' region considers electron and ion/proton transfer in electrochemical devices.
    Temperature and thermal properties will be mapped from and mapped to the main domain.
    The region holds the potential, current and conductivity fields and the
    dissolved-water model. The potential equations and the operating control
    (galvanostatic/potentiostatic) are solved by electricSystem
    (constant/electricProperties).

SourceFiles
    electric.C

\*---------------------------------------------------------------------------*/

#ifndef electric_H
#define electric_H

#include "dictionary.H"
#include "regionType.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "wordList.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fuelCellSystem;
class sigmaModelList;
class dissolvedModel;
class electricSystem;
class coupledMatrix;
class electricControl;

namespace regionTypes
{

/*---------------------------------------------------------------------------*\
                           Class electric Declaration
\*---------------------------------------------------------------------------*/

class electric
:
    public regionType
{
private:
        //-mesh
        const fvMesh& mesh_;

        //- Current density (A/m2)
        volVectorField i_;

        //- Potential (V)
        volScalarField phi_;

        //- Source term (A/m3)
        volScalarField j_;

        //- Electric conductivity
        volScalarField sigmaField_;

        //- Temperature
        volScalarField T_;

        //- Sigma model (electric conductivity)
        autoPtr<sigmaModelList> sigma_;

        //- Dissolved water
        //- In PEM fuel cell/electrolyzer, water can dissolve from one side
        //- to the other. This process is usually affected by the proton transfer.
        //- Refer to the region of CCM (carbon coated membrane)
        Switch dissolveOnOff_;

        autoPtr<dissolvedModel> dissolved_;

        //- Stop with a hint if a key moved to electricProperties is found
        void checkMovedKeys() const;

        //- The potential equations are assembled and solved here
        friend class Foam::electricSystem;
        friend class Foam::coupledMatrix;
        friend class Foam::electricControl;

public:
        TypeName("electric");

    // Constructors

        electric
        (
            const fvMesh& mesh,
            const word& regionName
        );


    //- Destructor
        virtual ~electric();


    // Member Functions

        //- Correct fields
        virtual void correct();

        //- Set reversed delta t
        virtual void setRDeltaT();

        //- Solve PDEs (done by electricSystem)
        virtual void solve();

        //- Map heat source/sink to the parent mesh
        virtual void mapToCell(fuelCellSystem&);

        //- Map temperature from the parent mesh
        virtual void mapFromCell(fuelCellSystem&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace regionTypes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

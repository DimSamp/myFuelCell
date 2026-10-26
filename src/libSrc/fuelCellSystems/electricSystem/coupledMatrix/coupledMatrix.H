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
    Foam::coupledMatrix

Description
    One matrix for the unknowns [phiE; phiI] of an electron and an ion
    region, with the Butler-Volmer kinetics linearised (inverse Newton
    update) in the paired catalyst cells:

        row e: (aP + r) phiE - r phiI = -V J* + r (eta* + E)
        row i: (aP + r) phiI - r phiE = +V J* - r (eta* + E)

    with r = omega V dJ/deta(eta*).

    The addressing (cell pairs, face order, processor interfaces) is built
    once. assemble() fills the coefficients and the source b0 of the current
    iteration; solve() may be called several times on the same matrix with
    different right-hand sides (e.g. a unit response for implicit
    galvanostatic control); realise() sets the current J from given
    potentials.

    In parallel all copies of a cell are on the same processor
    (decomposeParID), so the cell pairs are always local.

SourceFiles
    coupledMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef coupledMatrix_H
#define coupledMatrix_H

#include "lduPrimitiveMesh.H"
#include "lduMatrix.H"
#include "FieldField.H"
#include "electric.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class coupledMatrix Declaration
\*---------------------------------------------------------------------------*/

class coupledMatrix
{
    // Private data

        //- Electron and ion regions
        regionTypes::electric& eR_;
        regionTypes::electric& iR_;

        //- Newton damping of the coupling coefficient, (0, 1]
        const scalar omega_;

        //- Number of cells: electron, ion, both
        const label nE_;
        const label nI_;
        const label nCells_;

        //- Number of internal faces: electron, ion
        const label nFE_;
        const label nFI_;

        //- Cell pairs: ion cell pairI_[p] <-> electron cell pairE_[p]
        labelList pairE_;
        labelList pairI_;

        //- Combined face index -> upper-triangular face index
        labelList oldToNew_;

        //- Coupled (processor) patches, in interface order:
        //  (0 = electron / 1 = ion region, patch index)
        List<labelPair> procPatches_;

        //- Interface fields of the combined matrix (the phi patch fields)
        lduInterfaceFieldPtrsList interfaces_;

        //- Combined mesh, owns the processor interfaces
        autoPtr<lduPrimitiveMesh> lduMesh_;

        //- Combined matrix (symmetric: diag and upper)
        autoPtr<lduMatrix> matrix_;

        //- Interface coefficients of the current assembly
        FieldField<Field, scalar> interfaceBouCoeffs_;
        FieldField<Field, scalar> interfaceIntCoeffs_;

        //- Source b0 of the current assembly
        scalarField source_;

        //- Coupling coefficients r = omega V dJ/deta, per pair
        scalarField r_;


    // Private Member Functions

        //- Build the combined addressing
        void buildAddressing();

        //- No copy construct
        coupledMatrix(const coupledMatrix&) = delete;

        //- No copy assignment
        void operator=(const coupledMatrix&) = delete;


public:

    // Constructors

        //- Construct from the electron and ion regions
        coupledMatrix
        (
            regionTypes::electric& electron,
            regionTypes::electric& ion,
            const scalar omega
        );


    //- Destructor
    ~coupledMatrix() = default;


    // Member Functions

        // Access

            //- Number of unknowns
            label size() const
            {
                return nCells_;
            }

            //- Ion cells of the pairs
            const labelList& pairI() const
            {
                return pairI_;
            }

            //- Electron cells of the pairs
            const labelList& pairE() const
            {
                return pairE_;
            }

            //- Coupling coefficients of the current assembly
            const scalarField& r() const
            {
                return r_;
            }

            //- Source b0 of the current assembly
            const scalarField& source() const
            {
                return source_;
            }


        // Evolution

            //- Assemble matrix and source from the current fields, with
            //  dJdEta = G'(eta*) and etaE = eta* + E in the ion cells
            void assemble
            (
                const scalarField& dJdEta,
                const scalarField& etaE
            );

            //- Current potentials as one vector [phiE; phiI]
            tmp<scalarField> psi() const;

            //- Solve M psi = rhs, psi is the initial guess on entry
            solverPerformance solve
            (
                scalarField& psi,
                const scalarField& rhs,
                const dictionary& solverControls
            ) const;

            //- Set the internal potentials from [phiE; phiI] and correct the
            //  boundary conditions
            void scatter(const scalarField& psi);

            //- Realised current J = J* + (r/V)(eta_new - eta*) from the
            //  current potentials; it becomes the base point of the next
            //  iteration
            void realise(const scalarField& etaE);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

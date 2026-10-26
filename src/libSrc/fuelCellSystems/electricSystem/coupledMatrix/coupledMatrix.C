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

#include "coupledMatrix.H"
#include "fvCFD.H"
#include "lduPrimitiveProcessorInterface.H"
#include "processorFvPatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::coupledMatrix::buildAddressing()
{
    //- Cell pairs: ion cell -> electron cell (via the parent mesh)
    DynamicList<label> pairE;
    DynamicList<label> pairI;
    forAll(iR_.cellMapIO_, cellI)
    {
        const auto iter = eR_.cellMap_.cfind(iR_.cellMapIO_[cellI]);
        if (iter.good())
        {
            pairI.append(cellI);
            pairE.append(iter.val());
        }
    }
    pairE_.transfer(pairE);
    pairI_.transfer(pairI);
    const label nP = pairI_.size();

    //- Addressing: phiE faces, phiI faces (offset nE), one link per pair
    const labelUList& lE = eR_.lduAddr().lowerAddr();
    const labelUList& uE = eR_.lduAddr().upperAddr();
    const labelUList& lI = iR_.lduAddr().lowerAddr();
    const labelUList& uI = iR_.lduAddr().upperAddr();
    const label nF = nFE_ + nFI_ + nP;

    labelList lower(nF);
    labelList upper(nF);

    for (label f = 0; f < nFE_; ++f)
    {
        lower[f] = lE[f];
        upper[f] = uE[f];
    }
    for (label f = 0; f < nFI_; ++f)
    {
        lower[nFE_ + f] = nE_ + lI[f];
        upper[nFE_ + f] = nE_ + uI[f];
    }
    forAll(pairI_, p)
    {
        lower[nFE_ + nFI_ + p] = pairE_[p];
        upper[nFE_ + nFI_ + p] = nE_ + pairI_[p];
    }

    //- Upper-triangular face order
    oldToNew_ = lduPrimitiveMesh::upperTriOrder(nCells_, lower, upper);

    labelList lowerNew(nF);
    labelList upperNew(nF);
    forAll(oldToNew_, f)
    {
        lowerNew[oldToNew_[f]] = lower[f];
        upperNew[oldToNew_[f]] = upper[f];
    }

    //- Processor patches become interfaces of the combined matrix: their
    //- faceCells are offset like the cells, the communication is done by
    //- the processor patch fields of phi, which take the faceCells from
    //- the combined addressing (lduAddr.patchAddr).
    PtrList<const lduInterface> procInterfaces;

    auto addInterfaces = [&](const regionTypes::electric& r, label ri)
    {
        const label off = ri*nE_;
        const lduInterfaceFieldPtrsList phiInterfaces =
            r.phi_.boundaryField().scalarInterfaces();

        forAll(r.phi_.boundaryField(), patchi)
        {
            if (!r.phi_.boundaryField()[patchi].coupled())
            {
                continue;
            }

            if (!isA<processorFvPatch>(r.boundary()[patchi]))
            {
                FatalErrorInFunction
                    << "Coupled patch " << r.boundary()[patchi].name()
                    << " in region " << r.name() << " is not supported"
                    << exit(FatalError);
            }

            const processorFvPatch& pp =
                refCast<const processorFvPatch>(r.boundary()[patchi]);
            const labelUList& fc = r.boundary()[patchi].faceCells();

            labelList fcOff(fc.size());
            forAll(fc, k)
            {
                fcOff[k] = off + fc[k];
            }

            const label k = procInterfaces.size();

            procInterfaces.resize(k + 1);
            procInterfaces.set
            (
                k,
                new lduPrimitiveProcessorInterface
                (
                    fcOff,
                    pp.myProcNo(),
                    pp.neighbProcNo(),
                    pp.forwardT(),
                    pp.tag(),
                    pp.comm()
                )
            );

            interfaces_.resize(k + 1);
            interfaces_.set(k, phiInterfaces.get(patchi));

            procPatches_.append(labelPair(ri, patchi));
        }
    };

    addInterfaces(eR_, 0);
    addInterfaces(iR_, 1);

    interfaceBouCoeffs_.resize(procPatches_.size());
    interfaceIntCoeffs_.resize(procPatches_.size());

    lduInterfacePtrsList meshInterfaces(procInterfaces.size());
    forAll(procInterfaces, k)
    {
        meshInterfaces.set(k, &procInterfaces[k]);
    }

    const lduSchedule schedule
    (
        lduPrimitiveMesh::nonBlockingSchedule<processorLduInterface>
        (
            meshInterfaces
        )
    );

    //- The mesh takes ownership of procInterfaces
    lduMesh_.reset
    (
        new lduPrimitiveMesh
        (
            nCells_,
            lowerNew,
            upperNew,
            procInterfaces,
            schedule,
            UPstream::worldComm
        )
    );

    matrix_.reset(new lduMatrix(lduMesh_()));

    label nPairs = nP;
    reduce(nPairs, sumOp<label>());

    Info<< "Coupled matrix " << eR_.name() << " + " << iR_.name()
        << ": cells = " << returnReduce(nCells_, sumOp<label>())
        << ", pairs = " << nPairs
        << ", processor interfaces = "
        << returnReduce(procPatches_.size(), sumOp<label>()) << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coupledMatrix::coupledMatrix
(
    regionTypes::electric& electron,
    regionTypes::electric& ion,
    const scalar omega
)
:
    eR_(electron),
    iR_(ion),
    omega_(omega),
    nE_(electron.nCells()),
    nI_(ion.nCells()),
    nCells_(nE_ + nI_),
    nFE_(electron.lduAddr().lowerAddr().size()),
    nFI_(ion.lduAddr().lowerAddr().size()),
    source_(nCells_, Zero)
{
    buildAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::coupledMatrix::assemble
(
    const scalarField& dJdEta,
    const scalarField& etaE
)
{
    //- Segregated matrices: -laplacian(sigma, phi) - J
    fvScalarMatrix eEqn
    (
      - fvm::laplacian(eR_.sigmaField_, eR_.phi_, "laplacian(sigma,phi)")
      - eR_.j_
    );
    fvScalarMatrix iEqn
    (
      - fvm::laplacian(iR_.sigmaField_, iR_.phi_, "laplacian(sigma,phi)")
      - iR_.j_
    );

    const label nP = pairI_.size();
    const label nF = nFE_ + nFI_ + nP;

    scalarField upperCoeffs(nF);

    const scalarField& upperE = eEqn.upper();
    const scalarField& upperI = iEqn.upper();

    for (label f = 0; f < nFE_; ++f)
    {
        upperCoeffs[f] = upperE[f];
    }
    for (label f = 0; f < nFI_; ++f)
    {
        upperCoeffs[nFE_ + f] = upperI[f];
    }

    //- Diagonal and source, boundary contributions folded in
    //- (as fvMatrix::solve does with addBoundaryDiag/addBoundarySource)
    scalarField diag(nCells_);

    SubList<scalar>(diag, nE_, 0) = eEqn.diag();
    SubList<scalar>(diag, nI_, nE_) = iEqn.diag();
    SubList<scalar>(source_, nE_, 0) = eEqn.source();
    SubList<scalar>(source_, nI_, nE_) = iEqn.source();

    label k = 0;

    auto addBoundary = [&]
    (
        const regionTypes::electric& r,
        const fvScalarMatrix& eqn,
        label off
    )
    {
        forAll(r.phi_.boundaryField(), patchi)
        {
            const labelUList& fc = r.boundary()[patchi].faceCells();
            const scalarField& ic = eqn.internalCoeffs()[patchi];
            const scalarField& bc = eqn.boundaryCoeffs()[patchi];

            forAll(fc, f)
            {
                diag[off + fc[f]] += ic[f];
            }

            if (!r.phi_.boundaryField()[patchi].coupled())
            {
                forAll(fc, f)
                {
                    source_[off + fc[f]] += bc[f];
                }
                continue;
            }

            interfaceBouCoeffs_.set(k, new scalarField(bc));
            interfaceIntCoeffs_.set(k, new scalarField(ic));
            ++k;
        }
    };

    addBoundary(eR_, eEqn, 0);
    addBoundary(iR_, iEqn, nE_);

    //- Coupling terms, r = omega*V*G'(eta*):
    //  row e: (aP + r) phiE - r phiI = -V J* + r (eta* + E)
    //  row i: (aP + r) phiI - r phiE = +V J* - r (eta* + E)
    //  (the V J* parts are already in the source through "- J")
    const scalarField& V = iR_.V();

    r_.setSize(nP);

    forAll(pairI_, p)
    {
        const label e = pairE_[p];
        const label i = pairI_[p];

        r_[p] = omega_*V[i]*dJdEta[i];

        upperCoeffs[nFE_ + nFI_ + p] = -r_[p];

        diag[e] += r_[p];
        diag[nE_ + i] += r_[p];

        source_[e] += r_[p]*etaE[i];
        source_[nE_ + i] -= r_[p]*etaE[i];
    }

    //- Upper-triangular face order
    scalarField upperCoeffsNew(nF);
    forAll(oldToNew_, f)
    {
        upperCoeffsNew[oldToNew_[f]] = upperCoeffs[f];
    }

    //- Symmetric matrix: only diag and upper are set
    matrix_->diag() = diag;
    matrix_->upper() = upperCoeffsNew;
}


Foam::tmp<Foam::scalarField> Foam::coupledMatrix::psi() const
{
    tmp<scalarField> tpsi(new scalarField(nCells_));
    scalarField& psi = tpsi.ref();

    SubList<scalar>(psi, nE_, 0) = eR_.phi_.primitiveField();
    SubList<scalar>(psi, nI_, nE_) = iR_.phi_.primitiveField();

    return tpsi;
}


Foam::solverPerformance Foam::coupledMatrix::solve
(
    scalarField& psi,
    const scalarField& rhs,
    const dictionary& solverControls
) const
{
    return lduMatrix::solver::New
    (
        "phiCoupled",
        matrix_(),
        interfaceBouCoeffs_,
        interfaceIntCoeffs_,
        interfaces_,
        solverControls
    )->solve(psi, rhs);
}


void Foam::coupledMatrix::scatter(const scalarField& psi)
{
    eR_.phi_.primitiveFieldRef() = SubList<scalar>(psi, nE_, 0);
    iR_.phi_.primitiveFieldRef() = SubList<scalar>(psi, nI_, nE_);

    eR_.phi_.correctBoundaryConditions();
    iR_.phi_.correctBoundaryConditions();
}


void Foam::coupledMatrix::realise(const scalarField& etaE)
{
    const scalarField& V = iR_.V();

    scalar maxDelta = 0;
    scalar Ipos = 0;
    scalar Ineg = 0;

    forAll(pairI_, p)
    {
        const label e = pairE_[p];
        const label i = pairI_[p];

        const scalar delta = eR_.phi_[e] - iR_.phi_[i] - etaE[i];
        const scalar Jreal = iR_.j_[i] + r_[p]/V[i]*delta;

        iR_.j_[i] = Jreal;
        eR_.j_[e] = -Jreal;

        maxDelta = max(maxDelta, mag(delta));
        Ipos += max(Jreal, scalar(0))*V[i];
        Ineg += min(Jreal, scalar(0))*V[i];
    }

    eR_.j_.correctBoundaryConditions();
    iR_.j_.correctBoundaryConditions();

    label nPairs = pairI_.size();
    reduce(nPairs, sumOp<label>());
    reduce(maxDelta, maxOp<scalar>());
    reduce(Ipos, sumOp<scalar>());
    reduce(Ineg, sumOp<scalar>());

    Info<< "Coupled: pairs = " << nPairs
        << ", max|eta_new - eta*| = " << maxDelta
        << " V, I+ (A) = " << Ipos
        << ", I- (A) = " << Ineg
        << ", rel. imbalance = " << (Ipos + Ineg)/max(Ipos, VSMALL)
        << endl;
}


// ************************************************************************* //

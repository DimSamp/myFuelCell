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

Application
    regionProcAddressing

Description
    Write the processor addressing (cell, face, point and boundary
    ProcAddressing) of region meshes that were split from the parent mesh
    in parallel, after the decomposition (parallel.sh).

    Such a processor region mesh is not a decomposition of the serial region
    mesh, so decomposePar -region cannot provide its addressing, and the
    files it leaves behind belong to a different mesh. Without the right
    addressing, reconstructPar -region and decomposePar -region -fields
    scramble the fields or fail on the list sizes.

    The addressing is composed from three maps that are already on disk:

        processor region mesh -> processor parent mesh   (*RegionAddressing)
        processor parent mesh -> serial parent mesh      (*ProcAddressing)
        serial parent mesh    -> serial region mesh      (*RegionAddressing,
                                                          inverted)

    and written next to each processor region mesh. After that the standard
    tools work on the regions:

        reconstructPar -region <name>           (or -allRegions)
        decomposePar -region <name> -fields     (or -allRegions)

    The second one puts the fields of a serial time directory on processor
    region meshes that already exist, which is how a solution is moved to a
    different number of processors.

    Every map is checked before anything is written: each serial cell, face
    and point must be covered, faces must keep their patch, and only
    processor-patch faces may be flipped.

Usage
    regionProcAddressing -allRegions [-fileHandler collated]
    regionProcAddressing -region <name>

    Serial utility, like decomposePar and reconstructPar. Needs the serial
    parent and region meshes of the case (constant/polyMesh and
    constant/<region>/polyMesh) besides the processor directories.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "labelIOList.H"
#include "polyBoundaryMeshEntries.H"
#include "processorPolyPatch.H"
#include "processorCyclicPolyPatch.H"
#include "regionProperties.H"

using namespace Foam;

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Patch names, types and face ranges, read from a polyMesh/boundary file
struct patchTable
{
    wordList names;
    wordList types;
    labelList starts;
    labelList sizes;
    label nInternalFaces;

    patchTable
    (
        const Time& db,
        const word& instance,
        const fileName& meshDir,
        const label nFaces
    )
    :
        nInternalFaces(nFaces)
    {
        const polyBoundaryMeshEntries entries
        (
            IOobject
            (
                "boundary",
                instance,
                meshDir,
                db,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        );

        names.resize(entries.size());
        types.resize(entries.size());
        starts.resize(entries.size());
        sizes.resize(entries.size());

        forAll(entries, patchi)
        {
            const dictionary& dict = entries[patchi].dict();

            names[patchi] = entries[patchi].keyword();
            types[patchi] = dict.get<word>("type");
            starts[patchi] = dict.get<label>("startFace");
            sizes[patchi] = dict.get<label>("nFaces");

            nInternalFaces = min(nInternalFaces, starts[patchi]);
        }
    }

    bool processor(const label patchi) const
    {
        return
        (
            types[patchi] == processorPolyPatch::typeName
         || types[patchi] == processorCyclicPolyPatch::typeName
        );
    }
};


labelList readAddressing
(
    const Time& db,
    const word& instance,
    const fileName& meshDir,
    const word& name
)
{
    labelIOList addr
    (
        IOobject
        (
            name,
            instance,
            meshDir,
            db,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    return labelList(std::move(static_cast<labelList&>(addr)));
}


void writeAddressing
(
    const Time& db,
    const word& instance,
    const fileName& meshDir,
    const word& name,
    const labelList& addr
)
{
    labelIOList
    (
        IOobject
        (
            name,
            instance,
            meshDir,
            db,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        addr
    ).write();
}


// Invert a region -> parent map. The result holds, for every parent item,
// the region item + 1 (0: the parent item is not in the region).
// Cell and point maps hold the parent index (offset 0). Face maps hold the
// parent index + 1 (offset 1), negative when the region face is flipped,
// and the result then carries the same sign.
labelList invertAddressing(const labelUList& regionToParent, const label offset)
{
    label n = 0;
    for (const label val : regionToParent)
    {
        n = max(n, mag(val) - offset + 1);
    }

    labelList parentToRegion(n, Zero);

    forAll(regionToParent, i)
    {
        const label val = regionToParent[i];

        parentToRegion[mag(val) - offset] = (val < 0 ? -(i + 1) : (i + 1));
    }

    return parentToRegion;
}


// Entry of an inverted map, 0 when the parent item is beyond the region
inline label inverted(const labelUList& parentToRegion, const label parenti)
{
    return (parenti < parentToRegion.size() ? parentToRegion[parenti] : 0);
}


// A region -> parent map must point into the parent map it is combined with
void checkRange
(
    const word& name,
    const label proci,
    const labelUList& regionToParent,
    const label offset,
    const label parentSize
)
{
    for (const label val : regionToParent)
    {
        const label parenti = mag(val) - offset;

        if (parenti < 0 || parenti >= parentSize)
        {
            FatalErrorInFunction
                << "Processor " << proci << ": " << name
                << " points to item " << parenti
                << " of a processor parent mesh with " << parentSize
                << " items." << nl
                << "    The processor region mesh was not split from this"
                << " processor parent mesh."
                << exit(FatalError);
        }
    }
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Write the cell, face, point and boundary ProcAddressing of region"
        " meshes that were split from the parent mesh in parallel, so that"
        " reconstructPar -region and decomposePar -region -fields work"
    );

    argList::noParallel();

    #include "addAllRegionOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "getAllRegionOptions.H"

    if (regionNames.empty() || regionNames.found(polyMesh::defaultRegion))
    {
        FatalErrorInFunction
            << "Specify the regions with -region <name> or -allRegions." << nl
            << "The parent mesh already has its addressing from decomposePar."
            << exit(FatalError);
    }

    const label nProcs = fileHandler().nProcs(args.path());

    if (!nProcs)
    {
        FatalErrorInFunction
            << "No processor* directories found"
            << exit(FatalError);
    }

    // Warn fileHandler of number of processors
    const_cast<fileOperation&>(fileHandler()).nProcs(nProcs);

    Info<< "Number of processors: " << nProcs << nl << endl;

    // Processor databases
    PtrList<Time> databases(nProcs);

    forAll(databases, proci)
    {
        databases.set
        (
            proci,
            new Time
            (
                Time::controlDictName,
                args.rootPath(),
                args.caseName()/("processor" + Foam::name(proci)),
                false,  // no function objects
                false   // no libs
            )
        );

        databases[proci].setTime(runTime);
    }


    // Processor parent mesh -> serial parent mesh. Cells and points hold
    // the index, faces the index + 1 with the orientation in the sign.
    List<labelList> parentCellAddr(nProcs);
    List<labelList> parentFaceAddr(nProcs);
    List<labelList> parentPointAddr(nProcs);

    forAll(databases, proci)
    {
        const Time& db = databases[proci];
        const fileName meshDir(polyMesh::meshSubDir);
        const word instance(db.findInstance(meshDir, "faces"));

        parentCellAddr[proci] =
            readAddressing(db, instance, meshDir, "cellProcAddressing");
        parentFaceAddr[proci] =
            readAddressing(db, instance, meshDir, "faceProcAddressing");
        parentPointAddr[proci] =
            readAddressing(db, instance, meshDir, "pointProcAddressing");
    }


    for (const word& regionName : regionNames)
    {
        const fileName meshDir(regionName/polyMesh::meshSubDir);

        Info<< "Region " << regionName << endl;

        // Serial region mesh -> serial parent mesh, inverted
        const word instance(runTime.findInstance(meshDir, "faces"));

        const labelList serialCellAddr
        (
            readAddressing(runTime, instance, meshDir, "cellRegionAddressing")
        );
        const labelList serialFaceAddr
        (
            readAddressing(runTime, instance, meshDir, "faceRegionAddressing")
        );
        const labelList serialPointAddr
        (
            readAddressing(runTime, instance, meshDir, "pointRegionAddressing")
        );

        const label nCells = serialCellAddr.size();
        const label nFaces = serialFaceAddr.size();
        const label nPoints = serialPointAddr.size();

        // Parent item -> region item + 1 (faces: signed), 0 if absent
        const labelList parentToCell(invertAddressing(serialCellAddr, 0));
        const labelList parentToFace(invertAddressing(serialFaceAddr, 1));
        const labelList parentToPoint(invertAddressing(serialPointAddr, 0));

        const patchTable serialPatches(runTime, instance, meshDir, nFaces);


        // The new addressing, and how often each serial item is used
        List<labelList> cellAddr(nProcs);
        List<labelList> faceAddr(nProcs);
        List<labelList> pointAddr(nProcs);
        List<labelList> boundaryAddr(nProcs);
        wordList procInstances(nProcs);

        labelList cellCount(nCells, Zero);
        labelList faceCount(nFaces, Zero);
        labelList faceSignSum(nFaces, Zero);
        labelList pointCount(nPoints, Zero);

        label nErrors = 0;

        forAll(databases, proci)
        {
            const Time& db = databases[proci];

            procInstances[proci] = db.findInstance(meshDir, "faces");
            const word& procInstance = procInstances[proci];

            // Processor region mesh -> processor parent mesh
            const labelList regionCellAddr
            (
                readAddressing
                (
                    db, procInstance, meshDir, "cellRegionAddressing"
                )
            );
            const labelList regionFaceAddr
            (
                readAddressing
                (
                    db, procInstance, meshDir, "faceRegionAddressing"
                )
            );
            const labelList regionPointAddr
            (
                readAddressing
                (
                    db, procInstance, meshDir, "pointRegionAddressing"
                )
            );

            const patchTable procPatches
            (
                db, procInstance, meshDir, regionFaceAddr.size()
            );

            const labelList& pCellAddr = parentCellAddr[proci];
            const labelList& pFaceAddr = parentFaceAddr[proci];
            const labelList& pPointAddr = parentPointAddr[proci];

            checkRange
            (
                regionName/"cellRegionAddressing", proci,
                regionCellAddr, 0, pCellAddr.size()
            );
            checkRange
            (
                regionName/"faceRegionAddressing", proci,
                regionFaceAddr, 1, pFaceAddr.size()
            );
            checkRange
            (
                regionName/"pointRegionAddressing", proci,
                regionPointAddr, 0, pPointAddr.size()
            );

            // Cells
            labelList& cAddr = cellAddr[proci];
            cAddr.resize(regionCellAddr.size());

            forAll(cAddr, celli)
            {
                const label parentCelli =
                    pCellAddr[regionCellAddr[celli]];

                const label serialCelli =
                    inverted(parentToCell, parentCelli) - 1;

                cAddr[celli] = serialCelli;

                if (serialCelli < 0)
                {
                    ++nErrors;
                }
                else
                {
                    ++cellCount[serialCelli];
                }
            }

            // Points
            labelList& pAddr = pointAddr[proci];
            pAddr.resize(regionPointAddr.size());

            forAll(pAddr, pointi)
            {
                const label parentPointi =
                    pPointAddr[regionPointAddr[pointi]];

                const label serialPointi =
                    inverted(parentToPoint, parentPointi) - 1;

                pAddr[pointi] = serialPointi;

                if (serialPointi < 0)
                {
                    ++nErrors;
                }
                else
                {
                    ++pointCount[serialPointi];
                }
            }

            // Faces: index + 1, negative when the processor region face
            // points the other way than the serial region face
            labelList& fAddr = faceAddr[proci];
            fAddr.resize(regionFaceAddr.size());

            forAll(fAddr, facei)
            {
                const label toProcParent = regionFaceAddr[facei];
                const label toParent = pFaceAddr[mag(toProcParent) - 1];
                const label toSerial =
                    inverted(parentToFace, mag(toParent) - 1);

                if (!toSerial)
                {
                    fAddr[facei] = 0;
                    ++nErrors;
                    continue;
                }

                const label flip =
                    sign(toProcParent)*sign(toParent)*sign(toSerial);

                fAddr[facei] = flip*mag(toSerial);

                ++faceCount[mag(toSerial) - 1];
                faceSignSum[mag(toSerial) - 1] += flip;
            }

            // Internal faces stay internal and keep their orientation
            for (label facei = 0; facei < procPatches.nInternalFaces; ++facei)
            {
                if
                (
                    fAddr[facei] <= 0
                 || fAddr[facei] > serialPatches.nInternalFaces
                )
                {
                    ++nErrors;
                }
            }

            // Patches
            labelList& bAddr = boundaryAddr[proci];
            bAddr.resize(procPatches.names.size());

            forAll(bAddr, patchi)
            {
                const SubList<label> patchAddr
                (
                    fAddr,
                    procPatches.sizes[patchi],
                    procPatches.starts[patchi]
                );

                if (procPatches.processor(patchi))
                {
                    // Internal faces of the serial region mesh, either way
                    // round
                    bAddr[patchi] = -1;

                    for (const label val : patchAddr)
                    {
                        if (!val || mag(val) > serialPatches.nInternalFaces)
                        {
                            ++nErrors;
                        }
                    }

                    continue;
                }

                const label serialPatchi =
                    serialPatches.names.find(procPatches.names[patchi]);

                bAddr[patchi] = serialPatchi;

                if (serialPatchi < 0)
                {
                    if (patchAddr.size())
                    {
                        FatalErrorInFunction
                            << "Region " << regionName << ", processor "
                            << proci << ": patch "
                            << procPatches.names[patchi]
                            << " is not a patch of the serial region mesh "
                            << runTime.path()/instance/meshDir << nl
                            << "    Serial patches: "
                            << flatOutput(serialPatches.names)
                            << exit(FatalError);
                    }

                    continue;
                }

                // Same patch in the serial region mesh, same orientation
                const label first = serialPatches.starts[serialPatchi] + 1;
                const label last = first + serialPatches.sizes[serialPatchi];

                for (const label val : patchAddr)
                {
                    if (val < first || val >= last)
                    {
                        ++nErrors;
                    }
                }
            }
        }

        // Every serial cell once, every serial point at least once. Every
        // serial face once, or twice with opposite orientations when it
        // lies between two processors.
        for (const label count : cellCount)
        {
            if (count != 1)
            {
                ++nErrors;
            }
        }
        for (const label count : pointCount)
        {
            if (count < 1)
            {
                ++nErrors;
            }
        }
        forAll(faceCount, facei)
        {
            const bool once =
                (faceCount[facei] == 1 && faceSignSum[facei] == 1);

            const bool shared =
            (
                facei < serialPatches.nInternalFaces
             && faceCount[facei] == 2
             && faceSignSum[facei] == 0
            );

            if (!once && !shared)
            {
                ++nErrors;
            }
        }

        if (nErrors)
        {
            FatalErrorInFunction
                << "Region " << regionName << ": " << nErrors
                << " inconsistencies between the processor region meshes"
                << " and the serial region mesh "
                << runTime.path()/instance/meshDir << nl
                << "    Is this the serial mesh that was decomposed?" << nl
                << "    Nothing was written for this region."
                << exit(FatalError);
        }

        // Write, processor by processor: in serial the collated file
        // handler appends the processors to each file in this order
        forAll(databases, proci)
        {
            const Time& db = databases[proci];
            const word& procInstance = procInstances[proci];

            writeAddressing
            (
                db, procInstance, meshDir, "cellProcAddressing",
                cellAddr[proci]
            );
            writeAddressing
            (
                db, procInstance, meshDir, "faceProcAddressing",
                faceAddr[proci]
            );
            writeAddressing
            (
                db, procInstance, meshDir, "pointProcAddressing",
                pointAddr[proci]
            );
            writeAddressing
            (
                db, procInstance, meshDir, "boundaryProcAddressing",
                boundaryAddr[proci]
            );
        }

        fileHandler().flush();

        Info<< "    cells " << nCells
            << ", faces " << nFaces
            << ", points " << nPoints
            << ": addressing checked and written" << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

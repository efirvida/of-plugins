#include <map>
#include <tuple>

#include "Interface.H"
#include "Utilities.H"
#include "faceTriangulation.H"
#include "cellSet.H"
#include "Pstream.H"

#include <algorithm>


using namespace Foam;

preciceAdapter::Interface::Interface(
    precice::Participant& precice,
    const fvMesh& mesh,
    std::string meshName,
    std::string locationsType,
    std::vector<std::string> patchNames,
    std::vector<std::string> cellSetNames,
    bool meshConnectivity,
    bool restartFromDeformed,
    const std::string& namePointDisplacement,
    const std::string& nameCellDisplacement)
: precice_(precice),
  meshName_(meshName),
  patchNames_(patchNames),
  cellSetNames_(cellSetNames),
  meshConnectivity_(meshConnectivity),
  restartFromDeformed_(restartFromDeformed)
{
    dim_ = precice_.getMeshDimensions(meshName);

    if (dim_ == 2 && meshConnectivity_ == true)
    {
        DEBUG(adapterInfo("meshConnectivity is currently only supported for 3D cases. \n"
                          "You might set up a 3D case and restrict the 3rd dimension by z-dead = true. \n"
                          "Have a look in the adapter documentation for detailed information.",
                          "warning"));
    }

    if (locationsType == "faceCenters" || locationsType == "faceCentres")
    {
        locationType_ = LocationType::faceCenters;
    }
    else if (locationsType == "faceNodes")
    {
        locationType_ = LocationType::faceNodes;
    }
    else if (locationsType == "volumeCenters" || locationsType == "volumeCentres")
    {
        locationType_ = LocationType::volumeCenters;
    }
    else if (locationsType == "globalData" || locationsType == "global")
    {
        locationType_ = LocationType::globalData;
    }
    else
    {
        adapterInfo("Interface points location type \""
                    "locations = "
                        + locationsType + "\" is invalid.",
                    "error-deferred");
    }


    // For every patch that participates in the coupling
    for (uint j = 0; j < patchNames.size() && locationType_ != LocationType::globalData; j++)
    {
        // Get the patchID
        int patchID = mesh.boundaryMesh().findPatchID(patchNames.at(j));

        // Throw an error if the patch was not found
        if (patchID == -1)
        {
            adapterInfo("Patch \""
                            + patchNames.at(j) + "\" does not exist and therefore cannot be used as a coupling interface for mesh \""
                            + meshName + "\". Check the system/preciceDict.",
                        "error");
        }

        // Add the patch in the list
        patchIDs_.push_back(patchID);
    }

    // Configure the mesh (set the data locations)
    configureMesh(mesh, namePointDisplacement, nameCellDisplacement);
}

// ---------------------------------------------------------------------------
// Gather-scatter helper: each rank contributes its local vertex coordinates;
// rank 0 registers the full global set with preCICE and receives all global
// vertex IDs; the appropriate per-rank slice of IDs is scattered back to
// vertexIDs_ on every rank.
// In sequential (non-MPI) runs the function degenerates to a direct call.
// ---------------------------------------------------------------------------
void preciceAdapter::Interface::gatherRegisterScatterIDs(
    const std::vector<double>& localVertices)
{
    if (!Pstream::parRun())
    {
        vertexIDs_.resize(numDataLocations_);
        precice_.setMeshVertices(meshName_, localVertices, vertexIDs_);
        globalNumDataLocations_ = numDataLocations_;
        rankDataCount_.setSize(1, label(0));
        rankDataCount_[0] = label(numDataLocations_);
        rankDataOffset_.setSize(2, label(0));
        rankDataOffset_[1] = label(numDataLocations_);
        return;
    }

    // -- Step 1: exchange per-rank vertex counts --------------------------------
    rankDataCount_.setSize(Pstream::nProcs(), label(0));
    rankDataCount_[Pstream::myProcNo()] = label(numDataLocations_);
    Pstream::gatherList(rankDataCount_);
    Pstream::scatterList(rankDataCount_);

    // -- Step 2: prefix-sum offsets (identical on all ranks) -------------------
    rankDataOffset_.setSize(Pstream::nProcs() + 1, label(0));
    for (int p = 0; p < Pstream::nProcs(); ++p)
        rankDataOffset_[p + 1] = rankDataOffset_[p] + rankDataCount_[p];
    globalNumDataLocations_ = rankDataOffset_[Pstream::nProcs()];

    // -- Step 3: gather local coordinates to rank 0 ----------------------------
    List<List<double>> allVerts(Pstream::nProcs());
    {
        const label _nVerts = static_cast<label>(localVertices.size());
        List<double> localList(_nVerts);
        for (label k = 0; k < _nVerts; ++k)
            localList[k] = localVertices[k];
        allVerts[Pstream::myProcNo()] = std::move(localList);
    }
    Pstream::gatherList(allVerts);

    // -- Step 4: rank 0 registers all vertices with preCICE --------------------
    if (Pstream::master())
    {
        std::vector<double> globalVerts;
        globalVerts.reserve(
            static_cast<std::size_t>(dim_) *
            static_cast<std::size_t>(globalNumDataLocations_));
        for (int p = 0; p < Pstream::nProcs(); ++p)
            for (double d : allVerts[p])
                globalVerts.push_back(d);

        allVertexIDs_.resize(static_cast<std::size_t>(globalNumDataLocations_));
        precice_.setMeshVertices(meshName_, globalVerts, allVertexIDs_);

        DEBUG(adapterInfo(
            "Mesh \"" + meshName_ + "\": registered " +
            std::to_string(globalNumDataLocations_) +
            " vertices on rank 0 (gather-scatter mode).",
            "info"));
    }

    // -- Step 5: scatter preCICE vertex IDs back to each rank -----------------
    List<List<int>> allIDs(Pstream::nProcs());
    if (Pstream::master())
    {
        for (int p = 0; p < Pstream::nProcs(); ++p)
        {
            allIDs[p].setSize(rankDataCount_[p]);
            for (label i = 0; i < rankDataCount_[p]; ++i)
                allIDs[p][i] = allVertexIDs_[rankDataOffset_[p] + i];
        }
    }
    Pstream::scatterList(allIDs);

    const List<int>& myIDs = allIDs[Pstream::myProcNo()];
    vertexIDs_.assign(myIDs.begin(), myIDs.end());
}

// ---------------------------------------------------------------------------
void preciceAdapter::Interface::configureMesh(const fvMesh& mesh, const std::string& namePointDisplacement, const std::string& nameCellDisplacement)
{
    // Diagnostic lambda: log per-rank and global vertex counts; warn on sparse
    // interfaces (some ranks have 0 vertices while others do not).
    const auto reportVertexDistribution = [this](const std::string& locationLabel)
    {
        if (!Pstream::parRun())
            return;

        label localCount = numDataLocations_;
        label minCount   = localCount;
        label maxCount   = localCount;
        reduce(minCount, minOp<label>());
        reduce(maxCount, maxOp<label>());

        DEBUG(adapterInfo(
            "Mesh \"" + meshName_ + "\" (" + locationLabel + ") local vertices=" +
            std::to_string(localCount) +
            ", global min=" + std::to_string(minCount) +
            ", global max=" + std::to_string(maxCount),
            "info"));

        if (minCount == 0 && maxCount > 0)
        {
            DEBUG(adapterInfo(
                "Sparse coupling interface on mesh \"" + meshName_ +
                "\": some MPI ranks have no local interface entities. "
                "Handled via gather-scatter registration.",
                "warning"));
        }
    };

    // The way we configure the mesh differs between meshes based on face centers
    // and meshes based on face nodes.

    if (locationType_ == LocationType::faceCenters)
    {
        // Count the data locations for all the patches
        for (uint j = 0; j < patchIDs_.size(); j++)
        {
            numDataLocations_ +=
                mesh.boundaryMesh()[patchIDs_.at(j)].faceCentres().size();
        }
        DEBUG(adapterInfo("Number of face centres: " + std::to_string(numDataLocations_)));
        reportVertexDistribution("faceCenters");

        // In case we want to perform the reset later on, look-up the corresponding data field name
        Foam::volVectorField const* cellDisplacement = nullptr;
        if (mesh.foundObject<volVectorField>(nameCellDisplacement))
            cellDisplacement =
                &mesh.lookupObject<volVectorField>(nameCellDisplacement);

        // Array of the mesh vertices.
        // One mesh is used for all the patches and each vertex has 3D coordinates.
        std::vector<double> vertices(dim_ * numDataLocations_);

        // Array of the indices of the mesh vertices.
        // Each vertex has one index, but three coordinates.
        vertexIDs_.resize(numDataLocations_);

        // Initialize the index of the vertices array
        int verticesIndex = 0;

        // Get the locations of the mesh vertices (here: face centers)
        // for all the patches
        for (uint j = 0; j < patchIDs_.size(); j++)
        {
            // Get the face centers of the current patch
            vectorField faceCenters =
                mesh.boundaryMesh()[patchIDs_.at(j)].faceCentres();

            // Move the interface according to the current values of the cellDisplacement field,
            // to account for any displacements accumulated before restarting the simulation.
            // This is information that OpenFOAM reads from its result/restart files.
            // If the simulation is not restarted, the displacement should be zero and this line should have no effect.
            if (cellDisplacement != nullptr && !restartFromDeformed_)
                faceCenters -= cellDisplacement->boundaryField()[patchIDs_.at(j)];

            // Assign the (x,y,z) locations to the vertices
            for (int i = 0; i < faceCenters.size(); i++)
                for (unsigned int d = 0; d < dim_; ++d)
                    vertices[verticesIndex++] = faceCenters[i][d];

            // Check if we are in the right layer in case of preCICE dimension 2
            if (dim_ == 2)
            {
                const pointField faceNodes =
                    mesh.boundaryMesh()[patchIDs_.at(j)].localPoints();
                const auto faceNodesSize = faceNodes.size();
                std::array<double, 2> z_location({0, 0});
                constexpr unsigned int z_axis = 2;

                if (faceNodesSize > 0)
                    z_location[0] = faceNodes[0][z_axis];

                for (int i = 0; i < faceNodesSize; i++)
                {
                    if (z_location[0] == faceNodes[i][z_axis])
                        continue;
                    else
                    {
                        z_location[1] = faceNodes[i][z_axis];
                        break;
                    }
                }

                for (int i = 0; i < faceNodesSize; i++)
                {
                    if (z_location[0] == faceNodes[i][z_axis] || z_location[1] == faceNodes[i][z_axis])
                        continue;
                    else
                    {
                        adapterInfo("It seems like you are using preCICE in 2D and your geometry is not located int the xy-plane. "
                                    "The OpenFOAM adapter implementation supports preCICE 2D cases only with the z-axis as out-of-plane direction."
                                    "Please rotate your geometry so that the geometry is located in the xy-plane."
                                    "If you are running a 2D axisymmetric case just ignore this.",
                                    "warning");
                    }
                }
            }
        }

        // Gather all local vertices to rank 0, register with preCICE there,
        // then scatter the assigned IDs back to each rank's vertexIDs_.
        gatherRegisterScatterIDs(vertices);
    }
    else if (locationType_ == LocationType::faceNodes)
    {
        // Count the data locations for all the patches
        for (uint j = 0; j < patchIDs_.size(); j++)
        {
            numDataLocations_ +=
                mesh.boundaryMesh()[patchIDs_.at(j)].localPoints().size();
        }
        DEBUG(adapterInfo("Number of face nodes: " + std::to_string(numDataLocations_)));
        reportVertexDistribution("faceNodes");

        // In case we want to perform the reset later on, look-up the corresponding data field name
        Foam::pointVectorField const* pointDisplacement = nullptr;
        if (mesh.foundObject<pointVectorField>(namePointDisplacement))
            pointDisplacement =
                &mesh.lookupObject<pointVectorField>(namePointDisplacement);

        // Array of the mesh vertices.
        std::vector<double> vertices(dim_ * numDataLocations_);
        vertexIDs_.resize(numDataLocations_);
        int verticesIndex = 0;

        // Map between OpenFOAM vertices and preCICE vertex IDs
        std::map<std::tuple<double, double, double>, int> verticesMap;

        // Get the locations of the mesh vertices (here: face nodes) for all the patches
        for (uint j = 0; j < patchIDs_.size(); j++)
        {
            pointField faceNodes =
                mesh.boundaryMesh()[patchIDs_.at(j)].localPoints();

            if (pointDisplacement != nullptr && !restartFromDeformed_)
            {
                const vectorField& resetField = refCast<const vectorField>(
                    pointDisplacement->boundaryField()[patchIDs_.at(j)]);
                faceNodes -= resetField;
            }

            for (int i = 0; i < faceNodes.size(); i++)
                for (unsigned int d = 0; d < dim_; ++d)
                    vertices[verticesIndex++] = faceNodes[i][d];
        }

        // Gather all local vertices to rank 0, register with preCICE there,
        // then scatter the globally-assigned IDs back to each rank's vertexIDs_.
        gatherRegisterScatterIDs(vertices);

        if (meshConnectivity_)
        {
            // Rebuild verticesMap using the globally-assigned preCICE vertex IDs
            // (scattered back to this rank by gatherRegisterScatterIDs).
            for (std::size_t i = 0; i < vertexIDs_.size(); ++i)
            {
                verticesMap.emplace(
                    std::make_tuple(vertices[3 * i], vertices[3 * i + 1], vertices[3 * i + 2]),
                    vertexIDs_[i]);
            }

            // Accumulate triangle connectivity across all patches, then
            // gather to rank 0 and register with one preCICE call.
            std::vector<int> localAllTriVertIDs;
            for (uint j = 0; j < patchIDs_.size(); j++)
            {
                const int triaPerQuad = 2;
                const int nodesPerTria = 3;

                const List<face> faceField = mesh.boundaryMesh()[patchIDs_.at(j)].localFaces();
                Field<point> pointCoords   = mesh.boundaryMesh()[patchIDs_.at(j)].localPoints();

                if (pointDisplacement != nullptr && !restartFromDeformed_)
                {
                    const vectorField& resetField = refCast<const vectorField>(
                        pointDisplacement->boundaryField()[patchIDs_.at(j)]);
                    pointCoords -= resetField;
                }

                std::vector<int> triVertIDs;
                triVertIDs.reserve(faceField.size() * triaPerQuad * nodesPerTria);

                forAll(faceField, facei)
                {
                    const face& faceQuad = faceField[facei];
                    faceTriangulation faceTri(pointCoords, faceQuad, false);

                    for (uint triIndex = 0; triIndex < triaPerQuad; triIndex++)
                        for (uint nodeIndex = 0; nodeIndex < nodesPerTria; nodeIndex++)
                            triVertIDs.push_back(verticesMap.at(std::make_tuple(
                                pointCoords[faceTri[triIndex][nodeIndex]][0],
                                pointCoords[faceTri[triIndex][nodeIndex]][1],
                                pointCoords[faceTri[triIndex][nodeIndex]][2])));
                }

                DEBUG(adapterInfo("Number of triangles: " + std::to_string(faceField.size() * triaPerQuad)));

                localAllTriVertIDs.insert(localAllTriVertIDs.end(),
                    triVertIDs.begin(), triVertIDs.end());
            }

            // Register mesh triangles: gather all to rank 0, then one setMeshTriangles call.
            if (Pstream::parRun())
            {
                List<List<int>> allTriIDs(Pstream::nProcs());
                {
                    const label _nTri = static_cast<label>(localAllTriVertIDs.size());
                    List<int> localList(_nTri);
                    for (label k = 0; k < _nTri; ++k)
                        localList[k] = localAllTriVertIDs[k];
                    allTriIDs[Pstream::myProcNo()] = std::move(localList);
                }
                Pstream::gatherList(allTriIDs);
                if (Pstream::master())
                {
                    std::vector<int> globalTriIDs;
                    for (int p = 0; p < Pstream::nProcs(); ++p)
                        for (int id : allTriIDs[p])
                            globalTriIDs.push_back(id);
                    if (!globalTriIDs.empty())
                        precice_.setMeshTriangles(meshName_, globalTriIDs);
                }
            }
            else
            {
                if (!localAllTriVertIDs.empty())
                    precice_.setMeshTriangles(meshName_, localAllTriVertIDs);
            }
        }
    }
    else if (locationType_ == LocationType::volumeCenters)
    {
        // The volume coupling implementation considers the mesh points in the volume and
        // on the boundary patches in order to take the boundary conditions into account

        std::vector<labelList> overlapCells;

        if (!cellSetNames_.empty())
        {
            for (uint j = 0; j < cellSetNames_.size(); j++)
            {
                cellSet overlapRegion(mesh, cellSetNames_[j]);
                overlapCells.push_back(overlapRegion.toc());
                numDataLocations_ += overlapCells[j].size();
            }
        }
        else
        {
            numDataLocations_ = mesh.C().size();
        }

        for (uint j = 0; j < patchIDs_.size(); j++)
        {
            numDataLocations_ +=
                mesh.boundaryMesh()[patchIDs_.at(j)].faceCentres().size();
        }
        DEBUG(adapterInfo("Number of coupling volumes: " + std::to_string(numDataLocations_)));
        reportVertexDistribution("volumeCenters");

        std::vector<double> vertices(dim_ * numDataLocations_);
        vertexIDs_.resize(numDataLocations_);
        int verticesIndex = 0;

        if (!cellSetNames_.empty())
        {
            for (uint j = 0; j < cellSetNames_.size(); j++)
            {
                const labelList& cells = overlapCells.at(j);
                for (int i = 0; i < cells.size(); i++)
                {
                    vertices[verticesIndex++] = mesh.C().internalField()[cells[i]].x();
                    vertices[verticesIndex++] = mesh.C().internalField()[cells[i]].y();
                    if (dim_ == 3)
                        vertices[verticesIndex++] = mesh.C().internalField()[cells[i]].z();
                }
            }
        }
        else
        {
            const vectorField& CellCenters = mesh.C();
            for (int i = 0; i < CellCenters.size(); i++)
            {
                vertices[verticesIndex++] = CellCenters[i].x();
                vertices[verticesIndex++] = CellCenters[i].y();
                if (dim_ == 3)
                    vertices[verticesIndex++] = CellCenters[i].z();
            }
        }

        for (uint j = 0; j < patchIDs_.size(); j++)
        {
            const vectorField faceCenters =
                mesh.boundaryMesh()[patchIDs_.at(j)].faceCentres();
            for (int i = 0; i < faceCenters.size(); i++)
            {
                vertices[verticesIndex++] = faceCenters[i].x();
                vertices[verticesIndex++] = faceCenters[i].y();
                if (dim_ == 3)
                    vertices[verticesIndex++] = faceCenters[i].z();
            }
        }

        // Gather all local vertices to rank 0, register with preCICE there,
        // then scatter the assigned IDs back to each rank's vertexIDs_.
        gatherRegisterScatterIDs(vertices);
    }
    else if (locationType_ == LocationType::globalData)
    {
        numDataLocations_ = 1;
        vertexIDs_.clear();

        if (Pstream::master())
        {
            std::vector<double> vertices(dim_, 0.0);
            vertexIDs_.resize(1);
            precice_.setMeshVertices(meshName_, vertices, vertexIDs_);
        }
    }
}


void preciceAdapter::Interface::addCouplingDataWriter(
    const FieldConfig& fieldConfig,
    CouplingDataUser* couplingDataWriter)
{
    couplingDataWriter->setDataName(fieldConfig.name);
    couplingDataWriter->setFlipNormal(fieldConfig.flip_normal);
    couplingDataWriter->setPatchIDs(patchIDs_);
    couplingDataWriter->setCellSetNames(cellSetNames_);
    couplingDataWriter->setLocationsType(locationType_);
    couplingDataWriter->checkDataLocation(meshConnectivity_);
    couplingDataWriter->initialize();
    couplingDataWriters_.push_back(couplingDataWriter);
}


void preciceAdapter::Interface::addCouplingDataReader(
    const FieldConfig& fieldConfig,
    preciceAdapter::CouplingDataUser* couplingDataReader)
{
    couplingDataReader->setDataName(fieldConfig.name);
    couplingDataReader->setFlipNormal(fieldConfig.flip_normal);
    couplingDataReader->setPatchIDs(patchIDs_);
    couplingDataReader->setLocationsType(locationType_);
    couplingDataReader->setCellSetNames(cellSetNames_);
    couplingDataReader->checkDataLocation(meshConnectivity_);
    couplingDataReader->initialize();
    couplingDataReaders_.push_back(couplingDataReader);
}

void preciceAdapter::Interface::createBuffer()
{
    bool needsVectorData = false;
    int dataBufferSize = 0;

    for (uint i = 0; i < couplingDataReaders_.size(); i++)
        if (couplingDataReaders_.at(i)->hasVectorData())
            needsVectorData = true;

    for (uint i = 0; i < couplingDataWriters_.size(); i++)
        if (couplingDataWriters_.at(i)->hasVectorData())
            needsVectorData = true;

    dataBufferSize = needsVectorData ? dim_ * numDataLocations_ : numDataLocations_;

    // An interface has only one data buffer, shared between several CouplingDataUsers.
    dataBuffer_.resize(dataBufferSize);
}

// ---------------------------------------------------------------------------
// readCouplingData
// For parallel non-global data: rank 0 reads the full global dataset from
// preCICE using allVertexIDs_, then scatters per-rank slices to every rank.
// For globalData: unchanged master-read + broadcast semantics.
// For sequential: unchanged direct read.
// ---------------------------------------------------------------------------
void preciceAdapter::Interface::readCouplingData(double relativeReadTime)
{
    for (uint i = 0; i < couplingDataReaders_.size(); i++)
    {
        preciceAdapter::CouplingDataUser* couplingDataReader = couplingDataReaders_.at(i);
        const int dataDim =
            static_cast<int>(precice_.getDataDimensions(meshName_, couplingDataReader->dataName()));

        if (locationType_ == LocationType::globalData)
        {
            // Global data: master reads, broadcasts to all ranks
            const std::size_t nGlobal = static_cast<std::size_t>(dataDim);
            precice::span<double> dataSpan {dataBuffer_.data(), nGlobal};
            if (Pstream::master())
            {
                precice_.readData(meshName_, couplingDataReader->dataName(),
                    vertexIDs_, relativeReadTime, dataSpan);
            }
            if (Pstream::parRun())
            {
                for (std::size_t k = 0; k < nGlobal; ++k)
                    Pstream::broadcast(dataBuffer_[k]);
            }
            couplingDataReader->applyFlipNormal(dataSpan);
            couplingDataReader->read(dataBuffer_.data(), dim_);
        }
        else if (Pstream::parRun())
        {
            // Parallel non-global: rank 0 reads all, scatter to each rank
            const std::size_t nLocalRead =
                static_cast<std::size_t>(numDataLocations_) *
                static_cast<std::size_t>(dataDim);

            List<List<double>> allBufs(Pstream::nProcs());
            if (Pstream::master())
            {
                std::vector<double> globalData(
                    static_cast<std::size_t>(globalNumDataLocations_) *
                    static_cast<std::size_t>(dataDim));
                precice_.readData(
                    meshName_,
                    couplingDataReader->dataName(),
                    allVertexIDs_,
                    relativeReadTime,
                    {globalData.data(), globalData.size()});

                // Split into per-rank slices
                for (int p = 0; p < Pstream::nProcs(); ++p)
                {
                    const label start = static_cast<label>(dataDim) * rankDataOffset_[p];
                    const label sz    = static_cast<label>(dataDim) * rankDataCount_[p];
                    allBufs[p].setSize(sz);
                    for (label k = 0; k < sz; ++k)
                        allBufs[p][k] = globalData[start + k];
                }
            }
            Pstream::scatterList(allBufs);

            // Copy scattered slice into this rank's dataBuffer_
            const List<double>& myBuf = allBufs[Pstream::myProcNo()];
            std::copy(myBuf.begin(), myBuf.end(), dataBuffer_.begin());

            precice::span<double> localSpan {dataBuffer_.data(), nLocalRead};
            couplingDataReader->applyFlipNormal(localSpan);
            couplingDataReader->read(dataBuffer_.data(), dim_);
        }
        else
        {
            // Sequential run: direct read
            const std::size_t nRead =
                vertexIDs_.size() * static_cast<std::size_t>(dataDim);
            precice::span<double> dataSpan {dataBuffer_.data(), nRead};
            precice_.readData(meshName_, couplingDataReader->dataName(),
                vertexIDs_, relativeReadTime, dataSpan);
            couplingDataReader->applyFlipNormal(dataSpan);
            couplingDataReader->read(dataBuffer_.data(), dim_);
        }
    }
}

// ---------------------------------------------------------------------------
// writeCouplingData
// For parallel non-global data: each rank fills its local dataBuffer_ via
// write(), then all buffers are gathered to rank 0 which writes the assembled
// global dataset to preCICE via allVertexIDs_.
// For globalData: unchanged broadcast + consistency-check + master-write.
// For sequential: unchanged direct write.
// ---------------------------------------------------------------------------
void preciceAdapter::Interface::writeCouplingData()
{
    for (uint i = 0; i < couplingDataWriters_.size(); i++)
    {
        preciceAdapter::CouplingDataUser* couplingDataWriter = couplingDataWriters_.at(i);
        const int dataDim =
            static_cast<int>(precice_.getDataDimensions(meshName_, couplingDataWriter->dataName()));

        // Populate the local data buffer from OpenFOAM fields
        auto nWrittenData = couplingDataWriter->write(dataBuffer_.data(), meshConnectivity_, dim_);

        if (locationType_ == LocationType::globalData)
        {
            // Global data: check consistency, master writes
            if (Pstream::parRun())
            {
                std::vector<double> localBuffer(
                    dataBuffer_.begin(), dataBuffer_.begin() + nWrittenData);
                for (std::size_t k = 0; k < nWrittenData; ++k)
                    Pstream::broadcast(dataBuffer_[k]);
                scalar maxDiff = 0.0;
                for (std::size_t k = 0; k < nWrittenData; ++k)
                    maxDiff = max(maxDiff, mag(localBuffer[k] - dataBuffer_[k]));
                reduce(maxDiff, maxOp<scalar>());
                if (maxDiff > SMALL)
                {
                    adapterInfo(
                        "Global data \"" + couplingDataWriter->dataName() +
                        "\" differs across MPI ranks. Use a single shared value on all ranks "
                        "before writing to preCICE.",
                        "error");
                }
            }
            const std::size_t nWrite =
                vertexIDs_.size() * static_cast<std::size_t>(dataDim);
            precice::span<double> dataSpan {dataBuffer_.data(), nWrite};
            couplingDataWriter->applyFlipNormal(dataSpan);
            if (Pstream::master())
                precice_.writeData(meshName_, couplingDataWriter->dataName(),
                    vertexIDs_, dataSpan);
        }
        else if (Pstream::parRun())
        {
            // Parallel non-global: apply flip normal, gather, rank 0 writes all
            const std::size_t nLocalWrite =
                static_cast<std::size_t>(numDataLocations_) *
                static_cast<std::size_t>(dataDim);
            precice::span<double> localSpan {dataBuffer_.data(), nLocalWrite};
            couplingDataWriter->applyFlipNormal(localSpan);

            List<List<double>> allBufs(Pstream::nProcs());
            {
                const label _nCopy = static_cast<label>(nLocalWrite);
                List<double> localList(_nCopy);
                for (label k = 0; k < _nCopy; ++k)
                    localList[k] = dataBuffer_[k];
                allBufs[Pstream::myProcNo()] = std::move(localList);
            }
            Pstream::gatherList(allBufs);

            if (Pstream::master())
            {
                std::vector<double> globalData;
                globalData.reserve(
                    static_cast<std::size_t>(globalNumDataLocations_) *
                    static_cast<std::size_t>(dataDim));
                for (int p = 0; p < Pstream::nProcs(); ++p)
                    for (double v : allBufs[p])
                        globalData.push_back(v);
                precice_.writeData(
                    meshName_,
                    couplingDataWriter->dataName(),
                    allVertexIDs_,
                    {globalData.data(), globalData.size()});
            }
        }
        else
        {
            // Sequential run: direct write
            const std::size_t nWrite =
                vertexIDs_.size() * static_cast<std::size_t>(dataDim);
            precice::span<double> dataSpan {dataBuffer_.data(), nWrite};
            couplingDataWriter->applyFlipNormal(dataSpan);
            precice_.writeData(meshName_, couplingDataWriter->dataName(),
                vertexIDs_, dataSpan);
        }
    }
}

preciceAdapter::Interface::~Interface()
{
    for (uint i = 0; i < couplingDataReaders_.size(); i++)
        delete couplingDataReaders_.at(i);
    couplingDataReaders_.clear();

    for (uint i = 0; i < couplingDataWriters_.size(); i++)
        delete couplingDataWriters_.at(i);
    couplingDataWriters_.clear();
}

preciceAdapter::LocationType preciceAdapter::Interface::locationType() const
{
    return locationType_;
}

unsigned int preciceAdapter::Interface::dataDimensions(const std::string& dataName) const
{
    return precice_.getDataDimensions(meshName_, dataName);
}

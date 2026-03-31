#ifndef INTERFACE_H
#define INTERFACE_H

#include <memory>
#include <string>
#include <vector>
#include "fvCFD.H"
#include "CouplingDataUser.H"
#include <precice/precice.hpp>

#include "pointPatchField.H"

namespace preciceAdapter
{

//- field options
struct FieldConfig
{
    std::string name;
    std::string solver_name;
    std::string operation;
    bool flip_normal;
    Foam::dimensionSet dimensions;
    bool has_dimensions = false;
};

class Interface
{
protected:
    //- preCICE solver interface
    precice::Participant& precice_;

    //- Mesh name used in the preCICE configuration
    std::string meshName_;

    //- Point locations type (faceCenters vs faceNodes)
    LocationType locationType_ {LocationType::none};

    //- Mesh ID assigned by preCICE to the interface
    int meshID_;

    //- Names of the OpenFOAM patches that form the interface
    std::vector<std::string> patchNames_;

    //- OpenFOAM patches that form the interface
    std::vector<int> patchIDs_;

    //- Names of the OpenFOAM cell sets to be coupled (for volume coupling)
    std::vector<std::string> cellSetNames_;

    //- Number of data points (cell centers) on the interface (local to this rank)
    int numDataLocations_ = 0;

    //- Vertex IDs assigned by preCICE (local slice, scatter-received from rank 0)
    std::vector<int> vertexIDs_;

    //- Per-rank vertex count (all ranks, size = nProcs). Valid after configureMesh.
    Foam::labelList rankDataCount_;

    //- Prefix-sum offsets into the global vertex array (size = nProcs+1, all ranks).
    Foam::labelList rankDataOffset_;

    //- Total vertex count across all ranks. Valid after configureMesh.
    Foam::label globalNumDataLocations_ = 0;

    //- All preCICE vertex IDs for this mesh (rank 0 only, size = globalNumDataLocations_).
    std::vector<int> allVertexIDs_;

    //- Buffer for the coupling data
    std::vector<double> dataBuffer_;

    //- Vector of CouplingDataReaders
    std::vector<std::unique_ptr<CouplingDataUser>> couplingDataReaders_;

    //- Vector of CouplingDataWriters
    std::vector<std::unique_ptr<CouplingDataUser>> couplingDataWriters_;

    //Switch for faceTriangulation (nearest projection)
    bool meshConnectivity_;

    //- Reset the displacement during interface definition
    bool restartFromDeformed_;

    // Simulation dimension
    unsigned int dim_;

    //- Gather local mesh vertices to rank 0, call setMeshVertices there, then
    //  scatter the preCICE vertex IDs back to each rank's vertexIDs_.
    //  In sequential runs calls setMeshVertices directly.
    void gatherRegisterScatterIDs(const std::vector<double>& localVertices);

    //- Extracts the locations of the face centers or face nodes
    //  and exposes them to preCICE via gatherRegisterScatterIDs.
    void configureMesh(const Foam::fvMesh& mesh,
                       const std::string& namePointDisplacement,
                       const std::string& nameCellDisplacement);

public:
    //- Constructor
    Interface(
        precice::Participant& precice,
        const Foam::fvMesh& mesh,
        const std::string& meshName,
        const std::string& locationsType,
        const std::vector<std::string>& patchNames,
        const std::vector<std::string>& cellSetNames,
        bool meshConnectivity,
        bool restartFromDeformed,
        const std::string& namePointDisplacement,
        const std::string& nameCellDisplacement);

    //- Add a CouplingDataUser to read data from the interface
    void addCouplingDataReader(
        const FieldConfig& fieldConfig,
        CouplingDataUser* couplingDataReader);


    //- Add a CouplingDataUser to write data on the interface
    void addCouplingDataWriter(
        const FieldConfig& fieldConfig,
        CouplingDataUser* couplingDataWriter);

    //- Allocate an appropriate buffer for scalar or vector data.
    //  If at least one couplingDataUser has vector data, then
    //  define a buffer for 3D data. Otherwise, for 1D data.
    void createBuffer();

    //- Call read() on each registered couplingDataReader to read the coupling
    //  data from the buffer and apply the boundary conditions
    void readCouplingData(double relativeReadTime);

    //- Call write() on each registered couplingDataWriter to extract the boundary
    //  data and write them into the buffer
    void writeCouplingData();

    //- Return the configured interface location type
    LocationType locationType() const;

    //- Return the preCICE data dimensions for a configured data set
    unsigned int dataDimensions(const std::string& dataName) const noexcept;

    //- Destructor
    ~Interface();
};

}

#endif

#include "DisplacementDelta.H"

using namespace Foam;

preciceAdapter::FSI::DisplacementDelta::DisplacementDelta(
    const Foam::fvMesh& mesh,
    const std::string namePointDisplacement,
    const std::string nameCellDisplacement)
: pointDisplacement_(
    const_cast<pointVectorField*>(
        &mesh.lookupObject<pointVectorField>(namePointDisplacement))),
  cellDisplacement_(
      const_cast<volVectorField*>(
          &mesh.lookupObject<volVectorField>(nameCellDisplacement))),
  mesh_(mesh)
{
    dataType_ = vector;
}

// We cannot do this step in the constructor by design of the adapter since the information of the CouplingDataUser is
// defined later. Hence, we call this method after the CouplingDaaUser has been configured
void preciceAdapter::FSI::DisplacementDelta::initialize()
{
    // Initialize appropriate objects for each interface patch, namely the volField and the interpolation object
    // this is only necessary for face based FSI
    if (this->locationType_ == LocationType::faceCenters)
    {
        for (unsigned int j = 0; j < patchIDs_.size(); ++j)
        {
            const unsigned int patchID = patchIDs_.at(j);
            interpolationObjects_.emplace_back(new primitivePatchInterpolation(mesh_.boundaryMesh()[patchID]));
        }
    }
}


std::size_t preciceAdapter::FSI::DisplacementDelta::write(double* buffer, bool meshConnectivity, const unsigned int dim)
{
    /* TODO: Implement
     * We need two nested for-loops for each patch,
     * the outer for the locations and the inner for the dimensions.
     * See the preCICE writeBlockVectorData() implementation.
     */
    adapterInfo("Writing displacementDelta is not supported.", "error");
    return 0;
}

// return the displacement to use later in the velocity?
void preciceAdapter::FSI::DisplacementDelta::read(double* buffer, const unsigned int dim)
{
    int bufferIndex = 0;
    for (unsigned int j = 0; j < patchIDs_.size(); j++)
    {
        // Get the ID of the current patch
        const unsigned int patchID = patchIDs_.at(j);

        if (this->locationType_ == LocationType::faceCenters)
        {

            // the boundaryCellDisplacement is a vector and ordered according to the iterator j
            // and not according to the patchID
            // First, copy the buffer data into the center based vectorFields on each interface patch
            // For DisplacementDelta, set absolute values here and sum the interpolated values up to the point field
            // since the temporary field in this class is not reloaded in the implicit coupling
            forAll(cellDisplacement_->boundaryField()[patchID], i)
            {
                for (unsigned int d = 0; d < dim; ++d)
                    cellDisplacement_->boundaryFieldRef()[patchID][i][d] = buffer[bufferIndex++];
            }

            // FSI-DISP-LOG: log incoming cell displacement from preCICE buffer
            {
                const vectorField& cellBF = cellDisplacement_->boundaryField()[patchID];
                Info << "[FSI-DISP-LOG] DisplacementDelta::read"
                     << " T=" << mesh_.time().timeName()
                     << " patch=" << mesh_.boundary()[patchID].name()
                     << " cellDisp_buffer: min=" << gMin(cellBF)
                     << " max=" << gMax(cellBF)
                     << " maxMag=" << gMax(mag(cellBF))
                     << endl;
            }

            // Get a reference to the displacement on the point patch in order to overwrite it
            vectorField& pointDisplacementFluidPatch(
                refCast<vectorField>(
                    pointDisplacement_->boundaryFieldRef()[patchID]));

            // Overwrite the node based patch using the interpolation objects and the cell based vector field
            // Afterwards, continue as usual
            pointDisplacementFluidPatch += interpolationObjects_[j]->faceToPointInterpolate(cellDisplacement_->boundaryField()[patchID]);

            // FSI-DISP-LOG: log accumulated point displacement after += 
            {
                Info << "[FSI-DISP-LOG] DisplacementDelta::read"
                     << " T=" << mesh_.time().timeName()
                     << " patch=" << mesh_.boundary()[patchID].name()
                     << " pointDisp_accumulated: min=" << gMin(pointDisplacementFluidPatch)
                     << " max=" << gMax(pointDisplacementFluidPatch)
                     << " maxMag=" << gMax(mag(pointDisplacementFluidPatch))
                     << endl;
            }
        }
        else if (this->locationType_ == LocationType::faceNodes)
        {

            // Get the displacement on the patch
            fixedValuePointPatchVectorField& pointDisplacementFluidPatch(
                refCast<fixedValuePointPatchVectorField>(
                    pointDisplacement_->boundaryFieldRef()[patchID]));

            // Overwrite the nodes on the interface directly
            forAll(pointDisplacement_->boundaryFieldRef()[patchID], i)
            {
                for (unsigned int d = 0; d < dim; ++d)
                    pointDisplacementFluidPatch[i][d] += buffer[bufferIndex++];
            }

            // FSI-DISP-LOG: log accumulated point displacement (faceNodes)
            {
                Info << "[FSI-DISP-LOG] DisplacementDelta::read (faceNodes)"
                     << " T=" << mesh_.time().timeName()
                     << " patch=" << mesh_.boundary()[patchID].name()
                     << " pointDisp_accumulated: min=" << gMin(static_cast<const vectorField&>(pointDisplacementFluidPatch))
                     << " max=" << gMax(static_cast<const vectorField&>(pointDisplacementFluidPatch))
                     << " maxMag=" << gMax(mag(static_cast<const vectorField&>(pointDisplacementFluidPatch)))
                     << endl;
            }
        }
    }
}

bool preciceAdapter::FSI::DisplacementDelta::isLocationTypeSupported(const bool meshConnectivity) const
{
    // Solid solver *could* allow connectivity for writing displacement
    if (meshConnectivity)
    {
        return (this->locationType_ == LocationType::faceNodes);
    }
    else
    {
        return (this->locationType_ == LocationType::faceCenters || this->locationType_ == LocationType::faceNodes);
    }
}

std::string preciceAdapter::FSI::DisplacementDelta::getDataName() const
{
    return "DisplacementDelta";
}

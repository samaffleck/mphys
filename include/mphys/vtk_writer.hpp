#pragma once

#include <string>
#include <vector>

#include "mphys/mesh_model.hpp"
#include "mphys/topology.hpp"

namespace mphys {

// Write a structured mesh and its cell-centred fields to a VTK ImageData file
// (.vti, ASCII) that ParaView/VTK can open directly. Works for 1D/2D/3D
// structured meshes; the cell ordering (i fastest, then j, then k) already
// matches VTK's expectation, so no connectivity is emitted.
//
// `field_names[k]` labels `field_data[k]`, each a cell array of length
// mesh.NCells(). Throws std::invalid_argument if the mesh is not structured or
// the array sizes are inconsistent.
void WriteVtkImageData(const std::string& path, const Mesh& mesh,
                       const std::vector<std::string>& field_names,
                       const std::vector<const std::vector<double>*>& field_data);

// Convenience: dump every field of a MeshModel to `path` (a .vti file).
void WriteVtk(const std::string& path, const MeshModel& model);

}  // namespace mphys

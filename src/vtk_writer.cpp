#include "mphys/vtk_writer.hpp"

#include <fstream>
#include <stdexcept>

namespace mphys {

void WriteVtkImageData(
    const std::string& path, const Mesh& mesh,
    const std::vector<std::string>& field_names,
    const std::vector<const std::vector<double>*>& field_data) {
  const StructuredInfo& s = mesh.structured;
  if (!s.valid) {
    throw std::invalid_argument(
        "WriteVtkImageData: mesh is not a structured grid");
  }
  if (field_names.size() != field_data.size()) {
    throw std::invalid_argument("WriteVtkImageData: names/data size mismatch");
  }
  const int n_cells = mesh.NCells();
  for (const auto* d : field_data) {
    if (static_cast<int>(d->size()) != n_cells) {
      throw std::invalid_argument("WriteVtkImageData: field size != NCells");
    }
  }

  std::ofstream out(path);
  if (!out) throw std::runtime_error("WriteVtkImageData: cannot open " + path);

  // WholeExtent is given in points: nx/ny/nz cells -> nx/ny/nz point spans.
  out << "<?xml version=\"1.0\"?>\n";
  out << "<VTKFile type=\"ImageData\" version=\"1.0\" "
         "byte_order=\"LittleEndian\">\n";
  out << "  <ImageData WholeExtent=\"0 " << s.nx << " 0 " << s.ny << " 0 "
      << s.nz << "\" Origin=\"" << s.x0 << " " << s.y0 << " " << s.z0
      << "\" Spacing=\"" << s.dx << " " << s.dy << " " << s.dz << "\">\n";
  out << "    <Piece Extent=\"0 " << s.nx << " 0 " << s.ny << " 0 " << s.nz
      << "\">\n";

  out << "      <CellData";
  if (!field_names.empty()) out << " Scalars=\"" << field_names[0] << "\"";
  out << ">\n";
  for (std::size_t k = 0; k < field_names.size(); ++k) {
    out << "        <DataArray type=\"Float64\" Name=\"" << field_names[k]
        << "\" format=\"ascii\">\n          ";
    const auto& v = *field_data[k];
    for (std::size_t i = 0; i < v.size(); ++i) {
      out << v[i] << (i + 1 < v.size() ? ' ' : '\n');
    }
    out << "        </DataArray>\n";
  }
  out << "      </CellData>\n";

  out << "    </Piece>\n";
  out << "  </ImageData>\n";
  out << "</VTKFile>\n";
}

void WriteVtk(const std::string& path, const MeshModel& model) {
  std::vector<std::string> names;
  std::vector<const std::vector<double>*> data;
  names.reserve(model.NFields());
  data.reserve(model.NFields());
  for (int k = 0; k < model.NFields(); ++k) {
    names.push_back(model.field_name(k));
    data.push_back(&model.fields()[k]);
  }
  WriteVtkImageData(path, model.mesh(), names, data);
}

}  // namespace mphys

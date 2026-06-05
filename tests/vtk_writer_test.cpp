#include <gtest/gtest.h>

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "mphys/topology.hpp"
#include "mphys/vtk_writer.hpp"

namespace {

std::string ReadFile(const std::string& path) {
  std::ifstream in(path);
  std::stringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

std::string TempPath(const std::string& name) {
  return (std::filesystem::temp_directory_path() / name).string();
}

}  // namespace

// The VTI writer must emit a well-formed ImageData file with the grid extent,
// spacing, and a DataArray carrying every cell value in order.
TEST(VtkWriter, ImageData2D) {
  const mphys::Mesh mesh = mphys::MakeStructuredMesh2D(0.0, 2.0, 3, 0.0, 1.0, 2);
  ASSERT_TRUE(mesh.structured.valid);

  std::vector<double> u(mesh.NCells());
  for (int c = 0; c < mesh.NCells(); ++c) u[c] = static_cast<double>(c);

  const std::string path = TempPath("mphys_vtk_test.vti");
  mphys::WriteVtkImageData(path, mesh, {"u"}, {&u});

  const std::string text = ReadFile(path);
  EXPECT_NE(text.find("type=\"ImageData\""), std::string::npos);
  EXPECT_NE(text.find("WholeExtent=\"0 3 0 2 0 1\""), std::string::npos);
  EXPECT_NE(text.find("Spacing=\"" ), std::string::npos);
  EXPECT_NE(text.find("Name=\"u\""), std::string::npos);
  // Every cell value should appear (0..5).
  EXPECT_NE(text.find("0 1 2 3 4 5"), std::string::npos);

  std::filesystem::remove(path);
}

// A non-structured mesh has no grid metadata and must be rejected clearly.
TEST(VtkWriter, RejectsUnstructured) {
  mphys::Mesh mesh;  // default: structured.valid == false
  std::vector<double> u;
  EXPECT_THROW(mphys::WriteVtkImageData(TempPath("x.vti"), mesh, {"u"}, {&u}),
               std::invalid_argument);
}

// Copyright 2026 Spencer Evans-Cole

#include "aether/output/output.hpp"

#include <string>
#include <vector>
// Appends a <DataArray> element with whitespace-separated text content.
// ncomp <= 0 omits the NumberOfComponents attribute.
inline void append_data_array(pugi::xml_node parent, const char* type, const char* name, int ncomp,
                              const std::string& content) {
  pugi::xml_node da = parent.append_child("DataArray");
  da.append_attribute("type") = type;
  if (name != nullptr) da.append_attribute("Name") = name;
  if (ncomp > 0) da.append_attribute("NumberOfComponents") = ncomp;
  da.append_attribute("format") = "ascii";
  da.text().set(content.c_str());
}

namespace aether::output {
// Writes a serial VTK XML UnstructuredGrid (.vtu) for a 2D P1 triangular mesh.
//   points       : one (x, y) per node; z is written as 0.
//   triangles    : one triple of 0-based global node indices per element.
//   point_fields : named nodal scalar fields (e.g. {"u", solution}); each
//               must have the same size as points.
void Output::write_solution_vtk(
    const std::string& path, const std::vector<Vec2>& points,
    const std::vector<std::array<NodeIndex, 3>>& triangles,
    const std::vector<std::pair<std::string, Eigen::VectorXd>>& point_fields) const {
  pugi::xml_document doc;

  pugi::xml_node decl = doc.append_child(pugi::node_declaration);
  decl.append_attribute("version") = "1.0";

  pugi::xml_node vtkfile = doc.append_child("VTKFile");
  vtkfile.append_attribute("type") = "UnstructuredGrid";
  vtkfile.append_attribute("version") = "1.0";
  vtkfile.append_attribute("byte_order") = "LittleEndian";

  pugi::xml_node piece = vtkfile.append_child("UnstructuredGrid").append_child("Piece");
  piece.append_attribute("NumberOfPoints") =
      static_cast<unsigned long long>(points.size());  // NOLINT
  piece.append_attribute("NumberOfCells") =
      static_cast<unsigned long long>(triangles.size());  // NOLINT

  // Points -- VTK always expects 3 components.
  {
    std::ostringstream oss;
    oss.precision(17);
    for (const auto& p : points) oss << p[0] << ' ' << p[1] << " 0\n";
    append_data_array(piece.append_child("Points"), "Float64", nullptr, 3, oss.str());
  }

  // Cells -- connectivity, offsets, and per-cell type (5 = VTK_TRIANGLE).
  {
    pugi::xml_node cells = piece.append_child("Cells");

    std::ostringstream conn;
    for (const auto& t : triangles) {
      int indices[3] = {t[0].get(), t[1].get(), t[2].get()};
      conn << indices[0] << ' ' << indices[1] << ' ' << indices[2] << '\n';
    }
    append_data_array(cells, "Int32", "connectivity", 0, conn.str());

    std::ostringstream offsets;
    for (std::size_t c = 0; c < triangles.size(); ++c) offsets << 3 * (c + 1) << ' ';
    append_data_array(cells, "Int32", "offsets", 0, offsets.str());

    std::ostringstream types;
    for (std::size_t c = 0; c < triangles.size(); ++c) types << "5 ";
    append_data_array(cells, "UInt8", "types", 0, types.str());
  }

  // Nodal scalar fields.
  if (!point_fields.empty()) {
    pugi::xml_node pd = piece.append_child("PointData");
    pd.append_attribute("Scalars") = point_fields.front().first.c_str();
    for (const auto& [name, values] : point_fields) {
      if (static_cast<std::size_t>(values.size()) != points.size())
        throw std::runtime_error("write_vtu: field '" + name + "' size != number of points");
      std::ostringstream oss;
      oss.precision(17);
      for (Eigen::Index i = 0; i < values.size(); ++i) oss << values(i) << ' ';
      append_data_array(pd, "Float64", name.c_str(), 0, oss.str());
    }
  }

  if (!doc.save_file(path.c_str(), "  "))
    throw std::runtime_error("write_vtu: failed to write " + path);
}
}  // namespace aether::output

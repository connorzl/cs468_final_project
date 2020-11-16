#pragma once

#include "geometrycentral/numerical/linear_solvers.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/heat_method_distance.h"
#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/utilities/vector2.h"
#include "geometrycentral/utilities/vector3.h"

#include <complex>
#include <memory>
#include <tuple>
#include <vector>

namespace geometrycentral {
namespace surface {

class VectorHeatSolver {

public:
  VectorHeatSolver(IntrinsicGeometryInterface& geom);
  VertexData<double> extendScalar(const std::vector<SurfacePoint>& source_points, const std::vector<double>& source_values);
  VertexData<Vector2> transportTangentVectors(const std::vector<SurfacePoint>& source_points, const std::vector<Vector2>& source_vectors);

private:
  // Inputs.
  IntrinsicGeometryInterface& geom;
  SurfaceMesh& mesh;

  // Timestep and solvers.
  double t;
  PositiveDefiniteSolver<double>* scalarHeatFlowSolver;
  LinearSolver<std::complex<double>>* vectorHeatFlowSolver;
};


} // namespace surface
} // namespace geometrycentral

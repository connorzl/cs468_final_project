#pragma once

#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/vector_heat_method.h"

namespace geometrycentral {
namespace surface {
SurfacePoint run_weiszfeld(IntrinsicGeometryInterface& geom, VectorHeatMethodSolver& solver,
                            const std::vector<Vertex>& vertices);

Vector2 weiszfeld_energy(SurfacePoint aboutPoint, VectorHeatMethodSolver& solver, SurfaceMesh& mesh);

} // namespace surface
} // namespace geometrycentral

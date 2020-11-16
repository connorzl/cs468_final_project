#include "weiszfeld.h"

#include "geometrycentral/surface/trace_geodesic.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/vector_heat_method.h"

#include "geometrycentral/utilities/utilities.h"

namespace geometrycentral {
namespace surface {

int num_iters = 1000;
double thresh = 1 / 100.;   
double eps = 1e-6;

SurfacePoint run_weiszfeld(IntrinsicGeometryInterface& geom, VectorHeatMethodSolver& solver,
                            const std::vector<Vertex>& vertices) {
    geom.requireFaceAreas();
    geom.requireHalfedgeVectorsInVertex();
    geom.requireHalfedgeVectorsInFace();

    SurfaceMesh& mesh = geom.mesh;
    // randomly initialize starting point
    SurfacePoint result = SurfacePoint(mesh.vertex(randomIndex(mesh.nVertices())));         
    result = result.inSomeFace();
    bool converged = false;
    for (int i = 0; i < num_iters; i++) {
        Vector2 update = weiszfeld_energy(result, solver, mesh);
        result = traceGeodesic(geom, result, update).endPoint.inSomeFace();

        double faceScale = std::sqrt(geom.faceAreas[result.inSomeFace().face]);
        if (update.norm() < thresh * faceScale) {
            break;
        }
    }
    geom.unrequireFaceAreas();
    geom.unrequireHalfedgeVectorsInVertex();
    geom.unrequireHalfedgeVectorsInFace();
    return result;
}

Vector2 weiszfeld_energy(SurfacePoint point, VectorHeatMethodSolver& solver, SurfaceMesh& mesh) {
    // compute logmap
    VertexData<Vector2> logmap = solver.computeLogMap(point);

    // update step
    Vector2 update = Vector2::zero();
    double updateWeightSum = 0.;
    for (Vertex v : mesh.vertices()) {
        Vector2 pointCoord = logmap[v];
        double dist = std::sqrt(pointCoord.norm2());
        update += pointCoord / (dist + eps);
        updateWeightSum += 1.0 / (dist + eps);
    }
    update /= updateWeightSum;
    return update;
}

} // namespace surface
} // namespace geometrycentral
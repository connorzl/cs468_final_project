#include "vector_heat.h"

namespace geometrycentral {
namespace surface {

VectorHeatSolver::VectorHeatSolver(IntrinsicGeometryInterface& geom_) : geom(geom_), mesh(geom_.mesh) {
  geom.requireEdgeLengths();
  geom.requireVertexLumpedMassMatrix();
  geom.requireCotanLaplacian();
  geom.requireVertexConnectionLaplacian();

  // Set timestep as squared average edge length.
  t = 0;
  for (Edge e : mesh.edges()) {
    t += geom.edgeLengths[e];
  }
  t /= mesh.nEdges();
  t *= t;

  // Initialize solver for scalar heat flow.  
  SparseMatrix<double> LHS_scalar = geom.vertexLumpedMassMatrix + t * geom.cotanLaplacian;
  scalarHeatFlowSolver = new PositiveDefiniteSolver<double>(LHS_scalar);

  // Initialize solver for vector heat flow.
  SparseMatrix<std::complex<double>> LHS_vector = geom.vertexLumpedMassMatrix.cast<std::complex<double>>() + t * geom.vertexConnectionLaplacian;
  vectorHeatFlowSolver = new SquareSolver<std::complex<double>>(LHS_vector);

  geom.unrequireEdgeLengths();
  geom.unrequireCotanLaplacian();
  geom.unrequireVertexConnectionLaplacian();
  geom.unrequireVertexLumpedMassMatrix();
}

VertexData<double> VectorHeatSolver::extendScalar(const std::vector<SurfacePoint>& source_points, const std::vector<double>& source_values) {
  geom.requireVertexIndices();
  Vector<double> u_0 = Vector<double>::Zero(mesh.nVertices());
  Vector<double> phi_0 = Vector<double>::Zero(mesh.nVertices());
  for (size_t i = 0; i < source_points.size(); i++) {
    SurfacePoint facePoint = source_points[i].inSomeFace();
    Vector3 baryCoords = facePoint.faceCoords;
    Halfedge he = facePoint.face.halfedge();
    for (size_t j = 0; j < 3; j++) {
      size_t k = geom.vertexIndices[he.vertex()];
      u_0[k] += baryCoords[j] * source_values[i];
      phi_0[k] += baryCoords[j];
      he = he.next();
    }
  }
  Vector<double> u = scalarHeatFlowSolver->solve(u_0);
  Vector<double> phi = scalarHeatFlowSolver->solve(phi_0);
  VertexData<double> result(mesh, u.array() / phi.array());
  geom.unrequireVertexIndices();
  return result;
}

VertexData<Vector2> VectorHeatSolver::transportTangentVectors(const std::vector<SurfacePoint>& source_points, const std::vector<Vector2>& source_vectors) {
  geom.requireVertexIndices();
  Vector<std::complex<double>> Y_0 = Vector<std::complex<double>>::Zero(mesh.nVertices());

  // Compute directions.
  for (size_t i = 0; i < source_points.size(); i++) {
    std::complex<double> v = Vector2::fromComplex(source_vectors[i]).normalize();
    SurfacePoint facePoint = source_points[i].inSomeFace();
    Vector3 baryCoords = facePoint.faceCoords;
    Halfedge he = facePoint.face.halfedge();
    for (size_t j = 0; j < 3; j++) {
      size_t k = geom.vertexIndices[he.vertex()];
      Y_0[k] += baryCoords[j] * v;
      he = he.next();
    }
  }
  Vector<std::complex<double>> Y = vectorHeatFlowSolver->solve(Y_0);

  // Compute magnitudes.
  std::vector<double> source_magnitudes;
  for (size_t i = 0; i < source_vectors.size(); i++) {
    source_magnitudes.push_back(source_vectors[i].norm());
  }
  VertexData<double> r = extendScalar(source_points, source_magnitudes);

  VertexData<Vector2> result(mesh);
  for (Vertex v : mesh.vertices()) {
    Vector2 xy = Vector2::fromComplex(Y[geom.vertexIndices[v]]).normalize();
    result[v] = r[v] * xy;
  }
  geom.unrequireVertexIndices();
  return result;
}

} // namespace surface
} // namespace geometrycentral
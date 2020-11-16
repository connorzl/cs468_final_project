#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "vector_heat.h"
#include "weiszfeld.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/vector_heat_method.h"

#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

#include <sstream>

using namespace geometrycentral;
using namespace geometrycentral::surface;

std::unique_ptr<HalfedgeMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;

// Sources.
std::vector<Vertex> source_vertices;
std::vector<double> source_values;
std::vector<Vector2> source_vectors;

polyscope::SurfaceMesh* psMesh;
std::unique_ptr<VectorHeatSolver> solver;

void addSource(size_t v_ind, double val, double r, double phi) {
  Vertex v = mesh->vertex(v_ind);
  source_vertices.push_back(v);
  source_values.push_back(val);
  source_vectors.push_back(r * Vector2::fromAngle(phi));
}

void drawSources() {
  // Draw each source point and source vector.
  std::vector<std::pair<size_t, double>> sourcePairs;
  VertexData<Vector2> sourceVectors(*mesh, Vector2::zero());
  for (size_t i = 0; i < source_vertices.size(); i++) {
    size_t ind = geometry->vertexIndices[source_vertices[i]];
    sourcePairs.emplace_back(ind, source_values[i]);
    sourceVectors[source_vertices[i]] = source_vectors[i];
  }

  // Point at sources.
  auto scalarQ = polyscope::getSurfaceMesh()->addVertexIsolatedScalarQuantity("source scalars", sourcePairs);
  scalarQ->pointRadius *= 5;
  scalarQ->cMap = "blues";
  scalarQ->setEnabled(true);

  // Vector at sources.
  auto vectorQ = polyscope::getSurfaceMesh()->addVertexIntrinsicVectorQuantity("source vectors", sourceVectors);
  vectorQ->setVectorLengthScale(0.07);
  vectorQ->setVectorRadius(0.01);
  vectorQ->setVectorColor(glm::vec3{206 / 255., 129 / 255., 247 / 255.});
  vectorQ->setEnabled(true);
}

int main(int argc, char** argv) {
  args::ArgumentParser parser("Implementation of Vector Heat Method.");
  args::Positional<std::string> inputFilename(parser, "mesh", "Input mesh file.");
  try {
    parser.ParseCLI(argc, argv);
  } catch (args::Help&) {
    std::cout << parser;
    return 0;
  } catch (args::ParseError& e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return 1;
  }

  polyscope::init();
  std::tie(mesh, geometry) = loadMesh(args::get(inputFilename));
  psMesh = polyscope::registerSurfaceMesh(polyscope::guessNiceNameFromPath(args::get(inputFilename)),
                                          geometry->inputVertexPositions, mesh->getFaceVertexList(),
                                          polyscopePermutations(*mesh));

  // Compute tangent spaces.
  geometry->requireVertexTangentBasis();
  geometry->requireFaceTangentBasis();
  VertexData<Vector3> vBasisX(*mesh);
  FaceData<Vector3> fBasisX(*mesh);
  for (Vertex v : mesh->vertices()) {
    vBasisX[v] = geometry->vertexTangentBasis[v][0];
  }
  for (Face f : mesh->faces()) {
    fBasisX[f] = geometry->faceTangentBasis[f][0];
  }
  polyscope::getSurfaceMesh()->setVertexTangentBasisX(vBasisX);
  polyscope::getSurfaceMesh()->setFaceTangentBasisX(fBasisX);

  // To start, pick two vertices as sources
  geometry->requireVertexIndices();
  addSource(0, 1.0, 1.0, 0);
  addSource(mesh->nVertices() / 2, 3.0, 1.0, 0);

  /*
  addSource(12214, 1.0, 1.0, 0);
  addSource(2197, 2.0, 1.0, 0);
  addSource(28899, 3.0, 1.0, 0);
  addSource(37742, 4.0, 1.0, 0);
  */

  drawSources();
  solver.reset(new VectorHeatSolver(*geometry));

  // Scalar flow.
  std::vector<SurfacePoint> source_points;
  for (size_t i = 0; i < source_vertices.size(); i++) {
    source_points.push_back(source_vertices[i]);
  }
  VertexData<double> scalarFlow = solver->extendScalar(source_points, source_values);
  auto psScalar = psMesh->addVertexScalarQuantity("Flow", scalarFlow);
  psScalar->setEnabled(true);

  // Vector flow.
  VertexData<Vector2> vectorFlow = solver->transportTangentVectors(source_points, source_vectors);
  auto psVec = psMesh->addVertexIntrinsicVectorQuantity("Vector Flow", vectorFlow);
  psVec->setEnabled(true);

  // Geometric medians.
  std::unique_ptr<VectorHeatMethodSolver> gt_solver;
  gt_solver.reset(new VectorHeatMethodSolver(*geometry));
  SurfacePoint center = run_weiszfeld(*geometry, *gt_solver, source_vertices);
  std::vector<Vector3> point{center.interpolate(geometry->inputVertexPositions)};
  auto pointQ = polyscope::registerPointCloud("center", point);
  pointQ->setPointRadius(5.);

  polyscope::show();
  return EXIT_SUCCESS;
}
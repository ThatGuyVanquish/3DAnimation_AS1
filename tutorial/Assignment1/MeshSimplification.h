#pragma once

#include <utility>
#include "IglMeshLoader.h"
#include "igl_inline.h"
#include <igl/circulation.h>
#include <igl/collapse_edge.h>
#include <igl/edge_flaps.h>
#include <igl/parallel_for.h>
#include <igl/read_triangle_mesh.h>
#include <Eigen/Core>
#include <iostream>
#include "../engine/Mesh.h"
#include <Eigen/LU> 
#include <map>
#include <igl/per_face_normals.h>
#include <numeric>

/*

    Mesh Simplification object holds a simplified mesh object
    with

*/
class MeshSimplification
{
public:
    MeshSimplification(std::string filename, int _decimations);
    std::shared_ptr<cg3d::Mesh> getMesh();

private:
    std::shared_ptr<cg3d::Mesh> currentMesh;
    int decimations;
    int collapseCounter;
    int QResetInterval;
    std::set<int> verticesToUpdate;
    Eigen::MatrixXd faceNormals;
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::MatrixXi E;
    Eigen::VectorXi EMAP;
    Eigen::MatrixXi EF;
    Eigen::MatrixXi EI;
    igl::min_heap<std::tuple<double, int, int>> Q;
    Eigen::VectorXi EQ;
    Eigen::MatrixXd C;
    std::vector<Eigen::Matrix4d> verticesToQ;
    std::map<int, std::set<int>> verticesToFaces;
    
    // Helper functions
    Eigen::Vector4d ThreeDimVecToFourDim(Eigen::Vector3d vertex);
    Eigen::Vector3d FourDimVecToThreeDim(Eigen::Vector4d vertex);

    // Methods to calculate Q matrices
    Eigen::Vector4d calculatePlaneNormal(int face);
    Eigen::Vector4d equation_plane(Eigen::Vector3i triangle);
    Eigen::Matrix4d calculateKp(Eigen::Vector4d planeVector);
    double calculateCost(const int vertex);
    void buildVerticesToFaces();
    void calculateQ(int v);
    void calculateQs(std::vector<int> verticesToCalculate);

    // Edge collapsing methods
    void updateVerticesToFacesPostCollapse(const int e);
    void post_collapse(const int e);
    Eigen::Matrix4d calculateQDerive(Eigen::Matrix4d currentQ);
    void quadratic_error_simplification(const int e, double& cost, Eigen::RowVectorXd& p);
    bool collapse_edge();
    
    // Run methods
    void Init();
    void createDecimatedMesh();
    
    
};
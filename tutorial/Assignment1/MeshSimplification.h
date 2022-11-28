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
    std::map<int, std::vector<int>> verticesToFaces;

    void Init();
    void createDecimatedMesh();
    Eigen::Vector4d equation_plane(Eigen::Vector3i triangle, Eigen::MatrixXd& V);
    Eigen::Matrix4d calculateKp(Eigen::Vector4d planeVector);
    void calculateQs();
    double calculateCost(Eigen::Matrix4d QMatrix, Eigen::Vector3d vertex);
    Eigen::Matrix4d calculateQDerive(Eigen::Matrix4d currentQ);
    void quadratic_error_simplification(const int e, double& cost, Eigen::RowVectorXd& p);
    bool collapse_edge();
};
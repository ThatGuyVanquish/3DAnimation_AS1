#pragma once

#include "Scene.h"

#include "AutoMorphingModel.h"
#include <utility>
#include "ObjLoader.h"
#include "IglMeshLoader.h"
#include "igl_inline.h"
#include <igl/circulation.h>
#include <igl/collapse_edge.h>
#include <igl/edge_flaps.h>
#include <igl/decimate.h>
#include <igl/shortest_edge_and_midpoint.h>
#include <igl/parallel_for.h>
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>
#include <iostream>
#include <set>
#include "../engine/Mesh.h"
#include <Eigen/LU> 

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
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::MatrixXi E;
    Eigen::VectorXi EMAP;
    Eigen::MatrixXi EF;
    Eigen::MatrixXi EI;
    igl::min_heap<std::tuple<double, int>> Q;
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
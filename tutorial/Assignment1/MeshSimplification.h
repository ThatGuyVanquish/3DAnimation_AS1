#pragma once

#include <utility>
#include "IglMeshLoader.h"
#include "igl_inline.h"
#include <igl/collapse_edge.h>
#include <igl/edge_flaps.h>
#include <igl/parallel_for.h>
#include <igl/read_triangle_mesh.h>
#include <igl/circulation.h>
#include <Eigen/Core>
#include <iostream>
#include "../engine/Mesh.h"
#include <Eigen/LU> 
#include <map>
#include <igl/per_face_normals.h>
#include <igl/vertex_triangle_adjacency.h>
#include <numeric>
#include <cmath>

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
    
    // HELPER METHODS
    Eigen::Vector4d ThreeDimVecToFourDim(Eigen::Vector3d vertex);
    Eigen::Vector3d FourDimVecToThreeDim(Eigen::Vector4d vertex);

    // METHODS TO CALCULATE Q MATRICES
    Eigen::Vector4d planeEquation(Eigen::Vector3i triangle, const Eigen::MatrixXd& V);
    Eigen::Vector4d calculatePlaneNormal(const Eigen::MatrixXd& V, Eigen::Vector3d threeDimNormal, int vi);
    Eigen::Matrix4d calculateKp(Eigen::Vector4d planeVector);
    void calculateQs(const Eigen::MatrixXd& V, 
        const Eigen::MatrixXi &F, 
        Eigen::MatrixXd& faceNormals,
        std::vector<std::vector<int>> &verticesToFaces,
        std::vector<std::vector<int>> &facesBeforeIndex,
        std::vector<Eigen::Matrix4d> &verticesToQ);

    // METHODS TO CALCULATE COST AND PLACEMENT
    Eigen::Matrix4d calculateQDerive(Eigen::Matrix4d currentQ);
    double calculateCost(Eigen::Vector4d vertex, Eigen::Matrix4d Q);
    void preCalculateCostAndPos(
        const igl::decimate_cost_and_placement_callback& callback,
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& F,
        const Eigen::MatrixXi& E,
        const Eigen::VectorXi& EMAP,
        const Eigen::MatrixXi& EF,
        const Eigen::MatrixXi& EI,
        igl::min_heap<std::tuple<double, int, int>>& Q,
        Eigen::VectorXi& EQ,
        Eigen::MatrixXd& C
    );
    
    // RUN METHODS
    void createDecimatedMesh(std::string fileName);
};
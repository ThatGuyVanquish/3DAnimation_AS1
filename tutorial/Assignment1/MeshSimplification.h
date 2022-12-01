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
    int QResetInterval;
    
    // Helper functions
    Eigen::Vector4d ThreeDimVecToFourDim(Eigen::Vector3d vertex);
    Eigen::Vector3d FourDimVecToThreeDim(Eigen::Vector4d vertex);

    // Methods to calculate Q matrices
    Eigen::Vector4d calculatePlaneNormal(const Eigen::MatrixXd& V, Eigen::Vector3d threeDimNormal, int vi);
    Eigen::Matrix4d calculateKp(Eigen::Vector4d planeVector);

    // Edge collapsing methods
    Eigen::Matrix4d calculateQDerive(Eigen::Matrix4d currentQ);
    double calculateCost(Eigen::Vector4d vertex, Eigen::Matrix4d Q);
    
    // Run methods
    void createDecimatedMesh(std::string fileName);
};
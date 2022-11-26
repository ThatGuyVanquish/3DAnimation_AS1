#pragma once

#include "Scene.h"

#include "AutoMorphingModel.h"
#include <utility>
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
#include "eigen/Eigen/src/Core/PlainObjectBase.h"
#include "igl/per_vertex_normals.h"
#include "MeshSimplification.h"

class BasicScene : public cg3d::Scene
{
public:
    explicit BasicScene(std::string name, cg3d::Display* display) : Scene(std::move(name), display) {};
    void Init(float fov, int width, int height, float near, float far);
    void Update(const cg3d::Program& program, const Eigen::Matrix4f& proj, const Eigen::Matrix4f& view, const Eigen::Matrix4f& model) override;
    void KeyCallback(cg3d::Viewport* viewport, int x, int y, int key, int scancode, int action, int mods) override;

private:
    std::shared_ptr<Movable> root;
    std::shared_ptr<cg3d::Model> camel;
    std::shared_ptr<cg3d::AutoMorphingModel> autoCamel;
    std::shared_ptr<MeshSimplification> simplificationObject;
};

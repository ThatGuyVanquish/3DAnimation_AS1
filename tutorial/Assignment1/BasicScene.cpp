#include "BasicScene.h"
#include <utility>
#include "ObjLoader.h"
#include "IglMeshLoader.h"

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

// #include "AutoMorphingModel.h"

using namespace cg3d;

void BasicScene::Init(float fov, int width, int height, float near, float far)
{
    using namespace std;
    using namespace Eigen;
    using namespace igl;
    camera = Camera::Create( "camera", fov, float(width) / height, near, far);
    
    AddChild(root = Movable::Create("root")); // a common (invisible) parent object for all the shapes
    auto daylight{std::make_shared<Material>("daylight", "shaders/cubemapShader")}; 
    daylight->AddTexture(0, "textures/cubemaps/Daylight Box_", 3);
    auto background{Model::Create("background", Mesh::Cube(), daylight)};
    AddChild(background);
    background->Scale(120, Axis::XYZ);
    background->SetPickable(false);
    background->SetStatic();

    auto program = std::make_shared<Program>("shaders/basicShader");
    auto material{ std::make_shared<Material>("material", program)}; // empty material
 
    material->AddTexture(0, "textures/box0.bmp", 2);
    auto camelMesh{IglLoader::MeshFromFiles("cyl_igl","data/camel_b.obj")};
    camel = Model::Create( "camel", camelMesh, material);
    camel->Translate({3,0,0});
    camel->Scale(0.12f);
    camel->showWireframe = true;
    camera->Translate(20, Axis::Z);
    root->AddChild(camel);
    
    auto mesh = camel->GetMeshList();

    string filename("data/camel_b.obj");

    MatrixXd V, OV;
    MatrixXi F, OF;
    read_triangle_mesh(filename, OV, OF);

    // Prepare array-based edge data structures and priority queue
    VectorXi EMAP;
    MatrixXi E, EF, EI;
    min_heap< std::tuple<double, int, int> > Q;
    VectorXi EQ;
    // If an edge were collapsed, we'd collapse it to these points:
    MatrixXd C;
    int num_collapsed;

    // Function to reset original mesh and data structures
    const auto& reset = [&]()
    {
        F = OF;
        V = OV;
        edge_flaps(F, E, EMAP, EF, EI);
        C.resize(E.rows(), V.cols());
        VectorXd costs(E.rows());
        Q = {};
        EQ = VectorXi::Zero(E.rows());
        {
            VectorXd costs(E.rows());
            parallel_for(E.rows(), [&](const int e)
                {
                    double cost = e;
                    RowVectorXd p(1, 3);
                    shortest_edge_and_midpoint(e, V, F, E, EMAP, EF, EI, cost, p);
                    C.row(e) = p;
                    costs(e) = cost;
                }, 10000);
            for (int e = 0; e < E.rows(); e++)
            {
                Q.emplace(costs(e), e, 0);
            }
        }
        num_collapsed = 0;
    };
    
}

void BasicScene::Update(const Program& program, const Eigen::Matrix4f& proj, const Eigen::Matrix4f& view, const Eigen::Matrix4f& model)
{
    Scene::Update(program, proj, view, model);
    program.SetUniform4f("lightColor", 1.0f, 1.0f, 1.0f, 0.5f);
    program.SetUniform4f("Kai", 1.0f, 1.0f, 1.0f, 1.0f);
    camel->Rotate(0.01f, Axis::Y);
}

void BasicScene::KeyCallback(Viewport* viewport, int x, int y, int key, int scancode, int action, int mods)
{
    auto system = camera->GetRotation().transpose();

    if (action == GLFW_PRESS || action == GLFW_REPEAT) {
        switch (key) // NOLINT(hicpp-multiway-paths-covered)
        {
        case GLFW_KEY_ESCAPE:
            glfwSetWindowShouldClose(window, GLFW_TRUE);
            break;
        case GLFW_KEY_UP:
            camera->RotateInSystem(system, 0.1f, Axis::X);
            break;
        case GLFW_KEY_DOWN:
            camera->RotateInSystem(system, -0.1f, Axis::X);
            break;
        case GLFW_KEY_LEFT:
            camera->RotateInSystem(system, 0.1f, Axis::Y);
            break;
        case GLFW_KEY_RIGHT:
            camera->RotateInSystem(system, -0.1f, Axis::Y);
            break;
        case GLFW_KEY_W:
            camera->TranslateInSystem(system, { 0, 0.05f, 0 });
            break;
        case GLFW_KEY_S:
            camera->TranslateInSystem(system, { 0, -0.05f, 0 });
            break;
        case GLFW_KEY_A:
            camera->TranslateInSystem(system, { -0.05f, 0, 0 });
            break;
        case GLFW_KEY_D:
            camera->TranslateInSystem(system, { 0.05f, 0, 0 });
            break;
        case GLFW_KEY_B:
            camera->TranslateInSystem(system, { 0, 0, 0.05f });
            break;
        case GLFW_KEY_F:
            camera->TranslateInSystem(system, { 0, 0, -0.05f });
            break;
        case GLFW_KEY_SPACE:
            /* logic to decimate*/
        }
    }
}

void BasicScene::Simpilifcation() 
{

}
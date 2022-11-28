#include "BasicScene.h"

using namespace cg3d;

std::shared_ptr<cg3d::Mesh> BasicScene::createDecimatedMesh(std::string filename, int decimations)
{
    auto currentMesh{ IglLoader::MeshFromFiles("currentMesh", filename) };
    Eigen::MatrixXd V, OV;
    Eigen::MatrixXi F, OF;
    OV = currentMesh->data[0].vertices;
    OF = currentMesh->data[0].faces;

    // Prepare array-based edge data structures and priority queue
    Eigen::VectorXi EMAP;
    Eigen::MatrixXi E, EF, EI;
    igl::min_heap< std::tuple<double, int, int>> Q;
    Eigen::VectorXi EQ;
    // If an edge were collapsed, we'd collapse it to these points:
    Eigen::MatrixXd C;

    // Function to reset original mesh and data structures
    const auto& reset = [&]()
    {
        F = OF;
        V = OV;
        igl::edge_flaps(F, E, EMAP, EF, EI);
        C.resize(E.rows(), V.cols());
        Q = {};
        EQ = Eigen::VectorXi::Zero(E.rows());
        {
            Eigen::VectorXd costs(E.rows());
            igl::parallel_for(E.rows(), [&](const int e)
                {
                    double cost = e;
                    Eigen::RowVectorXd p(1, 3);
                    igl::shortest_edge_and_midpoint(e, V, F, E, EMAP, EF, EI, cost, p);
                    C.row(e) = p;
                    costs(e) = cost;
                }, 10000);
            for (int e = 0; e < E.rows(); e++)
            {
                Q.emplace(costs(e), e, 0);
            }
        }
    };

    reset();

    for (int i = 0; i < decimations; i++)
    {
        int collapsed_edges = 0;
        const int max_iter = std::ceil(0.01 * Q.size());
        for (int j = 0; j < max_iter; j++)
        {
            if (!igl::collapse_edge(igl::shortest_edge_and_midpoint, V, F, E, EMAP, EF, EI, Q, EQ, C))
            {
                break;
            }
            collapsed_edges++;
        }
        if (collapsed_edges > 0) {
            currentMesh->data.push_back({ V, F, currentMesh->data[0].vertexNormals, currentMesh->data[0].textureCoords });
        }
    }
    return currentMesh;
}

void BasicScene::Init(float fov, int width, int height, float near, float far)
{
    using namespace std;
    using namespace Eigen;
    using namespace igl;
    camera = Camera::Create("camera", fov, float(width) / height, near, far);

    AddChild(root = Movable::Create("root")); // a common (invisible) parent object for all the shapes
    auto daylight{ std::make_shared<Material>("daylight", "shaders/cubemapShader") };
    daylight->AddTexture(0, "textures/cubemaps/Daylight Box_", 3);
    auto background{ Model::Create("background", Mesh::Cube(), daylight) };
    AddChild(background);
    background->Scale(120, Axis::XYZ);
    background->SetPickable(false);
    background->SetStatic();

    auto program = std::make_shared<Program>("shaders/basicShader");
    auto material{ std::make_shared<Material>("material", program) }; // empty material
    material->AddTexture(0, "textures/box0.bmp", 2);

    std::string objFile = "data/cube.off";
    int decimations = 6;

    auto morphFunc = [](Model* model, cg3d::Visitor* visitor) {
        return model->meshIndex;
    };

    myAutoModel = cg3d::AutoMorphingModel::Create(
        *Model::Create("myModel", createDecimatedMesh(objFile, decimations), material),
        morphFunc
    );
    myAutoModel->Translate({ 1,-3,0 });
    myAutoModel->Scale(3.0f);
    myAutoModel->showWireframe = true;
    root->AddChild(myAutoModel);
    camera->Translate(10, Axis::Z);
}

void BasicScene::Update(const Program& program, const Eigen::Matrix4f& proj, const Eigen::Matrix4f& view, const Eigen::Matrix4f& model)
{
    Scene::Update(program, proj, view, model);
    program.SetUniform4f("lightColor", 1.0f, 1.0f, 1.0f, 0.5f);
    program.SetUniform4f("Kai", 1.0f, 1.0f, 1.0f, 1.0f);
    myAutoModel->Rotate(0.001f, Axis::Y);
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
            if (myAutoModel->meshIndex > 0)
                myAutoModel->meshIndex--;
            break;
        case GLFW_KEY_DOWN:
            if (myAutoModel->meshIndex < myAutoModel->GetMesh(0)->data.size())
                myAutoModel->meshIndex++;
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
        }
    }
}
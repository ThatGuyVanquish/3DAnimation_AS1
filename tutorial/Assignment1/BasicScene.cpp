#include "BasicScene.h"

using namespace cg3d;

enum keyPress { UP, DOWN, UNPRESSED };
keyPress lastKey = UNPRESSED;
int lastState = 0;

std::vector<MeshData> BasicScene::createDecimatedMesh(std::string filename)
{
    auto currentMesh{ IglLoader::MeshFromFiles("cyl_igl", filename) };
    Eigen::MatrixXd V, OV;
    Eigen::MatrixXi F, OF;

    igl::read_triangle_mesh(filename, OV, OF);

    // Prepare array-based edge data structures and priority queue
    Eigen::VectorXi EMAP;
    Eigen::MatrixXi E, EF, EI;
    igl::min_heap< std::tuple<double, int, int>> Q;
    Eigen::VectorXi EQ;
    // If an edge were collapsed, we'd collapse it to these points:
    Eigen::MatrixXd C;
    int num_collapsed;

    // Function to reset original mesh and data structures
    const auto& reset = [&]()
    {
        F = OF;
        V = OV;
        igl::edge_flaps(F, E, EMAP, EF, EI);
        C.resize(E.rows(), V.cols());
        Eigen::VectorXd costs(E.rows());
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
        num_collapsed = 0;
    };

    reset();

    std::vector<cg3d::MeshData> meshDataVector;
    meshDataVector.push_back(currentMesh->data[0]);

    for (int i = 0; i < 6; i++)
    {
        const int max_iter = std::ceil(0.01 * Q.size());
        for (int j = 0; j < max_iter; j++)
        {
            if (!igl::collapse_edge(igl::shortest_edge_and_midpoint, V, F, E, EMAP, EF, EI, Q, EQ, C))
            {
                break;
            }
        }
        std::cout << "this is a test\n";
        cg3d::Mesh temp("new mesh", V, F, currentMesh->data[0].vertexNormals, currentMesh->data[0].textureCoords);
        meshDataVector.push_back(temp.data[0]);
    }
    std::cout << "Decimated mesh size is: " << meshDataVector.size() << std::endl;
    return meshDataVector;
}

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

    std::shared_ptr<cg3d::Mesh> camelMesh(new cg3d::Mesh(std::string("camelWithDecimations"), createDecimatedMesh("data/bunny.off")));
    camel = Model::Create("camel", camelMesh, material);

    auto morphFunc = [](Model* model, cg3d::Visitor* visitor) {
        if (lastKey == UP && lastState > 0)
        {
            lastKey = UNPRESSED;
            return --lastState;
        }
            
        if (lastKey == DOWN && lastState < 6)
        {
            lastKey = UNPRESSED;
            return ++lastState;
        }

        lastKey = UNPRESSED;
        return lastState;
    };

    autoCamel = cg3d::AutoMorphingModel::Create(*camel, morphFunc);
    autoCamel->Translate({ 1,-3,0 });
    autoCamel->Scale(30.0f);
    autoCamel->showWireframe = true;
    root->AddChild(autoCamel);
    camera->Translate(10, Axis::Z);
}

void BasicScene::Update(const Program& program, const Eigen::Matrix4f& proj, const Eigen::Matrix4f& view, const Eigen::Matrix4f& model)
{
    Scene::Update(program, proj, view, model);
    program.SetUniform4f("lightColor", 1.0f, 1.0f, 1.0f, 0.5f);
    program.SetUniform4f("Kai", 1.0f, 1.0f, 1.0f, 1.0f);
    autoCamel->Rotate(0.001f, Axis::Y);
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
            lastKey = UP;
            break;
        case GLFW_KEY_DOWN:
            lastKey = DOWN;
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
        default:
            lastKey = UNPRESSED;
            break;
        }
    }
}

Eigen::Vector4d BasicScene::equation_plane(std::vector<int> triangle, Eigen::MatrixXd& V)
{
    auto x2MinusX1 = V.row(triangle[1]) - V.row(triangle[0]);
    double a1 = x2MinusX1[0];
    double b1 = x2MinusX1[1];
    double c1 = x2MinusX1[2];
    auto x3MinusX1 = V.row(triangle[2]) - V.row(triangle[0]);
    double a2 = x3MinusX1[0];
    double b2 = x3MinusX1[1];
    double c2 = x3MinusX1[2];
    double a = b1 * c2 - b2 * c1;
    double b = a2 * c1 - a1 * c2;
    double c = a1 * b2 - b1 * a2;
    double d = (-a * V.row(triangle[0])[0] - b * V.row(triangle[0])[1] - c * V.row(triangle[0])[2]);
    double normalizer = sqrt(pow(a, 2) + pow(b, 2) + pow(c, 2));
    Eigen::Vector4d ret;
    ret(0) = a / normalizer, ret(1) = b / normalizer, ret(2) = c / normalizer, ret(3) = d / normalizer;
    return ret;
}

Eigen::Matrix4d BasicScene::calculateQ(Eigen::MatrixXd planeMatrix)
{
    /* each row of planeMatrix is a plane equation*/
    Eigen::Matrix4d ret;
    ret.Zero();
    for (int i = 0; i < planeMatrix.rows(); i++)
    {
        ret = ret + planeMatrix.row(i).transpose() * planeMatrix.row(i);
    }
    return ret;
}

/*
    some function that runs equation_plane in a loop for every face and adds the last vector to a matrix
*/

double BasicScene::calculateCost(Eigen::Matrix4d QMatrix, Eigen::Vector4d vertex)
{
    return vertex.transpose() * QMatrix * vertex;
}

#include "BasicScene.h"
#include <Eigen/LU> 

using namespace cg3d;

void BasicScene::Init(float fov, int width, int height, float near, float far)
{
    //using namespace igl;
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
    std::shared_ptr<cg3d::Mesh> camelMesh(new cg3d::Mesh(std::string("camelWithDecimations"), simplificationObject->MeshSimplification::createDecimatedMesh("data/cube.off")));
    camel = Model::Create("camel", camelMesh, material);

    auto morphFunc = [](Model* model, cg3d::Visitor* visitor) {
        return model->meshIndex;
    };

    autoCamel = cg3d::AutoMorphingModel::Create(*camel, morphFunc);
    //autoCamel->Translate({ 1,-3,0 });
    autoCamel->Translate({ 1,0,0 });
    //autoCamel->Scale(30.0f);
    autoCamel->Scale(1.0f);
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
            if (autoCamel->meshIndex > 0)
                autoCamel->meshIndex--;
            break;
        case GLFW_KEY_DOWN:
            if (autoCamel->meshIndex < 10)
                autoCamel->meshIndex++;
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
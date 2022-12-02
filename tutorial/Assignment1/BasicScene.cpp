#include "BasicScene.h"

using namespace cg3d;

void BasicScene::Init(float fov, int width, int height, float near, float far)
{
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
    std::vector<std::string> objFiles{ "data/bunny.off", /* 0 */
        "data/sphere.obj", /* 1 */
        "data/cheburashka.off", /* 2 */
        "data/fertility.off" /* 3 */};
    int objIndex = 3;
    int decimations = 10;
    myMeshObj = std::make_shared<MeshSimplification>(MeshSimplification(objFiles[objIndex], decimations));
    auto morphFunc = [](Model* model, cg3d::Visitor* visitor) {
        return model->meshIndex;
    };

    myAutoModel = cg3d::AutoMorphingModel::Create(
        *cg3d::Model::Create("My Model", myMeshObj->getMesh(), material),
        morphFunc
    );
    float cameraTranslate = 0;
    
    switch (objIndex)
    {
    case 0: /* Bunny */
        myAutoModel->Scale(40.0f);
        myAutoModel->Translate({ 1,-4,0 });
        cameraTranslate = 10;
        break;
    case 1: /* Sphere */
        myAutoModel->Scale(5.0f);
        myAutoModel->Translate({ 0,0,0 });
        cameraTranslate = 10;
        break;
    case 2: /* Cheburashka */
        myAutoModel->Scale(25.0f);
        myAutoModel->Translate({ -12.45,-12,0 });
        cameraTranslate = 40;
        break;
    case 3: /* Fertility */
        myAutoModel->Scale(0.156f);
        myAutoModel->Translate({ -3.5,0,0 });
        cameraTranslate = 40;
        break;
    default: 
        myAutoModel->Scale(3.0f);
        myAutoModel->Translate({ 0,0,0 });
        cameraTranslate = 10;
    }
    myAutoModel->showWireframe = true;
    root->AddChild(myAutoModel);
    camera->Translate(cameraTranslate, Axis::Z);
}

void BasicScene::Update(const Program& program, const Eigen::Matrix4f& proj, const Eigen::Matrix4f& view, const Eigen::Matrix4f& model)
{
    Scene::Update(program, proj, view, model);
    program.SetUniform4f("lightColor", 1.0f, 1.0f, 1.0f, 0.5f);
    program.SetUniform4f("Kai", 1.0f, 1.0f, 1.0f, 1.0f);
    //myAutoModel->Rotate(0.01f, Axis::Y);
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
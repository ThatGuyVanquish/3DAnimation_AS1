#include "MeshSimplification.h"

MeshSimplification::MeshSimplification(std::string filename):
    currentMesh(cg3d::IglLoader::MeshFromFiles("cyl_igl", filename))
{
    V = currentMesh->data[0].vertices, F = currentMesh->data[0].faces;
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
                // verticesToQ = calculateQs(F, V);
                    //quadratic_error_simplification(e, V, F, E, EMAP, EF, EI, cost, p);
                C.row(e) = p;
                costs(e) = cost;
            }, 10000);
        for (int e = 0; e < E.rows(); e++)
        {
            Q.emplace(costs(e), e);
        }
    }
    Init();
    calculateQs();
}

void MeshSimplification::Init()
{
    for (int i = 0; i < F.rows(); i++)
    {
        Eigen::VectorXi currentRowInF = F.row(i);
        // prepare vertex to face map
        verticesToFaces.clear();
        for (int j = 0; j < currentRowInF.cols(); j++)
        {
            int vertexId = currentRowInF(j);
            verticesToFaces[vertexId].push_back(i);
        }
    }
}

Eigen::Vector4d MeshSimplification::equation_plane(Eigen::Vector3i triangle, Eigen::MatrixXd& V)
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

Eigen::Matrix4d MeshSimplification::calculateKp(Eigen::Vector4d planeVector)
{
    std::cout << planeVector * planeVector.transpose() << std::endl;
    return planeVector * planeVector.transpose();
}

double MeshSimplification::calculateCost(Eigen::Matrix4d QMatrix, Eigen::Vector3d vertex)
{
    Eigen::Vector4d convertedPoint;
    convertedPoint[0] = vertex[0], convertedPoint[1] = vertex[1], convertedPoint[2] = vertex[2], convertedPoint[3] = 1;
    return convertedPoint.transpose() * QMatrix * convertedPoint;
}

void MeshSimplification::calculateQs()
{
    std::vector <Eigen::Matrix4d> QforVertex(V.rows());
    for (Eigen::Matrix4d currentMatrix : QforVertex) // set every matrix to the zero matrix
        currentMatrix = Eigen::Matrix4d::Zero();

    for (int i = 0; i < F.rows(); i++)
    {
        Eigen::VectorXi currentRowInF = F.row(i);
        auto currentPlaneVector = equation_plane(currentRowInF, V);
        auto KpForPlaneI = calculateKp(currentPlaneVector);

        for (int j = 0; j < 3; j++)
        {
            int currentVertex = currentRowInF[j];
            QforVertex[currentVertex] += KpForPlaneI;
        }
    }
    verticesToQ = QforVertex;
}

Eigen::Matrix4d MeshSimplification::calculateQDerive(Eigen::Matrix4d currentQ)
{
    auto lastRow = currentQ.row(3);
    lastRow[0] = 0, lastRow[1] = 0, lastRow[2] = 0, lastRow[3] = 1;
    return currentQ;
}

IGL_INLINE void MeshSimplification::quadratic_error_simplification(
    const int e,
    double& cost,
    Eigen::RowVectorXd& p)
{
    // Calculate p
    int vertex1 = E(e, 0), vertex2 = E(e, 1);
    Eigen::Matrix4d Q1 = verticesToQ[vertex1], Q2 = verticesToQ[vertex2], Qtag = Q1 + Q2;
    Eigen::Matrix4d derivedQtag = calculateQDerive(Qtag);
    Eigen::Matrix4d derivedQtagInverse;
    bool invertible = false;
    double determinant;
    derivedQtag.computeInverseAndDetWithCheck(derivedQtagInverse, determinant, invertible);
    Eigen::Vector4d vtag;
    Eigen::Vector4d v0001;
    v0001[0] = 0, v0001[1] = 0, v0001[2] = 0, v0001[3] = 1;
    if (!invertible)
    {
        Eigen::Vector4d v1 = V.row(vertex1), v2 = V.row(vertex2), v3 = (v1 + v2) / 2;
        double costV1 = v1.transpose() * Qtag * v1, costV2 = v2.transpose() * Qtag * v2, costV3 = v3.transpose() * Qtag * v3;
        vtag = costV1 <= costV2 ?
            (costV1 <= costV3 ? v1 : v3) :
            costV2 <= costV3 ? v2 : v3;
    }
    else
    {
        vtag = derivedQtagInverse * v0001;
    }
    p = vtag;
    cost = vtag.transpose() * Qtag * vtag;
}

bool MeshSimplification::collapse_edge()
{   
    int e; // address of edge to collapse

    // a. takes out the lowest cost edge from queue
    if (Q.empty())
    {
        e = -1;
        return false;
    }

    std::tuple<double, int> currentEdge = Q.top();
    if (std::get<0>(currentEdge) == std::numeric_limits<double>::infinity())
    {
        e = -1;
        // min cost edge is infinite cost
        return false;
    }
    else e = std::get<1>(currentEdge);
    Q.pop();
    // b. add the new vertex into V
    int v1 = E(e, 0), v2 = E(e, 1);
    //double cost;
    Eigen::RowVectorXd p = C.row(e);
    //quadratic_error_simplification(e, cost, p);
    V.conservativeResize(V.rows() + 1, V.cols());
    int newVertexIndex = V.rows() - 1;
    V.row(newVertexIndex) = p;
    
    // c. remove old vertices from V, F so we could call edge_flaps to rebuild E, EF, EI, EMAP
    F.row(EF(e, 0)) = Eigen::Vector3i(0, 0, 0);
    F.row(EF(e, 1)) = Eigen::Vector3i(0, 0, 0);
    std::vector<int> facesOfV1 = verticesToFaces[v1];
    for (int face : facesOfV1)
    {
        if (face == EF(e, 0) || face == EF(e, 1))
        {
            continue;
        }
        verticesToFaces[newVertexIndex].push_back(face);
        for (int i = 0; i < 3; i++)
            if (F(face, i) == v1)
                F(face, i) = newVertexIndex;
    }
    std::vector<int> facesOfV2 = verticesToFaces[v2];
    for (int face : facesOfV2)
    {
        if (face == EF(e, 0) || face == EF(e, 1))
        {
            continue;
        }
        verticesToFaces[newVertexIndex].push_back(face);
        for (int i = 0; i < 3; i++)
            if (F(face, i) == v2)
                F(face, i) = newVertexIndex;
    }

    verticesToFaces.erase(v1);
    verticesToFaces.erase(v2);
    igl::edge_flaps(F, E, EMAP, EF, EI);

    return true;
}

std::vector<cg3d::MeshData> MeshSimplification::createDecimatedMesh(std::string filename)
{
    //std::cout << "EMAP:\n*********************\n" << EMAP << "\n*********************\n" << std::endl;
    //std::cout << "E:\n*********************\n" << E << "\n*********************\n" << std::endl;
    ////std::cout << "EF:\n*********************\n" << EF << "\n*********************\n" << std::endl;
    ////std::cout << "EI:\n*********************\n" << EI << "\n*********************\n" << std::endl;
    //std::cout << "EI:\n*********************\n" << EQ << "\n*********************\n" << std::endl;
    //std::cout << "V:" << V << std::endl;
    //std::cout << "F:\n*********************\n" << F << "\n*********************\n" << std::endl;

    ////std::cout << "normals before collapse\n" << currentMesh->data[0].vertexNormals << std::endl;
    //auto facezero = F.row(0);
    //std::cout << "face to calculate normal is:" << facezero << std::endl;
    //std::cout << "First vertex is " << V.row(facezero(0, 0)) << std::endl;
    //std::cout << "Second vertex is " << V.row(facezero(0, 1)) << std::endl;
    //std::cout << "Third vertex is " << V.row(facezero(0, 2)) << std::endl;
    //auto normalized = equation_plane(facezero, V);
    //std::cout << "Our normalization is " << normalized << std::endl;
    //std::cout << "other normalization is " << currentMesh->data[0].vertexNormals << std::endl;

    //std::vector<cg3d::MeshData> meshDataVector;
    //meshDataVector.push_back(currentMesh->data[0]);

    for (int i = 0; i < 10; i++)
    {
        const int max_iter = std::ceil(0.1 * Q.size());
        for (int j = 0; j < max_iter; j++)
        {
            if (!collapse_edge())
            {
                break;
            }
        }
        //igl::per_vertex_normals(V, F, N);
        //std::cout << "\nnormals after collapse number " << i << "\n" << N << std::endl;
        cg3d::Mesh temp("new mesh", V, F, currentMesh->data[0].vertexNormals, currentMesh->data[0].textureCoords);
        meshDataVector.push_back(temp.data[0]);

    }
    std::cout << "Decimated mesh size is: " << meshDataVector.size() << std::endl;
    return meshDataVector;
}

int main() {
    auto currentMesh{ cg3d::IglLoader::MeshFromFiles("cyl_igl", "data/cube.off")};
    Eigen::MatrixXd V, OV;
    Eigen::MatrixXi F, OF;

    igl::read_triangle_mesh("data/cube.off", OV, OF);
    // Prepare array-based edge data structures and priority queue
    Eigen::VectorXi EMAP;
    Eigen::MatrixXi E, EF, EI;
    igl::min_heap< std::tuple<double, int, int>> Q;
    Eigen::VectorXi EQ;
    // If an edge were collapsed, we'd collapse it to these points:
    Eigen::MatrixXd C;
    F = OV;
    V = OF;
    igl::edge_flaps(F, E, EMAP, EF, EI);

    std::cout << "EMAP:\n*********************\n" << EMAP << "\n*********************\n" << std::endl;
    std::cout << "E:\n*********************\n" << E << "\n*********************\n" << std::endl;
    std::cout << "EF:\n*********************\n" << EF << "\n*********************\n" << std::endl;
    std::cout << "EI:\n*********************\n" << EI << "\n*********************\n" << std::endl;
    std::cout << "EI:\n*********************\n" << EQ << "\n*********************\n" << std::endl;
    std::cout << "V:" << V << std::endl;
    std::cout << "F:\n*********************\n" << F << "\n*********************\n" << std::endl;
    return 0;
}
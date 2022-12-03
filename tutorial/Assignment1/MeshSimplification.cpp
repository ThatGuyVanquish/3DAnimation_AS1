#include "MeshSimplification.h"

MeshSimplification::MeshSimplification(std::string filename, int _decimations, const int _recalcQsRate) :
    currentMesh(cg3d::IglLoader::MeshFromFiles("Current Mesh", filename)),
    decimations(_decimations),
    collapseCounter(0),
    recalcQsRate(_recalcQsRate)
{
    createDecimatedMesh(filename);
}

/*

    HELPER METHODS

*/

std::shared_ptr<cg3d::Mesh> MeshSimplification::getMesh()
{
    return currentMesh;
}

Eigen::Vector4d MeshSimplification::ThreeDimVecToFourDim(Eigen::Vector3d vertex)
{
    Eigen::Vector4d newVertex;
    newVertex[0] = vertex[0], newVertex[1] = vertex[1], newVertex[2] = vertex[2], newVertex[3] = 1;
    return newVertex;
}

Eigen::Vector3d MeshSimplification::FourDimVecToThreeDim(Eigen::Vector4d vertex)
{
    Eigen::Vector3d newVertex;
    newVertex[0] = vertex[0], newVertex[1] = vertex[1], newVertex[2] = vertex[2];
    return newVertex;
}

/*

    METHODS TO CALCULATE Q MATRICES

*/

Eigen::Vector4d MeshSimplification::planeEquation(Eigen::Vector3i triangle, const Eigen::MatrixXd& V)
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

Eigen::Vector4d MeshSimplification::calculatePlaneNormal(const Eigen::MatrixXd& V, Eigen::Vector3d threeDimNormal, int vi)
{
    Eigen::Vector3d vertex = V.row(vi);
    double d = -(vertex[0] * threeDimNormal[0] + vertex[1] * threeDimNormal[1] + vertex[2] * threeDimNormal[2]);
    Eigen::Vector4d planeIn4D = ThreeDimVecToFourDim(threeDimNormal);
    planeIn4D[3] = d;
    return planeIn4D;
}

Eigen::Matrix4d MeshSimplification::calculateKp(Eigen::Vector4d planeVector)
{
    return planeVector * planeVector.transpose();
}

void MeshSimplification::calculateQs(const Eigen::MatrixXd& V, 
    const Eigen::MatrixXi &F, 
    Eigen::MatrixXd& faceNormals,
    std::vector<std::vector<int>>& verticesToFaces,
    std::vector<std::vector<int>>& facesBeforeIndex,
    std::vector<Eigen::Matrix4d>& verticesToQ)
{
    igl::per_face_normals(V, F, faceNormals); // Re-calculate the normals of all faces
    igl::vertex_triangle_adjacency(V.rows(), F, verticesToFaces, facesBeforeIndex);
    for (int i = 0; i < V.rows(); i++)
    {
        verticesToQ[i] = Eigen::Matrix4d::Zero();
        for (int j = 0; j < verticesToFaces[i].size(); j++)
        {
            int currentFace = verticesToFaces[i][j];
            Eigen::Vector3d currentFaceNormal = faceNormals.row(currentFace);
            Eigen::Vector4d normal = calculatePlaneNormal(V, currentFaceNormal, i);
            Eigen::Matrix4d Kp = calculateKp(normal);
            verticesToQ[i] += Kp;
        }
    }
}

/*
    
    METHODS TO CALCULATE COST AND PLACEMENT

*/

Eigen::Matrix4d MeshSimplification::calculateQDerive(Eigen::Matrix4d currentQ)
{
    auto lastRow = currentQ.row(3);
    lastRow[0] = 0, lastRow[1] = 0, lastRow[2] = 0, lastRow[3] = 1;
    return currentQ;
}

double MeshSimplification::calculateCost(Eigen::Vector4d vertex, Eigen::Matrix4d Q)
{
    double cost = vertex.transpose() * Q * vertex;
    return cost;
}

void MeshSimplification::preCalculateCostAndPos(
    const igl::decimate_cost_and_placement_callback& callback,
    const Eigen::MatrixXd &V,
    const Eigen::MatrixXi &F,
    const Eigen::MatrixXi &E,
    const Eigen::VectorXi &EMAP,
    const Eigen::MatrixXi &EF,
    const Eigen::MatrixXi &EI,
    igl::min_heap<std::tuple<double, int, int>>& Q,
    Eigen::VectorXi& EQ,
    Eigen::MatrixXd& C
)
{
    C.resize(E.rows(), V.cols());
    Q = {};
    EQ = Eigen::VectorXi::Zero(E.rows());
    {
        Eigen::VectorXd costs(E.rows());
        igl::parallel_for(E.rows(), [&](const int e)
            {
                double cost = e;
                Eigen::RowVectorXd p(1, 3);
                callback(e, V, F, E, EMAP, EF, EI, cost, p);
                C.row(e) = p;
                costs(e) = cost;
            }, 10000);
        for (int e = 0; e < E.rows(); e++)
        {
            Q.emplace(costs(e), e, 0);
        }
    }
}

void MeshSimplification::createDecimatedMesh(std::string fileName)
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::read_triangle_mesh(fileName, V, F);
    std::vector<Eigen::Matrix4d> verticesToQ(V.rows());
    Eigen::MatrixXi E;
    Eigen::VectorXi EMAP;
    Eigen::MatrixXi EF;
    Eigen::MatrixXi EI;
    igl::edge_flaps(F, E, EMAP, EF, EI);

    Eigen::MatrixXd faceNormals;
    std::vector<std::vector<int>> verticesToFaces;
    std::vector<std::vector<int>> facesBeforeIndex;

    calculateQs(V, F, faceNormals, verticesToFaces, facesBeforeIndex, verticesToQ);

    auto calculateCostAndPos = [&](
        const int e,
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi&/*F*/,
        const Eigen::MatrixXi& E,
        const Eigen::VectorXi&/*EMAP*/,
        const Eigen::MatrixXi&/*EF*/,
        const Eigen::MatrixXi&/*EI*/,
        double& cost,
        Eigen::RowVectorXd& p
        )
    {
        // Calculate the point p, the coordinates of the combined vertices
        int vertex1 = E(e, 0), vertex2 = E(e, 1);
        Eigen::Matrix4d Q1 = verticesToQ[vertex1], Q2 = verticesToQ[vertex2], Qtag = Q1 + Q2;

        Eigen::Matrix4d derivedQtag = calculateQDerive(Qtag); // Deriviated Q, 4th row is 0,0,0,1

        // Try to calculate Q' inverse, if possible the new point is ((Q')^-1) * Transpose(0,0,0,1)
        Eigen::Matrix4d derivedQtagInverse;
        bool invertible = false;
        double determinant;
        derivedQtag.computeInverseAndDetWithCheck(derivedQtagInverse, determinant, invertible);
        Eigen::Vector4d vtag;
        Eigen::Vector4d v0001;
        v0001[0] = 0, v0001[1] = 0, v0001[2] = 0, v0001[3] = 1;
        /*
            If Q' isn't invertible calculate the coordinates of the new vertex p by
            the minimum cost between the vertx v1, the vertex v2 and the middle of the edge e=(v1,v2), denoted by v3
        */
        if (!invertible)
        {
            Eigen::Vector4d v1 = ThreeDimVecToFourDim(V.row(vertex1)), v2 = ThreeDimVecToFourDim(V.row(vertex2));
            Eigen::Vector4d v3 = (v1 + v2) / 2;

            double costV1 = calculateCost(v1, Qtag), costV2 = calculateCost(v2, Qtag), costV3 = calculateCost(v3, Qtag);
            vtag = costV1 <= costV2 ?
                (costV1 <= costV3 ? v1 : v3) : // costV1 <= costV3 <= costV2, otherwise costV3 < costV1 <= costV2
                costV2 <= costV3 ? v2 : v3;    // costV2 <= costV1 <= costV3, otherwise costV3 < costV2 <= costV1
        }
        else
        {
            vtag = derivedQtagInverse * v0001;
        }
        p = FourDimVecToThreeDim(vtag);
        //p = (V.row(vertex1) + V.row(vertex2)) / 2;
        cost = calculateCost(vtag, Qtag);
    };

    igl::min_heap<std::tuple<double, int, int>> Q;
    Eigen::VectorXi EQ;
    Eigen::MatrixXd C;

    preCalculateCostAndPos(calculateCostAndPos, V, F, E, EMAP, EF, EI, Q, EQ, C);

    auto preCollapse = [&](
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi&/*F*/,
        const Eigen::MatrixXi&/*E*/,
        const Eigen::VectorXi&/*EMAP*/,
        const Eigen::MatrixXi&/*EF*/,
        const Eigen::MatrixXi&/*EI*/,
        const igl::min_heap<std::tuple<double, int, int>> /*Q*/,
        const Eigen::VectorXi /*EQ*/,
        const Eigen::MatrixXd&/*C*/,
        const int e
        )
    {
        return true;
    };

    auto postCollapse = [&](
        const Eigen::MatrixXd& V, 
        const Eigen::MatrixXi&/*F*/,
        const Eigen::MatrixXi &/*E*/,
        const Eigen::VectorXi &/*EMAP*/,
        const Eigen::MatrixXi &/*EF*/,
        const Eigen::MatrixXi &/*EI*/,
        const igl::min_heap<std::tuple<double,int,int>> /*Q*/, 
        const Eigen::VectorXi /*EQ*/, 
        const Eigen::MatrixXd &/*C*/, 
        const int e, 
        const int e1, 
        const int e2, 
        const int f1, 
        const int f2, 
        const bool collapsed
        )
    {
        if (collapsed)
        {
            int v1 = E(e, 0), v2 = E(e, 1);
            verticesToQ[v1] += verticesToQ[v2];
            verticesToQ[v2] = verticesToQ[v1];
        }
    };

    int currentNumOfEdges = E.rows();
    int counter = 0;
    int recalc = 0;
    for (int i = 0; i < decimations; i++)
    {
        int collapsed_edges = 0;
        const int max_iter = (std::ceil(0.1 * currentNumOfEdges));
        for (int j = 0; j < max_iter; j++)
        {
            if (counter == recalcQsRate)
            {
                calculateQs(V, F, faceNormals, verticesToFaces, facesBeforeIndex, verticesToQ);
                recalc++;
                counter = 0;
            }
            std::tuple<double, int, int> currentTop = Q.top();
            if (!igl::collapse_edge(calculateCostAndPos, preCollapse, postCollapse, V, F, E, EMAP, EF, EI, Q, EQ, C))
            {
                break;
            }
            counter++;
            collapsed_edges += 3;
            int e = std::get<1>(currentTop);
            int v1 = E(e, 0), v2 = E(e, 1);

            std::cout << "Edge: " << e << ", Cost = " << std::get<0>(currentTop) << ", New position: ("
                << C.row(e) << ")" << std::endl;
        }
        if (collapsed_edges > 0)
        {
            collapseCounter += collapsed_edges;
            currentNumOfEdges -= collapsed_edges;
            currentMesh->data.push_back(
                { V, F, currentMesh->data[0].vertexNormals, currentMesh->data[0].textureCoords }
            );
        }
    }
    std::cout << "Initial Number of edges: " << E.rows() << std::endl;
    std::cout << "Collapsed edges: " << collapseCounter << std::endl;
    std::cout << "Qs were recalculated " << recalc << " times" << std::endl;
}
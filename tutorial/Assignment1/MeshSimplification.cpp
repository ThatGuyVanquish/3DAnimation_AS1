#include "MeshSimplification.h"

MeshSimplification::MeshSimplification(std::string filename, int _decimations) :
    currentMesh(cg3d::IglLoader::MeshFromFiles("Current Mesh", filename)),
    decimations(_decimations),
    collapseCounter(0),
    verticesToUpdate(),
    QResetInterval(0)
{
    V = currentMesh->data[0].vertices, F = currentMesh->data[0].faces;
    igl::edge_flaps(F, E, EMAP, EF, EI);
    faceNormals = Eigen::MatrixXd::Zero(F.rows(), 3);
    //igl::per_face_normals(V, F, faceNormals);
    for (int i = 0; i < V.rows(); i++)
        verticesToQ.push_back(Eigen::Matrix4d::Zero());
    //std::cout << "Per face normals\n" << blyat;
    //std::cout << "\n";
    Init();
    createDecimatedMesh();
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
    newVertex[0] = vertex[0], newVertex[1] = vertex[1], newVertex[2] = vertex[2], newVertex[3] = d;
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

Eigen::Vector4d MeshSimplification::equation_plane(Eigen::Vector3i triangle)
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

Eigen::Vector4d MeshSimplification::calculatePlaneNormal(int face)
{
    Eigen::Vector3d planeIn3D = faceNormals.row(face);
    Eigen::Vector3d vertex = V.row(F(face,0));
    int d = -(vertex[0] * planeIn3D[0] + vertex[1] * planeIn3D[1] + vertex[2] * planeIn3D[2]);
    Eigen::Vector4d planeIn4D = ThreeDimVecToFourDim(planeIn3D);
    std::cout << "Current face normal is\n" << planeIn4D << std::endl;
    //std::cout << "Current normal based on algebra is\n" << equation_plane(F.row(face)) << std::endl;
    return planeIn4D;
}

Eigen::Matrix4d MeshSimplification::calculateKp(Eigen::Vector4d planeVector)
{
    return planeVector * planeVector.transpose();
}

double MeshSimplification::calculateCost(const int vertex)
{
    Eigen::Vector4d v = ThreeDimVecToFourDim(V.row(vertex));
    Eigen::Matrix4d QofV = verticesToQ.at(vertex);
    return v.transpose() * QofV * v;
}

void MeshSimplification::buildVerticesToFaces()
{
    for (int i = 0; i < F.rows(); i++)
    {
        Eigen::VectorXi currentRowInF = F.row(i);
        // prepare vertex to face map
        for (int j = 0; j < currentRowInF.cols(); j++)
        {
            int vertexId = currentRowInF(j);
            verticesToFaces[vertexId].insert(i);
        }
    }
}

void MeshSimplification::calculateQ(int v)
{
    std::set<int> faces = verticesToFaces[v];
    verticesToQ[v] = Eigen::Matrix4d::Zero();
    for (int face : faces)
    {
        verticesToQ[v] += calculateKp(calculatePlaneNormal(face));
    }
}

void MeshSimplification::calculateQs(std::vector<int> verticesToCalculate)
{
    igl::per_face_normals(V, F, faceNormals); // Re-calculate the normals of all faces
    for (int vertex : verticesToCalculate)
    {
        calculateQ(vertex);
    }
}

/*

    METHODS FOR COLLAPSING EDGES

*/

void MeshSimplification::updateVerticesToFacesPostCollapse(const int e)
{
    int v0 = E(e, 0), v1 = E(e, 1);
    int f0 = EF(e, 0), f1 = EF(e, 1);
    int n0 = EI(e, 0), n1 = EI(e, 1);
    verticesToFaces[v0].erase(f0);
    verticesToFaces[v0].erase(f1);
    verticesToFaces[v1].erase(f0);
    verticesToFaces[v1].erase(f1);
    verticesToFaces[n0].erase(f0);
    verticesToFaces[n1].erase(f1);

    // Copy all of the remaining faces from v0 to v1 and vice versa
    verticesToFaces[v0].merge(verticesToFaces[v1]);
    verticesToFaces[v1] = verticesToFaces[v0];
}

void MeshSimplification::post_collapse(const int e)
{
    updateVerticesToFacesPostCollapse(e);
    // update the Q matrix of v'
    verticesToQ[E(e, 0)] += verticesToQ[E(e, 1)];
    verticesToQ[E(e, 1)] = verticesToQ[(e, 0)];

    // if statement based on count to see if we should re calculate Qs
    if (collapseCounter >= QResetInterval)
    {
        calculateQs(std::vector<int>(verticesToUpdate.begin(), verticesToUpdate.end()));
        verticesToUpdate.clear();
        collapseCounter = 0;
    }
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
        double costV1 = v1.transpose() * Qtag * v1, costV2 = v2.transpose() * Qtag * v2, costV3 = v3.transpose() * Qtag * v3;
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
    cost = vtag.transpose() * Qtag * vtag;
}

bool MeshSimplification::collapse_edge()
{
    int e, e1, e2, f1, f2;

    std::tuple<double, int, int> costEdgeTimestamp;

    while (true)
    {
        // Check if Q is empty
        if (Q.empty())
        {
            // no edges to collapse
            e = -1;
            return false;
        }
        // Pop from Q
        costEdgeTimestamp = Q.top();
        if (std::get<0>(costEdgeTimestamp) == std::numeric_limits<double>::infinity())
        {
            e = -1;
            // Min cost edge is infinite cost
            return false;
        }
        Q.pop();
        e = std::get<1>(costEdgeTimestamp);
        // Check if matches timestamp
        if (std::get<2>(costEdgeTimestamp) == EQ(e))
        {
            break;
        }
        // must be stale or dead
        assert(std::get<2>(costEdgeTimestamp) < EQ(e) || EQ(e) == -1);
        // therefore try again
    }

    std::vector<int> /*Nse,*/Nsf, Nsv;
    igl::circulation(e, true, F, EMAP, EF, EI,/*Nse,*/Nsv, Nsf);
    std::vector<int> /*Nde,*/Ndf, Ndv;
    igl::circulation(e, false, F, EMAP, EF, EI,/*Nde,*/Ndv, Ndf);

    bool collapsed = igl::collapse_edge(
        e, C.row(e),
        Nsv, Nsf, Ndv, Ndf,
        V, F, E, EMAP, EF, EI, e1, e2, f1, f2);

    if (collapsed)
    {
        for (int i = 0; i < Nsv.size(); i++)
            verticesToUpdate.insert(Nsv.at(i));
        for (int i = 0; i < Ndv.size(); i++)
            verticesToUpdate.insert(Ndv.at(i));
        post_collapse(e);

        std::cout << "Edge: " << e << ", Cost = " << std::get<0>(costEdgeTimestamp) << ", New position: ("
            << C.row(e) << ")" << std::endl;

        // Erase the two, other collapsed edges by marking their timestamps as -1
        EQ(e1) = -1;
        EQ(e2) = -1;
        
        // update local neighbors
        // loop over original face neighbors

        std::vector<int> Nf;
        Nf.reserve(Nsf.size() + Ndf.size()); // preallocate memory
        Nf.insert(Nf.end(), Nsf.begin(), Nsf.end());
        Nf.insert(Nf.end(), Ndf.begin(), Ndf.end());
        
        std::sort(Nf.begin(), Nf.end());
        Nf.erase(std::unique(Nf.begin(), Nf.end()), Nf.end());

        // Collect all edges that must be updated
        std::vector<int> Ne;
        Ne.reserve(3 * Nf.size());
        for (auto& n : Nf)
        {
            if (F(n, 0) != IGL_COLLAPSE_EDGE_NULL ||
                F(n, 1) != IGL_COLLAPSE_EDGE_NULL ||
                F(n, 2) != IGL_COLLAPSE_EDGE_NULL)
            {
                for (int v = 0; v < 3; v++)
                {
                    // get edge id
                    const int ei = EMAP(v * F.rows() + n);
                    Ne.push_back(ei);
                }
            }
        }
        // Only process edge once
        std::sort(Ne.begin(), Ne.end());
        Ne.erase(std::unique(Ne.begin(), Ne.end()), Ne.end());
        for (auto& ei : Ne)
        {
            // Compute cost and potential placement
            double cost;
            Eigen::RowVectorXd newCoords = Eigen::RowVectorXd::Zero(3);
            quadratic_error_simplification(ei, cost, newCoords);
            // Increment timestamp
            EQ(ei)++;
            // Replace in queue
            Q.emplace(cost, ei, EQ(ei));
            C.row(ei) = newCoords;
        }
    }
    else
    {
        // Re-insert with infinite weight (the provided cost function must **not**
        // have given this un-collapsable edge inf cost already)
        // Increment timestamp
        EQ(e)++;
        // Replace in queue
        Q.emplace(std::numeric_limits<double>::infinity(), e, EQ(e));
    }
    return collapsed;
}

void MeshSimplification::Init()
{
    buildVerticesToFaces();
    C = Eigen::MatrixXd::Zero(E.rows(), 3);
    Eigen::VectorXd costs(E.rows());
    Q = {};
    EQ = Eigen::VectorXi::Zero(E.rows());

    /*
    *
    * BUILD COSTS HEAP:
    *
    */
    // Calculate Q matrices for every vertex
    std::vector<int> verticesToCalculate;
    for (int i = 0; i < V.rows(); i++)
    {
        verticesToCalculate.push_back(i);
    }
    
    calculateQs(verticesToCalculate);

    // Calculate cost of collapse and new vertex placement after collapse
    for (int e = 0; e < E.rows(); e++)
    {
        double cost = e;
        Eigen::RowVectorXd p(1, 3);
        quadratic_error_simplification(e, cost, p);
        C.row(e) = p;
        costs(e) = cost;
    }
    //igl::parallel_for(E.rows(), [&](const int e)
    //    {
    //        double cost = e;
    //        Eigen::RowVectorXd p(1, 3);
    //        quadratic_error_simplification(e, cost, p);
    //        C.row(e) = p;
    //        costs(e) = cost;
    //    }, 10000);

    for (int e = 0; e < E.rows(); e++)
    {
        Q.emplace(costs(e), e, 0);
    }
}

void MeshSimplification::createDecimatedMesh()
{
    for (int i = 0; i < decimations; i++)
    {
        int collapsed_edges = 0;
        const int max_iter = (std::ceil(0.1 * Q.size()));
        //std::cout << "Max iter is " << max_iter << " at decimation i = " << i << std::endl;
        QResetInterval = (max_iter / 10) > 0 ? max_iter / 10 : 1;
        for (int j = 0; j < max_iter; j++)
        {
            if (!collapse_edge())
            {
                break;
            }
            collapsed_edges++;
            collapseCounter++;
        }
        if (collapsed_edges > 0)
        {
            currentMesh->data.push_back(
                { V, F, currentMesh->data[0].vertexNormals, currentMesh->data[0].textureCoords }
            );
        }
    }
}
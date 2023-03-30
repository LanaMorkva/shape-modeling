#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui.h>
/*** insert any necessary libigl headers here ***/
#include <igl/per_face_normals.h>
#include <igl/copyleft/marching_cubes.h>
#include <igl/bounding_box_diagonal.h>

// added class to handle spatial indices and brute force method
#include "ConstrHelper.h"
#include "MLSHelper.h"
#include <queue>

using namespace std;
using Viewer = igl::opengl::glfw::Viewer;

// Input: imported points, #P x3
Eigen::MatrixXd P;

// Input: imported normals, #P x3
Eigen::MatrixXd N;

// Normals evaluated via PCA method, #P x3
Eigen::MatrixXd NP;

// Intermediate result: constrained points, #C x3
Eigen::MatrixXd constrained_points;

// Intermediate result: implicit function values at constrained points, #C x1
Eigen::VectorXd constrained_values;

// Parameter: degree of the polynomial
int polyDegree = 0;

// Parameter: Wendland weight function radius (make this relative to the size of the mesh)
double wendlandRadius = 0.1;

// Parameter: grid resolution
int resolution = 20;

// Diagonal size
double diag_size = 0;

// Intermediate result: grid points, at which the imlicit function will be evaluated, #G x3
Eigen::MatrixXd grid_points;

// Intermediate result: implicit function values at the grid points, #G x1
Eigen::VectorXd grid_values;

// Intermediate result: grid point colors, for display, #G x3
Eigen::MatrixXd grid_colors;

// Intermediate result: grid lines, for display, #L x6 (each row contains
// starting and ending point of line segment)
Eigen::MatrixXd grid_lines;

// Output: vertex array, #V x3
Eigen::MatrixXd V;

// Output: face array, #F x3
Eigen::MatrixXi F;

// Output: face normals of the reconstructed mesh, #F x3
Eigen::MatrixXd FN;

std::string filename_def;
std::string filename_par;
bool non_aligned = false;
int k_neighb = 5;
bool flipN = false;
bool autoNormFlipOn = false;
bool spatialIndexOn = true;

// Functions
void createGrid();
void nonAlignedGrid();
void alignedGrid();
void evaluateImplicitFunc();
void evaluateImplicitFunc_PolygonSoup();
void getLines();
void pcaNormal();
void autoNormFlip(const std::vector<std::vector<int>> &conn_graph);
bool callback_key_down(Viewer &viewer, unsigned char key, int modifiers);

// Creates a grid_points array for the simple sphere example. The points are
// stacked into a single matrix, ordered first in the x, then in the y and
// then in the z direction. If you find it necessary, replace this with your own
// function for creating the grid.
void createGrid()
{
    grid_points.resize(0, 3);
    grid_colors.resize(0, 3);
    grid_lines.resize(0, 6);
    grid_values.resize(0);
    V.resize(0, 3);
    F.resize(0, 3);
    FN.resize(0, 3);

    if (non_aligned) {
        nonAlignedGrid();
    } else {
        alignedGrid();
    }
}

void alignedGrid() {
    // Grid bounds: axis-aligned bounding box
    Eigen::Array3d bb_min = P.colwise().minCoeff();
    Eigen::Array3d bb_max = P.colwise().maxCoeff();

    // enlarge bound box a bit
    bb_min -= 0.01*diag_size;
    bb_max += 0.01*diag_size;
    // Bounding box dimensions
    Eigen::RowVector3d dim = bb_max - bb_min;

    // Grid spacing
    const double dx = dim[0] / (double)(resolution - 1);
    const double dy = dim[1] / (double)(resolution - 1);
    const double dz = dim[2] / (double)(resolution - 1);
    // 3D positions of the grid points -- see slides or marching_cubes.h for ordering
    grid_points.resize(resolution * resolution * resolution, 3);
    // Create each gridpoint
    for (unsigned int x = 0; x < resolution; ++x)
    {
        for (unsigned int y = 0; y < resolution; ++y)
        {
            for (unsigned int z = 0; z < resolution; ++z)
            {
                // Linear index of the point at (x,y,z)
                int index = x + resolution * (y + resolution * z);
                // 3D point at (x,y,z)
                grid_points.row(index) = bb_min.matrix() + Eigen::Vector3d(x * dx, y * dy, z * dz);
            }
        }
    }
}

void nonAlignedGrid() {
    Eigen::MatrixXd centered = P.rowwise() - P.colwise().mean();
    Eigen::MatrixXd cov = (centered.adjoint() * centered) / double(P.rows() - 1);
    Eigen::EigenSolver<Eigen::MatrixXd> es(cov);
    Eigen::MatrixXd eigenVec = es.pseudoEigenvectors();
    Eigen::MatrixXd egVecInv = eigenVec.inverse();

    auto P_rot = P;
    for (int i = 0; i < P.rows(); i++) {
        P_rot.row(i) = egVecInv*P.row(i).transpose();
    }
    Eigen::Array3d bb_min = P_rot.colwise().minCoeff();
    Eigen::Array3d bb_max = P_rot.colwise().maxCoeff();

    bb_min -= 0.01*diag_size;
    bb_max += 0.01*diag_size;

    Eigen::Vector3d dim = bb_max - bb_min;

    // Grid spacing
    const double dx = dim[0] / (double)(resolution - 1);
    const double dy = dim[1] / (double)(resolution - 1);
    const double dz = dim[2] / (double)(resolution - 1);
    // 3D positions of the grid points -- see slides or marching_cubes.h for ordering
    grid_points.resize(resolution * resolution * resolution, 3);
    // Create each gridpoint
    for (unsigned int x = 0; x < resolution; ++x)
    {
        for (unsigned int y = 0; y < resolution; ++y)
        {
            for (unsigned int z = 0; z < resolution; ++z)
            {
                // Linear index of the point at (x,y,z)
                int index = x + resolution * (y + resolution * z);
                Eigen::Vector3d baseVec = bb_min + Eigen::Array3d(x * dx, y * dy, z * dz);
                // 3D point at (x,y,z)
                grid_points.row(index) = eigenVec * baseVec;
            }
        }
    }
}

void evaluateImplicitFunc()
{
    // scale wendlandRadius to the object size
    auto radius = wendlandRadius * diag_size;

    MLSHelper mlsHelper(grid_points, constrained_points, constrained_values, resolution,
                        polyDegree, radius);
    auto start = chrono::high_resolution_clock::now();

    // pass parameter "spatialIndexOn" to use fast or slow implementation
    mlsHelper.calcGridValues(spatialIndexOn);
    auto end = chrono::high_resolution_clock::now();
    auto time = std::chrono::duration<double, std::milli>(end - start).count();
    grid_values = mlsHelper.getGridValues();

    printf("Execution time MLS %f (%s)\n", time, spatialIndexOn ? "fast" : "slow");
}

void evaluateImplicitFunc_PolygonSoup()
{
    // Replace with your code here, for "key == '5'"
    auto radius = wendlandRadius * diag_size;
    MLSHelper mlsHelper(grid_points, P, {}, resolution, polyDegree, radius);
    mlsHelper.calcGridValues_NormConstr(N);
    grid_values = mlsHelper.getGridValues();
}

// Code to display the grid lines given a grid structure of the given form.
// Assumes grid_points have been correctly assigned
// Replace with your own code for displaying lines if need be.
void getLines()
{
    int nnodes = grid_points.rows();
    grid_lines.resize(3 * nnodes, 6);
    int numLines = 0;

    for (unsigned int x = 0; x < resolution; ++x)
    {
        for (unsigned int y = 0; y < resolution; ++y)
        {
            for (unsigned int z = 0; z < resolution; ++z)
            {
                int index = x + resolution * (y + resolution * z);
                if (x < resolution - 1)
                {
                    int index1 = (x + 1) + y * resolution + z * resolution * resolution;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
                if (y < resolution - 1)
                {
                    int index1 = x + (y + 1) * resolution + z * resolution * resolution;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
                if (z < resolution - 1)
                {
                    int index1 = x + y * resolution + (z + 1) * resolution * resolution;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
            }
        }
    }

    grid_lines.conservativeResize(numLines, Eigen::NoChange);
}

// Estimation of the normals via PCA.
void pcaNormal()
{
    auto k_neighb_alg = [=](const Eigen::Vector3d &val) {
        std::vector<int> indices;
        std::vector<double> distV;
        int maxValInd = 0;
        double maxDist = DBL_MAX;
        for (int i = 0; i < k_neighb; i++) {
            indices.push_back(-1);
            distV.push_back(DBL_MAX);
        }
        for (int i = 0; i < P.rows(); i++) {
            Eigen::Vector3d p = P.row(i);
            auto dist = (p-val).norm();
            if (dist < maxDist) {
                distV[maxValInd] = dist;
                indices[maxValInd] = i;
                maxValInd = std::max_element(distV.begin(), distV.end()) - distV.begin();
                maxDist = distV[maxValInd];
            }
        }
        return indices;
    };
    std::vector<std::vector<int>> conn_graph;
    NP = Eigen::MatrixXd(P.rows(), 3);
    conn_graph.resize(P.rows());

    for (int i = 0; i < P.rows(); i++) {
        Eigen::Vector3d p = P.row(i);
        const auto &indices = k_neighb_alg(p);
        Eigen::MatrixXd PI = Eigen::MatrixXd(indices.size(), 3);
        for (int j = 0; j < indices.size(); j++) {
            PI.row(j) = P.row(indices[j]);
            conn_graph[i].push_back(indices[j]);
        }
        Eigen::MatrixXd centered = PI.rowwise() - PI.colwise().mean();
        Eigen::MatrixXd cov = (centered.adjoint() * centered) / double(PI.rows() - 1);
        Eigen::EigenSolver<Eigen::MatrixXd> es(cov);
        Eigen::MatrixXd eigenVec = es.pseudoEigenvectors();
        Eigen::MatrixXd eigenVal = es.pseudoEigenvalueMatrix();
        int minEValInd = -1; double minEVal = DBL_MAX;
        for (int j = 0; j < eigenVal.rows(); j++) {
            auto val = eigenVal.coeff(j,j);
            if (val < minEVal) {
                minEVal = val;
                minEValInd = j;
            }
        }
        auto n = eigenVec.col(minEValInd);
        if (!autoNormFlipOn && n.dot(N.row(i)) < 0) {
            NP.row(i) = -n;
        } else {
            NP.row(i) = n;
        }
    }
    if (autoNormFlipOn) {
        autoNormFlip(conn_graph);
    }
}

void autoNormFlip(const std::vector<std::vector<int>> &conn_graph){
    int visitedNum = 0;
    std::vector<bool> visited;
    std::queue<int> visit_order;
    visit_order.push(0);
    visited.resize(conn_graph.size(), false);
    if (flipN) {
        NP.row(0) = -NP.row(0);
    }
    while (visitedNum < conn_graph.size()) {
        int curr;
        if (visit_order.empty()) {
            for (int i = 0; i < conn_graph.size(); i++) {
                if (!visited[i]) {
                    curr = i;
                    break;
                }
            }
        } else {
            curr = visit_order.front();
            visit_order.pop();
        }
        const auto &neighb = conn_graph[curr];
        for (int index : neighb) {
            if (index!=curr && !visited[index]) {
                const auto &neighb_n = conn_graph[index];
                auto it = std::find(neighb_n.begin(), neighb_n.end(), curr);
                if (it != neighb_n.end()) {
                    auto n_curr = NP.row(curr);
                    auto n_next = NP.row(index);
                    if (n_curr.dot(n_next) < 0) {
                        NP.row(index) = -n_next;
                    }
                    visit_order.push(index);
                    visited[index] = true;
                }
            }
        }
        visited[curr] = true;
        visitedNum++;
    }
}

bool callback_key_down(Viewer &viewer, unsigned char key, int modifiers)
{
    if (key == '1')
    {
        // Show imported points
        viewer.data().clear();
        viewer.core().align_camera_center(P);
        viewer.data().point_size = 11;
        viewer.data().add_points(P, Eigen::RowVector3d(0, 0, 0));
    }

    if (key == '2')
    {
        // Show all constraints
        viewer.data().clear();
        viewer.core().align_camera_center(P);
        // Add your code for computing auxiliary constraint points here
        // Add code for displaying all points, as above

        auto start = chrono::high_resolution_clock::now();

        // pass parameter "spatialIndexOn" to use fast or slow implementation
        ConstrHelper helper(P, N, spatialIndexOn, resolution);
        helper.calcConstraints();
        auto end = chrono::high_resolution_clock::now();
        auto time = std::chrono::duration<double, std::milli>(end - start).count();

        printf("Execution time constraints %f (%s)\n", time, spatialIndexOn ? "fast" : "slow");

        constrained_points = helper.constrPoints();
        constrained_values = helper.constrValues();

        Eigen::MatrixXd constrained_colors = Eigen::MatrixXd::Zero(P.rows() * 3, 3);
        constrained_colors.block(0, 2, P.rows(), 1) =
                Eigen::VectorXd::Ones(P.rows());
        constrained_colors.block(P.rows(), 0, P.rows(), 1) =
                Eigen::VectorXd::Ones(P.rows());
        constrained_colors.block(2*P.rows(), 1, P.rows(), 1) =
                Eigen::VectorXd::Ones(P.rows());
        viewer.data().point_size = 8;
        viewer.data().add_points(constrained_points, constrained_colors);
    }

    if (key == '3')
    {
        // add constraints again in case if 2 wasn't pressed or need to rebuild
        callback_key_down(viewer, '2', modifiers);

        // Show grid points with colored nodes and connected with lines
        viewer.data().clear();
        viewer.core().align_camera_center(P);

        // Make grid
        createGrid();

        // Evaluate implicit function
        evaluateImplicitFunc();

        // get grid lines
        getLines();

        // Code for coloring and displaying the grid points and lines
        // Assumes that grid_values and grid_points have been correctly assigned.
        grid_colors.setZero(grid_points.rows(), 3);

        // Build color map
        for (int i = 0; i < grid_points.rows(); ++i)
        {
            double value = grid_values(i);
            if (value < 0)
            {
                grid_colors(i, 1) = 1;
            }
            else
            {
                if (value > 0)
                    grid_colors(i, 0) = 1;
            }
        }

        // Draw lines and points
        viewer.data().point_size = 8;
        viewer.data().add_points(grid_points, grid_colors);
        viewer.data().add_edges(grid_lines.block(0, 0, grid_lines.rows(), 3),
                                grid_lines.block(0, 3, grid_lines.rows(), 3),
                                Eigen::RowVector3d(0.8, 0.8, 0.8));
    }

    if (key == '4')
    {
        // In case keys were not pressed before
        callback_key_down(viewer, '3', modifiers);

        // Show reconstructed mesh
        viewer.data().clear();

        // Code for computing the mesh (V,F) from grid_points and grid_values
        if ((grid_points.rows() == 0) || (grid_values.rows() == 0))
        {
            cerr << "Not enough data for Marching Cubes !" << endl;
            return true;
        }
        // Run marching cubes
        igl::copyleft::marching_cubes(grid_values, grid_points, resolution, resolution, resolution, V, F);
        if (V.rows() == 0)
        {
            cerr << "Marching Cubes failed!" << endl;
            return true;
        }

        igl::per_face_normals(V, F, FN);
        viewer.data().set_mesh(V, F);
        viewer.data().show_lines = true;
        viewer.data().show_faces = true;
        viewer.data().set_normals(FN);
        igl::writeOFF(filename_def, V, F);
    }

    if (key == '5')
    {
        // Use the structure for key=='3' but replace the function evaluateImplicitFunc();
        // with a function performing the approximation of the implicit surface from polygon soup
        // Ref: Chen Shen, James F. O’Brien, and Jonathan Richard Shewchuk. Interpolating and approximating implicit surfaces from polygon soup.

        // Show grid points with colored nodes and connected with lines
        viewer.data().clear();
        viewer.core().align_camera_center(P);

        // Make grid
        createGrid();

        // Evaluate implicit function --> Function to be modified here
        evaluateImplicitFunc_PolygonSoup();

        // get grid lines
        getLines();

        // Display the reconstruction
        callback_key_down(viewer, '4', modifiers);
    }

    if (key == '6' || key == '7' || key == '8')
    {
        // Implement PCA Normal Estimation --> Function to be modified here
        pcaNormal();

        // To use the normals estimated via PCA instead of the input normals and then restaurate the input normals
        Eigen::MatrixXd N_tmp = N;
        N = NP;

        switch (key)
        {
        case '6':
            callback_key_down(viewer, '2', modifiers);
            break;
        case '7':
            callback_key_down(viewer, '3', modifiers);
            break;
        case '8':
            callback_key_down(viewer, '3', modifiers);
            callback_key_down(viewer, '4', modifiers);
            break;
        default:
            break;
        }

        // Restore input normals
        N = N_tmp;
    }

    return true;
}

bool callback_load_mesh(Viewer &viewer, string filename)
{
    igl::readOFF(filename, P, F, N);
    callback_key_down(viewer, '1', 0);
    return true;
}

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        cout << "Usage ex2_bin <mesh.off>" << endl;
        igl::readOFF("../data/sphere.off", P, F, N);
        filename_def = "../res/sphere.off";
    }
    else
    {
        // Read points and normals
        igl::readOFF(argv[1], P, F, N);
        filename_def = argv[1];
        auto it = filename_def.find("/data/");
        if (it != std::string::npos) {
            filename_def.replace(it, 6, "/res/");
        }
    }
    diag_size = igl::bounding_box_diagonal(P);

    Viewer viewer;
    igl::opengl::glfw::imgui::ImGuiPlugin plugin;
    viewer.plugins.push_back(&plugin);
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    plugin.widgets.push_back(&menu);

    viewer.callback_key_down = callback_key_down;

    menu.callback_draw_viewer_menu = [&]()
    {
        // Draw parent menu content
        menu.draw_viewer_menu();

        // Add new group
        if (ImGui::CollapsingHeader("Reconstruction Options", ImGuiTreeNodeFlags_DefaultOpen))
        {
            // Expose variable directly ...
            ImGui::Checkbox("Spatial index", &spatialIndexOn);
            ImGui::InputInt("Resolution", &resolution, 0, 0);
            ImGui::InputInt("polyDegree", &polyDegree, 0, 0);
            ImGui::InputDouble("wendlandRadius", &wendlandRadius, 0, 0);
            ImGui::Checkbox("Non-aligned grid", &non_aligned);

            if (ImGui::Button("Reset Grid", ImVec2(-1, 0)))
            {
                if (polyDegree < 0 || polyDegree > 2) {
                    polyDegree = 0;
                }
                if (wendlandRadius < 0) {
                    wendlandRadius = 0.1;
                }
                std::cout << "ResetGrid\n";
                // Recreate the grid
                createGrid();
                // Switch view to show the grid
                callback_key_down(viewer, '3', 0);
            }
        }
        if (ImGui::CollapsingHeader("Normals Estimation", ImGuiTreeNodeFlags_DefaultOpen)) {

            ImGui::InputInt("K-Nearest Neighbours", &k_neighb, 0, 0);
            ImGui::Checkbox("Auto normal flipping", &autoNormFlipOn);
            ImGui::Checkbox("Flip start normal", &flipN);
        }
    };

    callback_key_down(viewer, '1', 0);

    viewer.launch();
}

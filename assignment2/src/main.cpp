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
int resolution = 10;

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

// Functions
void createGrid();
void evaluateImplicitFunc();
void evaluateImplicitFunc_PolygonSoup();
void getLines();
void pcaNormal();
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

void evaluateImplicitFunc()
{
    // scale wendlandRadius to the object size
    auto radius = wendlandRadius * diag_size;

    // slow
    MLSHelper mlsHelper_sl(grid_points, constrained_points, constrained_values, resolution,
                        polyDegree, radius);
    auto start = chrono::high_resolution_clock::now();
    mlsHelper_sl.calcGridValues();
    auto end = chrono::high_resolution_clock::now();
    auto slow_grid_values = mlsHelper_sl.getGridValues();
    auto timeSlow = std::chrono::duration<double, std::milli>(end - start).count();

    // fast
    MLSHelper mlsHelper(grid_points, constrained_points, constrained_values, resolution,
                        polyDegree, radius);
    start = chrono::high_resolution_clock::now();
    mlsHelper.calcGridValues(true);
    end = chrono::high_resolution_clock::now();
    grid_values = mlsHelper.getGridValues();
    auto timeFast = std::chrono::duration<double, std::milli>(end - start).count();

    printf("Execution time MLS: slow %f; fast %f\n", timeSlow, timeFast);
    if (!grid_values.isApprox(slow_grid_values)) {
        printf("Error (MLS): values are not equal\n");
    } else {
        printf("Yay (MLS): they are equal\n");
    }
}

void evaluateImplicitFunc_PolygonSoup()
{
    // Replace with your code here, for "key == '5'"
    evaluateImplicitFunc();
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
    NP = -N; // to be replaced with your code
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
        ConstrHelper helper_fast(P, N, true, resolution);
        helper_fast.calcConstraints();
        auto end = chrono::high_resolution_clock::now();
        auto timeFast = std::chrono::duration<double, std::milli>(end - start).count();

        start = chrono::high_resolution_clock::now();
        ConstrHelper helper_slow(P, N);
        helper_slow.calcConstraints();
        end = chrono::high_resolution_clock::now();
        auto timeSlow = std::chrono::duration<double, std::milli>(end - start).count();

        printf("Execution time constraints: slow %f; fast %f\n", timeSlow, timeFast);

        constrained_points = helper_fast.constrPoints();
        constrained_values = helper_fast.constrValues();
        auto slow_constrp = helper_slow.constrPoints();

        if (!slow_constrp.isApprox(constrained_points)) {
            printf("Error (constraints): values are not equal\n");
        } else {
            printf("Yay (constraints): they are equal\n");
        }

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
        // Show grid points with colored nodes and connected with lines
        viewer.data().clear();
        viewer.core().align_camera_center(P);

        // add constraints again in case if 2 wasn't pressed or need to rebuild
        ConstrHelper helper_fast(P, N, true, resolution);
        helper_fast.calcConstraints();

        constrained_points = helper_fast.constrPoints();
        constrained_values = helper_fast.constrValues();

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
        /*** end: sphere example ***/
    }

    if (key == '4')
    {
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
    }
    else
    {
        // Read points and normals
        igl::readOFF(argv[1], P, F, N);
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
            ImGui::InputInt("Resolution", &resolution, 0, 0);
            ImGui::InputInt("polyDegree", &polyDegree, 0, 0);
            ImGui::InputDouble("wendlandRadius", &wendlandRadius, 0, 0);

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
    };

    callback_key_down(viewer, '1', 0);

    viewer.launch();
}

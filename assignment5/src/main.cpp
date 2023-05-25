#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui.h>
#include <igl/slice_into.h>
#include <igl/rotate_by_quat.h>
#include <igl/sum.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_face_normals.h>
#include <igl/adjacency_list.h>
#include <igl/heat_geodesics.h>

#include "Lasso.h"
#include "Colors.h"

//activate this for alternate UI (easier to debug but no interactive updates, turn this OFF for your report)
//#define UPDATE_ONLY_ON_UP

using namespace std;
using namespace Eigen;
using Viewer = igl::opengl::glfw::Viewer;

Viewer viewer;

//vertex array, #V x3
Eigen::MatrixXd V(0, 3), V_cp(0, 3);
//face array, #F x3
Eigen::MatrixXi F(0, 3);

//mouse interaction
enum MouseMode {
    SELECT, TRANSLATE, ROTATE, NONE
};
MouseMode mouse_mode = NONE;
bool doit = false;
int down_mouse_x = -1, down_mouse_y = -1;

//for selecting vertices
std::unique_ptr<Lasso> lasso;
//list of currently selected vertices
Eigen::VectorXi selected_v(0, 1);

//for saving constrained vertices
//vertex-to-handle index, #V x1 (-1 if vertex is free)
Eigen::VectorXi handle_id(0, 1);
//list of all vertices belonging to handles, #HV x1
Eigen::VectorXi handle_vertices(0, 1);
//centroids of handle regions, #H x1
Eigen::MatrixXd handle_centroids(0, 3);
//updated positions of handle vertices, #HV x3
Eigen::MatrixXd handle_vertex_positions(0, 3);
//index of handle being moved
int moving_handle = -1;
//rotation and translation for the handle being moved
Eigen::Vector3f translation(0, 0, 0);
Eigen::Vector4f rotation(0, 0, 0, 1.);
typedef Eigen::Triplet<double> T;
//per vertex color array, #V x3
Eigen::MatrixXd vertex_colors;

Eigen::SparseMatrix<double> A, mAfc, GDG, GD;
Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>, Eigen::RowMajor> solver;
Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>, Eigen::RowMajor> transferSolver;

Eigen::VectorXi free_indices;
std::vector<std::vector<long>> VV;
std::vector<long> longest_edge;
Eigen::MatrixXd S, D, N, B, B2, S2;
int showView = 0;

//When false, use standard displacement vectors for details, when true use Deformation Transfer from part 2
bool use_deformation_transfer = false;
// if true, solve will be called in the next pre-draw call
bool needs_solve = false;

//function declarations (see below for implementation)
bool solve(Viewer &viewer);

void get_new_handle_locations();

Eigen::Vector3f
computeTranslation(Viewer &viewer, int mouse_x, int from_x, int mouse_y, int from_y, Eigen::RowVector3d pt3D);

Eigen::Vector4f
computeRotation(Viewer &viewer, int mouse_x, int from_x, int mouse_y, int from_y, Eigen::RowVector3d pt3D);

void compute_handle_centroids();

Eigen::MatrixXd readMatrix(const char *filename);

bool callback_mouse_down(Viewer &viewer, int button, int modifier);

bool callback_mouse_move(Viewer &viewer, int mouse_x, int mouse_y);

bool callback_mouse_up(Viewer &viewer, int button, int modifier);

bool callback_pre_draw(Viewer &viewer);

bool callback_key_down(Viewer &viewer, unsigned char key, int modifiers);

void onNewHandleID();

void applySelection();

void calculateA();
void calculateTransferMatrices();

void findNewPositions();

void prefactorise();

void calcDisplacement();

void addLocalDetail();

bool solve(Viewer &viewer) {
    igl::slice_into(handle_vertex_positions, handle_vertices, 1, V);

    findNewPositions();
    B2 = V;
    addLocalDetail();
    S2 = V;

    needs_solve = false;
    return true;
};

void get_new_handle_locations() {
    int count = 0;
    for (long vi = 0; vi < V.rows(); ++vi) {
        if (handle_id[vi] >= 0) {
            Eigen::RowVector3f goalPosition = V.row(vi).cast<float>();
            if (handle_id[vi] == moving_handle) {
                if (mouse_mode == TRANSLATE)
                    goalPosition += translation.transpose();
                else if (mouse_mode == ROTATE) {
                    Eigen::RowVector3f goalPositionCopy = goalPosition;
                    goalPosition -= handle_centroids.row(moving_handle).cast<float>();
                    igl::rotate_by_quat(goalPosition.data(), rotation.data(), goalPositionCopy.data());
                    goalPosition = goalPositionCopy;
                    goalPosition += handle_centroids.row(moving_handle).cast<float>();
                }
            }
            handle_vertex_positions.row(count++) = goalPosition.cast<double>();
        }
    }
}

bool load_mesh(string filename) {
    igl::read_triangle_mesh(filename, V, F);
    viewer.data().clear();
    viewer.data().set_mesh(V, F);

    viewer.core().align_camera_center(V);
    V_cp = V;
    handle_id.setConstant(V.rows(), 1, -1);
    // Initialize selector
    lasso = std::unique_ptr<Lasso>(new Lasso(V, F, viewer));

    selected_v.resize(0, 1);
    calculateA();

    return true;
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        cout << "Usage assignment5 mesh.off>" << endl;
        load_mesh("../data/woody-lo.off");
    } else {
        load_mesh(argv[1]);
    }

    igl::opengl::glfw::imgui::ImGuiPlugin plugin;
    viewer.plugins.push_back(&plugin);
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    plugin.widgets.push_back(&menu);

    menu.callback_draw_viewer_menu = [&]() {
        // Draw parent menu content
        menu.draw_viewer_menu();

        // Add new group
        if (ImGui::CollapsingHeader("Deformation Controls", ImGuiTreeNodeFlags_DefaultOpen)) {
            int mouse_mode_type = static_cast<int>(mouse_mode);

            if (ImGui::Combo("Mouse Mode", &mouse_mode_type, "SELECT\0TRANSLATE\0ROTATE\0NONE\0")) {
                mouse_mode = static_cast<MouseMode>(mouse_mode_type);
            }

            if (ImGui::Button("Clear Selection", ImVec2(-1, 0))) {
                selected_v.resize(0, 1);
            }

            if (ImGui::Button("Apply Selection", ImVec2(-1, 0))) {
                applySelection();
            }

            if (ImGui::Button("Clear Constraints", ImVec2(-1, 0))) {
                handle_id.setConstant(V.rows(), 1, -1);
            }
            if (ImGui::Checkbox("Deformation Transfer", &use_deformation_transfer)) {
                calculateTransferMatrices();
            }

            std::vector<std::string> views = {"None", "S", "B", "B'", "S'"};
            ImGui::Combo("View", &showView, views);
        }
    };

    viewer.callback_key_down = callback_key_down;
    viewer.callback_mouse_down = callback_mouse_down;
    viewer.callback_mouse_move = callback_mouse_move;
    viewer.callback_mouse_up = callback_mouse_up;
    viewer.callback_pre_draw = callback_pre_draw;

    viewer.data().point_size = 10;
    viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);
    viewer.launch();
}


bool callback_mouse_down(Viewer &viewer, int button, int modifier) {
    if (button == (int) Viewer::MouseButton::Right)
        return false;

    down_mouse_x = viewer.current_mouse_x;
    down_mouse_y = viewer.current_mouse_y;

    if (mouse_mode == SELECT) {
        if (lasso->strokeAdd(viewer.current_mouse_x, viewer.current_mouse_y) >= 0)
            doit = true;
        else
            lasso->strokeReset();
    } else if ((mouse_mode == TRANSLATE) || (mouse_mode == ROTATE)) {
        int vi = lasso->pickVertex(viewer.current_mouse_x, viewer.current_mouse_y);
        if (vi >= 0 && handle_id[vi] >= 0)  //if a region was found, mark it for translation/rotation
        {
            moving_handle = handle_id[vi];
            get_new_handle_locations();
            doit = true;
        }
    }
    return doit;
}

bool callback_mouse_move(Viewer &viewer, int mouse_x, int mouse_y) {
    if (!doit)
        return false;
    if (mouse_mode == SELECT) {
        lasso->strokeAdd(mouse_x, mouse_y);
        return true;
    }
    if ((mouse_mode == TRANSLATE) || (mouse_mode == ROTATE)) {
        if (mouse_mode == TRANSLATE) {
            translation = computeTranslation(viewer,
                                             mouse_x,
                                             down_mouse_x,
                                             mouse_y,
                                             down_mouse_y,
                                             handle_centroids.row(moving_handle));
        } else {
            rotation = computeRotation(viewer,
                                       mouse_x,
                                       down_mouse_x,
                                       mouse_y,
                                       down_mouse_y,
                                       handle_centroids.row(moving_handle));
        }
        get_new_handle_locations();
#ifndef UPDATE_ONLY_ON_UP
        needs_solve = true;
        down_mouse_x = mouse_x;
        down_mouse_y = mouse_y;
#endif
        return true;

    }
    return false;
}

bool callback_mouse_up(Viewer &viewer, int button, int modifier) {
    if (!doit)
        return false;
    doit = false;
    if (mouse_mode == SELECT) {
        selected_v.resize(0, 1);
        lasso->strokeFinish(selected_v);
        return true;
    }

    if ((mouse_mode == TRANSLATE) || (mouse_mode == ROTATE)) {
#ifdef UPDATE_ONLY_ON_UP
        if(moving_handle>=0)
          solve(viewer);
#endif
        translation.setZero();
        rotation.setZero();
        rotation[3] = 1.;
        moving_handle = -1;

        compute_handle_centroids();

        return true;
    }

    return false;
};

bool callback_pre_draw(Viewer &viewer) {
#ifndef UPDATE_ONLY_ON_UP
    if (needs_solve) {
        solve(viewer);
    }
#endif
    // initialize vertex colors
    vertex_colors = Eigen::MatrixXd::Constant(V.rows(), 3, .9);

    // first, color constraints
    for (int i = 0; i < V.rows(); ++i)
        if (handle_id[i] != -1) {
            int r = handle_id[i] % MAXNUMREGIONS;
            vertex_colors.row(i) << regionColors[r][0], regionColors[r][1], regionColors[r][2];
        }
    // then, color selection
    for (int i = 0; i < selected_v.size(); ++i)
        vertex_colors.row(selected_v[i]) << 131. / 255, 131. / 255, 131. / 255.;

    viewer.data().set_colors(vertex_colors);

    //clear points and lines
    viewer.data().set_points(Eigen::MatrixXd::Zero(0, 3), Eigen::MatrixXd::Zero(0, 3));
    viewer.data().set_edges(Eigen::MatrixXd::Zero(0, 3), Eigen::MatrixXi::Zero(0, 3), Eigen::MatrixXd::Zero(0, 3));

    //draw the stroke of the selection
    for (unsigned int i = 0; i < lasso->strokePoints.size(); ++i) {
        viewer.data().add_points(lasso->strokePoints[i], Eigen::RowVector3d(0.4, 0.4, 0.4));
        if (i > 1)
            viewer.data().add_edges(lasso->strokePoints[i - 1], lasso->strokePoints[i],
                                    Eigen::RowVector3d(0.7, 0.7, 0.7));
    }

    //update the vertex position all the time
    switch (showView) {
        case 0: viewer.data().set_mesh(V, F); break;
        case 1: viewer.data().set_mesh(S, F); break;
        case 2: viewer.data().set_mesh(B, F); break;
        case 3: viewer.data().set_mesh(B2, F); break;
        case 4: viewer.data().set_mesh(S2, F); break;
    }

#ifdef UPDATE_ONLY_ON_UP
    //draw only the moving parts with a white line
    if (moving_handle>=0)
    {
      Eigen::MatrixXd edges(3*F.rows(),6);
      int num_edges = 0;
      for (int fi = 0; fi<F.rows(); ++fi)
      {
        int firstPickedVertex = -1;
        for(int vi = 0; vi<3 ; ++vi)
          if (handle_id[F(fi,vi)] == moving_handle)
          {
            firstPickedVertex = vi;
            break;
          }
        if(firstPickedVertex==-1)
          continue;


        Eigen::Matrix3d points;
        for(int vi = 0; vi<3; ++vi)
        {
          int vertex_id = F(fi,vi);
          if (handle_id[vertex_id] == moving_handle)
          {
            int index = -1;
            // if face is already constrained, find index in the constraints
            (handle_vertices.array()-vertex_id).cwiseAbs().minCoeff(&index);
            points.row(vi) = handle_vertex_positions.row(index);
          }
          else
            points.row(vi) =  V.row(vertex_id);

        }
        edges.row(num_edges++) << points.row(0), points.row(1);
        edges.row(num_edges++) << points.row(1), points.row(2);
        edges.row(num_edges++) << points.row(2), points.row(0);
      }
      edges.conservativeResize(num_edges, Eigen::NoChange);
      viewer.data().add_edges(edges.leftCols(3), edges.rightCols(3), Eigen::RowVector3d(0.9,0.9,0.9));
    }
#endif
    return false;

}

bool callback_key_down(Viewer &viewer, unsigned char key, int modifiers) {
    bool handled = false;
    if (key == 'S') {
        mouse_mode = SELECT;
        handled = true;
    }

    if ((key == 'T') && (modifiers == IGL_MOD_ALT)) {
        mouse_mode = TRANSLATE;
        handled = true;
    }

    if ((key == 'R') && (modifiers == IGL_MOD_ALT)) {
        mouse_mode = ROTATE;
        handled = true;
    }
    if (key == 'A') {
        applySelection();

        callback_key_down(viewer, '1', 0);
        handled = true;
    }

    return handled;
}

void onNewHandleID() {
    //store handle vertices too
    int numFree = (handle_id.array() == -1).cast<int>().sum();
    int num_handle_vertices = V.rows() - numFree;
    handle_vertices.setZero(num_handle_vertices);
    handle_vertex_positions.setZero(num_handle_vertices, 3);

    int count = 0;
    for (long vi = 0; vi < V.rows(); ++vi)
        if (handle_id[vi] >= 0)
            handle_vertices[count++] = vi;

    compute_handle_centroids();
}

void applySelection() {
    int index = handle_id.maxCoeff() + 1;
    for (int i = 0; i < selected_v.rows(); ++i) {
        const int selected_vertex = selected_v[i];
        if (handle_id[selected_vertex] == -1)
            handle_id[selected_vertex] = index;
    }
    selected_v.resize(0, 1);

    onNewHandleID();

    // precalculate everything that possible before transformation
    prefactorise();
    S = V;
    findNewPositions();
    B = V;
    calcDisplacement();
    V = S;
}

void compute_handle_centroids() {
    //compute centroids of handles
    int num_handles = handle_id.maxCoeff() + 1;
    handle_centroids.setZero(num_handles, 3);

    Eigen::VectorXi num;
    num.setZero(num_handles, 1);
    for (long vi = 0; vi < V.rows(); ++vi) {
        int r = handle_id[vi];
        if (r != -1) {
            handle_centroids.row(r) += V.row(vi);
            num[r]++;
        }
    }

    for (long i = 0; i < num_handles; ++i)
        handle_centroids.row(i) = handle_centroids.row(i).array() / num[i];

}

//computes translation for the vertices of the moving handle based on the mouse motion
Eigen::Vector3f computeTranslation(Viewer &viewer,
                                   int mouse_x,
                                   int from_x,
                                   int mouse_y,
                                   int from_y,
                                   Eigen::RowVector3d pt3D) {
    Eigen::Matrix4f modelview = viewer.core().view;// * viewer.data().model;
    //project the given point (typically the handle centroid) to get a screen space depth
    Eigen::Vector3f proj = igl::project(pt3D.transpose().cast<float>().eval(),
                                        modelview,
                                        viewer.core().proj,
                                        viewer.core().viewport);
    float depth = proj[2];

    double x, y;
    Eigen::Vector3f pos1, pos0;

    //unproject from- and to- points
    x = mouse_x;
    y = viewer.core().viewport(3) - mouse_y;
    pos1 = igl::unproject(Eigen::Vector3f(x, y, depth),
                          modelview,
                          viewer.core().proj,
                          viewer.core().viewport);


    x = from_x;
    y = viewer.core().viewport(3) - from_y;
    pos0 = igl::unproject(Eigen::Vector3f(x, y, depth),
                          modelview,
                          viewer.core().proj,
                          viewer.core().viewport);

    //translation is the vector connecting the two
    Eigen::Vector3f translation = pos1 - pos0;
    return translation;

}

//computes translation for the vertices of the moving handle based on the mouse motion
Eigen::Vector4f computeRotation(Viewer &viewer,
                                int mouse_x,
                                int from_x,
                                int mouse_y,
                                int from_y,
                                Eigen::RowVector3d pt3D) {

    Eigen::Vector4f rotation;
    rotation.setZero();
    rotation[3] = 1.;

    Eigen::Matrix4f modelview = viewer.core().view;// * viewer.data().model;

    //initialize a trackball around the handle that is being rotated
    //the trackball has (approximately) width w and height h
    double w = viewer.core().viewport[2] / 8;
    double h = viewer.core().viewport[3] / 8;

    //the mouse motion has to be expressed with respect to its center of mass
    //(i.e. it should approximately fall inside the region of the trackball)

    //project the given point on the handle(centroid)
    Eigen::Vector3f proj = igl::project(pt3D.transpose().cast<float>().eval(),
                                        modelview,
                                        viewer.core().proj,
                                        viewer.core().viewport);
    proj[1] = viewer.core().viewport[3] - proj[1];

    //express the mouse points w.r.t the centroid
    from_x -= proj[0];
    mouse_x -= proj[0];
    from_y -= proj[1];
    mouse_y -= proj[1];

    //shift so that the range is from 0-w and 0-h respectively (similarly to a standard viewport)
    from_x += w / 2;
    mouse_x += w / 2;
    from_y += h / 2;
    mouse_y += h / 2;

    //get rotation from trackball
    Eigen::Vector4f drot = viewer.core().trackball_angle.coeffs();
    Eigen::Vector4f drot_conj;
    igl::quat_conjugate(drot.data(), drot_conj.data());
    igl::trackball(w, h, float(1.), rotation.data(), from_x, from_y, mouse_x, mouse_y, rotation.data());

    //account for the modelview rotation: pre-rotate by modelview (place model back to the original
    //unrotated frame), post-rotate by inverse modelview
    Eigen::Vector4f out;
    igl::quat_mult(rotation.data(), drot.data(), out.data());
    igl::quat_mult(drot_conj.data(), out.data(), rotation.data());
    return rotation;
}

void calculateA() {
    igl::adjacency_list(F, VV);

    Eigen::SparseVector<double> Asum;
    Eigen::SparseMatrix<double> Adiag, L, M;
    igl::cotmatrix(V, F, L);
    igl::sum(L, 1, Asum);
    igl::diag(Asum, Adiag);
    L = Adiag - L;
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M);
    A =  L * M.cwiseInverse() * L;
}

void calculateTransferMatrices() {
    igl::HeatGeodesicsData<double> data;
    igl::heat_geodesics_precompute(V, F, data);
    GD = data.Div;
    GDG = GD * data.Grad;
    transferSolver.compute(GDG);

    prefactorise();
    S = V;
    findNewPositions();
    B = V;
    calcDisplacement();
    V = S;
}

void prefactorise() {
    long constr_num = handle_vertices.size();
    long free_num = V.rows() - constr_num;
    free_indices.setZero(free_num);
    long index = 0;
    for (long vi = 0; vi < V.rows(); vi++) {
        if (handle_id[vi] != -1) {
            continue;
        }
        free_indices[index] = vi;
        index++;
    }
    Eigen::SparseMatrix<double> Aff = Eigen::SparseMatrix<double>(free_num, free_num);
    mAfc = Eigen::SparseMatrix<double>(free_num, constr_num);
    Eigen::SparseMatrix<double> Af = Eigen::SparseMatrix<double>(free_num, V.rows());
    igl::slice(A, free_indices, 1, Af);
    igl::slice(Af, free_indices, 2, Aff);
    igl::slice(Af, handle_vertices, 2, mAfc);

    solver.compute(Aff);
}

void findNewPositions() {
    Eigen::MatrixXd vc = Eigen::MatrixXd(handle_vertices.size(), 3);
    igl::slice(V, handle_vertices, 1, vc);
    auto b = -1 * mAfc * vc;
    Eigen::MatrixXd Vf = solver.solve(b);
    for (long index = 0; index < free_indices.size(); index++) {
        long V_row = free_indices[index];
        V.row(V_row) = Vf.row(index);
    }
}

void calcDisplacement() {
    if (use_deformation_transfer) {
        igl::per_face_normals(V, F, N);
        D = Eigen::MatrixXd(F.rows() * 3, 3);
        for (long ti = 0; ti < F.rows(); ti++) {
            auto face = F.row(ti);
            auto ni = N.row(ti);
            auto q1 = V.row(face[0]);
            auto q2 = V.row(face[1]);
            auto q3 = V.row(face[2]);
            Eigen::Matrix3d displB, displBt;
            displB << q1-q3, q2-q3, ni;
            displBt = displB.transpose();
            D.row(ti*3) = displBt.row(0);
            D.row(ti*3+1) = displBt.row(1);
            D.row(ti*3+2) = displBt.row(2);
        }
        return;
    }
    longest_edge = {};
    long free_num = free_indices.size();
    igl::per_vertex_normals(V, F, N);
    D = Eigen::MatrixXd(free_num, 3);
    for (long i = 0; i < free_num; i++) {
        long V_row = free_indices[i];
        Eigen::Vector3d vb = V.row(V_row);
        Eigen::Vector3d si = S.row(V_row);
        Eigen::Vector3d ni = N.row(V_row).normalized();
        Eigen::Vector3d di = si - vb;
        double maxDist = 0;
        long maxIndex;
        Eigen::Vector3d xi;
        for (long nV: VV[V_row]) {
            Eigen::Vector3d edgeP = V.row(nV);
            Eigen::Vector3d edge = edgeP - vb;
            auto edgeProject = edge - ni.dot(edge) * ni;
            double dist = edgeProject.norm();
            if (dist > maxDist) {
                maxDist = dist;
                xi = edgeProject.normalized();
                maxIndex = nV;
            }
        }
        longest_edge.push_back(maxIndex);
        Eigen::Vector3d yi = ni.cross(xi).normalized();
        Eigen::Matrix3d localBasis;
        localBasis << xi, yi, ni;
        D.row(i) = localBasis.inverse() * di;
    }
}

void addLocalDetail() {
    if (use_deformation_transfer) {
        igl::per_face_normals(V, F, N);
        Eigen::MatrixXd D2 = Eigen::MatrixXd(F.rows() * 3, 3);
        for (long ti = 0; ti < F.rows(); ti++) {
            auto face = F.row(ti);
            auto ni = N.row(ti);
            auto q1 = V.row(face[0]);
            auto q2 = V.row(face[1]);
            auto q3 = V.row(face[2]);

            Eigen::Matrix3d displB, displB2, diff, diffT;
            displB2 << q1-q3, q2-q3, ni;
            displB << D.row(ti*3), D.row(ti*3+1), D.row(ti*3+2);

            diff = displB2.transpose() * displB.inverse();
            diffT = diff.transpose();
            D2.row(ti*3) = diffT.row(0);
            D2.row(ti*3+1) = diffT.row(1);
            D2.row(ti*3+2) = diffT.row(2);
        }

        auto b = GD * D2;
        V = transferSolver.solve(b);
        return;
    }
    long free_num = free_indices.size();
    igl::per_vertex_normals(V, F, N);
    Eigen::MatrixXd D2 = Eigen::MatrixXd(free_num, 3);
    for (long i = 0; i < free_num; i++) {
        long V_row = free_indices[i];
        Eigen::Vector3d vb = V.row(V_row);
        Eigen::Vector3d ni = N.row(V_row).normalized();
        Eigen::Vector3d edgeP = V.row(longest_edge[i]);
        Eigen::Vector3d edge = edgeP - vb;
        auto edgeProject = edge - ni.dot(edge) * ni;
        Eigen::Vector3d xi = edgeProject.normalized();
        Eigen::Vector3d yi = ni.cross(xi).normalized();
        Eigen::Vector3d di = D.row(i);
        Eigen::Matrix3d localBasis;
        localBasis << xi, yi, ni;
        D2.row(i) = localBasis * di;
    }

    for (long i = 0; i < free_num; i++) {
        long V_row = free_indices[i];
        V.row(V_row) += D2.row(i);
    }
}
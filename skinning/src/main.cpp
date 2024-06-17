#include <iostream>
#include <igl/readOFF.h>
#include <igl/readTGF.h>
#include <igl/readDMAT.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <unsupported/Eigen/MatrixFunctions>
#include <imgui.h>

#include <igl/forward_kinematics.h>
#include <igl/directed_edge_parents.h>
#include <igl/deform_skeleton.h>
#include <igl/cotmatrix.h>
#include <igl/diag.h>
#include <igl/sum.h>
#include <igl/cat.h>
#include <igl/dqs.h>
#include <igl/lbs_matrix.h>
#include <igl/grad.h>

using namespace std;
using Viewer = igl::opengl::glfw::Viewer;
Viewer mViewer;
const std::string folder = "../data/hand/";

// Mesh
// Vertex array, #V x3
Eigen::MatrixXd V;
// Face array, #F x3
Eigen::MatrixXi F;
// Normals
Eigen::MatrixXd N;

// Skeleton
// Joints, #J x3
Eigen::MatrixXd J;
// Bones, #B x3
Eigen::MatrixXi B;
// Loaded (given) and constructed handles
Eigen::VectorXi LH, CH;
// Handles for each bone (better format to visualize them)
std::vector<std::vector<int>> load_handles, constr_handles;

// Animations
// Rotation quaternions: RQ - loaded as quaternions; RM - loaded matrices and converted to quaternions
Eigen::MatrixXd RQ;
Eigen::MatrixXd rot_matrices;
std::vector<std::vector<Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond>>> RM;
// direct edge parent
Eigen::VectorXi PI;

// Global transformations
Eigen::MatrixXd GGT;
std::vector<Eigen::Vector3d> GT;
std::vector<Eigen::Quaterniond,Eigen::aligned_allocator<Eigen::Quaterniond>> GQ;
Eigen::MatrixXd load_weights, constr_weights, face_load_weights, face_constr_weights;
std::vector<Eigen::Matrix4d> logQ;

// colors
Eigen::MatrixXd colorP, colorE;
Eigen::RowVector3d red, blue, background;

// per-face solver
Eigen::SparseMatrix<double> GDG, GD, GDGfc;
Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>, Eigen::RowMajor> face_solver;
Eigen::VectorXi constrV, freeV;
Eigen::MatrixXd constrV_pos;

// Animations parameters
double anim_t = 0.0;
double anim_t_dir = 0.015;
double anim_seconds = 2;

int quaternion_mode = true;
bool show_handles = false;
bool show_weights = false;
bool given_weights = false;
bool without_stitching = false;
int skinning = 2;
int handle_num = 0;
float scale = 2.0;

void set_colors() {
    Eigen::RowVector3d green;
    red << 1, 0, 0;
    green << 0, 1, 0;
    blue << 0, 0, 1;
    background << 1, 1, 0;

    colorP = red.replicate(J.rows(), 1);
    colorE = green.replicate(B.rows(), 1);
}

bool load_mesh(const string &filename) {
    if (!igl::readOFF(filename,V,F)) {
        return false;
    }
    igl::per_face_normals(V, F, N);
    mViewer.data().clear();
    mViewer.data().set_mesh(V,F);
    mViewer.data().compute_normals();
    mViewer.core().align_camera_center(V, F);
    return true;
}

bool load_skeleton(const string &filename) {
    if (!igl::readTGF(filename, J, B)) {
        return false;
    }

    set_colors();
    mViewer.data().set_points(J, colorP);
    mViewer.data().set_edges(J, B, colorE);
    mViewer.core().align_camera_center(J, B);
    return true;
}

void quaternion_global_calc(double t) {
    Eigen::MatrixXd NJ = Eigen::MatrixXd(J.rows(), J.cols());

    // Interpolate
    if (quaternion_mode) {
        std::vector<Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond>> anim_pose(B.rows());
        for (int e = 0; e < B.rows(); e++) {
            Eigen::Quaterniond start = Eigen::Quaterniond::Identity();
            Eigen::Vector4d endVector = RQ.block(e * 4, 0, 4, 1);
            Eigen::Quaterniond end = Eigen::Quaterniond(endVector);
            anim_pose[e] = start.slerp(t, end);
        }

        igl::forward_kinematics(J, B, PI, anim_pose, GQ, GT);
    } else {
        long currentFrame = std::floor(t * RM.size());
        // Converted matrices to quaternions on the start
        igl::forward_kinematics(J, B, PI, RM[currentFrame], GQ, GT);
    }

    GGT = Eigen::MatrixXd(B.rows()*4, 3);
    logQ.resize(F.rows());
    for(int e = 0;e < B.rows(); e++) {
        Eigen::Affine3d a = Eigen::Affine3d::Identity();
        a.translate(GT[e]);
        a.rotate(GQ[e]);
        logQ[e] = a.matrix().log();
        GGT.block(e*4, 0, 4, 3) =
                a.matrix().transpose().block(0, 0, 4, 3);
    }

    Eigen::MatrixXi NB;
    igl::deform_skeleton(J, B, GGT, NJ, NB);

    mViewer.data().set_points(NJ, colorP);
    mViewer.data().set_edges(NJ, NB, colorE);
}

bool callback_pre_draw(Viewer &viewer) {
    const auto &weights = given_weights ? load_weights : constr_weights;
    const auto &face_weights = given_weights ? face_load_weights : face_constr_weights;
    const auto &handles = given_weights ? load_handles : constr_handles;
    if (viewer.core().is_animating) {
        // Find pose interval
        const double t = anim_t / anim_seconds;

        Eigen::MatrixXd NV = V;
        quaternion_global_calc(t);
        if (skinning == 2) {
            Eigen::MatrixXd ST = Eigen::MatrixXd(F.rows() * 3, 3);
            Eigen::MatrixXd verts = Eigen::MatrixXd(F.rows() * 3, 3);
            Eigen::MatrixXi new_faces = Eigen::MatrixXi(F.rows(), 3);
            Eigen::MatrixXd D = Eigen::MatrixXd(F.rows() * 3, 3);
            for (long ti = 0; ti < F.rows(); ti++) {
                auto face = F.row(ti);
                auto ni = N.row(ti);
                auto q1 = V.row(face[0]);
                auto q2 = V.row(face[1]);
                auto q3 = V.row(face[2]);
                Eigen::VectorXd fW = face_weights.row(ti);
                Eigen::Matrix4d sum = Eigen::Matrix4d::Zero();
                for(int e = 0; e < B.rows(); e++) {
                    sum += fW[e] * logQ[e];
                }
                sum = sum.exp();
                Eigen::Vector3d nv1 = sum.block(0, 0, 3, 3) * Eigen::Vector3d(q1);
                Eigen::Vector3d nv2 = sum.block(0, 0, 3, 3) * Eigen::Vector3d(q2);
                Eigen::Vector3d nv3 = sum.block(0, 0, 3, 3) * Eigen::Vector3d(q3);

                verts.row(ti*3) = nv1;
                verts.row(ti*3+1) = nv2;
                verts.row(ti*3+2) = nv3;
                new_faces.row(ti) << ti*3, ti*3+1, ti*3+2;

                Eigen::Matrix3d displB, displBt;
                displB << q1-q3, q2-q3, ni;
                displBt = displB.transpose();
                D.row(ti*3) = displBt.row(0);
                D.row(ti*3+1) = displBt.row(1);
                D.row(ti*3+2) = displBt.row(2);
            }
            Eigen::MatrixXd new_norms;
            igl::per_face_normals(verts, new_faces, new_norms);
            for (long ti = 0; ti < F.rows(); ti++) {
                auto face = F.row(ti);
                auto ni = new_norms.row(ti);
                auto v1 = verts.row(ti*3);
                auto v2 = verts.row(ti*3+1);
                auto v3 = verts.row(ti*3+2);

                Eigen::Matrix3d displB, displB2, diff, diffT;
                displB2 << v1-v3, v2-v3, ni;
                displB << D.row(ti*3), D.row(ti*3+1), D.row(ti*3+2);

                diff = displB2.transpose() * displB.inverse();
                diffT = diff.transpose();
                ST.row(ti) = diffT.row(0);
                ST.row(ti+F.rows()) = diffT.row(1);
                ST.row(ti+2*F.rows()) = diffT.row(2);
            }
            if (without_stitching) {
                mViewer.data().set_points(verts, colorP);
            }
            Eigen::SparseMatrix<double> b = (GD * ST).sparseView();
            Eigen::SparseMatrix<double> free_b;
            igl::slice(b, freeV, 1, free_b);
            free_b -= (GDGfc*constrV_pos).sparseView();

            Eigen::MatrixXd V_free = face_solver.solve(free_b);
            for (long i = 0; i < freeV.size(); i++) {
                long V_row = freeV[i];
                NV.row(V_row) = V_free.row(i);
            }
        } else if (skinning == 1) {
            igl::dqs(V, weights, GQ, GT, NV);
        } else {
            Eigen::MatrixXd M;
            igl::lbs_matrix(V, weights, M);
            NV = M * GGT;
        }

        mViewer.data().set_mesh(NV, F);
        mViewer.data().compute_normals();
        mViewer.core().align_camera_center(V, F);
        double nextTime = anim_t + anim_t_dir;
        if (nextTime > anim_seconds || nextTime < 0) {
            anim_t_dir *= -1;
        }
        anim_t += anim_t_dir;
    }
    Eigen::MatrixXd vertex_colors = background.replicate(V.rows(), 1);
    if (show_weights || show_handles) {
        int e = handle_num;
        if (e >= B.rows()) {
            e = 0;
        }
        if (show_weights) {
            Eigen::VectorXd w = weights.col(e);
            for (int v = 0; v < V.rows(); v++) {
                double vw = w[v];
                vertex_colors.row(v) << vw, 0., 1. - vw;
            }
        } else if (show_handles) {
            vertex_colors = blue.replicate(V.rows(), 1);
            for (long v: handles[e]) {
                vertex_colors.row(v) = red;
            }
        }
    }
    mViewer.data().set_colors(vertex_colors);
    return false;
}

bool callback_key_down(Viewer &viewer, unsigned char key, int)
{
    // Change modes: quaternions -> matrices; matrices -> quaternions
    if (key == 'M') {
        quaternion_mode = !bool(quaternion_mode);
    }

    return true;
}

double distancePointSegment(const Eigen::RowVector3d& vert, const Eigen::RowVector3d& j_start,
                            const Eigen::RowVector3d& j_end)
{
    Eigen::RowVector3d segment = j_end - j_start;
    Eigen::RowVector3d toPoint = vert - j_start;

    double t = toPoint.dot(segment) / segment.dot(segment);

    if (t < 0.0) {
        return (vert - j_start).norm();
    } else if (t > 1.0) {
        return (vert - j_end).norm();
    }
    Eigen::RowVector3d closestPoint = j_start + t * segment;
    return (vert - closestPoint).norm();
}

void construct_handles() {
    std::vector<Eigen::VectorXd> squaredDistances(B.rows());
    std::vector<long> closestVertexIndices(B.rows());
    std::set<long> alreadyAdded;
    CH = Eigen::VectorXi::Constant(V.rows(),-1);

    for (long e = 0; e < B.rows(); e++) {
        Eigen::RowVector3d j_start = J.row(B(e, 0));
        Eigen::RowVector3d j_end = J.row(B(e, 1));

        std::vector<double> distances;
        double minDist = DBL_MAX;
        for (long v = 0; v < V.rows(); v++) {
            Eigen::RowVector3d vert = V.row(v);
            double dist = distancePointSegment(vert, j_start, j_end);
            distances.push_back(dist);
            if (dist < minDist) {
                minDist = dist;
            }
        }
        double threshold = minDist * scale;
        for (long v = 0; v < V.rows(); v++) {
            double dist = distances[v];
            if (dist < threshold) {
                CH[v] = e;
            }
        }
    }
}

void prefactorise() {
    int root = -1;
    int constr_size = 0;
    for (int i = 0; i < PI.rows(); i++) {
        if (PI[i] == -1) {
            root = i;
            constr_size = load_handles[i].size();
            break;
        }
    }
    constrV = Eigen::VectorXi(constr_size);
    freeV = Eigen::VectorXi(V.rows() - constr_size);
    constrV_pos = Eigen::MatrixXd(constr_size, 3);
    int constr_Num = 0; int free_Num = 0;
    for (int v = 0; v < V.rows(); v++) {
        auto it = std::find(load_handles[root].begin(), load_handles[root].end(), v);
        if (it != load_handles[root].end()) {
            constrV[constr_Num] = v;
            constrV_pos.row(constr_Num++) = V.row(v);
        } else {
            freeV[free_Num++] = v;
        }
    }

    Eigen::VectorXd dblA;
    Eigen::SparseMatrix<double> Grad,Div;
    igl::doublearea(V, F, dblA);
    igl::grad(V,F,Grad);
    GD = Grad.transpose() * dblA.colwise().replicate(3).asDiagonal();
    GDG = GD * Grad;
    Eigen::SparseMatrix<double> GDGff = Eigen::SparseMatrix<double>(free_Num, free_Num);
    GDGfc = Eigen::SparseMatrix<double>(free_Num, constr_Num);
    igl::slice(GDG, freeV, freeV, GDGff);
    igl::slice(GDG, freeV, constrV, GDGfc);
    face_solver.compute(GDGff);
}

void calculate_weights(const Eigen::VectorXi &H, std::vector<std::vector<int>> &handles, Eigen::MatrixXd &weights,
                       Eigen::MatrixXd &face_weights) {
    handles = {};
    weights = Eigen::MatrixXd(V.rows(), B.rows());
    face_weights = Eigen::MatrixXd(F.rows(), B.rows());
    handles.resize(B.rows());

    Eigen::SparseVector<double> Asum;
    Eigen::SparseMatrix<double> A, Adiag, C;
    Eigen::VectorXd b, d;
    igl::cotmatrix(V, F, A);
    igl::sum(A, 1, Asum);
    igl::diag(Asum, Adiag);
    A = Adiag - A;
    b = Eigen::VectorXd::Zero(V.rows());

    for (int e = 0; e < B.rows(); e++) {
        std::vector<long> constr_zero;
        for (long i = 0; i < H.rows(); i++) {
            long v = H[i];
            if (v == e) {
                handles[e].push_back(i);
            } else if (v != -1) {
                constr_zero.push_back(i);
            }
        }
        long constr_len = handles[e].size() + constr_zero.size();
        C = Eigen::SparseMatrix<double>(constr_len, V.rows());
        d = Eigen::VectorXd::Zero(constr_len);

        for (long i = 0; i < handles[e].size(); i++) {
            auto vInd = handles[e][i];
            C.coeffRef(i, vInd) = 1.0;
            d[i] = 1.0;
        }
        for (long i = 0; i < constr_zero.size(); i++) {
            auto vInd = constr_zero[i];
            C.coeffRef(i+handles[e].size(), vInd) = 1.0;
        }

        Eigen::SparseMatrix<double> col1, col2, Ares;
        Eigen::VectorXd bres;
        igl::cat(1, b, d, bres);
        igl::cat(1, A, C, col1);
        igl::cat(1, Eigen::SparseMatrix<double>(C.transpose()),
                 Eigen::SparseMatrix<double>(C.rows(), C.rows()), col2);
        igl::cat(2, col1, col2, Ares);

        Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
        solver.analyzePattern(Ares);
        solver.factorize(Ares);

        if (solver.info() != Eigen::Success) {
            std::cout << solver.lastErrorMessage() << std::endl;
            continue;
        }
        weights.col(e) = solver.solve(bres).head(V.rows());
    }
    for (int f = 0; f < F.rows(); f++) {
        Eigen::Vector3i face = F.row(f);
        // face's vertices
        auto n1 = weights.row(face(0));
        auto n2 = weights.row(face(1));
        auto n3 = weights.row(face(2));
        face_weights.row(f) = (n1+n2+n3) / 3;
    }
}

void rotmat_to_quaternions() {
    int frameNums = rot_matrices.rows() / (3 * B.rows());
    RM.resize(frameNums);
    for (int f = 0; f < frameNums; f++) {
        RM[f].resize(B.rows());
        for (int e = 0; e < B.rows(); e++) {
            Eigen::Matrix3d rotMatrix = rot_matrices.block((f*B.rows() + e) * 3, 0, 3, 3);
            Eigen::Quaterniond quatRot(rotMatrix);
            RM[f][e] = quatRot;
        }
    }
}

int main(int argc, char *argv[]) {
    load_mesh(folder + "hand.off");
    load_skeleton(folder + "hand.tgf");
    const auto quat_filename = folder + "hand-pose_quat.dmat";
    const auto mat_filename = folder + "hand-pose_matrot.dmat";
    if (!igl::readDMAT(quat_filename, RQ)) {
        std::cout << "ERROR: can't load quaternions" << std::endl;
        return 1;
    }
    if (!igl::readDMAT(mat_filename, rot_matrices)) {
        std::cout << "ERROR: can't load matrices" << std::endl;
        return 1;
    }
    int frameRate = 20;
    int frameNums = std::floor(rot_matrices.rows() / (3 * B.rows()));
    anim_seconds = std::floor(frameNums / frameRate);
    rotmat_to_quaternions();

    igl::directed_edge_parents(B, PI);
    igl::readDMAT(folder + "hand-handles.dmat", LH);
    construct_handles();
    calculate_weights(LH, load_handles, load_weights, face_load_weights);
    calculate_weights(CH, constr_handles, constr_weights, face_constr_weights);

    prefactorise();

    // Attach a menu plugin
    igl::opengl::glfw::imgui::ImGuiPlugin plugin;
    mViewer.plugins.push_back(&plugin);
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    plugin.widgets.push_back(&menu);

    mViewer.callback_key_down = callback_key_down;

    menu.callback_draw_viewer_menu = [&]() {
        // Draw parent menu content
        menu.draw_viewer_menu();

        // Add new group
        if (ImGui::CollapsingHeader("Animation Options", ImGuiTreeNodeFlags_DefaultOpen)) {
            ImGui::Checkbox("Animation On/Off", &mViewer.core().is_animating);
            std::vector<std::string> modes = {"Matrix", "Quaternion"};
            ImGui::Combo("Rotation representation", &quaternion_mode, modes);
            ImGui::Checkbox("With given weights", &given_weights);
            ImGui::Checkbox("Show handles", &show_handles);
            ImGui::Checkbox("Show weights", &show_weights);
            ImGui::Checkbox("Show without stitching", &without_stitching);
            ImGui::InputInt("Handle number", &handle_num, 0, 0);
            if (ImGui::InputFloat("Scale for threshold (handles)", &scale, 0, 0)) {
                construct_handles();
                calculate_weights(CH, constr_handles, constr_weights, face_constr_weights);
                prefactorise();
            }
            std::vector<std::string> sk_types = {"Linear", "DQS", "Per-face"};
            ImGui::Combo("Skinning type", &skinning, sk_types);
        }
    };

    mViewer.data().point_size = 5;
    mViewer.core().is_animating = true;
    mViewer.data().show_overlay_depth = false;
    mViewer.callback_pre_draw = callback_pre_draw;
    mViewer.launch();
}

#include <iostream>
#include <igl/readOFF.h>
#include <igl/readTGF.h>
#include <igl/readDMAT.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui.h>

#include <igl/forward_kinematics.h>
#include <igl/directed_edge_parents.h>
#include <igl/deform_skeleton.h>
#include <igl/cotmatrix.h>
#include <igl/diag.h>
#include <igl/sum.h>
#include <igl/cat.h>

using namespace std;
using Viewer = igl::opengl::glfw::Viewer;
Viewer mViewer;
const std::string folder = "../data/hand/";

// Mesh
// Vertex array, #V x3
Eigen::MatrixXd V;
// Face array, #F x3
Eigen::MatrixXi F;

// Skeleton
// Joints, #J x3
Eigen::MatrixXd J;
// Bones, #B x3
Eigen::MatrixXi B;
Eigen::VectorXi handles;

// Animations
// Rotation quaternions/matrices
Eigen::MatrixXd RQ, RM, R;
// direct edge parent
Eigen::MatrixXi PI;
// Global transformations
Eigen::MatrixXd GP;
std::vector<Eigen::Vector3d> GT;
std::vector<Eigen::Quaterniond,Eigen::aligned_allocator<Eigen::Quaterniond>> GQ;
std::vector<Eigen::VectorXd> weights;

// colors
Eigen::MatrixXd colorP, colorE;

// Animations parameters
double anim_t = 0.0;
double anim_t_dir = 0.015;
double anim_seconds = 2;

int quaternion_mode = true;

bool load_mesh(const string &filename) {
    if (!igl::readOFF(filename,V,F)) {
        return false;
    }
    mViewer.data().clear();
    mViewer.data().set_mesh(V,F);
    mViewer.data().compute_normals();
    mViewer.core().align_camera_center(V, F);
}

bool load_skeleton(const string &filename) {
    if (!igl::readTGF(filename, J, B)) {
        return false;
    }
    Eigen::RowVector3d green;
    green << 1, 0, 0;
    colorP = green.replicate(J.rows(), 1);

    Eigen::RowVector3d red;
    red << 0, 1, 0;
    colorE = red.replicate(B.rows(), 1);

    mViewer.data().set_points(J, colorP);
    mViewer.data().set_edges(J, B, colorE);
    mViewer.core().align_camera_center(J, B);
}

void quaternion_global_calc(std::vector<Eigen::Affine3d> &gT, double t) {
    // Interpolate pose and identity
    std::vector<Eigen::Quaterniond,Eigen::aligned_allocator<Eigen::Quaterniond>> anim_pose(B.rows());
    for(int e = 0; e < B.rows(); e++) {
        Eigen::Quaterniond start = Eigen::Quaterniond::Identity();
        Eigen::Vector4d endVector = R.block(e*4, 0, 4, 1);
        Eigen::Quaterniond end = Eigen::Quaterniond(endVector);
        anim_pose[e] = start.slerp(t, end);
    }

    // Propagate relative rotations via FK to retrieve absolute transformations
    igl::forward_kinematics(J, B, PI, anim_pose, GQ, GT);

    for(int e = 0; e < B.rows(); e++) {
        Eigen::Affine3d a = Eigen::Affine3d::Identity();
        a.translate(GT[e]);
        a.rotate(GQ[e]);
        gT.push_back(a);
    }
}

void to_global(long currentFrame, long currentEdge, Eigen::Affine3d &a) {
    if (PI(currentEdge) == -1) {
        long startRow = (currentFrame*B.rows() + currentEdge) * 3;
        Eigen::Matrix3d rotation = R.block(startRow, 0, 3, 3);

        a.linear() = rotation;
        return;
    }
    long nextEdge = PI(currentEdge);
    to_global(currentFrame, nextEdge, a);

    long startRow = (currentFrame*B.rows() + currentEdge) * 3;
    Eigen::Matrix3d rotation = R.block(startRow, 0, 3, 3);
    Eigen::Affine3d localTransform = Eigen::Affine3d::Identity();

    localTransform.linear() = rotation;
    a = a * localTransform;
}

void matrices_global_calc(std::vector<Eigen::Affine3d> &gT, double t) {
    long frameNum = R.rows() / (3*B.rows());
    long currentFrame = std::floor(t * frameNum);
    for(int e = 0; e < B.rows(); e++) {
        Eigen::Affine3d a = Eigen::Affine3d::Identity();
        to_global(currentFrame, e, a);
        gT.push_back(a);
    }
}

bool callback_pre_draw(Viewer &viewer) {
    if (viewer.core().is_animating) {
        R = quaternion_mode ? RQ : RM;

        // Find pose interval
        const double t = anim_t / anim_seconds;

        // Global transformations
        std::vector<Eigen::Affine3d> GGT;

        if (quaternion_mode) {
            quaternion_global_calc(GGT, t);
        } else {
            matrices_global_calc(GGT, t);
        }

        // Also deform skeleton edges
        Eigen::MatrixXd NJ, NV;
        NV = Eigen::MatrixXd::Zero(V.rows(), V.cols());
        NJ = Eigen::MatrixXd(J.rows(), J.cols());
        for (int e = 0; e < B.rows(); e++) {
            Eigen::VectorXd w = weights[e];
            Eigen::Affine3d a = GGT[e];
            Eigen::Vector3d joint1 = J.row(B(e, 0));
            Eigen::Vector3d joint2 = J.row(B(e, 1));
            NJ.row(B(e, 0)) = a.linear() * joint1 + a.translation();
            NJ.row(B(e, 1)) = a.linear() * joint2 + a.translation();

            for (int v = 0; v < V.rows(); v++) {
                Eigen::Vector3d pos = V.row(v);
                NV.row(v) += w[v] * (a * pos);
            }
        }

        mViewer.data().set_points(NJ, colorP);
        mViewer.data().set_edges(NJ, B, colorE);

        mViewer.data().set_mesh(NV, F);
        mViewer.data().compute_normals();
        mViewer.core().align_camera_center(V, F);

        double nextTime = anim_t + anim_t_dir;
        if (nextTime > anim_seconds || nextTime < 0) {
            anim_t_dir *= -1;
        }
        anim_t += anim_t_dir;
    }
    return false;
}

bool callback_key_down(Viewer &viewer, unsigned char key, int modifiers)
{
    // Change modes: quaternions -> matrices; matrices -> quaternions
    if (key == 'M') {
        quaternion_mode = !bool(quaternion_mode);
    }

    return true;
}

void calculate_weights() {
    weights.resize(B.rows());

    Eigen::SparseVector<double> Asum;
    Eigen::SparseMatrix<double> A, Adiag, C;
    Eigen::VectorXd b, d;
    igl::cotmatrix(V, F, A);
    igl::sum(A, 1, Asum);
    igl::diag(Asum, Adiag);
    A = Adiag - A;
    b = Eigen::VectorXd::Zero(V.rows());

    for (int e = 0; e < B.rows(); e++) {
        std::vector<long> constr_one, constr_zero;
        for (long i = 0; i < handles.rows(); i++) {
            long v = handles[i];
            if (v == e) {
                constr_one.push_back(i);
            } else if (v != -1) {
                constr_zero.push_back(i);
            }
        }
        long constr_len = constr_one.size() + constr_zero.size();
        C = Eigen::SparseMatrix<double>(constr_len, V.rows());
        d = Eigen::VectorXd::Zero(constr_len);

        for (long i = 0; i < constr_one.size(); i++) {
            auto vInd = constr_one[i];
            C.coeffRef(i, vInd) = 1.0;
            d[i] = 1.0;
        }
        for (long i = 0; i < constr_zero.size(); i++) {
            auto vInd = constr_zero[i];
            C.coeffRef(i+constr_one.size(), vInd) = 1.0;
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
        weights[e] = solver.solve(bres);
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
    if (!igl::readDMAT(mat_filename, RM)) {
        std::cout << "ERROR: can't load quaternions" << std::endl;
        return 1;
    }
    int frameRate = 20;
    int frameNums = RM.rows() / (3 * B.rows());
    anim_seconds = frameNums / frameRate;

    igl::directed_edge_parents(B, PI);
    igl::readDMAT(folder + "hand-handles.dmat", handles);
    calculate_weights();

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
            ImGui::Checkbox("Animations On/Off", &mViewer.core().is_animating);

            std::vector<std::string> modes = {"Matrix", "Quaternion"};
            ImGui::Combo("Rotation representation", &quaternion_mode, modes);
        }
    };

    mViewer.data().point_size = 10;
    mViewer.core().is_animating = true;
    mViewer.data().show_overlay_depth = false;
    mViewer.callback_pre_draw = callback_pre_draw;
    mViewer.launch();
}

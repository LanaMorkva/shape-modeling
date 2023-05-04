#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui.h>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <igl/local_basis.h>
#include <igl/grad.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>


/*** insert any necessary libigl headers here ***/
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/harmonic.h>
#include <igl/lscm.h>
#include <igl/adjacency_matrix.h>
#include <igl/sum.h>
#include <igl/diag.h>
#include <igl/speye.h>
#include <igl/repdiag.h>
#include <igl/cat.h>
#include <igl/dijkstra.h>
#include <igl/adjacency_list.h>

using namespace std;
using namespace Eigen;
using Viewer = igl::opengl::glfw::Viewer;

Viewer viewer;

// vertex array, #V x3
Eigen::MatrixXd V;

// face array, #F x3
Eigen::MatrixXi F;

// UV coordinates, #V x2
Eigen::MatrixXd UV;

Eigen::MatrixXd F_colors;

VectorXi fixed_UV_indices;
MatrixXd fixed_UV_positions;

bool showingUV = false;
bool showingFixedVert = false;
bool showingDistortion = false;
bool freeBoundary = false;
int distType = 0;
double TextureResolution = 10;
igl::opengl::ViewerCore temp3D;
igl::opengl::ViewerCore temp2D;

void Redraw()
{
	viewer.data().clear();
    if (showingDistortion) {
        viewer.data().set_mesh(V, F);
        viewer.data().set_colors(F_colors);
        viewer.data().set_face_based(true);
        return;
    }

	if (!showingUV)
	{
		viewer.data().set_mesh(V, F);
		viewer.data().set_face_based(false);

        if(UV.size() != 0)
        {
          viewer.data().set_uv(TextureResolution*UV);
          viewer.data().show_texture = true;
        }
	}
	else
	{
		viewer.data().show_texture = false;
		viewer.data().set_mesh(UV, F);
	}

    if (showingFixedVert) {
        Eigen::MatrixXd colors, points;
        colors.setZero(fixed_UV_indices.rows(), 3);
        points.setZero(fixed_UV_indices.rows(), 3);

        // Build color map
        for (int i = 0; i < fixed_UV_indices.rows(); ++i)
        {
            points.row(i) = V.row(fixed_UV_indices[i]);
            colors(i, 0) = 1;
        }
        viewer.data().point_size = 10;
        viewer.data().add_points(points, colors);
    }
}

bool callback_mouse_move(Viewer &viewer, int mouse_x, int mouse_y)
{
	if (showingUV)
		viewer.mouse_mode = igl::opengl::glfw::Viewer::MouseMode::Translation;
	return false;
}

static void computeSurfaceGradientMatrix(SparseMatrix<double> & D1, SparseMatrix<double> & D2)
{
	MatrixXd F1, F2, F3;
	SparseMatrix<double> DD, Dx, Dy, Dz;

	igl::local_basis(V, F, F1, F2, F3);
	igl::grad(V, F, DD);

	Dx = DD.topLeftCorner(F.rows(), V.rows());
	Dy = DD.block(F.rows(), 0, F.rows(), V.rows());
	Dz = DD.bottomRightCorner(F.rows(), V.rows());

	D1 = F1.col(0).asDiagonal()*Dx + F1.col(1).asDiagonal()*Dy + F1.col(2).asDiagonal()*Dz;
	D2 = F2.col(0).asDiagonal()*Dx + F2.col(1).asDiagonal()*Dy + F2.col(2).asDiagonal()*Dz;
}
static inline void SSVD2x2(const Eigen::Matrix2d& J, Eigen::Matrix2d& U, Eigen::Matrix2d& S, Eigen::Matrix2d& V)
{
	double e = (J(0) + J(3))*0.5;
	double f = (J(0) - J(3))*0.5;
	double g = (J(1) + J(2))*0.5;
	double h = (J(1) - J(2))*0.5;
	double q = sqrt((e*e) + (h*h));
	double r = sqrt((f*f) + (g*g));
	double a1 = atan2(g, f);
	double a2 = atan2(h, e);
	double rho = (a2 - a1)*0.5;
	double phi = (a2 + a1)*0.5;

	S(0) = q + r;
	S(1) = 0;
	S(2) = 0;
	S(3) = q - r;

	double c = cos(phi);
	double s = sin(phi);
	U(0) = c;
	U(1) = s;
	U(2) = -s;
	U(3) = c;

	c = cos(rho);
	s = sin(rho);
	V(0) = c;
	V(1) = -s;
	V(2) = s;
	V(3) = c;
}

void calculateDistortion()
{
    if (UV.size() == 0) {
        std::cout << "No UV, can't calculate the distortion";
        return;
    }
    VectorXd Dist = VectorXd(F.rows());
    SparseMatrix<double> Dx, Dy;
    MatrixXd I = MatrixXd::Identity(2, 2);
    computeSurfaceGradientMatrix(Dx, Dy);
    for (int i = 0; i < F.rows(); i++) {
        Eigen::Matrix2d J;
        J << Dx.row(i)*UV.col(0), Dy.row(i)*UV.col(0),
             Dx.row(i)*UV.col(1), Dy.row(i)*UV.col(1);
        switch (distType) {
        case 0: Dist[i] = (J+J.transpose()-J.trace()*I).squaredNorm(); break;
        case 1: {
            Eigen::Matrix2d u, s, v, r;
            SSVD2x2(J, u, s, v);
            r = u * v.transpose();
            // only determinant > 0 is suitable; change sign by multiplying one row by -1
            if (r.determinant() < 0) {
                r.row(0) *= -1;
            }
            Dist[i] = (J-r).squaredNorm();
        } break;
        case 2: {
            auto val = J.determinant() - 1;
            Dist[i] = val * val;
        }
        break;
        }
    }
    double maxDist = Dist.maxCoeff();
    F_colors = MatrixXd::Ones(F.rows(), 3);
    F_colors.col(1) -= Dist / maxDist;
    F_colors.col(2) -= Dist / maxDist;
}

void ConvertConstraintsToMatrixForm(const VectorXi &indices, const MatrixXd &positions, Eigen::SparseMatrix<double> &C,
                                    VectorXd &d)
{
    long size = indices.size();
    long len = 2 * size;
    // C rows: #of constraints - fixed points; cols: #of vertices*2 (for u and v)
    C = Eigen::SparseMatrix<double>(len, 2*V.rows());
    d = VectorXd::Zero(len);

    for (long i=0; i<size; i++) {
        auto ind = indices[i];
        C.coeffRef(i, ind) = 1.0;
        C.coeffRef(i+size, ind+V.rows()) = 1.0;
        d[i] = positions.coeff(i, 0);
        d[i+size] = positions.coeff(i, 1);
    }
}

void computeParameterization(int type)
{
	SparseMatrix<double> A;
	VectorXd b;
	Eigen::SparseMatrix<double> C;
	VectorXd d;
	// Find the indices of the boundary vertices of the mesh and put them in fixed_UV_indices
	if (!freeBoundary || type < '3') {
        // need to recalculate (important to not recalculate for big meshes and especially for key '4')
        if (fixed_UV_indices.size() <= 2) {
            igl::boundary_loop(F, fixed_UV_indices);
            igl::map_vertices_to_circle(V, fixed_UV_indices, fixed_UV_positions);
        }
	} else {
        // need to recalculate (important to not recalculate for big meshes and especially for key '4')
        if (fixed_UV_indices.size() != 2) {
            VectorXi prev_pos;
            VectorXd minDist;
            std::vector<std::vector<int>> VV;
            igl::adjacency_list(F, VV);
            double maxDist = 0;
            // indices for 2 most distant points
            int ind1, ind2;
            for (int i = 0; i < V.rows(); i++) {
                igl::dijkstra(V, VV, i, {}, minDist, prev_pos);
                for (int j = 0; j < minDist.size(); j++) {
                    if (maxDist < minDist[j]) {
                        maxDist = minDist[j];
                        ind1 = i;
                        ind2 = j;
                    }
                }
            }

            fixed_UV_indices = VectorXi(2);
            fixed_UV_indices << ind1, ind2;
            igl::map_vertices_to_circle(V, fixed_UV_indices, fixed_UV_positions);
        }
    }

	ConvertConstraintsToMatrixForm(fixed_UV_indices, fixed_UV_positions, C, d);

	if (type == '1') {
        auto facesNum = F.rows();
        auto vNum = V.rows();
        b = VectorXd::Zero(2*vNum);

        // duplicate faces (with offset) to receive 2 times bigger matrix A (for u and v)
        Eigen::MatrixXi Fuv = Eigen::MatrixXi(facesNum*2, F.cols());
        for (int i = 0; i < facesNum; i++) {
            Eigen::VectorXi face = F.row(i);
            Fuv.row(i) = face;
            Fuv.row(i+facesNum) = face.array() + vNum;
        }

        SparseVector<double> Asum;
        SparseMatrix<double> Adiag;
        igl::adjacency_matrix(Fuv, A);
        igl::sum(A,1,Asum);
        igl::diag(Asum,Adiag);
        A = Adiag - A;
	}

	if (type == '2') {
        auto facesNum = F.rows();
        auto vNum = V.rows();
        b = VectorXd::Zero(2*vNum);

        // duplicate vertices matrix and faces (with offset) to receive 2 times bigger matrix A (for u and v)
        Eigen::MatrixXi Fuv = Eigen::MatrixXi(facesNum*2, F.cols());
        Eigen::MatrixXd Vuv = Eigen::MatrixXd(V.rows()*2, V.cols());
        for (int i = 0; i < facesNum; i++) {
            Eigen::VectorXi face = F.row(i);
            Fuv.row(i) = face;
            Fuv.row(i+facesNum) = face.array() + vNum;
        }
        for (int i = 0; i < V.rows(); i++) {
            Vuv.row(i) = V.row(i);
            Vuv.row(i+V.rows()) = V.row(i);
        }

        SparseVector<double> Asum;
        SparseMatrix<double> Adiag;
        igl::cotmatrix(Vuv, Fuv, A);
        igl::sum(A,1,Asum);
        igl::diag(Asum,Adiag);
        A = Adiag - A;
	}

    SparseMatrix<double> Dx, Dy;
    computeSurfaceGradientMatrix(Dx, Dy);

    Eigen::VectorXd dblA;
    igl::doublearea(V, F, dblA);
    SparseMatrix<double> dblADiag = SparseMatrix<double>(dblA.asDiagonal());

	if (type == '3') {
        SparseMatrix<double> a, b1, b2, col1, col2;

        // a - topLeft and bottomRight; b1 - topRight; b2 - bottomLeft
        a = Dx.transpose() * dblADiag * Dx + Dy.transpose() * dblADiag * Dy;
        b1 = -1 * Dx.transpose() * dblADiag * Dy + Dy.transpose() * dblADiag * Dx;
        b2 = -1 * Dy.transpose() * dblADiag * Dx + Dx.transpose() * dblADiag * Dy;
        igl::cat(1, a, b2, col1);
        igl::cat(1, b1, a, col2);
        igl::cat(2, col1, col2, A);
        b = VectorXd::Zero(A.rows());
	}

	if (type == '4') {
        if (UV.size() == 0) {
            computeParameterization('3');
        }
        SparseMatrix<double> area = (dblADiag * 0.5).cwiseSqrt();
        VectorXd R = VectorXd(4*F.rows());
        for (int i = 0; i < F.rows(); i++) {
            auto fArea = area.coeff(i, i);
            Eigen::Matrix2d J, U, S, v, r;
            J << Dx.row(i)*UV.col(0), Dy.row(i)*UV.col(0),
                 Dx.row(i)*UV.col(1), Dy.row(i)*UV.col(1);
            SSVD2x2(J, U, S, v);
            r = U * v.transpose();

            // only determinant > 0 is suitable; change sign by multiplying one row by -1
            if (r.determinant() < 0) {
                r.row(0) *= -1;
            }

            R[i] = fArea*r.coeff(0, 0);
            R[i+F.rows()] = fArea*r.coeff(0, 1);
            R[i+2*F.rows()] = fArea*r.coeff(1, 0);
            R[i+3*F.rows()] = fArea*r.coeff(1, 1);
        }
        SparseMatrix<double> areaDx = area*Dx;
        SparseMatrix<double> areaDy = area*Dy;
        SparseMatrix<double> a, row1, row2, row3, row4, row12, row34, zero;
        zero = SparseMatrix<double>(Dx.rows(), Dx.cols());
        // construct A (Dx 0; Dy 0; 0 Dx; 0 Dy) multiplied by sqrt(doublearea / 2)
        igl::cat(2, areaDx, zero, row1);
        igl::cat(2, areaDy, zero, row2);
        igl::cat(2, zero, areaDx, row3);
        igl::cat(2, zero, areaDy, row4);
        igl::cat(1, row1, row2, row12);
        igl::cat(1, row3, row4, row34);
        igl::cat(1, row12, row34, a);

        A = a.transpose() * a;

        // construct b equal to A.transpose() * R
        b = a.transpose() * R;
	}

    SparseMatrix<double> col1, col2, Ares;
    VectorXd bres;
    igl::cat(1, b, d, bres);
    igl::cat(1, A, C, col1);
    igl::cat(1, SparseMatrix<double>(C.transpose()),
            Eigen::SparseMatrix<double>(C.rows(), C.rows()), col2);
    igl::cat(2, col1, col2, Ares);

    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
    solver.analyzePattern(Ares);
    solver.factorize(Ares);

    if (solver.info() != Success) {
        std::cout << solver.lastErrorMessage() << std::endl;
        return;
    }
    VectorXd x = solver.solve(bres);

    // The solver will output a vector
	UV.resize(V.rows(), 2);
	UV.col(0) = x.head(V.rows());
	UV.col(1) = x.segment(V.rows(), V.rows());
}

bool callback_key_pressed(Viewer &viewer, unsigned char key, int modifiers) {
    showingDistortion = false;
	switch (key) {
	case '1':
	case '2':
	case '3':
	case '4':
		computeParameterization(key);
		break;

    // key for the distortion visualization (removed previous comments)
	case '5':
        showingDistortion = true;
        calculateDistortion();
		break;
	case '+':
		TextureResolution /= 2;
		break;
	case '-':
		TextureResolution *= 2;
		break;
	case ' ': // space bar -  switches view between mesh and parameterization
    if(showingUV)
    {
      temp2D = viewer.core();
      viewer.core() = temp3D;
      showingUV = false;
    }
    else
    {
      if(UV.rows() > 0)
      {
        temp3D = viewer.core();
        viewer.core() = temp2D;
        showingUV = true;
      }
      else { std::cout << "ERROR ! No valid parameterization\n"; }
    }
    break;
	}
	Redraw();
	return true;
}

bool load_mesh(string filename)
{
  igl::read_triangle_mesh(filename,V,F);
  Redraw();
  viewer.core().align_camera_center(V);
  showingUV = false;

  return true;
}

bool callback_init(Viewer &viewer)
{
	temp3D = viewer.core();
	temp2D = viewer.core();
	temp2D.orthographic = true;

	return false;
}

int main(int argc,char *argv[]) {
  if(argc != 2) {
    cout << "Usage ex4_bin <mesh.off/obj>" << endl;
    load_mesh("../data/cathead.obj");
  }
  else
  {
    // Read points and normals
    load_mesh(argv[1]);
  }

	igl::opengl::glfw::imgui::ImGuiPlugin plugin;
  viewer.plugins.push_back(&plugin);
  igl::opengl::glfw::imgui::ImGuiMenu menu;
  plugin.widgets.push_back(&menu);

	menu.callback_draw_viewer_menu = [&]()
	{
		// Draw parent menu content
		menu.draw_viewer_menu();

		// Add new group
		if (ImGui::CollapsingHeader("Parmaterization", ImGuiTreeNodeFlags_DefaultOpen))
		{
			// Expose variable directly ...
			ImGui::Checkbox("Free boundary", &freeBoundary);
            ImGui::Checkbox("Show fixed vertices", &showingFixedVert);
            std::vector<std::string> distortions = {"Conformal", "Isometric", "Authalic"};
            ImGui::Combo("Distortion type", &distType, distortions);
		}
	};

  viewer.callback_key_pressed = callback_key_pressed;
  viewer.callback_mouse_move = callback_mouse_move;
  viewer.callback_init = callback_init;

  viewer.launch();
}

#include <iostream>
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/adjacency_list.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_corner_normals.h>
#include <igl/facet_components.h>
#include <igl/jet.h>
#include <igl/barycenter.h>
#include <igl/edge_topology.h>

using namespace std;
using Viewer = igl::opengl::glfw::Viewer;

// Vertex array, #V x3
Eigen::MatrixXd V;
// Face array, #F x3
Eigen::MatrixXi F;
// Per-face normal array, #F x3
Eigen::MatrixXd FN;
// Per-vertex normal array, #V x3
Eigen::MatrixXd VN;
// Per-corner normal array, (3#F) x3
Eigen::MatrixXd CN;
// Vectors of indices for adjacency relations
std::vector<std::vector<int> > VF, VFi, VV;
// Integer vector of component IDs per face, #F x1
Eigen::VectorXi cid;
// Per-face color array, #F x3
Eigen::MatrixXd component_colors_per_face;

void subdivide_sqrt3(const Eigen::MatrixXd &V,
					 const Eigen::MatrixXi &F,
					 Eigen::MatrixXd &Vout,
					 Eigen::MatrixXi &Fout) {
    std::vector<std::vector<int>> vv;
    igl::adjacency_list(F, vv);

    Eigen::MatrixXd centers;
    igl::barycenter(V, F, centers);

    Eigen::MatrixXi FNew = Eigen::MatrixXi(centers.rows() * 3, F.cols());
    for (int i = 0; i < F.rows(); i++) {
        Eigen::Vector3i verts = F.row(i);
        for (int j = 0; j < 2; j++) {
            FNew.row(i * 3 + j) << verts[j], verts[j + 1], V.rows() + i;
        }
        FNew.row(i * 3 + 2) << verts[2], verts[0], V.rows() + i;
    }

    Eigen::MatrixXd P = Eigen::MatrixXd(V.rows(), V.cols());
    for (int i = 0; i < V.rows(); i++) {
        int n = vv[i].size();
        double an = (4.0f - 2.0f * std::cos((2 * M_PI) / n)) / 9.0f;
        Eigen::Vector3d sumV = Eigen::Vector3d::Zero();
        for (int index: vv[i]) {
            sumV += V.row(index);
        }
        Eigen::Vector3d prevVal = V.row(i);
        P.row(i) = (1.0f - an) * prevVal + (an / n) * sumV;
    }
    Vout = Eigen::MatrixXd(V.rows() + centers.rows(), V.cols());
    Vout << P, centers;
    Fout = FNew;

    Eigen::MatrixXi ev, fe, ef;
    igl::edge_topology(Vout, Fout, ev, fe, ef);
    for (int i = 0; i < ef.rows(); i++) {
        bool flip = (ev.row(i).array() < V.rows()).all();
        if (!flip) {
            continue;
        }
        std::array<int, 2> faces = {ef.coeff(i, 0), ef.coeff(i, 1)};
        if (faces[0] < 0 || faces[1] < 0) {
            continue;
        }
        std::vector<int> notOnEdge;
        for (int f: faces) {
            auto verts = FNew.row(f);
            for (int v = 0; v < 3; v++) {
                int vert = verts[v];
                if ((ev.row(i).array() != vert).all()) {
                    notOnEdge.push_back(vert);
                }
            }
        }
        if (notOnEdge.size() != 2) {
            continue;
        }
        Fout.row(faces[0]) << notOnEdge[0], ev.coeff(i, 0), notOnEdge[1];
        Fout.row(faces[1]) << notOnEdge[1], ev.coeff(i, 1), notOnEdge[0];
    }
}

bool callback_key_down(Viewer& viewer, unsigned char key, int modifiers) {
    if (key == '1') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);

        int n = V.rows();
        igl::vertex_triangle_adjacency(n, F, VF, VFi);
        printf("Vertex-to-Face Relationship\n");
        for (int i = 0; i < n; i++) {
            if (i >= VF.size()) {
                printf("ERROR: vertices num should be equal to VF size\n");
                continue;
            }
            printf("V %d: ", i);
            for (int face: VF[i]) {
                printf("%d ", face);
            }
            printf("\n");
        }
    }

    if (key == '2') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);

        igl::adjacency_list(F, VV);
        printf("Vertex-to-Vertex Relationship\n");
        for (int i = 0; i < V.rows(); i++) {
            if (i >= VV.size()) {
                printf("ERROR: vertices num should be equal to VV size\n");
                continue;
            }
            printf("V %d: ", i);
            for (int vert: VV[i]) {
                printf("%d ", vert);
            }
            printf("\n");
        }
    }

    if (key == '3') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        FN.setZero(F.rows(),3);

        igl::per_face_normals(V, F, FN);
        viewer.data().set_normals(FN);
    }

    if (key == '4') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);

        igl::per_vertex_normals(V, F, VN);
        viewer.data().set_normals(VN);
    }

    if (key == '5') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);

        igl::per_corner_normals(V, F, 80, CN);
        viewer.data().set_normals(CN);
    }

    if (key == '6') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        component_colors_per_face.setZero(F.rows(),3);

        igl::facet_components(F, cid);
        int max_value = cid.maxCoeff();
        std::vector<int> facePerComp(max_value + 1, 0);
        for (int i = 0; i < cid.size(); i++) {
            int comp = cid.coeff(i);
            double val = comp / (double)max_value;
            Eigen::Vector3d rgb;
            igl::jet(val, rgb.data());
            component_colors_per_face.row(i) = rgb;
            facePerComp[comp]++;
        }
        printf("Number of components %d\n", max_value + 1);
        for (int fpc : facePerComp) {
            printf("%d/", fpc);
        }
        printf("\n");
        viewer.data().set_colors(component_colors_per_face);
    }

    if (key == '7') {
		Eigen::MatrixXd Vout;
		Eigen::MatrixXi Fout;
        // Fill the subdivide_sqrt3() function with your code for sqrt(3) subdivision.
		subdivide_sqrt3(V,F,Vout,Fout);
        // Set up the viewer to display the new mesh
        V = Vout; F = Fout;
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
    }

    return true;
}

bool load_mesh(Viewer& viewer,string filename, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
  igl::readOFF(filename,V,F);
  viewer.data().clear();
  viewer.data().set_mesh(V,F);
  viewer.data().compute_normals();
  viewer.core().align_camera_center(V, F);
  return true;
}

int main(int argc, char *argv[]) {
    // Show the mesh
    Viewer viewer;
    viewer.callback_key_down = callback_key_down;
    
    std::string filename;
    if (argc == 2) {
        filename = std::string(argv[1]); // Mesh provided as command line argument
    }
    else {
        filename = std::string("../data/bunny.off"); // Default mesh
    }
	
    load_mesh(viewer,filename,V,F);

    callback_key_down(viewer, '1', 0);

    // Attach a menu plugin
    igl::opengl::glfw::imgui::ImGuiPlugin plugin;
	viewer.plugins.push_back(&plugin);
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    plugin.widgets.push_back(&menu);
    
    viewer.launch();
}

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
					 Eigen::MatrixXi &Fout){
	
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

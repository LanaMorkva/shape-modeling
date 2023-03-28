//
// Created by Sveta Morkva on 27.03.2023.
//

#include "MLSHelper.h"
#include <cfloat>
#include <Eigen/Dense>
#include <igl/bounding_box_diagonal.h>


MLSHelper::MLSHelper(const Eigen::MatrixXd &GP, const Eigen::MatrixXd &CP, const Eigen::VectorXd &CV, int res, int degree,
                     double radius) {
    m_GP = GP; m_CP = CP; m_CV = CV; m_res = res;
    m_degree = degree; m_radius = radius;
}

void MLSHelper::initGrid() {
    double eps_m = 0.01 * igl::bounding_box_diagonal(m_CP);
    bb_min = m_CP.colwise().minCoeff();
    bb_max = m_CP.colwise().maxCoeff();

    bb_min -= eps_m;
    bb_max += eps_m;

    // Bounding box dimensions
    dim = bb_max - bb_min;

    // Make cell size the same as radius
    Eigen::RowVector3d cell_sizef = dim / (double)(m_res - 1);
    const int dx = std::ceil(cell_sizef[0]);
    const int dy = std::ceil(cell_sizef[1]);
    const int dz = std::ceil(cell_sizef[2]);
    cell_size << dx, dy, dz;

    uni_grid.resize(m_res*m_res*m_res);
    for (int index = 0; index < m_CP.rows(); index++) {
        Eigen::Array3d p = m_CP.row(index);
        Eigen::Array3d cell_dist = (p - bb_min) / cell_size;
        int i = std::floor(cell_dist[0]);
        int j = std::floor(cell_dist[1]);
        int k = std::floor(cell_dist[2]);
        int cell_ind = i + j * m_res + k * m_res * m_res;
        uni_grid[cell_ind].push_back(index);
    }
}

void MLSHelper::calcGridValues(bool fast) {
    if (fast) {
        initGrid();
    }
    int size = m_res * m_res * m_res;
    m_GV.resize(size);
    for (int index = 0; index < size; index++) {
        auto gp = m_GP.row(index);
        const auto &points = fast ? withinDist(gp) : withinDist_Slow(gp);
        if (points.empty()) {
            // assign large positive value
            m_GV[index] = DBL_MAX;
            continue;
        }
        Eigen::VectorXd weights = Eigen::VectorXd(points.size());
        Eigen::VectorXd constrVal = Eigen::VectorXd(points.size());
        Eigen::MatrixXd A = Eigen::MatrixXd(points.size(), polynomVecSize());
        for (int p_it = 0; p_it < points.size(); p_it++) {
            int p_ind = points[p_it];
            auto cp = m_CP.row(p_ind);
            auto r = (cp - gp).norm();
            auto r_div_h = r/m_radius;

            weights[p_it] = std::pow(1.0 - r_div_h, 4) * (4*r_div_h+1);
            constrVal[p_it] = m_CV[p_ind];
            A.row(p_it) = getPolyVector(cp) * weights[p_it];
        }
        Eigen::VectorXd b = constrVal.cwiseProduct(weights);
        Eigen::RowVectorXd c = A.colPivHouseholderQr().solve(b);

        double gv = getPolyVector(gp).dot(c);
        m_GV[index] = gv;
    }
}

std::vector<int> MLSHelper::withinDist_Slow(const Eigen::RowVector3d &q) {
    std::vector<int> res;
    for (int i = 0; i < m_CP.rows(); i++) {
        auto dist = (m_CP.row(i) - q).norm();
        if (dist < m_radius) {
            res.push_back(i);
        }
    }
    return res;
}

std::vector<int> MLSHelper::withinDist(const Eigen::Array3d &q) {
    std::vector<int> res;
    int circlesDefInDist = std::floor(m_radius / cell_size.maxCoeff()) - 1;
    int maxCircleNum = std::ceil(m_radius / cell_size.minCoeff());

    Eigen::Array3d dist = (q - bb_min) / cell_size;
    Eigen::Array3d prev_cell = dist.floor();

    int i = prev_cell[0];
    int j = prev_cell[1];
    int k = prev_cell[2];

    for (int x = i-circlesDefInDist; x < i+circlesDefInDist+1; x++) {
        for (int y = j-circlesDefInDist; y < j+circlesDefInDist+1; y++) {
            for (int z = k-circlesDefInDist; z < k+circlesDefInDist+1; z++) {
                int cell = x + y * m_res + z * m_res * m_res;
                if (cell >= uni_grid.size()) {
                    continue;
                }
                for (int num: uni_grid[cell]) {
                    res.push_back(num);
                }
            }
        }
    }
    auto check_cell = [=, &res](int index) {
        if (index >= uni_grid.size()) {
            return;
        }
        for (int num: uni_grid[index]) {
            Eigen::Array3d p = m_CP.row(num);
            if ((p - q).matrix().norm() < m_radius) {
                res.push_back(num);
            }
        }
    };
    int leftCircles = maxCircleNum - circlesDefInDist;
    while (leftCircles > 0) {
        int circleNum = maxCircleNum - leftCircles + 1;
        int xL = i-circleNum;
        int xR = i+circleNum;
        for (int z = k-circleNum; z < k+circleNum+1; z++) {
            for (int y = j - circleNum; y < j + circleNum + 1; y++) {
                int base = y * m_res + z * m_res*m_res;
                int cellL = xL + base;
                int cellR = xR + base;
                check_cell(cellL);
                check_cell(cellR);
            }
        }
        int yD = j-circleNum;
        int yU = j+circleNum;
        for (int z = k-circleNum; z < k+circleNum+1; z++) {
            for (int x = xL+1; x < xR; x++) {
                int base = z * m_res*m_res;
                int cellD = x + yD*m_res + base;
                int cellU = x + yU*m_res + base;
                check_cell(cellD);
                check_cell(cellU);
            }
        }
        int zD = k-circleNum;
        int zU = k+circleNum;
        for (int y = yD+1; y < yU; y++) {
            for (int x = xL + 1; x < xR; x++) {
                int base = x + y*m_res;
                int cellD = base + zD*m_res*m_res;
                int cellU = base + zU*m_res*m_res;
                check_cell(cellD);
                check_cell(cellU);
            }
        }
        leftCircles--;
    }
    return res;
}

int MLSHelper::polynomVecSize() const {
    switch (m_degree) {
        case 1: return 4;
        case 2: return 10;
    }
    return 1;
}

Eigen::VectorXd MLSHelper::getPolyVector(const Eigen::RowVector3d &q) const {
    Eigen::VectorXd res = Eigen::VectorXd(polynomVecSize());
    double x = q[0], y = q[1], z = q[2];
    switch (m_degree) {
        case 1:
            res << 1, x, y, z;
            return res;
        case 2:
            res << 1, x, y, z, x*y, x*z, y*z, x*x, y*y, z*z;
            return res;
    }
    res << 1;
    return res;
}
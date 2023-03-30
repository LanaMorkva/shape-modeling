//
// Created by Sveta Morkva on 27.03.2023.
//

#include "MLSHelper.h"
#include <cfloat>
#include <iostream>
#include <Eigen/Dense>
#include <igl/bounding_box_diagonal.h>


MLSHelper::MLSHelper(const Eigen::MatrixXd &GP, const Eigen::MatrixXd &CP, const Eigen::VectorXd &CV, int res, int degree,
                     double radius) {
    m_GP = GP; m_CP = CP; m_CV = CV; m_res = res;
    m_degree = degree; m_radius = radius;
}

void MLSHelper::initGrid() {
    double eps_m = 0.05 * igl::bounding_box_diagonal(m_CP);
    bb_min = m_CP.colwise().minCoeff();
    bb_max = m_CP.colwise().maxCoeff();

    bb_min -= eps_m;
    bb_max += eps_m;

    // Bounding box dimensions
    dim = bb_max - bb_min;

    // Make cell size the same as radius
    cell_size << eps_m, eps_m, eps_m;
    Eigen::RowVector3d cell_num = dim / eps_m;
    m_gridx = std::ceil(cell_num[0]);
    m_gridy = std::ceil(cell_num[1]);
    m_gridz = std::ceil(cell_num[2]);
    m_gridxdy = m_gridx*m_gridy;

    uni_grid.resize(m_gridxdy*m_gridz);
    for (int index = 0; index < m_CP.rows(); index++) {
        Eigen::Array3d p = m_CP.row(index);
        Eigen::Array3d cell_dist = (p - bb_min) / cell_size;
        int i = std::floor(cell_dist[0]);
        int j = std::floor(cell_dist[1]);
        int k = std::floor(cell_dist[2]);
        int cell_ind = i + j * m_gridx + k * m_gridxdy;
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

void MLSHelper::calcGridValues_NormConstr(const Eigen::MatrixXd &N) {
    int size = m_res * m_res * m_res;
    m_GV.resize(size);
    for (int index = 0; index < size; index++) {
        auto gp = m_GP.row(index);
        const auto &points = withinDist_Slow(gp);
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
            Eigen::VectorXd pDiff = gp-cp;
            Eigen::VectorXd n = N.row(p_ind).normalized();
            auto r = pDiff.norm();
            auto r_div_h = r/m_radius;

            weights[p_it] = std::pow(1.0 - r_div_h, 4) * (4*r_div_h+1);
            constrVal[p_it] = pDiff.dot(n);
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

std::vector<int> MLSHelper::withinDist(const Eigen::Vector3d &q) {
    std::vector<int> res;
    Eigen::Array3d dist = (q.array() - bb_min) / cell_size;
    Eigen::Array3d prev_cell = dist.floor();
    Eigen::Array3d next_cell = prev_cell+1;

    // find minimum dist to cell's sides
    Eigen::Array3d prevCellDist = (dist - prev_cell);
    Eigen::Array3d nextCellDist = (next_cell - dist);
    Eigen::Array3d lNumCell = (m_radius/cell_size - prevCellDist).ceil();
    Eigen::Array3d rNumCell = (m_radius/cell_size - nextCellDist).ceil();

    int i = prev_cell[0];
    int j = prev_cell[1];
    int k = prev_cell[2];

    int minx = std::min(std::max(i-(int)lNumCell[0], 0), m_gridx-1);
    int miny = std::min(std::max(j-(int)lNumCell[1], 0), m_gridy-1);
    int minz = std::min(std::max(k-(int)lNumCell[2], 0), m_gridz-1);

    int maxx = std::min(std::max(i+(int)rNumCell[0], 0), m_gridx-1);
    int maxy = std::min(std::max(j+(int)rNumCell[1], 0), m_gridy-1);
    int maxz = std::min(std::max(k+(int)rNumCell[2], 0), m_gridz-1);

    int cell = minx + miny * m_gridx + minz * m_gridxdy;
    for (int z = minz; z < maxz+1; z++) {
        for (int y = miny; y < maxy+1; y++) {
            for (int x = minx; x < maxx+1; x++) {
                int distX = (x<i) ? x : x+1;
                int distY = (y<j) ? y : y+1;
                int distZ = (z<k) ? z : z+1;
                Eigen::Vector3d distCorn;
                distCorn << distX * cell_size[0], distY * cell_size[1], distZ * cell_size[2];
                auto maxDistInCell = (distCorn - q).norm();
                // full cell within distance
                if (maxDistInCell < m_radius) {
                    res.insert(res.end(), uni_grid[cell].begin(), uni_grid[cell].end());
                } else {
                    for (int num: uni_grid[cell]) {
                        Eigen::Vector3d p = m_CP.row(num);
                        auto norm = (p - q).norm();
                        if (norm < m_radius) {
                            res.push_back(num);
                        }
                    }
                }
                cell++;
            }
            cell += m_gridx - (maxx + 1 - minx);
        }
        cell += m_gridxdy - m_gridx * (maxy + 1 - miny);
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
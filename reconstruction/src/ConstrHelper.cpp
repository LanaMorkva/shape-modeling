//
// Created by Sveta Morkva on 27.03.2023.
//

#include "ConstrHelper.h"
#include <cfloat>
#include <igl/bounding_box_diagonal.h>

ConstrHelper::ConstrHelper(const Eigen::MatrixXd &P, const Eigen::MatrixXd &N, bool fast, int res) {
    m_P = P; m_N = N; m_fast = fast; m_res = res;
    if (m_fast) {
        initGrid();
    }
}

void ConstrHelper::initGrid() {
    double eps_m = 0.01 * igl::bounding_box_diagonal(m_P);
    bb_min = m_P.colwise().minCoeff();
    bb_max = m_P.colwise().maxCoeff();

    bb_min -= eps_m;
    bb_max += eps_m;

    // Bounding box dimensions
    dim = bb_max - bb_min;

    // Grid spacing
//    const int size = std::ceil(eps_m);
    cell_size << eps_m, eps_m, eps_m;
    Eigen::RowVector3d cell_num = dim / eps_m;
    m_gridx = std::ceil(cell_num[0]);
    m_gridy = std::ceil(cell_num[1]);
    m_gridz = std::ceil(cell_num[2]);

    uni_grid.resize(m_gridx*m_gridy*m_gridz);
    for (int index = 0; index < m_P.rows(); index++) {
        Eigen::Array3d p = m_P.row(index);
        Eigen::Array3d cell_dist = (p - bb_min) / cell_size;
        int i = std::floor(cell_dist[0]);
        int j = std::floor(cell_dist[1]);
        int k = std::floor(cell_dist[2]);
        int cell_ind = i + j * m_gridx + k * m_gridx * m_gridy;
        uni_grid[cell_ind].push_back(index);
    }
}

void ConstrHelper::calcConstraints() {
    m_CP = Eigen::MatrixXd(m_P.rows() * 3, m_P.cols());
    m_CV = Eigen::VectorXd(m_P.rows() * 3);

    // add points to constraints
    m_CP.block(0,0, m_P.rows(), 3) = m_P;
    m_CV.topRows(m_P.rows()) = Eigen::VectorXd::Zero(m_P.rows());

    double eps_m = 0.01 * igl::bounding_box_diagonal(m_P);
    // add + epsilon
    m_row_num = m_P.rows();
    for (int i = 0; i < m_P.rows(); i++) {
        addEpsConstr(i, eps_m, 1);
        m_row_num++;
    }
    // add - epsilon
    for (int i = 0; i < m_P.rows(); i++) {
        addEpsConstr(i, eps_m, -1);
        m_row_num++;
    }
}

void ConstrHelper::addEpsConstr(int p_ind, double eps, int sign) {
    Eigen::Vector3d p = m_P.row(p_ind);
    Eigen::Vector3d n = m_N.row(p_ind).normalized();
    while (true) {
        Eigen::Vector3d pEps = p + sign * eps * n;
        int minDistInd = m_fast ? closestToPoint(pEps) : closestToPoint_Slow(pEps);

        if (minDistInd == p_ind) {
            m_CP.row(m_row_num) = pEps;
            m_CV[m_row_num] = sign*eps;

            // add new point to grid
            if (m_fast) {
                Eigen::Array3d dist = (pEps.array() - bb_min) / cell_size;
                Eigen::Array3d prev_cell = dist.floor();
                int i = prev_cell[0];
                int j = prev_cell[1];
                int k = prev_cell[2];
                int q_cell = i + j * m_gridx + k * m_gridx * m_gridy;
                uni_grid[q_cell].push_back(m_row_num);
            }
            return;
        };

        eps /= 2.0;
    }
}

int ConstrHelper::closestToPoint_Slow(const Eigen::RowVector3d &q) {
    int minDistInd = -1;
    double minDist = DBL_MAX;
    for (int i = 0; i < m_row_num; i++) {
        auto dist = (m_CP.row(i) - q).norm();
        minDistInd = dist < minDist ? i : minDistInd;
        minDist = dist < minDist ? dist : minDist;
    }
    return minDistInd;
}

int ConstrHelper::closestToPoint(const Eigen::Array3d &q) {
    int minDistInd = -1;
    double minDist = DBL_MAX;

    // Grid spacing
    Eigen::Array3d dist = (q - bb_min) / cell_size;
    Eigen::Array3d prev_cell = dist.floor();
    Eigen::Array3d next_cell = prev_cell+1;

    // find minimum dist to cell's sides
    Eigen::Array3d prevCellDist = (dist - prev_cell) * cell_size;
    Eigen::Array3d nextCellDist = (next_cell - dist) * cell_size;

    int i = prev_cell[0];
    int j = prev_cell[1];
    int k = prev_cell[2];
    int q_cell = i + j * m_gridx + k * m_gridx * m_gridy;

    auto check_cell = [=, &minDistInd, &minDist](int index) {
        if (index >= uni_grid.size()) {
            return;
        }
        for (int num: uni_grid[index]) {
            Eigen::Array3d p = m_CP.row(num);
            auto norm = (p - q).matrix().norm();
            if (norm < minDist) {
                minDistInd = num;
                minDist = norm;
            }
        }
    };

    // check the cell with q
    check_cell(q_cell);
    for (int x = i-1; x < i+2; x++) {
        for (int y = j-1; y < j+2; y++) {
            for (int z = k-1; z < k+2; z++) {
                if (x==i && y==j && z==k) {
                    continue;
                }
                if ((x<i && prevCellDist[0]>minDist) || (x>i && nextCellDist[0]>minDist) ||
                (y<j && prevCellDist[1]>minDist) || (y>j && nextCellDist[1]>minDist) ||
                (z<k && prevCellDist[2]>minDist) || (z>k && nextCellDist[2]>minDist)) {
                    continue;
                }
                int cell = x + y * m_gridx + z * m_gridx * m_gridy;
                check_cell(cell);
            }
        }
    }
    return minDistInd;
}
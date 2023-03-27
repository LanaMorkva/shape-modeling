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
    bb_min = m_P.colwise().minCoeff();
    bb_max = m_P.colwise().maxCoeff();

    // Bounding box dimensions
    dim = bb_max - bb_min;

    // Grid spacing
    Eigen::RowVector3d cell_sizef = dim / (double)(m_res - 1);
    const int dx = std::ceil(cell_sizef[0]);
    const int dy = std::ceil(cell_sizef[1]);
    const int dz = std::ceil(cell_sizef[2]);
    cell_size << dx, dy, dz;

    uni_grid.resize(m_res*m_res*m_res);
    for (int index = 0; index < m_P.rows(); index++) {
        Eigen::Array3d p = m_P.row(index);
        Eigen::Array3d cell_dist = (p - bb_min) / cell_size;
        int i = std::floor(cell_dist[0]);
        int j = std::floor(cell_dist[1]);
        int k = std::floor(cell_dist[2]);
        int cell_ind = i + j * m_res + k * m_res * m_res;
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
                int q_cell = i + j * m_res + k * m_res * m_res;
                uni_grid[q_cell].push_back(m_row_num);
            }
            return;
        }
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
    double prevCellDist = ((dist - prev_cell) * cell_size).minCoeff();
    double nextCellDist = ((next_cell - dist) * cell_size).minCoeff();
    double borderDist = std::min(prevCellDist, nextCellDist);

    int i = prev_cell[0];
    int j = prev_cell[1];
    int k = prev_cell[2];
    int q_cell = i + j * m_res + k * m_res * m_res;

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
    // return if dist to sides bigger than found point
    if (minDist < borderDist) {
        return minDistInd;
    }

    // go in circles around the closest cells
    double minCellDist = cell_size.minCoeff();
    int circleNum = 1;
    borderDist += minCellDist;
    while (true) {
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
        if (minDist < borderDist) {
            return minDistInd;
        }
        borderDist += minCellDist;
        circleNum++;
    }
}
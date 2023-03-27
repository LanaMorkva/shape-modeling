//
// Created by Sveta Morkva on 27.03.2023.
//

#include "MLSHelper.h"
#include <cfloat>
#include <Eigen/Dense>


MLSHelper::MLSHelper(const Eigen::MatrixXd &GP, const Eigen::MatrixXd &CP, const Eigen::VectorXd &CV, int res, int degree,
                     double radius) {
    m_GP = GP; m_CP = CP; m_CV = CV; m_res = res;
    m_degree = degree; m_radius = radius;
}

void MLSHelper::calcGridValues() {
    m_GV.resize(m_res * m_res * m_res);
    for (int x = 0; x < m_res; x++) {
        for (int y = 0; y < m_res; y++) {
            for (int z = 0; z < m_res; z++) {
                int index = x + m_res * (y + m_res * z);
                auto gp = m_GP.row(index);
                const auto &points = withinDist_Slow(gp);
                // TODO: add fast method
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
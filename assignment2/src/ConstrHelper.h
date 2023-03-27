//
// Created by Sveta Morkva on 27.03.2023.
//

#ifndef ASSIGNMENT2_CONSTRHELPER_H
#define ASSIGNMENT2_CONSTRHELPER_H

#include <Eigen/Core>
#include <vector>

class ConstrHelper {
public:
    ConstrHelper(const Eigen::MatrixXd &P, const Eigen::MatrixXd &N, bool fast = false, int res = 0);
    void calcConstraints();
    const Eigen::MatrixXd &constrPoints() const { return m_CP; }
    const Eigen::VectorXd &constrValues() const { return m_CV; }

private:
    void initGrid();
    void addEpsConstr(int p_ind, double eps_m, int sign);
    int closestToPoint_Slow(const Eigen::RowVector3d &q);
    int closestToPoint(const Eigen::Array3d &q);

    Eigen::MatrixXd m_P, m_N, m_CP;
    Eigen::VectorXd m_CV;
    bool m_fast;
    int m_res, m_row_num;
    std::vector<std::vector<int>> uni_grid;
    Eigen::Array3d bb_min, bb_max, dim;
    Eigen::Array3d cell_size;
};


#endif //ASSIGNMENT2_CONSTRHELPER_H

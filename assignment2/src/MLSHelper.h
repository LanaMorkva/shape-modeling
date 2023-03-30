//
// Created by Sveta Morkva on 27.03.2023.
//

#ifndef ASSIGNMENT2_MLSHELPER_H
#define ASSIGNMENT2_MLSHELPER_H

#include <Eigen/Core>
#include <vector>

class MLSHelper {
public:
    MLSHelper(const Eigen::MatrixXd &GP, const Eigen::MatrixXd &CP, const Eigen::VectorXd &CV, int res, int degree,
              double radius);

    void calcGridValues(bool fast = false);

    // need to pass regular points, without added constraints to the constructor; CV are not needed
    void calcGridValues_NormConstr(const Eigen::MatrixXd &N);
    const Eigen::VectorXd &getGridValues() const { return m_GV; }

private:
    void initGrid();
    std::vector<int> withinDist_Slow(const Eigen::RowVector3d &q);
    std::vector<int> withinDist(const Eigen::Vector3d &q);
    int polynomVecSize() const;
    Eigen::VectorXd getPolyVector(const Eigen::RowVector3d &q) const;

    Eigen::MatrixXd m_GP, m_CP;
    Eigen::VectorXd m_GV, m_CV;
    int m_res, m_degree;
    double m_radius;

    std::vector<std::vector<int>> uni_grid;
    Eigen::Array3d bb_min, bb_max, dim;
    Eigen::Array3d cell_size;
    int m_gridx, m_gridy, m_gridz, m_gridxdy;
};


#endif //ASSIGNMENT2_MLSHELPER_H

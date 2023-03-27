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

    void calcGridValues();
    const Eigen::VectorXd &getGridValues() const { return m_GV; }

private:
    std::vector<int> withinDist_Slow(const Eigen::RowVector3d &q);
    int polynomVecSize() const;
    Eigen::VectorXd getPolyVector(const Eigen::RowVector3d &q) const;

    Eigen::MatrixXd m_GP, m_CP;
    Eigen::VectorXd m_GV, m_CV;
    int m_res, m_degree;
    double m_radius;
};


#endif //ASSIGNMENT2_MLSHELPER_H

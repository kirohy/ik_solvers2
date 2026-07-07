#ifndef IK_CONSTRAINT2_ORIENTATIONCONSTRAINT_H
#define IK_CONSTRAINT2_ORIENTATIONCONSTRAINT_H

#include <ik_constraint2/IKConstraint.h>
#include <cnoid/EigenUtil>
#include <cnoid/LinkPath>
#include <iostream>

namespace ik_constraint2{
  class OrientationConstraint : public IKConstraint
  {
  public:
    // A_link中のA_localRとB_link中のB_localRを一致させる.
    // リンクがnullptrならworld座標系を意味する.
    // maxError, weight, precisionはeval系の姿勢3軸に対する値.
    const cnoid::LinkPtr& A_link() const { return A_link_;}
    cnoid::LinkPtr& A_link() { return A_link_;}
    const cnoid::Matrix3d& A_localR() const { return A_localR_;}
    cnoid::Matrix3d& A_localR() { return A_localR_;}
    const cnoid::LinkPtr& B_link() const { return B_link_;}
    cnoid::LinkPtr& B_link() { return B_link_;}
    const cnoid::Matrix3d& B_localR() const { return B_localR_;}
    cnoid::Matrix3d& B_localR() { return B_localR_;}
    const cnoid::Vector3& maxError() const { return maxError_;}
    cnoid::Vector3& maxError() { return maxError_;}
    const double& precision() const { return precision_;}
    double& precision() { return precision_;}
    const cnoid::Vector3& weight() const { return weight_;}
    cnoid::Vector3& weight() { return weight_;}
    const cnoid::LinkPtr& eval_link() const { return eval_link_;}
    cnoid::LinkPtr& eval_link() { return eval_link_;}
    const cnoid::Matrix3d& eval_localR() const { return eval_localR_;}
    cnoid::Matrix3d& eval_localR() { return eval_localR_;}

    virtual void updateBounds () override;
    virtual void updateJacobian (const std::vector<cnoid::LinkPtr>& joints) override;
    virtual bool isSatisfied () const override;
    virtual double distance() const override;
    virtual std::vector<cnoid::SgNodePtr>& getDrawOnObjects() override;
    virtual std::shared_ptr<IKConstraint> clone(const std::map<cnoid::BodyPtr, cnoid::BodyPtr>& modelMap) const override;
    void copy(std::shared_ptr<OrientationConstraint> ret, const std::map<cnoid::BodyPtr, cnoid::BodyPtr>& modelMap) const;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  protected:
    cnoid::LinkPtr A_link_ = nullptr;
    cnoid::Matrix3d A_localR_ = cnoid::Matrix3d::Identity();
    cnoid::LinkPtr B_link_ = nullptr;
    cnoid::Matrix3d B_localR_ = cnoid::Matrix3d::Identity();
    cnoid::Vector3 maxError_ = cnoid::Vector3::Constant(0.05);
    double precision_ = 1e-3;
    cnoid::Vector3 weight_ = cnoid::Vector3::Ones();
    cnoid::LinkPtr eval_link_ = nullptr;
    cnoid::Matrix3d eval_localR_ = cnoid::Matrix3d::Identity();

    cnoid::SgLineSetPtr lines_;
    cnoid::Vector3 current_error_eval_ = cnoid::Vector3::Zero();

    std::vector<cnoid::LinkPtr> path_A_joints_;
    std::vector<cnoid::LinkPtr> path_B_joints_;
    std::vector<cnoid::LinkPtr> path_eval_joints_;
    Eigen::SparseMatrix<double,Eigen::RowMajor> jacobian_A_full_;
    Eigen::SparseMatrix<double,Eigen::RowMajor> jacobian_B_full_;
    Eigen::SparseMatrix<double,Eigen::RowMajor> jacobian_eval_full_;
    Eigen::SparseMatrix<double,Eigen::RowMajor> jacobian_full_local_;
    cnoid::LinkPtr jacobian_A_link_ = nullptr;
    cnoid::LinkPtr jacobian_B_link_ = nullptr;
    cnoid::LinkPtr jacobian_eval_link_ = nullptr;

    std::vector<cnoid::LinkPtr> jacobian_joints_;
    std::unordered_map<cnoid::LinkPtr,int> jacobianColMap_;
  };
}

#endif

#include <ik_constraint2/VelocityConstraint.h>
#include <ik_constraint2/Jacobian.h>

namespace ik_constraint2{
  void VelocityConstraint::updateBounds () {
    const cnoid::Vector3& A_v = (this->A_link_) ? this->A_link_->v() + this->A_link_->w().cross(A_localpos_.translation()) : this->A_link_->body()->rootLink()->v();
    const cnoid::Vector3& A_w = (this->A_link_) ? this->A_link_->w() : this->A_link_->body()->rootLink()->w();
    const cnoid::Vector3& B_v = (this->B_link_) ? this->B_link_->v() + this->B_link_->w().cross(B_localpos_.translation()) : cnoid::Vector3(base_velocity_.head<3>());
    const cnoid::Vector3& B_w = (this->B_link_) ? this->B_link_->w() : cnoid::Vector3(base_velocity_.tail<3>());

    cnoid::Matrix3d eval_R = (this->eval_link_) ? this->eval_link_->R() * this->eval_localR_ : this->eval_localR_;

    cnoid::Vector6 error; // world frame. A-B
    const cnoid::Vector3 v_error = A_v - B_v - this->target_velocity_.head<3>();
    cnoid::Vector3 w_error = A_w - B_w - this->target_velocity_.tail<3>();

    cnoid::Vector6 error_eval; // eval frame. A-B
    error_eval.head<3>() = eval_R.transpose() * v_error * this->dt_;
    error_eval.tail<3>() = eval_R.transpose() * w_error * this->dt_;

    // A-Bの目標変位を計算し、this->eq_に入れる
    if(this->eq_.rows()!=(this->weight_.array() > 0.0).count()) this->eq_ = Eigen::VectorXd((this->weight_.array() > 0.0).count());
    int idx=0;
    for(size_t i=0; i<6; i++){
      if(this->weight_[i]>0.0) {
        this->eq_[idx] = std::min(std::max(-error_eval[i],-this->maxError_[i]),this->maxError_[i]) * this->weight_[i];
        idx++;
      }
    }

    // distanceを計算する
    this->current_error_eval_ = error_eval;

    if(this->debugLevel_>=1){
      std::cerr << "VelocityConstraint " << (this->A_link_?this->A_link_->name():std::string("world")) << " : " << (this->B_link_?this->B_link_->name():std::string("world")) << std::endl;
      std::cerr << "A_v" << std::endl;
      std::cerr << A_v << std::endl;
      std::cerr << "A_w" << std::endl;
      std::cerr << A_w << std::endl;
      std::cerr << "B_v" << std::endl;
      std::cerr << B_v << std::endl;
      std::cerr << "B_w" << std::endl;
      std::cerr << B_w << std::endl;
      std::cerr << "error_eval" << std::endl;
      std::cerr << error_eval.transpose() << std::endl;
      std::cerr << "eq" << std::endl;
      std::cerr << this->eq_.transpose() << std::endl;
    }
  }

  void VelocityConstraint::updateJacobian (const std::vector<cnoid::LinkPtr>& joints) {
    // this->jacobian_を計算する
    // 行列の初期化. 前回とcol形状が変わっていないなら再利用
    if(!IKConstraint::isJointsSame(joints,this->jacobian_joints_)
       || this->A_link_ != this->jacobian_A_link_
       || this->B_link_ != this->jacobian_B_link_
       || this->eval_link_ != this->jacobian_eval_link_){
      this->jacobian_joints_ = joints;
      this->jacobian_A_link_ = this->A_link_;
      this->jacobian_B_link_ = this->B_link_;
      this->jacobian_eval_link_ = this->eval_link_;

      ik_constraint2::calc6DofJacobianShape(this->jacobian_joints_,//input
                                            this->jacobian_A_link_,//input
                                            this->jacobian_A_full_,
                                            this->jacobianColMap_,
                                            this->path_A_joints_
                                            );
      ik_constraint2::calc6DofJacobianShape(this->jacobian_joints_,//input
                                            this->jacobian_B_link_,//input
                                            this->jacobian_B_full_,
                                            this->jacobianColMap_,
                                            this->path_B_joints_
                                            );
      ik_constraint2::calc6DofJacobianShape(this->jacobian_joints_,//input
                                            this->jacobian_eval_link_,//input
                                            this->jacobian_eval_full_,
                                            this->jacobianColMap_,
                                            this->path_eval_joints_
                                            );
    }

    ik_constraint2::calc6DofJacobianCoef(this->jacobian_joints_,//input
                                         this->jacobian_A_link_,//input
                                         this->A_localpos_.translation(),//input
                                         this->jacobianColMap_,//input
                                         this->path_A_joints_,//input
                                         this->jacobian_A_full_
                                         );
    ik_constraint2::calc6DofJacobianCoef(this->jacobian_joints_,//input
                                         this->jacobian_B_link_,//input
                                         this->B_localpos_.translation(),//input
                                         this->jacobianColMap_,//input
                                         this->path_B_joints_,//input
                                         this->jacobian_B_full_
                                         );
    ik_constraint2::calc6DofJacobianCoef(this->jacobian_joints_,//input
                                         this->jacobian_eval_link_,//input
                                         cnoid::Vector3::Zero(),//input
                                         this->jacobianColMap_,//input
                                         this->path_eval_joints_,//input
                                         this->jacobian_eval_full_
                                         );

    cnoid::Matrix3d eval_R = (this->eval_link_) ? this->eval_link_->R() * this->eval_localR_ : this->eval_localR_;
    Eigen::SparseMatrix<double,Eigen::RowMajor> eval_R_sparse(3,3);
    for(int i=0;i<3;i++) for(int j=0;j<3;j++) eval_R_sparse.insert(i,j) = eval_R(i,j);
    this->jacobian_full_local_.resize(6, this->jacobian_A_full_.cols());

    this->jacobian_full_local_.topRows<3>() = eval_R_sparse.transpose() * this->jacobian_A_full_.topRows<3>();
    this->jacobian_full_local_.bottomRows<3>() = eval_R_sparse.transpose() * this->jacobian_A_full_.bottomRows<3>();
    this->jacobian_full_local_.topRows<3>() -= Eigen::SparseMatrix<double,Eigen::RowMajor>(eval_R_sparse.transpose() * this->jacobian_B_full_.topRows<3>());
    this->jacobian_full_local_.bottomRows<3>() -= Eigen::SparseMatrix<double,Eigen::RowMajor>(eval_R_sparse.transpose() * this->jacobian_B_full_.bottomRows<3>());
    this->jacobian_full_local_.topRows<3>() += IKConstraint::cross(this->current_error_eval_.head<3>()) * eval_R_sparse.transpose() * this->jacobian_eval_full_.bottomRows<3>();
    this->jacobian_full_local_.bottomRows<3>() += IKConstraint::cross(this->current_error_eval_.tail<3>()) * eval_R_sparse.transpose() * this->jacobian_eval_full_.bottomRows<3>();

    this->jacobian_.resize((this->weight_.array() > 0.0).count(),this->jacobian_full_local_.cols());
    for(size_t i=0, idx=0;i<6;i++){
      if(this->weight_[i]>0.0) {
        this->jacobian_.row(idx) = this->weight_[i] * this->jacobian_full_local_.row(i);
        idx++;
      }
    }

    // this->jacobianIneq_のサイズだけそろえる
    this->jacobianIneq_.resize(0,this->jacobian_.cols());

    if(this->debugLevel_>=1){
      std::cerr << "VelocityConstraint" << std::endl;
      std::cerr << "jacobian" << std::endl;
      std::cerr << this->jacobian_ << std::endl;
    }

    return;
  }

  bool VelocityConstraint::isSatisfied () const {
    return this->distance() <= this->precision_;
  }

  double VelocityConstraint::distance () const {
    return this->current_error_eval_.cwiseProduct(this->weight_).norm();
  }


  std::vector<cnoid::SgNodePtr>& VelocityConstraint::getDrawOnObjects(){
    if(!this->lines_){
      this->lines_ = new cnoid::SgLineSet;
      this->lines_->setLineWidth(1.0);
      this->lines_->getOrCreateColors()->resize(4);
      this->lines_->getOrCreateColors()->at(0) = cnoid::Vector3f(1.0,1.0,1.0);
      this->lines_->getOrCreateColors()->at(1) = cnoid::Vector3f(1.0,0.0,0.0);
      this->lines_->getOrCreateColors()->at(2) = cnoid::Vector3f(0.0,1.0,0.0);
      this->lines_->getOrCreateColors()->at(3) = cnoid::Vector3f(0.0,0.0,1.0);
      // A, A_x, A_y, A_z, B, B_x, B_y, B_z
      this->lines_->getOrCreateVertices()->resize(8);
      this->lines_->colorIndices().resize(0);
      this->lines_->addLine(0,1); this->lines_->colorIndices().push_back(1); this->lines_->colorIndices().push_back(1);
      this->lines_->addLine(0,2); this->lines_->colorIndices().push_back(2); this->lines_->colorIndices().push_back(2);
      this->lines_->addLine(0,3); this->lines_->colorIndices().push_back(3); this->lines_->colorIndices().push_back(3);
      this->lines_->addLine(4,5); this->lines_->colorIndices().push_back(1); this->lines_->colorIndices().push_back(1);
      this->lines_->addLine(4,6); this->lines_->colorIndices().push_back(2); this->lines_->colorIndices().push_back(2);
      this->lines_->addLine(4,7); this->lines_->colorIndices().push_back(3); this->lines_->colorIndices().push_back(3);
      this->lines_->addLine(0,4); this->lines_->colorIndices().push_back(0); this->lines_->colorIndices().push_back(0);

      this->drawOnObjects_ = std::vector<cnoid::SgNodePtr>{this->lines_};
    }

    const cnoid::Isometry3& A_pos = (this->A_link_) ? this->A_link_->T() * this->A_localpos_ : this->A_localpos_;
    const cnoid::Isometry3& B_pos = (this->B_link_) ? this->B_link_->T() * this->B_localpos_ : this->B_localpos_;

    this->lines_->getOrCreateVertices()->at(0) = A_pos.translation().cast<cnoid::Vector3f::Scalar>();
    this->lines_->getOrCreateVertices()->at(1) = (A_pos * (0.05 * cnoid::Vector3::UnitX())).cast<cnoid::Vector3f::Scalar>();
    this->lines_->getOrCreateVertices()->at(2) = (A_pos * (0.05 * cnoid::Vector3::UnitY())).cast<cnoid::Vector3f::Scalar>();
    this->lines_->getOrCreateVertices()->at(3) = (A_pos * (0.05 * cnoid::Vector3::UnitZ())).cast<cnoid::Vector3f::Scalar>();
    this->lines_->getOrCreateVertices()->at(4) = B_pos.translation().cast<cnoid::Vector3f::Scalar>();
    this->lines_->getOrCreateVertices()->at(5) = (B_pos * (0.05 * cnoid::Vector3::UnitX())).cast<cnoid::Vector3f::Scalar>();
    this->lines_->getOrCreateVertices()->at(6) = (B_pos * (0.05 * cnoid::Vector3::UnitY())).cast<cnoid::Vector3f::Scalar>();
    this->lines_->getOrCreateVertices()->at(7) = (B_pos * (0.05 * cnoid::Vector3::UnitZ())).cast<cnoid::Vector3f::Scalar>();

    return this->drawOnObjects_;
  }

  std::shared_ptr<IKConstraint> VelocityConstraint::clone(const std::map<cnoid::BodyPtr, cnoid::BodyPtr>& modelMap) const {
    std::shared_ptr<VelocityConstraint> ret = std::make_shared<VelocityConstraint>(*this);
    this->copy(ret, modelMap);
    return ret;
  }

  void VelocityConstraint::copy(std::shared_ptr<VelocityConstraint> ret, const std::map<cnoid::BodyPtr, cnoid::BodyPtr>& modelMap) const {
    if(this->A_link_ && modelMap.find(this->A_link_->body()) != modelMap.end()) ret->A_link() = modelMap.find(this->A_link_->body())->second->link(this->A_link_->index());
    if(this->B_link_ && modelMap.find(this->B_link_->body()) != modelMap.end()) ret->B_link() = modelMap.find(this->B_link_->body())->second->link(this->B_link_->index());
    if(this->eval_link_ && modelMap.find(this->eval_link_->body()) != modelMap.end()) ret->eval_link() = modelMap.find(this->eval_link_->body())->second->link(this->eval_link_->index());
  }
}

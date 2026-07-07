#include <ik_constraint2/OrientationConstraint.h>
#include <algorithm>
#include <cmath>

namespace ik_constraint2{
  namespace {
    inline Eigen::Matrix3d orientCoordToAxis(const Eigen::Matrix3d& m, const Eigen::Vector3d& axis, const Eigen::Vector3d& localaxis){
      // axisとlocalaxisはノルムが1, mは回転行列でなければならない.
      // axisとlocalaxisがピッタリ180反対向きの場合、回転方向が定まらないので不安定.
      Eigen::AngleAxisd m_ = Eigen::AngleAxisd(m);
      Eigen::Vector3d localaxisdir = m_ * localaxis;
      Eigen::Vector3d cross = localaxisdir.cross(axis);
      double dot = std::min(1.0,std::max(-1.0,localaxisdir.dot(axis)));
      if(cross.norm()==0){
        if(dot == -1) return Eigen::Matrix3d(-m);
        else return Eigen::Matrix3d(m_);
      }else{
        double angle = std::acos(dot);
        Eigen::Vector3d axis = cross.normalized();
        return Eigen::Matrix3d(Eigen::AngleAxisd(angle, axis) * m_);
      }
    }

    void pushBackAngularTripletList(std::vector<Eigen::Triplet<double> >& tripletList, const cnoid::LinkPtr& joint, int idx){
      if(joint->isFreeJoint()){
        for(size_t d=0;d<3;d++) tripletList.push_back(Eigen::Triplet<double>(d,idx+3+d,1));
      } else if(joint->isRevoluteJoint()){
        for(size_t d=0;d<3;d++) tripletList.push_back(Eigen::Triplet<double>(d,idx,1));
      }
    }

    void fillAngularJacobian(Eigen::SparseMatrix<double,Eigen::RowMajor>& jacobian, const cnoid::LinkPtr& joint, int idx, int sign){
      if(joint->isFreeJoint()){
        for(size_t d=0;d<3;d++) jacobian.coeffRef(d,idx+3+d) = sign;
      } else if(joint->isRevoluteJoint()){
        cnoid::Vector3 omega = joint->R() * joint->a();
        for(size_t d=0;d<3;d++) jacobian.coeffRef(d,idx) = sign * omega[d];
      }
    }

    void calcAngularJacobianShape(const std::vector<cnoid::LinkPtr>& joints,
                                  const cnoid::LinkPtr& A_link,
                                  Eigen::SparseMatrix<double,Eigen::RowMajor>& jacobian,
                                  std::unordered_map<cnoid::LinkPtr,int>& jacobianColMap,
                                  std::vector<cnoid::LinkPtr>& path_A_joints){
      jacobianColMap.clear();
      int num_variables = 0;
      for(size_t i=0;i<joints.size();i++){
        jacobianColMap[joints[i]] = num_variables;
        num_variables += IKConstraint::getJointDOF(joints[i]);
      }

      std::vector<Eigen::Triplet<double> > tripletList;
      if(A_link){
        tripletList.reserve(100);
        path_A_joints.clear();
        cnoid::LinkPath path(A_link);
        for(size_t j=0;j<path.size();j++) path_A_joints.push_back(path[j]);
        for(size_t j=0;j<path_A_joints.size();j++){
          cnoid::LinkPtr joint = path_A_joints[j];
          if(jacobianColMap.find(joint)==jacobianColMap.end()) continue;
          pushBackAngularTripletList(tripletList,joint,jacobianColMap[joint]);
        }
      }

      jacobian = Eigen::SparseMatrix<double,Eigen::RowMajor>(3,num_variables);
      jacobian.setFromTriplets(tripletList.begin(), tripletList.end());
    }

    void calcAngularJacobianCoef(const cnoid::LinkPtr& A_link,
                                 std::unordered_map<cnoid::LinkPtr,int>& jacobianColMap,
                                 const std::vector<cnoid::LinkPtr>& path_A_joints,
                                 Eigen::SparseMatrix<double,Eigen::RowMajor>& jacobian,
                                 int sign){
      if(!A_link) return;
      for(size_t j=0;j<path_A_joints.size();j++){
        cnoid::LinkPtr joint = path_A_joints[j];
        if(jacobianColMap.find(joint)==jacobianColMap.end()) continue;
        fillAngularJacobian(jacobian,joint,jacobianColMap[joint],sign);
      }
    }
  }

  void OrientationConstraint::updateBounds () {
    const cnoid::Matrix3d A_R = (this->A_link_) ? this->A_link_->R() * this->A_localR_ : this->A_localR_;
    const cnoid::Matrix3d B_R = (this->B_link_) ? this->B_link_->R() * this->B_localR_ : this->B_localR_;
    cnoid::Matrix3d eval_R = (this->eval_link_) ? this->eval_link_->R() * this->eval_localR_ : this->eval_localR_;

    cnoid::Vector3 rot_error = cnoid::Vector3::Zero();
    if((this->weight_.array() > 0.0).count() == 2 &&
       ((this->eval_link_ == this->A_link_ && this->A_localR_ == this->eval_localR_) || (this->eval_link_ == this->B_link_ && this->B_localR_ == this->eval_localR_)) ) {
      cnoid::Vector3 axis;
      if(this->weight_[0] == 0.0) axis = cnoid::Vector3::UnitX();
      else if(this->weight_[1] == 0.0) axis = cnoid::Vector3::UnitY();
      else if(this->weight_[2] == 0.0) axis = cnoid::Vector3::UnitZ();

      cnoid::Vector3 A_axis;
      cnoid::Vector3 B_axis;
      if(this->eval_link_ == this->A_link_){
        A_axis = eval_R * axis;
        B_axis = B_R * A_R.transpose() * A_axis;
      }else{
        B_axis = eval_R * axis;
        A_axis = A_R * B_R.transpose() * B_axis;
      }
      Eigen::Vector3d cross = B_axis.cross(A_axis);
      double dot = std::min(1.0,std::max(-1.0,B_axis.dot(A_axis)));
      if(cross.norm()==0){
        if(dot == -1){
          if(this->weight_[0] == 0.0) rot_error = eval_R * M_PI * cnoid::Vector3::UnitY();
          else if(this->weight_[1] == 0.0) rot_error = eval_R * M_PI * cnoid::Vector3::UnitZ();
          else if(this->weight_[2] == 0.0) rot_error = eval_R * M_PI * cnoid::Vector3::UnitX();
        }
      }else{
        double angle = std::acos(dot);
        Eigen::Vector3d axis_ = cross.normalized();
        rot_error = angle * axis_;
      }
    }else if((this->weight_.array() > 0.0).count() == 1 &&
             ((this->eval_link_ == this->A_link_ && this->A_localR_ == this->eval_localR_) || (this->eval_link_ == this->B_link_ && this->B_localR_ == this->eval_localR_)) ) {
      cnoid::Vector3 axis;
      if(this->weight_[0] > 0.0) axis = cnoid::Vector3::UnitX();
      else if(this->weight_[1] > 0.0) axis = cnoid::Vector3::UnitY();
      else if(this->weight_[2] > 0.0) axis = cnoid::Vector3::UnitZ();

      cnoid::Matrix3 A_R_aligned;
      cnoid::Matrix3 B_R_aligned;
      if(this->eval_link_ == this->A_link_ && this->A_localR_ == this->eval_localR_){
        A_R_aligned = A_R;
        B_R_aligned = orientCoordToAxis(B_R, A_R * axis, axis);
      }else{
        A_R_aligned = orientCoordToAxis(A_R, B_R * axis, axis);
        B_R_aligned = B_R;
      }
      const cnoid::AngleAxis angleAxis = cnoid::AngleAxis(A_R_aligned * B_R_aligned.transpose());
      rot_error = angleAxis.angle()*angleAxis.axis();
    }else{
      const cnoid::AngleAxis angleAxis = cnoid::AngleAxis(A_R * B_R.transpose());
      rot_error = angleAxis.angle()*angleAxis.axis();
    }

    cnoid::Vector3 error_eval = eval_R.transpose() * rot_error;
    if(this->eq_.rows()!=(this->weight_.array() > 0.0).count()) this->eq_ = Eigen::VectorXd((this->weight_.array() > 0.0).count());
    int idx=0;
    for(size_t i=0; i<3; i++){
      if(this->weight_[i]>0.0) {
        this->eq_[idx] = std::min(std::max(-error_eval[i],-this->maxError_[i]),this->maxError_[i]) * this->weight_[i];
        idx++;
      }
    }

    this->current_error_eval_ = error_eval;

    if(this->debugLevel_>=1){
      std::cerr << "OrientationConstraint " << (this->A_link_?this->A_link_->name():std::string("world")) << " : " << (this->B_link_?this->B_link_->name():std::string("world")) << std::endl;
      std::cerr << "error_eval" << std::endl;
      std::cerr << error_eval.transpose() << std::endl;
      std::cerr << "eq" << std::endl;
      std::cerr << this->eq_.transpose() << std::endl;
    }
  }

  void OrientationConstraint::updateJacobian (const std::vector<cnoid::LinkPtr>& joints) {
    if(!IKConstraint::isJointsSame(joints,this->jacobian_joints_)
       || this->A_link_ != this->jacobian_A_link_
       || this->B_link_ != this->jacobian_B_link_
       || this->eval_link_ != this->jacobian_eval_link_){
      this->jacobian_joints_ = joints;
      this->jacobian_A_link_ = this->A_link_;
      this->jacobian_B_link_ = this->B_link_;
      this->jacobian_eval_link_ = this->eval_link_;

      ik_constraint2::calcAngularJacobianShape(this->jacobian_joints_,
                                               this->jacobian_A_link_,
                                               this->jacobian_A_full_,
                                               this->jacobianColMap_,
                                               this->path_A_joints_);
      std::unordered_map<cnoid::LinkPtr,int> unusedColMap;
      ik_constraint2::calcAngularJacobianShape(this->jacobian_joints_,
                                               this->jacobian_B_link_,
                                               this->jacobian_B_full_,
                                               unusedColMap,
                                               this->path_B_joints_);
      ik_constraint2::calcAngularJacobianShape(this->jacobian_joints_,
                                               this->jacobian_eval_link_,
                                               this->jacobian_eval_full_,
                                               unusedColMap,
                                               this->path_eval_joints_);
    }

    for(int k=0;k<this->jacobian_A_full_.outerSize();++k) for(Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator it(this->jacobian_A_full_,k); it; ++it) it.valueRef() = 0.0;
    for(int k=0;k<this->jacobian_B_full_.outerSize();++k) for(Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator it(this->jacobian_B_full_,k); it; ++it) it.valueRef() = 0.0;
    for(int k=0;k<this->jacobian_eval_full_.outerSize();++k) for(Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator it(this->jacobian_eval_full_,k); it; ++it) it.valueRef() = 0.0;

    ik_constraint2::calcAngularJacobianCoef(this->jacobian_A_link_,
                                            this->jacobianColMap_,
                                            this->path_A_joints_,
                                            this->jacobian_A_full_,
                                            1);
    ik_constraint2::calcAngularJacobianCoef(this->jacobian_B_link_,
                                            this->jacobianColMap_,
                                            this->path_B_joints_,
                                            this->jacobian_B_full_,
                                            1);
    ik_constraint2::calcAngularJacobianCoef(this->jacobian_eval_link_,
                                            this->jacobianColMap_,
                                            this->path_eval_joints_,
                                            this->jacobian_eval_full_,
                                            1);

    cnoid::Matrix3d eval_R = (this->eval_link_) ? this->eval_link_->R() * this->eval_localR_ : this->eval_localR_;
    Eigen::SparseMatrix<double,Eigen::RowMajor> eval_R_sparse(3,3);
    for(int i=0;i<3;i++) for(int j=0;j<3;j++) eval_R_sparse.insert(i,j) = eval_R(i,j);
    this->jacobian_full_local_.resize(3, this->jacobian_A_full_.cols());

    this->jacobian_full_local_ = eval_R_sparse.transpose() * this->jacobian_A_full_;
    this->jacobian_full_local_ -= Eigen::SparseMatrix<double,Eigen::RowMajor>(eval_R_sparse.transpose() * this->jacobian_B_full_);
    this->jacobian_full_local_ += IKConstraint::cross(this->current_error_eval_) * eval_R_sparse.transpose() * this->jacobian_eval_full_;

    this->jacobian_.resize((this->weight_.array() > 0.0).count(),this->jacobian_full_local_.cols());
    for(size_t i=0, idx=0;i<3;i++){
      if(this->weight_[i]>0.0) {
        this->jacobian_.row(idx) = this->weight_[i] * this->jacobian_full_local_.row(i);
        idx++;
      }
    }

    this->jacobianIneq_.resize(0,this->jacobian_.cols());

    if(this->debugLevel_>=1){
      std::cerr << "OrientationConstraint" << std::endl;
      std::cerr << "jacobian" << std::endl;
      std::cerr << this->jacobian_ << std::endl;
    }
  }

  bool OrientationConstraint::isSatisfied () const {
    return this->distance() <= this->precision_;
  }

  double OrientationConstraint::distance () const {
    return this->current_error_eval_.cwiseProduct(this->weight_).norm();
  }

  std::vector<cnoid::SgNodePtr>& OrientationConstraint::getDrawOnObjects(){
    if(!this->lines_){
      this->lines_ = new cnoid::SgLineSet;
      this->lines_->setLineWidth(1.0);
      this->lines_->getOrCreateColors()->resize(4);
      this->lines_->getOrCreateColors()->at(0) = cnoid::Vector3f(1.0,1.0,1.0);
      this->lines_->getOrCreateColors()->at(1) = cnoid::Vector3f(1.0,0.0,0.0);
      this->lines_->getOrCreateColors()->at(2) = cnoid::Vector3f(0.0,1.0,0.0);
      this->lines_->getOrCreateColors()->at(3) = cnoid::Vector3f(0.0,0.0,1.0);
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

    const cnoid::Matrix3d A_R = (this->A_link_) ? this->A_link_->R() * this->A_localR_ : this->A_localR_;
    const cnoid::Matrix3d B_R = (this->B_link_) ? this->B_link_->R() * this->B_localR_ : this->B_localR_;
    cnoid::Vector3 A_p = cnoid::Vector3::Zero();
    cnoid::Vector3 B_p = cnoid::Vector3::Zero();
    if(this->A_link_) A_p = this->A_link_->p();
    if(this->B_link_) B_p = this->B_link_->p();

    this->lines_->getOrCreateVertices()->at(0) = A_p.cast<cnoid::Vector3f::Scalar>();
    this->lines_->getOrCreateVertices()->at(1) = (A_p + A_R * (0.05 * cnoid::Vector3::UnitX())).cast<cnoid::Vector3f::Scalar>();
    this->lines_->getOrCreateVertices()->at(2) = (A_p + A_R * (0.05 * cnoid::Vector3::UnitY())).cast<cnoid::Vector3f::Scalar>();
    this->lines_->getOrCreateVertices()->at(3) = (A_p + A_R * (0.05 * cnoid::Vector3::UnitZ())).cast<cnoid::Vector3f::Scalar>();
    this->lines_->getOrCreateVertices()->at(4) = B_p.cast<cnoid::Vector3f::Scalar>();
    this->lines_->getOrCreateVertices()->at(5) = (B_p + B_R * (0.05 * cnoid::Vector3::UnitX())).cast<cnoid::Vector3f::Scalar>();
    this->lines_->getOrCreateVertices()->at(6) = (B_p + B_R * (0.05 * cnoid::Vector3::UnitY())).cast<cnoid::Vector3f::Scalar>();
    this->lines_->getOrCreateVertices()->at(7) = (B_p + B_R * (0.05 * cnoid::Vector3::UnitZ())).cast<cnoid::Vector3f::Scalar>();

    return this->drawOnObjects_;
  }

  std::shared_ptr<IKConstraint> OrientationConstraint::clone(const std::map<cnoid::BodyPtr, cnoid::BodyPtr>& modelMap) const {
    std::shared_ptr<OrientationConstraint> ret = std::make_shared<OrientationConstraint>(*this);
    this->copy(ret, modelMap);
    return ret;
  }

  void OrientationConstraint::copy(std::shared_ptr<OrientationConstraint> ret, const std::map<cnoid::BodyPtr, cnoid::BodyPtr>& modelMap) const {
    ret->A_link() = applyModelMap(this->A_link_, modelMap);
    ret->B_link() = applyModelMap(this->B_link_, modelMap);
    ret->eval_link() = applyModelMap(this->eval_link_, modelMap);
    ret->lines_ = nullptr;
  }
}

#pragma once

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <memory>

// constraint class
struct Dirichlet
{
  int flag[3];   ///< 1=constrainted, 0=not constrainted
  double val[3]; ///< constrained value
};

/// node class
struct Node
{
  double x[3];                       ///< coordinate
  std::shared_ptr<Dirichlet> dirich; ///< Dirichlet constraint
};

// element class
class Element
{
public:
  void stiffnessMatrix(int ndim, int voigt, const Eigen::MatrixXd &D,
                       std::vector<Eigen::Triplet<double>> &triplets,
                       std::vector<Node> &nodes);

  double misesStress(int ndim, int voigt, const Eigen::MatrixXd &D,
                     std::vector<Node> &nodes, Eigen::VectorXd &displacement);

  const int ne = 3;    /// nodes in this element
  const int ipmax = 1; ///< number of gauss point

  int numdof;              ///< dof in this element
  std::vector<int> nodeID; ///< node ID in this element
};

// Create a Bmatrix of triangular elements.
void bmatrixTri3(Eigen::MatrixXd &B, Eigen::MatrixXd &X, double &jac,
                 double &weight, int ip);

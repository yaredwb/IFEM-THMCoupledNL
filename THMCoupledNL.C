// $Id$
//==============================================================================
//!
//! \file THMCoupledNL.C
//!
//! \date
//!
//! \author Yared Bekele
//!
//! \brief Integrand implementations for THM coupled problems in ground freezing
//!
//==============================================================================

#include "THMCoupledNL.h"
#include "ASMbase.h"
#include "ASMmxBase.h"
#include "FiniteElement.h"
#include "TimeDomain.h"
#include "Utilities.h"
#include "Tensor.h"
#include "ElmMats.h"
#include "ElmNorm.h"
#include "Vec3Oper.h"
#include "VTF.h"
#include "StabilizationUtils.h"

//! \brief Enum for element level solution vectors
enum SolutionVectors
{
  Uo = 0,             //!< Previous displacement
  Po = 1,             //!< Previous pore pressure
  To = 2,             //!< Previous temperature
  Uc = 3,             //!< Current displacement
  Pc = 4,             //!< Current pore pressure
  Tc = 5,             //!< Current temperature
  NSOL = 6            //!< Number of solution vectors
};


//! \brief Enum for element level right-hand-side vectors
enum ResidualVectors
{
  Ru = 0,             //!< RHS vector for displacement
  Rp = 1,             //!< RHS vector for pore pressure
  RT = 2,             //!< RHS vector for temperature
  Rprev = 3,          //!< Internal vector from previous step
  Rnow = 4,           //!< Internal vector for current step
  Rres = 5,           //!< Residual RHS vector
  Rc = 6,             //!< Weak Dirichlet boundary contribution to RHS
  NVEC = 7            //!< Number of RHS vectors
};


//! \brief Enum for element level left-hand-side matrices
enum TangentMatrices
{
  uu = 0,             //!< Stiffness matrix
  up = 1,             //!< Mechanical-hydraulic coupling matrix
  uT = 2,             //!< Mechanical-thermal coupling matrix
  pu = 3,             //!< Hydro-mechanical coupling matrix
  pp = 4,             //!< Hydraulic matrix
  pT = 5,             //!< Hydro-thermal coupling matrix
  Tp = 6,             //!< Thermo-hydraulic coupling matrix
  TT = 7,             //!< Thermal matrix
  Kprev = 8,          //!< Fully coupled matrix from previous step
  Know = 9,           //!< Fully coupled matrix at current step
  Ktan = 10,          //!< Fully coupled Jacobian matrix
  NMAT = 11           //!< Number of LHS element matrices
};


THMCoupledNL::MixedElmMats::MixedElmMats()
{
  this->resize(NMAT,NVEC);
}


const Matrix& THMCoupledNL::MixedElmMats::getNewtonMatrix() const
{
  Matrix& N = const_cast<Matrix&>(A[Ktan]);

  size_t ru = A[uu].rows();
  size_t rp = A[pp].rows();

  for (size_t i = 1; i <= ru; i++)
  {
    for (size_t j = 1; j <= ru; j++)
    {
      N(i,j) = A[uu](i,j);
    }
    for (size_t j = 1; j <= rp; j++)
    {
      size_t k = ru+2*j-1;
      N(i,k) = A[up](i,j);
      N(k,i) = A[pu](j,i);
      size_t l = ru+2*j;
      N(i,l) = A[uT](i,j);
    }
  }

  for (size_t i = 1; i <= rp; i++)
  {
    for (size_t j = 1; j <= rp; j++)
    {
      size_t ki = ru+2*i-1;
      size_t kj = ru+2*j-1;
      size_t li = ru+2*i;
      size_t lj = ru+2*j;
      N(ki,kj) = A[pp](i,j);
      N(ki,lj) = A[pT](i,j);
      N(li,lj) = A[TT](i,j);
      N(li,kj) = A[Tp](i,j);
    }
  }

  return A[Ktan];
}


const Vector& THMCoupledNL::MixedElmMats::getRHSVector() const
{
  Vector& R = const_cast<Vector&>(b[Rres]);

  size_t ru = b[Ru].size();
  size_t rp = b[Rp].size();

  for (size_t i = 1; i <= ru; i++)
    R(i) = b[Ru](i);

  for (size_t i = 1; i <= rp; i++)
  {
    R(ru+2*i-1) = b[Rp](i);
    R(ru+2*i  ) = b[RT](i);
  }

  R += b[Rprev];
  R -= b[Rnow];

  // Robin boundary contribution to the residual
  // TO DO: This should be done properly by calling finalizeElementBou after
  //        integration of the Robin boundary terms.
  R -= A[Know]*b[Rc];

  return b[Rres];
}


THMCoupledNL::THMCoupledNL(unsigned short int n, int order, bool stab) :
  nsd(n), gacc(9.81), mat(nullptr), SUPG(stab)
{
  primsol.resize(1+order);
  tracFld = nullptr;
  fluxFld = nullptr;
}


Vec3 THMCoupledNL::getTraction(const Vec3& X, const Vec3& n) const
{
  if (fluxFld)
    return (*fluxFld)(X);
  else if (tracFld)
    return (*tracFld)(X,n);
  else
    return Vec3();
}


LocalIntegral* THMCoupledNL::getLocalIntegral(const std::vector<size_t>& nen,
                                              size_t, bool neumann) const
{
  const size_t nedof1 = nsd*nen[0];           //!< Number of DOFs on basis 1
  const size_t nedof  = nedof1 + 2*nen[1];    //!< Total number of DOFs

  ElmMats* result = new MixedElmMats();

  result->rhsOnly = neumann;
  result->withLHS = !neumann;
  result->b[Ru].resize(nedof1);
  result->b[Rp].resize(nen[1]);
  result->b[RT].resize(nen[1]);
  result->b[Rprev].resize(nedof);
  result->b[Rnow].resize(nedof);
  result->b[Rres].resize(nedof);
  result->b[Rc].resize(nedof);

  if(!neumann)
  {
    result->A[uu].resize(nedof1,nedof1);
    result->A[up].resize(nedof1,nen[1]);
    result->A[uT].resize(nedof1,nen[1]);
    result->A[pu].resize(nen[1],nedof1);
    result->A[pp].resize(nen[1],nen[1]);
    result->A[pT].resize(nen[1],nen[1]);
    result->A[Tp].resize(nen[1],nen[1]);
    result->A[TT].resize(nen[1],nen[1]);
    result->A[Kprev].resize(nedof,nedof);
    result->A[Know].resize(nedof,nedof);
    result->A[Ktan].resize(nedof,nedof);
  }

  return result;
}


bool THMCoupledNL::initElement(const std::vector<int>& MNPC,
                               const std::vector<size_t>& elem_sizes,
                               const std::vector<size_t>& basis_sizes,
                               LocalIntegral& elmInt)
{
  if(primsol.front().empty()) return true;

  // Extract element level solution vectors
  elmInt.vec.resize(NSOL);
  std::vector<int>::const_iterator fstart = MNPC.begin() + elem_sizes[0];
  int ierr = utl::gather(IntVec(MNPC.begin(), fstart),nsd,primsol[0],elmInt.vec[Uc])
           + utl::gather(IntVec(fstart,MNPC.end()),0,2,primsol[0],elmInt.vec[Pc],nsd*basis_sizes[0],basis_sizes[0])
           + utl::gather(IntVec(fstart,MNPC.end()),1,2,primsol[0],elmInt.vec[Tc],nsd*basis_sizes[0],basis_sizes[0])
           + utl::gather(IntVec(MNPC.begin(), fstart),nsd,primsol[1],elmInt.vec[Uo])
           + utl::gather(IntVec(fstart,MNPC.end()),0,2,primsol[1],elmInt.vec[Po],nsd*basis_sizes[0],basis_sizes[0])
           + utl::gather(IntVec(fstart,MNPC.end()),1,2,primsol[1],elmInt.vec[To],nsd*basis_sizes[0],basis_sizes[0]);

  if (ierr == 0) return true;

  std::cerr << " *** THMCoupledNL::initElement: Detected " << ierr/3
            << " node numbers out of range." << std::endl;

  return false;
}


bool THMCoupledNL::initElementBou(const std::vector<int>& MNPC,
                                  const std::vector<size_t>& elem_sizes,
                                  const std::vector<size_t>& basis_sizes,
                                  LocalIntegral& elmInt)
{
  return this->IntegrandBase::initElementBou(MNPC,elem_sizes,basis_sizes,elmInt);
}


bool THMCoupledNL::evalIntMx(LocalIntegral& elmInt, const MxFiniteElement& fe,
                             const TimeDomain& time, const Vec3& X) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  if (!mat)
  {
    std::cerr << __FUNCTION__ << ": No material data." << std::endl;
    return false;
  }

  size_t i,j,k;

  Matrix Bmat, Cmat, CB;

  // Get the strain-displacement matrix
  if(!mat->formBmatrix(Bmat,fe.grad(1),nsd))
    return false;

  // Get the updated pressure and temperature values at the current point
  double p = elMat.vec[Pc].dot(fe.basis(2));
  double T = elMat.vec[Tc].dot(fe.basis(2));
  // Rescale the pressure and temperature values
  p *= scl1;
  T *= scl2;

  // Evaluate the updated material tangent stiffness
  if(!mat->formElasticMatrix(Cmat,X,nsd,p,T))
    return false;

  // Integration of the stiffness matrix
  CB.multiply(Cmat,Bmat,false,false);
  CB *= 1.0*fe.detJxW;
  elMat.A[uu].multiply(Bmat,CB,true,false,true);

  // Get the densities of water and ice and the porosity
  double rhow = mat->getWaterDensity(X);
  double rhoi = mat->getIceDensity(X);
  double poro = mat->getPorosity(X);
  // Evaluate the ice pressure at the current point
  double pi = mat->getIcePressure(X,p,T);
  // Evaluate the degrees of saturation (Sw,Si) and water capacities (Sp,ST)
  double Sw = mat->getWaterSaturation(X,p,T);
  double Si = mat->getIceSaturation(X,p,T);
  double Sp = mat->getWaterCapacity(X,p,T,true);
  double ST = mat->getWaterCapacity(X,p,T,false);
  // Evaluate the relative permeability coefficient and the unfrozen permeability
  double kr = mat->getRelPermCoeff(X,p,T);
  Vec3 perm = mat->getPermeability(X);

#if SP_DEBUG > 3
  std::cout << "Degree of water saturation,     Sw = " << Sw << std::endl;
  std::cout << "Degree of ice saturation,       Si = " << Si << std::endl;
  std::cout << "Isothermal water capacity,      Sp = " << Sp << std::endl;
  std::cout << "Nonisothermal water capacity,   ST = " << ST << std::endl;
  std::cout << "Relative permeability coeff.,   kr = " << kr << std::endl;
#endif

  // Define the unit Voigt vector
  Vector m, Cm;
  m.resize(Cmat.rows());
  for (i = 1; i <= Cmat.rows(); i++)
    if (i <= nsd)
      m(i) = 1;

  // Evaluate Cmat*m
  Cm.resize(Cmat.rows());
  Cm = Cmat*m;

  // Integration of the mechanical-hydraulic matrix
  Matrix Cuptmp;
  const size_t nstrc = nsd*(nsd+1)/2;
  Cuptmp.resize(nstrc,fe.basis(2).size());
  // Evaluate the derivatives of total stress and phase change strain wrt pressure
  double dsdp = Sp * (pi - p) - (Sw + Si * (rhoi/rhow));
  double dedp = poro * (rhoi - rhow) * Sp / (rhow*Sw + rhoi*Si) / 3.0;

#if SP_DEBUG > 3
  std::cout << "Change of stress with pore pressure,       dsdp = " << dsdp << std::endl;
  std::cout << "Change of phase strain with pore pressure, dedp = " << dedp << std::endl;
#endif

  for (i = 1; i <= nstrc; i++)
    for (j = 1; j <= fe.basis(2).size(); j++)
      Cuptmp(i,j) += scl1*(dsdp*m(i) - dedp*Cm(i))*fe.basis(2)(j)*fe.detJxW;

  elMat.A[up].multiply(Bmat,Cuptmp,true,false,true);

  // Integration of the thermo-mechanical matrix
  // Get the latent heat of fusion
  double Lf = mat->getLatentHeat(X);
  // Evaluate the derivatives of total stress and phase change strain wrt temperature
  double dsdT = ST * (pi - p) + Si * rhoi * Lf/T;
  double dedT = poro * (rhoi - rhow) * ST / (rhow*Sw + rhoi*Si) / 3.0;

#if SP_DEBUG > 3
  std::cout << "Change of stress with temperature,       dsdT = " << dsdT << std::endl;
  std::cout << "Change of phase strain with temperature, dedT = " << dedT << std::endl;
#endif

  Matrix CuTtmp;
  CuTtmp.resize(nstrc,fe.basis(2).size());
  for (i = 1; i <= nstrc; i++)
    for (j = 1; j <= fe.basis(2).size(); j++)
      CuTtmp(i,j) += scl2*(dsdT*m(i) - dedT*Cm(i))*fe.basis(2)(j)*fe.detJxW;

  elMat.A[uT].multiply(Bmat,CuTtmp,true,false,true);

  // Integration of the hydraulic matrix (coefficient to constant)
  Matrix Kpp;
  Kpp.resize(fe.basis(2).size(),fe.basis(2).size());
  for (i = 1; i <= fe.basis(2).size(); i++)
    for (j= 1; j <= fe.basis(2).size(); j++)
      for (k = 1; k <= nsd; k++)
        Kpp(i,j) += scl1*scl1*fe.grad(2)(i,k)*(kr/(rhow*gacc))*perm[k-1]*fe.grad(2)(j,k)*fe.detJxW;

  // Integration of the hydro-mechanical matrix
  Matrix Cputmp;
  Cputmp.resize(fe.basis(2).size(),nstrc);
  double theta = Si/Sw;
  for (i = 1; i <= fe.basis(2).size(); i++)
    for (j = 1; j <= nstrc; j++)
      Cputmp(i,j) += scl1*((rhow*Sw + rhoi*Si)/(rhow + theta*rhoi))*fe.basis(2)(i)*m(j)*fe.detJxW;

  elMat.A[pu].multiply(Cputmp,Bmat,false,false,true);

  // Integration of the hydraulic matrix
  Matrix Cpp;
  Cpp.resize(fe.basis(2).size(),fe.basis(2).size());
  double tmp = poro*(rhow - rhoi) / (rhow + theta*rhoi);
  for (i = 1; i <= fe.basis(2).size(); i++)
    for (j = 1; j <= fe.basis(2).size(); j++)
      Cpp(i,j) += scl1*scl1*fe.basis(2)(i)*tmp*Sp*fe.basis(2)(j)*fe.detJxW;

  elMat.A[pp] += Cpp;
  elMat.A[pp].add(Kpp,time.dt);

  // Integration of the hydro-thermal matrix
  Matrix CpT;
  CpT.resize(fe.basis(2).size(),fe.basis(2).size());
  // Introduce a second scaling parameter for temperature coupling
  for (i = 1; i <= fe.basis(2).size(); i++)
    for (j = 1; j <= fe.basis(2).size(); j++)
      CpT(i,j) += scl1*scl2*fe.basis(2)(i)*tmp*ST*fe.basis(2)(j)*fe.detJxW;

  elMat.A[pT] += CpT;

  // Integration of the thermo-hdraulic matrix
  double xi = (poro * rhoi) / (Sw + (rhoi/rhow) * Si);
  Matrix CTp;
  CTp.resize(fe.basis(2).size(),fe.basis(2).size());
  for (i = 1; i <= fe.basis(2).size(); i++)
    for (j = 1; j <= fe.basis(2).size(); j++)
      CTp(i,j) += scl1*scl2*fe.basis(2)(i)*Lf*xi*Sp*fe.basis(2)(j)*fe.detJxW;

  elMat.A[Tp] += CTp;

  // Integration of the thermal matrix (coefficient to constant)
  Matrix KTT;
  KTT.resize(fe.basis(2).size(),fe.basis(2).size());
  double cw = mat->getWaterHeatCapacity(T);
  double ci = mat->getIceHeatCapacity(T);
  double rhoc_adv = rhow * cw + theta * rhoi * ci;
  Vector gradP;
  fe.grad(2).multiply(elMat.vec[Pc],gradP,true);
  Vec3 grav = this->getGravity();
  Vec3 vel;
  for (k = 1; k <= nsd; k++)
    vel[k-1] = -1.0*scl1*(kr / (rhow*gacc)) * perm[k-1] * (gradP(k) - rhow*grav[k-1]);

#if SP_DEBUG > 3
  std::cout << "Darcy velocity at current point, vel = " << vel << std::endl;
#endif

  // Evaluate the overall thermal conductivity
  double lambda = mat->getThermalConductivity(X,p,T);
  for (i = 1; i <= fe.basis(2).size(); i++)
  {
    for (j = 1; j <= fe.basis(2).size(); j++)
    {
      double laplace = 0.0, convection = 0.0;
      for (k = 1; k <= nsd; k++)
      {
        laplace += fe.grad(2)(i,k)*fe.grad(2)(j,k);
        convection += fe.basis(2)(i)*rhoc_adv*vel[k-1]*fe.grad(2)(j,k);
      }
      KTT(i,j) += scl2*scl2*(laplace*lambda + convection)*fe.detJxW;
    }
  }

  utl::zero_print_tol = 1e-30;

  // Integration of the thermal matrix
  Matrix CTT;
  CTT.resize(fe.basis(2).size(),fe.basis(2).size());
  // Evaluate the effective heat capacity
  double rhoc_eff = mat->getEffHeatCapacity(X,p,T);
  for (i = 1; i <= fe.basis(2).size(); i++)
    for (j = 1; j <= fe.basis(2).size(); j++)
      CTT(i,j) += scl2*scl2*fe.basis(2)(i)*(rhoc_eff + Lf*xi*ST)*fe.basis(2)(j)*fe.detJxW;

  elMat.A[TT] += CTT;
  elMat.A[TT].add(KTT,time.dt);

  // Add flow driving body forces vector to the RHS
  for (size_t i = 1; i <= fe.basis(2).size(); i++)
    for (size_t k = 1; k <= nsd; k++)
      elMat.b[Rp](i) += time.dt*fe.grad(2)(i,k)*(kr / (rhow*gacc))*perm[k-1]*rhow*grav[k-1]*fe.detJxW;

#if SP_DEBUG > 3
  std::cout << "Cup = " << elMat.A[up] << std::endl;
  std::cout << "CuT = " << elMat.A[uT] << std::endl;
  std::cout << "Cuu = " << elMat.A[uu] << std::endl;
  std::cout << "Kpp = " << Kpp << std::endl;
  std::cout << "Cpp = " << Cpp << std::endl;
  std::cout << "App = " << elMat.A[pp] << std::endl;
  std::cout << "CTp = " << CTp << std::endl;
  std::cout << "KTT = " << KTT << std::endl;
  std::cout << "CTT = " << CTT << std::endl;
  std::cout << "ATT = " << elMat.A[TT] << std::endl;
#endif

  // Evaluate elements of Kprev on basis 2
  size_t ru = elMat.A[uu].rows();
  size_t rp = elMat.A[pp].rows();

  for (i = 1; i <= rp; i++)
  {
    for (j = 1; j <= rp; j++)
    {
      size_t ki = ru+2*i-1;
      size_t kj = ru+2*j-1;
      elMat.A[Kprev](ki,kj) += Cpp(i,j);
      size_t li = ru+2*i;
      size_t lj = ru+2*j;
      elMat.A[Kprev](li,lj) += CTT(i,j);
    }
  }

  // If SUPG stabilization is required, evaluate the stabilizing matrices
  if (SUPG)
  {
    Matrix KTTs, CTTs, CTps;
    KTTs.resize(fe.basis(2).size(),fe.basis(2).size());
    CTTs.resize(fe.basis(2).size(),fe.basis(2).size());
    CTps.resize(fe.basis(2).size(),fe.basis(2).size());
    double h = StabilizationUtils::getElementSize(fe.XC,nsd);
    double alpha_e = h * fabs(nsd*vel[0]) / (2*lambda);
    double taue;
    if (alpha_e > 1)
      taue = h / (2 * fabs(nsd*vel[0])) * (1/tanh(alpha_e) - 1/alpha_e);
    else
      taue = h * h / (12*lambda);

#if SP_DEBUG > 3
    std::cout << "Element size,                      h = " << h << std::endl;
    std::cout << "Peclet like number,          alpha_e = " << alpha_e << std::endl;
    std::cout << "SUPG stabilization coefficient, taue = " << taue << std::endl;
#endif

    for (i = 1; i <= fe.basis(2).size(); i++) {
      for (j = 1; j <= fe.basis(2).size(); j++) {
        double a = 0.0, b = 0.0, c = 0.0;
        for (k = 1; k <= nsd; k++) {
          a += fe.grad(2)(i,k)*vel[k-1]*vel[k-1]*fe.grad(2)(j,k);
          b += fe.grad(2)(i,k)*vel[k-1]*fe.hess(2)(i,k,k);
          c += fe.grad(2)(i,k)*vel[k-1]*fe.basis(2)(j);
        }
        KTTs(i,j) += scl2 * scl2 * taue * rhoc_adv * (rhoc_adv * a + lambda * b) * fe.detJxW;
        CTTs(i,j) += scl2 * scl2 * taue * rhoc_adv * (rhoc_eff * c + Lf * xi * ST * c) * fe.detJxW;
        CTps(i,j) += scl1 * scl2 * taue * rhoc_adv * Lf * xi * Sp * c * fe.detJxW;
      }
    }
    elMat.A[TT] += CTTs;
    elMat.A[TT].add(KTTs,time.dt);
    for (i = 1; i <= rp; i++) {
      for (j = 1; j <= rp; j++) {
        size_t li = ru+2*i;
        size_t lj = ru+2*j;
        elMat.A[Kprev](li,lj) += CTTs(i,j);
      }
    }
  }

  return true;
}


bool THMCoupledNL::evalBouMx(LocalIntegral& elmInt,
                             const MxFiniteElement& fe,
                             const TimeDomain& time,
                             const Vec3& X, const Vec3& normal) const
{
  if (!tracFld && !fluxFld)
  {
    std::cerr << "*** THMCoupledNL::evalBouMx: No fluxes/tractions." << std::endl;
    return false;
  }

  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  //! Evaluate the surface traction
  Vec4 Xt = static_cast<const Vec4&>(X);
  Xt.t = time.t;
  Vec3 tr2 = this->getTraction(Xt,normal);
  Xt.t -= time.dt;
  Vec3 tr1 = this->getTraction(Xt,normal);
  Vec3 dtr;
  dtr = tr2 - tr1;

  // Integration of Ru
  for (size_t i = 1; i <= fe.basis(1).size(); i++)
    for (unsigned short int j = 1; j <= nsd; j++)
      elMat.b[Ru](nsd*(i-1)+j) += 1.0*dtr[j-1]*fe.basis(1)(i)*fe.detJxW;

  // Integration of Rp
  for (size_t i = 1; i <= fe.basis(2).size(); i++)
    for (size_t k = 1; k <= nsd; k++)
      elMat.b[Rp](i) += 0.0;    // TO DO

  // Integration of RT
  for (size_t i = 1; i <= fe.basis(2).size(); i++)
    for (size_t k = 1; k <= nsd; k++)
      elMat.b[RT](i) += 0.0;    // TO DO

  return true;
}


bool THMCoupledNL::finalizeElement(LocalIntegral& elmInt,
                                   const TimeDomain&, size_t)
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  size_t ru = elMat.A[uu].rows();
  size_t rp = elMat.A[pp].rows();

  // Evaluate the previous and updated internal matrices
  for (size_t i = 1; i <= ru; i++)
  {
    for (size_t j = 1; j <= ru; j++)
    {
      elMat.A[Kprev](i,j) = elMat.A[uu](i,j);
      elMat.A[Know](i,j)  = elMat.A[uu](i,j);
    }
    for (size_t j = 1; j <= rp; j++)
    {
      size_t k = ru+2*j-1;
      elMat.A[Kprev](i,k) = elMat.A[up](i,j);
      elMat.A[Know](i,k)  = elMat.A[up](i,j);
      elMat.A[Kprev](k,i) = elMat.A[pu](j,i);
      elMat.A[Know](k,i)  = elMat.A[pu](j,i);
      size_t l = ru+2*j;
      elMat.A[Kprev](i,l) = elMat.A[uT](i,j);
      elMat.A[Know](i,l)  = elMat.A[uT](i,j);
    }
  }

  for (size_t i = 1; i <= rp; i++)
  {
    for (size_t j = 1; j <= rp; j++)
    {
      size_t ki = ru+2*i-1;
      size_t kj = ru+2*j-1;
      size_t li = ru+2*i;
      size_t lj = ru+2*j;
      elMat.A[Know](ki,kj)  = elMat.A[pp](i,j);
      elMat.A[Know](li,lj)  = elMat.A[TT](i,j);
      elMat.A[Kprev](ki,lj) = elMat.A[pT](i,j);
      elMat.A[Kprev](li,kj) = elMat.A[Tp](i,j);
      elMat.A[Know](ki,lj)  = elMat.A[pT](i,j);
      elMat.A[Know](li,kj)  = elMat.A[Tp](i,j);
    }
  }

  // Get the previous and current solution vectors on basis 2
  Vector PoTo, PcTc;
  PoTo.resize(2*elMat.vec[Po].size());
  PcTc.resize(2*elMat.vec[Pc].size());
  for (size_t i = 1; i <= elMat.vec[Po].size(); i++)
  {
    PoTo(2*i-1) = elMat.vec[Po](i);
    PoTo(2*i  ) = elMat.vec[To](i);
    PcTc(2*i-1) = elMat.vec[Pc](i);
    PcTc(2*i  ) = elMat.vec[Tc](i);
  }

  // Get the previous and current solution vectors
  Vector prevSol, currSol;
  prevSol = elMat.vec[Uo];
  prevSol.insert(prevSol.end(),PoTo.begin(),PoTo.end());

  currSol = elMat.vec[Uc];
  currSol.insert(currSol.end(),PcTc.begin(),PcTc.end());

  elMat.b[Rprev] = elMat.A[Kprev]*prevSol;
  elMat.b[Rnow]  = elMat.A[Know]*currSol;

  return true;
}


bool THMCoupledNL::evalSol(Vector& s, const MxFiniteElement& fe,
                           const Vec3& X, const std::vector<int>& MNPC,
                           const std::vector<size_t>& elem_sizes,
                           const std::vector<size_t>& basis_sizes) const
{
  std::vector<int>::const_iterator fstart = MNPC.begin() + elem_sizes[0];
  IntVec MNPC2(fstart,MNPC.end());
  for (size_t i = 0; i < MNPC2.size(); i++)
    MNPC2[i] += basis_sizes[0];

  Vector disp, pres, temp;
  int ierr = utl::gather(IntVec(MNPC.begin(), fstart),nsd,primsol[0],disp)
           + utl::gather(MNPC2,0,2,primsol[0],pres,nsd*basis_sizes[0],basis_sizes[0])
           + utl::gather(MNPC2,1,2,primsol[0],temp,nsd*basis_sizes[0],basis_sizes[0]);

  if (ierr != 0)
    std::cerr << " *** THMCoupledNL::initElement: Detected " << ierr/3
              << " node numbers out of range." << std::endl;

  double p = pres.dot(fe.basis(2));
  p *= scl1;
  double T = temp.dot(fe.basis(2));
  T *= scl2;

  double pi = mat->getIcePressure(X,p,T);
  double Si = mat->getIceSaturation(X,p,T);
  double Sw = mat->getWaterSaturation(X,p,T);

  s.resize(9);
  s(1) = pi;
  s(2) = Si;
  s(3) = Sw;

  Matrix Bmat, Cmat;
  if(!mat->formBmatrix(Bmat,fe.grad(1),nsd))
    return false;

  if(!mat->formElasticMatrix(Cmat,X,nsd,p,T))
    return false;

  Vector strain, stress;
  strain = Bmat*disp;
  stress = Cmat*strain;

  for (size_t i = 1; i <= strain.size(); i++) {
    s(3+i) = strain(i);
    s(6+i) = stress(i);
  }

  return true;
}


size_t THMCoupledNL::getNoFields(int fld) const
{
  size_t noSecSol = 9;
  if (fld < 2)
    return nsd+2;
  else
    return noSecSol;
}


const char* THMCoupledNL::getField1Name(size_t i, const char* prefix) const
{
  static const char* s[4] = { "u_x", "u_y", "p^w", "T" };

  if(!prefix) return s[i];

  static std::string name;
  name = prefix + std::string(" ") + s[i];

  return name.c_str();
}


const char* THMCoupledNL::getField2Name(size_t i, const char* prefix) const
{
  size_t noSecSol = 9;
  if (i >= noSecSol) return 0;

  static const char* s2[] = {"p^i","Si","Sw","eps_x","eps_y","eps_xy",
                             "sig_x", "sig_y","sig_xy"};
  if (!prefix) return s2[i];

  static std::string name;
  name = prefix + std::string(" ") + s2[i];

  return name.c_str();
}


THMCoupledNL::WeakDirichlet::WeakDirichlet(unsigned short int n) :
    nsd(n), flux(nullptr)
{
  primsol.resize(2);
}

LocalIntegral* THMCoupledNL::WeakDirichlet::getLocalIntegral(const std::vector<size_t>& nen,
                                                             size_t, bool neumann) const
{
  const size_t nedof1 = nsd*nen[0];           //!< Number of DOFs on basis 1
  const size_t nedof  = nedof1 + 2*nen[1];    //!< Total number of DOFs

  ElmMats* result = new MixedElmMats();

  //result->rhsOnly = false;
  result->withLHS = true;
  result->b[Ru].resize(nedof1);
  result->b[Rp].resize(nen[1]);
  result->b[RT].resize(nen[1]);
  result->b[Rprev].resize(nedof);
  result->b[Rnow].resize(nedof);
  result->b[Rres].resize(nedof);
  result->b[Rc].resize(nedof);

  result->A[uu].resize(nedof1,nedof1);
  result->A[up].resize(nedof1,nen[1]);
  result->A[uT].resize(nedof1,nen[1]);
  result->A[pu].resize(nen[1],nedof1);
  result->A[pp].resize(nen[1],nen[1]);
  result->A[pT].resize(nen[1],nen[1]);
  result->A[Tp].resize(nen[1],nen[1]);
  result->A[TT].resize(nen[1],nen[1]);
  result->A[Kprev].resize(nedof,nedof);
  result->A[Know].resize(nedof,nedof);
  result->A[Ktan].resize(nedof,nedof);

  return result;
}


bool THMCoupledNL::WeakDirichlet::initElementBou(const std::vector<int>& MNPC,
                                                 const std::vector<size_t>& elem_sizes,
                                                 const std::vector<size_t>& basis_sizes,
                                                 LocalIntegral& elmInt)
{
  if(primsol.front().empty()) return true;

  // Extract element level solution vectors
  elmInt.vec.resize(NSOL);
  std::vector<int>::const_iterator fstart = MNPC.begin() + elem_sizes[0];
  int ierr = utl::gather(IntVec(MNPC.begin(), fstart),nsd,primsol[0],elmInt.vec[Uc])
           + utl::gather(IntVec(fstart,MNPC.end()),0,2,primsol[0],elmInt.vec[Pc],nsd*basis_sizes[0],basis_sizes[0])
           + utl::gather(IntVec(fstart,MNPC.end()),1,2,primsol[0],elmInt.vec[Tc],nsd*basis_sizes[0],basis_sizes[0])
           + utl::gather(IntVec(MNPC.begin(), fstart),nsd,primsol[1],elmInt.vec[Uo])
           + utl::gather(IntVec(fstart,MNPC.end()),0,2,primsol[1],elmInt.vec[Po],nsd*basis_sizes[0],basis_sizes[0])
           + utl::gather(IntVec(fstart,MNPC.end()),1,2,primsol[1],elmInt.vec[To],nsd*basis_sizes[0],basis_sizes[0]);

  if (ierr == 0) return true;

  std::cerr << " *** THMCoupledNL::initElement: Detected " << ierr/3
            << " node numbers out of range." << std::endl;

  return false;
}


bool THMCoupledNL::WeakDirichlet::evalBouMx(LocalIntegral& elmInt,
                                            const MxFiniteElement& fe,
                                            const TimeDomain& time,
                                            const Vec3& X,
                                            const Vec3& normal) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  // Evaluate the Neumann heat flux on the boundary
  double qT = 0.0;
  if (flux)
    qT = (*flux)(X);

  size_t ru = elMat.A[uu].rows();

  for (size_t i = 1; i <= fe.basis(2).size(); i++)
  {
    for (size_t j = 1; j <= fe.basis(2).size(); j++)
    {
      size_t li = ru+2*i;
      size_t lj = ru+2*j;
      elMat.A[TT](i,j) += time.dt *fe.basis(2)(i)*lambdae*fe.basis(2)(j)*fe.detJxW;
      elMat.A[Know](li,lj) += time.dt *fe.basis(2)(i)*lambdae*fe.basis(2)(j)*fe.detJxW;
    }
    elMat.b[RT](i) += time.dt * (qT*fe.basis(2)(i) + lambdae*Te*fe.basis(2)(i)) * fe.detJxW;
  }

  // Get the previous and current solution vectors on basis 2
  Vector PcTc;
  PcTc.resize(2*elMat.vec[Pc].size());
  for (size_t i = 1; i <= elMat.vec[Pc].size(); i++)
  {
    PcTc(2*i-1) = elMat.vec[Pc](i);
    PcTc(2*i  ) = elMat.vec[Tc](i);
  }

  elMat.b[Rc] = elMat.vec[Uc];
  elMat.b[Rc].insert(elMat.b[Rc].end(),PcTc.begin(),PcTc.end());

  return true;
}

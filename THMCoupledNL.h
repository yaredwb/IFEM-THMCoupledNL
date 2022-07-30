// $Id$
//==============================================================================
//!
//! \file THMCoupledNL.h
//!
//! \date
//!
//! \author Yared Bekele
//!
//! \brief Integrand implementations for THM coupled problems in ground freezing
//!
//==============================================================================

#ifndef _THM_COUPLED_NL_H
#define _THM_COUPLED_NL_H

#include "Vec3.h"
#include "ElmMats.h"
#include "IntegrandBase.h"
#include "Poro3Material.h"

/*!
 * \brief Class representing the integrands of fully coupled thermo-hydro-mechanical
    processes during ground freezing. A general class name 'THMCoupledNL' is chosen
    with a view to extending the implementation for all practical temperature ranges
    of interest.
*/

class THMCoupledNL : public IntegrandBase
{
  /*!
   * \brief Class representing an element matrix for fully coupled THM processes
     during ground freezing.
  */
  class MixedElmMats : public ElmMats
  {
  public:
    //! \brief Default constructor
    MixedElmMats();
    //! \brief Empty destructor
    virtual ~MixedElmMats() {}
    //! \brief Returns the element level Newton matrix
    virtual const Matrix& getNewtonMatrix() const;
    //! \brief Returns the element level RHS vector
    virtual const Vector& getRHSVector() const;
  };

public:
  /*!
   * \brief Class for integrating Robin boundary conditions
  */
  class WeakDirichlet : public IntegrandBase
  {
  public:
    //! \brief Default constructor.
    //! \param[in] n Number of spatial dimensions
    WeakDirichlet(unsigned short int n);

    //! \brief Empty destructor.
    virtual ~WeakDirichlet() {}

    //! \brief Defines the flux function
    void setFlux(RealFunc* f) { flux = f; }
    //! \brief Defines the temperature of the environment.
    void setEnvtTemperature(double envT) { Te = envT; }
    //! \brief Defines the thermal conductivity of the environment.
    void setEnvtConductivity(double envCond) { lambdae = envCond; }

    //! \brief Returns that this integrand has no interior contributions.
    virtual bool hasInteriorTerms() const { return false; }

    //! \brief Returns a local integral container for the given element.
    //! \param[in] nen Number of nodes on element
    //! \param[in] neumann Whether or not we are assembling Neumann BCs
    virtual LocalIntegral* getLocalIntegral(const std::vector<size_t>& nen,
                                            size_t, bool neumann) const;

    //! \brief Initializes current element boundary for numerical integration.
    //! \param[in] MNPC Matrix of nodal point correspondance for current element
    //! \param[in] elem_sizes Size of each basis on the element
    //! \param[in] basis_sizes Size of each basis on the patch
    //! \param elmInt The local integral object for the current element
    virtual bool initElementBou(const std::vector<int>& MNPC,
                                const std::vector<size_t>& basis_sizes,
                                const std::vector<size_t>& elem_sizes,
                                LocalIntegral& elmInt);

    //! \brief Evaluates the integrands at a boundary point.
    //! \param elmInt The local integral object to receive the contributions
    //! \param[in] fe Finite element data of current integration point
    //! \param[in] time Parameters for nonlinear and time-dependent simulations
    //! \param[in] X Cartesian coordinates of current integration point
    //! \param[in] normal Boundary normal vector at integration point
    virtual bool evalBouMx(LocalIntegral& elmInt, const MxFiniteElement& fe,
                           const TimeDomain& time, const Vec3& X,
                           const Vec3& normal) const;

    //! \brief Finalizes the element quantities after the numerical integration.
    //! \details This method is invoked once for each element, afte the numerical
    //! integration loop over interior points is finished and before the resulting
    //! element quantities are assembled into their system level equivalents
    //virtual bool finalizeElement(LocalIntegral&, const TimeDomain&, size_t);

  protected:
    unsigned short int nsd;       //!< Number of space dimensions
    RealFunc* flux;               //!< Flux function
    double Te;                    //!< Temperature of environment
    double lambdae;               //!< Thermal conductivity of environment
  };

  //! \brief The default constructor initializes all pointers to zero.
  //! \param[in] n Number of spatial dimensions
  //! \param[in] order The order of the time integration
  //! \param[in] stab Boolean flag SUPG stabilization
  THMCoupledNL(unsigned short int n = 2, int order = 1, bool stab = false);

  //! \brief The destructor frees the dynamically allocated data objects.
  virtual ~THMCoupledNL() {}

  //! \brief Defines the traction field to use in Neumann boundary conditions.
  void setTraction(TractionFunc* tf) { tracFld = tf; }

  //! \brief Defines the traction field to use in Neumann boundary conditions.
  void setTraction(VecFunc* tf) { fluxFld = tf; }

  //! \brief Defines scaling values for U and P and U and T
  void setScalingValues(double scl_UP, double scl_UT)
  { scl1 = scl_UP; scl2 = scl_UT; }

  //! \brief Defines the graviation vector
  //! \param[in] grav Gravity vector
  void setGravity(const Vec3& gravity) { grav = gravity; }

  //! \brief Evaluates the boundary traction field (if any) at a specified point.
  //! \param[in] X Cartesian coordinate of the current integration point
  //! \param[in] n Outward-directed unit normal vector at current point
  virtual Vec3 getTraction(const Vec3& X, const Vec3& n) const;

  //! Get the gravity vector
  const Vec3 getGravity() const { return grav; }

  //! \brief Defines the material properties.
  virtual void setMaterial(Poro3Material* material) { mat = material; }

  virtual int getIntegrandType() const
  { return SUPG ? ELEMENT_CORNERS | SECOND_DERIVATIVES : STANDARD; }

  //! \brief Returns a local integral container for the given element.
  //! \param[in] nen Number of nodes on element
  //! \param[in] neumann Whether or not we are assembling Neumann BCs
  virtual LocalIntegral* getLocalIntegral(const std::vector<size_t>& nen,
                                          size_t, bool neumann) const;

  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param[in] elem_sizes Size of each basis on the element
  //! \param[in] basis_sizes Size of each basis on the patch
  //! \param elmInt The local integral object for current element
  virtual bool initElement(const std::vector<int>& MNPC,
                           const std::vector<size_t>& basis_sizes,
                           const std::vector<size_t>& elem_sizes,
                           LocalIntegral& elmInt);

  //! \brief Initializes current element boundary for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param[in] elem_sizes Size of each basis on the element
  //! \param[in] basis_sizes Size of each basis on the patch
  //! \param elmInt The local integral object for current element
  virtual bool initElementBou(const std::vector<int>& MNPC,
                              const std::vector<size_t>& basis_sizes,
                              const std::vector<size_t>& elem_sizes,
                              LocalIntegral& elmInt);

  //! \brief Evaluates the integrands at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integraion point
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalIntMx(LocalIntegral& elmInt, const MxFiniteElement& fe,
                         const TimeDomain& time, const Vec3& X) const;

  //! \brief Evaluates the integrands at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] time Parameters for nonlineat and time-dependent simulations
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at integration point
  virtual bool evalBouMx(LocalIntegral& elmInt, const MxFiniteElement& fe,
                         const TimeDomain& time, const Vec3& X,
                         const Vec3& normal) const;

  //! \brief Finalizes the element quantities after the numerical integration.
  //! \details This method is invoked once for each element, afte the numerical
  //! integration loop over interior points is finished and before the resulting
  //! element quantities are assembled into their system level equivalents
  virtual bool finalizeElement(LocalIntegral&, const TimeDomain&, size_t);

  //! \brief Evalautes the secondary solution at a result point (mixed).
  //! \param[out] s The solution field values at current point
  //! \param[in] fe Mixed finte element data at current point
  //! \param[in] X Cartesian coordiantes of current point
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param[in] elem_sizes Size of each basis on the element
  //! \param[in] basis_sizes Size of each basis on the patch
  virtual bool evalSol(Vector& s, const MxFiniteElement& fe, const Vec3& X,
                       const std::vector<int>& MNPC,
                       const std::vector<size_t>& elem_sizes,
                       const std::vector<size_t>& basis_sizes) const;

  //! \brief Returns whether a mixed formulation is used.
  virtual bool mixedFormulation() const { return true; }

  //! \brief Returns the number of primary/secondary solution field components.
  //! \param[in] fld Which field set to consider (1 = Primary, 2 = Secondary)
  virtual size_t getNoFields(int fld = 2) const;

  //! \brief Returns the name of a primary solution field component.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  virtual const char* getField1Name(size_t i, const char* prefix = 0) const;

  //! \brief Returns the name of a secondary solution field component.
  //! \param[in] i Field componet index
  //! \param[in] prefix Name prefix for all components
  virtual const char* getField2Name(size_t i, const char* prefix = 0) const;

private:
  Vec3 grav;                    //!< Gravitation vector

protected:
  unsigned short int nsd;       //!< Number of space dimensions
  double gacc;                  //!< Gravitational acceleration
  TractionFunc* tracFld;        //!< Pointer to implicit boundary traction field
  VecFunc* fluxFld;             //!< Pointer to explicit boundary traction field
  Poro3Material* mat;           //!< Material data
  double scl1;                  //!< Scaling between displacement and pressure
  double scl2;                  //!< Scaling between displacement and temperature
  bool SUPG;                    //!< Boolean flag for SUPG stabilization
};

#endif // _THM_COUPLED_NL_H

#ifndef OMEGA_FRAZIL_H
#define OMEGA_FRAZIL_H
//===-- ocn/Frazil.h - Frazil Ice Formation -------------------*- C++ -*-===//
//
// This header defines a scaffold for frazil-related tendencies and
// accumulators. Physics implementations are intentionally left empty.
//
//===----------------------------------------------------------------------===//

#include "Config.h"
// #include "DataTypes.h"
#include "GlobalConstants.h"
#include "HorzMesh.h"
#include "OmegaKokkos.h"
#include "VertCoord.h"

#include <map>
#include <memory>
#include <string>

#include <gswteos-10.h>

namespace OMEGA {

enum class FrazilType {
   BasicFrazil,  ///< Placeholder basic frazil option
   SimpleFrazil, ///< Placeholder simple frazil option
   TeosFrazil    ///< Placeholder TEOS frazil option
};

class FrazilMelt {
 public:
   /// constructor declaration
   FrazilMelt();

   //   The functor takes the full arrays of specific volume (inout),
   //   the indices ICell and KChunk, and the ocean tracers (conservative)
   //   temperature, and (absolute) salinity as inputs, and outputs the
   //   specific volume according to the Roquet et al. 2015 75 term expansion.
   KOKKOS_FUNCTION void operator()(Real &AccMIce, Real &AccMLiq, Real &AccMSalt,
                                   Real &AccELiq, Real &AccEIce, Real &HTend,
                                   Real &TTend, Real &STend) const {
      (void)AccMIce;
      (void)AccMLiq;
      (void)AccMSalt;
      (void)AccELiq;
      (void)AccEIce;
      HTend = 0.0_Real;
      TTend = 0.0_Real;
      STend = 0.0_Real;
   }
};

class FrazilFormation {
 public:
   /// constructor declaration
   FrazilFormation();

   /// Parameters for FrazilFormation (overwritten by config file if set there)
   Real Phi = 0.75_Real; ///< liquid mass fraction of frazil for export
   Real MassLimit =
       0.1_Real; ///<  layer mass fraction limit for thickness tendency

   //   The functor takes the full arrays of specific volume (inout),
   //   the indices ICell and KChunk, and the ocean tracers (conservative)
   //   temperature, and (absolute) salinity as inputs, and outputs the
   //   specific volume according to the Roquet et al. 2015 75 term expansion.
   KOKKOS_FUNCTION void operator()(const Real SA, const Real CT, const Real P,
                                   const Real h, const Real Dt, Real &AccMIce,
                                   Real &AccMLiq, Real &AccMSalt, Real &AccELiq,
                                   Real &AccEIce, Real &HTend, Real &TTend,
                                   Real &STend) const {

      Real CTnew;
      Real SAnew;
      Real wIh        = 0.0_Real;
      Real solidMass  = 0.0_Real;
      Real liquidMass = 0.0_Real;

      // type casting for now because gsw expects double
      // and Real can be either float or double depending on build configuration
      double SAnew_d = 0.0;
      double CTnew_d = 0.0;
      double wIh_d   = 0.0;

      gsw_frazil_properties_potential_poly(
          static_cast<double>(SA), static_cast<double>(Cp0Sw * CT),
          static_cast<double>(P), &SAnew_d, &CTnew_d, &wIh_d);

      SAnew = static_cast<Real>(SAnew_d);
      CTnew = static_cast<Real>(CTnew_d);
      wIh   = static_cast<Real>(wIh_d);

      const Real OneMinusPhi = Kokkos::max(1.0e-12_Real, 1.0_Real - Phi);
      solidMass              = h * Kokkos::min(wIh, OneMinusPhi * MassLimit);
      liquidMass             = (Phi / OneMinusPhi) * solidMass;

      HTend = -(solidMass + liquidMass) / Dt;
      TTend = (h - solidMass - liquidMass) * CTnew - h * CT;
      STend = (h - solidMass - liquidMass) * SAnew - h * SA;

      // these are not currently mass or energy, they all need a RhoSw factor
      AccMIce += solidMass;
      AccMLiq += liquidMass;
      AccMSalt +=
          liquidMass * SAnew * PPt2Salt; // the PPt2Salt scaling could be moved
                                         // to the coupler interaction
      AccELiq += liquidMass * Cp0Sw * CTnew;
      AccEIce += solidMass * gsw_pot_enthalpy_from_pt_ice_poly(CTnew);
   }
};

class Frazil {
 public:
   static void init();
   /// Creates a new frazil object and stores it in the AllFrazil map.
   static Frazil *create(const std::string &Name);

   /// Retrieve frazil object by name.
   static Frazil *get(const std::string &Name);

   /// Retrieve default frazil object.
   static Frazil *getDefault();

   /// Destructor
   ~Frazil();

   /// Deallocates arrays
   static void clear();

   /// Remove frazil object by name.
   static void erase(std::string InName); ///< [in] name to remove

   Array2DReal FrazilTTend;
   Array2DReal FrazilSTend;
   Array2DReal FrazilHTend;
   Array1DReal AccMIce;
   Array1DReal AccEIce;
   Array1DReal AccMLiq;
   Array1DReal AccELiq;
   Array1DReal AccMSalt;

   void computeFrazil(const Array2DReal &CT, const Array2DReal &SA,
                      const Array2DReal &P, const Array2DReal &H);

 private:
   static Frazil *DefaultFrazil;
   static std::map<std::string, std::unique_ptr<Frazil>> AllFrazil;

   Frazil(const HorzMesh *Mesh, const VertCoord *VCoord);

   // Forbid copy and move construction/assignment.
   Frazil(const Frazil &)            = delete;
   Frazil &operator=(const Frazil &) = delete;
   Frazil(Frazil &&)                 = delete;
   Frazil &operator=(Frazil &&)      = delete;

   FrazilType frazilChoice;
   FrazilFormation computeFrazilFormation;
   FrazilMelt computeFrazilMelt;
   Real massLimit;
   Real phi;
   //    Real dt;
   I4 NCellsAll;
   I4 NChunks;

   const HorzMesh *MeshPtr;
   const VertCoord *VCoordPtr;
};

} // namespace OMEGA

#endif

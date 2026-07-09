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
   BasicFrazil,  ///< MPAS-O style basic frazil option
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
   KOKKOS_FUNCTION void operator()(const Real SA, const Real CT, const Real P,
                                   const Real h, Real &AccMIce, Real &AccMLiq,
                                   Real &AccMSalt, Real &AccELiq, Real &AccEIce,
                                   Real &HTend, Real &TTend,
                                   Real &STend) const {

      constexpr Real Eps = 1.0e-12_Real;

      if (AccMIce <= Eps) { // potential leak if we dont redistribute
         AccMIce = 0.0_Real;
         HTend   = 0.0_Real;
         TTend   = 0.0_Real;
         STend   = 0.0_Real;
         return;
      }

      const Real safeAccMLiq       = Kokkos::max(AccMLiq, Eps);
      const Real safeAccMIce       = Kokkos::max(AccMIce, Eps);
      const Real frazilIceFraction = Kokkos::min(
          1.0_Real,
          Kokkos::max(0.0_Real, AccMIce / (safeAccMLiq + safeAccMIce)));
      const Real brineSalinity  = AccMSalt / safeAccMLiq;
      const Real brineEnthalpy  = AccELiq / safeAccMLiq;
      const Real potEnthalpyIce = AccEIce / safeAccMIce;

      const Real layerMass     = h + AccMIce;
      const Real safeLayerMass = Kokkos::max(layerMass, Eps);

      const Real layerIceFraction =
          Kokkos::min(1.0_Real, Kokkos::max(0.0_Real, AccMIce / safeLayerMass));

      // typecasting for now but will be simplified
      const double SA_d    = static_cast<double>(SA);
      const double CT_d    = static_cast<double>(CT);
      const double P_d     = static_cast<double>(P);
      const double wIhIn_d = static_cast<double>(layerIceFraction);
      const double pt0Ice_d =
          gsw_pt_from_pot_enthalpy_ice(static_cast<double>(potEnthalpyIce));
      const double tIce_d = gsw_t_from_pt0_ice(pt0Ice_d, P_d);
      double SAnew_d      = SA_d;
      double CTnew_d      = CT_d;
      double wIhOut_d     = wIhIn_d;

      gsw_melting_ice_into_seawater(SA_d, CT_d, P_d, wIhIn_d, tIce_d, &SAnew_d,
                                    &CTnew_d, &wIhOut_d);

      const Real wIhOut = Kokkos::min(
          1.0_Real, Kokkos::max(0.0_Real, static_cast<Real>(wIhOut_d)));

      const Real finalSolidMass =
          Kokkos::min(AccMIce, Kokkos::max(0.0_Real, wIhOut * safeLayerMass));
      const Real solidMass = Kokkos::max(0.0_Real, AccMIce - finalSolidMass);

      if (solidMass <= Eps) {
         return;
      }

      const Real liquidMass =
          Kokkos::min(AccMLiq, solidMass * (1.0_Real - frazilIceFraction) /
                                   Kokkos::max(frazilIceFraction, Eps));
      const Real solidEnthalpy  = solidMass * potEnthalpyIce;
      const Real liquidEnthalpy = liquidMass * brineEnthalpy;

      HTend = +(solidMass + liquidMass);
      TTend = +(liquidEnthalpy + solidEnthalpy) / Cp0Sw;
      STend = +(liquidMass * brineSalinity);

      AccMIce  = Kokkos::max(0.0_Real, AccMIce - solidMass);
      AccMLiq  = Kokkos::max(0.0_Real, AccMLiq - liquidMass);
      AccMSalt = Kokkos::max(0.0_Real, AccMSalt - liquidMass * brineSalinity);
      AccELiq -= liquidEnthalpy;
      AccEIce -= solidEnthalpy;
      if (AccMIce <= Eps) {
         AccEIce = 0.0_Real;
      }
      if (AccMLiq <= Eps) {
         AccELiq = 0.0_Real;
      }
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
                                   const Real h, Real &AccMIce, Real &AccMLiq,
                                   Real &AccMSalt, Real &AccELiq, Real &AccEIce,
                                   Real &HTend, Real &TTend,
                                   Real &STend) const {

      Real CTnew;
      Real SAnew;
      Real wIh            = 0.0_Real;
      Real solidMass      = 0.0_Real;
      Real liquidMass     = 0.0_Real;
      Real solidEnthalpy  = 0.0_Real;
      Real liquidEnthalpy = 0.0_Real;

      // type casting for now because gsw expects double
      // and Real can be either float or double depending on build configuration
      double SAnew_d = 0.0;
      double CTnew_d = 0.0;
      double wIh_d   = 0.0;
      // double PTnew_d = 0.0;

      gsw_frazil_properties_potential_poly(
          static_cast<double>(SA), static_cast<double>(Cp0Sw * CT),
          static_cast<double>(P), &SAnew_d, &CTnew_d, &wIh_d);
      double PTnew_d =
          gsw_pt_from_ct(SAnew_d, CTnew_d); // convert to potential temperature

      SAnew = static_cast<Real>(SAnew_d);
      CTnew = static_cast<Real>(CTnew_d);
      wIh   = static_cast<Real>(wIh_d);

      const Real OneMinusPhi = Kokkos::max(1.0e-12_Real, 1.0_Real - Phi);
      // anything called mass below is in pseudo-thickness units (m) and needs
      // to be scaled by RhoSw for coupling
      solidMass      = h * Kokkos::min(wIh, OneMinusPhi * MassLimit);
      liquidMass     = (Phi / OneMinusPhi) * solidMass;
      solidEnthalpy  = solidMass * gsw_pot_enthalpy_from_pt_ice_poly(PTnew_d);
      liquidEnthalpy = liquidMass * Cp0Sw * CTnew;

      // per timestep (not scaled by dt here)
      HTend = -(solidMass + liquidMass);
      // TTend = (h - solidMass - liquidMass) * CTnew - h * CT;
      // STend = (h - solidMass - liquidMass) * SAnew - h * SA;
      TTend = -(liquidEnthalpy + solidEnthalpy) /
              Cp0Sw; // convert back to CT tendency
      STend = -(liquidMass * SAnew);

      // Local unit of mass is pseudo thickness (m)
      // these all need a RhoSw factor before coupling
      AccMIce += solidMass;           // m
      AccMLiq += liquidMass;          // m
      AccMSalt += liquidMass * SAnew; // (m)(g/kg)
      AccELiq += liquidEnthalpy;      // (m)(J/kg)
      AccEIce += solidEnthalpy;       // (m)(J/kg)
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
   bool conservationCheck = false;
   Real depthLimit        = -1.0_Real;

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
   BasicFrazilFormation computeBasicFrazilFormation;
   BasicFrazilMelt computeBasicFrazilMelt;
   FrazilFormation computeFrazilFormation;
   FrazilMelt computeFrazilMelt;
   Real massLimit;
   Real phi;
   I4 NCellsAll;
   I4 NChunks;

   const HorzMesh *MeshPtr;
   const VertCoord *VCoordPtr;

   void checkColumnConservation() const;
};

} // namespace OMEGA

#endif

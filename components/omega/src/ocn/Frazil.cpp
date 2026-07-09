//===-- ocn/Frazil.cpp - Frazil Ice Formation -----------------*- C++ -*-===//
//
// The Frazil class manages frazil tendencies and accumulators.
// This initial implementation only has a teos-10 configuration.
//
//===----------------------------------------------------------------------===//

#include "Frazil.h"
#include "Eos.h"
#include "Error.h"
#include "Logging.h"

namespace OMEGA {

namespace {

KOKKOS_INLINE_FUNCTION bool isApprox(const Real X, const Real Y,
                                     const Real RTol, const Real ATol = 0) {
   if (Kokkos::isnan(X) || Kokkos::isnan(Y) || Kokkos::isinf(X) ||
       Kokkos::isinf(Y)) {
      return false;
   }

   return Kokkos::abs(X - Y) <=
          Kokkos::max(ATol, RTol * Kokkos::max(Kokkos::abs(X), Kokkos::abs(Y)));
}

} // namespace

Frazil *Frazil::DefaultFrazil = nullptr;
std::map<std::string, std::unique_ptr<Frazil>> Frazil::AllFrazil;

/// Constructor for FrazilFormation
FrazilFormation::FrazilFormation() {}

/// Constructor for FrazilMelt
FrazilMelt::FrazilMelt() {}

/// Constructor for BasicFrazilFormation
BasicFrazilFormation::BasicFrazilFormation() {}

/// Constructor for BasicFrazilMelt
BasicFrazilMelt::BasicFrazilMelt() {}

void Frazil::init() {

   if (!HorzMesh::getDefault() or !VertCoord::getDefault()) {
      ABORT_ERROR("Frazil::init: HorzMesh and VertCoord must be initialized");
   }

   // Frazil freezing-temperature calculations depend on EOS configuration.
   Eos::init();

   if (!DefaultFrazil) {
      Error Err;
      bool FrazilTendencyEnable = false;
      Config *OmegaConfig       = Config::getOmegaConfig();
      Config TendConfig("Tendencies");

      Err += OmegaConfig->get(TendConfig);
      CHECK_ERROR_ABORT(Err,
                        "Frazil::init: Tendencies group not found in Config");

      Err += TendConfig.get("FrazilTendencyEnable", FrazilTendencyEnable);
      CHECK_ERROR_ABORT(
          Err, "Frazil::init: FrazilTendencyEnable not found in Tendencies");

      if (!FrazilTendencyEnable) {
         LOG_INFO("Frazil::init: Frazil tendency disabled; skipping default "
                  "frazil object creation");
         LOG_INFO("All frazil is off - frazil parameters will be ignored");
         return;
      }

      DefaultFrazil = create("Default");
   }
}

Frazil::Frazil(const HorzMesh *Mesh, const VertCoord *VCoord)
    : frazilChoice(FrazilType::TeosFrazil), computeFrazilFormation(),
      computeFrazilMelt(), NCellsAll(Mesh->NCellsAll),
      NChunks((VCoord->NVertLayers + VecLength - 1) / VecLength), MeshPtr(Mesh),
      VCoordPtr(VCoord) {

   FrazilTTend =
       Array2DReal("FrazilTTend", Mesh->NCellsSize, VCoord->NVertLayers);
   FrazilSTend =
       Array2DReal("FrazilSTend", Mesh->NCellsSize, VCoord->NVertLayers);
   FrazilHTend =
       Array2DReal("FrazilHTend", Mesh->NCellsSize, VCoord->NVertLayers);

   AccMIce  = Array1DReal("AccMIce", Mesh->NCellsSize);
   AccEIce  = Array1DReal("AccEIce", Mesh->NCellsSize);
   AccMLiq  = Array1DReal("AccMLiq", Mesh->NCellsSize);
   AccELiq  = Array1DReal("AccELiq", Mesh->NCellsSize);
   AccMSalt = Array1DReal("AccMSalt", Mesh->NCellsSize);

   deepCopy(FrazilTTend, 0.0_Real);
   deepCopy(FrazilSTend, 0.0_Real);
   deepCopy(FrazilHTend, 0.0_Real);
   deepCopy(AccMIce, 0.0_Real);
   deepCopy(AccEIce, 0.0_Real);
   deepCopy(AccMLiq, 0.0_Real);
   deepCopy(AccELiq, 0.0_Real);
   deepCopy(AccMSalt, 0.0_Real);
}

Frazil::~Frazil() {}

Frazil *Frazil::create(const std::string &Name) {
   if (AllFrazil.find(Name) != AllFrazil.end()) {
      LOG_ERROR("Attempted to create Frazil {} but it already exists", Name);
      return nullptr;
   }

   auto *NewFrazil =
       new Frazil(HorzMesh::getDefault(), VertCoord::getDefault());
   AllFrazil.emplace(Name, NewFrazil);

   Error Err;
   Config *OmegaConfig = Config::getOmegaConfig();
   Config FrazilConfig("Frazil");
   Err += OmegaConfig->get(FrazilConfig);
   CHECK_ERROR_ABORT(Err, "Frazil::create: Frazil group not found in Config");

   std::string FrazilTypeStr;
   Err += FrazilConfig.get("FrazilType", FrazilTypeStr);
   CHECK_ERROR_ABORT(Err,
                     "Frazil::create: FrazilType not found in Frazil config");

   if ((FrazilTypeStr == "Basic") or (FrazilTypeStr == "basic") or
       (FrazilTypeStr == "BasicFrazil")) {
      NewFrazil->frazilChoice = FrazilType::BasicFrazil;
   } else if ((FrazilTypeStr == "Simple") or (FrazilTypeStr == "simple") or
              (FrazilTypeStr == "SimpleFrazil")) {
      NewFrazil->frazilChoice = FrazilType::SimpleFrazil;
      ABORT_ERROR("Frazil::create: SimpleFrazil not supported yet");
   } else if ((FrazilTypeStr == "teos") or (FrazilTypeStr == "Teos") or
              (FrazilTypeStr == "TEOS") or (FrazilTypeStr == "Teos10") or
              (FrazilTypeStr == "teos10") or (FrazilTypeStr == "TEOS10")) {
      NewFrazil->frazilChoice = FrazilType::TeosFrazil;
   } else {
      ABORT_ERROR("Frazil::create: Unknown FrazilType requested");
   }

   Err += FrazilConfig.get("MassLimit",
                           NewFrazil->computeFrazilFormation.massLimit);
   Err += FrazilConfig.get("MassLimit", NewFrazil->computeFrazilMelt.massLimit);
   CHECK_ERROR_ABORT(Err,
                     "Frazil::create: MassLimit not found in Frazil config");

   NewFrazil->computeBasicFrazilFormation.FractionalThicknessLimit =
       NewFrazil->computeFrazilFormation.MassLimit;
   NewFrazil->computeBasicFrazilMelt.FractionalThicknessLimit =
       NewFrazil->computeFrazilFormation.MassLimit;

   Err += FrazilConfig.get("Phi", NewFrazil->computeFrazilFormation.Phi);
   CHECK_ERROR_ABORT(Err, "Frazil::create: Phi not found in Frazil config");

   Err += FrazilConfig.get("ConservationCheck", NewFrazil->conservationCheck);
   CHECK_ERROR_ABORT(
       Err, "Frazil::create: ConservationCheck not found in Frazil config");

   Err += FrazilConfig.get("DepthLimit", NewFrazil->depthLimit);
   CHECK_ERROR_ABORT(Err,
                     "Frazil::create: DepthLimit not found in Frazil config");

   if (Name == "Default") {
      DefaultFrazil = NewFrazil;
   }

   return NewFrazil;
}

Frazil *Frazil::getDefault() { return DefaultFrazil; }

Frazil *Frazil::get(const std::string &Name) {
   auto it = AllFrazil.find(Name);
   if (it != AllFrazil.end()) {
      return it->second.get();
   }

   LOG_ERROR("Frazil::get: Attempted to retrieve non-existent Frazil {}", Name);
   return nullptr;
}

void Frazil::erase(std::string InName) {
   auto *ToErase = get(InName);
   AllFrazil.erase(InName);
   if (ToErase == DefaultFrazil) {
      DefaultFrazil = nullptr;
   }
}

void Frazil::clear() {
   AllFrazil.clear();
   DefaultFrazil = nullptr;
}

void Frazil::checkColumnConservation() const {
   auto MinLayerCellH = createHostMirrorCopy(VCoordPtr->MinLayerCell);
   auto MaxLayerCellH = createHostMirrorCopy(VCoordPtr->MaxLayerCell);
   auto FrazilHTendH  = createHostMirrorCopy(FrazilHTend);
   auto FrazilTTendH  = createHostMirrorCopy(FrazilTTend);
   auto FrazilSTendH  = createHostMirrorCopy(FrazilSTend);
   auto AccMIceH      = createHostMirrorCopy(AccMIce);
   auto AccMLiqH      = createHostMirrorCopy(AccMLiq);
   auto AccMSaltH     = createHostMirrorCopy(AccMSalt);
   auto AccELiqH      = createHostMirrorCopy(AccELiq);
   auto AccEIceH      = createHostMirrorCopy(AccEIce);

   constexpr Real RTol = 1.0e-10_Real;

   for (I4 ICell = 0; ICell < NCellsAll; ++ICell) {
      const I4 KMin = MinLayerCellH(ICell);
      const I4 KMax = MaxLayerCellH(ICell);

      Real MassTend   = 0.0_Real;
      Real EnergyTend = 0.0_Real;
      Real SaltTend   = 0.0_Real;

      for (I4 K = KMin; K <= KMax; ++K) {
         MassTend += FrazilHTendH(ICell, K);
         EnergyTend += FrazilTTendH(ICell, K);
         SaltTend += FrazilSTendH(ICell, K);
      }

      const Real MassTotal   = AccMIceH(ICell) + AccMLiqH(ICell);
      const Real EnergyTotal = AccELiqH(ICell) + AccEIceH(ICell);
      const Real SaltTotal   = AccMSaltH(ICell);
      if (ICell == 0) {
         LOG_INFO("Frazil column conservation check: cell {} MassTend={} "
                  "MassTotal={} "
                  "EnergyTend={} EnergyTotal={} SaltTend={} SaltTotal={}",
                  ICell, MassTend * RhoSw, MassTotal,
                  EnergyTend * Cp0Sw * RhoSw, EnergyTotal,
                  SaltTend * RhoSw * PPt2Salt, SaltTotal);
         LOG_INFO("Frazil column conservation check: cell {} EpsMass={} "
                  "EpsE={} EpsS={} ",
                  ICell, MassTend * RhoSw + MassTotal,
                  EnergyTend * Cp0Sw * RhoSw + EnergyTotal,
                  SaltTend * RhoSw * PPt2Salt + SaltTotal);
      }

      if (!isApprox(-MassTend * RhoSw, MassTotal, RTol)) {
         LOG_INFO(
             "Frazil column mass check failed: cell {} tendency={} total={}",
             ICell, -MassTend * RhoSw, MassTotal);
      }
      if (!isApprox(-EnergyTend * Cp0Sw * RhoSw, EnergyTotal, RTol)) {
         LOG_INFO(
             "Frazil column energy check failed: cell {} tendency={} total={}",
             ICell, -EnergyTend * Cp0Sw * RhoSw, EnergyTotal);
      }
      if (!isApprox(-SaltTend * RhoSw * PPt2Salt, SaltTotal, RTol)) {
         LOG_INFO(
             "Frazil column salt check failed: cell {} tendency={} total={}",
             ICell, -SaltTend * RhoSw * PPt2Salt, SaltTotal);
      }
   }
}

void Frazil::computeFrazilBasicImpl(const Array2DReal &CT,
                                    const Array2DReal &SA, const Array2DReal &P,
                                    const Array2DReal &LayerH) {
   const EosType LocEosChoice = Eos::getInstance()->EosChoice;
   const Real LocDepthLimit   = depthLimit;

   OMEGA_SCOPE(MinLayerCell, VCoordPtr->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoordPtr->MaxLayerCell);
   OMEGA_SCOPE(LocGeomZMid, VCoordPtr->GeomZMid);

   OMEGA_SCOPE(LocComputeBasicFrazilFormation, computeBasicFrazilFormation);
   OMEGA_SCOPE(LocComputeBasicFrazilMelt, computeBasicFrazilMelt);
   OMEGA_SCOPE(LocFrazilTTend, FrazilTTend);
   OMEGA_SCOPE(LocFrazilSTend, FrazilSTend);
   OMEGA_SCOPE(LocFrazilHTend, FrazilHTend);
   OMEGA_SCOPE(LocAccMIce, AccMIce);
   OMEGA_SCOPE(LocAccEIce, AccEIce);
   OMEGA_SCOPE(LocAccMLiq, AccMLiq);
   OMEGA_SCOPE(LocAccELiq, AccELiq);
   OMEGA_SCOPE(LocAccMSalt, AccMSalt);
   OMEGA_SCOPE(LocIceRefSal, IceRefSal);
   OMEGA_SCOPE(LocLatIce, LatIce);

   parallelFor(
       {NCellsAll}, KOKKOS_LAMBDA(I4 ICell) {
          const I4 KMin = MinLayerCell(ICell);
          const I4 KMax = MaxLayerCell(ICell);

          I4 Klim          = KMax;
          bool HasKlim     = true;
          const bool Limit = (LocDepthLimit >= 0.0_Real);

          if (Limit) {
             HasKlim = false;
             for (I4 K = KMax; K >= KMin; --K) {
                if (Kokkos::abs(LocGeomZMid(ICell, K)) <= LocDepthLimit) {
                   Klim    = K;
                   HasKlim = true;
                   break;
                }
             }
          }

          // Explicit accumulation order: bottom layer to top layer.
          for (I4 K = KMax; K >= KMin; --K) {
             if (!HasKlim || K > Klim) {
                LocFrazilHTend(ICell, K) = 0.0_Real;
                LocFrazilTTend(ICell, K) = 0.0_Real;
                LocFrazilSTend(ICell, K) = 0.0_Real;
                continue;
             }

             const Real SAIn = SA(ICell, K);
             const Real CTIn = CT(ICell, K);
             const Real PIn  = P(ICell, K);
             const Real PDb  = PIn * Pa2Db;
             const Real H    = LayerH(ICell, K);

             const Real Tfrz =
                 Eos::calcCtFreezing(SAIn, PDb, 0.0_Real, LocEosChoice);

             Real HTend = 0.0_Real;
             Real TTend = 0.0_Real;
             Real STend = 0.0_Real;

             if (CTIn < Tfrz) {
                LocComputeBasicFrazilFormation(
                    SAIn, CTIn, PDb, H, LocAccMIce(ICell), LocAccMSalt(ICell),
                    LocAccEIce(ICell), HTend, TTend, STend, Tfrz);
             } else if (LocAccMIce(ICell) > 0.0_Real) {
                LocComputeBasicFrazilMelt(SAIn, CTIn, PDb, H, LocAccMIce(ICell),
                                          LocAccMSalt(ICell), LocAccEIce(ICell),
                                          HTend, TTend, STend, Tfrz);
             }

             // temporary log -- TBRemoved
             // if  (ICell == 0) {
             // LOG_INFO("computeFrazil cell = {}, SAIn = {}, CTIn = {}, PIn = "
             //          "{}, H = {}, Tfrz = {}",
             //          ICell, SAIn, CTIn, PIn, H, Tfrz);
             // LOG_INFO("computeFrazil cell={} K={} (cold={}) AccMIce={} "
             //          "AccMLiq={} AccMSalt={} AccELiq={} AccEIce={}",
             //          ICell, K, (CTIn < Tfrz), LocAccMIce(ICell),
             //          LocAccMLiq(ICell), LocAccMSalt(ICell),
             //          LocAccELiq(ICell), LocAccEIce(ICell));
             // LOG_INFO("                                     HTend= {} "
             //          "TTend= {} STend= {}",
             //          HTend, TTend, STend);
             //}

             LocFrazilHTend(ICell, K) = HTend; // not scaled by dt
             LocFrazilTTend(ICell, K) = TTend;
             LocFrazilSTend(ICell, K) = STend;
          } // end of vertical loop

          //  // temporary log -- TBRemoved
          //  if (ICell == 0) {
          //     LOG_INFO(
          //         "Frazil::computeFrazil: (basic={}), cell={} AccMIce={} "
          //         "AccMSaltCalc={} AccMSaltCpl={} "
          //         "AccEIceCalc={} AccEIceCpl={}",
          //         true, ICell, LocAccMIce(ICell),
          //         LocAccMSalt(ICell), LocAccMIce(ICell) * LocIceRefSal,
          //         LocAccEIce(ICell), -LocAccMIce(ICell) * LocLatIce);
          //  }

          // Redistribute excess salt at the surface. No treatment of low
          // salinity frazil for now.
          LocFrazilSTend(ICell, KMin) += Kokkos::max(
              0.0_Real, LocAccMIce(ICell) * LocIceRefSal - LocAccMSalt(ICell));
          // hijack total terms before the coupling
          LocAccMSalt(ICell) = LocAccMIce(ICell) * LocIceRefSal;
          LocAccEIce(ICell)  = -LocAccMIce(ICell) * LocLatIce;

          // Convert to coupler units
          LocAccMIce(ICell)  = LocAccMIce(ICell) * RhoSw;
          LocAccMLiq(ICell)  = LocAccMLiq(ICell) * RhoSw;
          LocAccMSalt(ICell) = LocAccMSalt(ICell) * RhoSw * PPt2Salt;
          LocAccELiq(ICell)  = LocAccELiq(ICell) * RhoSw;
          LocAccEIce(ICell)  = LocAccEIce(ICell) * RhoSw;
       }); // end of NCells loop
}

void Frazil::computeFrazilTeosImpl(const Array2DReal &CT, const Array2DReal &SA,
                                   const Array2DReal &P,
                                   const Array2DReal &LayerH) {
   const EosType LocEosChoice = Eos::getInstance()->EosChoice;
   const Real LocDepthLimit   = depthLimit;

   OMEGA_SCOPE(MinLayerCell, VCoordPtr->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoordPtr->MaxLayerCell);
   OMEGA_SCOPE(LocGeomZMid, VCoordPtr->GeomZMid);

   OMEGA_SCOPE(LocComputeFrazilFormation, computeFrazilFormation);
   OMEGA_SCOPE(LocComputeFrazilMelt, computeFrazilMelt);
   OMEGA_SCOPE(LocFrazilTTend, FrazilTTend);
   OMEGA_SCOPE(LocFrazilSTend, FrazilSTend);
   OMEGA_SCOPE(LocFrazilHTend, FrazilHTend);
   OMEGA_SCOPE(LocAccMIce, AccMIce);
   OMEGA_SCOPE(LocAccEIce, AccEIce);
   OMEGA_SCOPE(LocAccMLiq, AccMLiq);
   OMEGA_SCOPE(LocAccELiq, AccELiq);
   OMEGA_SCOPE(LocAccMSalt, AccMSalt);

   parallelFor(
       {NCellsAll}, KOKKOS_LAMBDA(I4 ICell) {
          const I4 KMin = MinLayerCell(ICell);
          const I4 KMax = MaxLayerCell(ICell);

          I4 Klim          = KMax;
          bool HasKlim     = true;
          const bool Limit = (LocDepthLimit >= 0.0_Real);

          if (Limit) {
             HasKlim = false;
             for (I4 K = KMax; K >= KMin; --K) {
                if (Kokkos::abs(LocGeomZMid(ICell, K)) <= LocDepthLimit) {
                   Klim    = K;
                   HasKlim = true;
                   break;
                }
             }
          }

          // Explicit accumulation order: bottom layer to top layer.
          for (I4 K = KMax; K >= KMin; --K) {
             if (!HasKlim || K > Klim) {
                LocFrazilHTend(ICell, K) = 0.0_Real;
                LocFrazilTTend(ICell, K) = 0.0_Real;
                LocFrazilSTend(ICell, K) = 0.0_Real;
                continue;
             }

             const Real SAIn = SA(ICell, K);
             const Real CTIn = CT(ICell, K);
             const Real PIn  = P(ICell, K);
             const Real PDb  = PIn * Pa2Db;
             const Real H    = LayerH(ICell, K);

             const Real Tfrz =
                 Eos::calcCtFreezing(SAIn, PDb, 0.0_Real, LocEosChoice);

             Real HTend = 0.0_Real;
             Real TTend = 0.0_Real;
             Real STend = 0.0_Real;

             if (CTIn < Tfrz) {
                LocComputeFrazilFormation(SAIn, CTIn, PDb, H, LocAccMIce(ICell),
                                          LocAccMLiq(ICell), LocAccMSalt(ICell),
                                          LocAccELiq(ICell), LocAccEIce(ICell),
                                          HTend, TTend, STend);
             } else if (LocAccMIce(ICell) > 0.0_Real) {
                LocComputeFrazilMelt(SAIn, CTIn, PDb, H, LocAccMIce(ICell),
                                     LocAccMLiq(ICell), LocAccMSalt(ICell),
                                     LocAccELiq(ICell), LocAccEIce(ICell),
                                     HTend, TTend, STend);
             }

             //  // temporary kernel logging -- TBRemoved
             //  if (ICell == 0) {
             //     const Real Hf  = H + HTend;
             //     const Real SAf = (H * SAIn + STend) / Hf;
             //     const Real Tf  = (H * CTIn + TTend) / Hf;
             //     // LOG_INFO("computeFrazil cell = {}, SAIn = {}, CTIn = {},
             //     PIn
             //     // = "
             //     //          "{}, H = {}, Tfrz = {}",
             //     //          ICell, SAIn, CTIn, PIn, H, Tfrz);
             //     LOG_INFO("computeFrazil cell={} K={} (cold={}) AccMIce={} "
             //              "AccMLiq={} AccMSalt={} AccELiq={} AccEIce={}",
             //              ICell, K, (CTIn < Tfrz), LocAccMIce(ICell),
             //              LocAccMLiq(ICell), LocAccMSalt(ICell),
             //              LocAccELiq(ICell), LocAccEIce(ICell));
             //     LOG_INFO("                                     HTend= {} "
             //              "TTend= {} STend= {}",
             //              HTend, TTend, STend);
             //     LOG_INFO("                                    "
             //              "H= {}-->{}, T= {}-->{}, S= {}-->{}",
             //              H, Hf, CTIn, Tf, SAIn, SAf);
             //  }

             LocFrazilHTend(ICell, K) = HTend; // not scaled by dt
             LocFrazilTTend(ICell, K) = TTend;
             LocFrazilSTend(ICell, K) = STend;
          } // end of vertical loop

          // Convert to coupler units
          LocAccMIce(ICell)  = LocAccMIce(ICell) * RhoSw;
          LocAccMLiq(ICell)  = LocAccMLiq(ICell) * RhoSw;
          LocAccMSalt(ICell) = LocAccMSalt(ICell) * RhoSw * PPt2Salt;
          LocAccELiq(ICell)  = LocAccELiq(ICell) * RhoSw;
          LocAccEIce(ICell)  = LocAccEIce(ICell) * RhoSw;
       }); // end of NCells loop
}

void Frazil::computeFrazil(const Array2DReal &CT, const Array2DReal &SA,
                           const Array2DReal &P, const Array2DReal &LayerH) {
   Eos *DefEos = Eos::getInstance();
   if (!DefEos) {
      ABORT_ERROR("Frazil::computeFrazil: Eos must be initialized before "
                  "computeFrazil");
   }

   const EosType LocEosChoice = DefEos->EosChoice;
   if (LocEosChoice != EosType::Teos10Eos) {
      ABORT_ERROR("Frazil::computeFrazil: CtFreezing not implemented for "
                  "non-TEOS-10 EOS");
   }

   switch (frazilChoice) {
   case FrazilType::BasicFrazil:
      computeFrazilBasicImpl(CT, SA, P, LayerH);
      break;
   case FrazilType::TeosFrazil:
      computeFrazilTeosImpl(CT, SA, P, LayerH);
      break;
   case FrazilType::SimpleFrazil:
      ABORT_ERROR("Frazil::computeFrazil: SimpleFrazil not supported yet");
      break;
   default:
      ABORT_ERROR("Frazil::computeFrazil: Unknown frazilChoice");
      break;
   }

   if (conservationCheck) {
      checkColumnConservation();
   }
} // end of computeFrazil

} // namespace OMEGA

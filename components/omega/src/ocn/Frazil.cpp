//===-- ocn/Frazil.cpp - Frazil Ice Formation -----------------*- C++ -*-===//
//
// The Frazil class manages frazil tendencies and accumulators.
// This initial scaffold wires allocation and configuration only.
//
//===----------------------------------------------------------------------===//

#include "Frazil.h"
#include "Error.h"
#include "Logging.h"

namespace OMEGA {

Frazil *Frazil::DefaultFrazil = nullptr;
std::map<std::string, std::unique_ptr<Frazil>> Frazil::AllFrazil;

/// Constructor for FrazilFormation
FrazilFormation::FrazilFormation() {}

/// Constructor for FrazilMelt
FrazilMelt::FrazilMelt() {}

void Frazil::init() {

   if (!HorzMesh::getDefault() or !VertCoord::getDefault()) {
      ABORT_ERROR("Frazil::init: HorzMesh and VertCoord must be initialized");
   }

   if (!DefaultFrazil) {
      DefaultFrazil = create("Default");
   }
}

Frazil::Frazil(const HorzMesh *Mesh, const VertCoord *VCoord)
    : frazilChoice(FrazilType::TeosFrazil), massLimit(0.1_Real), phi(0.75_Real),
      NCellsAll(Mesh->NCellsAll),
      NChunks((VCoord->NVertLayers + VecLength - 1) / VecLength), MeshPtr(Mesh),
      VCoordPtr(VCoord), computeFrazilFormation(), computeFrazilMelt() {

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
   } else if ((FrazilTypeStr == "teos") or (FrazilTypeStr == "Teos") or
              (FrazilTypeStr == "TEOS") or (FrazilTypeStr == "Teos10") or
              (FrazilTypeStr == "teos10") or (FrazilTypeStr == "TEOS10")) {
      NewFrazil->frazilChoice = FrazilType::TeosFrazil;
   } else {
      ABORT_ERROR("Frazil::create: Unknown FrazilType requested");
   }

   Err += FrazilConfig.get("massLimit",
                           NewFrazil->computeFrazilFormation.MassLimit);
   CHECK_ERROR_ABORT(Err,
                     "Frazil::create: massLimit not found in Frazil config");

   Err += FrazilConfig.get("phi", NewFrazil->computeFrazilFormation.Phi);
   CHECK_ERROR_ABORT(Err, "Frazil::create: phi not found in Frazil config");

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

void Frazil::computeFrazil(const Array2DReal &CT, const Array2DReal &SA,
                           const Array2DReal &P, const Array2DReal &LayerH) {
   OMEGA_SCOPE(MinLayerCell, VCoordPtr->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoordPtr->MaxLayerCell);

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

          // Explicit accumulation order: bottom layer to top layer.
          for (I4 K = KMax; K >= KMin; --K) {
             const Real SAIn = SA(ICell, K);
             const Real CTIn = CT(ICell, K);
             const Real PIn  = P(ICell, K);
             const Real H    = LayerH(ICell, K);

             const Real Tfrz = gsw_ct_freezing_poly(SAIn, PIn, 0.0_Real);

             Real HTend = 0.0_Real;
             Real TTend = 0.0_Real;
             Real STend = 0.0_Real;
             Real Dt    = 1800.0_Real; // hard-coded for now (30min in s)

             if (CTIn < Tfrz) {
                LocComputeFrazilFormation(
                    SAIn, CTIn, PIn, H, Dt, LocAccMIce(ICell),
                    LocAccMLiq(ICell), LocAccMSalt(ICell), LocAccELiq(ICell),
                    LocAccEIce(ICell), HTend, TTend, STend);
             } else {
                LocComputeFrazilMelt(LocAccMIce(ICell), LocAccMLiq(ICell),
                                     LocAccMSalt(ICell), LocAccELiq(ICell),
                                     LocAccEIce(ICell), HTend, TTend, STend);
             }

             // LocAccMIce(ICell) += solidMass;
             // LocAccMLiq(ICell) += liquidMass;
             // LocAccMSalt(ICell) += liquidMass * SAnew;
             // LocAccELiq(ICell) += liquidMass * Cp0Sw * CTnew;
             // LocAccEIce(ICell) +=
             //     solidMass * gsw_pot_enthalpy_from_pt_ice_poly(CTnew);

             LocFrazilHTend(ICell, K) = HTend;
             LocFrazilTTend(ICell, K) = TTend;
             LocFrazilSTend(ICell, K) = STend;
          }
       });
}

} // namespace OMEGA

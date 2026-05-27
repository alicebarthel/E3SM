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
      NChunks((VCoord->NVertLayers + VecLength - 1) / VecLength) {

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

   Err += FrazilConfig.get("massLimit", NewFrazil->massLimit);
   CHECK_ERROR_ABORT(Err,
                     "Frazil::create: massLimit not found in Frazil config");

   Err += FrazilConfig.get("phi", NewFrazil->phi);
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

void Frazil::checkFrazil() {
   // Placeholder for frazil consistency checks.
}

void Frazil::computeFrazilFormation() {
   // Placeholder for frazil formation implementation.
}

void Frazil::computeFrazilMelt() {
   // Placeholder for frazil melt implementation.
}

} // namespace OMEGA

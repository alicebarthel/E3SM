#include "Forcing.h"
#include "Field.h"
#include "Logging.h"
#include "Pacer.h"

namespace OMEGA {

Forcing *Forcing::DefaultForcing = nullptr;
std::map<std::string, std::unique_ptr<Forcing>> Forcing::AllForcing;

static std::string stripDefault(const std::string &Name) {
   return Name != "Default" ? Name : "";
}

Forcing::Forcing(const std::string &Name, const HorzMesh *Mesh, Halo *MeshHalo)
    : Name(stripDefault(Name)), MomForcingAux(stripDefault(Name), Mesh),
      Mesh(Mesh), MeshHalo(MeshHalo) {}

Forcing::~Forcing() { unregisterFields(); }

void Forcing::registerFields(const std::string &MeshName) const {
   MomForcingAux.registerFields(MeshName);
}

void Forcing::unregisterFields() const { MomForcingAux.unregisterFields(); }

Forcing *Forcing::create(const std::string &Name, const HorzMesh *Mesh,
                         Halo *MeshHalo) {
   if (AllForcing.find(Name) != AllForcing.end()) {
      LOG_ERROR("Attempted to create new Forcing with name {} but it already "
                "exists",
                Name);
      return nullptr;
   }

   auto *NewForcing = new Forcing(Name, Mesh, MeshHalo);
   AllForcing.emplace(Name, NewForcing);

   return NewForcing;
}

void Forcing::init() {
   if (DefaultForcing != nullptr) {
      return;
   }

   FieldGroup::create("Forcing");

   const HorzMesh *DefMesh = HorzMesh::getDefault();
   Halo *DefHalo           = Halo::getDefault();

   DefaultForcing = Forcing::create("Default", DefMesh, DefHalo);

   if (DefaultForcing == nullptr) {
      ABORT_ERROR("Forcing: failed to initialize default forcing state");
   }

   DefaultForcing->registerFields(DefMesh->MeshName);

   Config *OmegaConfig = Config::getOmegaConfig();
   DefaultForcing->readConfigOptions(OmegaConfig);
}

Forcing *Forcing::getDefault() { return DefaultForcing; }

Forcing *Forcing::get(const std::string &Name) {
   auto it = AllForcing.find(Name);
   if (it != AllForcing.end()) {
      return it->second.get();
   }

   LOG_ERROR("Forcing::get: Attempt to retrieve non-existent forcing state:");
   LOG_ERROR("{} has not been defined or has been removed", Name);
   return nullptr;
}

bool Forcing::exists(const std::string &Name) {
   return AllForcing.find(Name) != AllForcing.end();
}

void Forcing::erase(const std::string &Name) { AllForcing.erase(Name); }

void Forcing::clear() {
   AllForcing.clear();
   DefaultForcing = nullptr;
   if (FieldGroup::exists("Forcing")) {
      FieldGroup::destroy("Forcing");
   }
}

void Forcing::readConfigOptions(Config *OmegaConfig) {
   Error Err;

   Config SrfStressConfig("SrfStress");
   Err += OmegaConfig->get(SrfStressConfig);

   std::string SrfStressInterpTypeStr;
   Err += SrfStressConfig.get("InterpType", SrfStressInterpTypeStr);
   CHECK_ERROR_ABORT(Err, "Forcing: InterpType not found in SrfStressConfig");

   if (SrfStressInterpTypeStr == "Isotropic") {
      this->MomForcingAux.InterpChoice = InterpCellToEdgeOption::Isotropic;
   } else if (SrfStressInterpTypeStr == "Anisotropic") {
      this->MomForcingAux.InterpChoice = InterpCellToEdgeOption::Anisotropic;
   } else {
      ABORT_ERROR("Forcing: Unknown InterpType requested");
   }
}

void Forcing::computeSrfStressForcingOnEdge() const {
   OMEGA_SCOPE(LocMomForcingAux, MomForcingAux);

   Pacer::start("Forcing:edgeAuxState1", 2);
   parallelFor(
       "Forcing:edgeAuxState1", {Mesh->NEdgesAll},
       KOKKOS_LAMBDA(int IEdge) { LocMomForcingAux.computeVarsOnEdge(IEdge); });
   Pacer::stop("Forcing:edgeAuxState1", 2);
}

I4 Forcing::exchangeHalo() const {
   I4 Err = 0;

   Err +=
       MeshHalo->exchangeFullArrayHalo(MomForcingAux.ZonalStressCell, OnCell);
   Err +=
       MeshHalo->exchangeFullArrayHalo(MomForcingAux.MeridStressCell, OnCell);

   return Err;
}

} // namespace OMEGA

#include "Forcing.h"
#include "Eos.h"
#include "Logging.h"
#include "Pacer.h"
#include "Tracers.h"

namespace OMEGA {

Forcing *Forcing::DefaultForcing = nullptr;
std::map<std::string, std::unique_ptr<Forcing>> Forcing::AllForcing;

static std::string stripDefault(const std::string &Name) {
   return Name != "Default" ? Name : "";
}

Forcing::Forcing(const std::string &Name, const HorzMesh *Mesh, Halo *MeshHalo,
                 const VertCoord *VCoord, int NTracers)
    : Name(stripDefault(Name)), CplForcingAux(stripDefault(Name), Mesh),
      SurfTracerRestAux(stripDefault(Name), Mesh, NTracers),
      WindForcingAux(stripDefault(Name), Mesh), Mesh(Mesh), MeshHalo(MeshHalo),
      VCoord(VCoord) {}

Forcing::~Forcing() { unregisterFields(); }

void Forcing::registerFields(const std::string &AuxGroupName,
                             const std::string &MeshName) const {
   WindForcingAux.registerFields(AuxGroupName, MeshName);
   CplForcingAux.registerFields(AuxGroupName, MeshName);
   SurfTracerRestAux.registerFields(AuxGroupName, MeshName);
}

void Forcing::unregisterFields() const {
   WindForcingAux.unregisterFields();
   CplForcingAux.unregisterFields();
   SurfTracerRestAux.unregisterFields();
}

Forcing *Forcing::create(const std::string &Name, const HorzMesh *Mesh,
                         Halo *MeshHalo, const VertCoord *VCoord,
                         int NTracers) {
   if (AllForcing.find(Name) != AllForcing.end()) {
      LOG_ERROR("Attempted to create new Forcing with name {} but it "
                "already exists",
                Name);
      return nullptr;
   }

   auto *NewForcing = new Forcing(Name, Mesh, MeshHalo, VCoord, NTracers);
   AllForcing.emplace(Name, NewForcing);

   return NewForcing;
}

void Forcing::init() {
   if (DefaultForcing != nullptr) {
      return;
   }

   const HorzMesh *DefMesh    = HorzMesh::getDefault();
   Halo *DefHalo              = Halo::getDefault();
   const VertCoord *DefVCoord = VertCoord::getDefault();
   int NTracers               = Tracers::getNumTracers();

   DefaultForcing =
       Forcing::create("Default", DefMesh, DefHalo, DefVCoord, NTracers);

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
}

void Forcing::readConfigOptions(Config *OmegaConfig) {
   Error Err;

   Config WindStressConfig("WindStress");
   Err += OmegaConfig->get(WindStressConfig);

   std::string WindStressInterpTypeStr;
   Err += WindStressConfig.get("InterpType", WindStressInterpTypeStr);
   CHECK_ERROR_ABORT(Err, "Forcing: InterpType not found in WindStressConfig");

   if (WindStressInterpTypeStr == "Isotropic") {
      this->WindForcingAux.InterpChoice = InterpCellToEdgeOption::Isotropic;
   } else if (WindStressInterpTypeStr == "Anisotropic") {
      this->WindForcingAux.InterpChoice = InterpCellToEdgeOption::Anisotropic;
   } else {
      ABORT_ERROR("Forcing: Unknown InterpType requested");
   }
}

void Forcing::computeWindForcingOnEdge() const {
   OMEGA_SCOPE(LocWindForcingAux, WindForcingAux);

   Pacer::start("Forcing:edgeAuxState1", 2);
   parallelFor(
       "Forcing:edgeAuxState1", {Mesh->NEdgesAll}, KOKKOS_LAMBDA(int IEdge) {
          LocWindForcingAux.computeVarsOnEdge(IEdge);
       });
   Pacer::stop("Forcing:edgeAuxState1", 2);
}

void Forcing::computeSurfInsituTemp(const Array3DReal &TracerArray) const {
   CplForcingAux.computeSurfInsituTemp(TracerArray, VCoord, Eos::getInstance());
}

I4 Forcing::exchangeHalo() const {
   I4 Err = 0;

   Err +=
       MeshHalo->exchangeFullArrayHalo(WindForcingAux.ZonalStressCell, OnCell);
   Err +=
       MeshHalo->exchangeFullArrayHalo(WindForcingAux.MeridStressCell, OnCell);

   const I4 NTracers =
       SurfTracerRestAux.TracersMonthlySurfClimoCell.extent_int(0);
   for (I4 LTracer = 0; LTracer < NTracers; ++LTracer) {
      auto TracerSurfClimoCell = Kokkos::subview(
          SurfTracerRestAux.TracersMonthlySurfClimoCell, LTracer, Kokkos::ALL);
      Err += MeshHalo->exchangeFullArrayHalo(TracerSurfClimoCell, OnCell);
   }

   Err += MeshHalo->exchangeFullArrayHalo(CplForcingAux.SnowFluxCell, OnCell);
   Err += MeshHalo->exchangeFullArrayHalo(CplForcingAux.RainFluxCell, OnCell);
   Err += MeshHalo->exchangeFullArrayHalo(CplForcingAux.EvaporationFluxCell,
                                          OnCell);
   Err += MeshHalo->exchangeFullArrayHalo(
       CplForcingAux.SeaIceFreshWaterFluxCell, OnCell);
   Err +=
       MeshHalo->exchangeFullArrayHalo(CplForcingAux.IceRunoffFluxCell, OnCell);
   Err += MeshHalo->exchangeFullArrayHalo(CplForcingAux.RiverRunoffFluxCell,
                                          OnCell);
   Err += MeshHalo->exchangeFullArrayHalo(CplForcingAux.LatentHeatFluxCell,
                                          OnCell);
   Err += MeshHalo->exchangeFullArrayHalo(CplForcingAux.SensibleHeatFluxCell,
                                          OnCell);
   Err += MeshHalo->exchangeFullArrayHalo(CplForcingAux.LongWaveHeatFluxUpCell,
                                          OnCell);
   Err += MeshHalo->exchangeFullArrayHalo(
       CplForcingAux.LongWaveHeatFluxDownCell, OnCell);
   Err += MeshHalo->exchangeFullArrayHalo(CplForcingAux.SeaIceHeatFluxCell,
                                          OnCell);
   Err += MeshHalo->exchangeFullArrayHalo(CplForcingAux.ShortWaveHeatFluxCell,
                                          OnCell);
   Err += MeshHalo->exchangeFullArrayHalo(CplForcingAux.SeaIceSaltFluxCell,
                                          OnCell);

   return Err;
}

} // namespace OMEGA

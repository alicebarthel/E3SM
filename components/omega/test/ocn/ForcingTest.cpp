//===-- ocn/ForcingTest.cpp - Forcing Unit Test ---------------*- C++ -*-===//
//
/// \file
/// \brief Scoped tests for the Forcing class (SfcStress + SfcStressForcingAux)
//
//===----------------------------------------------------------------------===//

#include "Forcing.h"
#include "Config.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Dimension.h"
#include "Error.h"
#include "Field.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "IOStream.h"
#include "Logging.h"
#include "MachEnv.h"
#include "Pacer.h"
#include "TimeStepper.h"
#include "mpi.h"

#include <limits>
#include <string>

using namespace OMEGA;

namespace {

const std::string DefaultMeshFile = "OmegaMesh.nc";

int initForcingTest(const std::string &MeshFile) {
   int Err = 0;

   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   initLogging(DefEnv);

   Config("Omega");
   Config::readAll("omega.yml");

   TimeStepper::init1();
   TimeStepper *DefStepper = TimeStepper::getDefault();
   Clock *ModelClock       = DefStepper->getClock();

   IO::init(DefComm);
   Decomp::init(MeshFile);

   Field::init(ModelClock);
   IOStream::init(ModelClock);

   Err = Halo::init();
   if (Err != 0) {
      ABORT_ERROR("ForcingTest: error initializing default halo");
   }

   HorzMesh::init(ModelClock);
   Forcing::init();

   return 0;
}

void finalizeForcingTest() {
   Forcing::clear();
   HorzMesh::clear();
   Field::clear();
   Dimension::clear();
   IOStream::finalize();
   TimeStepper::clear();
   Halo::clear();
   Decomp::clear();
   MachEnv::removeAll();
}

int testForcingInitAndConfig() {
   int Err = 0;

   Forcing *DefForcing = Forcing::getDefault();
   if (DefForcing == nullptr) {
      LOG_ERROR("ForcingTest: default forcing instance is null");
      return 1;
   }

   Config *OmegaConfig = Config::getOmegaConfig();
   if (OmegaConfig == nullptr) {
      LOG_ERROR("ForcingTest: Omega config unavailable");
      return 1;
   }

   Error CfgErr;
   Config SfcStressConfig("SfcStress");
   CfgErr = OmegaConfig->get(SfcStressConfig);
   CHECK_ERROR_ABORT(CfgErr, "ForcingTest: missing Omega.SfcStress config");

   std::string InterpType;
   CfgErr = SfcStressConfig.get("InterpType", InterpType);
   CHECK_ERROR_ABORT(CfgErr, "ForcingTest: missing SfcStress.InterpType");

   InterpCellToEdgeOption ExpectedChoice;
   if (InterpType == "Isotropic") {
      ExpectedChoice = InterpCellToEdgeOption::Isotropic;
   } else if (InterpType == "Anisotropic") {
      ExpectedChoice = InterpCellToEdgeOption::Anisotropic;
   } else {
      LOG_ERROR("ForcingTest: unknown InterpType in config: {}", InterpType);
      return 1;
   }

   if (DefForcing->SfcStressForcingAux.InterpChoice != ExpectedChoice) {
      LOG_ERROR("ForcingTest: InterpChoice mismatch after Forcing::init");
      Err++;
   }

   if (Err == 0) {
      LOG_INFO("ForcingTest: Init/config PASS");
   }

   return Err;
}

int testForcingComputeAll() {
   int Err = 0;

   const HorzMesh *Mesh = HorzMesh::getDefault();
   Forcing *DefForcing  = Forcing::getDefault();
   if (Mesh == nullptr || DefForcing == nullptr) {
      LOG_ERROR("ForcingTest: missing mesh or forcing for compute test");
      return 1;
   }

   auto &ZonalStressCell = DefForcing->SfcStressForcingAux.ZonalStressCell;
   auto &MeridStressCell = DefForcing->SfcStressForcingAux.MeridStressCell;
   auto &NormalStress    = DefForcing->SfcStressForcingAux.NormalStressEdge;
   const auto AngleEdge  = Mesh->AngleEdge;

   // First pass: pure zonal stress should project to cos(AngleEdge).
   deepCopy(ZonalStressCell, 1._Real);
   deepCopy(MeridStressCell, 0._Real);
   deepCopy(NormalStress, 0._Real);

   DefForcing->computeAll();

   Real MaxErrCos = 0._Real;
   parallelReduce(
       {Mesh->NEdgesOwned},
       KOKKOS_LAMBDA(int IEdge, Real &LocalMax) {
          const Real Expected = Kokkos::cos(AngleEdge(IEdge));
          const Real AbsErr   = Kokkos::abs(NormalStress(IEdge) - Expected);
          if (AbsErr > LocalMax) {
             LocalMax = AbsErr;
          }
       },
       Kokkos::Max<Real>(MaxErrCos));

   Real GlobalMaxErrCos = 0._Real;
   MPI_Datatype MpiReal =
       sizeof(Real) == sizeof(double) ? MPI_DOUBLE : MPI_FLOAT;
   MPI_Allreduce(&MaxErrCos, &GlobalMaxErrCos, 1, MpiReal, MPI_MAX,
                 MachEnv::getDefault()->getComm());

   // Second pass: pure meridional stress should project to sin(AngleEdge).
   deepCopy(ZonalStressCell, 0._Real);
   deepCopy(MeridStressCell, 1._Real);
   deepCopy(NormalStress, 0._Real);

   DefForcing->computeAll();

   Real MaxErrSin = 0._Real;
   parallelReduce(
       {Mesh->NEdgesOwned},
       KOKKOS_LAMBDA(int IEdge, Real &LocalMax) {
          const Real Expected = Kokkos::sin(AngleEdge(IEdge));
          const Real AbsErr   = Kokkos::abs(NormalStress(IEdge) - Expected);
          if (AbsErr > LocalMax) {
             LocalMax = AbsErr;
          }
       },
       Kokkos::Max<Real>(MaxErrSin));

   Real GlobalMaxErrSin = 0._Real;
   MPI_Allreduce(&MaxErrSin, &GlobalMaxErrSin, 1, MpiReal, MPI_MAX,
                 MachEnv::getDefault()->getComm());

   const Real Tol       = 1e-11;
   const Real GlobalErr = Kokkos::max(GlobalMaxErrCos, GlobalMaxErrSin);
   if (GlobalErr > Tol) {
      LOG_ERROR(
          "ForcingTest: normal stress mismatch in cos/sin checks, max error "
          "{} > tol {}",
          GlobalErr, Tol);
      Err++;
   }

   if (Err == 0) {
      LOG_INFO("ForcingTest: computeAll PASS");
   }

   return Err;
}

int testForcingAPISmoke() {
   int Err = 0;

   const std::string Name = "UnitTestForcing";

   if (Forcing::exists(Name)) {
      Forcing::erase(Name);
   }

   Forcing *Named =
       Forcing::create(Name, HorzMesh::getDefault(), Halo::getDefault());
   if (Named == nullptr) {
      LOG_ERROR("ForcingTest: failed creating named forcing instance");
      Err++;
   }

   if (!Forcing::exists(Name)) {
      LOG_ERROR("ForcingTest: created forcing instance not found in map");
      Err++;
   }

   if (Forcing::get(Name) == nullptr) {
      LOG_ERROR("ForcingTest: get() failed for named forcing instance");
      Err++;
   }

   Forcing::erase(Name);
   if (Forcing::exists(Name)) {
      LOG_ERROR("ForcingTest: erase() failed for named forcing instance");
      Err++;
   }

   if (Err == 0) {
      LOG_INFO("ForcingTest: API smoke PASS");
   }

   return Err;
}

int forcingTest() {
   int Err = 0;

   Err += initForcingTest(DefaultMeshFile);
   Err += testForcingInitAndConfig();
   Err += testForcingComputeAll();
   Err += testForcingAPISmoke();

   if (Err == 0) {
      LOG_INFO("ForcingTest: Successful completion");
   }

   finalizeForcingTest();
   return Err;
}

} // namespace

int main(int argc, char *argv[]) {
   int RetErr = 0;

   MPI_Init(&argc, &argv);
   Kokkos::initialize(argc, argv);
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");

   RetErr = forcingTest();

   Pacer::finalize();
   Kokkos::finalize();
   MPI_Finalize();

   return RetErr;
}

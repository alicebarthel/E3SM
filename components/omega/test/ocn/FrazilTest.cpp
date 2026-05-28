//===-- Test driver for OMEGA FrazilFormation ---------------------------*- C++
//-*-===//
//
/// \file
/// \brief Minimal test driver for OMEGA frazil formation functor
//
//===-----------------------------------------------------------------------===/

#include "Frazil.h"
#include "Config.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Dimension.h"
#include "Field.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "IOStream.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanTestCommon.h"
#include "OmegaKokkos.h"
#include "Pacer.h"
#include "VertCoord.h"
#include "mpi.h"

using namespace OMEGA;

constexpr int NVertLayers = 60;

void initFrazilTest(const std::string &mesh) {
   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   initLogging(DefEnv);
   LOG_INFO("------ Frazil Unit Tests ------");

   Config("Omega");
   Config::readAll("omega.yml");

   IO::init(DefComm);
   IOStream::init();
   Decomp::init(mesh);
   Halo::init();
   HorzMesh::init();
   VertCoord::init(false);
   Frazil::init();
}

void finalizeFrazilTest() {
   Frazil::clear();
   VertCoord::clear();
   HorzMesh::clear();
   Halo::clear();
   Decomp::clear();
   Field::clear();
   Dimension::clear();
   MachEnv::removeAll();
}

void testFrazilFormationCold() {
   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();

   VCoord->NVertLayers = NVertLayers;

   const Real SAIn = 35.0_Real;
   const Real CTIn = -2.0_Real;
   const Real PIn  = 100.0_Real;
   const Real h    = 10.0_Real;
   const Real Dt   = 1800.0_Real;
   const Real RTol = 1e-10_Real;

   (void)Mesh;

   FrazilFormation ComputeFrazilFormation;

   Real AccMIce  = 0.0_Real;
   Real AccMLiq  = 0.0_Real;
   Real AccMSalt = 0.0_Real;
   Real AccELiq  = 0.0_Real;
   Real AccEIce  = 0.0_Real;

   Real HTend = 0.0_Real;
   Real TTend = 0.0_Real;
   Real STend = 0.0_Real;

   ComputeFrazilFormation(SAIn, CTIn, PIn, h, Dt, AccMIce, AccMLiq, AccMSalt,
                          AccELiq, AccEIce, HTend, TTend, STend);

   if (AccMIce <= 0.0_Real) {
      ABORT_ERROR("FrazilTestCold: accumulated ice mass is non-positive: {}",
                  AccMIce);
   }
   if (AccMLiq <= 0.0_Real) {
      ABORT_ERROR("FrazilTestCold: accumulated liquid mass is non-positive: {}",
                  AccMLiq);
   }
   if (AccMSalt <= 0.0_Real) {
      ABORT_ERROR("FrazilTestCold: accumulated salt mass is non-positive: {}",
                  AccMSalt);
   }

   if (isApprox(AccELiq, 0.0_Real, RTol)) {
      ABORT_ERROR(
          "FrazilTestCold: accumulated liquid energy is effectively zero: {}",
          AccELiq);
   }
   if (isApprox(AccEIce, 0.0_Real, RTol)) {
      ABORT_ERROR(
          "FrazilTestCold: accumulated ice energy is effectively zero: {}",
          AccEIce);
   }
   if (isApprox(HTend, 0.0_Real, RTol)) {
      ABORT_ERROR("FrazilTestCold: HTend is effectively zero: {}", HTend);
   }

   if (isApprox(TTend, 0.0_Real, RTol)) {
      ABORT_ERROR("FrazilTestCold: TTend is zero: {}", TTend);
   }

   if (isApprox(STend, 0.0_Real, RTol)) {
      ABORT_ERROR("FrazilTestCold: STend is effectively zero: {}", STend);
   }
   LOG_INFO("FrazilTestCold: AccMIce = {}, AccMLiq = {}, AccMSalt = {}, "
            "AccELiq = {}, AccEIce = {}, HTend = {}, TTend = {}, STend = {}",
            AccMIce, AccMLiq, AccMSalt, AccELiq, AccEIce, HTend, TTend, STend);
}

void testFrazilFormationWarm() {
   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();

   VCoord->NVertLayers = NVertLayers;

   const Real SAIn = 35.0_Real;
   const Real CTIn = 10.0_Real;
   const Real PIn  = 100.0_Real;
   const Real h    = 10.0_Real;
   const Real Dt   = 1800.0_Real;
   const Real RTol = 1e-10_Real;

   (void)Mesh;

   FrazilFormation ComputeFrazilFormation;

   Real AccMIce  = 0.0_Real;
   Real AccMLiq  = 0.0_Real;
   Real AccMSalt = 0.0_Real;
   Real AccELiq  = 0.0_Real;
   Real AccEIce  = 0.0_Real;

   Real HTend = 0.0_Real;
   Real TTend = 0.0_Real;
   Real STend = 0.0_Real;

   ComputeFrazilFormation(SAIn, CTIn, PIn, h, Dt, AccMIce, AccMLiq, AccMSalt,
                          AccELiq, AccEIce, HTend, TTend, STend);

   if (!isApprox(AccMIce, 0.0_Real, RTol)) {
      ABORT_ERROR("FrazilTest warm: expected zero AccMIce, got {}", AccMIce);
   }

   if (!isApprox(AccMSalt, 0.0_Real, RTol)) {
      ABORT_ERROR("FrazilTest warm: expected zero AccMSalt, got {}", AccMSalt);
   }
   if (!isApprox(AccMLiq, 0.0_Real, RTol)) {
      ABORT_ERROR("FrazilTest warm: expected zero AccMLiq, got {}", AccMLiq);
   }

   if (!isApprox(AccELiq, 0.0_Real, RTol)) {
      ABORT_ERROR("FrazilTest warm: expected zero AccELiq, got {}", AccELiq);
   }
   if (!isApprox(AccEIce, 0.0_Real, RTol)) {
      ABORT_ERROR("FrazilTest warm: expected zero AccEIce, got {}", AccEIce);
   }
   if (!isApprox(HTend, 0.0_Real, RTol)) {
      ABORT_ERROR("FrazilTest warm: expected zero HTend, got {}", HTend);
   }

   if (!isApprox(TTend, 0.0_Real, RTol)) {
      ABORT_ERROR("FrazilTest warm: expected zero TTend, got {}", TTend);
   }

   if (!isApprox(STend, 0.0_Real, RTol)) {
      ABORT_ERROR("FrazilTest warm: expected zero STend, got {}", STend);
   }
   LOG_INFO("FrazilTestWarm: AccMIce = {}, AccMLiq = {}, AccMSalt = {}, "
            "AccELiq = {}, AccEIce = {}, HTend = {}, TTend = {}, STend = {}",
            AccMIce, AccMLiq, AccMSalt, AccELiq, AccEIce, HTend, TTend, STend);
}

void frazilTest(const std::string &MeshFile = "OmegaMesh.nc") {
   initFrazilTest(MeshFile);
   testFrazilFormationCold();
   testFrazilFormationWarm();
   finalizeFrazilTest();
}

int main(int argc, char *argv[]) {
   MPI_Init(&argc, &argv);
   Kokkos::initialize(argc, argv);
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");

   frazilTest();

   LOG_INFO("------ Frazil Unit Tests Successful ------");

   Pacer::finalize();
   Kokkos::finalize();
   MPI_Finalize();

   return 0;
}

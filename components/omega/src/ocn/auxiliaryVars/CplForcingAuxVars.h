#ifndef OMEGA_AUX_CPL_FORCING_H
#define OMEGA_AUX_CPL_FORCING_H

#include "DataTypes.h"
#include "HorzMesh.h"

#include <string>

namespace OMEGA {

class CplForcingAuxVars {
 public:
   Array1DReal SnowFluxCell;
   Array1DReal RainFluxCell;
   Array1DReal EvaporationFluxCell;
   Array1DReal SeaIceFreshWaterFluxCell;
   Array1DReal IceRunoffFluxCell;
   Array1DReal RiverRunoffFluxCell;

   Array1DReal LatentHeatFluxCell;
   Array1DReal SensibleHeatFluxCell;
   Array1DReal LongWaveHeatFluxUpCell;
   Array1DReal LongWaveHeatFluxDownCell;
   Array1DReal SeaIceHeatFluxCell;
   Array1DReal ShortWaveHeatFluxCell;

   Array1DReal SeaIceSaltFluxCell;

   CplForcingAuxVars(const std::string &AuxStateSuffix, const HorzMesh *Mesh);

   void registerFields(const std::string &AuxGroupName,
                       const std::string &MeshName) const;
   void unregisterFields() const;

 private:
   std::string ForcingGroupName;
};

} // namespace OMEGA

#endif

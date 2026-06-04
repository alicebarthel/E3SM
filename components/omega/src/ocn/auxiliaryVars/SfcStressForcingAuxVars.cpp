#include "SfcStressForcingAuxVars.h"
#include "DataTypes.h"
#include "Field.h"

#include <limits>

namespace OMEGA {

SfcStressForcingAuxVars::SfcStressForcingAuxVars(
    const std::string &AuxStateSuffix, const HorzMesh *Mesh)
    : NormalStressEdge("NormalStressEdge" + AuxStateSuffix, Mesh->NEdgesSize),
      ZonalStressCell("SfcStressZonal" + AuxStateSuffix, Mesh->NCellsSize),
      MeridStressCell("SfcStressMeridional" + AuxStateSuffix, Mesh->NCellsSize),
      CellsOnEdge(Mesh->CellsOnEdge), AngleEdge(Mesh->AngleEdge), Interp(Mesh) {
}

void SfcStressForcingAuxVars::registerFields(
    const std::string &MeshName // name of horizontal mesh
) const {

   const Real FillValue = -9.99e30;
   int NDims            = 1;
   std::vector<std::string> DimNames(NDims);
   std::string DimSuffix;
   if (MeshName == "Default") {
      DimSuffix = "";
   } else {
      DimSuffix = MeshName;
   }

   DimNames[0] = "NCells" + DimSuffix;
   auto ZonalStressCellField =
       Field::create(ZonalStressCell.label(),          // field name
                     "zonal wind stress",              // long name/describe
                     "N m^{-2}",                       // units
                     "",                               // CF standard Name
                     std::numeric_limits<Real>::min(), // min valid value
                     std::numeric_limits<Real>::max(), // max valid value
                     FillValue,                        // scalar for undefined
                     NDims,                            // number of dimensions
                     DimNames                          // dim names
       );

   auto MeridStressCellField =
       Field::create(MeridStressCell.label(),  // field name
                     "meridional wind stress", // long Name or description
                     "N m^{-2}",               // units
                     "",                       // CF standard Name
                     std::numeric_limits<Real>::min(), // min valid value
                     std::numeric_limits<Real>::max(), // max valid value
                     FillValue,                        // scalar used undefined
                     NDims,                            // number of dimensions
                     DimNames                          // dimension names
       );

   FieldGroup::addFieldToGroup(ZonalStressCell.label(), "Forcing");
   FieldGroup::addFieldToGroup(MeridStressCell.label(), "Forcing");

   ZonalStressCellField->attachData<Array1DReal>(ZonalStressCell);
   MeridStressCellField->attachData<Array1DReal>(MeridStressCell);
}

void SfcStressForcingAuxVars::unregisterFields() const {
   Field::destroy(ZonalStressCell.label());
   Field::destroy(MeridStressCell.label());
}

} // namespace OMEGA

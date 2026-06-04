#ifndef OMEGA_FORCING_H
#define OMEGA_FORCING_H

#include "Config.h"
#include "DataTypes.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "auxiliaryVars/SfcStressForcingAuxVars.h"

#include <map>
#include <memory>
#include <string>

namespace OMEGA {

class Forcing {
 public:
   std::string Name;

   SfcStressForcingAuxVars SfcStressForcingAux;

   ~Forcing();

   static void init();

   static Forcing *create(const std::string &Name, const HorzMesh *Mesh,
                          Halo *MeshHalo);

   static Forcing *getDefault();

   static Forcing *get(const std::string &Name);

   static bool exists(const std::string &Name);

   static void erase(const std::string &Name);

   static void clear();

   void readConfigOptions(Config *OmegaConfig);

   void registerFields(const std::string &MeshName) const;

   void unregisterFields() const;

   void computeAll() const;

   void computeSfcStressForcingOnEdge() const;

   I4 exchangeHalo() const;

 private:
   Forcing(const std::string &Name, const HorzMesh *Mesh, Halo *MeshHalo);

   Forcing(const Forcing &) = delete;
   Forcing(Forcing &&)      = delete;

   const HorzMesh *Mesh;
   Halo *MeshHalo;

   static Forcing *DefaultForcing;
   static std::map<std::string, std::unique_ptr<Forcing>> AllForcing;
};

} // namespace OMEGA

#endif

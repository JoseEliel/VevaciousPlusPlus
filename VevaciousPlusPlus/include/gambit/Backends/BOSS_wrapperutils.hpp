#ifndef __BOSS_wrapperutils_vevacious_1_0_hpp__
#define __BOSS_wrapperutils_vevacious_1_0_hpp__

#include "backend_types/vevacious_1_0/wrapper_VevaciousPlusPlus.hpp"
#include "gambit/Backends/abstracttypedefs.hpp"
#include "gambit/Backends/wrappertypedefs.hpp"

VevaciousPlusPlus::Wrapper_VevaciousPlusPlus* wrapper_creator(VevaciousPlusPlus::Abstract_VevaciousPlusPlus*);

void wrapper_deleter(VevaciousPlusPlus::Wrapper_VevaciousPlusPlus*);

void set_delete_BEptr(VevaciousPlusPlus::Wrapper_VevaciousPlusPlus*, bool);

#endif /* __BOSS_wrapperutils_vevacious_1_0_hpp__ */

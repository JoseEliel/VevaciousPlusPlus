#include "gambit/Backends/BOSS_wrapperutils.hpp"

VevaciousPlusPlus::Wrapper_VevaciousPlusPlus* wrapper_creator(VevaciousPlusPlus::Abstract_VevaciousPlusPlus* abs_ptr)
{
    return new VevaciousPlusPlus::Wrapper_VevaciousPlusPlus(abs_ptr);
}

void wrapper_deleter(VevaciousPlusPlus::Wrapper_VevaciousPlusPlus* wptr)
{
    wptr->set_delete_BEptr(false);
    delete wptr;
}

void set_delete_BEptr(VevaciousPlusPlus::Wrapper_VevaciousPlusPlus* wptr, bool setting)
{
    wptr->set_delete_BEptr(setting);
}

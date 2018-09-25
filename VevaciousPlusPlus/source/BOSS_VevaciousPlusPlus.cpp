#include <string>
#include <vector>
#include "gambit/Backends/abstracttypedefs.hpp"
#include "gambit/Backends/wrappertypedefs.hpp"
#include "VevaciousPlusPlus.hpp"

void VevaciousPlusPlus::VevaciousPlusPlus::AppendResultsToLhaFile__BOSS(const ::std::basic_string<char, std::char_traits<char>, std::allocator<char> >& lhaFilename)
{
    AppendResultsToLhaFile(lhaFilename);
}




#include "backend_types/VevaciousPlusPlus_1_0/identification.hpp"

VevaciousPlusPlus::Abstract_VevaciousPlusPlus* VevaciousPlusPlus::VevaciousPlusPlus::pointer_copy__BOSS()
{
    VevaciousPlusPlus::Abstract_VevaciousPlusPlus* new_ptr = new VevaciousPlusPlus(*this);
    return new_ptr;
}

void VevaciousPlusPlus::VevaciousPlusPlus::pointer_assign__BOSS(VevaciousPlusPlus::Abstract_VevaciousPlusPlus* in)
{
    CAT_3(BACKENDNAME,_,SAFE_VERSION)::VevaciousPlusPlus::VevaciousPlusPlus* wptr_temp = Abstract_VevaciousPlusPlus::get_wptr();
    *this = *dynamic_cast<VevaciousPlusPlus*>(in);
    Abstract_VevaciousPlusPlus::set_wptr(wptr_temp);
}

#include "gambit/Backends/backend_undefs.hpp"

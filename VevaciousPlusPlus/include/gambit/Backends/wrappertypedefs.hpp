#include "backend_types/vevacious_1_0/forward_decls_wrapper_classes.hpp"
#include "backend_types/vevacious_1_0/identification.hpp"

namespace VevaciousPlusPlus
{
    typedef CAT_3(BACKENDNAME,_,SAFE_VERSION)::VevaciousPlusPlus::VevaciousPlusPlus Wrapper_VevaciousPlusPlus;
}

#include "gambit/Backends/backend_undefs.hpp"

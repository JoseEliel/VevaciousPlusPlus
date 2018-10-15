#include "VevaciousPlusPlus.hpp"
#include "backend_types/vevacious_1_0/wrapper_VevaciousPlusPlus.hpp"
#include <string>
#include "gambit/Backends/abstracttypedefs.hpp"
#include "gambit/Backends/wrappertypedefs.hpp"

extern "C"
{
namespace VevaciousPlusPlus
{
    Abstract_VevaciousPlusPlus* Factory_VevaciousPlusPlus_0__BOSS_1(const ::std::basic_string<char, std::char_traits<char>, std::allocator<char> >& initializationFileName)
    {
        return new VevaciousPlusPlus(initializationFileName);
    }
    
}
}


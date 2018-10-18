#ifndef __wrapper_VevaciousPlusPlus_def_vevacious_1_0_hpp__
#define __wrapper_VevaciousPlusPlus_def_vevacious_1_0_hpp__

#include <string>
#include <vector>

#include "identification.hpp"

namespace CAT_3(BACKENDNAME,_,SAFE_VERSION)
{
    
    namespace VevaciousPlusPlus
    {
        
        // Member functions: 
        inline void VevaciousPlusPlus::RunPoint(const ::std::basic_string<char, std::char_traits<char>, std::allocator<char> >& newInput)
        {
            get_BEptr()->RunPoint(newInput);
        }
        
        inline void VevaciousPlusPlus::ReadLhaBlock(const ::std::basic_string<char, std::char_traits<char>, std::allocator<char> >& uppercaseBlockName, const double scale, const ::std::vector<std::pair<int, double>, std::allocator<std::pair<int, double> > >& parameters, const int dimension)
        {
            get_BEptr()->ReadLhaBlock(uppercaseBlockName, scale, parameters, dimension);
        }
        
        inline void VevaciousPlusPlus::WriteResultsAsXmlFile(const ::std::basic_string<char, std::char_traits<char>, std::allocator<char> >& xmlFilename)
        {
            get_BEptr()->WriteResultsAsXmlFile(xmlFilename);
        }
        
        inline ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > VevaciousPlusPlus::GetResultsAsString()
        {
            return get_BEptr()->GetResultsAsString();
        }

        inline double VevaciousPlusPlus::GetLifetimeInSeconds()
        {
            return get_BEptr()->GetLifetimeInSeconds();
        }

        inline double VevaciousPlusPlus::GetThermalProbability()
        {
            return get_BEptr()->GetThermalProbability();
        }

        inline void VevaciousPlusPlus::AppendResultsToLhaFile(const ::std::basic_string<char, std::char_traits<char>, std::allocator<char> >& lhaFilename, const bool writeWarnings)
        {
            get_BEptr()->AppendResultsToLhaFile(lhaFilename, writeWarnings);
        }
        
        inline void VevaciousPlusPlus::AppendResultsToLhaFile(const ::std::basic_string<char, std::char_traits<char>, std::allocator<char> >& lhaFilename)
        {
            get_BEptr()->AppendResultsToLhaFile__BOSS(lhaFilename);
        }
        
        
        // Wrappers for original constructors: 
        inline VevaciousPlusPlus::VevaciousPlusPlus::VevaciousPlusPlus(const ::std::basic_string<char, std::char_traits<char>, std::allocator<char> >& initializationFileName) :
            WrapperBase(__factory0(initializationFileName))
        {
            get_BEptr()->set_wptr(this);
            get_BEptr()->set_delete_wrapper(false);
        }
        
        // Special pointer-based constructor: 
        inline VevaciousPlusPlus::VevaciousPlusPlus::VevaciousPlusPlus(Abstract_VevaciousPlusPlus* in) :
            WrapperBase(in)
        {
            get_BEptr()->set_wptr(this);
            get_BEptr()->set_delete_wrapper(false);
        }
        
        // Copy constructor: 
        inline VevaciousPlusPlus::VevaciousPlusPlus::VevaciousPlusPlus(const VevaciousPlusPlus& in) :
            WrapperBase(in.get_BEptr()->pointer_copy__BOSS())
        {
            get_BEptr()->set_wptr(this);
            get_BEptr()->set_delete_wrapper(false);
        }
        
        // Assignment operator: 
        inline VevaciousPlusPlus& VevaciousPlusPlus::operator=(const VevaciousPlusPlus& in)
        {
            if (this != &in)
            {
                get_BEptr()->pointer_assign__BOSS(in.get_BEptr());
            }
            return *this;
        }
        
        
        // Destructor: 
        inline VevaciousPlusPlus::VevaciousPlusPlus::~VevaciousPlusPlus()
        {
            if (get_BEptr() != 0)
            {
                get_BEptr()->set_delete_wrapper(false);
                if (can_delete_BEptr())
                {
                    delete BEptr;
                    BEptr = 0;
                }
            }
            set_delete_BEptr(false);
        }
        
        // Returns correctly casted pointer to Abstract class: 
        inline Abstract_VevaciousPlusPlus* VevaciousPlusPlus::VevaciousPlusPlus::get_BEptr() const
        {
            return dynamic_cast<Abstract_VevaciousPlusPlus*>(BEptr);
        }
    }
    
}


#include "gambit/Backends/backend_undefs.hpp"

#endif /* __wrapper_VevaciousPlusPlus_def_vevacious_1_0_hpp__ */

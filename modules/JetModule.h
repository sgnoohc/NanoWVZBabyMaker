#ifndef JetModule_h
#define JetModule_h

#include "rooutil.h"
#include "analysis.h"
#include "Config.h"

namespace wvzModule
{
    //__________________________________________________________________
    // MET module
    class JetModule: public RooUtil::Module
    {
        public:
            JetModule() {}
            virtual void AddOutput();
            virtual void FillOutput();
    };
}

#endif

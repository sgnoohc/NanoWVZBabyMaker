#ifndef AnalysisVariableModule_h
#define AnalysisVariableModule_h

#include "rooutil.h"
#include "analysis.h"
#include "Config.h"
#include "LeptonModule.h"

namespace wvzModule
{
    //__________________________________________________________________
    // AnalysisVariable module
    class AnalysisVariableModule: public RooUtil::Module
    {
        public:
            AnalysisVariableModule() {}
            virtual void AddOutput();
            virtual void FillOutput();
    };
}

#endif

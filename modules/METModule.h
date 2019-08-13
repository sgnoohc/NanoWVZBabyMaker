#ifndef METModule_h
#define METModule_h

#include "rooutil.h"
#include "METCorrectionHandler.h"
#include "analysis.h"
#include "Config.h"
#include "MetSelections.h"

namespace wvzModule
{
    //__________________________________________________________________
    // MET module
    class METModule: public RooUtil::Module
    {
        public:
            METCorrectionHandler metcorrector;
            METObject metobj;
            METObject metobj_corrected;
            METModule() { metcorrector.setup(gconf.year, TString::Format("%d", gconf.year), "StopAnalysis/StopCORE/METCorr/METSFs/"); }
            virtual void AddOutput();
            virtual void FillOutput();
            virtual void correctMET();
    };
}

#endif

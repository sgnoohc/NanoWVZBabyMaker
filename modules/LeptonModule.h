#ifndef LeptonModule_h
#define LeptonModule_h

#include "rooutil.h"
#include "analysis.h"
#include "Config.h"

namespace wvzModule
{
    //__________________________________________________________________
    // Lepton module
    class LeptonModule: public RooUtil::Module
    {
        public:
            LeptonModule() {}
            virtual void AddOutput();
            virtual void FillOutput();
            static bool isPt10AnalysisVetoMuon(int idx);
            static bool isPt10AnalysisVetoElectron(int idx);
            static bool isPt10AnalysisNominalMuon(int idx);
            static bool isPt10AnalysisNominalElectron(int idx);
            static int passCount(const std::vector<int>& v);

    };
}

#endif

#include "BDTInputVariableModule.h"

//==============================================================================================================
//
// BDTInputVariable Module
//
//==============================================================================================================

//##############################################################################################################
void wvzModule::BDTInputVariableModule::AddOutput()
{

    tx->createBranch<int>("lep_Z_idx0");
    tx->createBranch<int>("lep_Z_idx1");
    tx->createBranch<int>("lep_N_idx0");
    tx->createBranch<int>("lep_N_idx1");
    //needed for wwz vs ZZ BDT
    //lep1: Z cand lepton, id > 0
    //lep2: Z cand lepton, id < 0
    //lep3: W cand lepton, leading pt
    //lep4: W cand lepton, subleading pt
    tx->createBranch<float>("MllN");
    tx->createBranch<float>("MllZ");
    tx->createBranch<float>("ZPt");
    tx->createBranch<float>("lep1Pt");
    tx->createBranch<float>("lep2Pt");
    tx->createBranch<float>("lep3Pt");
    tx->createBranch<float>("lep4Pt");
    tx->createBranch<float>("lep3Id");
    tx->createBranch<float>("lep4Id");
    tx->createBranch<float>("lep3MT");
    tx->createBranch<float>("lep4MT");
    tx->createBranch<float>("lep34MT");
    tx->createBranch<float>("lep1dZ");
    tx->createBranch<float>("lep2dZ");
    tx->createBranch<float>("lep3dZ");
    tx->createBranch<float>("lep4dZ");
    tx->createBranch<float>("pt_zeta");
    tx->createBranch<float>("pt_zeta_vis");
    tx->createBranch<float>("phi0");
    tx->createBranch<float>("phi");
    tx->createBranch<float>("phiH");
    tx->createBranch<float>("theta0");
    tx->createBranch<float>("theta1");
    tx->createBranch<float>("theta2");
    //needed for wwz vs ttZ BDT
    //for jet1/2/3/4, if not enough jets, fill with -999
    tx->createBranch<float>("minDRJetToLep3");
    tx->createBranch<float>("minDRJetToLep4");
    tx->createBranch<float>("jet1Pt");
    tx->createBranch<float>("jet2Pt");
    tx->createBranch<float>("jet3Pt");
    tx->createBranch<float>("jet4Pt");
    tx->createBranch<float>("jet1BtagScore");
    tx->createBranch<float>("jet2BtagScore");
    tx->createBranch<float>("jet3BtagScore");
    tx->createBranch<float>("jet4BtagScore");

}

//##############################################################################################################
void wvzModule::BDTInputVariableModule::FillOutput()
{

    // BDT input variables are calculated only when event passes nVlep == 4
    const int& nVlep = tx->getBranch<int>("nVlep");
    const int& nNlep = tx->getBranch<int>("nNlep", false);
    if (nVlep >= 4 and nNlep >= 2)
    {

        const int& lep_ZCand_idx1 = tx->getBranch<int>("lep_ZCand_idx1");
        const int& lep_ZCand_idx2 = tx->getBranch<int>("lep_ZCand_idx2");
        const int& lep_Nom_idx1 = tx->getBranch<int>("lep_Nom_idx1");
        const int& lep_Nom_idx2 = tx->getBranch<int>("lep_Nom_idx2");
        const std::vector<int>& lep_id = tx->getBranch<vector<int>>("lep_id");
        const std::vector<float>& lep_pt = tx->getBranch<vector<float>>("lep_pt");
        const std::vector<float>& lep_eta = tx->getBranch<vector<float>>("lep_eta");
        const std::vector<float>& lep_phi = tx->getBranch<vector<float>>("lep_phi");
        const std::vector<LV>& lep_p4 = tx->getBranch<vector<LV>>("lep_p4");
        const std::vector<LV>& jets_p4 = tx->getBranch<vector<LV>>("jets_p4", false);
        const std::vector<float>& jets_btag_score = tx->getBranch<vector<float>>("jets_btag_score", false);
        const std::vector<float>& lep_dz = tx->getBranch<vector<float>>("lep_dz");
        const float& met_pt = tx->getBranch<float>("met_pt");
        const float& met_phi = tx->getBranch<float>("met_phi");

        int lep1_idx = lep_ZCand_idx1;
        int lep2_idx = lep_ZCand_idx2;
        int lep3_idx = lep_Nom_idx1;
        int lep4_idx = lep_Nom_idx2;

        //Z: lep1 is id>0 lepton
        if(lep_id.at(lep_ZCand_idx1) < 0) 
        {
            lep1_idx = lep_ZCand_idx2;
            lep2_idx = lep_ZCand_idx1;
        } 
        //W: lep3 is higher pt lepton
        if(lep_pt.at(lep_Nom_idx1) < lep_pt.at(lep_Nom_idx2)) 
        {
            lep3_idx = lep_Nom_idx2;
            lep4_idx = lep_Nom_idx1;
        } 

        // For MT calculation
        auto VarMT = [&](LV system)
        {
            return sqrt(2 * met_pt * system.Et() * (1.0 - cos(system.Phi() - met_phi)));
        };

        // X->\tautau reconstruction related variable
        auto VarPtZeta = [&](bool onlyvis)
        {
            TLorentzVector metv3_, lep1_, lep2_;
            metv3_.SetPtEtaPhiM(met_pt, 0., met_phi, 0);
            lep1_.SetPtEtaPhiM(lep_pt.at(lep_Nom_idx1), 0,  lep_phi.at(lep_Nom_idx1), 0);
            lep2_.SetPtEtaPhiM(lep_pt.at(lep_Nom_idx2), 0,  lep_phi.at(lep_Nom_idx2), 0);
            LV metv3 = RooUtil::Calc::getLV(metv3_);
            LV lep1 = RooUtil::Calc::getLV(lep1_);
            LV lep2 = RooUtil::Calc::getLV(lep2_);
            float lep1mag = sqrt(lep1.Px()*lep1.Px() + lep1.Py()*lep1.Py() + lep1.Pz()*lep1.Pz());
            float lep2mag = sqrt(lep2.Px()*lep2.Px() + lep2.Py()*lep2.Py() + lep2.Pz()*lep2.Pz());
            LV zeta = lep1 * lep2mag + lep2 * lep1mag; // find bisector
            float zetamag = sqrt(zeta.Px()*zeta.Px() + zeta.Py()*zeta.Py() + zeta.Pz()*zeta.Pz());
            LV sum_vis = lep1 + lep2;
            LV sum = onlyvis ? sum_vis : sum_vis + metv3;
            return  (sum.Px()*zeta.Px() + sum.Py()*zeta.Py() + sum.Pz()*zeta.Pz()) / zetamag;
        };

        auto VarMinDRJetsToLep = [&](int idx)
        {
            float minDR = 9999;
            for (unsigned int jj = 0 ; jj < jets_p4.size(); jj++)
            {
                float dphi = jets_p4[jj].phi() - lep_phi.at(idx);
                float deta = jets_p4[jj].eta() - lep_eta.at(idx);
                float dR_this = sqrt(dphi * dphi + deta * deta);
                if (dR_this < minDR) minDR = dR_this;
            }
            if (minDR > 9998)
                minDR = -999; //default value: -999
            return minDR;
        };

        LeptonVectors vLeptons = GetLeptonVectors();
        HZZ4lEventParameters eventParameters = ConvertVectorsToAngles(vLeptons);
        //start filling...
        tx->setBranch<int>("lep_Z_idx0", lep1_idx);
        tx->setBranch<int>("lep_Z_idx1", lep2_idx);
        tx->setBranch<int>("lep_N_idx0", lep3_idx);
        tx->setBranch<int>("lep_N_idx1", lep4_idx);
        //tx->setBranch<float>("eventweight", this->EventWeight());
        //tx->setBranch<float>("lepsf", this->LeptonScaleFactor());
        tx->setBranch<float>("MllN", (lep_p4.at(lep_Nom_idx1) + lep_p4.at(lep_Nom_idx2)).M());
        tx->setBranch<float>("MllZ", (lep_p4.at(lep_ZCand_idx1) + lep_p4.at(lep_ZCand_idx2)).M());
        tx->setBranch<float>("ZPt", (lep_p4.at(lep_ZCand_idx1) + lep_p4.at(lep_ZCand_idx2)).Pt());
        tx->setBranch<float>("lep1Pt", lep_pt.at(lep1_idx));
        tx->setBranch<float>("lep2Pt", lep_pt.at(lep2_idx));
        tx->setBranch<float>("lep3Pt", lep_pt.at(lep3_idx));
        tx->setBranch<float>("lep4Pt", lep_pt.at(lep4_idx));
        tx->setBranch<float>("lep3Id", lep_id.at(lep3_idx));
        tx->setBranch<float>("lep4Id", lep_id.at(lep4_idx));
        tx->setBranch<float>("lep3MT", VarMT(lep_p4.at(lep3_idx)));
        tx->setBranch<float>("lep4MT", VarMT(lep_p4.at(lep4_idx)));
        tx->setBranch<float>("lep34MT", VarMT((lep_p4.at(lep3_idx) + lep_p4.at(lep4_idx))));
        tx->setBranch<float>("lep1dZ", lep_dz.at(lep1_idx));
        tx->setBranch<float>("lep2dZ", lep_dz.at(lep2_idx));
        tx->setBranch<float>("lep3dZ", lep_dz.at(lep3_idx));
        tx->setBranch<float>("lep4dZ", lep_dz.at(lep4_idx));
        tx->setBranch<float>("pt_zeta", VarPtZeta(false));
        tx->setBranch<float>("pt_zeta_vis", VarPtZeta(true));
        tx->setBranch<float>("phi0", eventParameters.Phi0);
        tx->setBranch<float>("phi", eventParameters.Phi);
        tx->setBranch<float>("phiH", eventParameters.PhiH);
        tx->setBranch<float>("theta0", eventParameters.Theta0);
        tx->setBranch<float>("theta1", eventParameters.Theta1);
        tx->setBranch<float>("theta2", eventParameters.Theta2);
        tx->setBranch<float>("minDRJetToLep3", VarMinDRJetsToLep(lep3_idx));
        tx->setBranch<float>("minDRJetToLep4", VarMinDRJetsToLep(lep4_idx));
        tx->setBranch<float>("jet1Pt", jets_p4.size() > 0 ? jets_p4[0].pt() : -999);
        tx->setBranch<float>("jet2Pt", jets_p4.size() > 1 ? jets_p4[1].pt() : -999);
        tx->setBranch<float>("jet3Pt", jets_p4.size() > 2 ? jets_p4[2].pt() : -999);
        tx->setBranch<float>("jet4Pt", jets_p4.size() > 3 ? jets_p4[3].pt() : -999);
        tx->setBranch<float>("jet1BtagScore", jets_btag_score.size() > 0 ? jets_btag_score[0] : -999);
        tx->setBranch<float>("jet2BtagScore", jets_btag_score.size() > 1 ? jets_btag_score[1] : -999);
        tx->setBranch<float>("jet3BtagScore", jets_btag_score.size() > 2 ? jets_btag_score[2] : -999);
        tx->setBranch<float>("jet4BtagScore", jets_btag_score.size() > 3 ? jets_btag_score[3] : -999);
    }

}


//______________________________________________________________________________________________
LeptonVectors wvzModule::BDTInputVariableModule::GetLeptonVectors()
{
    const int& lep_ZCand_idx1 = tx->getBranch<int>("lep_ZCand_idx1");
    const int& lep_ZCand_idx2 = tx->getBranch<int>("lep_ZCand_idx2");
    const int& lep_Nom_idx1 = tx->getBranch<int>("lep_Nom_idx1");
    const int& lep_Nom_idx2 = tx->getBranch<int>("lep_Nom_idx2");
    const std::vector<int>& lep_id = tx->getBranch<vector<int>>("lep_id");
    const std::vector<float>& lep_pt = tx->getBranch<vector<float>>("lep_pt");
    const std::vector<float>& lep_eta = tx->getBranch<vector<float>>("lep_eta");
    const std::vector<float>& lep_phi = tx->getBranch<vector<float>>("lep_phi");

    int lep11_idx = lep_ZCand_idx1;
    int lep12_idx = lep_ZCand_idx2;
    int lep21_idx = lep_Nom_idx1;
    int lep22_idx = lep_Nom_idx2;

    //Z: lep11 is id>0 lepton
    if (lep_id.at(lep_ZCand_idx1) < 0)
    {
        lep11_idx = lep_ZCand_idx2;
        lep12_idx = lep_ZCand_idx1;
    }
    //W: lep21 id>0 lepton
    if (lep_id.at(lep_Nom_idx1) < 0)
    {
        lep21_idx = lep_Nom_idx2;
        lep22_idx = lep_Nom_idx1;
    }

    double lep11Mass = 0, lep12Mass = 0, lep21Mass = 0, lep22Mass = 0;
    if (abs(lep_id.at(lep11_idx)) == 11) lep11Mass = 0.000511;
    else if (abs(lep_id.at(lep11_idx)) == 13) lep11Mass = 0.1057;
    if (abs(lep_id.at(lep12_idx)) == 11) lep12Mass = 0.000511;
    else if (abs(lep_id.at(lep12_idx)) == 13) lep12Mass = 0.1057;
    if (abs(lep_id.at(lep21_idx)) == 11) lep21Mass = 0.000511;
    else if (abs(lep_id.at(lep21_idx)) == 13) lep21Mass = 0.1057;
    if (abs(lep_id.at(lep22_idx)) == 11) lep22Mass = 0.000511;
    else if (abs(lep_id.at(lep22_idx)) == 13) lep22Mass = 0.1057;
    LeptonVectors vLeptons;
    vLeptons.Lepton11.SetPtEtaPhiMass(lep_pt.at(lep11_idx), lep_eta.at(lep11_idx), lep_phi.at(lep11_idx), lep11Mass);
    vLeptons.Lepton12.SetPtEtaPhiMass(lep_pt.at(lep12_idx), lep_eta.at(lep12_idx), lep_phi.at(lep12_idx), lep12Mass);
    vLeptons.Lepton21.SetPtEtaPhiMass(lep_pt.at(lep21_idx), lep_eta.at(lep21_idx), lep_phi.at(lep21_idx), lep21Mass);
    vLeptons.Lepton22.SetPtEtaPhiMass(lep_pt.at(lep22_idx), lep_eta.at(lep22_idx), lep_phi.at(lep22_idx), lep22Mass);
    return vLeptons;
}

//----------------------------------------------------------------------------
FourVector::FourVector()
{
    P[0] = 0;
    P[1] = 0;
    P[2] = 0;
    P[3] = 0;
}
//----------------------------------------------------------------------------
FourVector::FourVector(double p[4])
{
    P[0] = p[0];
    P[1] = p[1];
    P[2] = p[2];
    P[3] = p[3];
}
//----------------------------------------------------------------------------
FourVector::FourVector(double e, double px, double py, double pz)
{
    P[0] = e;
    P[1] = px;
    P[2] = py;
    P[3] = pz;
}
//----------------------------------------------------------------------------
FourVector::~FourVector()
{
}

//----------------------------------------------------------------------------
void FourVector::SetPtEtaPhi(double pt, double eta, double phi)
{
    SetPtEtaPhiMass(pt, eta, phi, 0);
}
//----------------------------------------------------------------------------
void FourVector::SetPtEtaPhiMass(double pt, double eta, double phi, double mass)
{
    P[1] = pt * cos(phi);
    P[2] = pt * sin(phi);
    P[3] = pt * sinh(eta);
    P[0] = sqrt(mass * mass + SpatialDot(*this));
}
//----------------------------------------------------------------------------
void FourVector::SetPtYPhi(double pt, double y, double phi)
{
    SetPtYPhiMass(pt, y, phi, 0);
}
//----------------------------------------------------------------------------
void FourVector::SetPtYPhiMass(double pt, double y, double phi, double mass)
{
    P[1] = pt * cos(phi);
    P[2] = pt * sin(phi);
    P[0] = sqrt(pt * pt +  mass * mass) * cosh(y);
    P[3] = P[0] * tanh(y);
}
//----------------------------------------------------------------------------
void FourVector::SetSizeEtaPhi(double size, double eta, double phi)
{
    SetSizeEtaPhiMass(size, eta, phi, 0);
}
//----------------------------------------------------------------------------
void FourVector::SetSizeEtaPhiMass(double size, double eta, double phi, double mass)
{
    P[0] = sqrt(size * size + mass * mass);
    P[1] = size / cosh(eta) * cos(phi);
    P[2] = size / cosh(eta) * sin(phi);
    P[3] = size * tanh(eta);
}
//----------------------------------------------------------------------------
void FourVector::SetSizeEtaPhiEnergy(double size, double eta, double phi, double energy)
{
    P[0] = energy;
    P[1] = size / cosh(eta) * cos(phi);
    P[2] = size / cosh(eta) * sin(phi);
    P[3] = size * tanh(eta);
}
//----------------------------------------------------------------------------
void FourVector::SetSizeThetaPhi(double size, double theta, double phi)
{
    SetSizeThetaPhiMass(size, theta, phi, 0);
}
//----------------------------------------------------------------------------
void FourVector::SetSizeThetaPhiMass(double size, double theta, double phi, double mass)
{
    P[0] = sqrt(size * size + mass * mass);
    P[1] = size * sin(theta) * cos(phi);
    P[2] = size * sin(theta) * sin(phi);
    P[3] = size * cos(theta);
}
//----------------------------------------------------------------------------
double &FourVector::operator [](int index)
{
    if (index >= 0 && index <= 3)
        return P[index];
    return P[0];
}
//----------------------------------------------------------------------------
double FourVector::operator [](int index) const
{
    if (index >= 0 && index <= 3)
        return P[index];
    return 0;
}
//----------------------------------------------------------------------------
FourVector &FourVector::operator =(const FourVector &Other)
{
    P[0] = Other.P[0];
    P[1] = Other.P[1];
    P[2] = Other.P[2];
    P[3] = Other.P[3];
    return *this;
}
//----------------------------------------------------------------------------
FourVector FourVector::operator +(const FourVector &Other) const
{
    FourVector Out;
    Out.P[0] = P[0] + Other.P[0];
    Out.P[1] = P[1] + Other.P[1];
    Out.P[2] = P[2] + Other.P[2];
    Out.P[3] = P[3] + Other.P[3];
    return Out;
}
//----------------------------------------------------------------------------
FourVector FourVector::operator -() const
{
    FourVector Out;
    Out.P[0] = -P[0];
    Out.P[1] = -P[1];
    Out.P[2] = -P[2];
    Out.P[3] = -P[3];
    return Out;
}
//----------------------------------------------------------------------------
FourVector FourVector::operator -(const FourVector &Other) const
{
    FourVector Out;
    Out.P[0] = P[0] - Other.P[0];
    Out.P[1] = P[1] - Other.P[1];
    Out.P[2] = P[2] - Other.P[2];
    Out.P[3] = P[3] - Other.P[3];
    return Out;
}
//----------------------------------------------------------------------------
FourVector FourVector::operator *(double Scale) const
{
    FourVector Out;
    Out.P[0] = P[0] * Scale;
    Out.P[1] = P[1] * Scale;
    Out.P[2] = P[2] * Scale;
    Out.P[3] = P[3] * Scale;
    return Out;
}
//----------------------------------------------------------------------------
FourVector FourVector::operator /(double Scale) const
{
    FourVector Out;
    Out.P[0] = P[0] / Scale;
    Out.P[1] = P[1] / Scale;
    Out.P[2] = P[2] / Scale;
    Out.P[3] = P[3] / Scale;
    return Out;
}
//----------------------------------------------------------------------------
double FourVector::GetMass() const
{
    double Mass2 = GetMass2();
    if (Mass2 >= 0)
        return sqrt(Mass2);
    return 0;
}
//----------------------------------------------------------------------------
double FourVector::GetMass2() const
{
    return MetricDot(*this);
}

//----------------------------------------------------------------------------
double FourVector::GetP() const
{
    return sqrt(SpatialDot(*this));
}
//----------------------------------------------------------------------------
double FourVector::GetP2() const
{
    return SpatialDot(*this);
}
//----------------------------------------------------------------------------
double FourVector::GetPT() const
{
    return sqrt(P[1] * P[1] + P[2] * P[2]);
}
//----------------------------------------------------------------------------
double FourVector::GetPT2() const
{
    return (P[1] * P[1] + P[2] * P[2]);
}
//----------------------------------------------------------------------------
double FourVector::GetEta() const
{
    double Momentum = GetP();
    return 0.5 * log((Momentum + P[3]) / (Momentum - P[3]));
}
//----------------------------------------------------------------------------
double FourVector::GetAbsEta() const
{
    double Eta = GetEta();
    if (Eta < 0)
        Eta = -Eta;
    return Eta;
}
//----------------------------------------------------------------------------
double FourVector::GetRapidity() const
{
    return 0.5 * log((P[0] + P[3]) / (P[0] - P[3]));
}
//----------------------------------------------------------------------------
double FourVector::GetY() const
{
    return GetRapidity();
}
//----------------------------------------------------------------------------
double FourVector::GetPhi() const
{
    double PT = GetPT();
    double Angle = acos(P[1] / PT);
    if (P[2] < 0)
        Angle = -Angle;
    return Angle;
}
//----------------------------------------------------------------------------
double FourVector::GetTheta() const
{
    return acos(P[3] / GetP());
}
//----------------------------------------------------------------------------
double FourVector::GetBeta() const
{
    double Gamma = GetGamma();
    return sqrt(1 - 1 / (Gamma * Gamma));
}
//----------------------------------------------------------------------------
double FourVector::GetGamma() const
{
    return P[0] / GetMass();
}
//----------------------------------------------------------------------------
FourVector FourVector::RotateX(double Angle) const
{
    FourVector Out;
    Out.P[0] = P[0];
    Out.P[1] = P[1];
    Out.P[2] = cos(Angle) * P[2] - sin(Angle) * P[3];
    Out.P[3] = sin(Angle) * P[2] + cos(Angle) * P[3];
    return Out;
}

//----------------------------------------------------------------------------
FourVector FourVector::RotateY(double Angle) const
{
    FourVector Out;
    Out.P[0] = P[0];
    Out.P[1] = sin(Angle) * P[3] + cos(Angle) * P[1];
    Out.P[2] = P[2];
    Out.P[3] = cos(Angle) * P[3] - sin(Angle) * P[1];
    return Out;
}

//----------------------------------------------------------------------------
FourVector FourVector::RotateZ(double Angle) const
{
    FourVector Out;
    Out.P[0] = P[0];
    Out.P[1] = cos(Angle) * P[1] - sin(Angle) * P[2];
    Out.P[2] = sin(Angle) * P[1] + cos(Angle) * P[2];
    Out.P[3] = P[3];
    return Out;
}
//----------------------------------------------------------------------------
FourVector FourVector::Rotate(const FourVector Axis, double Angle) const
{
    // rotate "axis" and input to y-z plane, then rotate "axis" to z axis,
    //    rotate input with respect to z axis
    //    and then rotate back
    double Psi = PI / 2 - Axis.GetPhi();
    double Theta = acos(Axis[3] / Axis.GetP());
    return RotateZ(Psi).RotateX(Theta).RotateZ(Angle).RotateX(-Theta).RotateZ(-Psi);
}
//----------------------------------------------------------------------------
FourVector FourVector::BoostX(double Beta) const
{
    double Gamma = BetaToGamma(Beta);
    FourVector Out;
    Out.P[0] = Gamma * P[0] + Beta * Gamma * P[1];
    Out.P[1] = Beta * Gamma * P[0] + Gamma * P[1];
    Out.P[2] = P[2];
    Out.P[3] = P[3];
    return Out;
}
//----------------------------------------------------------------------------
FourVector FourVector::BoostY(double Beta) const
{
    double Gamma = BetaToGamma(Beta);
    FourVector Out;
    Out.P[0] = Gamma * P[0] + Beta * Gamma * P[2];
    Out.P[1] = P[1];
    Out.P[2] = Beta * Gamma * P[0] + Gamma * P[2];
    Out.P[3] = P[3];
    return Out;
}
//----------------------------------------------------------------------------
FourVector FourVector::BoostZ(double Beta) const
{
    double Gamma = BetaToGamma(Beta);
    FourVector Out;
    Out.P[0] = Gamma * P[0] + Beta * Gamma * P[3];
    Out.P[1] = P[1];
    Out.P[2] = P[2];
    Out.P[3] = Beta * Gamma * P[0] + Gamma * P[3];
    return Out;
}
//----------------------------------------------------------------------------
FourVector FourVector::Boost(const FourVector Axis, double Beta) const
{
    if (Axis.GetPT() < 1e-8)  // axis along z direction
    {
        if (Axis[3] > 0)
            return BoostZ(Beta);
        else
            return BoostZ(-Beta);
    }
    double Psi = PI / 2 - Axis.GetPhi();
    double Theta = acos(Axis[3] / Axis.GetP());
    return RotateZ(Psi).RotateX(Theta).BoostZ(Beta).RotateX(-Theta).RotateZ(-Psi);
}
//----------------------------------------------------------------------------
FourVector FourVector::SpatialCross(const FourVector Other) const
{
    FourVector Out;
    Out.P[0] = P[0];
    Out.P[1] = P[2] * Other.P[3] - P[3] * Other.P[2];
    Out.P[2] = P[3] * Other.P[1] - P[1] * Other.P[3];
    Out.P[3] = P[1] * Other.P[2] - P[2] * Other.P[1];
    return Out;
}
//----------------------------------------------------------------------------
FourVector FourVector::SpatialNormalize() const
{
    FourVector Out;
    Out = (*this) / GetP();
    return Out;
}
//----------------------------------------------------------------------------
double FourVector::SpatialDot(const FourVector &Other) const
{
    return P[1] * Other.P[1] + P[2] * Other.P[2] + P[3] * Other.P[3];
}
//----------------------------------------------------------------------------
double FourVector::MetricDot(const FourVector &Other) const
{
    return P[0] * Other.P[0] - SpatialDot(Other);
}
//----------------------------------------------------------------------------
FourVector operator *(double Scale, FourVector P)
{
    return P * Scale;
}
//----------------------------------------------------------------------------
double GetAngle(const FourVector P1, const FourVector P2)
{
    return acos(P1.SpatialDot(P2) / P1.GetP() / P2.GetP());
}

//----------------------------------------------------------------------------
double GetDR(const FourVector P1, const FourVector P2)
{
    double DEta = P1.GetEta() - P2.GetEta();
    double DPhi = GetDPhi(P1, P2);
    return sqrt(DPhi * DPhi + DEta * DEta);
}
//----------------------------------------------------------------------------
double GetDPhi(const FourVector P1, const FourVector P2)
{
    double DPhi = P1.GetPhi() - P2.GetPhi();
    if (DPhi > PI)
        DPhi = 2 * PI - DPhi;
    if (DPhi < -PI)
        DPhi = DPhi + 2 * PI;
    return DPhi;
}
//----------------------------------------------------------------------------
double GetMT(const FourVector P1, const FourVector P2)
{
    double PT1 = P1.GetPT();
    double PT2 = P2.GetPT();
    return sqrt(2 * (PT1 * PT2 - P1[1] * P2[1] - P1[2] * P2[2]));
}
//----------------------------------------------------------------------------
double GetMinRadius(const FourVector P1, const FourVector P2, const FourVector P3)   // in eta-phi space
{
    double Eta1 = P1.GetEta();
    double Phi1 = P1.GetPhi();
    double Eta2 = P2.GetEta();
    double Phi2 = P2.GetPhi();
    double Eta3 = P3.GetEta();
    double Phi3 = P3.GetPhi();
    double BestResidual2 = 9999999;
    double Best1 = 0;
    double Best2 = 0;
    double Best3 = 0;
    for (int i1 = 0; i1 <= 1; i1++)
    {
        for (int i2 = 0; i2 <= 1; i2++)
        {
            for (int i3 = 0; i3 <= 1; i3++)
            {
                double AveragePhi = (Phi1 + Phi2 + Phi3 + (i1 + i2 + i3) * 2 * PI) / 3;
                double Residual2 = (AveragePhi - Phi1 - i1 * 2 * PI) * (AveragePhi - Phi1 - i1 * 2 * PI)
                    + (AveragePhi - Phi2 - i2 * 2 * PI) * (AveragePhi - Phi2 - i2 * 2 * PI)
                    + (AveragePhi - Phi3 - i3 * 2 * PI) * (AveragePhi - Phi3 - i3 * 2 * PI);
                if (Residual2 < BestResidual2)
                {
                    Best1 = i1;
                    Best2 = i2;
                    Best3 = i3;
                    BestResidual2 = Residual2;
                }
            }
        }
    }
    return GetMinRadius(Eta1, Phi1 + Best1 * 2 * PI, Eta2, Phi2 + Best2 * 2 * PI, Eta3, Phi3 + Best3 * 2 * PI);
}
//----------------------------------------------------------------------------
double GetMinRadius(const double X1, const double Y1, const double X2, const double Y2,
        const double X3, const double Y3)
{
    // compare two radii:
    //    - common circle radius
    //    - maximum of edge length (divided by two)
    // return the smaller of the two
    // calculate common circle radius
    double C1X0 = 2 * (X1 - X2);
    double C1Y0 = 2 * (Y1 - Y2);
    double C1 = X1 * X1 + Y1 * Y1 - X2 * X2 - Y2 * Y2;
    double C2X0 = 2 * (X1 - X3);
    double C2Y0 = 2 * (Y1 - Y3);
    double C2 = X1 * X1 + Y1 * Y1 - X3 * X3 - Y3 * Y3;
    double Distance2 = 99999999;
    if (fabs(C1X0 * C2Y0 - C2X0 * C1Y0) > 1e-8)  // Otherwise three points too close to a straight line
    {
        double M = C1X0 * C2Y0 - C2X0 * C1Y0;
        double MX = C1 * C2Y0 - C2 * C1Y0;
        double MY = C1X0 * C2 - C2X0 * C1;
        double X0 = MX / M;
        double Y0 = MY / M;
        Distance2 = (X1 - X0) * (X1 - X0) + (Y1 - Y0) * (Y1 - Y0);
    }
    // calculate max of edge
    double MaxEdge2 = 0;
    if (MaxEdge2 < (X1 - X2) * (X1 - X2) + (Y1 - Y2) * (Y1 - Y2))
        MaxEdge2 = (X1 - X2) * (X1 - X2) + (Y1 - Y2) * (Y1 - Y2);
    if (MaxEdge2 < (X2 - X3) * (X2 - X3) + (Y2 - Y3) * (Y2 - Y3))
        MaxEdge2 = (X2 - X3) * (X2 - X3) + (Y2 - Y3) * (Y2 - Y3);
    if (MaxEdge2 < (X3 - X1) * (X3 - X1) + (Y3 - Y1) * (Y3 - Y1))
        MaxEdge2 = (X3 - X1) * (X3 - X1) + (Y3 - Y1) * (Y3 - Y1);
    MaxEdge2 = MaxEdge2 / 4;
    // minimum of the two
    return sqrt(std::min(MaxEdge2, Distance2));
}
//----------------------------------------------------------------------------
double GetMR(const FourVector P1, const FourVector P2)
{
    double Temp1 = P1[0] * P2[3] - P1[3] * P2[0];
    double Temp2 = P1[3] - P2[3];
    double Temp3 = P1[0] - P2[0];
    return 2 * sqrt(Temp1 * Temp1 / (Temp2 * Temp2 - Temp3 * Temp3));
}
//----------------------------------------------------------------------------
double GetMRStar(const FourVector P1, const FourVector P2)
{
    double Temp1 = P1[0] + P2[0];
    double Temp2 = P1[3] + P2[3];
    double Temp3 = P1.GetPT() * P1.GetPT() - P2.GetPT() * P2.GetPT();
    double Temp4 = (P1 + P2).GetPT();
    return sqrt((Temp1 * Temp1) - (Temp2 * Temp2) - (Temp3 * Temp3) / (Temp4 * Temp4));
}
//----------------------------------------------------------------------------
double Get2011MR(const FourVector P1, const FourVector P2)
{
    return GetMRStar(P1, P2) * GetGammaRStar(P1, P2);
}
//----------------------------------------------------------------------------
double GetISRRemovedMR(const FourVector P1, const FourVector P2, const FourVector POther, double ME3Assumption)
{
    FourVector ME = -(P1 + P2 + POther);
    ME[3] = ME3Assumption;
    ME[0] = ME.GetP();
    FourVector Total = P1 + P2 + ME;
    double Beta = Total.GetPT() / Total[0];
    FourVector Direction = -Total;
    FourVector NewP1 = P1.Boost(Direction, Beta);
    FourVector NewP2 = P2.Boost(Direction, Beta);
    return GetMR(NewP1, NewP2);
}
//----------------------------------------------------------------------------
double GetISRRemoved2011MR(const FourVector P1, const FourVector P2, const FourVector POther, double ME3Assumption)
{
    FourVector ME = -(P1 + P2 + POther);
    ME[3] = ME3Assumption;
    ME[0] = ME.GetP();
    FourVector Total = P1 + P2 + ME;
    double Beta = Total.GetPT() / Total[0];
    FourVector Direction = -Total;
    FourVector NewP1 = P1.Boost(Direction, Beta);
    FourVector NewP2 = P2.Boost(Direction, Beta);
    return Get2011MR(NewP1, NewP2);
}
//----------------------------------------------------------------------------
double GetISR2011MR(const FourVector P1, const FourVector P2, const FourVector ME, const FourVector ISR, int Assumption)
{
    if (Assumption == 1)
    {
        FourVector ME1 = ME;
        ME1[3] = 0;
        ME1[0] = ME1.GetP();
        FourVector Total1 = P1 + P2 + ME1;
        double Beta1 = Total1.GetPT() / Total1[0];
        Total1[3] = 0;
        FourVector BP1 = P1.Boost(-Total1, Beta1);
        FourVector BP2 = P2.Boost(-Total1, Beta1);
        ME1 = ME1.Boost(-Total1, Beta1);
        return Get2011MR(BP1, BP2);
    }
    if (Assumption == 2)
    {
        FourVector ME2 = ME;
        ME2[3] = 0;
        double Bottom = P1.SpatialDot(P2) * P1.SpatialDot(P2) - P1.GetP2() * P2.GetP2();
        double Projection1Top = P2.SpatialDot(ME2) * P2.SpatialDot(P1) - P1.SpatialDot(ME2) * P2.GetP2();
        double Projection2Top = P1.SpatialDot(ME2) * P1.SpatialDot(P2) - P2.SpatialDot(ME2) * P1.GetP2();
        FourVector ME2P1rojection = P1 * Projection1Top / Bottom;
        FourVector ME2P2rojection = P2 * Projection2Top / Bottom;
        ME2[0] = ME2P1rojection.GetP() + ME2P2rojection.GetP();
        FourVector Total2 = P1 + P2 + ME2;
        double Beta2 = Total2.GetPT() / Total2[0];
        Total2[3] = 0;
        FourVector P12 = P1.Boost(-Total2, Beta2);
        FourVector P22 = P2.Boost(-Total2, Beta2);
        ME2 = ME2.Boost(-Total2, Beta2);
        return Get2011MR(P12, P22);
    }
    if (Assumption == 3)
    {
        FourVector ME3 = ME;
        ME3[3] = -(P1[3] + P2[3]);
        ME3[0] = ME3.GetP();
        FourVector Total3 = P1 + P2 + ME3;
        double Beta3 = Total3.GetPT() / Total3[0];
        Total3[3] = 0;
        FourVector P13 = P1.Boost(-Total3, Beta3);
        FourVector P23 = P2.Boost(-Total3, Beta3);
        ME3 = ME3.Boost(-Total3, Beta3);
        return Get2011MR(P13, P23);
    }
    if (Assumption == 4)
    {
        FourVector ME4 = ME;
        ME4[3] = -(P1[3] + P2[3]);
        double Bottom4 = P1.SpatialDot(P2) * P1.SpatialDot(P2) - P1.GetP2() * P2.GetP2();
        double Projection1Top4 = P2.SpatialDot(ME4) * P2.SpatialDot(P1) - P1.SpatialDot(ME4) * P2.GetP2();
        double Projection2Top4 = P1.SpatialDot(ME4) * P1.SpatialDot(P2) - P2.SpatialDot(ME4) * P1.GetP2();
        FourVector ME2P1rojection4 = P1 * Projection1Top4 / Bottom4;
        FourVector ME2P2rojection4 = P2 * Projection2Top4 / Bottom4;
        ME4[0] = ME2P1rojection4.GetP() + ME2P2rojection4.GetP();
        FourVector Total4 = P1 + P2 + ME4;
        double Beta4 = Total4.GetPT() / Total4[0];
        Total4[3] = 0;
        FourVector P14 = P1.Boost(-Total4, Beta4);
        FourVector P24 = P2.Boost(-Total4, Beta4);
        ME4 = ME4.Boost(-Total4, Beta4);
        return Get2011MR(P14, P24);
    }
    if (Assumption == 5)
    {
        FourVector P1Temp5 = P1;
        FourVector P2Temp5 = P2;
        FourVector METemp5 = ME;
        FourVector Total5 = P1 + P2 + ME;
        P1Temp5[3] = 0;
        P1Temp5[0] = P1Temp5.GetPT();
        P2Temp5[3] = 0;
        P2Temp5[0] = P2Temp5.GetPT();
        METemp5[3] = 0;
        Total5[3] = 0;
        P1Temp5 = P1Temp5.RotateZ(-Total5.GetPhi());
        P2Temp5 = P2Temp5.RotateZ(-Total5.GetPhi());
        METemp5 = METemp5.RotateZ(-Total5.GetPhi());
        // visible: particle 1 (B1) and 3 (B2), partners: particle 2 and 4
        // rotated so that boost is in the x direction
        // y direction is then trivial
        double P2XPrime5 = (P1Temp5[1] + P2Temp5[1] + METemp5[1]) / 2 - P1Temp5[1];
        double P2YPrime5 = -P1Temp5[2];
        double P4XPrime5 = (P1Temp5[1] + P2Temp5[1] + METemp5[1]) / 2 - P2Temp5[1];
        double P4YPrime5 = -P2Temp5[2];
        double E2 = sqrt(P2XPrime5 * P2XPrime5 + P2YPrime5 * P2YPrime5);
        double E4 = sqrt(P4XPrime5 * P4XPrime5 + P4YPrime5 * P4YPrime5);
        double Beta5 = (P1Temp5[1] + P2XPrime5 + P2Temp5[1] + P4XPrime5) / (P1Temp5[0] + E2 + P2Temp5[0] + E4);
        METemp5[0] = E2 + E4;
        FourVector ME5 = ME;
        ME5[3] = 0;
        ME5[0] = sqrt((E2 + E4) * (E2 + E4) + ME5[3] * ME5[3]);
        FourVector P15 = P1.Boost(-Total5, Beta5);
        FourVector P25 = P2.Boost(-Total5, Beta5);
        ME5 = ME5.Boost(-Total5, Beta5);
        return Get2011MR(P15, P25);
    }
    if (Assumption == 6)
    {
        FourVector P1Temp6 = P1;
        FourVector P2Temp6 = P2;
        FourVector METemp6 = ME;
        FourVector Total6 = P1 + P2 + ME;
        P1Temp6[3] = 0;
        P1Temp6[0] = P1Temp6.GetPT();
        P2Temp6[3] = 0;
        P2Temp6[0] = P2Temp6.GetPT();
        METemp6[3] = 0;
        Total6[3] = 0;
        P1Temp6 = P1Temp6.RotateZ(-Total6.GetPhi());
        P2Temp6 = P2Temp6.RotateZ(-Total6.GetPhi());
        METemp6 = METemp6.RotateZ(-Total6.GetPhi());
        double E13_6 = P1Temp6[0] + P2Temp6[0];
        double P13x_6 = P1Temp6[1] + P2Temp6[1];
        double P24x_6 = METemp6[1];
        double P1234x_6 = P13x_6 + P24x_6;
        double Beta6 = (-E13_6 + sqrt(E13_6 * E13_6 + P24x_6 * P24x_6 - P13x_6 * P13x_6)) / P1234x_6;
        double E24_6 = (-P1234x_6 - Beta6 * E13_6) / Beta6;
        FourVector ME6 = ME;
        ME6[3] = 0;
        ME6[0] = sqrt(E24_6 * E24_6 + ME6[3] * ME6[3]);
        Total6 = P1Temp6 + P2Temp6 + ME6;
        FourVector P16 = P1.Boost(-Total6, fabs(Beta6));
        FourVector P26 = P2.Boost(-Total6, fabs(Beta6));
        ME6 = ME6.Boost(-Total6, fabs(Beta6));
        return Get2011MR(P16, P26);
    }
    if (Assumption == 7)
    {
        FourVector P1Temp7 = P1;
        FourVector P2Temp7 = P2;
        FourVector METemp7 = ME;
        FourVector Total7 = P1 + P2 + ME;
        P1Temp7[3] = 0;
        P1Temp7[0] = P1Temp7.GetPT();
        P2Temp7[3] = 0;
        P2Temp7[0] = P2Temp7.GetPT();
        METemp7[3] = 0;
        Total7[3] = 0;
        P1Temp7 = P1Temp7.RotateZ(-Total7.GetPhi());
        P2Temp7 = P2Temp7.RotateZ(-Total7.GetPhi());
        METemp7 = METemp7.RotateZ(-Total7.GetPhi());
        double PxJet7 = P1Temp7[1] + P2Temp7[1];
        double EJet7 = P1Temp7[0] + P2Temp7[0];
        double Beta7 = (EJet7 - sqrt(EJet7 * EJet7 - PxJet7 * PxJet7 + METemp7[1] * METemp7[1]))
            / (PxJet7 - METemp7[1]);
        METemp7[3] = 0;
        METemp7[0] = EJet7 + Beta7 * (METemp7[1] - PxJet7);
        Total7 = P1Temp7 + P2Temp7 + METemp7;
        FourVector P17 = P1Temp7.Boost(-Total7, -Beta7);
        FourVector P27 = P2Temp7.Boost(-Total7, -Beta7);
        FourVector ME7 = METemp7.Boost(-Total7, -Beta7);
        return Get2011MR(P17, P27);
    }
    if (Assumption == 8)
    {
        FourVector B1PTemp8 = P1;
        FourVector B2PTemp8 = P2;
        FourVector METemp8 = ME;
        FourVector Total8 = P1 + P2 + ME;
        B1PTemp8 = B1PTemp8.RotateZ(-Total8.GetPhi());
        B2PTemp8 = B2PTemp8.RotateZ(-Total8.GetPhi());
        METemp8 = METemp8.RotateZ(-Total8.GetPhi());
        double MinimumBetaX = 1000;
        double MinimumDifference = 999999;
        double SearchCenter = 0;
        double SearchStep = 0.1;
        for (int i = 0; i < 5; i++)
        {
            for (double BetaX = SearchCenter - SearchStep * 10; BetaX < SearchCenter + SearchStep * 10;
                    BetaX = BetaX + SearchStep)
            {
                double Difference = GetDifference8(B1PTemp8, B2PTemp8, METemp8, BetaX);
                if (fabs(Difference) < MinimumDifference)
                {
                    MinimumBetaX = BetaX;
                    MinimumDifference = fabs(Difference);
                }
            }
            SearchCenter = MinimumBetaX;
            SearchStep = SearchStep / 10;
        }
        double Beta8 = MinimumBetaX;
        double GammaX8 = BetaToGamma(Beta8);
        double BetaZ8 = GammaX8 * ((B1PTemp8[0] - B2PTemp8[0]) - Beta8 * (B1PTemp8[1] - B2PTemp8[1]))
            / (P1[3] - P2[3]);
        double GammaZ8 = BetaToGamma(BetaZ8);
        // double Beta8E8 = METemp8[1] + B1PTemp8[1] + B2PTemp8[1] - Beta8 * (B1PTemp8[0] + B2PTemp8[0]);
        // double FinalMEx = GammaX8 * METemp8[1] - GammaX8 * Beta8E8;
        // double FinalMEy = METemp8[2];
        // double FinalMET = sqrt(FinalMEx * FinalMEx + FinalMEy * FinalMEy);
        // double FinalJx = GammaX8 * B1PTemp8[1] - GammaX8 * Beta8 * B1PTemp8[0]
        //    + GammaX8 * B2PTemp8[1] - GammaX8 * Beta8 * B2PTemp8[0];
        // double FinalJy = B1PTemp8[2] + B2PTemp8[2];
        // double FinalJT = sqrt(FinalJx * FinalJx + FinalJy * FinalJy);
        double MR8 = GammaZ8 * ((GammaX8 * (B1PTemp8[0] + B2PTemp8[0])
                    - GammaX8 * Beta8 * (B1PTemp8[1] + B2PTemp8[1]))
                - BetaZ8 * (P1[3] + P2[3]));
        // double MT8 = sqrt(2 * (FinalMET * FinalJT - FinalMEx * FinalJx - FinalMEy * FinalJy));
        return MR8;
    }
    if (Assumption == 9)
    {
        FourVector B1PTemp9 = P1;
        FourVector B2PTemp9 = P2;
        FourVector METemp9 = ME;
        FourVector Total9 = P1 + P2 + ME;
        B1PTemp9 = B1PTemp9.RotateZ(-Total9.GetPhi());
        B2PTemp9 = B2PTemp9.RotateZ(-Total9.GetPhi());
        METemp9 = METemp9.RotateZ(-Total9.GetPhi());
        double MinimumBetaZ = 1000;
        double MinimumDifference9 = 99999999;
        double SearchCenter9 = 0;
        double SearchStep9 = 0.05;
        for (int i = 0; i < 5; i++)
        {
            for (double BetaZ = SearchCenter9 - SearchStep9 * 20; BetaZ < SearchCenter9 + SearchStep9 * 20;
                    BetaZ = BetaZ + SearchStep9)
            {
                double Difference = GetDifference9(B1PTemp9, B2PTemp9, METemp9, BetaZ);
                if (fabs(Difference) < MinimumDifference9)
                {
                    MinimumBetaZ = BetaZ;
                    MinimumDifference9 = fabs(Difference);
                }
            }
            SearchCenter9 = MinimumBetaZ;
            SearchStep9 = SearchStep9 / 10;
        }
        double GammaZ9 = BetaToGamma(MinimumBetaZ);
        double BetaX9 = (GammaZ9 * (B1PTemp9[0] - B2PTemp9[0])
                - GammaZ9 * MinimumBetaZ * (B1PTemp9[3] - B2PTemp9[3]))
            / (B1PTemp9[1] - B2PTemp9[1]);
        double GammaX9 = BetaToGamma(BetaX9);
        double FinalJE = GammaX9 * (GammaZ9 * (B1PTemp9[0] + B2PTemp9[0])
                - GammaZ9 * MinimumBetaZ * (B1PTemp9[3] - B2PTemp9[3]))
            - GammaX9 * BetaX9 * (B1PTemp9[1] + B2PTemp9[0]);
        //double FinalJx9 = GammaX9 * (B1PTemp9[1] + B2PTemp9[1]) - GammaX9 * BetaX9
        //   * (GammaZ9 * (B1PTemp9[0] + B2PTemp9[0]) - GammaZ9 * MinimumBetaZ * (B1PTemp9[3] - B2PTemp9[3]));
        //double FinalJy9 = B1PTemp9[2] + B2PTemp9[2];
        // double FinalJT9 = sqrt(FinalJx9 * FinalJx9 + FinalJy9 * FinalJy9);
        double MR9 = FinalJE;
        // double MT9 = 2 * FinalJT9;
        // double R9 = MT9 / MR9 / 2;
        return MR9;
    }
    if (Assumption == 11)
    {
        double JJMass2 = (P1 + P2).GetMass2();
        FourVector METemp11 = ME;
        METemp11[3] = FindMR11MinimumPz(P1, P2, METemp11, ISR);
        METemp11[0] = sqrt(JJMass2 + METemp11.GetP2());
        return EstimateMass11(P1, P2, METemp11, ISR);
    }
    if (Assumption == -11)
    {
        double JJMass2 = (P1 + P2).GetMass2();
        FourVector METemp11 = ME;
        METemp11[3] = FindMR11MinimumPz(P1, P2, METemp11, ISR);
        METemp11[0] = sqrt(JJMass2 + METemp11.GetP2());
        return EstimateMass11(P1, P2, METemp11, ISR);
    }
    return 0;
}
//----------------------------------------------------------------------------
double GetMRT(const FourVector P1, const FourVector P2, const FourVector ME)
{
    double Temp1 = ME.GetPT() * (P1.GetPT() + P2.GetPT());
    double Temp2 = ME[1] * (P1[1] + P2[1]) + ME[2] * (P1[2] + P2[2]);
    return sqrt((Temp1 - Temp2) / 2);
}
//----------------------------------------------------------------------------
double GetR(const FourVector P1, const FourVector P2, const FourVector ME)
{
    return GetMRT(P1, P2, ME) / GetMR(P1, P2);
}
//----------------------------------------------------------------------------
double GetRStar(const FourVector P1, const FourVector P2, const FourVector ME)
{
    return GetMRT(P1, P2, ME) / GetMRStar(P1, P2) / GetGammaRStar(P1, P2);
}
//----------------------------------------------------------------------------
double Get2011R(const FourVector P1, const FourVector P2, const FourVector ME)
{
    return GetMRT(P1, P2, ME) / GetMRStar(P1, P2) / GetGammaRStar(P1, P2);
}
//----------------------------------------------------------------------------
double GetISRRemovedR(const FourVector P1, const FourVector P2, const FourVector POther, double ME3Assumption)
{
    FourVector ME = -(P1 + P2 + POther);
    ME[3] = ME3Assumption;
    ME[0] = ME.GetP();
    FourVector Total = P1 + P2 + ME;
    double Beta = Total.GetPT() / Total[0];
    FourVector Direction = -Total;
    FourVector NewP1 = P1.Boost(Direction, Beta);
    FourVector NewP2 = P2.Boost(Direction, Beta);
    FourVector NewME = ME.Boost(Direction, Beta);
    return GetR(NewP1, NewP2, NewME);
}
//----------------------------------------------------------------------------
double GetISRRemoved2011R(const FourVector P1, const FourVector P2, const FourVector POther, double ME3Assumption)
{
    FourVector ME = -(P1 + P2 + POther);
    ME[3] = ME3Assumption;
    ME[0] = ME.GetP();
    FourVector Total = P1 + P2 + ME;
    double Beta = Total.GetPT() / Total[0];
    FourVector Direction = -Total;
    FourVector NewP1 = P1.Boost(Direction, Beta);
    FourVector NewP2 = P2.Boost(Direction, Beta);
    FourVector NewME = ME.Boost(Direction, Beta);
    return Get2011R(NewP1, NewP2, NewME);
}
//----------------------------------------------------------------------------
double GetISR2011R(const FourVector P1, const FourVector P2, const FourVector ME, const FourVector ISR, int Assumption, char AdditionalVariant)
{
    if (Assumption == 1)
    {
        FourVector ME1 = ME;
        ME1[3] = 0;
        ME1[0] = ME1.GetP();
        FourVector Total1 = P1 + P2 + ME1;
        double Beta1 = Total1.GetPT() / Total1[0];
        Total1[3] = 0;
        FourVector BP1 = P1.Boost(-Total1, Beta1);
        FourVector BP2 = P2.Boost(-Total1, Beta1);
        ME1 = ME1.Boost(-Total1, Beta1);
        return Get2011R(BP1, BP2, ME1);
    }
    if (Assumption == 2)
    {
        FourVector ME2 = ME;
        ME2[3] = 0;
        double Bottom = P1.SpatialDot(P2) * P1.SpatialDot(P2) - P1.GetP2() * P2.GetP2();
        double Projection1Top = P2.SpatialDot(ME2) * P2.SpatialDot(P1) - P1.SpatialDot(ME2) * P2.GetP2();
        double Projection2Top = P1.SpatialDot(ME2) * P1.SpatialDot(P2) - P2.SpatialDot(ME2) * P1.GetP2();
        FourVector ME2P1rojection = P1 * Projection1Top / Bottom;
        FourVector ME2P2rojection = P2 * Projection2Top / Bottom;
        ME2[0] = ME2P1rojection.GetP() + ME2P2rojection.GetP();
        FourVector Total2 = P1 + P2 + ME2;
        double Beta2 = Total2.GetPT() / Total2[0];
        Total2[3] = 0;
        FourVector P12 = P1.Boost(-Total2, Beta2);
        FourVector P22 = P2.Boost(-Total2, Beta2);
        ME2 = ME2.Boost(-Total2, Beta2);
        return Get2011R(P12, P22, ME2);
    }
    if (Assumption == 3)
    {
        FourVector ME3 = ME;
        ME3[3] = -(P1[3] + P2[3]);
        ME3[0] = ME3.GetP();
        FourVector Total3 = P1 + P2 + ME3;
        double Beta3 = Total3.GetPT() / Total3[0];
        Total3[3] = 0;
        FourVector P13 = P1.Boost(-Total3, Beta3);
        FourVector P23 = P2.Boost(-Total3, Beta3);
        ME3 = ME3.Boost(-Total3, Beta3);
        return Get2011R(P13, P23, ME3);
    }
    if (Assumption == 4)
    {
        FourVector ME4 = ME;
        ME4[3] = -(P1[3] + P2[3]);
        double Bottom4 = P1.SpatialDot(P2) * P1.SpatialDot(P2) - P1.GetP2() * P2.GetP2();
        double Projection1Top4 = P2.SpatialDot(ME4) * P2.SpatialDot(P1) - P1.SpatialDot(ME4) * P2.GetP2();
        double Projection2Top4 = P1.SpatialDot(ME4) * P1.SpatialDot(P2) - P2.SpatialDot(ME4) * P1.GetP2();
        FourVector ME2P1rojection4 = P1 * Projection1Top4 / Bottom4;
        FourVector ME2P2rojection4 = P2 * Projection2Top4 / Bottom4;
        ME4[0] = ME2P1rojection4.GetP() + ME2P2rojection4.GetP();
        FourVector Total4 = P1 + P2 + ME4;
        double Beta4 = Total4.GetPT() / Total4[0];
        Total4[3] = 0;
        FourVector P14 = P1.Boost(-Total4, Beta4);
        FourVector P24 = P2.Boost(-Total4, Beta4);
        ME4 = ME4.Boost(-Total4, Beta4);
        return Get2011R(P14, P24, ME4);
    }
    if (Assumption == 5)
    {
        FourVector P1Temp5 = P1;
        FourVector P2Temp5 = P2;
        FourVector METemp5 = ME;
        FourVector Total5 = P1 + P2 + ME;
        P1Temp5[3] = 0;
        P1Temp5[0] = P1Temp5.GetPT();
        P2Temp5[3] = 0;
        P2Temp5[0] = P2Temp5.GetPT();
        METemp5[3] = 0;
        Total5[3] = 0;
        P1Temp5 = P1Temp5.RotateZ(-Total5.GetPhi());
        P2Temp5 = P2Temp5.RotateZ(-Total5.GetPhi());
        METemp5 = METemp5.RotateZ(-Total5.GetPhi());
        // visible: particle 1 (B1) and 3 (B2), partners: particle 2 and 4
        // rotated so that boost is in the x direction
        // y direction is then trivial
        double P2XPrime5 = (P1Temp5[1] + P2Temp5[1] + METemp5[1]) / 2 - P1Temp5[1];
        double P2YPrime5 = -P1Temp5[2];
        double P4XPrime5 = (P1Temp5[1] + P2Temp5[1] + METemp5[1]) / 2 - P2Temp5[1];
        double P4YPrime5 = -P2Temp5[2];
        double E2 = sqrt(P2XPrime5 * P2XPrime5 + P2YPrime5 * P2YPrime5);
        double E4 = sqrt(P4XPrime5 * P4XPrime5 + P4YPrime5 * P4YPrime5);
        double Beta5 = (P1Temp5[1] + P2XPrime5 + P2Temp5[1] + P4XPrime5) / (P1Temp5[0] + E2 + P2Temp5[0] + E4);
        METemp5[0] = E2 + E4;
        FourVector ME5 = ME;
        ME5[3] = 0;
        ME5[0] = sqrt((E2 + E4) * (E2 + E4) + ME5[3] * ME5[3]);
        FourVector P15 = P1.Boost(-Total5, Beta5);
        FourVector P25 = P2.Boost(-Total5, Beta5);
        ME5 = ME5.Boost(-Total5, Beta5);
        return Get2011R(P15, P25, ME5);
    }
    if (Assumption == 6)
    {
        FourVector P1Temp6 = P1;
        FourVector P2Temp6 = P2;
        FourVector METemp6 = ME;
        FourVector Total6 = P1 + P2 + ME;
        P1Temp6[3] = 0;
        P1Temp6[0] = P1Temp6.GetPT();
        P2Temp6[3] = 0;
        P2Temp6[0] = P2Temp6.GetPT();
        METemp6[3] = 0;
        Total6[3] = 0;
        P1Temp6 = P1Temp6.RotateZ(-Total6.GetPhi());
        P2Temp6 = P2Temp6.RotateZ(-Total6.GetPhi());
        METemp6 = METemp6.RotateZ(-Total6.GetPhi());
        double E13_6 = P1Temp6[0] + P2Temp6[0];
        double P13x_6 = P1Temp6[1] + P2Temp6[1];
        double P24x_6 = METemp6[1];
        double P1234x_6 = P13x_6 + P24x_6;
        double Beta6 = (-E13_6 + sqrt(E13_6 * E13_6 + P24x_6 * P24x_6 - P13x_6 * P13x_6)) / P1234x_6;
        double E24_6 = (-P1234x_6 - Beta6 * E13_6) / Beta6;
        FourVector ME6 = ME;
        ME6[3] = 0;
        ME6[0] = sqrt(E24_6 * E24_6 + ME6[3] * ME6[3]);
        Total6 = P1Temp6 + P2Temp6 + ME6;
        FourVector P16 = P1.Boost(-Total6, fabs(Beta6));
        FourVector P26 = P2.Boost(-Total6, fabs(Beta6));
        ME6 = ME6.Boost(-Total6, fabs(Beta6));
        return Get2011R(P16, P26, ME6);
    }
    if (Assumption == 7)
    {
        FourVector P1Temp7 = P1;
        FourVector P2Temp7 = P2;
        FourVector METemp7 = ME;
        FourVector Total7 = P1 + P2 + ME;
        P1Temp7[3] = 0;
        P1Temp7[0] = P1Temp7.GetPT();
        P2Temp7[3] = 0;
        P2Temp7[0] = P2Temp7.GetPT();
        METemp7[3] = 0;
        Total7[3] = 0;
        P1Temp7 = P1Temp7.RotateZ(-Total7.GetPhi());
        P2Temp7 = P2Temp7.RotateZ(-Total7.GetPhi());
        METemp7 = METemp7.RotateZ(-Total7.GetPhi());
        double PxJet7 = P1Temp7[1] + P2Temp7[1];
        double EJet7 = P1Temp7[0] + P2Temp7[0];
        double Beta7 = (EJet7 - sqrt(EJet7 * EJet7 - PxJet7 * PxJet7 + METemp7[1] * METemp7[1]))
            / (PxJet7 - METemp7[1]);
        METemp7[3] = 0;
        METemp7[0] = EJet7 + Beta7 * (METemp7[1] - PxJet7);
        Total7 = P1Temp7 + P2Temp7 + METemp7;
        FourVector P17 = P1Temp7.Boost(-Total7, -Beta7);
        FourVector P27 = P2Temp7.Boost(-Total7, -Beta7);
        FourVector ME7 = METemp7.Boost(-Total7, -Beta7);
        return Get2011R(P17, P27, ME7);
    }
    if (Assumption == 8)
    {
        FourVector B1PTemp8 = P1;
        FourVector B2PTemp8 = P2;
        FourVector METemp8 = ME;
        FourVector Total8 = P1 + P2 + ME;
        B1PTemp8 = B1PTemp8.RotateZ(-Total8.GetPhi());
        B2PTemp8 = B2PTemp8.RotateZ(-Total8.GetPhi());
        METemp8 = METemp8.RotateZ(-Total8.GetPhi());
        double MinimumBetaX = 1000;
        double MinimumDifference = 999999;
        double SearchCenter = 0;
        double SearchStep = 0.1;
        for (int i = 0; i < 5; i++)
        {
            for (double BetaX = SearchCenter - SearchStep * 10; BetaX < SearchCenter + SearchStep * 10;
                    BetaX = BetaX + SearchStep)
            {
                double Difference = GetDifference8(B1PTemp8, B2PTemp8, METemp8, BetaX);
                if (fabs(Difference) < MinimumDifference)
                {
                    MinimumBetaX = BetaX;
                    MinimumDifference = fabs(Difference);
                }
            }
            SearchCenter = MinimumBetaX;
            SearchStep = SearchStep / 10;
        }
        double Beta8 = MinimumBetaX;
        double GammaX8 = BetaToGamma(Beta8);
        double BetaZ8 = GammaX8 * ((B1PTemp8[0] - B2PTemp8[0]) - Beta8 * (B1PTemp8[1] - B2PTemp8[1]))
            / (P1[3] - P2[3]);
        double GammaZ8 = BetaToGamma(BetaZ8);
        double Beta8E8 = METemp8[1] + B1PTemp8[1] + B2PTemp8[1] - Beta8 * (B1PTemp8[0] + B2PTemp8[0]);
        double FinalMEx = GammaX8 * METemp8[1] - GammaX8 * Beta8E8;
        double FinalMEy = METemp8[2];
        double FinalMET = sqrt(FinalMEx * FinalMEx + FinalMEy * FinalMEy);
        double FinalJx = GammaX8 * B1PTemp8[1] - GammaX8 * Beta8 * B1PTemp8[0]
            + GammaX8 * B2PTemp8[1] - GammaX8 * Beta8 * B2PTemp8[0];
        double FinalJy = B1PTemp8[2] + B2PTemp8[2];
        double FinalJT = sqrt(FinalJx * FinalJx + FinalJy * FinalJy);
        double MR8 = GammaZ8 * ((GammaX8 * (B1PTemp8[0] + B2PTemp8[0])
                    - GammaX8 * Beta8 * (B1PTemp8[1] + B2PTemp8[1]))
                - BetaZ8 * (P1[3] + P2[3]));
        double MT8 = sqrt(2 * (FinalMET * FinalJT - FinalMEx * FinalJx - FinalMEy * FinalJy));
        return MT8 / MR8 / 2;
    }
    if (Assumption == 9)
    {
        FourVector B1PTemp9 = P1;
        FourVector B2PTemp9 = P2;
        FourVector METemp9 = ME;
        FourVector Total9 = P1 + P2 + ME;
        B1PTemp9 = B1PTemp9.RotateZ(-Total9.GetPhi());
        B2PTemp9 = B2PTemp9.RotateZ(-Total9.GetPhi());
        METemp9 = METemp9.RotateZ(-Total9.GetPhi());
        double MinimumBetaZ = 1000;
        double MinimumDifference9 = 99999999;
        double SearchCenter9 = 0;
        double SearchStep9 = 0.05;
        for (int i = 0; i < 5; i++)
        {
            for (double BetaZ = SearchCenter9 - SearchStep9 * 20; BetaZ < SearchCenter9 + SearchStep9 * 20;
                    BetaZ = BetaZ + SearchStep9)
            {
                double Difference = GetDifference9(B1PTemp9, B2PTemp9, METemp9, BetaZ);
                if (fabs(Difference) < MinimumDifference9)
                {
                    MinimumBetaZ = BetaZ;
                    MinimumDifference9 = fabs(Difference);
                }
            }
            SearchCenter9 = MinimumBetaZ;
            SearchStep9 = SearchStep9 / 10;
        }
        double GammaZ9 = BetaToGamma(MinimumBetaZ);
        double BetaX9 = (GammaZ9 * (B1PTemp9[0] - B2PTemp9[0])
                - GammaZ9 * MinimumBetaZ * (B1PTemp9[3] - B2PTemp9[3]))
            / (B1PTemp9[1] - B2PTemp9[1]);
        double GammaX9 = BetaToGamma(BetaX9);
        double FinalJE = GammaX9 * (GammaZ9 * (B1PTemp9[0] + B2PTemp9[0])
                - GammaZ9 * MinimumBetaZ * (B1PTemp9[3] - B2PTemp9[3]))
            - GammaX9 * BetaX9 * (B1PTemp9[1] + B2PTemp9[0]);
        double FinalJx9 = GammaX9 * (B1PTemp9[1] + B2PTemp9[1]) - GammaX9 * BetaX9
            * (GammaZ9 * (B1PTemp9[0] + B2PTemp9[0]) - GammaZ9 * MinimumBetaZ * (B1PTemp9[3] - B2PTemp9[3]));
        double FinalJy9 = B1PTemp9[2] + B2PTemp9[2];
        double FinalJT9 = sqrt(FinalJx9 * FinalJx9 + FinalJy9 * FinalJy9);
        double MR9 = FinalJE;
        double MT9 = 2 * FinalJT9;
        double R9 = MT9 / MR9 / 2;
        return R9;
    }
    if (Assumption == 11)
    {
        double JJMass2 = (P1 + P2).GetMass2();
        FourVector METemp11 = ME;
        METemp11[3] = FindMR11MinimumPz(P1, P2, METemp11, ISR);
        METemp11[0] = sqrt(JJMass2 + METemp11.GetP2());
        double MR11 = EstimateMass11(P1, P2, METemp11, ISR);
        double MRT11 = EstimateTransverseMass11(P1, P2, METemp11, ISR, AdditionalVariant);
        return MRT11 / MR11;
    }
    return 0;
}
//----------------------------------------------------------------------------
double GetGammaRStar(const FourVector P1, const FourVector P2)
{
    double Temp1 = P1[0] + P2[0];
    double Temp2 = P1[3] + P2[3];
    double Temp3 = P1.GetPT() * P1.GetPT() - P2.GetPT() * P2.GetPT();
    double Temp4 = (P1 + P2).GetPT();
    double Upper = Temp1 * Temp1 - Temp2 * Temp2;
    double Lower = Temp1 * Temp1 - Temp2 * Temp2 - Temp3 * Temp3 / Temp4 / Temp4;
    return sqrt(Upper / Lower);
    /*
       double Beta2 = (Temp3 * Temp3) / (Temp4 * Temp4 * (Temp1 * Temp1 - Temp2 * Temp2));
       double Gamma = 1.0 / sqrt(1.0 - Beta2);
       return Gamma;
       */
}
//----------------------------------------------------------------------------
double BetaToGamma(double Beta)
{
    return 1 / sqrt(1 - Beta * Beta);
}
//----------------------------------------------------------------------------
double GammaToBeta(double Gamma)
{
    return sqrt(1 - 1 / (Gamma * Gamma));
}
//----------------------------------------------------------------------------
std::vector<FourVector> SplitIntoGroups(std::vector<FourVector> &Input, bool ZeroMass)
{
    std::vector<FourVector> Result;
    if (Input.size() == 0)
    {
        Result.push_back(FourVector());
        Result.push_back(FourVector());
        return Result;
    }
    if (Input.size() == 1)
    {
        Result.push_back(Input[0]);
        Result.push_back(FourVector());
        return Result;
    }
    if (Input.size() == 2)
        return Input;
    // let's start with easy (potentially slow) way: try out all possibilities
    // if speed becomes a problem then we'll come up with something else
    int InputSize = Input.size();
    std::vector<int> Groups(InputSize);
    for (int i = 0; i < InputSize; i++)
        Groups[i] = 0;
    Groups[0] = 1;
    FourVector Group1;
    FourVector Group2;
    double MinMass2 = -1;
    while (Groups[InputSize - 1] == 0) // last one is always in group "0"
    {
        FourVector Vector1Temp;
        FourVector Vector2Temp;
        for (int i = 0; i < InputSize; i++)
        {
            if (Groups[i] == 0)
                Vector1Temp = Vector1Temp + Input[i];
            else
                Vector2Temp = Vector2Temp + Input[i];
        }
        double MinMass2Temp = Vector1Temp.GetMass2() + Vector2Temp.GetMass2();
        if (MinMass2 < 0 || MinMass2Temp < MinMass2)
        {
            MinMass2 = MinMass2Temp;
            Group1 = Vector1Temp;
            Group2 = Vector2Temp;
        }
        Groups[0] = Groups[0] + 1;
        for (int i = 0; i < InputSize - 1; i++)
        {
            while (Groups[i] >= 2)  // just in case something went wrong....an if statement should be enough
            {
                Groups[i] = Groups[i] - 2;
                Groups[i + 1] = Groups[i + 1] + 1;
            }
        }
    }
    if (ZeroMass == true)
    {
        Group1[0] = Group1.GetP();
        Group2[0] = Group2.GetP();
    }
    Result.push_back(Group1);
    Result.push_back(Group2);
    return Result;
}
//----------------------------------------------------------------------------
double GetDifference8(FourVector &P1, FourVector &P2, FourVector &ME, double BetaX)
{
    double DeltaPx = P1[1] - P2[1];
    double Pxj = P1[1] + P2[1];
    double Px = ME[1];
    double DeltaPz = P1[3] - P2[3];
    double Pzj = P1[3] + P2[3];
    double DeltaE = P1[0] - P2[0];
    double Ej = P1[0] + P2[0];
    double GammaX = BetaToGamma(BetaX);
    double DeltaPzBetaZ = GammaX * DeltaE - GammaX * BetaX * DeltaPx;
    double BetaXE = Px + Pxj - BetaX * Ej;
    double BetaXDeltaPzPz = -BetaX * DeltaPz * Pzj
        + DeltaPzBetaZ * GammaX * (BetaXE + BetaX * Ej - BetaX * BetaX * (Px + Pxj));
    double Left = DeltaPz * DeltaPz * BetaX * GammaX * Ej
        - GammaX * BetaX * BetaX * DeltaPz * Pxj * DeltaPz
        - Pzj * BetaX * DeltaPzBetaZ * DeltaPz;
    double Right = GammaX * BetaXE * DeltaPz * DeltaPz
        - GammaX * BetaX * BetaX * Px * DeltaPz * DeltaPz
        - DeltaPzBetaZ * BetaXDeltaPzPz;
    return Left - Right;
}
//----------------------------------------------------------------------------
double GetDifference9(FourVector &P1, FourVector &P2, FourVector &ME, double BetaZ)
{
    double DeltaPx = P1[1] - P2[1];
    double Pxj = P1[1] + P2[1];
    double Px = ME[1];
    double DeltaPz = P1[3] - P2[3];
    double Pzj = P1[3] + P2[3];
    double DeltaE = P1[0] - P2[0];
    double Ej = P1[0] + P2[0];
    double GammaZ = BetaToGamma(BetaZ);
    double InvGammaZ = sqrt(1 - BetaZ * BetaZ);
    double DeltaPxBetaX = GammaZ * DeltaE - GammaZ * BetaZ * DeltaPz;
    double DeltaPxInvGammaE = GammaZ * Ej * DeltaPx * (1 + BetaZ * BetaZ)
        - 2 * GammaZ * BetaZ * Pzj * DeltaPx
        + DeltaPxBetaX * (Px - Pxj);
    double DeltaPxInvGammaPz = DeltaPx * InvGammaZ * Ej * BetaZ
        - DeltaPx * InvGammaZ * Pzj + BetaZ * DeltaPxInvGammaE;
    double Left = DeltaPx * DeltaPx * InvGammaZ * (Pxj + Px);
    double Right = DeltaPxBetaX * (GammaZ * (DeltaPxInvGammaE + DeltaPx * InvGammaZ * Ej)
            - GammaZ * BetaZ * (DeltaPxInvGammaPz + DeltaPx * InvGammaZ * Pzj));
    return fabs(Left - Right);
}
//----------------------------------------------------------------------------
double FindMR11MinimumPz(FourVector J1, FourVector J2, FourVector ME, FourVector ISR)
{
    // do some basic caching to save time repeating search for MR and R
    static FourVector PreviousJ1 = FourVector(0, 0, 0, 0);
    static FourVector PreviousJ2 = FourVector(0, 0, 0, 0);
    static FourVector PreviousME = FourVector(0, 0, 0, 0);
    static FourVector PreviousISR = FourVector(0, 0, 0, 0);
    static double PreviousPz = 0;
    if ((J1 - PreviousJ1).GetP() < 0.1
            && (J2 - PreviousJ2).GetP() < 0.1
            && (ME - PreviousME).GetP() < 0.1
            && (ISR - PreviousISR).GetP() < 0.1)
        return PreviousPz;
    // start calculation
    double JJMass2 = (J1 + J2).GetMass2();
    int InitialStep = 400;
    double InitialStepSize = 5;
    double InitialCenter = 0;
    int SearchStep = 10;
    double SearchStepSize = 1;
    std::vector<double> Masses;
    for (int i = 0; i <= InitialStep; i++)
    {
        ME[3] = InitialCenter - InitialStep / 2 * InitialStepSize + i * InitialStepSize;
        ME[0] = sqrt(JJMass2 + ME.GetP2());
        double Mass = EstimateMass11(J1, J2, ME, ISR);
        Masses.push_back(Mass);
    }
    std::vector<double> LocalMinima;
    for (int i = 1; i < InitialStep; i++)
        if (Masses[i] <= Masses[i - 1] && Masses[i] <= Masses[i + 1])
            LocalMinima.push_back(InitialCenter - InitialStep / 2 * InitialStepSize + i * InitialStepSize);
    if (Masses[0] <= Masses[1])
        LocalMinima.push_back(InitialCenter - InitialStep / 2 * InitialStepSize);
    if (Masses[InitialStep] <= Masses[InitialStep - 1])
        LocalMinima.push_back(InitialCenter + InitialStep / 2 * InitialStepSize);
    for (int i = 0; i <= 5; i++)
    {
        std::vector<double> NewMinima;
        for (int j = 0; j < (int)LocalMinima.size(); j++)
        {
            double SearchCenter = LocalMinima[j];
            Masses.clear();
            for (int k = 0; k <= SearchStep; k++)
            {
                ME[3] = SearchCenter - SearchStep / 2 * SearchStepSize + k * SearchStepSize;
                ME[0] = sqrt(JJMass2 + ME.GetP2());
                double Mass = EstimateMass11(J1, J2, ME, ISR);
                Masses.push_back(Mass);
            }
            for (int k = 1; k < SearchStep; k++)
                if (Masses[k] <= Masses[k - 1] && Masses[k] <= Masses[k + 1])
                    NewMinima.push_back(SearchCenter - SearchStep / 2 * SearchStepSize + k * SearchStepSize);
            if (Masses[0] <= Masses[1])
                NewMinima.push_back(SearchCenter - SearchStep / 2 * SearchStepSize);
            if (Masses[InitialStep] <= Masses[InitialStep - 1])
                NewMinima.push_back(SearchCenter + SearchStep / 2 * SearchStepSize);
        }
        SearchStepSize = SearchStepSize / 5;
        LocalMinima = NewMinima;
    }
    double BestPz = -1;
    double BestMass = -1;
    for (int i = 0; i < (int)LocalMinima.size(); i++)
    {
        ME[3] = LocalMinima[i];
        ME[0] = sqrt(JJMass2 + ME.GetP2());
        double Mass = EstimateMass11(J1, J2, ME, ISR);
        if (Mass < BestMass || BestMass < 0)
        {
            BestPz = LocalMinima[i];
            BestMass = Mass;
        }
    }
    PreviousJ1 = J1;
    PreviousJ2 = J2;
    PreviousME = ME;
    PreviousISR = ISR;
    PreviousPz = BestPz;
    return BestPz;
}
//----------------------------------------------------------------------------
double EstimateMass11(FourVector J1, FourVector J2, FourVector ME, FourVector ISR, bool Reversal)
{
    Reversal = false;
    FourVector TempTotal = ME + J1 + J2 + ISR;
    double TempBetaZ = TempTotal[3] / TempTotal[0];
    FourVector TempJ1 = J1.Boost(FourVector(1, 0, 0, 1), TempBetaZ);
    FourVector TempJ2 = J2.Boost(FourVector(1, 0, 0, 1), TempBetaZ);
    FourVector TempME = ME.Boost(FourVector(1, 0, 0, 1), TempBetaZ);
    FourVector TempISR = ISR.Boost(FourVector(1, 0, 0, 1), TempBetaZ);
    TempTotal = TempJ1 + TempJ2 + TempME;
    double TempBeta = TempTotal.GetP() / TempTotal[0];
    TempJ1 = TempJ1.Boost(TempTotal, TempBeta);
    TempJ2 = TempJ2.Boost(TempTotal, TempBeta);
    TempME = TempME.Boost(TempTotal, TempBeta);
    double EMET = TempME[0];
    double SumE = TempJ1[0] + TempJ2[0];
    double DeltaE = TempJ1[0] - TempJ2[0];
    double ES = (TempJ1 + TempJ2).GetP();
    double ED = (TempJ1 - TempJ2).GetP();
    double A = SumE * DeltaE / ES;
    double B = sqrt(ED * ED - A * A);
    double m0 = (SumE - EMET) / (2 * DeltaE) - SumE * DeltaE / ES / ES;
    double m1 = DeltaE * (SumE + EMET) / 2;
    double EQA = m0 * m0 + B * B / ES / ES;
    double EQB = 2 * m0 * m1 - B * B;
    double EQC = m1 * m1;
    double X2Max = (-EQB + sqrt(EQB * EQB - 4 * EQA * EQC)) / (2 * EQA);
    double X2Min = (-EQB - sqrt(EQB * EQB - 4 * EQA * EQC)) / (2 * EQA);
    if (X2Min > ES * ES)
        return -1;
    if (X2Max < 0)
        return -1;
    if (X2Max > ES * ES)
        X2Max = ES * ES;
    double X = sqrt(X2Max);
    double Y = ((X * X + DeltaE * DeltaE) * SumE - (X * X - DeltaE * DeltaE) * EMET) / (2 * X * DeltaE);
    double M2 = (SumE * X - DeltaE * Y) * (SumE * X - DeltaE * Y) / (X * X - DeltaE * DeltaE) / 4;
    if (M2 < 0)
        return -1;
    return sqrt(M2) * 2;
}
//----------------------------------------------------------------------------
double EstimateTransverseMass11(FourVector J1, FourVector J2, FourVector ME, FourVector ISR, char Variant, bool Reversal)
{
    Reversal = false;
    FourVector TempTotal = ME + J1 + J2 + ISR;
    double TempBetaZ = TempTotal[3] / TempTotal[0];
    FourVector TempJ1 = J1.Boost(FourVector(1, 0, 0, 1), TempBetaZ);
    FourVector TempJ2 = J2.Boost(FourVector(1, 0, 0, 1), TempBetaZ);
    FourVector TempME = ME.Boost(FourVector(1, 0, 0, 1), TempBetaZ);
    FourVector TempISR = ISR.Boost(FourVector(1, 0, 0, 1), TempBetaZ);
    TempTotal = TempJ1 + TempJ2 + TempME;
    double TempBeta = TempTotal.GetP() / TempTotal[0];
    TempJ1 = TempJ1.Boost(TempTotal, TempBeta);
    TempJ2 = TempJ2.Boost(TempTotal, TempBeta);
    TempME = TempME.Boost(TempTotal, TempBeta);
    double MT = 0;
    double EMET = TempME[0];
    double SumE = TempJ1[0] + TempJ2[0];
    double DeltaE = TempJ1[0] - TempJ2[0];
    double ES = (TempJ1 + TempJ2).GetP();
    double ED = (TempJ1 - TempJ2).GetP();
    double A = SumE * DeltaE / ES;
    double B = sqrt(ED * ED - A * A);
    double P1P1 = TempJ1.SpatialDot(TempJ1);
    double P1P2 = TempJ1.SpatialDot(TempJ2);
    double m0 = (SumE - EMET) / (2 * DeltaE) - SumE * DeltaE / ES / ES;
    double m1 = DeltaE * (SumE + EMET) / 2;
    double EQA = m0 * m0 + B * B / ES / ES;
    double EQB = 2 * m0 * m1 - B * B;
    double EQC = m1 * m1;
    double X2Max = (-EQB + sqrt(EQB * EQB - 4 * EQA * EQC)) / (2 * EQA);
    // double X2Min = (-EQB - sqrt(EQB * EQB - 4 * EQA * EQC)) / (2 * EQA);
    double X = sqrt(X2Max);
    double Y = ((X * X + DeltaE * DeltaE) * SumE - (X * X - DeltaE * DeltaE) * EMET) / (2 * X * DeltaE);
    double M2 = (SumE * X - DeltaE * Y) * (SumE * X - DeltaE * Y) / (X * X - DeltaE * DeltaE) / 4;
    double A0 = X / ES;
    double B0 = (Y - A * A0) / B;
    double C0 = sqrt(1 - A0 * A0 - B0 * B0);
    if (A0 * A0 + B0 * B0 >= 1)
        C0 = 0;
    FourVector DirectionA = (TempJ1 + TempJ2).SpatialNormalize();
    FourVector DirectionB = ((TempJ1 - TempJ2) - DirectionA * A).SpatialNormalize();
    FourVector DirectionC = DirectionA.SpatialCross(DirectionB);   // ambiguity in direction...
    double BetaCMSize = DeltaE / X;
    FourVector BetaCM = (DirectionA * A0 + DirectionB * B0 + DirectionC * C0).SpatialNormalize() * BetaCMSize;
    FourVector J1Boosted = TempJ1.Boost(BetaCM, -BetaCMSize);
    FourVector J2Boosted = TempJ2.Boost(BetaCM, BetaCMSize);
    if (fabs(J1Boosted.GetP() - J2Boosted.GetP()) > 0.1)  // wrong boost direction!
    {
        J1Boosted = TempJ1.Boost(BetaCM, BetaCMSize);
        J2Boosted = TempJ2.Boost(BetaCM, -BetaCMSize);
    }
    if (Variant == 'a')
        MT = GetMRT(TempJ1, TempJ2, TempME);
    if (Variant == 'b')
    {
        FourVector TransverseVector = (TempJ1 - (TempJ1 + TempJ2) * (P1P1 + P1P2) / (ES * ES)) * 2;
        MT = sqrt(M2 * 4 - TransverseVector.GetP2());
    }
    if (Variant == 'c')
    {
        double MT1 = J1Boosted.GetPT();
        double MT2 = J2Boosted.GetPT();
        MT = sqrt(M2) * 2 - (MT1 + MT2);
        J2Boosted = TempJ2.Boost(BetaCM, -BetaCMSize);
    }
    if (Variant == 'a')
        MT = GetMRT(TempJ1, TempJ2, TempME);
    if (Variant == 'b')
    {
        FourVector TransverseVector = (TempJ1 - (TempJ1 + TempJ2) * (P1P1 + P1P2) / (ES * ES)) * 2;
        MT = sqrt(M2 * 4 - TransverseVector.GetP2());
    }
    if (Variant == 'c')
    {
        double MT1 = J1Boosted.GetPT();
        double MT2 = J2Boosted.GetPT();
        MT = sqrt(M2) * 2 - (MT1 + MT2);
    }
    if (Variant == 'd')
    {
        double MT1 = J1Boosted.GetPT();
        double MT2 = J2Boosted.GetPT();
        MT = sqrt(M2) * 2 - sqrt(MT1 * MT1 + MT2 * MT2) * sqrt(2);
    }
    if (Variant == 'e')
        MT = (J1Boosted + J2Boosted).GetPT();
    if (Variant == 'f')
    {
        double MT1 = J1Boosted.GetPT();
        double MT2 = J2Boosted.GetPT();
        MT = sqrt(2 * M2 - MT1 * MT1 - MT2 * MT2);
    }
    if (Variant == 'g')
        MT = GetMRT(TempJ1, TempJ2, TempME) / BetaToGamma(BetaCMSize);
    return MT;
}

//----------------------------------------------------------------------------
LeptonVectors ConvertAnglesToVectors(const HZZ4lEventParameters &Angles, double HiggsPT, double HiggsEta)
{
    // Collect information
    double HMass = Angles.HMass;
    double ZMass = Angles.ZMass;
    double Z2Mass = Angles.Z2Mass;
    double PhiOffset = Angles.PhiH;
    double Phi0 = Angles.Phi0;
    double Theta0 = Angles.Theta0;
    double Phi = Angles.Phi;
    double Theta1 = Angles.Theta1;
    double Theta2 = Angles.Theta2;
    double Gamma1 = HMass / (2 * ZMass) * (1 + (ZMass * ZMass - Z2Mass * Z2Mass) / (HMass * HMass));
    double Gamma2 = HMass / (2 * Z2Mass) * (1 - (ZMass * ZMass - Z2Mass * Z2Mass) / (HMass * HMass));
    double Beta1 = GammaToBeta(Gamma1);
    double Beta2 = GammaToBeta(Gamma2);
    if (HMass < ZMass + Z2Mass)
        std::cerr << "[ConvertAnglesToVectors] WTF!  HMass < ZMass + Z2Mass.  Expect errors" << std::endl;
    // Calculate Higgs direction and boost in the ZZ frame (in that famous illustration used everywhere)
    FourVector BoostDirection;
    BoostDirection[0] = 1;
    BoostDirection[1] = -sin(Theta0) * sin(Phi0 - PI);
    BoostDirection[2] = sin(Theta0) * cos(Phi0 - PI);
    BoostDirection[3] = cos(Theta0);
    double HiggsBoostGamma =
        sqrt(HiggsPT * HiggsPT / HMass / HMass * cosh(HiggsEta) * cosh(HiggsEta) + 1);
    double HiggsBoostBeta = GammaToBeta(HiggsBoostGamma);
    // Write down the four-vectors in the ZZ frame
    FourVector Lepton11, Lepton12, Lepton21, Lepton22;
    Lepton11[0] = Gamma1 * (1 + Beta1 * cos(Theta1));
    Lepton11[1] = -sin(Phi) * sin(Theta1);
    Lepton11[2] = -cos(Phi) * sin(Theta1);
    Lepton11[3] = -Gamma1 * (Beta1 + cos(Theta1));
    Lepton11 = Lepton11 * ZMass / 2;
    Lepton12[0] = Gamma1 * (1 - Beta1 * cos(Theta1));
    Lepton12[1] = sin(Phi) * sin(Theta1);
    Lepton12[2] = cos(Phi) * sin(Theta1);
    Lepton12[3] = -Gamma1 * (Beta1 - cos(Theta1));
    Lepton12 = Lepton12 * ZMass / 2;
    Lepton21[0] = Gamma2 * (1 + Beta2 * cos(Theta2));
    Lepton21[1] = 0;
    Lepton21[2] = sin(Theta2);
    Lepton21[3] = Gamma2 * (Beta2 + cos(Theta2));
    Lepton21 = Lepton21 * Z2Mass / 2;
    Lepton22[0] = Gamma2 * (1 - Beta2 * cos(Theta2));
    Lepton22[1] = 0;
    Lepton22[2] = -sin(Theta2);
    Lepton22[3] = Gamma2 * (Beta2 - cos(Theta2));
    Lepton22 = Lepton22 * Z2Mass / 2;
    // Write down final higgs direction in the lab frame
    FourVector HiggsDirection;
    HiggsDirection.SetPtEtaPhiMass(HiggsPT, HiggsEta, 0, HMass);
    if (HiggsPT < 1e-5)
        HiggsDirection = FourVector(1, 0, 0, 1);
    // Calculate what is the rotation axis needed (and amount) to go from boosted-frame to lab frame
    FourVector RotationAxis = BoostDirection.SpatialCross(HiggsDirection).SpatialNormalize();
    double RotationAngle = BoostDirection.SpatialDot(HiggsDirection)
        / HiggsDirection.GetP() / BoostDirection.GetP();
    RotationAngle = acos(RotationAngle);
    // Boost to rotated lab frame
    Lepton11 = Lepton11.Boost(BoostDirection, HiggsBoostBeta);
    Lepton12 = Lepton12.Boost(BoostDirection, HiggsBoostBeta);
    Lepton21 = Lepton21.Boost(BoostDirection, HiggsBoostBeta);
    Lepton22 = Lepton22.Boost(BoostDirection, HiggsBoostBeta);
    // Rotate to lab frame
    Lepton11 = Lepton11.Rotate(RotationAxis, RotationAngle);
    Lepton12 = Lepton12.Rotate(RotationAxis, RotationAngle);
    Lepton21 = Lepton21.Rotate(RotationAxis, RotationAngle);
    Lepton22 = Lepton22.Rotate(RotationAxis, RotationAngle);
    // The extra degree of freedom that is not discussed anywhere - will become relevant for S != 0
    Lepton11 = Lepton11.Rotate(HiggsDirection, PhiOffset);
    Lepton12 = Lepton12.Rotate(HiggsDirection, PhiOffset);
    Lepton21 = Lepton21.Rotate(HiggsDirection, PhiOffset);
    Lepton22 = Lepton22.Rotate(HiggsDirection, PhiOffset);
    // Collect output
    LeptonVectors Result;
    Result.Lepton11 = Lepton11;
    Result.Lepton12 = Lepton12;
    Result.Lepton21 = Lepton21;
    Result.Lepton22 = Lepton22;
    return Result;
}

//----------------------------------------------------------------------------
HZZ4lEventParameters ConvertVectorsToAngles(const LeptonVectors &Leptons)
{
    // Some short-hand notation
    FourVector L11 = Leptons.Lepton11;
    FourVector L12 = Leptons.Lepton12;
    FourVector L21 = Leptons.Lepton21;
    FourVector L22 = Leptons.Lepton22;
    // Calculate basic information
    double HMass = (L11 + L12 + L21 + L22).GetMass();
    double ZMass = (L11 + L12).GetMass();
    double Z2Mass = (L21 + L22).GetMass();
    // Calculate higgs direction in lab frame, and boost to go to CM frame
    FourVector HiggsLab = L11 + L12 + L21 + L22;
    double HiggsBoostGamma = HiggsLab[0] / HMass;
    double HiggsBoostBeta = GammaToBeta(HiggsBoostGamma);
    // Boost everything back to the higgs rest frame
    L11 = L11.Boost(HiggsLab, -HiggsBoostBeta);
    L12 = L12.Boost(HiggsLab, -HiggsBoostBeta);
    L21 = L21.Boost(HiggsLab, -HiggsBoostBeta);
    L22 = L22.Boost(HiggsLab, -HiggsBoostBeta);
    double Theta0 = PI - GetAngle(L11 + L12, HiggsLab);
    // Z directions
    FourVector Z1 = L11 + L12;
    FourVector Z2 = L21 + L22;
    // Subtract out Z1-projection from L11, and same thing for L21 and Higgs
    // From these we can cook up the phi's
    FourVector L11Perpendicular = L11 - Z1 * (L11.SpatialDot(Z1)) / Z1.GetP2();
    FourVector L21Perpendicular = L21 - Z1 * (L21.SpatialDot(Z1)) / Z1.GetP2();
    FourVector HPerpendicular = HiggsLab - Z1 * (HiggsLab.SpatialDot(Z1)) / Z1.GetP2();
    double Phi0 = GetAngle(-HPerpendicular, L21Perpendicular);
    if (Z2.SpatialDot(HPerpendicular.SpatialCross(L21Perpendicular)) < 0)
        Phi0 = 2 * PI - Phi0;
    double Phi = GetAngle(-L21Perpendicular, L11Perpendicular);
    if (Z2.SpatialDot(L21Perpendicular.SpatialCross(L11Perpendicular)) < 0)
        Phi = 2 * PI - Phi;
    while (Phi < 0)         Phi = Phi + 2 * PI;
    while (Phi >= 2 * PI)   Phi = Phi - 2 * PI;
    // now, finally, the remaining theta's
    double Z1BoostBeta = GammaToBeta(Z1[0] / ZMass);
    double Z2BoostBeta = GammaToBeta(Z2[0] / Z2Mass);
    double Theta1 = GetAngle(Z1, L11.Boost(Z1, -Z1BoostBeta));
    double Theta2 = GetAngle(Z2, L21.Boost(Z2, -Z2BoostBeta));
    // collect output
    HZZ4lEventParameters Result;
    Result.Phi0 = Phi0;
    Result.Theta0 = Theta0;
    Result.Phi = Phi;
    Result.Theta1 = Theta1;
    Result.Theta2 = Theta2;
    Result.HMass = HMass;
    Result.ZMass = ZMass;
    Result.Z2Mass = Z2Mass;
    Result.PhiH = 0;
    return Result;
}

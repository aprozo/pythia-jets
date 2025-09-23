// main.cc
#include "Pythia8/Pythia.h"
#include <iostream>
#include <vector>
#include <cmath>

#include "fastjet/ClusterSequence.hh"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"

using namespace Pythia8;

struct DijetPair {
   int lead; // index in jet array
   int sub;

   double dphi;      // in [0, M_PI]
   double closeness; // = M_PI - dphi (smaller is better / closer to back-to-back)
};

double deltaPhi(double phi1, double phi2) // return value in (-PI, PI]
{
   double dphi = phi1 - phi2;
   while (dphi > M_PI)
      dphi -= 2 * M_PI;
   while (dphi <= -M_PI)
      dphi += 2 * M_PI;
   return dphi;
}

double deltaR(double eta1, double phi1, double eta2, double phi2)
{
   const double dphi = deltaPhi(phi1, phi2);
   const double deta = eta1 - eta2;
   return std::sqrt(deta * deta + dphi * dphi);
}

// Count charged final-state particles in a cone of radius R around (eta0, phi0)
template <class PartContainer>
int countInCone(const PartContainer &parts, double eta0, double phi0, double R, double partPtMin, double partEtaMax)
{
   if (std::abs(eta0) > partEtaMax - R)
      return 0; // require cone fully inside
   int n = 0;
   for (const auto &p : parts) {
      const double pt = p.pt();
      if (pt < partPtMin)
         continue;
      const double eta = p.eta();
      if (std::abs(eta) > partEtaMax)
         continue;
      const double phi = p.phi();
      if (deltaR(eta, phi, eta0, phi0) < R)
         ++n;
   }
   return n;
}

static std::string trim_trailing_zeros(double x)
{
   std::ostringstream os;
   os << std::fixed << std::setprecision(3) << x;
   std::string s = os.str();
   // strip trailing zeros and possibly trailing dot
   while (!s.empty() && s.back() == '0')
      s.pop_back();
   if (!s.empty() && s.back() == '.')
      s.pop_back();
   return s.empty() ? "0" : s;
}

int main(int argc, char *argv[])
{
   if (argc < 3) {
      std::cerr << "Usage: " << argv[0]
                << " pTHatMin pTHatMax|inf [nEvents=50000] [SEED=12345]"
                   " [OUTPREFIX=pp200_HardQCD]\n";
      return 1;
   }

   // Required: pTHatMin
   const double ptHatMin = std::stod(argv[1]);

   // Required: pTHatMax (can be "inf" or negative for open upper bound)
   double ptHatMax;
   std::string s = argv[2];
   if (s == "inf" || s == "Inf" || s == "INF") {
      ptHatMax = -1.0;
   } else {
      ptHatMax = std::stod(s);
   }

   if (ptHatMax > 0.0 && ptHatMax < ptHatMin) {
      std::cerr << "[error] pTHatMax < pTHatMin\n";
      return 1;
   }

   int nEvents = (argc > 3) ? std::atoi(argv[3]) : 50000;
   int seed = (argc > 4) ? std::stoi(argv[4]) : 12345;
   std::string out = (argc > 5) ? argv[5] : "pp200";

   // jet parameter
   const double jetRadius = 0.4;
   const double jetEtaMax = 1.0 - jetRadius;
   const double dPhiMin = 0.75 * M_PI; // back-to-back requirement
   const double jetPtMin = 3.0;
   // particle parameters
   const double partPtMin = 0.15;
   const double partEtaMax = 1.0;

   // Nice label for filenames
   const std::string labMin = trim_trailing_zeros(ptHatMin);
   const std::string labMax = (ptHatMax > 0.0) ? trim_trailing_zeros(ptHatMax) : "-1";
   const std::string outFile = out + "_pThat_" + labMin + "_" + labMax + ".root";

   // --- Pythia setup ---
   Pythia8::Pythia pythia8;
   pythia8.readString("Beams:idA = 2212");
   pythia8.readString("Beams:idB = 2212");
   pythia8.readString("Beams:eCM = 200.");

   pythia8.readString("HardQCD:all = on");

   // mdcy(106, 1) = 0; // PI+ 211
   // mdcy(116, 1) = 0; // K+ 321
   // mdcy(112, 1) = 0; // K_SHORT 310
   // mdcy(105, 1) = 0; // K_LONG 130
   // mdcy(164, 1) = 0; // LAMBDA0 3122
   // mdcy(162, 1) = 0; // SIGMA- 3112
   // mdcy(169, 1) = 0; // SIGMA+ 3222
   // mdcy(172, 1) = 0; // Xi- 3312
   // mdcy(174, 1) = 0; // Xi0 3322
   // mdcy(176, 1) = 0; // OMEGA- 3334
   // mdcy(102, 1) = 0; // PI0 111
   // mdcy(109, 1) = 0; // ETA 221
   // mdcy(167, 1) = 0; // SIGMA0 3212

   // pythia8.readString(
   //    "211:mayDecay = off; 321:mayDecay = off; 310:mayDecay = off; 130:mayDecay = off;"
   //    "3122:mayDecay =  off; 3112:mayDecay = off; 3222:mayDecay = off; 3312:mayDecay = off; 3322:mayDecay = off; "
   //    "3334:mayDecay = off; 111:mayDecay = off; 221:mayDecay = off; 3212:mayDecay = off");

   // Phase space cuts
   {
      std::ostringstream s1;
      s1 << "PhaseSpace:pTHatMin = " << ptHatMin;
      pythia8.readString(s1.str());
      if (ptHatMax > 0.0) {
         std::ostringstream s2;
         s2 << "PhaseSpace:pTHatMax = " << ptHatMax;
         pythia8.readString(s2.str());
      } else {
         pythia8.readString("PhaseSpace:pTHatMax = -1"); // no upper bound
      }
   }

   // Random seed
   if (seed != 0) {
      pythia8.readString("Random:setSeed = on");
      pythia8.readString(("Random:seed = " + std::to_string(seed)).c_str());
   }

   // Init
   if (!pythia8.init()) {
      std::cerr << "[error] PYTHIA init() failed.\n";
      return 2;
   }

   // --- ROOT output ---
   TFile *fout = new TFile(outFile.c_str(), "RECREATE");
   TTree *t = new TTree("events", "dijet events");

   // Jet branches (store up to 10 dijets)

   int lead_n_charged, sub_n_charged;
   double lead_pt, sub_pt, lead_eta, sub_eta, lead_phi, sub_phi, closeness, background_mult_A, background_mult_B;

   t->Branch("lead_pt", &lead_pt, "lead_pt/D");
   t->Branch("sub_pt", &sub_pt, "sub_pt/D");
   t->Branch("lead_eta", &lead_eta, "lead_eta/D");
   t->Branch("sub_eta", &sub_eta, "sub_eta/D");
   t->Branch("lead_phi", &lead_phi, "lead_phi/D");
   t->Branch("sub_phi", &sub_phi, "sub_phi/D");
   t->Branch("lead_n_charged", &lead_n_charged, "lead_n_charged/I");
   t->Branch("sub_n_charged", &sub_n_charged, "sub_n_charged/I");
   t->Branch("background_mult_A", &background_mult_A, "background_mult_A/D");
   t->Branch("background_mult_B", &background_mult_B, "background_mult_B/D");
   t->Branch("closeness", &closeness, "closeness/D");

   fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, jetRadius);

   long long accepted = 0;

   // Event loop
   for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
      if (!pythia8.next())
         continue;

      // Read pTHat
      // double pthat = pythia8.info.pTHat();

      // Build input particles for jet finding
      std::vector<fastjet::PseudoJet> parts;
      parts.reserve(2000);

      for (int i = 0; i < pythia8.event.size(); ++i) {
         const auto &p = pythia8.event[i];
         // final-state, visible (no neutrinos), basic kinematic filter
         if (!p.isFinal() || !p.isVisible())
            continue;
         if (p.idAbs() == 12 || p.idAbs() == 14 || p.idAbs() == 16)
            continue;
         if (std::abs(p.eta()) > partEtaMax)
            continue; // wide acceptance for clustering
         if (p.pT() < partPtMin)
            continue;
         fastjet::PseudoJet pj(p.px(), p.py(), p.pz(), p.e());
         pj.set_user_index(i); // <— keep Pythia index to recover charge later
         parts.push_back(pj);
      }

      // Cluster
      fastjet::ClusterSequence cs(parts, jetDef);
      fastjet::Selector select_eta = fastjet::SelectorAbsEtaMax(jetEtaMax);
      fastjet::Selector select_pt = fastjet::SelectorPtMin(jetPtMin);
      fastjet::Selector select_both = select_pt && select_eta;

      auto all_jets = fastjet::sorted_by_pt(cs.inclusive_jets());
      auto jets = select_both(all_jets);
      // Need at least two jets
      if (jets.size() != 2)
         continue;

      accepted++;

      std::vector<DijetPair> myPairs;

      for (size_t i = 0; i < jets.size(); ++i) {
         double phi1 = jets[i].phi_std();
         for (size_t j = i + 1; j < jets.size(); ++j) {
            double phi2 = jets[j].phi_std();
            double dphi12 = deltaPhi(phi1, phi2);
            dphi12 = std::abs(dphi12); // make positive

            if (dphi12 < dPhiMin)
               continue;
            // make ordered pair with leading first
            int index_lead = i, index_sub = j;
            if (jets[j].pt() > jets[i].pt())
               std::swap(index_lead, index_sub);

            DijetPair myPair;
            myPair.lead = index_lead;
            myPair.sub = index_sub;
            myPair.closeness = M_PI - dphi12;

            myPairs.push_back(myPair);
         }
      }
      // sort pairs by closeness back-to-back
      std::sort(myPairs.begin(), myPairs.end(),
                [&](const DijetPair &A, const DijetPair &B) { return A.closeness <= B.closeness; });

      std::vector<int> used(jets.size(), 0);
      std::vector<DijetPair> chosenPairs;
      chosenPairs.reserve(myPairs.size());

      for (const auto &p : myPairs) {
         if (!used[p.lead] && !used[p.sub]) {
            chosenPairs.push_back(p);
            used[p.lead] = used[p.sub] = 1;
         }
      }

      auto countCharged = [&](const fastjet::PseudoJet &j) {
         int n = 0;
         std::vector<fastjet::PseudoJet> consts = j.constituents();
         for (const auto &c : consts) {
            int idx = c.user_index();
            // Safety: user_index() is -1 if not set; skip those
            if (idx >= 0 && idx < pythia8.event.size()) {
               if (pythia8.event[idx].isCharged())
                  ++n;
            }
         }
         return n;
      };

      for (const auto &pair : chosenPairs) {

         auto leadJet = jets[pair.lead];
         auto subJet = jets[pair.sub];

         lead_n_charged = countCharged(leadJet);
         sub_n_charged = countCharged(subJet);

         lead_pt = leadJet.pt();
         sub_pt = subJet.pt();
         lead_eta = leadJet.eta();
         sub_eta = subJet.eta();
         lead_phi = leadJet.phi_std();
         sub_phi = subJet.phi_std();
         closeness = pair.closeness;

         double phiA = deltaPhi(lead_phi, M_PI / 2);
         double phiB = deltaPhi(lead_phi, -M_PI / 2);

         background_mult_A = countInCone(parts, lead_eta, phiB, jetRadius, partPtMin, partEtaMax);
         background_mult_B = countInCone(parts, lead_eta, phiA, jetRadius, partPtMin, partEtaMax);

         t->Fill();
      }
   }

   // Cross sections (mb)
   const double sigmaGen = pythia8.info.sigmaGen();
   const double sigmaErr = pythia8.info.sigmaErr();

   TH1D *stats = new TH1D("stats", "stats", 6, 0, 6);
   vector<TString> statNames = {"nEvents", "nAccepted", "ptHatMin", "ptHatMax", "sigmaGen_mb", "sigmaErr_mb"};
   for (size_t i = 0; i < statNames.size(); ++i)
      stats->GetXaxis()->SetBinLabel(i + 1, statNames[i]);

   stats->SetBinContent(1, nEvents);
   stats->SetBinContent(2, accepted);
   // stats->SetBinContent(3, ptHatMin);
   // stats->SetBinContent(4, ptHatMax);
   stats->SetBinContent(5, sigmaGen);
   stats->SetBinError(5, sigmaErr);
   // stats->SetBinContent(6, sigmaErr);

   // Print and record
   std::cout << "[done] Wrote " << outFile << "\n"
             << "       N_accepted = " << accepted << "\n"
             << "       sigmaGen   = " << sigmaGen << " mb  (± " << sigmaErr << ")\n";

   pythia8.stat();

   fout->Write();
   fout->Close();
   delete fout;

   std::cout << "Accepted dijet-like events: " << accepted << " / " << nEvents << std::endl;
   return 0;
}

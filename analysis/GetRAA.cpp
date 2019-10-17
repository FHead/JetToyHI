#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

#include "CommandLine.h"
#include "PlotHelper4.h"
#include "SetStyle.h"

void PrintYiFile(string OutFileName, TH1D histogram[], int NRadii);

int main(int argc, char *argv[])
{
   SetThesisStyle();

   CommandLine CL(argc, argv);

   // Input one file for PP, and one for each PbPb type.
   vector<string> InputFileNames = CL.GetStringVector("input", ' ');
   vector<double> bins = CL.GetDoubleVector("bins", "200,250,300,400,500,1000");
   string OutputBase = CL.Get("output");
   string tag = CL.Get("tag");   // Should be "Jewel" or "Pyquen".

   PdfFileHelper PdfFile(OutputBase + ".pdf");
   PdfFile.AddTextPage(tag + " Spectra & RAA");

   TFile *outf = new TFile(Form("%s.root",OutputBase.c_str()), "recreate");

   const int NRadii = 6;
   string radius[NRadii] = {"02", "03", "04", "06", "08", "10"};
   
   const int NMaxType = 4;
   const int NMaxCentrality = 4;
   string CentralityTag[NMaxCentrality] = {"0to10", "10to30", "30to50", "50to90"};  // Used to name histograms.
   string centrality[NMaxCentrality] = {"0-10", "10-30", "30-50", "50-90"};         // Used to open files.

   int NCentrality;
   
   int NType;
   vector<string> type;
   
   // Define the different kinds of PbPb data depending on the theory type.
   if(tag == "Jewel")
   {
      NCentrality = 4;
      NType = 2;
      type.push_back("PbPbNoRecoil");
      type.push_back("PbPb");
   }
   else if(tag == "Pyquen")
   {
      NCentrality = 2;
      centrality[0] = "";
      centrality[1] = "0-10";
      CentralityTag[0] = "b=0";
      CentralityTag[1] = "0to10";
      NType = 2;
      type.push_back("PbPbWide");
      type.push_back("PbPb");
   }
   else
   {
      cout << "Input Error: Type should be \"Jewel\" or \"Pyquen\"" << endl;
      return -1;
   }

   // Create histograms.
   TH1D PPJetPT[NRadii];
   TH1D PbPbJetPT[NMaxType][NMaxCentrality][NRadii];
   TH1D JetRAA[NMaxType][NMaxCentrality][NRadii];
   TH1D JetRRAA[NMaxType][NMaxCentrality][NRadii];

   const int NBins = 160;
   
   for(int iR = 0; iR < NRadii; iR++)
   {
      PPJetPT[iR] = TH1D(Form("PP_R%s_JetPT",radius[iR].c_str()), Form("PP R = %s Weighted & Scaled Jet PT",radius[iR].c_str()), NBins, 200, 1000);
      PPJetPT[iR].Sumw2();

      for(int iT = 0; iT < NType; iT++)
      {
         for(int iC = 0; iC < NCentrality; iC++)
         {
            string HistogramName;
            string HistogramTitle;
            
            if(tag == "Jewel")
            {
               HistogramName = type[iT] + "_R" + radius[iR] + "_C" + CentralityTag[iC];
               HistogramTitle = type[iT] + " " + centrality[iC] + "% R = " + radius[iR];
            }
            else if(tag == "Pyquen")
            {
               HistogramName = type[iT] + "_R" + radius[iR] + "_C" + CentralityTag[iC];
               HistogramTitle = type[iT] + " " + centrality[iC] + " R = " + radius[iR];
            }
            
            PbPbJetPT[iT][iC][iR] = TH1D(Form("%s_JetPT",HistogramName.c_str()), Form("%s Weighted & Scaled Jet PT",HistogramTitle.c_str()), NBins, 200, 1000);
            PbPbJetPT[iT][iC][iR].Sumw2();

            JetRAA[iT][iC][iR] = TH1D(Form("%s_JetRAA",HistogramName.c_str()), Form("%s Jet RAA",HistogramTitle.c_str()), bins.size()-1, &bins[0]);
            JetRAA[iT][iC][iR].SetMaximum(1.5);
            JetRAA[iT][iC][iR].SetMinimum(0);
            JetRAA[iT][iC][iR].SetStats(0);
            
            JetRRAA[iT][iC][iR] = TH1D(Form("%s_JetRRAA",HistogramName.c_str()), Form("%s Jet RRAA",HistogramTitle.c_str()), bins.size()-1, &bins[0]);
            JetRRAA[iT][iC][iR].SetMaximum(2.5);
            JetRRAA[iT][iC][iR].SetMinimum(0);
            JetRRAA[iT][iC][iR].SetStats(0);
         }
      }
   }

   // Get spectra.
   for(string FileName : InputFileNames)
   {
      if (InputFileNames.size() != NType*NCentrality + 1)
      {
         cout << "Input Error: Please input " << NType*NCentrality + 1 << " files for PP";
         for (int iT = 0; iT < NType; iT++)
         {
            if (iT == NType - 1) cout << ", and ";
            else                 cout << ", ";
            cout << type[iT];
         }
      
         if(tag == "Jewel") cout << " for centrality 0-10%, 10-30%, 30-50%, and 50-90%" << endl;
         else if(tag == "Pyquen") cout << " for b = 0 and centrality 0-10%" << endl;
         
         return -1;
      }

      TFile File(FileName.c_str());
      cout << FileName.c_str() << endl;
      TTree *Tree = (TTree *)File.Get("JetTree");
      
      if(Tree == nullptr)
      {
         File.Close();
         cout << "Input Error: Null tree in file " << FileName << endl;
         return -1;
      }

      vector<double> *SignalJetPt[NRadii];
      vector<double> *EventWeight = nullptr;
      vector<double> *EventWeightR03 = nullptr;

      for(int iR = 0; iR < NRadii; iR++)
      {
         SignalJetPt[iR] = nullptr;
         Tree->SetBranchAddress(Form("SignalJet%s%sPt",radius[iR].c_str(),tag.c_str()), &SignalJetPt[iR]);
      }

      Tree->SetBranchAddress("EventWeight", &EventWeight);
      
      // Get a pointer to fill the correct histogram.
      // R = 3 for Jewel PP and PbPb (aka not PbPbNoRecoil) needs a different EventWeight because it was generated separately.
      bool R3Correction = false;
      TH1D *JetPTPtr;
      if(FileName.find("PP.root") != string::npos)
      {
         JetPTPtr = PPJetPT;
         if(tag == "Jewel") // R = 3 bug.
         {
            R3Correction = true;
            Tree->SetBranchAddress("EventWeightR03", &EventWeightR03);
         }
      }
      else for(int iT = 0; iT < NType; iT++)
      {
         for(int iC = 0; iC < NCentrality; iC++)
         {
            string CheckFileName;
            
            if(tag == "Jewel")         CheckFileName = type[iT] + "-" + centrality[iC] + ".root";
            else if (tag == "Pyquen")
            {
               if (centrality[iC] == "0-10")
                  CheckFileName = type[iT] + "-" + centrality[iC] + ".root";
               else
                  CheckFileName = type[iT] + ".root";
            }
            
            if (FileName.find(CheckFileName) != string::npos)
            {
               JetPTPtr = PbPbJetPT[iT][iC];
               
               if(tag == "Jewel" && type[iT] == "PbPb" && (centrality[iC] == "0-10" || centrality[iC] == "10-30")) // R = 3 bug.
               {
                  R3Correction = true;
                  Tree->SetBranchAddress("EventWeightR03", &EventWeightR03);
               }
            }
         }
      }

      // Fill PT histogram.
      int EntryCount = Tree->GetEntries();
      for(int iE = 0; iE < EntryCount; iE++)
      {
         Tree->GetEntry(iE);

         for(int iR = 0; iR < NRadii; iR++)
         {
            if(SignalJetPt[iR] == nullptr)
            {
               cout << "Input Error: Null branch in file " << FileName << " for radius " << radius[iR] << endl;
               return -1;
            }

            int NJet = SignalJetPt[iR]->size();
            
            // Need different event weights for R = 3 due to bug when generating.
            if (R3Correction & radius[iR] == "03")
            {
               for(int iJ = 0; iJ < NJet; iJ++)
                  JetPTPtr[iR].Fill((*SignalJetPt[iR])[iJ], (*EventWeightR03)[0]);
            }
            else
            {
               for(int iJ = 0; iJ < NJet; iJ++)
                  JetPTPtr[iR].Fill((*SignalJetPt[iR])[iJ], (*EventWeight)[0]);
            }
         }
      }

      for(int iR = 0; iR < NRadii; iR++)
         JetPTPtr[iR].Scale(1.0/(EntryCount));

      File.Close();
   }
   
   // Rebin spectra.
   TH1 *PPJetPTRebin[NRadii];
   TH1 *PbPbJetPTRebin[NMaxType][NMaxCentrality][NRadii];
   for(int iR = 0; iR < NRadii; iR++)
   {
      PPJetPTRebin[iR] = PPJetPT[iR].Rebin(bins.size()-1,Form("PP_%s_JetPT_Rebin",radius[iR].c_str()),&bins[0]);

      for(int iT = 0; iT < NType; iT++)
      {
         for(int iC = 0; iC < NCentrality; iC++)
         {
            string HistogramName;
            if(tag == "Jewel")         HistogramName = type[iT] + "_R" + radius[iR] + "_C" + CentralityTag[iC];
            else if(tag == "Pyquen")   HistogramName = type[iT] + "_R" + radius[iR];
            
            PbPbJetPTRebin[iT][iC][iR] = PbPbJetPT[iT][iC][iR].Rebin(bins.size()-1,Form("%s_JetRAA_Rebin",HistogramName.c_str(),type[iT].c_str()),&bins[0]);
         }
      }
   }
   
   // Produce jet RAA histograms.
   for(int iR = 0; iR < NRadii; iR++)
   {
      for(int iT = 0; iT < NType; iT++)
      {
         for(int iC = 0; iC < NCentrality; iC++)
         {
            for(int iB = 1; iB <= bins.size()-1; iB++)
            {
               double PbPbValue = PbPbJetPTRebin[iT][iC][iR]->GetBinContent(iB);
               double PbPbRelError = PbPbJetPTRebin[iT][iC][iR]->GetBinError(iB)/PbPbValue;
               double PPValue = PPJetPTRebin[iR]->GetBinContent(iB);
               double PPRelError = PPJetPTRebin[iR]->GetBinError(iB)/PPValue;
               double RAAValue = PbPbValue/PPValue;
               double RAAError = sqrt(pow(PbPbRelError,2) + pow(PPRelError,2))*RAAValue;

               JetRAA[iT][iC][iR].SetBinContent(iB,RAAValue);
               JetRAA[iT][iC][iR].SetBinError(iB,RAAError);
            }
         }
      }
   }
   
   // Produce RRAA histograms.
   for(int iR = 0; iR < NRadii; iR++) // Assuming iR = 0 is R = 02.
   {
      for(int iT = 0; iT < NType; iT++)
      {
         for(int iC = 0; iC < NCentrality; iC++)
         {
            for(int iB = 1; iB <= bins.size()-1; iB++)
            {
               if(iR == 0)
               {
                  JetRRAA[iT][iC][0].SetBinContent(iB,1);
                  JetRRAA[iT][iC][0].SetBinError(iB,0);
               }
               else
               {
                  double ThisRAAValue = JetRAA[iT][iC][iR].GetBinContent(iB);
                  double ThisRAARelError = JetRAA[iT][iC][iR].GetBinError(iB)/ThisRAAValue;
                  double R02RAAValue = JetRAA[iT][iC][0].GetBinContent(iB);
                  double R02RAARelError = JetRAA[iT][iC][0].GetBinError(iB)/R02RAAValue;
                  double RRAAValue = ThisRAAValue/R02RAAValue;
                  double RRAAError = sqrt(pow(ThisRAARelError,2) + pow(R02RAARelError,2))*RRAAValue;

                  JetRRAA[iT][iC][iR].SetBinContent(iB,RRAAValue);
                  JetRRAA[iT][iC][iR].SetBinError(iB,RRAAError);
               }
            }
         }
      }
   }
   
   // Write histograms to pdf and root files.
   for(int iR = 0; iR < NRadii; iR++)
   {
      PdfFile.AddPlot(PPJetPT[iR], "", true);
      for(int iT = 0; iT < NType; iT++)
         for(int iC = 0; iC < NCentrality; iC++)
            PdfFile.AddPlot(PbPbJetPT[iT][iC][iR], "", true);
   }
   
   for(int iT = 0; iT < NType; iT++)
      for(int iR = 0; iR < NRadii; iR++)
         for(int iC = 0; iC < NCentrality; iC++)
            PdfFile.AddPlot(JetRAA[iT][iC][iR], "", false);

   for(int iT = 0; iT < NType; iT++)
      for(int iR = 0; iR < NRadii; iR++)
         for(int iC = 0; iC < NCentrality; iC++)
            PdfFile.AddPlot(JetRRAA[iT][iC][iR], "", false);
   
   //PdfFile.AddTimeStampPage();
   PdfFile.Close();

   outf->cd();
   for(int iR = 0; iR < NRadii; iR++)
   {
      PPJetPT[iR].Write();
      for(int iT = 0; iT < NType; iT++)
         for(int iC = 0; iC < NCentrality; iC++)
            PbPbJetPT[iT][iC][iR].Write();
   }
   
   for(int iT = 0; iT < NType; iT++)
      for(int iR = 0; iR < NRadii; iR++)
         for(int iC = 0; iC < NCentrality; iC++)
            JetRAA[iT][iC][iR].Write();
   
   for(int iT = 0; iT < NType; iT++)
      for(int iR = 0; iR < NRadii; iR++)
         for(int iC = 0; iC < NCentrality; iC++)
            JetRRAA[iT][iC][iR].Write();

   // Write to text file for Yi's analysis code
   for(int iT = 0; iT < NType; iT++)
   {
      for (int iC = 0; iC < NCentrality; iC++)
      {
         string OutFileBase;
         if(tag == "Jewel") OutFileBase = OutputBase + "_" + type[iT] + "_" + CentralityTag[iC];
         else if(tag == "Pyquen")
         {
            if (centrality[iC] == "0-10")
               OutFileBase = OutputBase + "_" + type[iT] + "_" + CentralityTag[iC];
            else
               OutFileBase = OutputBase + "_" + type[iT];
         }
         
         PrintYiFile(OutFileBase + "_RAA.txt", JetRAA[iT][iC], NRadii);
         PrintYiFile(OutFileBase + "_RRAA.txt", JetRRAA[iT][iC], NRadii);
      }
   }

   return 0;
}

// All histograms must have the same binning.
void PrintYiFile(string OutFileName, TH1D histogram[], int NRadii)
{
   ofstream OutFile;
   OutFile.open(OutFileName.c_str());
   
   TAxis *axis = histogram[0].GetXaxis();
   
   for(int iB = 1; iB <= axis->GetNbins(); iB++)
   {
      OutFile << axis->GetBinLowEdge(iB) << " ";
      OutFile << axis->GetBinUpEdge(iB) << " ";
      
      for(int iR = 0; iR < NRadii; iR++)
      {
         double BinContent = histogram[iR].GetBinContent(iB);
         double BinError = histogram[iR].GetBinError(iB);
         OutFile << BinContent << " ";
         OutFile << BinContent + BinError << " ";
         OutFile << BinContent - BinError << " ";
      }
         
      OutFile << endl;
   }
   
   OutFile.close();
}

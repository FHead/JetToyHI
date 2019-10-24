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
   
   int NCentrality;
   vector<string> CentralityTag; // Used to name histograms.
   vector<string> centrality;    // Used to open files.
   
   int NType; // Holds different types of PbPb data produced by the same generator with different settings.
   vector<string> type;
   
   int NStat;  // For Pyquen which has a normal pthat sample and a pthat = 250 sample.
   vector<string> stat;
   
   // Define the different kinds of PbPb data depending on the theory type.
   if(tag == "Jewel")
   {
      NCentrality = 4;
      CentralityTag.push_back("0to10");
      CentralityTag.push_back("10to30");
      CentralityTag.push_back("30to50");
      CentralityTag.push_back("50to90");
      centrality.push_back("0-10");
      centrality.push_back("10-30");
      centrality.push_back("30-50");
      centrality.push_back("50-90");
      NType = 2;
      type.push_back("PbPbNoRecoil");
      type.push_back("PbPb");
      NStat = 1;
      stat.push_back("");
   }
   else if(tag == "Pyquen")
   {
      NCentrality = 1;
      CentralityTag.push_back("0to10");
      centrality.push_back("0-10");
      NType = 2;
      type.push_back("PbPbWide");
      type.push_back("PbPb");
      NStat = 2;
      stat.push_back("150");
      stat.push_back("250");
   }
   else
   {
      cout << "Input Error: Type should be \"Jewel\" or \"Pyquen\"" << endl;
      return -1;
   }

   // Create histograms.
   TH1D PPJetPT[NStat][NRadii];
   TH1D PbPbJetPT[NType][NStat][NCentrality][NRadii];
   TH1D JetRAA[NType][NStat][NCentrality][NRadii];
   TH1D JetRRAA[NType][NStat][NCentrality][NRadii];

   const int NBins = 160;
   
   for(int iS = 0; iS < NStat; iS++)
   {
      for(int iR = 0; iR < NRadii; iR++)
      {
         PPJetPT[iS][iR] = TH1D(Form("PP%s_R%s_JetPT",stat[iS].c_str(),radius[iR].c_str()), Form("PP R = %s Weighted & Scaled Jet PT",radius[iR].c_str()), NBins, 200, 1000);
         PPJetPT[iS][iR].Sumw2();

         for(int iT = 0; iT < NType; iT++)
         {
            for(int iC = 0; iC < NCentrality; iC++)
            {
               string HistogramName;
               string HistogramTitle;
               
               if(tag == "Jewel")
               {
                  HistogramName = type[iT] + stat[iS] + "_R" + radius[iR] + "_C" + CentralityTag[iC];
                  HistogramTitle = type[iT] + stat[iS] + " " + centrality[iC] + "% R = " + radius[iR];
               }
               else if(tag == "Pyquen")
               {
                  HistogramName = type[iT] + stat[iS] + "_R" + radius[iR] + "_C" + CentralityTag[iC];
                  HistogramTitle = type[iT] + stat[iS] + " " + centrality[iC] + " R = " + radius[iR];
               }
               
               PbPbJetPT[iS][iT][iC][iR] = TH1D(Form("%s_JetPT",HistogramName.c_str()), Form("%s Weighted & Scaled Jet PT",HistogramTitle.c_str()), NBins, 200, 1000);
               PbPbJetPT[iS][iT][iC][iR].Sumw2();

               JetRAA[iS][iT][iC][iR] = TH1D(Form("%s_JetRAA",HistogramName.c_str()), Form("%s Jet RAA",HistogramTitle.c_str()), bins.size()-1, &bins[0]);
               JetRAA[iS][iT][iC][iR].SetMaximum(1.5);
               JetRAA[iS][iT][iC][iR].SetMinimum(0);
               JetRAA[iS][iT][iC][iR].SetStats(0);
               
               JetRRAA[iS][iT][iC][iR] = TH1D(Form("%s_JetRRAA",HistogramName.c_str()), Form("%s Jet RRAA",HistogramTitle.c_str()), bins.size()-1, &bins[0]);
               JetRRAA[iS][iT][iC][iR].SetMaximum(2.5);
               JetRRAA[iS][iT][iC][iR].SetMinimum(0);
               JetRRAA[iS][iT][iC][iR].SetStats(0);
            }
         }
      }
   }

   // Get spectra.
   for(string FileName : InputFileNames)
   {
      if(tag == "Jewel" && InputFileNames.size() != 9)
      {
         cout << "Input Error: Please input 9 files for PP, PbPb, and PbPbNoRecoil centrality 0-10%, 10-30%, 30-50%, and 50-90%" << endl;
         return -1;
      }
      if(tag == "Pyquen" && InputFileNames.size() != 6)
      {
         cout << "Input Error: Please input 6 files for PP, PbPb, and PbPbNoRecoil for normal pthat and pthat = 250" << endl;
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
      
      for(int iS = 0; iS < NStat; iS++)
      {
         if(FileName.find("PP" + stat[iS] + ".root") != string::npos)
         {
            JetPTPtr = PPJetPT[iS];

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
               string CheckFileName = type[iT] + stat[iS] + "-" + centrality[iC] + ".root";
               
               if (FileName.find(CheckFileName) != string::npos)
               {
                  JetPTPtr = PbPbJetPT[iS][iT][iC];
                  
                  if(tag == "Jewel" && type[iT] == "PbPb" && (centrality[iC] == "0-10" || centrality[iC] == "10-30")) // R = 3 bug.
                  {
                     R3Correction = true;
                     Tree->SetBranchAddress("EventWeightR03", &EventWeightR03);
                  }
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
   TH1 *PPJetPTRebin[NStat][NRadii];
   TH1 *PbPbJetPTRebin[NStat][NType][NCentrality][NRadii];
   for(int iS = 0; iS < NStat; iS++)
   {
      for(int iR = 0; iR < NRadii; iR++)
      {
         PPJetPTRebin[iS][iR] = PPJetPT[iS][iR].Rebin(bins.size()-1,Form("PP%s_%s_JetPT_Rebin",stat[iS].c_str(),radius[iR].c_str()),&bins[0]);

         for(int iT = 0; iT < NType; iT++)
         {
            for(int iC = 0; iC < NCentrality; iC++)
            {
               string HistogramName = type[iT] + stat[iS] + "_R" + radius[iR] + "_C" + CentralityTag[iC];
               
               PbPbJetPTRebin[iS][iT][iC][iR] = PbPbJetPT[iS][iT][iC][iR].Rebin(bins.size()-1,Form("%s_JetRAA_Rebin",HistogramName.c_str(),type[iT].c_str()),&bins[0]);
            }
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
               double PbPbValue = PbPbJetPTRebin[0][iT][iC][iR]->GetBinContent(iB);
               double PbPbRelError = PbPbJetPTRebin[0][iT][iC][iR]->GetBinError(iB)/PbPbValue;
               double PPValue = PPJetPTRebin[0][iR]->GetBinContent(iB);
               double PPRelError = PPJetPTRebin[0][iR]->GetBinError(iB)/PPValue;
               double RAAValue = PbPbValue/PPValue;
               double RAAError = sqrt(pow(PbPbRelError,2) + pow(PPRelError,2))*RAAValue;

               JetRAA[0][iT][iC][iR].SetBinContent(iB,RAAValue);
               JetRAA[0][iT][iC][iR].SetBinError(iB,RAAError);
               
               if(tag == "Pyquen")
               {
                  if(bins[iB-1] >= 300) // Only want to use pthat = 250 for bins greater than 300.
                  {
                     PbPbValue = PbPbJetPTRebin[1][iT][iC][iR]->GetBinContent(iB);
                     PbPbRelError = PbPbJetPTRebin[1][iT][iC][iR]->GetBinError(iB)/PbPbValue;
                     PPValue = PPJetPTRebin[1][iR]->GetBinContent(iB);
                     PPRelError = PPJetPTRebin[1][iR]->GetBinError(iB)/PPValue;
                     RAAValue = PbPbValue/PPValue;
                     RAAError = sqrt(pow(PbPbRelError,2) + pow(PPRelError,2))*RAAValue;
                  }
                  JetRAA[1][iT][iC][iR].SetBinContent(iB,RAAValue);
                  JetRAA[1][iT][iC][iR].SetBinError(iB,RAAError);
               }
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
                  for(int iS = 0; iS < NStat; iS++)
                  {
                     JetRRAA[iS][iT][iC][0].SetBinContent(iB,1);
                     JetRRAA[iS][iT][iC][0].SetBinError(iB,0);
                  }
               }
               else
               {
                  double ThisRAAValue = JetRAA[0][iT][iC][iR].GetBinContent(iB);
                  double ThisRAARelError = JetRAA[0][iT][iC][iR].GetBinError(iB)/ThisRAAValue;
                  double R02RAAValue = JetRAA[0][iT][iC][0].GetBinContent(iB);
                  double R02RAARelError = JetRAA[0][iT][iC][0].GetBinError(iB)/R02RAAValue;
                  double RRAAValue = ThisRAAValue/R02RAAValue;
                  double RRAAError = sqrt(pow(ThisRAARelError,2) + pow(R02RAARelError,2))*RRAAValue;

                  JetRRAA[0][iT][iC][iR].SetBinContent(iB,RRAAValue);
                  JetRRAA[0][iT][iC][iR].SetBinError(iB,RRAAError);
                  
                  if(tag == "Pyquen")
                  {
                     if(bins[iB-1] >= 300) // Only want to use pthat = 250 for bins greater than 300.
                     {
                        ThisRAAValue = JetRAA[1][iT][iC][iR].GetBinContent(iB);
                        ThisRAARelError = JetRAA[1][iT][iC][iR].GetBinError(iB)/ThisRAAValue;
                        R02RAAValue = JetRAA[1][iT][iC][0].GetBinContent(iB);
                        R02RAARelError = JetRAA[1][iT][iC][0].GetBinError(iB)/R02RAAValue;
                        RRAAValue = ThisRAAValue/R02RAAValue;
                        RRAAError = sqrt(pow(ThisRAARelError,2) + pow(R02RAARelError,2))*RRAAValue;
                     }
                     JetRRAA[1][iT][iC][iR].SetBinContent(iB,RRAAValue);
                     JetRRAA[1][iT][iC][iR].SetBinError(iB,RRAAError);
                  }
               }
            }
         }
      }
   }
   
   // Write histograms to pdf and root files.
//   for(int iR = 0; iR < NRadii; iR++)
//   {
//      PdfFile.AddPlot(PPJetPT[iR], "", true);
//      for(int iT = 0; iT < NType; iT++)
//         for(int iC = 0; iC < NCentrality; iC++)
//            PdfFile.AddPlot(PbPbJetPT[iT][iC][iR], "", true);
//   }
//
//   for(int iT = 0; iT < NType; iT++)
//      for(int iR = 0; iR < NRadii; iR++)
//         for(int iC = 0; iC < NCentrality; iC++)
//            PdfFile.AddPlot(JetRAA[iT][iC][iR], "", false);
//
//   for(int iT = 0; iT < NType; iT++)
//      for(int iR = 0; iR < NRadii; iR++)
//         for(int iC = 0; iC < NCentrality; iC++)
//            PdfFile.AddPlot(JetRRAA[iT][iC][iR], "", false);
//
//   //PdfFile.AddTimeStampPage();
//   PdfFile.Close();

   outf->cd();
   for(int iS = 0; iS < NStat; iS++)
      for(int iR = 0; iR < NRadii; iR++)
      {
         PPJetPT[iS][iR].Write();
            for(int iT = 0; iT < NType; iT++)
               for(int iC = 0; iC < NCentrality; iC++)
                  PbPbJetPT[iS][iT][iC][iR].Write();
      }
   
   for(int iS = 0; iS < NStat; iS++)
      for(int iT = 0; iT < NType; iT++)
         for(int iR = 0; iR < NRadii; iR++)
            for(int iC = 0; iC < NCentrality; iC++)
               JetRAA[iS][iT][iC][iR].Write();
   
   for(int iS = 0; iS < NStat; iS++)
      for(int iT = 0; iT < NType; iT++)
         for(int iR = 0; iR < NRadii; iR++)
            for(int iC = 0; iC < NCentrality; iC++)
               JetRRAA[iS][iT][iC][iR].Write();

   // Write to text file for Yi's analysis code
   for(int iS = 0; iS < NStat; iS++)
   {
      for(int iT = 0; iT < NType; iT++)
      {
         for (int iC = 0; iC < NCentrality; iC++)
         {
            string OutFileBase = OutputBase + "_" + type[iT] + stat[iS] + "_" + CentralityTag[iC];
            
            PrintYiFile(OutFileBase + "_RAA.txt", JetRAA[iS][iT][iC], NRadii);
            PrintYiFile(OutFileBase + "_RRAA.txt", JetRRAA[iS][iT][iC], NRadii);
         }
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

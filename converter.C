class RhoCandList;
class RhoCandidate;
class PndAnaPidSelector;
class PndAnaPidCombiner;
class PndAnalysis;




void converter(int nevts = 0, TString prefix = "signal")
{
 	// *** some variables
	int i=0;
	gStyle->SetOptFit(1011);

	
	// *** the output file for FairRunAna
	TString OutFile="out_dummy.root";  
	
	// output 'CSV' file
	ofstream myfile;
	myfile.open ("outputPions.csv");   // name it	
					
	// *** the files coming from the simulation
	TString inPidFile  = prefix+"_pid.root";    // this file contains the PndPidCandidates and McTruth
	TString inParFile  = prefix+"_par.root";
	
	// *** PID table with selection thresholds; can be modified by the user
	TString pidParFile = TString(gSystem->Getenv("VMCWORKDIR"))+"/macro/params/all.par";	
	
	// *** initialization
	FairLogger::GetLogger()->SetLogToFile(kFALSE);
	FairRunAna* fRun = new FairRunAna();
	FairRuntimeDb* rtdb = fRun->GetRuntimeDb();
	fRun->SetSource(new FairFileSource(inPidFile));
	
	// *** setup parameter database 	
	FairParRootFileIo* parIO = new FairParRootFileIo();
	parIO->open(inParFile);
	FairParAsciiFileIo* parIOPid = new FairParAsciiFileIo();
	parIOPid->open(pidParFile.Data(),"in");
	
	rtdb->setFirstInput(parIO);
	rtdb->setSecondInput(parIOPid);
	rtdb->setOutput(parIO);  
	
	fRun->SetOutputFile(OutFile);
	fRun->Init(); 

	// *** the data reader object
	PndAnalysis* theAnalysis = new PndAnalysis();
	if (nevts==0) nevts= theAnalysis->GetEntries();
	
	// *** RhoCandLists for the analysis
	RhoCandList eminus, piminus, muminus, kminus;
	
//  ====================================================== Variables Begin ====================================================== //	
	// Basic Variables
	vector<double_t> momentumx, momentumy, momentumz, momentum, energy;
	vector<float_t>  positionx, positiony, positionz, position;
	vector<float_t>  fhitx, fhity, fhitz, lhitx, lhity,lhitz, fhit, lhit;
	vector<int>      charge, MCindex, Trackindex;
	
	// Detector Variables
	
    // == MVD ==
    vector<float_t> MvdDEDX;
    vector<int>     MvdHits;
    
    // == STT ==
    vector<float_t>  SttMeanDEDX;
    vector<int>      SttHits;
    
    // == GEM ==
    vector<int>     GemHits;
    
    // == TOF == 
    vector<float_t>  TofStopTime;
    vector<float_t>  TofM2;
    vector<float_t>  TofTrackLength;
    vector<float_t>  TofQuality;
    vector<int>      TofIndex; 
    vector<float>    TofBeta; 
    

    // == Barrel DIRC ==
    vector<float_t>  DrcThetaC;
    vector<float_t>  DrcThetaCErr;
    vector<float_t>  DrcQuality;
    vector<int>      DrcNumberOfPhotons;
    vector<int>      DrcIndex;    
    
    // == Disc DIRC ==
    vector<float_t>  DiscThetaC;
    vector<float_t>  DiscThetaCErr;
    vector<float_t>  DiscQuality;
    vector<int>      DiscNumberOfPhotons;
    vector<int>      DiscIndex;
    
    // == RICH ==
    vector<float_t>  RichThetaC;
    vector<float_t>  RichThetaCErr;
    vector<float_t>  RichQuality;
    vector<int>      RichNumberOfPhotons;
    vector<int>      RichIndex;
    
    // == EMC ==
    vector<float_t>  EmcRawEnergy;
    vector<float_t>  EmcCalEnergy;
    vector<float_t>  EmcQuality;
    vector<int>      EmcNumberOfCrystals;
    vector<int>      EmcNumberOfBumps;
    vector<int>      EmcModule;
    vector<int>      EmcIndex;
    
       // EMC Cluster properties
       vector<double_t>   EmcZ20;
       vector<double_t>   EmcZ53;
       vector<double_t>   EmcLat;
       vector<double_t>   EmcE1;
       vector<double_t>   EmcE9;
       vector<double_t>   EmcE25;
    
    // == Muons ==
	vector<float_t>  MuoProbability;
	vector<float_t>  MuoQuality;
	vector<float_t>  MuoIron;
	vector<float_t>  MuoMomentumIn;
	vector<int>      MuoNumberOfLayers;
	vector<int>      MuoModule;
	vector<int>      MuoHits;
	vector<int>      MuoIndex;
	
	// == Tracking ==
	vector<int>      DegreesOfFreedom;
	vector<int>      FitStatus;
	vector<float_t>  ChiSquared;	
//  ====================================================== Variables end   ====================================================== //	
	
	// ***
	// the event loop
	// ***
	int I=0;
	while (theAnalysis->GetEvent() && i++<nevts)
	{
		if ((i%100)==0) cout<<"evt " << i << endl;
		
		// *** Select with no PID info ('All'); type and mass are set
		theAnalysis->FillList( eminus, "ElectronAllMinus");
		theAnalysis->FillList( piminus, "PionAllMinus");
		theAnalysis->FillList( muminus, "MuonAllMinus");
		theAnalysis->FillList( kminus, "KaonAllMinus");
		    
		for (int j=0; j<piminus.GetLength(); ++j)
		{
			if (theAnalysis->McTruthMatch(piminus[j]))
			{
				PndPidCandidate *myCand = (PndPidCandidate*) piminus[j]->GetRecoCandidate();

		        momentumx.push_back(myCand->GetMomentum().X());
		        momentumy.push_back(myCand->GetMomentum().Y());
		        momentumz.push_back(myCand->GetMomentum().Z());
		        
		        double mom = ((myCand->GetMomentum().X()*myCand->GetMomentum().X()) + (myCand->GetMomentum().Y()*myCand->GetMomentum().Y()) + (myCand->GetMomentum().Z()*myCand->GetMomentum().Z()));		        
		        momentum.push_back(mom);
		        
		        energy.push_back(myCand->GetEnergy());
		        fhitx.push_back(myCand->GetFirstHit().X());
		        fhity.push_back(myCand->GetFirstHit().Y());
		        fhitz.push_back(myCand->GetFirstHit().Z());
		        lhitx.push_back(myCand->GetLastHit().X());
		        lhity.push_back(myCand->GetLastHit().Y());
		        lhitz.push_back(myCand->GetLastHit().Z());
		        positionx.push_back(myCand->GetPosition().X());	
		        positiony.push_back(myCand->GetPosition().Y());
		        positionz.push_back(myCand->GetPosition().Z());	
		        
		        float pos = ((myCand->GetPosition().X()*myCand->GetPosition().X()) + (myCand->GetPosition().Y()*myCand->GetPosition().Y()) + (myCand->GetPosition().Z()*myCand->GetPosition().Z()));				        
		        position.push_back(pos);
		        
		        charge.push_back(myCand->GetCharge());
		        MCindex.push_back(myCand->GetMcIndex());
		        Trackindex.push_back(myCand->GetTrackIndex());
		        // == MVD ==
		        MvdDEDX.push_back(myCand->GetMvdDEDX());
		        MvdHits.push_back(myCand->GetMvdHits());		  
		        // == STT ==
		        SttMeanDEDX.push_back(myCand->GetSttMeanDEDX());
		        SttHits.push_back(myCand->GetSttHits());
		        // == GEM ==
		        GemHits.push_back(myCand->GetGemHits());   
		        // == TOF ==
		        TofStopTime.push_back(myCand->GetTofStopTime());
		        TofM2.push_back(myCand->GetTofM2());
		        TofTrackLength.push_back(myCand->GetTofTrackLength());
		        TofQuality.push_back(myCand->GetTofQuality());
		        TofIndex.push_back(myCand->GetTofIndex());
		        TofBeta.push_back(myCand->GetTofBeta());
		        // == Barrel DIRC ==
                DrcThetaC.push_back(myCand->GetDrcThetaC());
                DrcThetaCErr.push_back(myCand->GetDrcThetaCErr());
                DrcQuality.push_back(myCand->GetDrcQuality());
                DrcNumberOfPhotons.push_back(myCand->GetDrcNumberOfPhotons());
                DrcIndex.push_back(myCand->GetDrcIndex());    
                // == Disc DIRC ==
                DiscThetaC.push_back(myCand->GetDiscThetaC());
                DiscThetaCErr.push_back(myCand->GetDiscThetaCErr());
                DiscQuality.push_back(myCand->GetDiscQuality());
                DiscNumberOfPhotons.push_back(myCand->GetDiscNumberOfPhotons());
                DiscIndex.push_back(myCand->GetDiscIndex());
                // == RICH ==
                RichThetaC.push_back(myCand->GetRichThetaC());
                RichThetaCErr.push_back(myCand->GetRichThetaCErr());
                RichQuality.push_back(myCand->GetRichQuality());
                RichNumberOfPhotons.push_back(myCand->GetRichNumberOfPhotons());
                RichIndex.push_back(myCand->GetRichIndex());
                // == EMC ==
                EmcRawEnergy.push_back(myCand->GetEmcRawEnergy());
                EmcCalEnergy.push_back(myCand->GetEmcCalEnergy());
                EmcQuality.push_back(myCand->GetEmcQuality());
                EmcNumberOfCrystals.push_back(myCand->GetEmcNumberOfCrystals());
                EmcNumberOfBumps.push_back(myCand->GetEmcNumberOfBumps());
                EmcModule.push_back(myCand->GetEmcModule());
                EmcIndex.push_back(myCand->GetEmcIndex());
                  // EMC Cluster properties
                  EmcZ20.push_back(myCand->GetEmcClusterZ20());
                  EmcZ53.push_back(myCand->GetEmcClusterZ53());
                  EmcLat.push_back(myCand->GetEmcClusterLat());
                  EmcE1.push_back(myCand->GetEmcClusterE1());
                  EmcE9.push_back(myCand->GetEmcClusterE9());
                  EmcE25.push_back(myCand->GetEmcClusterE25());
               // == Muons ==
               MuoProbability.push_back(myCand->GetMuoProbability());
               MuoQuality.push_back(myCand->GetMuoQuality());
               MuoIron.push_back(myCand->GetMuoIron());
               MuoMomentumIn.push_back(myCand->GetMuoMomentumIn());
	           MuoNumberOfLayers.push_back(myCand->GetMuoNumberOfLayers());
	           MuoModule.push_back(myCand->GetMuoModule());
	           MuoHits.push_back(myCand->GetMuoHits());
	           MuoIndex.push_back(myCand->GetMuoIndex());
	           // == Tracking ==
               DegreesOfFreedom.push_back(myCand->GetDegreesOfFreedom());
               FitStatus.push_back(myCand->GetFitStatus());
               ChiSquared.push_back(myCand->GetChiSquared());	                  
                
			}
		}
			    
		    
	}  // event loop
	
// ++++++++++++++++++++++++++++++++++++++ begin retrieve (output.csv) ++++++++++++++++++++++++++++++++++++++

   int label = 1;  // labels:  electron = 1,  muon = 2,  pion = 3,  kaon = 4,  proton = 5
   // --------------------------------------------------------------- //
   myfile<<"label,";
   for(unsigned int i = 0; i < momentumx.size(); i++){myfile<<label<<",";}
   myfile<<"\n";   
   // --------------------------------------------------------------- //	
	
   myfile<<"momentumx,";
   for(unsigned int i = 0; i < momentumx.size(); i++){myfile<<momentumx[i]<<",";}
   myfile<<"\n";
   
   myfile<<"momentumy,";
   for(unsigned int i = 0; i < momentumy.size(); i++){myfile<<momentumy[i]<<",";}
   myfile<<"\n";
   
   myfile<<"momentumz,";
   for(unsigned int i = 0; i < momentumz.size(); i++){myfile<<momentumz[i]<<",";}
   myfile<<"\n";

   myfile<<"momentum,";
   for(unsigned int i = 0; i < momentum.size(); i++){myfile<<momentum[i]<<",";}
   myfile<<"\n";
   
   myfile<<"energy,";
   for(unsigned int i = 0; i < energy.size(); i++){myfile<<energy[i]<<",";}
   myfile<<"\n";
   
   myfile<<"fhitx,";
   for(unsigned int i = 0; i < fhitx.size(); i++){myfile<<fhitx[i]<<",";}
   myfile<<"\n";
   
   myfile<<"fhity,";
   for(unsigned int i = 0; i < fhity.size(); i++){myfile<<fhity[i]<<",";}
   myfile<<"\n";
   
   myfile<<"fhitz,";
   for(unsigned int i = 0; i < fhitz.size(); i++){myfile<<fhitz[i]<<",";}
   myfile<<"\n";
      
   myfile<<"lhitx,";
   for(unsigned int i = 0; i < lhitx.size(); i++){myfile<<lhitx[i]<<",";}
   myfile<<"\n";
   
   myfile<<"lhity,";
   for(unsigned int i = 0; i < lhity.size(); i++){myfile<<lhity[i]<<",";}
   myfile<<"\n";
   
   myfile<<"lhitz,";
   for(unsigned int i = 0; i < lhitz.size(); i++){myfile<<lhitz[i]<<",";}
   myfile<<"\n";
   
   myfile<<"positionx,";
   for(unsigned int i = 0; i < positionx.size(); i++){myfile<<positionx[i]<<",";}
   myfile<<"\n";
   
   myfile<<"positiony,";
   for(unsigned int i = 0; i < positiony.size(); i++){myfile<<positiony[i]<<",";}
   myfile<<"\n";
   
   myfile<<"positionz,";
   for(unsigned int i = 0; i < positionz.size(); i++){myfile<<positionz[i]<<",";}
   myfile<<"\n";
   
   myfile<<"position,";
   for(unsigned int i = 0; i < position.size(); i++){myfile<<position[i]<<",";}
   myfile<<"\n";
   
   myfile<<"charge,";
   for(unsigned int i = 0; i < charge.size(); i++){myfile<<charge[i]<<",";}
   myfile<<"\n";
   
   myfile<<"MCindex,";
   for(unsigned int i = 0; i < MCindex.size(); i++){myfile<<MCindex[i]<<",";}
   myfile<<"\n";     
   
   myfile<<"Trackindex,";
   for(unsigned int i = 0; i < Trackindex.size(); i++){myfile<<Trackindex[i]<<",";}
   myfile<<"\n";
   
   myfile<<"MvdDEDX,";
   for(unsigned int i = 0; i < MvdDEDX.size(); i++){myfile<<MvdDEDX[i]<<",";}
   myfile<<"\n";  
   
   myfile<<"MvdHits,";
   for(unsigned int i = 0; i < MvdHits.size(); i++){myfile<<MvdHits[i]<<",";}
   myfile<<"\n";  
   
   myfile<<"SttMeanDEDX,";
   for(unsigned int i = 0; i < SttMeanDEDX.size(); i++){myfile<<SttMeanDEDX[i]<<",";}
   myfile<<"\n";  
   
   myfile<<"SttHits,";
   for(unsigned int i = 0; i < SttHits.size(); i++){myfile<<SttHits[i]<<",";}
   myfile<<"\n";  
   
   myfile<<"GemHits,";
   for(unsigned int i = 0; i < GemHits.size(); i++){myfile<<GemHits[i]<<",";}
   myfile<<"\n";  
   
   myfile<<"TofStopTime,";
   for(unsigned int i = 0; i < TofStopTime.size(); i++){myfile<<TofStopTime[i]<<",";}
   myfile<<"\n";  
   
   myfile<<"TofM2,";
   for(unsigned int i = 0; i < TofM2.size(); i++){myfile<<TofM2[i]<<",";}
   myfile<<"\n";  
   
   myfile<<"TofTrackLength,";
   for(unsigned int i = 0; i < TofTrackLength.size(); i++){myfile<<TofTrackLength[i]<<",";}
   myfile<<"\n";  
   
   myfile<<"TofQuality,";
   for(unsigned int i = 0; i < TofQuality.size(); i++){myfile<<TofQuality[i]<<",";}
   myfile<<"\n";  
   
   myfile<<"TofIndex,";
   for(unsigned int i = 0; i < TofIndex.size(); i++){myfile<<TofIndex[i]<<",";}
   myfile<<"\n";  
   
   myfile<<"TofBeta,";
   for(unsigned int i = 0; i < TofBeta.size(); i++){myfile<<TofBeta[i]<<",";}
   myfile<<"\n"; 
   
   myfile<<"DrcThetaC,";
   for(unsigned int i = 0; i < DrcThetaC.size(); i++){myfile<<DrcThetaC[i]<<",";}
   myfile<<"\n";
   
   myfile<<"DrcThetaCErr,";
   for(unsigned int i = 0; i < DrcThetaCErr.size(); i++){myfile<<DrcThetaCErr[i]<<",";}
   myfile<<"\n"; 
   
   myfile<<"DrcQuality,";
   for(unsigned int i = 0; i < DrcQuality.size(); i++){myfile<<DrcQuality[i]<<",";}
   myfile<<"\n"; 
   
   myfile<<"DrcNumberOfPhotons,";
   for(unsigned int i = 0; i < DrcNumberOfPhotons.size(); i++){myfile<<DrcNumberOfPhotons[i]<<",";}
   myfile<<"\n"; 
   
   myfile<<"DrcIndex,";
   for(unsigned int i = 0; i < DrcIndex.size(); i++){myfile<<DrcIndex[i]<<",";}
   myfile<<"\n"; 
   
   myfile<<"DiscThetaC,";
   for(unsigned int i = 0; i < DiscThetaC.size(); i++){myfile<<DiscThetaC[i]<<",";}
   myfile<<"\n"; 
   
   myfile<<"DiscThetaCErr,";
   for(unsigned int i = 0; i < DiscThetaCErr.size(); i++){myfile<<DiscThetaCErr[i]<<",";}
   myfile<<"\n"; 
   
   myfile<<"DiscQuality,";
   for(unsigned int i = 0; i < DiscQuality.size(); i++){myfile<<DiscQuality[i]<<",";}
   myfile<<"\n"; 
   
   myfile<<"DiscNumberOfPhotons,";
   for(unsigned int i = 0; i < DiscNumberOfPhotons.size(); i++){myfile<<DiscNumberOfPhotons[i]<<",";}
   myfile<<"\n"; 
   
   myfile<<"DiscIndex,";
   for(unsigned int i = 0; i < DiscIndex.size(); i++){myfile<<DiscIndex[i]<<",";}
   myfile<<"\n"; 
   
   myfile<<"RichThetaC,";
   for(unsigned int i = 0; i < RichThetaC.size(); i++){myfile<<RichThetaC[i]<<",";}
   myfile<<"\n"; 
   
   myfile<<"RichThetaCErr,";
   for(unsigned int i = 0; i < RichThetaCErr.size(); i++){myfile<<RichThetaCErr[i]<<",";}
   myfile<<"\n";
   
   myfile<<"RichQuality,";
   for(unsigned int i = 0; i < RichQuality.size(); i++){myfile<<RichQuality[i]<<",";}
   myfile<<"\n";
   
   myfile<<"RichNumberOfPhotons,";
   for(unsigned int i = 0; i < RichNumberOfPhotons.size(); i++){myfile<<RichNumberOfPhotons[i]<<",";}
   myfile<<"\n";
   
   myfile<<"RichIndex,";
   for(unsigned int i = 0; i < RichIndex.size(); i++){myfile<<RichIndex[i]<<",";}
   myfile<<"\n";
   
   myfile<<"EmcRawEnergy,";
   for(unsigned int i = 0; i < EmcRawEnergy.size(); i++){myfile<<EmcRawEnergy[i]<<",";}
   myfile<<"\n";
   
   myfile<<"EmcCalEnergy,";
   for(unsigned int i = 0; i < EmcCalEnergy.size(); i++){myfile<<EmcCalEnergy[i]<<",";}
   myfile<<"\n";
   
   myfile<<"EmcQuality,";
   for(unsigned int i = 0; i < EmcQuality.size(); i++){myfile<<EmcQuality[i]<<",";}
   myfile<<"\n";
   
   myfile<<"EmcNumberOfCrystals,";
   for(unsigned int i = 0; i < EmcNumberOfCrystals.size(); i++){myfile<<EmcNumberOfCrystals[i]<<",";}
   myfile<<"\n";
   
   myfile<<"EmcNumberOfBumps,";
   for(unsigned int i = 0; i < EmcNumberOfBumps.size(); i++){myfile<<EmcNumberOfBumps[i]<<",";}
   myfile<<"\n";
   
   myfile<<"EmcModule,";
   for(unsigned int i = 0; i < EmcModule.size(); i++){myfile<<EmcModule[i]<<",";}
   myfile<<"\n";
   
   myfile<<"EmcIndex,";
   for(unsigned int i = 0; i < EmcIndex.size(); i++){myfile<<EmcIndex[i]<<",";}
   myfile<<"\n";
   
   myfile<<"EmcZ20,";
   for(unsigned int i = 0; i < EmcZ20.size(); i++){myfile<<EmcZ20[i]<<",";}
   myfile<<"\n";
   
   myfile<<"EmcZ53,";
   for(unsigned int i = 0; i < EmcZ53.size(); i++){myfile<<EmcZ53[i]<<",";}
   myfile<<"\n";                        
   
   myfile<<"EmcLat,";
   for(unsigned int i = 0; i < EmcLat.size(); i++){myfile<<EmcLat[i]<<",";}
   myfile<<"\n";

   myfile<<"EmcE1,";
   for(unsigned int i = 0; i < EmcE1.size(); i++){myfile<<EmcE1[i]<<",";}
   myfile<<"\n"; 
   
   myfile<<"EmcE9,";
   for(unsigned int i = 0; i < EmcE9.size(); i++){myfile<<EmcE9[i]<<",";}
   myfile<<"\n";
   
   myfile<<"EmcE25,";
   for(unsigned int i = 0; i < EmcE25.size(); i++){myfile<<EmcE25[i]<<",";}
   myfile<<"\n";      
   
   myfile<<"MuoProbability,";
   for(unsigned int i = 0; i < MuoProbability.size(); i++){myfile<<MuoProbability[i]<<",";}
   myfile<<"\n";   
   
   myfile<<"MuoQuality,";
   for(unsigned int i = 0; i < MuoQuality.size(); i++){myfile<<MuoQuality[i]<<",";}
   myfile<<"\n";   
   
   myfile<<"MuoIron,";
   for(unsigned int i = 0; i < MuoIron.size(); i++){myfile<<MuoIron[i]<<",";}
   myfile<<"\n";   
   
   myfile<<"MuoMomentumIn,";
   for(unsigned int i = 0; i < MuoMomentumIn.size(); i++){myfile<<MuoMomentumIn[i]<<",";}
   myfile<<"\n";   
   
   myfile<<"MuoNumberOfLayers,";
   for(unsigned int i = 0; i < MuoNumberOfLayers.size(); i++){myfile<<MuoNumberOfLayers[i]<<",";}
   myfile<<"\n";   
   
   myfile<<"MuoModule,";
   for(unsigned int i = 0; i < MuoModule.size(); i++){myfile<<MuoModule[i]<<",";}
   myfile<<"\n";   
   
   myfile<<"MuoHits,";
   for(unsigned int i = 0; i < MuoHits.size(); i++){myfile<<MuoHits[i]<<",";}
   myfile<<"\n";   
   
   myfile<<"MuoIndex,";
   for(unsigned int i = 0; i < MuoIndex.size(); i++){myfile<<MuoIndex[i]<<",";}
   myfile<<"\n";     
   
   myfile<<"DegreesOfFreedom,";
   for(unsigned int i = 0; i < DegreesOfFreedom.size(); i++){myfile<<DegreesOfFreedom[i]<<",";}
   myfile<<"\n"; 
   
   myfile<<"FitStatus,";
   for(unsigned int i = 0; i < FitStatus.size(); i++){myfile<<FitStatus[i]<<",";}
   myfile<<"\n";   
   
   myfile<<"ChiSquared,";
   for(unsigned int i = 0; i < ChiSquared.size(); i++){myfile<<ChiSquared[i]<<",";}
   myfile<<"\n";     
            
// ++++++++++++++++++++++++++++++++++++++ end retrieve (output.csv) ++++++++++++++++++++++++++++++++++++++++
	
} // macro end


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
            
// ++++++++++++++++++++++++++++++++++++++ end retrieve (output.csv) ++++++++++++++++++++++++++++++++++++++++
	
} // macro end


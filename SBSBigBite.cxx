//////////////////////////////////////////////////////////////////////////
//
// Bare-bones SBB BigBite spectrometer class 
// Used for testing.
//
//////////////////////////////////////////////////////////////////////////

#include "SBSBigBite.h"
#include "THaTrack.h"
#include "THaPIDinfo.h"
#include "TList.h"
#include "SBSTimingHodoscope.h"
#include "SBSGEMSpectrometerTracker.h"
#include "SBSRasteredBeam.h"
#include "THaTrackingDetector.h"
#include "TH2D.h"

using namespace std;

ClassImp(SBSBigBite)

//_____________________________________________________________________________
SBSBigBite::SBSBigBite( const char* name, const char* description ) :
THaSpectrometer( name, description )
{
  SetMultiTracks(false);
  SetTrSorting(false);
  //fTrackerPitchAngle = 10.0;
  //fDetectorStackYaw = 0;
  //fDetectorStackRoll = 0;
  fOpticsOrder = -1;
  fFrontConstraintWidthX = 1.5;
  fFrontConstraintWidthY = 1.5;
  fBackConstraintWidthX = 1.5;
  fBackConstraintWidthY = 1.5;
  fFrontConstraintX0 = 0.0;
  fFrontConstraintY0 = 0.0;
  fBackConstraintX0 = 0.0;
  fBackConstraintY0 = 0.0;
  fTrackGrinchClusCorr_0 = 0.0;
  fTrackGrinchClusCorr_1 = 0.0;
  fTrackGrinchClusCorr_Sigma = 1.5;
  fPtheta_00000 = 0.0;
  fPtheta_10000 = 0.0;
  fPtheta_00100 = 0.0;
  fXptar_10000 = 0.0;
  fXptar_00100 = 0.0;
  fYtar_01000 = 0.0;
  fYtar_00010 = 0.0;
  fECaloFudgeFactor = 1.0;

  //Default ideal optics angle (central bend angle) to 10 deg:
    
  //This should match the g4sbs simulation coordinate system, and user should not ordinarily need to modify:
  //This is relative to target center along spectrometer axis:
  fMagDist = 1.80; //this will eventually be overridden by the DB, and will become a required parameter for BigBite
    
  fOpticsAngle = 10.0*TMath::DegToRad();
  fOpticsOrigin.SetXYZ( -0.1701, 0.0, 1.1087 + fMagDist );
    
  InitOpticsAxes( fOpticsAngle );
    
  //These will be the real GEM coordinates from the zero-field alignment:
  //Default these to match the simulation coordinate system:
  fGEMtheta = fOpticsAngle;
  fGEMphi   = 180.0*TMath::DegToRad();
  fGEMorigin = fOpticsOrigin;
    
  InitGEMAxes( fGEMtheta, fGEMphi, fGEMorigin );
    
  // Define our own PID parameters (undo what THaSpectrometer has set up)
  SBSBigBite::DefinePidParticles();
    
  // Enable PID calculations (call CalcPID)
  SetPID(true);
    
  fDownBendingMode = false;
    
  //Default
  fPrecon_flag = 0;
    
  fA_pth1 = 0.28615 * 0.97;
  fB_pth1 = 0.1976;
  fC_pth1 = 0.4764;

  fA_vy = 0.;
  fB_vy = 0.;

  fUseBeamPosInOptics = false;

  fIsMC = false;

  //Default to block size/sqrt(12) for shower and preshower:
  fSigmaX_shower = 0.085/sqrt(12.0); //2.45 cm
  fSigmaY_shower = 0.085/sqrt(12.0);
  fSigmaX_preshower = 0.09/sqrt(12.0); //2.6 cm
  fSigmaY_preshower = 0.37/sqrt(12.0); //10.7 cm
    
  //The hodoscope position resolutions are based on a fit of the hodo-track differences along X and Y:
  //NOTE: for X it is worse than bar size over sqrt(12), presumably the averaging of several bars does NOT improve the vertical position resolution of the hodoscope much
  // for Y it is based on the timing resolution
    
  fSigmaX_hodo = 0.021;
  fSigmaY_hodo = 0.041;
    
  // Constructor. Defines standard detectors
  //The standard BigBite detector package in the 12 GeV/SBS era will include:
  // pre-shower + shower calorimeters (inherit from THaNonTrackingDetector OR THaPidDetector)
  // Timing hodoscope (inherit from THaNonTrackingDetector)
  // GRINCH (inherit from THaPidDetector)
  // GEMs (five-layer) (inherit from THaTrackingDetector)
  /*
    h1_yVx_bcp = new TH2D("h1_yVx_bcp", ";x_{bcp} (m);y_{bcp} (m)", 300, -1.5, 1.5, 100, -0.5, 0.5);
    h1_x_fcpVbcp = new TH2D("h1_x_fcpVbcp", ";x_{bcp} (m);x_{fcp} (m)", 300, -1.5, 1.5, 300, -1.5, 1.5);
    h1_yVx_fcp = new TH2D("h1_yVx_fcp", ";x_{fcp} (m);y_{fcp} (m)", 300, -1.5, 1.5, 100, -0.5, 0.5);
    h1_y_fcpVbcp = new TH2D("h1_y_fcpVbcp", ";y_{bcp} (m);y_{fcp} (m)", 100, -0.5, 0.5, 100, -0.5, 0.5);
    h1_dyVdx = new TH2D("h1_dyVdx",";dx/dz;dy/dz", 100, -0.5, 0.5, 50, -0.25, 0.25);
  */
    
  fUseForwardOptics = false;
  fForwardOpticsOrder = -1;

}

//_____________________________________________________________________________
SBSBigBite::~SBSBigBite()
{
  // Destructor
}

void SBSBigBite::Clear( Option_t *opt )
{
  THaSpectrometer::Clear(opt);
  //  f_xtg_exp.clear();
  fEpsEtotRatio.clear();
  fEtot.clear();
  fEtotPratio.clear();
  fEpsEtotRatio.clear();
  fFrontConstraintX.clear();
  fFrontConstraintY.clear();
  fFrontConstraintZ.clear();
  fBackConstraintX.clear();
  fBackConstraintY.clear();
  fBackConstraintZ.clear();
  fProbaE.clear();
  fProbaPi.clear();
}

//_____________________________________________________________________________
void SBSBigBite::DefinePidParticles()
{
  // Define the default set of PID particles:
  //  electron, pion
    
  fPidParticles->Delete();    //make sure array is empty
    
  AddPidParticle( "e",  "electron",   0.511e-3, -1 );
  AddPidParticle( "pi", "pion",   0.139, -1 );
}

//_____________________________________________________________________________
Int_t SBSBigBite::ReadRunDatabase( const TDatime &date ){
  //const char* const here = "SBSBigBite::ReadRunDatabase()";
    
  Int_t err = THaSpectrometer::ReadRunDatabase( date );
  if( err ) return err;
    
  FILE* file = OpenRunDBFile( date );
  if( !file ) return kFileError;
    
  //Require magdist:
  const DBRequest req[] = {
    { "magdist", &fMagDist, kDouble, 0, 0, 1 },
    { nullptr }
  };
  err = LoadDB( file, date, req );
  fclose(file);
  if( err )
    return kInitError;
    
  fOpticsOrigin.SetXYZ( -0.1701, 0.0, 1.1087+fMagDist );
    
  //Default GEM origin to the same as optics origin; this will be overridden by ReadDatabase
  fGEMorigin = fOpticsOrigin;
    
  return kOK;
}

//_____________________________________________________________________________
Int_t SBSBigBite::ReadDatabase( const TDatime& date )
{
  // Hack from THaVDC::ReadDatabase()
  const char* const here = "SBSBigBite::ReadDatabase()";
    
  //THaSpectrometer::ReadRunDatabase();
    
  FILE* file = OpenFile( date );
  if( !file ){
    std::cerr << here << "(): database not found!"<< std::endl;
    return kFileError;
  }
    
  int pidflag = fPID ? 1 : 0;
  int downbend = fDownBendingMode ? 1 : 0;
    
  std::vector<Double_t> firstgem_offset;
    
  std::vector<Double_t> optics_param;
  std::vector<Double_t> foptics_param;
  std::vector<Double_t> doptics_param; //down bending optics
    
  std::vector<Double_t> pssh_pidproba;
  std::vector<Double_t> pcal_pidproba;
  std::vector<Double_t> grinch_pidproba;
    
  std::vector<Double_t> optics_origin;
    
  double gemthetadeg = fGEMtheta * TMath::RadToDeg();
  double gemphideg   = fGEMphi * TMath::RadToDeg();
  double opticsthetadeg = fOpticsAngle * TMath::RadToDeg();

  int mc_flag = fIsMC ? 1 : 0;
  int use_beampos = fUseBeamPosInOptics ? 1 : 0;
    
  const DBRequest request[] = {
    { "gemtheta", &gemthetadeg, kDouble, 0, 1, 1},
    { "gemphi", &gemphideg, kDouble, 0, 1, 1},
    { "gemorigin_xyz",    &firstgem_offset, kDoubleV,  0, 1, 1},
    { "opticstheta", &opticsthetadeg, kDouble, 0, 1, 1},
    { "optics_origin", &optics_origin, kDoubleV, 0, 1, 1},
    { "optics_order",    &fOpticsOrder, kUInt,  0, 1, 1},
    { "optics_parameters", &optics_param, kDoubleV, 0, 1, 1},
    { "ecalo_fudgefactor", &fECaloFudgeFactor, kDouble, 0, 1, 1},
    { "do_pid",    &pidflag, kInt,  0, 1, 1},
    { "frontconstraintwidth_x", &fFrontConstraintWidthX, kDouble, 0, 1, 0},
    { "frontconstraintwidth_y", &fFrontConstraintWidthY, kDouble, 0, 1, 0},
    { "backconstraintwidth_x", &fBackConstraintWidthX, kDouble, 0, 1, 0},
    { "backconstraintwidth_y", &fBackConstraintWidthY, kDouble, 0, 1, 0},
    { "frontconstraint_x0", &fFrontConstraintX0, kDouble, 0, 1, 0},
    { "frontconstraint_y0", &fFrontConstraintY0, kDouble, 0, 1, 0},
    { "backconstraint_x0", &fBackConstraintX0, kDouble, 0, 1, 0},
    { "backconstraint_y0", &fBackConstraintY0, kDouble, 0, 1, 0},
    { "trackgrinchcorr_const", &fTrackGrinchClusCorr_0, kDouble, 0, 1, 0},
    { "trackgrinchcorr_slope", &fTrackGrinchClusCorr_1, kDouble, 0, 1, 0},
    { "trackgrinchcorr_sigma", &fTrackGrinchClusCorr_Sigma, kDouble, 0, 1, 0},
    { "psshPIDprobatable",    &pssh_pidproba, kDoubleV,  0, 1, 0},
    { "pcalPIDprobatable",    &pcal_pidproba, kDoubleV,  0, 1, 0},
    { "grinchPIDpbins",    &fP_table, kDoubleV,  0, 1, 0},
    { "grinchPIDprobatable",    &grinch_pidproba, kDoubleV,  0, 1, 0},
    { "preconflag", &fPrecon_flag, kUInt, 0, 1, 1 },
    { "A_pth1", &fA_pth1, kDouble, 0, 1, 1 },
    { "B_pth1", &fB_pth1, kDouble, 0, 1, 1 },
    { "C_pth1", &fC_pth1, kDouble, 0, 1, 1 },
    { "A_pvy", &fA_vy, kDouble, 0, 1, 1 },
    { "B_pvy", &fB_vy, kDouble, 0, 1, 1 },
    { "xsigma_shower", &fSigmaX_shower, kDouble, 0, 1, 1 },
    { "ysigma_shower", &fSigmaY_shower, kDouble, 0, 1, 1 },
    { "xsigma_preshower", &fSigmaX_preshower, kDouble, 0, 1, 1 },
    { "ysigma_preshower", &fSigmaY_preshower, kDouble, 0, 1, 1 },
    { "xsigma_hodo", &fSigmaX_hodo, kDouble, 0, 1, 1 },
    { "ysigma_hodo", &fSigmaY_hodo, kDouble, 0, 1, 1 },
    { "forwardoptics_order", &fForwardOpticsOrder, kInt, 0, 1, 1 },
    { "forwardoptics_parameters", &foptics_param, kDoubleV, 0, 1, 1 },
    { "downbendoptics_order", &fOpticsOrderDownbend, kInt, 0, 1, 1 },
    { "downbendoptics_parameters", &doptics_param, kDoubleV, 0, 1, 1 },
    { "downbending_mode", &downbend, kInt, 0, 1, 1 },
    { "use_beampos", &use_beampos, kInt, 0, 1, 1 },
    { "is_mc",        &mc_flag,    kInt, 0, 1, 1 },
    {0}
  };
    
  Int_t status = LoadDB( file, date, request, fPrefix, 1 ); //The "1" after fPrefix means search up the tree
  fclose(file);
  if( status != 0 ){
    return status;
  }

  fUseBeamPosInOptics = ( use_beampos != 0 ) ? true : false;
  
  fOpticsAngle = opticsthetadeg * TMath::DegToRad();
  if( optics_origin.size() == 3 ){ //database overrides default values:
    fOpticsOrigin.SetXYZ( optics_origin[0],
			  optics_origin[1],
			  optics_origin[2] );
  }
    
  InitOpticsAxes( fOpticsAngle );
    
  //Database values for these angles are assumed to be given in degrees:
  // These will have been initialized to default values in degrees above, or if they were loaded from the database they will have been given in degrees:
  fGEMtheta = gemthetadeg * TMath::DegToRad();
  fGEMphi = gemphideg * TMath::DegToRad();
  if( firstgem_offset.size() == 3 ){ //database overrides default values:
    fGEMorigin.SetXYZ( firstgem_offset[0],
		       firstgem_offset[1],
		       firstgem_offset[2] );
  }
    
  InitGEMAxes( fGEMtheta, fGEMphi );
    
  if(fECaloFudgeFactor!=1.0)cout << "Setting the calorimeter energy fudge factor to " << fECaloFudgeFactor << endl;
    
  fPID = ( pidflag != 0 );
    
  //do we have non tracking detectors
  bool nontrackdet = false;
  TIter next( fNonTrackingDetectors );
  while( auto* theNonTrackDetector =
	 static_cast<THaNonTrackingDetector*>( next() )) {
    cout << theNonTrackDetector->GetName() << endl;
    nontrackdet = true;
  }//if we do not find any non tracking detectors, we force fPID to be false.
  if(!nontrackdet)fPID = false;
    
  if(fOpticsOrder>=0){
    unsigned int nterms = 0;
        
    for(int i = 0; i<=fOpticsOrder; i++){ //x
      for( int j=0; j<=fOpticsOrder-i; j++){ //y
	for( int k=0; k<=fOpticsOrder-i-j; k++){
	  for( int l=0; l<=fOpticsOrder-i-j-k; l++){
	    for( int m=0; m<=fOpticsOrder-i-j-k-l; m++ ){
	      nterms++;
	    }
	  }
	}
      }
    }
    cout << nterms << " lines of parameters expected for optics of order " << fOpticsOrder << endl;
        
    //int n_elem = TMath::FloorNint(optics_param.size()/nparam);
        
    //we expect 9 parameters per line: four coefficients plus five exponents:
        
    unsigned int nparams = 9*nterms;
        
    if(nparams!=optics_param.size()){
      std::cerr << "Warning: mismatch between " << optics_param.size()
		<< " optics parameters provided and " << nparams*9
		<< " optics parameters expected!" << std::endl;
      std::cerr << " Fix database! " << std::endl;
      return kInitError;
    }
        
    //int o_i, o_j, o_k, o_l, o_m;// shall we use those???
    fb_xptar.resize(nterms);
    fb_yptar.resize(nterms);
    fb_ytar.resize(nterms);
    fb_pinv.resize(nterms);
    f_oi.resize(nterms);
    f_oj.resize(nterms);
    f_ok.resize(nterms);
    f_ol.resize(nterms);
    f_om.resize(nterms);
        
    for(unsigned int i=0; i<nterms; i++){
      fb_xptar[i] = optics_param[9*i];
      fb_yptar[i] = optics_param[9*i+1];
      fb_ytar[i] = optics_param[9*i+2];
      fb_pinv[i] = optics_param[9*i+3];
      f_om[i] = int(optics_param[9*i+4]);
      f_ol[i] = int(optics_param[9*i+5]);
      f_ok[i] = int(optics_param[9*i+6]);
      f_oj[i] = int(optics_param[9*i+7]);
      f_oi[i] = int(optics_param[9*i+8]);
      if(f_om[i]==0 && f_ol[i]==0 && f_ok[i]==0 && f_oj[i]==0 && f_oi[i]==0){
	fPtheta_00000 = fb_pinv[i];
      }
      if(f_om[i]==1 && f_ol[i]==0 && f_ok[i]==0 && f_oj[i]==0 && f_oi[i]==0){
	fPtheta_10000 = fb_pinv[i];
	fXptar_10000 = fb_xptar[i];
      }
      if(f_om[i]==0 && f_ol[i]==0 && f_ok[i]==1 && f_oj[i]==0 && f_oi[i]==0){
	fPtheta_00100 = fb_pinv[i];
	fXptar_00100 = fb_xptar[i];
      }
      if(f_om[i]==0 && f_ol[i]==1 && f_ok[i]==0 && f_oj[i]==0 && f_oi[i]==0)
	fYtar_01000 = fb_ytar[i];
      if(f_om[i]==0 && f_ol[i]==0 && f_ok[i]==0 && f_oj[i]==1 && f_oi[i]==0)
	fYtar_00010 = fb_ytar[i];
    }
  }
    
  if( fForwardOpticsOrder >= 0 ){
    unsigned int nterms=0;
    for(int i = 0; i<=fForwardOpticsOrder; i++){ //x
      for( int j=0; j<=fForwardOpticsOrder-i; j++){ //y
	for( int k=0; k<=fForwardOpticsOrder-i-j; k++){
	  for( int l=0; l<=fForwardOpticsOrder-i-j-k; l++){
	    for( int m=0; m<=fForwardOpticsOrder-i-j-k-l; m++ ){
	      nterms++;
	    }
	  }
	}
      }
    }
    cout << nterms << " lines of parameters expected for optics of order " << fForwardOpticsOrder << endl;
        
    unsigned int nparams = 9*nterms;
    if(nparams!=foptics_param.size()){
      std::cerr << "Warning: mismatch between " << foptics_param.size()
		<< " forward optics parameters provided and " << nparams*9
		<< " forward optics parameters expected!" << std::endl;
      std::cerr << " Fix database! " << std::endl;
      return kInitError;
    }
        
    //I seem to be able to write anything except a std::cout
        
    fb_xfp.resize( nterms );
    fb_yfp.resize( nterms );
    fb_xpfp.resize( nterms );
    fb_ypfp.resize( nterms );
        
    f_foi.resize(nterms);
    f_foj.resize(nterms);
    f_fok.resize(nterms);
    f_fol.resize(nterms);
    f_fom.resize(nterms);
        
    for(unsigned int i=0; i<nterms; i++){
      fb_xfp[i] = foptics_param[9*i];
      fb_yfp[i] = foptics_param[9*i+1];
      fb_xpfp[i] = foptics_param[9*i+2];
      fb_ypfp[i] = foptics_param[9*i+3];
      f_fom[i] = int(foptics_param[9*i+4]);
      f_fol[i] = int(foptics_param[9*i+5]);
      f_fok[i] = int(foptics_param[9*i+6]);
      f_foj[i] = int(foptics_param[9*i+7]);
      f_foi[i] = int(foptics_param[9*i+8]);
    }
  }
    
  fDownBendingMode = downbend != 0 ? true : false;
  if( fDownBendingMode ){
    if( fOpticsOrderDownbend >= 0 ){
      unsigned int nterms=0;
      for(int i = 0; i<=fOpticsOrderDownbend; i++){ //x
	for( int j=0; j<=fOpticsOrderDownbend-i; j++){ //y
	  for( int k=0; k<=fOpticsOrderDownbend-i-j; k++){
	    for( int l=0; l<=fOpticsOrderDownbend-i-j-k; l++){
	      for( int m=0; m<=fOpticsOrderDownbend-i-j-k-l; m++ ){
		nterms++;
	      }
	    }
	  }
	}
      }
      cout << nterms << " lines of parameters expected for optics of order " << fOpticsOrderDownbend << endl;
            
      unsigned int nparams = 9*nterms;
      if(nparams!=doptics_param.size()){
	std::cerr << "Warning: mismatch between " << foptics_param.size()
		  << " downbending optics parameters provided and " << nparams*9
		  << " downbending optics parameters expected!" << std::endl;
	std::cerr << " Fix database! " << std::endl;
	return kInitError;
      }
            
      fb_xptar_downbend.resize(nterms);
      fb_yptar_downbend.resize(nterms);
      fb_ytar_downbend.resize(nterms);
      fb_pinv_downbend.resize(nterms);
      f_oi_downbend.resize(nterms);
      f_oj_downbend.resize(nterms);
      f_ok_downbend.resize(nterms);
      f_ol_downbend.resize(nterms);
      f_om_downbend.resize(nterms);
	    
      for(unsigned int i=0; i<nterms; i++){
	fb_xptar_downbend[i] = doptics_param[9*i];
	fb_yptar_downbend[i] = doptics_param[9*i+1];
	fb_ytar_downbend[i] = doptics_param[9*i+2];
	fb_pinv_downbend[i] = doptics_param[9*i+3];
	f_om_downbend[i] = int(doptics_param[9*i+4]);
	f_ol_downbend[i] = int(doptics_param[9*i+5]);
	f_ok_downbend[i] = int(doptics_param[9*i+6]);
	f_oj_downbend[i] = int(doptics_param[9*i+7]);
	f_oi_downbend[i] = int(doptics_param[9*i+8]);
	if(f_om_downbend[i]==0 && f_ol_downbend[i]==0 && f_ok_downbend[i]==0 && f_oj_downbend[i]==0 && f_oi_downbend[i]==0){
	  fPtheta_00000 = fb_pinv_downbend[i];
	}
	if(f_om_downbend[i]==1 && f_ol_downbend[i]==0 && f_ok_downbend[i]==0 && f_oj_downbend[i]==0 && f_oi_downbend[i]==0){
	  fPtheta_10000 = fb_pinv_downbend[i];
	  fXptar_10000 = fb_xptar_downbend[i];
	}
	if(f_om_downbend[i]==0 && f_ol_downbend[i]==0 && f_ok_downbend[i]==1 && f_oj_downbend[i]==0 && f_oi_downbend[i]==0){
	  fPtheta_00100 = fb_pinv_downbend[i];
	  fXptar_00100 = fb_xptar_downbend[i];
	}
	if(f_om_downbend[i]==0 && f_ol_downbend[i]==1 && f_ok_downbend[i]==0 && f_oj_downbend[i]==0 && f_oi_downbend[i]==0)
	  fYtar_01000 = fb_ytar_downbend[i];
	if(f_om_downbend[i]==0 && f_ol_downbend[i]==0 && f_ok_downbend[i]==0 && f_oj_downbend[i]==1 && f_oi_downbend[i]==0)
	  fYtar_00010 = fb_ytar_downbend[i];
                
      }
            
    } else {
      std::cerr << "Warning: downbending mode specified but no downbending optics parameters provided; fix database" << std::endl;
      return kInitError;
    }
  }
    
  //PID stuff
  fEpsEtotRatio_table.clear();
  fProba_e_PSSH_table.clear();
  fProba_pi_PSSH_table.clear();
    
  if(!pssh_pidproba.empty()){
    unsigned int npts = pssh_pidproba.size()/3;
    fEpsEtotRatio_table.resize(npts);
    fProba_e_PSSH_table.resize(npts);
    fProba_pi_PSSH_table.resize(npts);
        
    for(unsigned int i = 0; i<npts; i++){
      fEpsEtotRatio_table[i] = pssh_pidproba[3*i];
      fProba_e_PSSH_table[i] = pssh_pidproba[3*i+1];
      fProba_pi_PSSH_table[i] = pssh_pidproba[3*i+2];
    }
  }
    
  //PID stuff
  fEtotPratio_table.clear();
  fProba_e_PCAL_table.clear();
  fProba_pi_PCAL_table.clear();
    
  if(!pcal_pidproba.empty()){
    unsigned int npts = pcal_pidproba.size()/3;
    fEtotPratio_table.resize(npts);
    fProba_e_PCAL_table.resize(npts);
    fProba_pi_PCAL_table.resize(npts);
        
    for(unsigned int i = 0; i<npts; i++){
      fEtotPratio_table[i] = pcal_pidproba[3*i];
      fProba_e_PCAL_table[i] = pcal_pidproba[3*i+1];
      fProba_pi_PCAL_table[i] = pcal_pidproba[3*i+2];
    }
  }
    
  fNGRINCHPMTs_table.clear();
  fProba_e_GRINCH_table.clear();
  fProba_pi_GRINCH_table.clear();
    
  if(!grinch_pidproba.empty() && !fP_table.empty()){
    unsigned int n_ppts = 2+fP_table.size();
    unsigned int npts = grinch_pidproba.size()/(n_ppts);
    fNGRINCHPMTs_table.resize(npts);
    fProba_e_GRINCH_table.resize(npts);
    fProba_pi_GRINCH_table.resize(fP_table.size());
    for(unsigned int j = 0; j<fP_table.size(); j++){
      fProba_pi_GRINCH_table[j].resize(npts);
    }
    for(unsigned int i = 0; i<npts; i++){
      fNGRINCHPMTs_table[i] = grinch_pidproba[n_ppts*i];
      fProba_e_GRINCH_table[i] = grinch_pidproba[n_ppts*i+1];
      for(unsigned int j = 0; j<fP_table.size(); j++){
	fProba_pi_GRINCH_table[j][i] = grinch_pidproba[n_ppts*i+2+j];
      }
    }
  }

  fIsMC = (mc_flag != 0);

  fIsInit = true;
  return kOK;
}

Int_t SBSBigBite::DefineVariables( EMode mode ){
  THaSpectrometer::DefineVariables(mode);
    
  if( mode == kDefine and fIsSetup ) return kOK;
  fIsSetup = ( mode == kDefine );
    
    
  // removing all that stuff since apparently I'm not able to code properly...
    
  RVarDef beamtrackvars[] = {
    { "tr.tg_x", "target x", "fTracks.THaTrack.fTX" },
    { nullptr }
  };
  DefineVarsFromList( beamtrackvars, mode );
    
  RVarDef constraintvars[] = {
    { "x_fcp", "front track constraint x", "fFrontConstraintX" },
    { "y_fcp", "front track constraint y", "fFrontConstraintY" },
    { "z_fcp", "front track constraint z", "fFrontConstraintZ" },
    { "x_bcp", "back track constraint x", "fBackConstraintX" },
    { "y_bcp", "back track constraint y", "fBackConstraintY" },
    { "z_bcp", "back track constraing z", "fBackConstraintZ" },
    { nullptr }
  };
  DefineVarsFromList( constraintvars, mode );
    
  RVarDef pidvars[] = {
    { "eps_over_etot", "electron probability", "fEpsEtotRatio" },
    { "etot_over_p", "electron probability", "fEtotPratio" },
    { "prob_e", "electron probability", "fProbaE" },
    { "prob_pi", "pion probability", "fProbaPi" },
    { nullptr }
  };
  DefineVarsFromList( pidvars, mode );
    
  return 0;
}

//_____________________________________________________________________________
// Int_t SBSBigBite::End( THaRunBase* run )
// {
//   h1_yVx_bcp->Write();
//   h1_x_fcpVbcp->Write();
//   h1_yVx_fcp->Write();
//   h1_y_fcpVbcp->Write();
//   h1_dyVdx->Write();
// } 

//_____________________________________________________________________________
Int_t SBSBigBite::CoarseTrack()
{
  // Coarse track Reconstruction
  //std::cout << " SBSBigBite::CoarseTrack()...";
  THaSpectrometer::CoarseTrack();
  // TODO
  //std::cout << " call SBSBigBite::CoarseTrack" << std::endl;
  //std::cout << "done" << std::endl;
  return 0;
}

//_____________________________________________________________________________
Int_t SBSBigBite::CoarseReconstruct()
{
    
  //std::cout << "SBSBigBite::CoarseReconstruct()...";
  // Coarse Reconstruction of particles in spectrometer
  THaSpectrometer::CoarseReconstruct();
  // TODO
  // fetch the clusters from SBSBBShower detectors
  // FOR NOW: fetch the highest clusters from SBSBBShower detectors
  double x_fcp = 0, y_fcp = 0, z_fcp = 0;
  double x_bcp = 0, y_bcp = 0, z_bcp = 0;
  double sumweights_x = 0, sumweights_y = 0, sumweights_z = 0.0;
  double Etot = 0;
  //npts is incremented only if there are clusters in the preshower and shower
  int npts = 0;
  double EpsEtotRatio = 0;
  TIter next( fNonTrackingDetectors );
    
  return 0;
}

//_____________________________________________________________________________
Int_t SBSBigBite::Track()
{
  // Fine track Reconstruction
  THaSpectrometer::Track();
  // TODO
    
  return 0;
    
}

//_____________________________________________________________________________
Int_t SBSBigBite::Reconstruct()
{
  // Fine Reconstruction of particles in spectrometer
  //std::cout << "SBSBigBite::Reconstruct()..." << std::endl;
    
  THaSpectrometer::Reconstruct();
  // TODO
    
  //std::cout << "Done..." << std::endl;
    
  return 0;
    
}

//_____________________________________________________________________________
Int_t SBSBigBite::FindVertices( TClonesArray& tracks )
{
    
  //std::cout << "SBSBigBite::FindVertices()...";
  // Reconstruct target coordinates for all tracks found.
  Int_t n_trk = tracks.GetLast()+1;
  for( Int_t t = 0; t < n_trk; t++ ) {
    auto* theTrack = static_cast<THaTrack*>( tracks.At(t) );
    CalcOpticsCoords(theTrack);


    if(fOpticsOrder>=0){
      CalcTargetCoords(theTrack);
      if(fForwardOpticsOrder>=0){ //also calculate forward from target as consistency check:
	CalcFpCoords( theTrack );
      }
    }
  }
    
  //sort on other criteria than chi^2
  if( GetTrSorting() ) {
    fTracks->Sort();
    // Reassign track indexes. Sorting may have changed the order
    for( int i = 0; i < fTracks->GetLast()+1; i++ ) {
      auto* theTrack = static_cast<THaTrack*>( fTracks->At(i) );
      assert( theTrack );
      theTrack->SetIndex(i);
    }
  }
    
  if( GetNTracks() > 0 ) {
    // Select first track in the array. If there is more than one track
    // and track sorting is enabled, then this is the best fit track
    // (smallest chi2/ndof).  Otherwise, it is the track with the best
    // geometrical match (smallest residuals) between the U/V clusters
    // in the upper and lower VDCs (old behavior).
    //
    // Chi2/dof is a well-defined quantity, and the track selected in this
    // way is immediately physically meaningful. The geometrical match
    // criterion is mathematically less well defined and not usually used
    // in track reconstruction. Hence, chi2 sorting is preferable, albeit
    // obviously slower.


    fGoldenTrack = static_cast<THaTrack*>( fTracks->At(0) );
    fTrkIfo      = *fGoldenTrack;
    fTrk         = fGoldenTrack;
  } else
    fGoldenTrack = nullptr;
    
  //std::cout << "done." << std::endl;
    
  return 0;
}


void SBSBigBite::CalcOpticsCoords( THaTrack* track )
{
  Double_t x_fp, y_fp, xp_fp, yp_fp;
  //Double_t px, py, pz;// NB: not the actual momentum!
    
  x_fp = track->GetX();
  y_fp = track->GetY();
  xp_fp = track->GetTheta();
  yp_fp = track->GetPhi();
    
  TVector3 TrackPosLocal_GEM( x_fp, y_fp, 0.0 );
    
  //std::cout << "calculating optics coordinates: (xfp,yfp,xpfp,ypfp)=(" << x_fp << ", " << y_fp << ", " << xp_fp << ", " << yp_fp << ")" << std::endl;
    
  TVector3 TrackPosGlobal_GEM = fGEMorigin + TrackPosLocal_GEM.X() * fGEMxaxis_global + TrackPosLocal_GEM.Y() * fGEMyaxis_global + TrackPosLocal_GEM.Z() * fGEMzaxis_global;
    
  //std::cout << "Track pos global = " << endl;
  //TrackPosGlobal_GEM.Print();
    
  TVector3 TrackDirLocal_GEM( xp_fp, yp_fp, 1.0 );
  TrackDirLocal_GEM = TrackDirLocal_GEM.Unit();
    
  //  std::cout << "Track direction local = " << endl;
    
  // TrackDirLocal_GEM.Print();
    
  TVector3 TrackDirGlobal_GEM = TrackDirLocal_GEM.X() * fGEMxaxis_global + TrackDirLocal_GEM.Y() * fGEMyaxis_global + TrackDirLocal_GEM.Z() * fGEMzaxis_global;
  TrackDirGlobal_GEM = TrackDirGlobal_GEM.Unit(); //Likely unnecessary, but safer (I guess)
    
  //  std::cout << "Track direction global = " << endl;
  //TrackDirGlobal_GEM.Print();
    
  //Now project track to the z = 0 plane of the ideal optics system:
  // recall the formula to intersect a ray with a plane:
  // (x + s * trackdir - x0) dot planedir = 0.0
    
  double sintersect = (fOpticsOrigin - TrackPosGlobal_GEM).Dot(fOpticsZaxis_global)/ (TrackDirGlobal_GEM.Dot( fOpticsZaxis_global ) );
    
  TVector3 TrackIntersect_global = TrackPosGlobal_GEM + sintersect * TrackDirGlobal_GEM;
    
  //  std::cout << "Track intersection point, global = " << endl;
  //TrackIntersect_global.Print();
    
  //rather than modifying the x, y, theta, phi directly, let's use the RX, RY, RTheta, and RPhi coordinates:
  //TVector3 TrackIntersect_ = TrackIntersect_global - fOpticsOrigin;
    
  double xoptics = (TrackIntersect_global - fOpticsOrigin).Dot( fOpticsXaxis_global );
  double yoptics = (TrackIntersect_global - fOpticsOrigin).Dot( fOpticsYaxis_global );
    
  double xpoptics = TrackDirGlobal_GEM.Dot( fOpticsXaxis_global )/TrackDirGlobal_GEM.Dot( fOpticsZaxis_global );
  double ypoptics = TrackDirGlobal_GEM.Dot( fOpticsYaxis_global )/TrackDirGlobal_GEM.Dot( fOpticsZaxis_global );
    
  //std::cout << "GEM origin = " << std::endl;
  //fGEMorigin.Print();
  //std::cout << "Optics origin = " << std::endl;
  //fOpticsOrigin.Print();
    
  //std::cout << "GEM z axis global = " << std::endl;
  //fGEMzaxis_global.Print();
  //std::cout << "Optics z axis global = " << std::endl;
  //fOpticsZaxis_global.Print();
    
  //std::cout << "GEM (x,y,xp,yp) = " << x_fp << ", " << y_fp << ", " << xp_fp << ", " << yp_fp << std::endl;
  // std::cout << "Optics (x,y,xp,yp) = " << xoptics << ", " << yoptics << ", " << xpoptics << ", " << ypoptics << endl;
    
  track->SetR( xoptics, yoptics, xpoptics, ypoptics );
    
  //The following line is no longer necessary as we are using the "Rotated TRANSPORT coordinates" to store the track parameters in ideal optics system
  //track->Set(x_fp, y_fp, xp_fp, yp_fp);
    
    
}

void SBSBigBite::CalcTargetCoords( THaTrack* track )
{
  //std::cout << "SBSBigBite::CalcTargetCoords()...";
    
  //const double //make it configurable
  const double th_bb = GetThetaGeo();//retrieve the actual angle
    
  TVector3 BB_zaxis( sin(th_bb), 0, cos(th_bb) );
  TVector3 BB_xaxis(0,-1,0);
  TVector3 BB_yaxis = (BB_zaxis.Cross(BB_xaxis)).Unit();
    
  TVector3 spec_xaxis_fp,spec_yaxis_fp, spec_zaxis_fp;
  TVector3 spec_xaxis_tgt,spec_yaxis_tgt, spec_zaxis_tgt;
    
  spec_xaxis_tgt = BB_xaxis;
  spec_yaxis_tgt = BB_yaxis;
  spec_zaxis_tgt = BB_zaxis;
    
  spec_zaxis_fp = BB_zaxis;
  spec_yaxis_fp = BB_yaxis;
  spec_zaxis_fp.Rotate(-fOpticsAngle, spec_yaxis_fp);
  spec_xaxis_fp = spec_yaxis_fp.Cross(spec_zaxis_fp).Unit();
    
  Double_t x_fp, y_fp, xp_fp, yp_fp;
  //if( fCoordType == kTransport ) {
    
  if( track->HasRot() ){
    //    std::cout << "using rotated track coordinates for optics: " << endl;
    x_fp = track->GetRX();
    y_fp = track->GetRY();
    xp_fp = track->GetRTheta();
    yp_fp = track->GetRPhi();
  } else {
    //std::cout << "using non-rotated track coordinates for optics: " << endl;
    x_fp = track->GetX();
    y_fp = track->GetY();
    xp_fp = track->GetTheta();
    yp_fp = track->GetPhi();
  }
  //}
  //cout << x_fp << " " << y_fp << " " << xp_fp << " " << yp_fp << endl;
    
  //Beam position from BPMs and Rasters
  double ybeam=0.0,xbeam=0.0;

  //retrieve beam position from BPMs and Rasters
  TIter aiter(gHaApps);
  THaApparatus* app = 0;
  
  while( (app=(THaApparatus*)aiter()) ){
    if(app->InheritsFrom("SBSRasteredBeam")){
      //We want the Lrb variables
      if(app->GetName() != std::string("Lrb")) continue;
      
      SBSRasteredBeam* RasterBeam = reinterpret_cast<SBSRasteredBeam*>(app);
      ybeam = RasterBeam->GetBeamPosition().Y(); 
      xbeam = RasterBeam->GetBeamPosition().X();
    }
    //cout << app->GetName() << endl;
  }


  double /*vx, vy, vz, */px, py, pz;
  double p_fit, xptar_fit, yptar_fit, ytar_fit, xtar;
  double thetabend_fit;
  double pthetabend_fit;
  double vz_fit;

  //The fUseBeamPosInOptics flag protects against unintentionally using uncalibrated beam position in optics calculations;
  // or mixing beam position corrections with optics calibrated without them in an inconsistent way
  
  xtar = fUseBeamPosInOptics ? -ybeam : 0.0;
  //  if( fUseBeamPosInOptics ) xtar = -ybeam;
 
  //Three iterations needed to converge xtar reconstruction
  for(int iter = 0; iter < 2; iter++){
    int ipar = 0;
    xptar_fit = 0.0;
    yptar_fit = 0.0;
    ytar_fit = 0.0;
    pthetabend_fit = 0.0;
    for(int i=0; i<=fOpticsOrder; i++){
      for(int j=0; j<=fOpticsOrder-i; j++){
	for(int k=0; k<=fOpticsOrder-i-j; k++){
	  for(int l=0; l<=fOpticsOrder-i-j-k; l++){
	    for(int m=0; m<=fOpticsOrder-i-j-k-l; m++){
	      double term = pow(x_fp,m)*pow(y_fp,l)*pow(xp_fp,k)*pow(yp_fp,j)*pow(xtar,i);
	      
	      if( !fDownBendingMode ){
		xptar_fit += fb_xptar[ipar]*term;
		yptar_fit += fb_yptar[ipar]*term;
		ytar_fit += fb_ytar[ipar]*term;
		pthetabend_fit += fb_pinv[ipar]*term;
	      } else {
		xptar_fit += fb_xptar_downbend[ipar]*term;
		yptar_fit += fb_yptar_downbend[ipar]*term;
		ytar_fit += fb_ytar_downbend[ipar]*term;
		pthetabend_fit += fb_pinv_downbend[ipar]*term;
	      }

	      //pinv_fit += b_pinv(ipar)*term;
	      // cout << ipar << " " << term << " "
	      //      << b_xptar(ipar) << " " << b_yptar(ipar) << " "
	      //      << b_ytar(ipar) << " " << b_pinv(ipar) << endl;
	      ipar++;
	    }
	  }
	}
      }
    } // End loop over matrix
  
    vz_fit = -ytar_fit / (sin(th_bb) + cos(th_bb)*yptar_fit);
    
    xtar = -cos(GetThetaGeo()) * vz_fit * xptar_fit;    
    if( fUseBeamPosInOptics ) xtar += -ybeam;
  }
    
  //Let's simplify the bend angle reconstruction to avoid double-counting, even though
  //this calculation is almost certainly correct:
    
  TVector3 phat_tgt_fit(xptar_fit, yptar_fit, 1.0 );
  phat_tgt_fit = phat_tgt_fit.Unit();
    
  TVector3 phat_fp(xp_fp, yp_fp, 1.0 );
  phat_fp = phat_fp.Unit();
    
  TVector3 phat_fp_rot = phat_fp.X() * fOpticsXaxis_global +
    phat_fp.Y() * fOpticsYaxis_global +
    phat_fp.Z() * fOpticsZaxis_global;
    
  thetabend_fit = acos( phat_fp_rot.Dot( phat_tgt_fit ) );
    
  // TVector3 phat_tgt_fit_global = phat_tgt_fit.X() * spec_xaxis_tgt +
  //   phat_tgt_fit.Y() * spec_yaxis_tgt +
  //   phat_tgt_fit.Z() * spec_zaxis_tgt;
    
  // TVector3 phat_fp_fit(xp_fp, yp_fp, 1.0 );
  // phat_fp_fit = phat_fp_fit.Unit();
    
  // TVector3 phat_fp_fit_global = phat_fp_fit.X() * spec_xaxis_fp +
  //   phat_fp_fit.Y() * spec_yaxis_fp +
  //   phat_fp_fit.Z() * spec_zaxis_fp;
    
  //thetabend_fit = acos( phat_fp_fit_global.Dot( phat_tgt_fit_global ) );
  
  if( fPrecon_flag != 1 ){
    p_fit = pthetabend_fit/thetabend_fit;
  } else {
    double delta = pthetabend_fit;
    double p_firstorder = fA_pth1 * ( 1.0 + (fB_pth1 + fC_pth1*fMagDist)*xptar_fit ) / thetabend_fit;

    //If A_vy and B_vy parameters were defined AND the use of beam position corrections in optics is enabled,
    //then correct the momentum:
    p_fit = p_firstorder * (1.0 + delta);
    if( fUseBeamPosInOptics ) p_fit += -(fA_vy + fB_vy * ybeam);
    //p_fit = p_firstorder * (1.0 + delta) - (fA_vy + fB_vy * ybeam);
  }

  pz = p_fit*sqrt( 1.0/(xptar_fit*xptar_fit+yptar_fit*yptar_fit+1.) );
  px = xptar_fit * pz;
  py = yptar_fit * pz; 
 
  TVector3 pvect_BB = TVector3(px, py, pz);
    
  px = +pvect_BB.Z()*sin(th_bb)+pvect_BB.Y()*cos(th_bb);
  py = -pvect_BB.X();
  pz = pvect_BB.Z()*cos(th_bb)-pvect_BB.Y()*sin(th_bb);
    
  //  f_xtg_exp.push_back(xtar);
  
  track->SetTarget(xtar, ytar_fit, xptar_fit, yptar_fit);
  track->SetMomentum(p_fit);
  track->SetPvect(TVector3(px, py, pz));
  track->SetVertex(TVector3(xbeam, ybeam, vz_fit));
    
  //cout << px << " " << py << " " << pz << "   " << vz_fit << endl;
  //cout << track->GetLabPx() << " " << track->GetLabPy() << " " << track->GetLabPz()
  //   << "   " << track->GetVertexZ() << endl;
    
  //std::cout << "Done." << std::endl;
}

//_____________________________________________________________________________
void SBSBigBite::CalcFpCoords( THaTrack *track ){

  if( track->HasTarget() ){
    double xtar = track->GetTX();
    double ytar = track->GetTY();
    double xptar = track->GetTTheta();
    double yptar = track->GetTPhi();
    double p = track->GetP();
        
    double xfp_fit = 0.0;
    double yfp_fit = 0.0;
    double xpfp_fit = 0.0;
    double ypfp_fit = 0.0;

    //xtar = 0.0;
    
    int ipar=0;
    //forward optics expansion is (xfp, yfp, xpfp, ypfp) = sum_ijklm C_(xyxpyp)^ijklm * xptar^i yptar^j ytar^k (1/p)^l (xtar)^m:
    for( int i=0; i<=fForwardOpticsOrder; i++ ){
      for( int j=0; j<=fForwardOpticsOrder-i; j++ ){
	for( int k=0; k<=fForwardOpticsOrder-i-j; k++ ){
	  for( int l=0; l<=fForwardOpticsOrder-i-j-k; l++ ){
	    for( int m=0; m<=fForwardOpticsOrder-i-j-k-l; m++ ){
	      double term = pow( xptar, m )*pow( yptar, l )*pow( ytar, k ) * pow( 1.0/p, j ) * pow( xtar, i );
	      xfp_fit += fb_xfp[ipar] * term;
	      yfp_fit += fb_yfp[ipar] * term;
	      xpfp_fit += fb_xpfp[ipar] * term;
	      ypfp_fit += fb_ypfp[ipar] * term;
	      ipar++;
	    }
	  }
	}
      }
    }
        
    // std::cout << "(xptar, yptar, xtar, ytar, p (GeV) ) = (" << xptar << ", " << yptar << ", " << xtar << ", " << ytar << ", " << p << ")" << std::endl;
        
    // std::cout << "Fitted (xfp, yfp, xpfp, ypfp) = (" << xfp_fit << ", " << yfp_fit << ", " << xpfp_fit << ", " << ypfp_fit << ")" << std::endl;
        
        
    //now the "fit" focal plane coordinates are in the ideal optics system. Should we translate to the "detector" system? YES! because here is the best place to do so:
    TVector3 pos_optics( xfp_fit, yfp_fit, 0.0 );
    TVector3 dir_optics( xpfp_fit, ypfp_fit, 1.0 );
    dir_optics = dir_optics.Unit();
        
    TVector3 posglobal_optics = fOpticsOrigin + pos_optics.X() * fOpticsXaxis_global + pos_optics.Y() * fOpticsYaxis_global + pos_optics.Z() * fOpticsZaxis_global;
        
    TVector3 dirglobal_optics = dir_optics.X() * fOpticsXaxis_global + dir_optics.Y() * fOpticsYaxis_global + dir_optics.Z() * fOpticsZaxis_global;
        
    double sintersect = ( fGEMorigin - posglobal_optics ).Dot( fGEMzaxis_global )/( dirglobal_optics.Dot( fGEMzaxis_global ) );
        
    TVector3 TrackGEMintersect_global = posglobal_optics + sintersect * dirglobal_optics;
        
    double xGEM = (TrackGEMintersect_global - fGEMorigin).Dot( fGEMxaxis_global );
    double yGEM = (TrackGEMintersect_global - fGEMorigin).Dot( fGEMyaxis_global );
        
    double xpGEM = dirglobal_optics.Dot( fGEMxaxis_global )/dirglobal_optics.Dot( fGEMzaxis_global );
    double ypGEM = dirglobal_optics.Dot( fGEMyaxis_global )/dirglobal_optics.Dot( fGEMzaxis_global );
        
    track->SetD( xGEM, yGEM, xpGEM, ypGEM );
        
  }
}

//_____________________________________________________________________________
Int_t SBSBigBite::TrackCalc()
{
  // Additioal track calculations
  // Timing calculation goes here
  // PID calculation goes here!
  for( Int_t t = 0; t < fTracks->GetLast()+1; t++ ) {
    auto* theTrack = static_cast<THaTrack*>( fTracks->At(t) );
    CalcTrackTiming(theTrack);
  }
    
  return 0;
}

//_____________________________________________________________________________
Int_t SBSBigBite::CalcPID()
{
  // PID calculation goes here!
  for( Int_t t = 0; t < fTracks->GetLast()+1; t++ ) {
    auto* theTrack = static_cast<THaTrack*>( fTracks->At(t) );
    CalcTrackPID(theTrack);
  }
  return 0;
}

//_____________________________________________________________________________
void SBSBigBite::CalcTrackTiming(THaTrack* the_track)
{
  TIter next( fNonTrackingDetectors );
  while( auto* theNonTrackDetector =
	 static_cast<THaNonTrackingDetector*>( next() )) {
    // match a hodoscope cluster to a track:
    //the hodoscope has to be found for anything to be done.
    if(theNonTrackDetector->InheritsFrom("SBSTimingHodoscope")){
      SBSTimingHodoscope* TH = reinterpret_cast<SBSTimingHodoscope*>(theNonTrackDetector);
            
      //x, y of track at z = Z_hodoscope
      double x_track = the_track->GetX(TH->GetOrigin().Z());
      //double y_track = the_track->GetY(TH->GetOrigin().Z());
      //cout << x_track << " " << y_track << endl;
      // Not sure what to use for the hodoscope.
      // Perhaps we'd have to complete the class with a cluster
            
      for(int i=0; i<TH->GetNClusters(); i++){
	SBSTimingHodoscopeCluster* clus = TH->GetCluster(i);
	if(clus->GetXmean()-clus->GetSize()*TH->SizeRow()/2<x_track &&
	   x_track<clus->GetXmean()+clus->GetSize()*TH->SizeRow()/2){
	  //std::cout << clus->GetSize() << " " << clus->GetTmean() << " " << clus->GetXmean() << " " << clus->GetYmean() << " " << clus->GetToTmean() << std::endl;
	  the_track->SetTime(clus->GetTmean());
	}
      }
    }//end if inherits from hodoscope
        
        
  }
    
}

//_____________________________________________________________________________
void SBSBigBite::CalcTrackPID(THaTrack* the_track)
{
  if(fEpsEtotRatio.size()==0 || fEtot.size()==0 ||
     fFrontConstraintX.size()==0 || fFrontConstraintY.size()==0 ||
     fBackConstraintX.size()==0 || fBackConstraintY.size()==0)return;
    
  //particles: 0: electron, 1: pion
  //detectors: 0: calo, 1: GRINCH
  auto* pidinfo = the_track->GetPIDinfo();
    
  double eproba, piproba;
    
  TIter next( fNonTrackingDetectors );
}


Int_t SBSBigBite::proba_pssh(Double_t eps_etot_ratio, 
                             Double_t& proba_e, Double_t& proba_pi)
{
  if(fProba_e_PSSH_table.size()==0)return -1;
  proba_e = fProba_e_PSSH_table[fProba_e_PSSH_table.size()-1];
  proba_pi = fProba_pi_PSSH_table[fProba_pi_PSSH_table.size()-1];
  for(size_t i = 0; i<fEpsEtotRatio_table.size()-1; i++){
    if(fEpsEtotRatio_table[i]<=eps_etot_ratio &&
       eps_etot_ratio<fEpsEtotRatio_table[i+1]){
      proba_e = fProba_e_PSSH_table[i]+
	(fProba_e_PSSH_table[i+1]-fProba_e_PSSH_table[i])/
	(fEpsEtotRatio_table[i+1]-fEpsEtotRatio_table[i])*
	(eps_etot_ratio-fEpsEtotRatio_table[i]);
      proba_pi = fProba_pi_PSSH_table[i]+
	(fProba_pi_PSSH_table[i+1]-fProba_pi_PSSH_table[i])/
	(fEpsEtotRatio_table[i+1]-fEpsEtotRatio_table[i])*
	(eps_etot_ratio-fEpsEtotRatio_table[i]);
    }
  }
  return 0;
}

Int_t SBSBigBite::proba_pcal(Double_t etot_p_ratio, 
                             Double_t& proba_e, Double_t& proba_pi)
{
  if(fEtotPratio_table.size()==0)return -1;
  proba_e = fProba_e_PCAL_table[fProba_e_PCAL_table.size()-1];
  proba_pi = fProba_e_PCAL_table[fProba_pi_PCAL_table.size()-1];
  for(size_t i = 0; i<fEtotPratio_table.size()-1; i++){
    if(fEtotPratio_table[i]<=etot_p_ratio && etot_p_ratio<fEtotPratio_table[i+1]){
      proba_e = fProba_e_PCAL_table[i]+
	(fProba_e_PCAL_table[i+1]-fProba_e_PCAL_table[i])/
	(fEtotPratio_table[i+1]-fEtotPratio_table[i])*
	(etot_p_ratio-fEtotPratio_table[i]);
      proba_pi = fProba_pi_PCAL_table[i]+
	(fProba_pi_PCAL_table[i+1]-fProba_pi_PCAL_table[i])/
	(fEtotPratio_table[i+1]-fEtotPratio_table[i])*
	(etot_p_ratio-fEtotPratio_table[i]);
    }
  }
  return 0;
}

Int_t SBSBigBite::proba_grinch(Int_t npmt, Double_t p,
                               Double_t& proba_e, Double_t& proba_pi)
{
  if(fProba_e_GRINCH_table.size()==0)return -1;
  int j = fP_table.size()-1;
  if(j==-1)return -1;
  for(size_t i = 0; i<fP_table.size()-1;i++){
    if(fP_table[i]<p && p<fP_table[i])j = i;
  }
  proba_e = fProba_e_GRINCH_table[fProba_e_GRINCH_table.size()-1];
  proba_pi = fProba_e_GRINCH_table[fProba_pi_GRINCH_table.size()-1];
  for(size_t i = 0; i<fNGRINCHPMTs_table.size()-1; i++){
    if(fNGRINCHPMTs_table[i]<=npmt && npmt<fNGRINCHPMTs_table[i+1]){
      proba_e = fProba_e_GRINCH_table[i]+
	(fProba_e_GRINCH_table[i+1]-fProba_e_GRINCH_table[i])/
	(fNGRINCHPMTs_table[i+1]-fNGRINCHPMTs_table[i])*
	(npmt-fNGRINCHPMTs_table[i]);
      proba_pi = fProba_pi_GRINCH_table[j][i]+
	(fProba_pi_GRINCH_table[j][i+1]-fProba_pi_GRINCH_table[j][i])/
	(fNGRINCHPMTs_table[i+1]-fNGRINCHPMTs_table[i])*
	(npmt-fNGRINCHPMTs_table[i]);
    }
  }
  return 0;
}

//_____________________________________________________________________________
Bool_t SBSBigBite::SetTrSorting( Bool_t set )
{
  Bool_t oldset = TestBit(kSortTracks);
  SetBit( kSortTracks, set );
  return oldset;
}

//_____________________________________________________________________________
Bool_t SBSBigBite::GetTrSorting() const
{
  return TestBit(kSortTracks);
}

//_____________________________________________________________________________
Bool_t SBSBigBite::SetMultiTracks( Bool_t set )
{
  Bool_t oldset = TestBit(kMultiTracks);
  SetBit( kMultiTracks, set );
  return oldset;
}

//_____________________________________________________________________________
Bool_t SBSBigBite::GetMultiTracks() const
{
  return TestBit(kMultiTracks);
}

void SBSBigBite::InitOpticsAxes(double BendAngle, const TVector3 &Origin ){
  fOpticsOrigin = Origin;
  fOpticsYaxis_global.SetXYZ(0,1,0);
  fOpticsZaxis_global.SetXYZ(-sin(BendAngle), 0, cos(BendAngle) );
  fOpticsXaxis_global.SetXYZ(cos(BendAngle), 0, sin(BendAngle) );
}

void SBSBigBite::InitOpticsAxes(double BendAngle ){
  // fOpticsOrigin = Origin;
  fOpticsYaxis_global.SetXYZ(0,1,0);
  fOpticsZaxis_global.SetXYZ(-sin(BendAngle), 0, cos(BendAngle) );
  fOpticsXaxis_global.SetXYZ(cos(BendAngle), 0, sin(BendAngle) );
}

void SBSBigBite::InitGEMAxes( double theta, double phi, const TVector3 &Origin ){
  fGEMorigin = Origin;
  fGEMzaxis_global.SetXYZ( sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta) );
  fGEMxaxis_global = (fOpticsYaxis_global.Cross( fGEMzaxis_global) ).Unit(); // check to make sure this is consistent with definition in the zero-field alignment code
  fGEMyaxis_global = (fGEMzaxis_global.Cross(fGEMxaxis_global)).Unit();
}

void SBSBigBite::InitGEMAxes( double theta, double phi ){
  //fGEMorigin = Origin;
  fGEMzaxis_global.SetXYZ( sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta) );
  fGEMxaxis_global = (fOpticsYaxis_global.Cross( fGEMzaxis_global) ).Unit(); // check to make sure this is consistent with definition in the zero-field alignment code
  fGEMyaxis_global = (fGEMzaxis_global.Cross(fGEMxaxis_global)).Unit();
}



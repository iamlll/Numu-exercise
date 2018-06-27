//Making a contour plot at 90% confidence, 3sigma, and 5sigma for a nu_muon disappearance simulation with Gaussian smearing.

#include <iostream>
#include <fstream>
#include <math.h>
#include <random>
using namespace std;

Double_t iPOT = 6e20; //protons-on-target in Icarus
Double_t sPOT = 6e20;//protons-on-target in SBND

const Int_t SBND_dist = 100; //m, dist from SBND to the source
const Int_t Icarus_dist = 600; //m, distance from Icarus to the source

const Int_t activeSBND = 112; //112 tons active mass of LAr in SBND
const Int_t activeIc = 476; //476 tons active mass of LAr in Icarus

TCanvas c1("c1", "Plots", 1000, 500);
TLegend SBNDLeg(.7, 0.6, .9, 0.8, "Legend"); //legend for SBND events
TLegend legend(0.1,0.2,.3,.4);

TH1D* hflux;
Double_t mean, SD; //mean and standard deviation of flux histogram
TGraph* graph;
Int_t nbinsx, size; //# bins in provided histogram + # points in provided graph

void readFiles(){
   //Divide canvas  
   c1.Divide(2,1);
 
   //retrieving histogram from file
   TFile* flux = new TFile("flux.root");
   flux -> ls();
   hflux = (TH1D*)flux -> Get("numu_CV_AV_TPC");

   //find # bins in hist
   nbinsx = hflux->GetXaxis()->GetNbins();
   
   //Read file, pull graph from file, and draw TGraph
   TFile* cross = new TFile("cross_section.root");
   cross -> ls();
   graph = (TGraph*)cross -> Get("qel_cc_n");
}

/**
   Find the flux data given by the histogram in flux.root file; store data in array.
   @return an array containing the flux data from a histogram.
*/ 
Double_t* findFlux(){ 
   Double_t* fluxes = new Double_t[nbinsx];

   //transfer event counts into array
   for(int i=0; i<nbinsx; i++){
      Double_t n = hflux->GetBinContent(i+1); //not including underflow bin 0
      fluxes[i]=n;
   }
   return fluxes;
}

/**
   Finds energies corresponding to each bin, in GeV.
   @return array of energy values 
*/
Double_t* findEnergies(){
   Double_t* nrgs = new Double_t[nbinsx];
   for(int i=0; i<nbinsx; i++){
      Double_t en = hflux->GetXaxis()->GetBinCenter(i+1); //again, ignoring the underflow bin 0
      nrgs[i] = en;
   }
   return nrgs;
}

//Find cross-sections of neutrino events
Double_t* findCrossSections(){
   Double_t* areas = new Double_t[nbinsx];
   //Get total # points in graph
   size = graph->GetN();

   //Copy array of all y-values (AKA graph->GetY()) in plot to array 'allareas'
   //3rd param of 'memcpy' is total # bytes to copy, so do (# elements)*(bytes per element of the array given by GetY())
   Double_t allareas[size];
   memcpy(allareas, graph->GetY(), size*sizeof(*(graph->GetY())));

   //Take every 5th data point to cut down info to same # 'bins' as flux hist
   //start on 3rd data pt to take the data from the middle of each 5-data point cluster ("bin")
   for(int i = 2; i < size; i+=5){
      areas[(i-2)/5]=allareas[i];
   }
   return areas;
}

//Find "actual" event counts by multiplying data by event rate
TH1D* countEvents(Double_t energy[], Double_t fluxes[], Double_t areas[], Int_t distance, Int_t activeMass, Double_t PoT, Int_t draw){
   c1.cd(1);
   TH1D* hEvents = new TH1D("hEvents","Events", nbinsx, energy[0], energy[nbinsx-1]);
   Double_t* events = new Double_t[nbinsx];
   //units of cross-section: 10e-38 cm^2/Ar = 10e-38/100/100 m^2/Ar
   //units of flux: N/[m^2*(50MeV)*(10^6 PoT)] 
   
   //total # Argon atoms in detector: 
   auto totAr = activeMass * 1e6/39.948 * 6.02e23; //39.948 g/mol Ar; 6.02e23 atoms/mol 
   auto scalingfactor = pow(470./distance,2) * 1e-42 * totAr * (PoT*1e-6); //need scaling factor for distance (data normalized at 470 m from source); 1e-4 converting cm^2 to m^2; PoT=protons on-target     

   for(int i=0; i < nbinsx; i++){
      events[i] = fluxes[i]*areas[i] * scalingfactor ; //left with units of: (# neutrinos interacting)/50 MeV bin
      hEvents->SetBinContent(i+1, events[i]);
   }

   if(draw!=121){
      hEvents->Draw();
      SBNDLeg.AddEntry(hEvents, "original", "p");
   }

   return hEvents;
}

/**
   Convolve the neutrino event distribution with a Gaussian smearing (resolution function, resolution of user's choice). Save the convoluted values into an array of smeared event rates.

   @param resolution the energy resolution of muon neutrinos in the simulation.
   @return a smeared Gaussian with the given energy resolution.
*/
Double_t* convolution(Double_t resolution, Int_t color, Int_t option, TH1D* hEvents){
   TH1D *hconv = (TH1D*)hEvents->Clone("hconv");
   Double_t* smearedEvents = new Double_t[nbinsx];
   hconv->Reset();

   if(option != 121){
      if(color==1){
         hconv->SetLineColor(kRed);
         SBNDLeg.AddEntry(hconv, "5% res", "lp");
      }
      else if(color==2){
         hconv->SetLineColor(kYellow);
         SBNDLeg.AddEntry(hconv, "10% res", "lp");
      }
      else if(color==3){
         hconv->SetLineColor(kBlue);
         SBNDLeg.AddEntry(hconv, "25% res", "lp");
         SBNDLeg.Draw();
      }
   }

   auto binWidth = hEvents->GetXaxis()->GetBinWidth(1);
   for(Int_t i=1; i<nbinsx+1; i++){
      Double_t total=0.0;
      for(Int_t j=1; j<nbinsx+1; j++){
         auto xi=hEvents->GetXaxis()->GetBinCenter(i);
         auto xj=hEvents->GetXaxis()->GetBinCenter(j);
         auto yj=hEvents->GetBinContent(j);
         auto sigma = resolution*sqrt(xj);
         total = total + yj*ROOT::Math::gaussian_pdf(xi, sigma, xj);
      }  
      auto value = total*binWidth;
      hconv->SetBinContent(i, value);
      smearedEvents[i]=value;
   }
   if(option!=121) hconv->Draw("SAME");
   return smearedEvents;
}   


//Draw contours for 90% conf level, 3sigma, and 5sigma in same plot
void drawContours(Double_t **Chi, Int_t numX, Int_t numY, Double_t xs[], Double_t ys[], Int_t option){
   //Contour value for 5*sigma = 23.40
   c1.cd(2);
   legend.SetHeader("Legend", "C");   

   TH2D *hCont = new TH2D("hCont","Smeared 5*sigma Contour",numX-1, xs, numY-1, ys); 
   hCont->SetTitle("Smeared 5*sigma Contour; sin(2theta_mumu)^2; dm_41^2 (eV^2)");
   for(Int_t i=0; i< numX; i++){
      for(Int_t j=0; j<numY; j++){
         if(i!=0 && i!=numX-1 && j!=0 && j!=numY-1){
/*
             if((Chi[i][j+1]>23.4 && Chi[i][j-1]<23.4)||(Chi[i+1][j]>23.4 && Chi[i-1][j]<23.4)){ 
                if(option==1){
                   hCont->Fill(xs[i],ys[j],1); hCont->SetMarkerColor(46);
                }
                else if(option==2){ 
                   hCont->Fill(xs[i],ys[j],1); hCont->SetMarkerColor(42);
                }
                else if(option==3){
                   hCont->Fill(xs[i],ys[j],1); hCont->SetMarkerColor(38);
                }
             }
*/
             if(Chi[i][j]>=23.0 & Chi[i][j] <=24.0){
                if(option==1){
                   hCont->Fill(xs[i],ys[j],1); hCont->SetMarkerColor(46);
                }
                else if(option==2){
                   hCont->Fill(xs[i],ys[j],1); hCont->SetMarkerColor(42);
                }
                else if(option==3){
                   hCont->Fill(xs[i],ys[j],1); hCont->SetMarkerColor(38);
                }

             }
         }         
      }
   }
   gPad->SetLogy(); gPad->SetLogx();
   gStyle->SetPalette(1);

   if(option==1){
      legend.AddEntry(hCont, "5% res", "lp");
      hCont->Draw("COLZ");
   }
   else if(option==2){
      legend.AddEntry(hCont, "10% res", "lp");
      hCont->Draw("SAME");
   }   
   else if(option==3){
      legend.AddEntry(hCont, "25% res", "lp");
      hCont->Draw("SAME");
      legend.Draw();
   }
}

//Find and store all Chi-squared values for Icarus ONLY
Double_t** chiSq(Double_t nrgs[], Double_t farEvents[], Double_t nearEvents[], Int_t option){

   Double_t spaceM = 0.005;
   Double_t spaceA = 0.005;
   Int_t No_mass = Int_t((2.0-(-2.0))/spaceM); cout << "# mass entries: " << No_mass << endl;
   Int_t No_ang = Int_t(2.0/spaceA); cout << "# angle entries: " << No_ang << endl;
   Double_t** Chi2 = 0;
   Chi2 = new Double_t*[No_ang];
 
   Double_t xs[No_ang];
   Double_t ys[No_mass];
   
   //Loop through all plausible dm^2 and sin(2theta_mumu)^2 values
   
   Int_t xCount = 0; Int_t yCount = 0;   

   for(Double_t ang2 = pow(10,-2); ang2 < 1.0; ang2= ang2*pow(10,spaceA)){
      xs[xCount] = ang2; 
      xCount++;
   }
   for(Double_t dm2 = pow(10, -2); dm2 < pow(10,2); dm2= dm2*pow(10,spaceM)){
      ys[yCount] = dm2;  
      yCount++;
   }
   Double_t minChi=10.0; //minimum Chi-squared value   
 
   for(Int_t i=0; i < No_ang; i++){
      Chi2[i] = new Double_t[No_mass];
      for(Int_t j=0; j<No_mass;j++){ 
         auto Chi2sum = 0.0;

         for(int n=0; n < nbinsx; n++){
            if(farEvents[n]==0 || nearEvents[n]==0) {continue;}

/*taking hbar and c into account, from Boris Kayser's paper,

* survival probability P_{numu->numu} = 1-sin(2theta_mumu)^2*sin(1.27*dm^2(eV^2)*L(km)/E_nu(GeV))^2
*/
            auto prob = 1 - xs[i]*pow(sin(1.27*ys[j]*0.001*Icarus_dist/nrgs[n]),2); //for Icarus

/* We have E_flux=(0.3*N)^2, but this is an uncertainty that shifts the data from SBND and Icarus 
* by the same amount, 30% up and down (and so "cancels"). What we really care about now, is the 
* counting error propagated from SBND to Icarus - assuming 100% correlation between the near and
* far detector, we only know our measurements from the far detector as well as we knew it in the
* near detector, i.e. sqrt(N_SBND). However, this error relies not on the exact number of events
* seen from SBND, but rather the proportion of the error per bin: sqrt(N_100)/N_100. 
* Scaling this value to that of Icarus:
* E_propagated = [sqrt(N_100)/N_100]^2*(N_600/N_100). Thus,
* E_tot (matrix values) = [sqrt(N_100)/N_100]^2*(N_600/N_100) + sqrt(N_600)^2,
* since Icarus also has its own counting error to deal with.
* Chi^2 = Sum(N_i^{null} - N_i^{osc})^2*1/E_tot
*/
            auto totErr = pow((sqrt(nearEvents[n])/nearEvents[n]),2)*(farEvents[n]/nearEvents[n]) + farEvents[n];
            auto t1 = pow(totErr, -1);
            auto t2 = pow(farEvents[n]-prob*farEvents[n], 2);      
            Chi2sum += t1*t2;
         }
         Chi2[i][j] = Chi2sum;
         if(i==0 && j==0){ minChi = Chi2[i][j]; }
         if(Chi2sum < minChi) minChi = Chi2sum; //finding minimum Chi-sq value
      } 
   }

//update array of Chi-sq values and fill the histogram
   for(int nSine=0; nSine < No_ang; nSine++){
      for(int nMass=0; nMass < No_mass; nMass++){
         Chi2[nSine][nMass] = Chi2[nSine][nMass]-minChi;
      }
   }

   drawContours(Chi2, No_ang, No_mass, xs, ys, option);   

   return Chi2;
}

Int_t goose(){
   readFiles();
   Double_t* flux = findFlux();
   Double_t* energy = findEnergies();    
   Double_t* cross = findCrossSections();
   TH1D* nullSBND = countEvents(energy, flux, cross, SBND_dist, activeSBND, sPOT, 1);
   TH1D* nullIcarus = countEvents(energy, flux, cross, Icarus_dist, activeIc, iPOT, 121);

//   Double_t* SEvents5 = convolution(0.05, 1, 1, nullSBND);
//   Double_t* SEvents10 = convolution(0.1, 2, 1, nullSBND);
   Double_t* SEvents25 = convolution(0.25, 3, 1, nullSBND);
//   Double_t* IEvents5 = convolution(0.05, 1, 121, nullIcarus);
//   Double_t* IEvents10 = convolution(0.10, 2, 121, nullIcarus);
   Double_t* IEvents25 = convolution(0.25, 3, 121, nullIcarus);

//   chiSq(energy, IEvents5, SEvents5, 1);
//   chiSq(energy, IEvents10, SEvents10, 2);
   chiSq(energy, IEvents25, SEvents25, 3);

   return 0;
}

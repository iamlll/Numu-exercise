#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

const Int_t Icarus_dist = 600; //m, distance from Icarus to the source
TCanvas c1("c1", "Plots", 900, 450);
TH1D* h1;
TGraph* graph;
Int_t nbinsx, size; //# bins in provided histogram + # points in provided graph

ofstream myfile2a;

void readFiles(){
   //Split canvas
   c1.Divide(2,1);
   
   //retrieving histogram from file
   TFile* flux = new TFile("flux.root");
   flux -> ls();
   h1 = (TH1D*)flux -> Get("numu_CV_AV_TPC");
   
   //find # bins in hist
   nbinsx = h1->GetXaxis()->GetNbins();

   //Read file, pull graph from file, and draw TGraph
   TFile* cross = new TFile("cross_section.root");
   cross -> ls();
   graph = (TGraph*)cross -> Get("qel_cc_n");
}

//Find the flux data given by the histogram in flux.root file; store data in arrays and write results to new file.

Double_t* findFlux(){ 
   Double_t* fluxes = new Double_t[nbinsx];
   
//   cout << "Bin #\t Flux (N/[m^2*(50MeV)*(10^6 PoT)])" << "Energy (GeV)\n";
   
   //transfer event counts into array 
   for(int i=0; i<nbinsx; i++){
      Double_t n = h1->GetBinContent(i+1); //not including underflow bin 0
      fluxes[i]=n;  
//      cout << i << "\t" << fluxes[i] << "\n";
   }
   return fluxes; 
}

//finds energies corresponding to each bin, in GeV.
Double_t* findEnergies(){
   Double_t* nrgs = new Double_t[nbinsx];
   for(int i=0; i<nbinsx; i++){
      Double_t en = h1->GetXaxis()->GetBinCenter(i+1); //again, ignoring the underflow bin 0
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
Double_t* countEvents(Double_t fluxes[], Double_t areas[]){
   Double_t* events = new Double_t[nbinsx];
   //units of cross-section: 10e-38 cm^2/Ar = 10e-38/100/100 m^2/Ar
   //units of flux: N/[m^2*(50MeV)*(10^6 PoT)] 
   
   //total # Argon atoms in detector: 
   auto totAr = 476 * 1e6/39.948 * 6.02e23; //476 tons active mass of LAr in Icarus; 39.948 g/mol Ar; 6.02e23 atoms/mol 
   auto scalingfactor = (476./89)*pow(470./Icarus_dist,2) * 1e-42 * totAr * (6e20*1e-6); 
   //(476/89) = mass normalization factor; (470/600)^2 = scaling factor for distance (data normalized at 470 m from source); 1e-4 converting cm^2 to m^2; 6*10^20 protons on-target     

   for(int i=0; i < nbinsx; i++){
      events[i] = fluxes[i]*areas[i] * scalingfactor ; //left with units of: (# neutrinos interacting)/50 MeV
   }
   return events;
}

//Draw contours for 90% conf level, 3sigma, and 5sigma in separate plot
void drawContours(Double_t **Chi, Int_t numX, Int_t numY, Double_t xs[], Double_t ys[]){
   //Contour values = {1.64, 7.74, 23.40}, corresponding to dX^2_90, dX^2_3sigma, and dX^2_5sigma, respectively.
   
   TH2D *hCont = new TH2D("hCont","ChiSquared ContLines",numX-1, xs, numY-1, ys);
   hCont->SetTitle("ChiSquared Contour; sin(2theta_mumu)^2; dm_41^2 (eV^2)");
   for(Int_t i=0; i< numX; i++){
      for(Int_t j=0; j<numY; j++){
         if(i!=0 && i!=numX-1 && j!=0 && j!=numY-1){
             if((Chi[i][j+1]>1.64 && Chi[i][j-1]<1.64)||(Chi[i+1][j]>1.64 && Chi[i-1][j]<1.64)){ hCont->Fill(xs[i],ys[j],1.64); }
             if((Chi[i][j+1]>7.74 && Chi[i][j-1]<7.74)||(Chi[i+1][j]>7.74 && Chi[i-1][j]<7.74)){ hCont->Fill(xs[i], ys[j], 7.74); }
             if((Chi[i][j+1]>23.4 && Chi[i][j-1]<23.4)||(Chi[i+1][j]>23.4 && Chi[i-1][j]<23.4)){ hCont->Fill(xs[i],ys[j],23.4); }            
         }         
      }
   }
   c1.cd(2);
   gPad->SetLogy(); gPad->SetLogx(); gPad->SetLogz();
   gStyle->SetPalette(1);
   hCont->Draw("COLZ");      
}

//Find and store all Chi-squared values from data
Double_t** chiSq(Double_t nrgs[], Double_t events[]){
   c1.cd(1);

   Double_t spaceM = 0.005;
   Double_t spaceA = 0.005;
   Int_t No_mass = Int_t((2.0-(-2.0))/spaceM); cout << "# mass entries: " << No_mass << endl;
   Int_t No_ang = Int_t(3.0/spaceA); cout << "# angle entries: " << No_ang << endl;
   Double_t** Chi2 = 0;
   Chi2 = new Double_t*[No_ang];
 
   Double_t xs[No_ang];
   Double_t ys[No_mass];

   myfile2a.open("ChiSquareData2a.txt");
   
   //Loop through all plausible dm^2 and sin(2theta_mumu)^2 values
   
   myfile2a << "Entry #:\tMass^2:\tAngle^2:\tChi^2:\n";
   
   Int_t xCount = 0; Int_t yCount = 0;   

   for(Double_t ang2 = pow(10,-3); ang2 < 1.0; ang2= ang2*pow(10,spaceA)){
      xs[xCount] = ang2; 
      xCount++;
   }
   for(Double_t dm2 = pow(10, -2); dm2 < pow(10,2); dm2= dm2*pow(10,spaceM)){
      ys[yCount] = dm2;  
      yCount++;
   }

   TH2D *h2 = new TH2D("h2","ChiSquared RegCont",No_ang-1, xs, No_mass-1, ys);
   h2->SetTitle("ChiSquared Contour; sin(2theta_mumu)^2; dm_41^2 (eV^2)");
  
   Double_t minChi=10.0; //minimum Chi-squared value   
 
   for(Int_t i=0; i < No_ang; i++){
      Chi2[i] = new Double_t[No_mass];
      for(Int_t j=0; j<No_mass;j++){ 
         auto Chi2sum = 0.0;

         for(int n=0; n < nbinsx; n++){
            if(events[n]==0) {continue;}

/*taking hbar and c into account, from Boris Kayser's paper,
* probability P_{numu->numu} = 1-sin(2theta_mumu)^2*sin(1.27*dm^2(eV^2)*L(km)/E_nu(GeV))^2
*/
            auto prob = 1 - xs[i]*pow(sin(1.27*ys[j]*0.001*Icarus_dist/nrgs[n]),2);

/* We have E_flux=(0.3*N)^2, and still retain E_stat=(sqrt(N))^2=N; thus,
* E_tot = N+(0.3N)^2)
* New Chi^2 = Sum(N_i^{null} - N_i^{osc})^2*1/E_tot
*/
            auto totErr = events[n]+pow(0.3*events[n],2);
            auto t1 = pow(totErr, -1);
            auto t2 = pow(events[n]-prob*events[n], 2);      
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
         h2->Fill(xs[nSine], ys[nMass], Chi2[nSine][nMass]);
         myfile2a << "[" << nSine << "][" << nMass << "]\t" << xs[nSine] << "\t" << ys[nMass] << "\t" << Chi2[nSine][nMass] << endl;
      }
   }

   myfile2a.close();
   
   gPad->SetLogy(); gPad->SetLogx(); gPad->SetLogz(); 
   gStyle->SetPalette(1);
   h2->Draw("COLZ");   
   
   drawContours(Chi2, No_ang, No_mass, xs, ys);   

   return Chi2;
}

int uncertainty(){
   readFiles();
   Double_t* flux = findFlux();
   Double_t* energy = findEnergies();    
   Double_t* cross = findCrossSections();
   Double_t* events = countEvents(flux, cross);
   chiSq(energy, events);
   return 0;
}

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>

#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>

#include <string>
#include <fstream>
#include <iostream>
#include <ostream>
#include <sstream>

#include "Evento.h"

using namespace std;

int main(int argc, char *argv[]) 
{

  int status = 0, iL=0, NumProcessed=0, numfiles = 0, nread=0;
  TString outname("Out.root"), inname("In.txt");
  vector<TString> infiles; TString tmpSin;

  ifstream lista_file;
  char fname[200],linea[200];
  FILE *pfile;
  bool ReadError = false;
  int kSelWord=0, kSelDc = 0, kSelLyso = 0, maxevpro = 1000000000;
  int majority=0, iflag_analog=0;
  double Ethreshold = 0;

  static TTree *RootTree = 0;
  static Evento *pEv = 0;

  EVENTO_STRUCT eve;

  ifstream infile;


  for (int i = 0; i < argc; i++){
    if(strcmp(argv[i],"-in") == 0) {
      inname = TString(argv[++i]);
    }
    if(strcmp(argv[i],"-out") == 0) {
      outname = TString(argv[++i]);
    }
    if(strcmp(argv[i],"-nev") == 0) {
      maxevpro = atoi(argv[++i]);
    }
    if(strcmp(argv[i],"-sel") == 0) {
      kSelWord = atoi(argv[++i]);
      kSelLyso = (int)(((float)kSelWord)/10);
      kSelDc = (int)(kSelWord - 10*kSelLyso);
      cout<<endl<<"Txt2Root: Selection word = "<<kSelWord<<" SelLyso= "<<kSelLyso<<" SelDc= "<<kSelDc<<endl<<endl;
    }
    if(strcmp(argv[i],"-iL") == 0) { iL = 1; }
    if(strcmp(argv[i],"-help") == 0) {
      cout<<"Conversion of fluka TXT file : usage -> Txt2Root [opts] "<<endl;
      cout<<" possible opts are:"<<endl;
      cout<<"   -in  file    : [def=In.txt] Root input file"<<endl;
      cout<<"   -out  file   : [def=Out.root] Root output file"<<endl;
      cout<<"   -sel selw    : [def=0] select ev: 1*dc + 10*lyso "<<endl;
      cout<<"   -iL        : [def=none] input file is a list of files"<<endl;
      cout<<"   -nev        : [def=Inf] Max no. of events to process"<<endl;
      return 1;
    }
  }

  if(iL==0) {
    infiles.push_back(inname);
  } else {
    lista_file.open(inname.Data());
    cout<<"Processing data from "<<inname.Data()<<" LIST file!"<<endl;
    while (lista_file.getline(linea, 200, '\n')) {
      sscanf(linea,"%s",fname);
      cout<<"Adding "<<fname<<" to the list of files to be processed!"<<endl;
      tmpSin= fname;
      infiles.push_back(tmpSin);
    }
    lista_file.close();
  }
  numfiles = infiles.size();

  TFile *f_out = new TFile(outname,"RECREATE");
  f_out->cd();


  RootTree = new TTree("EventTree","gsimay");

  RootTree->Branch("nevent",&eve.num_event,"nevent/I");
  RootTree->Branch("nump",&eve.nump,"nump/I");
  RootTree->Branch("idpa",&eve.idpa,"idpa[nump]/I");
  RootTree->Branch("igen",&eve.igen,"igen[nump]/I");
  RootTree->Branch("icha",&eve.icha,"icha[nump]/I");
  RootTree->Branch("numreg",&eve.numreg,"numreg[nump]/I");
  RootTree->Branch("iba",&eve.iba,"iba[nump]/I");
  RootTree->Branch("idead",&eve.idead,"idead[nump]/I");
  RootTree->Branch("jpa",&eve.jpa,"jpa[nump]/I");
  RootTree->Branch("vxi",&eve.vxi,"vxi[nump]/D");
  RootTree->Branch("vyi",&eve.vyi,"vyi[nump]/D");
  RootTree->Branch("vzi",&eve.vzi,"vzi[nump]/D");
  RootTree->Branch("vxf",&eve.vxf,"vxf[nump]/D");
  RootTree->Branch("vyf",&eve.vyf,"vyf[nump]/D");
  RootTree->Branch("vzf",&eve.vzf,"vzf[nump]/D");
  RootTree->Branch("px",&eve.px,"px[nump]/D");
  RootTree->Branch("py",&eve.py,"py[nump]/D");
  RootTree->Branch("pz",&eve.pz,"pz[nump]/D");
  RootTree->Branch("pxf",&eve.pxf,"pxf[nump]/D");
  RootTree->Branch("pyf",&eve.pyf,"pyf[nump]/D");
  RootTree->Branch("pzf",&eve.pzf,"pzf[nump]/D");
  RootTree->Branch("amass",&eve.amass,"amass[nump]/D");
  RootTree->Branch("tempo",&eve.tempo,"tempo[nump]/D");
  RootTree->Branch("tof",&eve.tof,"tof[nump]/D");
  RootTree->Branch("trlen",&eve.trlen,"trlen[nump]/D");

  
  RootTree->Branch("ncross",&eve.ncross,"ncross/I");
  RootTree->Branch("idcross",&eve.idcross,"idcross[ncross]/I");
  RootTree->Branch("nregcross",&eve.nregcross,"nregcross[ncross]/I");
  RootTree->Branch("nregold",&eve.nregold,"nregold[ncross]/I");
  RootTree->Branch("xcross",&eve.xcross,"xcross[ncross]/D");
  RootTree->Branch("ycross",&eve.ycross,"ycross[ncross]/D");
  RootTree->Branch("zcross",&eve.zcross,"zcross[ncross]/D");
  RootTree->Branch("pxcross",&eve.pxcross,"pxcross[ncross]/D");
  RootTree->Branch("pycross",&eve.pycross,"pycross[ncross]/D");
  RootTree->Branch("pzcross",&eve.pzcross,"pzcross[ncross]/D");
  RootTree->Branch("mcross",&eve.mcross,"mcross[ncross]/D");
  RootTree->Branch("chcross",&eve.chcross,"chcross[ncross]/D");
  RootTree->Branch("tcross",&eve.tcross,"tcross[ncross]/D");

  pEv = new Evento();

  //    loop sui file della lista ( if any)


  for(int idfile=0; idfile<numfiles;idfile++){
    cout<<endl<<"Now processing data from "<<infiles.at(idfile)<<" file!"<<endl;
    ReadError = false;
    
    pfile = fopen(infiles.at(idfile),"r");
    /*    
    nread= fscanf(pfile,"%d %lf %d\n",&majority,&Ethreshold,&iflag_analog);        
    if(nread!=3){
      ReadError= true;
      cout<<"Wrong run header: read Error!! exiting  "<<endl;
      return 1;
    }
    */
    //	  loop sugli eventi del file  

    while((!feof(pfile))&&(!ReadError)){
      NumProcessed++;
      if (maxevpro>0) {
	if (NumProcessed>maxevpro) break;
      }
      eve.num_event = 0;
      eve.nump = 0 ;
      eve.ncross = 0 ;

      if(NumProcessed%10000==0){ cout<<"# event = "<<NumProcessed<<endl;}
      status = pEv->Clean();

      //	leggo l'header

      nread= fscanf(pfile,"%d %d %d \n",&eve.num_event,&eve.nump,
		   &eve.ncross);
      //    cout<<"header = "<<nread<<" "<<eve.num_event<<" "<<eve.nump<<" "<<eve.ncross<<endl;
      if(nread!=3){
	cout<<"ReadError in ev header section: nread = "<<nread<<
	  " instead of 6; ev= "<<NumProcessed<<endl;
	ReadError = true;
      }
      
      //	leggo il common della cinematica

      if(!ReadError){
	for(int jj =0;jj<eve.nump;jj++){
	  nread = fscanf(pfile,
			 "%d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n ",
			 &eve.idpa[jj],&eve.igen[jj],&eve.icha[jj],
			 &eve.numreg[jj],&eve.iba[jj],&eve.idead[jj],
			 &eve.jpa[jj],&eve.vxi[jj],
			 &eve.vyi[jj],&eve.vzi[jj],&eve.vxf[jj],&eve.vyf[jj],
			 &eve.vzf[jj],&eve.px[jj],&eve.py[jj],&eve.pz[jj],
			 &eve.pxf[jj],&eve.pyf[jj],&eve.pzf[jj],&eve.amass[jj],
			 &eve.tempo[jj],&eve.tof[jj],&eve.trlen[jj]);
	  /*
	  printf("%d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n ",
			 eve.idpa[jj],eve.igen[jj],eve.icha[jj],
 			 eve.numreg[jj],eve.iba[jj],eve.idead[jj],
			 eve.jpa[jj],eve.vxi[jj],
			 eve.vyi[jj],eve.vzi[jj],eve.vxf[jj],eve.vyf[jj],
			 eve.vzf[jj],eve.px[jj],eve.py[jj],eve.pz[jj],
			 eve.pxf[jj],eve.pyf[jj],eve.pzf[jj],eve.amass[jj],
			 eve.tempo[jj],eve.tof[jj],eve.trlen[jj]);
	  */
	  if(nread!=23){
	    cout<<"ReadError in kine section: nread = "<<nread<<
	      " instead of 23; ev= "<<NumProcessed<<endl;
	    
	    ReadError= true;
	    break;
	  }
	}
      }


      //	leggo i boundary crossing

      if(!ReadError){
	for(int jj=0; jj<eve.ncross;jj++){
	  nread = fscanf(pfile,"%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf \n",
			 &eve.idcross[jj],&eve.nregcross[jj],
			 &eve.nregold[jj],&eve.xcross[jj],
			 &eve.ycross[jj],&eve.zcross[jj],&eve.pxcross[jj],
			 &eve.pycross[jj],&eve.pzcross[jj],&eve.mcross[jj],
			 &eve.chcross[jj],&eve.tcross[jj]);
	  /*
	  printf("%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf \n",
			 eve.idcross[jj],eve.nregcross[jj],eve.xcross[jj],
			 eve.ycross[jj],eve.zcross[jj],eve.pxcross[jj],
			 eve.pycross[jj],eve.pzcross[jj],eve.mcross[jj],
			 eve.chcross[jj],eve.tcross[jj]);
	  */
	  if(nread!=12){
	    cout<<"ReadError in cross section: nread = "<<nread<<
	      " instead of 12; ev= "<<NumProcessed<<endl;
	    ReadError = true;
	    break;
	  }
	}
      }
      
      // end of event data structure reading
      
      if((!ReadError) &&(eve.nump<=MAXNUMP)&&(eve.ncross<=MAXCROSS)){
	RootTree->Fill() ;
      }
      else{
	cout<<"ReadError= "<<ReadError<<" um_event= "<<eve.num_event<<" nump= "<<eve.nump<<" ncross= "<<eve.ncross<<endl;  
      }
    }
     
    fclose(pfile);
  }
  RootTree->Write();
  
  f_out->Close();
  cout<<" total number of event safely converted= "<<NumProcessed<<endl;

  return 0;
  
}

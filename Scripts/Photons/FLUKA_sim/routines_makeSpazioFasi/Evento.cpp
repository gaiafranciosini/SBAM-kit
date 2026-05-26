#include <stdio.h>
#include <iostream>
#include <fstream>
#include "Evento.h"

using namespace std;

ClassImp(Evento);

/*-----------------------------------------------------------------*/
  
Evento::Evento()
{
  eve.num_event = 0;
  eve.nump = 0;
  eve.ncross = 0;

  for(int kk=0;kk<MAXNUMP;kk++){
    eve.idpa[kk] = 0;
    eve.igen[kk] = 0;
    eve.icha[kk] = -100;
    eve.numreg[kk] = 0;
    eve.iba[kk] = 0;
    eve.idead[kk] = 0;
    eve.jpa[kk] = -100;
    eve.vxi[kk] = 0.;
    eve.vyi[kk] = 0.;
    eve.vzi[kk] = 0.;
    eve.vxf[kk] = 0.;
    eve.vyf[kk] = 0.;
    eve.vzf[kk] = 0.;
    eve.px[kk] = 0.;
    eve.py[kk] = 0.;
    eve.pz[kk] = 0.;
    eve.pxf[kk] = 0.;
    eve.pyf[kk] = 0.;
    eve.pzf[kk] = 0.;
    eve.amass[kk] = 0.;
    eve.tempo[kk] = 0.;
    eve.tof[kk] = 0.;
    eve.trlen[kk] = 0.;
  }

  for(int kk=0;kk<MAXCROSS;kk++){
    eve.idcross[kk]=0;
    eve.nregcross[kk]=0;
    eve.nregold[kk]=0;
    eve.xcross[kk]=0.;
    eve.ycross[kk]=0.;
    eve.zcross[kk]=0.;
    eve.pxcross[kk]=0.;
    eve.pycross[kk]=0.;
    eve.pzcross[kk]=0.;
    eve.mcross[kk]=0;
    eve.chcross[kk]=0.;
    eve.tcross[kk]=0.;
  }
  
}

/*-----------------------------------------------------------------*/

Int_t Evento::SetEvent(Int_t fnumero_evento){
  //  cout <<" Entro in SetEvent"<<endl;
  eve.num_event = fnumero_evento;
  return 0;
}

/*-----------------------------------------------------------------*/

Int_t Evento::Clean(){
  //  cout <<" Entro in clean"<<endl;

  for(int kk=0;kk<eve.nump;kk++){
    eve.idpa[kk] = 0;
    eve.igen[kk] = 0;
    eve.icha[kk] = -100;
    eve.numreg[kk] = 0;
    eve.iba[kk] = 0;
    eve.idead[kk] = 0;
    eve.jpa[kk] = -100;
    eve.vxi[kk] = 0.;
    eve.vyi[kk] = 0.;
    eve.vzi[kk] = 0.;
    eve.vxf[kk] = 0.;
    eve.vyf[kk] = 0.;
    eve.vzf[kk] = 0.;
    eve.px[kk] = 0.;
    eve.py[kk] = 0.;
    eve.pz[kk] = 0.;
    eve.pxf[kk] = 0.;
    eve.pyf[kk] = 0.;
    eve.pzf[kk] = 0.;
    eve.amass[kk] = 0.;
    eve.tempo[kk] = 0.;
    eve.tof[kk] = 0.;
    eve.trlen[kk] = 0.;
  }


  for(int kk=0;kk<eve.ncross;kk++){
    eve.idcross[kk]=0;
    eve.nregcross[kk]=0;
    eve.nregold[kk]=0;
    eve.xcross[kk]=0.;
    eve.ycross[kk]=0.;
    eve.zcross[kk]=0.;
    eve.pxcross[kk]=0.;
    eve.pycross[kk]=0.;
    eve.pzcross[kk]=0.;
    eve.mcross[kk]=0;
    eve.chcross[kk]=0.;
    eve.tcross[kk]=0.;
  }

  eve.num_event = 0;
  eve.nump = 0;
  eve.ncross = 0;
  return 0;
}

/*-----------------------------------------------------------------*/
/*   twin function to be called before and together AddPx_After */

Int_t Evento::AddPart(Int_t fidpa, Int_t figen, 
		      Int_t ficha, Int_t fnumreg, Int_t fiba,
		      Int_t fidead, Int_t fjpa,
		      Double_t fvxi, Double_t fvyi, Double_t fvzi,
		      Double_t fvxf, Double_t fvyf, Double_t fvzf,
		      Double_t fpx, Double_t fpy, Double_t fpz,
		      Double_t fpxf, Double_t fpyf, Double_t fpzf,
		      Double_t famass, Double_t ftempo,
		      Double_t ftof, Double_t ftrlen){
 //  cout<<"Entro in Evento::AddPart "<<endl;
  if(eve.nump<MAXNUMP)
    {
      eve.nump ++;
      eve.idpa[eve.nump-1] = fidpa;
      eve.igen[eve.nump-1] = figen;
      eve.icha[eve.nump-1] = ficha;
      eve.numreg[eve.nump-1] = fnumreg;
      eve.iba[eve.nump-1] = fiba;
      eve.idead[eve.nump-1] = fidead;
      eve.jpa[eve.nump-1] = fjpa;
      eve.vxi[eve.nump-1] = fvxi;
      eve.vyi[eve.nump-1] = fvyi;
      eve.vzi[eve.nump-1] = fvzi;
      eve.vxf[eve.nump-1] = fvxf;
      eve.vyf[eve.nump-1] = fvyf;
      eve.vzf[eve.nump-1] = fvzf;
      eve.px[eve.nump-1] = fpx;
      eve.py[eve.nump-1] = fpy;
      eve.pz[eve.nump-1] = fpz;
      eve.pxf[eve.nump-1] = fpxf;
      eve.pyf[eve.nump-1] = fpyf;
      eve.pzf[eve.nump-1] = fpzf;
      eve.amass[eve.nump-1] = famass;
      eve.tempo[eve.nump-1] = ftempo;
      eve.tof[eve.nump-1] = ftof;
      eve.trlen[eve.nump-1] = ftrlen;
      return 0;
    }
  else
    {
      return -1;
    }
}


/*-----------------------------------------------------------------*/

Int_t Evento::AddCross(Int_t fidcross, Int_t fnregcross, Int_t fnregold,
		       Double_t fxcross, Double_t fycross, 
		       Double_t fzcross, Double_t fpxcross, 
		       Double_t fpycross, Double_t fpzcross, 
		       Double_t fmcross, Double_t fchcross, 
		       Double_t ftcross){
  //  cout<<"Entro in Evento::AddCross "<<endl;
  if(eve.ncross<MAXCROSS) {
    eve.ncross ++;
    eve.idcross[eve.ncross-1] = fidcross;
    eve.nregcross[eve.ncross-1] = fnregcross;
    eve.nregold[eve.ncross-1] = fnregold;
    eve.xcross[eve.ncross-1] = fxcross;
    eve.ycross[eve.ncross-1] = fycross;
    eve.zcross[eve.ncross-1] = fzcross;
    eve.pxcross[eve.ncross-1] = fpxcross;
    eve.pycross[eve.ncross-1] = fpycross;
    eve.pzcross[eve.ncross-1] = fpzcross;
    eve.mcross[eve.ncross-1] = fmcross;
    eve.chcross[eve.ncross-1] = fchcross;
    eve.tcross[eve.ncross-1] = ftcross;
    return 0;
  }
  else{
    return -1;
  }
}


/*-----------------------------------------------------------------*/
int Evento::FindBranches(TTree *RootTree, EVENTO_STRUCT *eve){

  RootTree->SetBranchAddress("nevent",&(eve->num_event));

  RootTree->SetBranchAddress("nump",&(eve->nump));
  RootTree->SetBranchAddress("idpa",&(eve->idpa));
  RootTree->SetBranchAddress("igen",&(eve->igen));
  RootTree->SetBranchAddress("icha",&(eve->icha));
  RootTree->SetBranchAddress("numreg",&(eve->numreg));
  RootTree->SetBranchAddress("iba",&(eve->iba));
  RootTree->SetBranchAddress("idead",&(eve->idead));
  RootTree->SetBranchAddress("jpa",&(eve->jpa));
  RootTree->SetBranchAddress("vxi",&(eve->vxi));
  RootTree->SetBranchAddress("vyi",&(eve->vyi));
  RootTree->SetBranchAddress("vzi",&(eve->vzi));
  RootTree->SetBranchAddress("vxf",&(eve->vxf));
  RootTree->SetBranchAddress("vyf",&(eve->vyf));
  RootTree->SetBranchAddress("vzf",&(eve->vzf));
  RootTree->SetBranchAddress("px",&(eve->px));
  RootTree->SetBranchAddress("py",&(eve->py));
  RootTree->SetBranchAddress("pz",&(eve->pz));
  RootTree->SetBranchAddress("pxf",&(eve->pxf));
  RootTree->SetBranchAddress("pyf",&(eve->pyf));
  RootTree->SetBranchAddress("pzf",&(eve->pzf));
  RootTree->SetBranchAddress("amass",&(eve->amass));
  RootTree->SetBranchAddress("tempo",&(eve->tempo));
  RootTree->SetBranchAddress("tof",&(eve->tof));
  RootTree->SetBranchAddress("trlen",&(eve->trlen));
  /*    */
  /**/
  RootTree->SetBranchAddress("ncross",&(eve->ncross));
  RootTree->SetBranchAddress("idcross",&(eve->idcross));
  RootTree->SetBranchAddress("nregcross",&(eve->nregcross));
  RootTree->SetBranchAddress("nregold",&(eve->nregold));
  RootTree->SetBranchAddress("xcross",&(eve->xcross));
  RootTree->SetBranchAddress("ycross",&(eve->ycross));
  RootTree->SetBranchAddress("zcross",&(eve->zcross));
  RootTree->SetBranchAddress("pxcross",&(eve->pxcross));
  RootTree->SetBranchAddress("pycross",&(eve->pycross));
  RootTree->SetBranchAddress("pzcross",&(eve->pzcross));
  RootTree->SetBranchAddress("mcross",&(eve->mcross));
  RootTree->SetBranchAddress("chcross",&(eve->chcross));
  RootTree->SetBranchAddress("tcross",&(eve->tcross));
                   /*    */
                   /*    */
    return 0;
}


/*-----------------------------------------------------------------*/


Int_t Evento::Dump(){
  return 0;

}

/*-----------------------------------------------------------------*/


EVENTO_STRUCT Evento::Output(){
  return eve;
}


/*-----------------------------------------------------------------*/

Evento::~Evento()
{
}



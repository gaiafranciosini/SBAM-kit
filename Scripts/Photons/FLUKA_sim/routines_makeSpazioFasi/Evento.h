#ifndef ROOT_Evento
#define ROOT_Evento

#include <iostream>
#include <fstream>
#include "TObject.h"
#include "TTree.h"

#include "EventStruct.h"

using namespace std;


class Evento : public TObject {

 public:

  Evento();
  virtual ~Evento();

  Int_t Clean();

  Int_t SetEvent(Int_t fnumero_evento);

  Int_t AddPart(Int_t fidpa, Int_t figen,
  Int_t ficha, Int_t fnumreg, Int_t fiba,
  Int_t fidead, Int_t fjpa,
  Double_t fvxi, Double_t fvyi, Double_t fvzi,
  Double_t fvxf, Double_t fvyf, Double_t fvzf,
  Double_t fpx, Double_t fpy, Double_t fpz,
  Double_t fpxf, Double_t fpyf, Double_t fpzf,
  Double_t famass, Double_t ftempo,
  Double_t ftof, Double_t ftrlen);


  Int_t AddCross(Int_t fidcross, Int_t fnreggcross, Int_t fnregold,
   Double_t fxcross, Double_t fycross, Double_t fzcross,
   Double_t fpxcross, Double_t fpycross, Double_t fpzcross,
   Double_t fmcross, Double_t fchcross,Double_t ftcross);

  int FindBranches(TTree *RootTree,EVENTO_STRUCT *eve);

  Int_t Dump();


  EVENTO_STRUCT Output();

 private:

  EVENTO_STRUCT eve;

 protected:

  ClassDef(Evento,1);

};

#endif

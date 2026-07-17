#ifndef _EVENT_STRUCT_HH
#define _EVENT_STRUCT_HH


const int MAXCROSS = 10000;
const int MAXNUMP = 1000;

typedef struct {
  Int_t num_event;
  Int_t nump;
  Int_t idpa[MAXNUMP];
  Int_t igen[MAXNUMP];
  Int_t icha[MAXNUMP];
  Int_t numreg[MAXNUMP];
  Int_t iba[MAXNUMP];
  Int_t idead[MAXNUMP];
  Int_t jpa[MAXNUMP];
  Double_t vxi[MAXNUMP];
  Double_t vyi[MAXNUMP];
  Double_t vzi[MAXNUMP];
  Double_t vxf[MAXNUMP];
  Double_t vyf[MAXNUMP];
  Double_t vzf[MAXNUMP];
  Double_t px[MAXNUMP];
  Double_t py[MAXNUMP];
  Double_t pz[MAXNUMP];
  Double_t pxf[MAXNUMP];
  Double_t pyf[MAXNUMP];
  Double_t pzf[MAXNUMP];
  Double_t amass[MAXNUMP];
  Double_t tempo[MAXNUMP];
  Double_t tof[MAXNUMP];
  Double_t trlen[MAXNUMP];


  Int_t ncross;
  Int_t idcross[MAXCROSS];
  Int_t nregcross[MAXCROSS];
  Int_t nregold[MAXCROSS];
  Double_t xcross[MAXCROSS];
  Double_t ycross[MAXCROSS];
  Double_t zcross[MAXCROSS];
  Double_t pxcross[MAXCROSS];
  Double_t pycross[MAXCROSS];
  Double_t pzcross[MAXCROSS];
  Double_t mcross[MAXCROSS];
  Double_t chcross[MAXCROSS];
  Double_t tcross[MAXCROSS];

} EVENTO_STRUCT;


#endif

// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME AnaFlukaDict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "ROOT/RConfig.hxx"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Header files passed as explicit arguments
#include "Evento.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_Evento(void *p = nullptr);
   static void *newArray_Evento(Long_t size, void *p);
   static void delete_Evento(void *p);
   static void deleteArray_Evento(void *p);
   static void destruct_Evento(void *p);
   static void streamer_Evento(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Evento*)
   {
      ::Evento *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::Evento >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("Evento", ::Evento::Class_Version(), "Evento.h", 14,
                  typeid(::Evento), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::Evento::Dictionary, isa_proxy, 16,
                  sizeof(::Evento) );
      instance.SetNew(&new_Evento);
      instance.SetNewArray(&newArray_Evento);
      instance.SetDelete(&delete_Evento);
      instance.SetDeleteArray(&deleteArray_Evento);
      instance.SetDestructor(&destruct_Evento);
      instance.SetStreamerFunc(&streamer_Evento);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Evento*)
   {
      return GenerateInitInstanceLocal(static_cast<::Evento*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Evento*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr Evento::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *Evento::Class_Name()
{
   return "Evento";
}

//______________________________________________________________________________
const char *Evento::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Evento*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int Evento::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Evento*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Evento::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Evento*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Evento::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Evento*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void Evento::Streamer(TBuffer &R__b)
{
   // Stream an object of class Evento.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      R__b.StreamObject(&(eve),typeid(eve));
      R__b.CheckByteCount(R__s, R__c, Evento::IsA());
   } else {
      R__c = R__b.WriteVersion(Evento::IsA(), kTRUE);
      TObject::Streamer(R__b);
      R__b.StreamObject(&(eve),typeid(eve));
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_Evento(void *p) {
      return  p ? new(p) ::Evento : new ::Evento;
   }
   static void *newArray_Evento(Long_t nElements, void *p) {
      return p ? new(p) ::Evento[nElements] : new ::Evento[nElements];
   }
   // Wrapper around operator delete
   static void delete_Evento(void *p) {
      delete (static_cast<::Evento*>(p));
   }
   static void deleteArray_Evento(void *p) {
      delete [] (static_cast<::Evento*>(p));
   }
   static void destruct_Evento(void *p) {
      typedef ::Evento current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_Evento(TBuffer &buf, void *obj) {
      ((::Evento*)obj)->::Evento::Streamer(buf);
   }
} // end of namespace ROOT for class ::Evento

namespace {
  void TriggerDictionaryInitialization_AnaFlukaDict_Impl() {
    static const char* headers[] = {
"Evento.h",
nullptr
    };
    static const char* includePaths[] = {
"/NFS_homes/software/root_latest/include/",
"/loc_work2/PHOTONS_FLASH/Fasi/routines/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "AnaFlukaDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$Evento.h")))  Evento;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "AnaFlukaDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "Evento.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"Evento", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("AnaFlukaDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_AnaFlukaDict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_AnaFlukaDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_AnaFlukaDict() {
  TriggerDictionaryInitialization_AnaFlukaDict_Impl();
}

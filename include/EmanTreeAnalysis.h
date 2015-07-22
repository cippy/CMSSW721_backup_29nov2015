#ifndef EmanTree_Analysis_h
#define EmanTree_Analysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

#include <iostream>
#include <cstdlib>
#include <vector>

#include "edimarcoTree.h"
//#include "edimarcoTreeFriend.h"
#include "edimarcoTree_noSkim.h"
#include "edimarcoTree_v2.h"

namespace myAnalyzerTEman {

  class znunujetsAna : public edimarcoTree_v2 /*,public edimarcoTreeFriend*/ {
  public:

    znunujetsAna(TTree *tree);
    virtual ~znunujetsAna() { std::cout<<"~znunujetsAna() called"<<std::endl; }
  
    void loop(const char* configFileName);

  };

  class zmumujetsAna : public edimarcoTree /*,public edimarcoTreeFriend*/ {
  public:

    zmumujetsAna(TTree *tree);
    virtual ~zmumujetsAna() { std::cout<<"~zmumujetsAna() called"<<std::endl; }
  
    void loop();

  };

  class zmumujetsAna_LepSk : public edimarcoTree /*,public edimarcoTreeFriend*/ {
  public:

    zmumujetsAna_LepSk(TTree *tree);
    virtual ~zmumujetsAna_LepSk() { std::cout<<"~zmumujetsAna_LepSk() called"<<std::endl; }
  
    void loop();

  };

  class zeejetsAna : public edimarcoTree /*,public edimarcoTreeFriend*/ {
  public:

    zeejetsAna(TTree *tree);
    virtual ~zeejetsAna() { std::cout<<"~zeejetsAna() called"<<std::endl; }
  
    void loop();

  };

  class zmumujetsControlSample : public edimarcoTree /*,public edimarcoTreeFriend*/ {
  public:

    zmumujetsControlSample(TTree *tree, const char*);
    virtual ~zmumujetsControlSample() { std::cout<<"~zmumujetsControlSample() called"<<std::endl; }
  
    void loop(std::vector< Double_t > &, std::vector< Double_t > &);
    std::string suffix;
  };

  class zeejetsControlSample : public edimarcoTree /*,public edimarcoTreeFriend*/ {
  public:

    zeejetsControlSample(TTree *tree, const char*);
    virtual ~zeejetsControlSample() { std::cout<<"~zeejetsControlSample() called"<<std::endl; }
  
    void loop(std::vector< Double_t > &, std::vector< Double_t > &);
    std::string suffix;
  };
  
  class zlljetsAna : public edimarcoTree /*,public edimarcoTreeFriend*/ {
  public:

    zlljetsAna(TTree *tree);
    virtual ~zlljetsAna() { std::cout<<"~zlljetsAna() called"<<std::endl; }
  
    void loop();

  };

  class zlljetsAna_new : public edimarcoTree_v2 /*,public edimarcoTreeFriend*/ {
  public:

    zlljetsAna_new(TTree *tree);
    virtual ~zlljetsAna_new() { std::cout<<"~zlljetsAna_new() called"<<std::endl; }
  
    void loop(const char* configFileName, const Int_t ISDATA_FLAG);

  };

  class zlljets_Axe_noSkim : public edimarcoTree_v2 /*,public edimarcoTreeFriend*/ {
  public:

    zlljets_Axe_noSkim(TTree *tree);
    virtual ~zlljets_Axe_noSkim() { std::cout<<"~zlljets_Axe_noSkim() called"<<std::endl; }
  
    void loop(const char* configFileName);

  };

  class zlljets_Axe_noSkim_light : public edimarcoTree_v2 /*,public edimarcoTreeFriend*/ {
  public:

    // same as previous, but only computes acceptance and efficiency in order to be faster
    zlljets_Axe_noSkim_light(TTree *tree);
    virtual ~zlljets_Axe_noSkim_light() { std::cout<<"~zlljets_Axe_noSkim_light() called"<<std::endl; }
  
    void loop(const char* configFileName);

  };

  class zlljets_resoResp : public edimarcoTree_v2 /*,public edimarcoTreeFriend*/ {
  public:

    // same as previous, but only computes acceptance and efficiency in order to be faster
    zlljets_resoResp(TTree *tree);
    virtual ~zlljets_resoResp() { std::cout<<"~zlljets_resoResp() called"<<std::endl; }
  
    void loop(const char* configFileName, const Int_t ISDATA_FLAG);

  };

}
#endif





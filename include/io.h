//--------------------------------------------------------------------------
/**
 * @file io.h
 * @brief data in out class (text and binary)
 * @ingroup YAMADA
 * @author  M.Yamada
 * @date    Sat Jun 11 22:09:45 2013
 */
//--------------------------------------------------------------------------
#ifndef IO_H_YAMADA_2013_11
#define IO_H_YAMADA_2013_11

#include <string>
#include <fstream>
#include <iostream>
#include <complex>
#include <stdlib.h>

class IODATA {

public:
IODATA();
template <typename DATA>
void callData(DATA in, int line, std::string path, std::string staticInfo, std::string physInfo, int iconf, int it);
template <typename DATA>  // when you use complex, you can use direct darraySize(number of complex elements)
void outData(DATA out, std::string path, std::string staticInfo, std::string physInfo, int iconf, int it ,int arraySize);
 template <typename XDATA ,typename DATA>   // when you use complex, you can use direct darraySize(number of complex elements)
void outErr(XDATA ix, DATA ave, DATA err,std::string path ,std::string staticInfo ,std::string physInfo ,int iconf ,int it ,int arraySize);
void setReadBinaryMode(bool swich);
void setWriteBinaryMode(bool swich);
void setConfSize(int confSize);
void setOutAve(bool swich);
void setD(std::string dim);
  
 private:
  template<typename DATA> int readData(DATA in ,char* inFileName);
  int readData(std::complex<double>* data ,char* inFileName);   //In text mode, complex<double> a.real() a.imag() <==can read ! ||| complex<double> (real,imag) <==type data can't read
  template<typename DATA> int writeData(DATA out,char* outFileName);
  template<typename XDATA ,typename DATA> int writeData1D(XDATA xd, DATA ave,DATA err,char* outFileName);
  template<typename DATA> int writeData3D(DATA ave,DATA err,char* outFileName);
    int writeData(std::complex<double>* data,char* outFileName);
  template<typename DATA> int writeData1D(DATA xd, std::complex<double>* ave,std::complex<double>* err,char* outFileName);
    int writeData3D( std::complex<double>* ave,std::complex<double>* err,char* outFileName);
    void checkConfsize();
    
    
    bool r_bSwitch_;
    bool w_bSwitch_;
    bool wAveSwitch_;
    char inFileName[500];
    char outFileName[500];
    std::string flg_dim_;
    int dataSizeOut_;
    int binnumber_;
    int line_;
};

IODATA::IODATA(){
  binnumber_ = -9999;
  dataSizeOut_ = 0;
  r_bSwitch_ = true;
  w_bSwitch_ = true;
  wAveSwitch_ = false;
  flg_dim_ = "1D";
  line_ = -9999;
}




void IODATA::setReadBinaryMode(bool swich){
  r_bSwitch_ = swich;
}
void IODATA::setWriteBinaryMode(bool swich){
  w_bSwitch_ = swich;
}
void IODATA::setConfSize(int confSize){
  binnumber_ = confSize;
}
void IODATA::setOutAve(bool swich){
  wAveSwitch_ = swich;  
}
void IODATA::setD(std::string flg_dim){
  flg_dim_ = flg_dim;
}

void IODATA::checkConfsize(){
  if(binnumber_ == -9999){std::cerr<<"ERROR confSize is not define :please set confSize using setConfSize(int confSize) "<<std::endl;}
}

template<typename DATA> void IODATA::callData(DATA in , int line, std::string inPath ,std::string staticInfo ,std::string physInfo , int iconf, int it){
  checkConfsize();
  line_ = line -1;
  if(r_bSwitch_)sprintf(inFileName,"%s/%s.%06d-%06d.%s.it%02d",inPath.c_str(),staticInfo.c_str(),binnumber_,iconf,physInfo.c_str(),it);
  else if(!r_bSwitch_)sprintf(inFileName,"%s/%s.%06d-%06d.%s.it%02d",inPath.c_str(),staticInfo.c_str(),binnumber_,iconf,physInfo.c_str(),it);
  else {std::cerr << "ERROR binary swich is warang::"<<std::endl;
    exit(1);}
  readData(in,inFileName);
}

template<typename DATA> void IODATA::outData(DATA out ,std::string outPath ,std::string staticInfo ,std::string physInfo ,int iconf ,int it,int arraySize){
  checkConfsize();
  dataSizeOut_ = arraySize;
  if(w_bSwitch_)sprintf(outFileName,"%s/%s.%06d-%06d.%s.it%02d",outPath.c_str(),staticInfo.c_str(),binnumber_,iconf,physInfo.c_str(),it);
  else if(!w_bSwitch_)sprintf(outFileName,"%s/%s.%06d-%06d.%s.it%02d",outPath.c_str(),staticInfo.c_str(),binnumber_,iconf,physInfo.c_str(),it);
    else {std::cerr << "ERROR binary swich is warang::"<<std::endl;
      exit(1);}
  writeData(out,outFileName);
}

template<typename XDATA,typename DATA> void IODATA::outErr(XDATA ix, DATA out ,DATA err ,std::string outPath ,std::string staticInfo ,std::string physInfo ,int iconf ,int it,int arraySize){
  checkConfsize();
  dataSizeOut_ = arraySize;
  if(w_bSwitch_) {sprintf(outFileName,"%s/%s.%06d-%06d.%s.it%02d.bin",outPath.c_str(),staticInfo.c_str(),binnumber_,iconf,physInfo.c_str(),it);
    if(flg_dim_ == "1D"){ writeData1D(ix,out,err,outFileName);}
    else if(flg_dim_ == "3D"){ writeData3D(out,err,outFileName);}
    else {
      std::cerr<<"ERROR should set Dimension setD(string 1D or 3D)"<<std::endl;
      exit(1);}
  }
  else if(!w_bSwitch_) {sprintf(outFileName,"%s/%s.%06d-%06d.%s.it%02d",outPath.c_str(),staticInfo.c_str(),binnumber_,iconf,physInfo.c_str(),it);
    if(flg_dim_ == "1D"){ writeData1D(ix,out,err,outFileName);}
    else if(flg_dim_ == "3D"){ writeData3D(out,err,outFileName);}
    else {
      std::cerr<<"ERROR should set Dimension setD(string 1D or 3D)"<<std::endl;
      exit(1);}
  }
  
  else {
    std::cerr << "ERROR binary swich is warang::"<<std::endl;
    exit(1);
}
}

template<typename DATA> int IODATA::readData(DATA data ,char* inFileName){
  std::fstream infile;
  if(r_bSwitch_){
    infile.open(inFileName,std::ios::in|std::ios::binary);
    if (!infile.is_open()) {
      std::cerr << "ERROR file can't open (no exist)::"<<inFileName<<std::endl;
      exit(1);
      return EXIT_FAILURE;
    }
    if (infile.fail()){
      std::cerr << "ERROR file size is 0 (can open)::"<<inFileName<<std::endl;
      exit(1);
      return EXIT_FAILURE;
    }
    int id = 0;
    while(!infile.eof()){
      infile.read( ( char * ) &(data[id]), sizeof( DATA ) );
      id=id+1;
    }
    static int tmp = 0;
    if (tmp==0) {
      std::cerr <<"reading binary data size is_open;;"<<id<<std::endl;
      tmp=tmp+1;
    }
    infile.close();
  }
  
  if(!r_bSwitch_){
    int id = 0;
    infile.open(inFileName,std::ios::in);
    if (!infile.is_open()) {
      std::cerr << "ERROR file can't open (no exist)::"<<inFileName<<std::endl;
      exit(1);
      return EXIT_FAILURE;
    }
    if (infile.fail()){
      std::cerr << "ERROR file size is 0 (can open)::"<<inFileName<<std::endl;
      exit(1);
      return EXIT_FAILURE;
    }
    std::string tmpdata;
    do{
      for (int tmpid =0;tmpid <line_;tmpid++){
	infile>>tmpdata;
      }
      infile>>data[id];
      int iChar = 0;
      while ((iChar = infile.get()) != '\n' && iChar != EOF) {}
      if (id>10000000) {
	std::cerr << "ERROR maybe imput file all number 0::"<<inFileName<<std::endl;
	exit(1);
      }
      id = id + 1;
    }
    while ( !infile.eof() ); // important for all 0 file                                                                                                                                            
    id = id -1;
    static int tmp = 0;
    if (tmp==0) {
      std::cout <<"reading text type data size is_open;;"<<id<<" line;;"<<line_ + 1<<std::endl;
      tmp=tmp+1;
    }
    infile.close();
  }
  return 0;
}

int IODATA::readData(std::complex<double>* data ,char* inFileName){
  std::ifstream infile;
  if(r_bSwitch_){
    infile.open(inFileName,std::ios::in|std::ios::binary);
    if (!infile.is_open()) {
      std::cerr << "ERROR file can't open (no exist)::"<<inFileName<<std::endl;
      exit(1);
      return EXIT_FAILURE;
    }
    if (infile.fail()){
      std::cerr << "ERROR file size is 0 (can open)::"<<inFileName<<std::endl;
      exit(1);
      return EXIT_FAILURE;
    }
    int id = 0;
    while(infile.read( ( char * ) &(data[id]), sizeof( double )*2 )){ // important for all 0 file
      id=id+1;
    }
    static int tmp = 0;
    if (tmp==0) {
      std::cout <<"reading binary type data size is_open;;"<<id<<std::endl;
      tmp=tmp+1;
    }
    infile.close();
  }
  
  if(!r_bSwitch_){
    int id = 0;
    infile.open(inFileName,std::ios::in);
    if (!infile.is_open()) {
      std::cerr << "ERROR file can't open (no exist)::"<<inFileName<<std::endl;
      exit(1);
      return EXIT_FAILURE;
    }
    if (infile.fail()){
      std::cerr << "ERROR file size is 0 (can open)::"<<inFileName<<std::endl;
      exit(1);
      return EXIT_FAILURE;
    }
    std::string tmpdata;
    do{
      double real;
      double imag;
      for (int tmpid =0;tmpid <line_;tmpid++){
	infile>>tmpdata;
      }
      
      infile >> real >> imag;
      data[id] = std::complex< double >(real,imag);
      int iChar = 0;
      while ((iChar = infile.get()) != '\n' && iChar != EOF) {}
      if (id>10000000) {
	std::cerr << "ERROR maybe imput file all number 0::"<<inFileName<<std::endl;
	exit(1);
      }
      id = id + 1;
    }
    while ( !infile.eof() ); // important for all 0 file                                                                                                                                            
    id = id -1;
    static int tmp = 0;
    if (tmp==0) {
      std::cout <<"reading text type data size is_open;;"<<id<<" line;;"<<line_ + 1<<std::endl;
      tmp=tmp+1;
    }
    infile.close();
  }
  return 0;
}

template<typename DATA> int IODATA::writeData(DATA data,char* outFileName)
{
  std::fstream ofs;
  if(w_bSwitch_){
    ofs.setf(std::ios::scientific);
    ofs.open(outFileName,std::ios::out|std::ios::binary|std::ios::trunc);
    if (!ofs.is_open()) {
      std::cout << "ERROR output file can't open (no exist)::"<<outFileName<<std::endl;
      exit(1);
      return EXIT_FAILURE;
    }
    for(int id = 0 ;id<dataSizeOut_;id++){
      ofs.write((const char*) &(data[id]),sizeof( DATA ));
    }
    ofs.close();
  }
  if(!w_bSwitch_){
    ofs.setf(std::ios::scientific);
    ofs.open(outFileName,std::ios::out);
    if (!ofs.is_open()) {
      std::cerr << "ERROR output file can't open (no exist)::"<<outFileName<<std::endl;
      exit(1);
      return EXIT_FAILURE;
    }
    for(int id = 0 ;id<dataSizeOut_;id++){
      ofs<<id<<" "<<data[id]<<std::endl;
    }
    ofs.close();}
  return 0;
}

int IODATA::writeData(std::complex<double>* data,char* outFileName)
{
  std::fstream ofs;
  if(w_bSwitch_){
    ofs.setf(std::ios::scientific);
    ofs.open(outFileName,std::ios::out|std::ios::binary|std::ios::trunc);
    if (!ofs.is_open()) {
      std::cout << "ERROR output file can't open (no exist)::"<<outFileName<<std::endl;
      exit(1);
      return EXIT_FAILURE;
    }
    for(int id = 0 ;id<dataSizeOut_;id++){
      ofs.write((const char*) &(data[id]),sizeof( std::complex<double> ));
    }
    ofs.close();
  }
  if(!w_bSwitch_){
    ofs.setf(std::ios::scientific);
    ofs.open(outFileName,std::ios::out);
    if (!ofs.is_open()) {
      std::cerr << "ERROR output file can't open (no exist)::"<<outFileName<<std::endl;
      exit(1);
      return EXIT_FAILURE;
    }
    for(int id = 0 ;id<dataSizeOut_;id++){
      ofs<<id<<" "<<data[id].real()<<" "<<data[id].imag()<<std::endl;
      
    }
    ofs.close();}
  return 0;
}





template<typename XDATA ,typename DATA> int IODATA::writeData1D(XDATA xd ,DATA ave,DATA err,char* outFileName)
{
  std::fstream ofs;
  if(w_bSwitch_){
    ofs.close();
    std::cout << "error output should use text mode "<<std::endl;
    exit(1);
    return EXIT_FAILURE;
  }
    
  if(!w_bSwitch_){
    ofs.setf(std::ios::scientific);
    ofs.open(outFileName,std::ios::out);
    if (!ofs.is_open()) {
      std::cerr << "ERROR output file can't open (no exist)::"<<outFileName<<std::endl;
      exit(1);
      return EXIT_FAILURE;
    }
    for(int id = 0; id<dataSizeOut_; id++){
      ofs<<xd[id]<<" "<< ave[id]<<" "<<err[id]<<std::endl;	  	    
    }     
    ofs.close();
  }
  return 0;
}

template<typename DATA>int IODATA::writeData1D(DATA xd, std::complex<double>* ave,std::complex<double>* err,char* outFileName)
{
  std::fstream ofs;
  if(w_bSwitch_){
    ofs.close();
    std::cout << "error output should use text mode "<<std::endl;
    exit(1);
    return EXIT_FAILURE;
  }
  if(!w_bSwitch_){
    ofs.setf(std::ios::scientific);
    ofs.open(outFileName,std::ios::out);
    if (!ofs.is_open()) {
      std::cerr << "ERROR output file can't open (no exist)::"<<outFileName<<std::endl;
      exit(1);
      return EXIT_FAILURE;
    }
    for(int id = 0; id<dataSizeOut_; id++){
      ofs<<xd[id]<<" "<< ave[id].real()<<" "<<err[id].real()<<std::endl;	  	    
    }
    ofs.close();
  }
  return 0;
}




template<typename DATA> int IODATA::writeData3D(DATA ave,DATA err,char* outFileName)
{
  std::fstream ofs;
  if(w_bSwitch_){
    ofs.close();
    std::cout << "error output should use text mode "<<std::endl;
    exit(1);
    return EXIT_FAILURE;
  }
    
  if(!w_bSwitch_){
    ofs.setf(std::ios::scientific);
    ofs.open(outFileName,std::ios::out);
    if (!ofs.is_open()) {
      std::cerr << "ERROR output file can't open (no exist)::"<<outFileName<<std::endl;
      exit(1);
      return EXIT_FAILURE;
    }
#define radius_sq(x,y,z) ((x)*(x) + (y)*(y) + (z)*(z))
#define min(a,b) (((a) < (b)) ? (a) : (b))
    double tmp = pow(dataSizeOut_,1.0/3.0);
    int Xsize = (int)(2*tmp+dataSizeOut_/tmp/tmp)/3;
    int Ysize = Xsize;
    int Zsize = Xsize;
    int r_sq_max = radius_sq(Xsize,Ysize,Zsize);
    int r_sq = 0;
    if(wAveSwitch_){	
    int* r_sq_count = new int[r_sq_max]();
    double* aveR = new double[r_sq_max]();
    double* errR = new double[r_sq_max]();
    for(int iz = 0; iz<Zsize; iz++){
      for(int iy = 0; iy<Ysize; iy++){
	for(int ix = 0; ix<Xsize; ix++){
	  r_sq = radius_sq( min(ix,Xsize-ix), min(iy,Ysize-iy), min(iz,Zsize-iz) );
	  r_sq_count[r_sq] = r_sq_count[r_sq] +1;
	  aveR[r_sq] =aveR[r_sq] + ave[(ix) +Xsize*((iy) + Ysize*((iz)))];
	  errR[r_sq] =errR[r_sq] + err[(ix) +Xsize*((iy) + Ysize*((iz)))];

	}
      }
    }
    for (int r_sq=0; r_sq<r_sq_max; r_sq++) {
      
      if ( r_sq_count[r_sq] == 0 ) continue;
      
      aveR[r_sq] =aveR[r_sq]/ ((double) r_sq_count[r_sq]);
      errR[r_sq] =errR[r_sq]/ ((double) r_sq_count[r_sq]);
      float rad = sqrt((float)r_sq);
      ofs<<rad<<" "<< aveR[r_sq]<<" "<<errR[r_sq]<<std::endl;
    }
    delete [] r_sq_count;
    delete [] aveR;
    delete [] errR;
      }
    else {
      for(int iz = 0; iz<Zsize; iz++){
	for(int iy = 0; iy<Ysize; iy++){
	  for(int ix = 0; ix<Xsize; ix++){
	    int r_sq = radius_sq( min(ix,Xsize-ix), min(iy,Ysize-iy), min(iz,Zsize-iz) );
	    float rad = (float)r_sq;
	    ofs<<rad<<" "<< ave[(ix) +Xsize*((iy) + Ysize*((iz)))]<<" "<<err[(ix) +Xsize*((iy) + Ysize*((iz)))]<<std::endl;	  	    
	}
      }
    } 
    }
     
    ofs.close();
  }
  return 0;
}

int IODATA::writeData3D(std::complex<double>* ave,std::complex<double>* err,char* outFileName)
{
  std::fstream ofs;
  if(w_bSwitch_){
    ofs.close();
    std::cout << "error output should use text mode "<<std::endl;
    exit(1);
    return EXIT_FAILURE;
  }
  if(!w_bSwitch_){
    ofs.setf(std::ios::scientific);
    ofs.open(outFileName,std::ios::out);
    if (!ofs.is_open()) {
      std::cerr << "ERROR output file can't open (no exist)::"<<outFileName<<std::endl;
      exit(1);
      return EXIT_FAILURE;
    }
#define radius_sq(x,y,z) ((x)*(x) + (y)*(y) + (z)*(z))
#define min(a,b) (((a) < (b)) ? (a) : (b))
    double tmp = pow(dataSizeOut_,1.0/3.0);
    int Xsize = (int)(2*tmp+dataSizeOut_/tmp/tmp)/3;
    int Ysize = Xsize;
    int Zsize = Xsize;
    int r_sq_max = radius_sq(Xsize,Ysize,Zsize);
    int r_sq = 0;
    if(wAveSwitch_)  {
    
    int* r_sq_count = new int[r_sq_max]();
    double* aveR = new double[r_sq_max]();
    double* errR = new double[r_sq_max]();
    for(int iz = 0; iz<Zsize; iz++){
      for(int iy = 0; iy<Ysize; iy++){
	for(int ix = 0; ix<Xsize; ix++){
	  r_sq = radius_sq( min(ix,Xsize-ix), min(iy,Ysize-iy), min(iz,Zsize-iz) );
	  r_sq_count[r_sq] = r_sq_count[r_sq] +1;
	  aveR[r_sq] =aveR[r_sq] + ave[(ix) +Xsize*((iy) + Ysize*((iz)))].real();
	  errR[r_sq] =errR[r_sq] + err[(ix) +Xsize*((iy) + Ysize*((iz)))].real();

	}
      }
    }
    for (int r_sq=0; r_sq<r_sq_max; r_sq++) {
      
      if ( r_sq_count[r_sq] == 0 ) continue;
      
      aveR[r_sq] =aveR[r_sq]/ ((double) r_sq_count[r_sq]);
      errR[r_sq] =errR[r_sq]/ ((double) r_sq_count[r_sq]);
      float rad = sqrt((float)r_sq);
      ofs<<rad<<" "<< aveR[r_sq]<<" "<<errR[r_sq]<<std::endl;
    }
    delete [] r_sq_count;
    delete [] aveR;
    delete [] errR;
    }
    else {
  for(int iz = 0; iz<Zsize; iz++){
    for(int iy = 0; iy<Ysize; iy++){
      for(int ix = 0; ix<Xsize; ix++){
	int r_sq = radius_sq( min(ix,Xsize-ix), min(iy,Ysize-iy), min(iz,Zsize-iz) );
	float rad = (float)r_sq;
	ofs<<rad<<" "<< ave[(ix) +Xsize*((iy) + Ysize*((iz)))].real()<<" "<<err[(ix) +Xsize*((iy) + Ysize*((iz)))].real()<<std::endl;	  	    
      }
    }
  }
    }
    
    ofs.close();
  }
  return 0;
}

#endif

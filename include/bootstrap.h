//--------------------------------------------------------------------------
/**
 * @File bootstrap.h
 * @brief calc bootstrap
 * @ingroup YAMADA
 * @author  M.YAMADA
 * @date    Thu Feb 3  2014
 */
//--------------------------------------------------------------------------
#ifndef BOOTSTRAP_CALC_YAMADA_20140203
#define BOOTSTRAP_CALC_YAMADA_20140203

#include<iostream>
#include<complex>
#include<memory.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

class BOOTSTRAP
{
public:
  BOOTSTRAP();
  ~BOOTSTRAP();

  void set(int DataSize, int ConfSize, int ResampleSize);
  template<typename DATA>
  void setData(DATA data, int iconf);
  void setData(std::complex<double>*  data, int iconf);
  template<typename DATA>
  void setResampleData(DATA data, int iconf);
  void setResampleData(std::complex<double>*  data, int iconf);
  double* calcResample(int b);
  void calcAve(double*);
  void calcErr(double*);

private:
  void resampling();
  void Ave();
  void Err();
  inline double confData_(int id, int j)
  {
    return ConfData_[(id) + DataSize_ *(j)];
  }
  inline double resamleData_(int id, int b)
  {
    return ResampleData_[(id) + DataSize_ *(b)];
  }
  inline int random(int n)
{
    return rand()%n;
}

  double* ConfData_;
  double* ResampleData_;
  double* DoubleResampleData_;
  double* reData_;
  double* ave_;
  double* err_;
  double* doubleAve_;

  int DataSize_;
  int ConfSize_;
  int ResampleSize_;
  bool aveInit_;
};





BOOTSTRAP::BOOTSTRAP()
{

  DataSize_ = 0;
  ConfSize_ = 0;
  ResampleSize_ = 0;
  aveInit_ = false;
  srand(time(NULL));
}
BOOTSTRAP::~BOOTSTRAP()
{
  delete [] ConfData_; ConfData_ = NULL;
  delete [] ResampleData_; ResampleData_ = NULL;
  delete [] DoubleResampleData_; DoubleResampleData_ = NULL;
  delete [] reData_; reData_ = NULL;
  delete [] ave_; ave_ = NULL;
  delete [] err_; err_ = NULL;
  delete [] doubleAve_; doubleAve_ = NULL;
}
void BOOTSTRAP::calcAve(double* outAve)
{
  Ave();
  memcpy(outAve,ave_,sizeof(ave_) * DataSize_);
}
void BOOTSTRAP::calcErr(double* outErr)
{
  if (!aveInit_)  Ave();
  Err();
  memcpy(outErr,err_,sizeof(err_) * DataSize_);
}
void BOOTSTRAP::set(int DataSize, int ConfSize, int ResampleSize)
{
  DataSize_ = DataSize;
  ConfSize_ = ConfSize;
  ResampleSize_ = ResampleSize;
  ConfData_ = new double[DataSize_*ConfSize_]();
  ResampleData_ = new double[DataSize_*ResampleSize_]();
  DoubleResampleData_ = new double[DataSize_*ResampleSize_]();
  reData_ = new double[DataSize_]();
  ave_ = new double[DataSize_]();
  err_ = new double[DataSize_]();
  doubleAve_    = new double[DataSize]();

}
template <typename DATA>     void BOOTSTRAP::setData(DATA in, int iconf){
  double* tmp =new double[DataSize_]();
  for(int id = 0; id <DataSize_; id++){
    tmp[id] = (double)in[id];
  }
  memcpy(ConfData_ + iconf*DataSize_,tmp,sizeof(tmp) * DataSize_);
  delete[] tmp;
  if(ConfSize_ == iconf + 1)  resampling();
}
void BOOTSTRAP::setData(std::complex<double>* in, int iconf)
{
  double* tmp =new double[DataSize_]();
  for(int id = 0; id <DataSize_; id++){
    tmp[id] = (double)in[id].real();
  }
  memcpy(ConfData_ + iconf*DataSize_,tmp,sizeof(tmp) * DataSize_);
  delete[] tmp;
  if(ConfSize_ == iconf + 1)  resampling();
}
template <typename DATA> void BOOTSTRAP::setResampleData(DATA in, int iconf)
{
  double*       tmp = new double[DataSize_]();
  for(int id = 0; id <DataSize_; id++){
    tmp[id] = (double)in[id];
  }
  memcpy(ResampleData_ + iconf*DataSize_,tmp,sizeof(tmp)*DataSize_);
  delete[] tmp;
}
void BOOTSTRAP::setResampleData(std::complex<double>* in, int iconf)
{
  double*       tmp = new double[DataSize_]();
  for(int id = 0; id <DataSize_; id++){
    tmp[id] = (double)in[id].real();
  }
  memcpy(ResampleData_ + iconf*DataSize_,tmp,sizeof(tmp)*DataSize_);
  delete[] tmp;
}

void BOOTSTRAP::resampling()
{
  for (int id = 0 ; id < DataSize_ ; id ++)
    {
      for (int b = 0 ; b < ResampleSize_ ; b++)
	{
	  double tmp = 0.0;
	  for (int j = 0 ; j < ConfSize_ ; j++)
	    {
	      tmp = tmp + confData_(id , random(ConfSize_));
	    }
	  ResampleData_[(id) + DataSize_ *(b)] = tmp / (double)ConfSize_;
	}
    }
}

double* BOOTSTRAP::calcResample(int b)
{
  for (int id =0 ; id <DataSize_ ; id ++)
    {
      reData_[id]  = ResampleData_[(id) + DataSize_ *(b)];
    }
  return reData_;
}

void BOOTSTRAP::Ave()
{
  for (int id = 0 ; id < DataSize_ ; id ++)
    {
      double tmp = 0.0;
      double tmp1 = 0.0;
      for (int b = 0 ; b < ResampleSize_ ; b++)
	{
	  DoubleResampleData_[id + DataSize_*b] =  resamleData_(id, b)*  resamleData_(id, b);
	  tmp = tmp + resamleData_(id, b);
	}
      for (int b = 0 ; b < ResampleSize_ ; b++)
	{
	  tmp1 = tmp1 + DoubleResampleData_[id + DataSize_*b];
	}
      doubleAve_[id] = tmp1/(double)ResampleSize_;
      ave_[id]  = tmp/(double)ResampleSize_;
    }
  aveInit_ = true;
}
void BOOTSTRAP::Err()
{
  for (int id = 0 ; id < DataSize_ ; id ++)
    {
      for (int b = 0 ; b < ResampleSize_ ; b++)
	{
	  err_[id] = sqrt((doubleAve_[id]-(ave_[id]*ave_[id])) * (double)ConfSize_ /((double)ConfSize_ -1.0));
	}
    }
}


#endif

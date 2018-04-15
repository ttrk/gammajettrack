#ifndef SCALEERRORTOOL_H
#define SCALEERRORTOOL_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <fstream>

#include "TMath.h"

//#include "include/doGlobalDebug.h"

class scaleErrorTool{
 public:
  scaleErrorTool();
  scaleErrorTool(std::string inErrorFile);
  ~scaleErrorTool(){};

  bool Init();
  bool Init(std::string inErrorFile);

  std::string getValidString(std::string inFullLine, std::vector<std::string> inCompVect);
  unsigned int getValidPos(std::string compStr, std::vector<std::string> inCompVect);
  std::string getValidStringFromPos(unsigned int pos, std::vector<std::string> inCompVect);
  unsigned int getKey(std::string fullLine);
  unsigned int getKey(std::string centStr, std::string absEtaStr, std::string rStr, std::string flowStr);

  std::string getStringFromKey(unsigned int key);

  double getMuDataMinusMC(std::string centStr, std::string absEtaStr, std::string rStr, std::string flowStr);
  double getMuDataMinusMCErr(std::string centStr, std::string absEtaStr, std::string rStr, std::string flowStr);
  double getSigDataMinusMC(std::string centStr, std::string absEtaStr, std::string rStr, std::string flowStr);
  double getSigDataMinusMCErr(std::string centStr, std::string absEtaStr, std::string rStr, std::string flowStr);

  double getMuDataMinusMC(int centVal, double absEtaVal, int rVal, std::string flowStr);
  double getMuDataMinusMCErr(int centVal, double absEtaVal, int rVal, std::string flowStr);
  double getSigDataMinusMC(int centVal, double absEtaVal, int rVal, std::string flowStr);
  double getSigDataMinusMCErr(int centVal, double absEtaVal, int rVal, std::string flowStr);

  double getMuDataMinusMC(std::string fullLine);
  double getMuDataMinusMCErr(std::string fullLine);
  double getSigDataMinusMC(std::string fullLine);
  double getSigDataMinusMCErr(std::string fullLine);

  std::string errorFile_;

  std::vector<std::string> centStrV = {"0-10%", "10-30%", "30-50%", "50-100%"};//0 thru 3
  std::vector<std::string> absEtaStrV = {"AbsEta0p0to1p0", "AbsEta1p0to2p0"};//0 thru 1
  std::vector<std::string> rStrV = {"R0p3", "R0p4", "R0p6", "R0p8", "R1p0"}; //0 thru 4
  std::vector<std::string> flowStrV = {"NoFlow", "FlowDefaultInRho"};//0 thru 1

  std::vector<int> centValLowV = {0, 10, 30, 50};//0 thru 3
  std::vector<int> centValHiV = {10, 30, 50, 90};//0 thru 3
  std::vector<double> absEtaValLowV = {0., 1.};//0 thru 1
  std::vector<double> absEtaValHiV = {1., 2.};//0 thru 1
  std::vector<int> rValV = {3, 4, 6, 8, 10}; //0 thru 4

  std::map<unsigned int, double> muDataMinusMC;
  std::map<unsigned int, double> sigDataMinusMC;
  std::map<unsigned int, double> muDataMinusMCErr;
  std::map<unsigned int, double> sigDataMinusMCErr;
};

scaleErrorTool::scaleErrorTool()
{
  errorFile_ = "";
  return;
}

scaleErrorTool::scaleErrorTool(std::string inErrorFile)
{
  errorFile_ = inErrorFile;
  return;
}

bool scaleErrorTool::Init()
{
  if(errorFile_.size() == 0  || errorFile_.find(".txt") == std::string::npos){
    std::cout << "SCALEERRORTOOL: Error file string \'" << errorFile_ << "\' is invalid. return false" << std::endl;
    return false;
  }

  std::ifstream file(errorFile_.c_str());
  std::string tempStr;
  while(std::getline(file, tempStr)){
    if(tempStr.find("AbsEta0p0to2p0") != std::string::npos) continue;
    if(tempStr.size() == 0) continue;
    
    unsigned int keyVal = getKey(tempStr);

    bool mcMinusData = true;
    if(tempStr.find("Data,MC") != std::string::npos) mcMinusData = false;

    tempStr.replace(0, tempStr.find("=")+1, "");
    double nominalMu = std::stod(tempStr.substr(0,tempStr.find("#")));
    if(mcMinusData) nominalMu *= -1;
    muDataMinusMC[keyVal] = nominalMu;
    tempStr.replace(0, tempStr.find("#pm")+3, "");
    muDataMinusMCErr[keyVal] = std::stod(tempStr.substr(0,tempStr.find(",")));

    tempStr.replace(0, tempStr.find("=")+1, "");
    double nominalSig = std::stod(tempStr.substr(0,tempStr.find("#")));
    if(mcMinusData) nominalSig *= -1;
    sigDataMinusMC[keyVal] = nominalSig;
    tempStr.replace(0, tempStr.find("#pm")+3, "");
    sigDataMinusMCErr[keyVal] = std::stod(tempStr);
  }

  std::cout << "SCALE TOOL ERROR INIT: " << std::endl;
  for(std::map<unsigned int, double>::iterator mI = muDataMinusMC.begin(); mI != muDataMinusMC.end(); ++mI){
    std::cout << " " << mI->first << ", " << getStringFromKey(mI->first) << ", " << mI->second << std::endl;
  }

  file.close();

  return true;
}


std::string scaleErrorTool::getValidString(std::string inFullLine, std::vector<std::string> inCompVect)
{
  std::string retStr = "";

  for(unsigned int i = 0; i < inCompVect.size(); ++i){
    if(inFullLine.find(inCompVect.at(i)) != std::string::npos){
      retStr = inCompVect.at(i);
      break;
    }
  }

  return retStr;
}

unsigned int scaleErrorTool::getValidPos(std::string compStr, std::vector<std::string> inCompVect)
{
  unsigned int retVal = 9;

  for(unsigned int i = 0; i < inCompVect.size(); ++i){
    if(compStr.find(inCompVect.at(i)) != std::string::npos && compStr.size() == inCompVect.at(i).size()){
      retVal = i;
      break;
    }
  }

  return retVal;
}


std::string scaleErrorTool::getValidStringFromPos(unsigned int pos, std::vector<std::string> inCompVect)
{
  return inCompVect.at(pos);
}


unsigned int scaleErrorTool::getKey(std::string fullLine)
{
  std::string centStr = getValidString(fullLine, centStrV);
  std::string absEtaStr = getValidString(fullLine, absEtaStrV);
  std::string rStr = getValidString(fullLine, rStrV);
  std::string flowStr = getValidString(fullLine, flowStrV);
  
  unsigned int retKey = 99999;
  
  if(centStr.size() == 0 || absEtaStr.size() == 0 || rStr.size() == 0 || flowStr.size() == 0){
    std::cout << "Missing key in line \'" << fullLine << "\'. Init return 9999 key" << std::endl;
    std::cout << " Cent: " << centStr << std::endl;
    std::cout << " AbsEta: " << absEtaStr << std::endl;
    std::cout << " R: " << rStr << std::endl;
    std::cout << " Flow: " << flowStr << std::endl;
  }
  else retKey = getKey(centStr, absEtaStr, rStr, flowStr);

  return retKey;
}

unsigned int scaleErrorTool::getKey(std::string centStr, std::string absEtaStr, std::string rStr, std::string flowStr)
{
  unsigned int retKey = 0;
  retKey += getValidPos(centStr, centStrV);
  retKey += 10*getValidPos(absEtaStr, absEtaStrV);
  retKey += 100*getValidPos(rStr, rStrV);
  retKey += 1000*getValidPos(flowStr, flowStrV);

  return retKey;
}


std::string scaleErrorTool::getStringFromKey(unsigned int key)
{
  std::string keyStr = getValidStringFromPos(key%10, centStrV);
  key /= 10;
  keyStr = keyStr + ", " + getValidStringFromPos(key%10, absEtaStrV);
  key /= 10;
  keyStr = keyStr + ", " + getValidStringFromPos(key%10, rStrV);
  key /= 10;
  keyStr = keyStr + ", " + getValidStringFromPos(key%10, flowStrV);

  return keyStr;
}


double scaleErrorTool::getMuDataMinusMC(std::string centStr, std::string absEtaStr, std::string rStr, std::string flowStr)
{
  unsigned int key = getKey(centStr, absEtaStr, rStr, flowStr);
  return muDataMinusMC[key];
}

double scaleErrorTool::getMuDataMinusMCErr(std::string centStr, std::string absEtaStr, std::string rStr, std::string flowStr)
{
  unsigned int key = getKey(centStr, absEtaStr, rStr, flowStr);
  return muDataMinusMCErr[key];
}

double scaleErrorTool::getSigDataMinusMC(std::string centStr, std::string absEtaStr, std::string rStr, std::string flowStr)
{
  unsigned int key = getKey(centStr, absEtaStr, rStr, flowStr);
  return sigDataMinusMC[key];
}

double scaleErrorTool::getSigDataMinusMCErr(std::string centStr, std::string absEtaStr, std::string rStr, std::string flowStr)
{
  unsigned int key = getKey(centStr, absEtaStr, rStr, flowStr);
  return sigDataMinusMCErr[key];
}

double scaleErrorTool::getMuDataMinusMC(int centVal, double absEtaVal, int rVal, std::string flowStr)
{
  std::string centStr = "";
  std::string absEtaStr = "";
  std::string rStr = "";

  absEtaVal = TMath::Abs(absEtaVal);

  for(unsigned int cI = 0; cI < centValLowV.size(); ++cI){
    if(centValLowV.at(cI) <= centVal && centVal < centValHiV.at(cI)){
      centStr = centStrV.at(cI);
      break;
    }
  }

  for(unsigned int cI = 0; cI < absEtaValLowV.size(); ++cI){
    if(absEtaValLowV.at(cI) <= absEtaVal && absEtaVal < absEtaValHiV.at(cI)){
      absEtaStr = absEtaStrV.at(cI);
      break;
    }
  }

  for(unsigned int cI = 0; cI < rValV.size(); ++cI){
    if(rValV.at(cI) == rVal){
      rStr = rStrV.at(cI);
      break;
    }
  }

  return getMuDataMinusMC(centStr, absEtaStr, rStr, flowStr);
}

double scaleErrorTool::getMuDataMinusMCErr(int centVal, double absEtaVal, int rVal, std::string flowStr)
{
  std::string centStr = "";
  std::string absEtaStr = "";
  std::string rStr = "";

  for(unsigned int cI = 0; cI < centValLowV.size(); ++cI){
    if(centValLowV.at(cI) <= centVal && centVal < centValHiV.at(cI)){
      centStr = centStrV.at(cI);
      break;
    }
  }

  for(unsigned int cI = 0; cI < absEtaValLowV.size(); ++cI){
    if(absEtaValLowV.at(cI) <= absEtaVal && absEtaVal < absEtaValHiV.at(cI)){
      absEtaStr = absEtaStrV.at(cI);
      break;
    }
  }

  for(unsigned int cI = 0; cI < rValV.size(); ++cI){
    if(rValV.at(cI) == rVal){
      rStr = rStrV.at(cI);
      break;
    }
  }

  return getMuDataMinusMCErr(centStr, absEtaStr, rStr, flowStr);
}

double scaleErrorTool::getSigDataMinusMC(int centVal, double absEtaVal, int rVal, std::string flowStr)
{
  std::string centStr = "";
  std::string absEtaStr = "";
  std::string rStr = "";

  for(unsigned int cI = 0; cI < centValLowV.size(); ++cI){
    if(centValLowV.at(cI) <= centVal && centVal < centValHiV.at(cI)){
      centStr = centStrV.at(cI);
      break;
    }
  }

  for(unsigned int cI = 0; cI < absEtaValLowV.size(); ++cI){
    if(absEtaValLowV.at(cI) <= absEtaVal && absEtaVal < absEtaValHiV.at(cI)){
      absEtaStr = absEtaStrV.at(cI);
      break;
    }
  }

  for(unsigned int cI = 0; cI < rValV.size(); ++cI){
    if(rValV.at(cI) == rVal){
      rStr = rStrV.at(cI);
      break;
    }
  }

  return getSigDataMinusMC(centStr, absEtaStr, rStr, flowStr);
}

double scaleErrorTool::getSigDataMinusMCErr(int centVal, double absEtaVal, int rVal, std::string flowStr)
{
  std::string centStr = "";
  std::string absEtaStr = "";
  std::string rStr = "";

  for(unsigned int cI = 0; cI < centValLowV.size(); ++cI){
    if(centValLowV.at(cI) <= centVal && centVal < centValHiV.at(cI)){
      centStr = centStrV.at(cI);
      break;
    }
  }

  for(unsigned int cI = 0; cI < absEtaValLowV.size(); ++cI){
    if(absEtaValLowV.at(cI) <= absEtaVal && absEtaVal < absEtaValHiV.at(cI)){
      absEtaStr = absEtaStrV.at(cI);
      break;
    }
  }

  for(unsigned int cI = 0; cI < rValV.size(); ++cI){
    if(rValV.at(cI) == rVal){
      rStr = rStrV.at(cI);
      break;
    }
  }

  return getSigDataMinusMCErr(centStr, absEtaStr, rStr, flowStr);
}


double scaleErrorTool::getMuDataMinusMC(std::string fullLine)
{
  std::string centStr = getValidString(fullLine, centStrV);
  std::string absEtaStr = getValidString(fullLine, absEtaStrV);
  std::string rStr = getValidString(fullLine, rStrV);
  std::string flowStr = getValidString(fullLine, flowStrV);

  unsigned int key = getKey(centStr, absEtaStr, rStr, flowStr);
  return muDataMinusMC[key];
}

double scaleErrorTool::getMuDataMinusMCErr(std::string fullLine)
{
  std::string centStr = getValidString(fullLine, centStrV);
  std::string absEtaStr = getValidString(fullLine, absEtaStrV);
  std::string rStr = getValidString(fullLine, rStrV);
  std::string flowStr = getValidString(fullLine, flowStrV);

  unsigned int key = getKey(centStr, absEtaStr, rStr, flowStr);
  return muDataMinusMCErr[key];
}

double scaleErrorTool::getSigDataMinusMC(std::string fullLine)
{
  std::string centStr = getValidString(fullLine, centStrV);
  std::string absEtaStr = getValidString(fullLine, absEtaStrV);
  std::string rStr = getValidString(fullLine, rStrV);
  std::string flowStr = getValidString(fullLine, flowStrV);

  unsigned int key = getKey(centStr, absEtaStr, rStr, flowStr);
  return sigDataMinusMC[key];
}

double scaleErrorTool::getSigDataMinusMCErr(std::string fullLine)
{
  std::string centStr = getValidString(fullLine, centStrV);
  std::string absEtaStr = getValidString(fullLine, absEtaStrV);
  std::string rStr = getValidString(fullLine, rStrV);
  std::string flowStr = getValidString(fullLine, flowStrV);

  unsigned int key = getKey(centStr, absEtaStr, rStr, flowStr);
  return sigDataMinusMCErr[key];
}


#endif

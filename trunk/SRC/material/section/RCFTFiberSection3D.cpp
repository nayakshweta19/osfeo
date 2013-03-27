/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Cenk Tort (tort0008@umn.edu)				      **
**   Jerome F. Hajjar (hajjar@struc.ce.umn.edu)                       **
**							              **
**   University of Minnesota - Twin Cities		   	      **
**								      **
** ****************************************************************** */

// $Revision: 1.3 $
// $Date: 2008-09-10 18:32:14 $
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/material/section/RCFTFiberSection3D.cpp,v $

// Written: cenk tort
// Created: 09/04
//
// Description: This file contains the class implementation of RCFTFiberSection3D.cpp

#include <stdlib.h>

#include <Channel.h>
#include <Vector.h>
#include <Matrix.h>
#include <MatrixUtil.h>
#include <Fiber.h>
#include <classTags.h>
#include <RCFTFiberSection3D.h>
#include <ID.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <MaterialResponse.h>
#include <UniaxialMaterial.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace::std;


ID RCFTFiberSection3D::code(12);

//using std::endl;
//using std::ios;
//using std::ifstream;
//using std::ofstream;

// constructors:
RCFTFiberSection3D::RCFTFiberSection3D(int tag, int num, Fiber **fibers, double gj, double d, double b, double t):
  SectionForceDeformation(tag, SEC_TAG_RCFTFiberSection3D), numFibers(num), theMaterials(0), matData(0),
  yBar(0.0), zBar(0.0), e(6), eCommit(6), s(0), ks(0), ksCommit(12,12), GJ(gj), D(d), B(b), T(t)
{
#ifdef COMPOSITE_DEBUG
    ofstream output;
    output.open("fibers.dat",ios::app);
#endif

    double EAc = 0.0;
    double EQzc = 0.0;
    double EQyc = 0.0;
    double EIzc = 0.0;
    double EIyc = 0.0;
    double EIyzc = 0.0;
    double EAs = 0.0;
    double EQzs = 0.0;
    double EQys = 0.0;
    double EIzs = 0.0;
    double EIys = 0.0;
    double EIyzs = 0.0;
	if (numFibers != 0){
		theMaterials = new UniaxialMaterial *[numFibers];

		if (theMaterials == 0) {
			opserr << "RCFTFiberSection3D::RCFTFiberSection3D -- failed to allocate Material pointers\n";
			exit(-1);
		}

		matData = new double [numFibers*3];

		if (matData == 0) {
			opserr << "RCFTFiberSection3D::RCFTFiberSection3D -- failed to allocate double array for material data\n";
			exit(-1);
		}

		double Qz = 0.0;
		double Qy = 0.0;
		double A  = 0.0;

		double Qzc = 0.0;
		double Qyc = 0.0;
		double Ac  = 0.0;

		double Qzs = 0.0;
		double Qys = 0.0;
		double As  = 0.0;

		for (int i = 0; i < numFibers; i++) {
			Fiber *theFiber = fibers[i];
			double Es, Ec, yLoc, zLoc, Area;
			theFiber->getFiberLocation(yLoc, zLoc);
			Area = theFiber->getArea();
			UniaxialMaterial *theMat = theFiber->getMaterial();
			if(theMat->getInitialTangent() < 25000.0){
				Ec = theMat->getInitialTangent();  
#ifdef COMPOSITE_DEBUG
				output<<i<<"  "<<yLoc<<"   "<<zLoc<<"   "<<Area<<"    "<<Ec<<endl;
#endif
				Qzc += yLoc*Area;
				Qyc += zLoc*Area;
				Ac  += Area;
				EQzc += yLoc*Area*Ec;
				EQyc += zLoc*Area*Ec;
				EAc  += Area*Ec;
				EIzc += yLoc*yLoc*Area*Ec;
				EIyc += zLoc*zLoc*Area*Ec;
				EIyzc += yLoc*zLoc*Area*Ec;
			}
			else if(theMat->getInitialTangent() > 25000.0){
				Es = theMat->getInitialTangent(); 
#ifdef COMPOSITE_DEBUG
				output<<i<<"  "<<yLoc<<"   "<<zLoc<<"   "<<Area<<"   "<<Es<<endl;
#endif
				Qzs += yLoc*Area;
				Qys += zLoc*Area;
				As  += Area;
				EQzs += yLoc*Area*Es;
				EQys += zLoc*Area*Es;
				EAs  += Area*Es;
				EIzs += yLoc*yLoc*Area*Es;
				EIys += zLoc*zLoc*Area*Es;
				EIyzs += yLoc*zLoc*Area*Es;
			}

			Qz += yLoc*Area;
			Qy += zLoc*Area;
			A  += Area;

			matData[i*3] = yLoc;
			matData[i*3+1] = zLoc;
			matData[i*3+2] = Area;

			theMaterials[i] = theMat->getCopy();

			if (theMaterials[i] == 0) {
				opserr << "RCFTFiberSection3D::RCFTFiberSection3D -- failed to get copy of a Material\n";
				exit(-1);
			}
		}

		yBar = Qz/A;
		zBar = Qy/A;

	}

    s = new Vector(sData,12);
    ks = new Matrix(kData,12,12);

    sData[0] = 0.0;
    sData[1] = 0.0;
    sData[2] = 0.0;
    sData[3] = 0.0;
    sData[4] = 0.0;
    sData[5] = 0.0;
    sData[6] = 0.0;
    sData[7] = 0.0;
    sData[8] = 0.0;
    sData[9] = 0.0;
    sData[10]= 0.0;
    sData[11]= 0.0;

    for (int i=0; i<144; i++)
      kData[i] = 0.0;

    kData[0] = EAc; 
    kData[13] = EIzc; 
    kData[26] = EIyc;
    kData[39] = EAs; 
    kData[52] = EIzs; 
    kData[65] = EIys;
   
    //kData[0] = 224801.03;
    //kData[13] = 3484423.0;
    //kData[26] = 3484423.0;
    //kData[39] = 224801.03;
    //kData[52] = 3484423.0;
    //kData[65] = 3484423.0;
    
    //kData[0] = 500.0;
    //kData[13] = 0.5;
    //kData[26] = 0.5;
    //kData[39] = 500.0;
    //kData[52] = 0.5;
    //kData[65] = 0.5;

    //kData[0]  = 500000.0;
    //kData[13] = 500.0;
    //kData[26] = 500.0;
    //kData[39] = 500000.0;
    //kData[52] = 500.0;
    //kData[65] = 500.0;
    
    //kData[0] = 145000000;
    //kData[13] =1450000;
    //kData[26] =1450000;
    //kData[39] =145000000;
    //kData[52] =1450000;
    //kData[65] =1450000;
    
    //kData[0] = 150000;
    //kData[13] = 1500000;
    //kData[26] = 1500000;
    //kData[39] = 150000;
    //kData[52] = 1500000;
    //kData[65] = 1500000;
    
    //kData[0] = 21600000;
    //kData[13] = 7200000;
    //kData[26] = 7200000;
    //kData[39] = 21600000;
    //kData[52] = 7200000;
    //kData[65] = 7200000;
  
    kData[78] = GJ; kData[91] = D; kData[104] = B; 
    kData[117] = T; kData[130] = 1; kData[143] = 1;

    ksCommit = *ks;

    code(0) = SECTION_RESPONSE_Pc;
    code(1) = SECTION_RESPONSE_MZc;
    code(2) = SECTION_RESPONSE_MYc;
    code(3) = SECTION_RESPONSE_LYYc;
    code(4) = SECTION_RESPONSE_LZZc;
    code(5) = SECTION_RESPONSE_LYZc;
    code(6) = SECTION_RESPONSE_Ps;
    code(7) = SECTION_RESPONSE_MZs;
    code(8) = SECTION_RESPONSE_MYs;
    code(9) = SECTION_RESPONSE_LYYs;
    code(10) = SECTION_RESPONSE_LZZs;
    code(11) = SECTION_RESPONSE_LYZs;

}

// constructor for blank object that recvSelf needs to be invoked upon
RCFTFiberSection3D::RCFTFiberSection3D():
  SectionForceDeformation(0, SEC_TAG_RCFTFiberSection3D),
  numFibers(0), theMaterials(0), matData(0),
  yBar(0.0), zBar(0.0), e(6), eCommit(6), s(0), ks(0), ksCommit(12,12), GJ(1.0), D(1.0), B(1.0), T(1.0)
{

  s = new Vector(sData, 12);
  ks = new Matrix(kData, 12, 12);

  sData[0] = 0.0;
  sData[1] = 0.0;
  sData[2] = 0.0;
  sData[3] = 0.0;
  sData[4] = 0.0;
  sData[5] = 0.0;
  sData[6] = 0.0;
  sData[7] = 0.0;
  sData[8] = 0.0;
  sData[9] = 0.0;
  sData[10] = 0.0;
  sData[11] = 0.0;

  for (int i=0; i<144; i++)
    kData[i] = 0.0;

  kData[78] = 1.0;
  kData[91] = 1.0;
  kData[104] = 1.0;
  kData[117] = 1.0;
  kData[130] = 1.0;
  kData[143] = 1.0;

  ksCommit = (*ks);

  code(0) = SECTION_RESPONSE_Pc;
  code(1) = SECTION_RESPONSE_MZc;
  code(2) = SECTION_RESPONSE_MYc;
  code(3) = SECTION_RESPONSE_LYYc;
  code(4) = SECTION_RESPONSE_LZZc;
  code(5) = SECTION_RESPONSE_LYZc;
  code(6) = SECTION_RESPONSE_Ps;
  code(7) = SECTION_RESPONSE_MZs;
  code(8) = SECTION_RESPONSE_MYs;
  code(9) = SECTION_RESPONSE_LYYs;
  code(10) = SECTION_RESPONSE_LZZs;
  code(11) = SECTION_RESPONSE_LYZs;

}

int
RCFTFiberSection3D::addFiber(Fiber &newFiber)
{
  // need to create a larger array
  int newSize = numFibers+1;

  UniaxialMaterial **newArray = new UniaxialMaterial *[newSize];

  double *newMatData = new double [3 * newSize];

  if (newArray == 0 || newMatData == 0) {
   opserr << "RCFTFiberSection3D::addFiber -- failed to allocate Fiber pointers\n";
   exit(-1);
 }
  // copy the old pointers
  int i;
  for (i = 0; i < numFibers; i++) {
    newArray[i] = theMaterials[i];
    newMatData[3*i] = matData[3*i];
    newMatData[3*i+1] = matData[3*i+1];
    newMatData[3*i+2] = matData[3*i+2];
  }


  // set the new pointers
  double yLoc, zLoc, Area;
  newFiber.getFiberLocation(yLoc, zLoc);
  Area = newFiber.getArea();
  newMatData[numFibers*3] = -yLoc;
  newMatData[numFibers*3+1] = zLoc;
  newMatData[numFibers*3+2] = Area;
  UniaxialMaterial *theMat = newFiber.getMaterial();
  newArray[numFibers] = theMat->getCopy();

  if (newArray[numFibers] == 0) {
    opserr << "RCFTFiberSection3D::addFiber -- failed to get copy of a Material\n";
    exit(-1);
    delete [] newArray;
    delete [] newMatData;
    return -1;
  }

  numFibers++;

  if (theMaterials != 0) {
    delete [] theMaterials;
    delete [] matData;
  }

  theMaterials = newArray;
  matData = newMatData;

  double Qz = 0.0;
  double Qy = 0.0;
  double A  = 0.0;

  double Qzc = 0.0;
  double Qyc = 0.0;
  double Ac  = 0.0;

  double Qzs = 0.0;
  double Qys = 0.0;
  double As  = 0.0;

  // Recompute centroid
  for (i = 0; i < numFibers; i++) {
    yLoc = matData[3*i];
    zLoc = matData[3*i+1];
    Area = matData[3*i+2];
    if( theMat->getInitialTangent() < 25000.0 )
    {
	Ac  += Area;
	Qzc += yLoc*Area;
	Qyc += zLoc*Area;
    }
    else if( theMat->getInitialTangent() < 25000.0 )
    {
	As  += Area;
	Qzs += yLoc*Area;
	Qys += zLoc*Area;
    }
    A  += Area;
    Qz += yLoc*Area;
    Qy += zLoc*Area;
  }

  yBar = Qz/A;
  zBar = Qy/A;

  return 0;
}



// destructor:
RCFTFiberSection3D::~RCFTFiberSection3D()
{
  if (theMaterials != 0) {
    for (int i = 0; i < numFibers; i++)
      if (theMaterials[i] != 0)
	delete theMaterials[i];

    delete [] theMaterials;
  }

  if (matData != 0)
    delete [] matData;

  if (s != 0)
    delete s;

  if (ks != 0)
    delete ks;

}

int
RCFTFiberSection3D::setTrialSectionDeformation (const Vector &deforms)
{

  int res = 0;
  
  e = deforms;

  kData[0] = 0.0;  kData[1] = 0.0;  kData[2] = 0.0;  kData[3] = 0.0;
  kData[4] = 0.0;  kData[5] = 0.0;  kData[6] = 0.0;  kData[7] = 0.0;
  kData[8] = 0.0;  kData[9] = 0.0;  kData[10] = 0.0; kData[11] = 0.0;
  kData[12] = 0.0; kData[13] = 0.0; kData[14] = 0.0; kData[15] = 0.0;
  kData[16] = 0.0; kData[17] = 0.0; kData[18] = 0.0; kData[19] = 0.0;
  kData[20] = 0.0; kData[21] = 0.0; kData[22] = 0.0; kData[23] = 0.0;
  kData[24] = 0.0; kData[25] = 0.0; kData[26] = 0.0; kData[27] = 0.0;
  kData[28] = 0.0; kData[29] = 0.0; kData[30] = 0.0; kData[31] = 0.0;
  kData[32] = 0.0; kData[33] = 0.0; kData[34] = 0.0; kData[35] = 0.0;
  kData[36] = 0.0; kData[37] = 0.0; kData[38] = 0.0; kData[39] = 0.0;
  kData[40] = 0.0; kData[41] = 0.0; kData[42] = 0.0; kData[43] = 0.0;
  kData[44] = 0.0; kData[45] = 0.0; kData[46] = 0.0; kData[47] = 0.0; 
  kData[48] = 0.0; kData[49] = 0.0; kData[50] = 0.0; kData[51] = 0.0;
  kData[52] = 0.0; kData[53] = 0.0; kData[54] = 0.0; kData[55] = 0.0;
  kData[56] = 0.0; kData[57] = 0.0; kData[58] = 0.0; kData[59] = 0.0;
  kData[60] = 0.0; kData[61] = 0.0; kData[62] = 0.0; kData[63] = 0.0;
  kData[64] = 0.0; kData[65] = 0.0; kData[66] = 0.0; kData[67] = 0.0;
  kData[68] = 0.0; kData[69] = 0.0; kData[70] = 0.0; kData[71] = 0.0;
  kData[72] = 0.0; kData[73] = 0.0; kData[74] = 0.0; kData[75] = 0.0;
  kData[76] = 0.0; kData[77] = 0.0; kData[78] = 0.0; kData[79] = 0.0;
  kData[80] = 0.0; kData[81] = 0.0; kData[82] = 0.0; kData[83] = 0.0;
  kData[84] = 0.0; kData[85] = 0.0; kData[86] = 0.0; kData[87] = 0.0;
  kData[88] = 0.0; kData[89] = 0.0; kData[90] = 0.0; kData[91] = 0.0;
  kData[92] = 0.0; kData[93] = 0.0; kData[94] = 0.0; kData[95] = 0.0;
  kData[96] = 0.0; kData[97] = 0.0; kData[98] = 0.0; kData[99] = 0.0;
  kData[100] = 0.0; kData[101] = 0.0; kData[102] = 0.0; kData[103] = 0.0;
  kData[104] = 0.0; kData[105] = 0.0; kData[106] = 0.0; kData[107] = 0.0;
  kData[108] = 0.0; kData[109] = 0.0; kData[110] = 0.0; kData[111] = 0.0;
  kData[112] = 0.0; kData[113] = 0.0; kData[114] = 0.0; kData[115] = 0.0;
  kData[116] = 0.0; kData[117] = 0.0; kData[118] = 0.0; kData[119] = 0.0;
  kData[120] = 0.0; kData[121] = 0.0; kData[122] = 0.0; kData[123] = 0.0;
  kData[124] = 0.0; kData[125] = 0.0; kData[126] = 0.0; kData[127] = 0.0;
  kData[128] = 0.0; kData[129] = 0.0; kData[130] = 0.0; kData[131] = 0.0;
  kData[132] = 0.0; kData[133] = 0.0; kData[134] = 0.0; kData[135] = 0.0;
  kData[136] = 0.0; kData[137] = 0.0; kData[138] = 0.0; kData[139] = 0.0;
  kData[140] = 0.0; kData[141] = 0.0; kData[142] = 0.0; kData[143] = 0.0;

  sData[0] = 0.0; sData[1] = 0.0;  sData[2] = 0.0;
  sData[3] = 0.0; sData[4] = 0.0;  sData[5] = 0.0;
  sData[6] = 0.0; sData[7] = 0.0;  sData[8] = 0.0;
  sData[9] = 0.0; sData[10] = 0.0; sData[11] = 0.0;

  int loc = 0;

  double d0 = deforms(0);
  double d1 = deforms(1);
  double d2 = deforms(2);
  double d3 = deforms(3);
  double d4 = deforms(4);
  double d5 = deforms(5);

  for (int i = 0; i < numFibers; i++) {
    UniaxialMaterial *theMat = theMaterials[i];
    double y = matData[loc++] - yBar; 
    double z = matData[loc++] - zBar; 
    double A = matData[loc++];

    // determine material strain and set it
    double strain; 
    double tangent, stress;
    double value;
    double vas1;
    double vas2;
    double vas1as1;
    double vas1as2;
    double vas2as2;
 
    if( theMat->getInitialTangent() < 25000.0 ){
    
	strain = d0 - z * d2 - y * d1;
	res = theMat->setTrial(strain, stress, tangent);
        value = tangent * A;
	vas1 = y*value;
	vas2 = z*value;
	vas1as1 = vas1*y;
	vas1as2 = vas1*z;
	vas2as2 = vas2*z;

	kData[0] += value;
	kData[1] -= vas1;
	kData[2] -= vas2;

	kData[12] -= vas1;
	kData[13] += vas1as1;
	kData[14] += vas1as2;

	kData[24] -= vas2;
	kData[25] += vas1as2;
	kData[26] += vas2as2;

	double fs0 = stress * A;
	sData[0] += fs0;
	sData[1] -= fs0 * y;
	sData[2] -= fs0 * z;
	sData[3] += fs0 * z * z;
	sData[4] += fs0 * y * y;
	sData[5] += fs0 * y * z;
    }
    if( theMat->getInitialTangent() > 25000.0 ){

        strain = d3 - z * d5 - y * d4;
	res = theMat->setTrial(strain, stress, tangent);
        value = tangent * A;
	vas1 = y*value;
	vas2 = z*value;
	vas1as1 = vas1*y;
	vas1as2 = vas1*z;
	vas2as2 = vas2*z;
	
	kData[39] += value;
	kData[40] -= vas1;
	kData[41] -= vas2;

	kData[51] -= vas1;
	kData[52] += vas1as1;
	kData[53] += vas1as2;

	kData[63] -= vas2;
	kData[64] += vas1as2;
	kData[65] += vas2as2;

	double fs0 = stress * A;
	sData[6] += fs0;
	sData[7] -= fs0 * y;
	sData[8] -= fs0 * z;
	sData[9] += fs0 * z * z;
	sData[10] += fs0 * y * y;
	sData[11] += fs0 * y * z;
    }
  }

  //check7<<sData[0]<<"  "<<sData[1]<<"  "<<sData[2]<<"  "<<sData[6]<<"   "<<sData[7]<<"   "<<sData[8]<<endl;
  //check7<<kData[0]<<"   "<<kData[1]<<"   "<<kData[2]<<"   "<<kData[3]<<"   "<<kData[4]<<"   "<<kData[5]<<endl; 
  //kData[39]<<"   "<<kData[40]<<"   "<<kData[41]<<endl;
  //kData[1]  = 0.0;
  //kData[2]  = 0.0;

  //kData[12] = 0.0;
  //kData[14] = 0.0;

  //kData[24] = 0.0;
  //kData[25] = 0.0;

  //kData[40] = 0.0;
  //kData[41] = 0.0;

  //kData[51] = 0.0;
  //kData[53] = 0.0;

  //kData[63] = 0.0;
  //kData[64] = 0.0;

  //kData[0] = 224801.0;
  //kData[13] = 3484423.0;
  //kData[26] = 3484423.0;

  //kData[39] = 224801.0;
  //kData[52] = 3484423.0;
  //kData[65] = 3484423.0;

  //kData[0] = 224801.0;
  //kData[13] = 3484423.0;
  //kData[26] = 3484423.0;

  //kData[39] = 224801.0;
  //kData[52] = 3484423.0;
  //kData[65] = 3484423.0;

  //kData[0] = 500.0;
  //kData[13] = 0.5;
  //kData[26] = 0.5;

  //kData[39] = 500.0;
  //kData[52] = 0.5;
  //kData[65] = 0.5;
  
  //kData[0] = 150000.0;
  //kData[13] = 1500000.0;
  //kData[26] = 1500000.0;
 
  //kData[39] = 150000.0;
  //kData[52] = 1500000.0;
  //kData[65] = 1500000.0;
  
  //kData[0] =  21600000.0;
  //kData[13] = 7200000.0;
  //kData[26] = 7200000.0;
  
  //kData[39] = 21600000.0;
  //kData[52] = 7200000.0;
  //kData[65] = 7200000.0;

  kData[3] = 0.0;
  kData[4] = 0.0;  kData[5] = 0.0;  kData[6] = 0.0;  kData[7] = 0.0;
  kData[8] = 0.0;  kData[9] = 0.0;  kData[10] = 0.0; kData[11] = 0.0;
  kData[15] = 0.0;
  kData[16] = 0.0; kData[17] = 0.0; kData[18] = 0.0; kData[19] = 0.0;
  kData[20] = 0.0; kData[21] = 0.0; kData[22] = 0.0; kData[23] = 0.0;
  kData[27] = 0.0;
  kData[28] = 0.0; kData[29] = 0.0; kData[30] = 0.0; kData[31] = 0.0;
  kData[32] = 0.0; kData[33] = 0.0; kData[34] = 0.0; kData[35] = 0.0;
  kData[36] = 0.0; kData[37] = 0.0; kData[38] = 0.0; 
  kData[42] = 0.0; kData[43] = 0.0;
  kData[44] = 0.0; kData[45] = 0.0; kData[46] = 0.0; kData[47] = 0.0;
  kData[48] = 0.0; kData[49] = 0.0; kData[50] = 0.0; 
  kData[54] = 0.0; kData[55] = 0.0;
  kData[56] = 0.0; kData[57] = 0.0; kData[58] = 0.0; kData[59] = 0.0;
  kData[60] = 0.0; kData[61] = 0.0; kData[62] = 0.0; 
  kData[66] = 0.0; kData[67] = 0.0;
  kData[68] = 0.0; kData[69] = 0.0; kData[70] = 0.0; kData[71] = 0.0;
  kData[72] = 0.0; kData[73] = 0.0; kData[74] = 0.0; kData[75] = 0.0;
  kData[76] = 0.0; kData[77] = 0.0; 
  kData[78] = GJ; 
  kData[79] = 0.0;
  kData[80] = 0.0; kData[81] = 0.0; kData[82] = 0.0; kData[83] = 0.0;
  kData[84] = 0.0; kData[85] = 0.0; kData[86] = 0.0; kData[87] = 0.0;
  kData[88] = 0.0; kData[89] = 0.0; kData[90] = 0.0; 
  kData[91] = D;
  kData[92] = 0.0; kData[93] = 0.0; kData[94] = 0.0; kData[95] = 0.0;
  kData[96] = 0.0; kData[97] = 0.0; kData[98] = 0.0; kData[99] = 0.0;
  kData[100] = 0.0; kData[101] = 0.0; kData[102] = 0.0; kData[103] = 0.0;
  kData[104] = B; 
  kData[105] = 0.0; kData[106] = 0.0; kData[107] = 0.0;
  kData[108] = 0.0; kData[109] = 0.0; kData[110] = 0.0; kData[111] = 0.0;
  kData[112] = 0.0; kData[113] = 0.0; kData[114] = 0.0; kData[115] = 0.0;
  kData[116] = 0.0; 
  kData[117] = T; 
  kData[118] = 0.0; kData[119] = 0.0;
  kData[120] = 0.0; kData[121] = 0.0; kData[122] = 0.0; kData[123] = 0.0;
  kData[124] = 0.0; kData[125] = 0.0; kData[126] = 0.0; kData[127] = 0.0;
  kData[128] = 0.0; kData[129] = 0.0; 
  kData[130] = 1.0; 
  kData[131] = 0.0;
  kData[132] = 0.0; kData[133] = 0.0; kData[134] = 0.0; kData[135] = 0.0;
  kData[136] = 0.0; kData[137] = 0.0; kData[138] = 0.0; kData[139] = 0.0;
  kData[140] = 0.0; kData[141] = 0.0; kData[142] = 0.0; 
  kData[143] = 1.0;

  for( int i = 0; i < 144; i++ ){
    if( fabs(kData[i]) < 1.0e-10 ){
        kData[i] = 0.0;   
    }
  }

  return res;
}

const Matrix&
RCFTFiberSection3D::getInitialTangent(void)
{
  kData[0] = 0.0;  kData[1] = 0.0;  kData[2] = 0.0;  kData[3] = 0.0;
  kData[4] = 0.0;  kData[5] = 0.0;  kData[6] = 0.0;  kData[7] = 0.0;
  kData[8] = 0.0;  kData[9] = 0.0;  kData[10] = 0.0; kData[11] = 0.0;
  kData[12] = 0.0; kData[13] = 0.0; kData[14] = 0.0; kData[15] = 0.0;
  kData[16] = 0.0; kData[17] = 0.0; kData[18] = 0.0; kData[19] = 0.0;
  kData[20] = 0.0; kData[21] = 0.0; kData[22] = 0.0; kData[23] = 0.0;
  kData[24] = 0.0; kData[25] = 0.0; kData[26] = 0.0; kData[27] = 0.0;
  kData[28] = 0.0; kData[29] = 0.0; kData[30] = 0.0; kData[31] = 0.0;
  kData[32] = 0.0; kData[33] = 0.0; kData[34] = 0.0; kData[35] = 0.0;
  kData[36] = 0.0; kData[37] = 0.0; kData[38] = 0.0; kData[39] = 0.0;
  kData[40] = 0.0; kData[41] = 0.0; kData[42] = 0.0; kData[43] = 0.0;
  kData[44] = 0.0; kData[45] = 0.0; kData[46] = 0.0; kData[47] = 0.0;
  kData[48] = 0.0; kData[49] = 0.0; kData[50] = 0.0; kData[51] = 0.0;
  kData[52] = 0.0; kData[53] = 0.0; kData[54] = 0.0; kData[55] = 0.0;
  kData[56] = 0.0; kData[57] = 0.0; kData[58] = 0.0; kData[59] = 0.0;
  kData[60] = 0.0; kData[61] = 0.0; kData[62] = 0.0; kData[63] = 0.0;
  kData[64] = 0.0; kData[65] = 0.0; kData[66] = 0.0; kData[67] = 0.0;
  kData[68] = 0.0; kData[69] = 0.0; kData[70] = 0.0; kData[71] = 0.0;
  kData[72] = 0.0; kData[73] = 0.0; kData[74] = 0.0; kData[75] = 0.0;
  kData[76] = 0.0; kData[77] = 0.0; kData[78] = 0.0; kData[79] = 0.0;
  kData[80] = 0.0; kData[81] = 0.0; kData[82] = 0.0; kData[83] = 0.0;
  kData[84] = 0.0; kData[85] = 0.0; kData[86] = 0.0; kData[87] = 0.0;
  kData[88] = 0.0; kData[89] = 0.0; kData[90] = 0.0; kData[91] = 0.0;
  kData[92] = 0.0; kData[93] = 0.0; kData[94] = 0.0; kData[95] = 0.0;
  kData[96] = 0.0; kData[97] = 0.0; kData[98] = 0.0; kData[99] = 0.0;
  kData[100] = 0.0; kData[101] = 0.0; kData[102] = 0.0; kData[103] = 0.0;
  kData[104] = 0.0; kData[105] = 0.0; kData[106] = 0.0; kData[107] = 0.0;
  kData[108] = 0.0; kData[109] = 0.0; kData[110] = 0.0; kData[111] = 0.0;
  kData[112] = 0.0; kData[113] = 0.0; kData[114] = 0.0; kData[115] = 0.0;
  kData[116] = 0.0; kData[117] = 0.0; kData[118] = 0.0; kData[119] = 0.0;
  kData[120] = 0.0; kData[121] = 0.0; kData[122] = 0.0; kData[123] = 0.0;
  kData[124] = 0.0; kData[125] = 0.0; kData[126] = 0.0; kData[127] = 0.0;
  kData[128] = 0.0; kData[129] = 0.0; kData[130] = 0.0; kData[131] = 0.0;
  kData[132] = 0.0; kData[133] = 0.0; kData[134] = 0.0; kData[135] = 0.0;
  kData[136] = 0.0; kData[137] = 0.0; kData[138] = 0.0; kData[139] = 0.0;
  kData[140] = 0.0; kData[141] = 0.0; kData[142] = 0.0; kData[143] = 0.0;
									  
  for (int i = 0; i < numFibers; i++) {
    UniaxialMaterial *theMat = theMaterials[i];
    double y = matData[3*i] - yBar;
    double z = matData[3*i+1] - zBar;
    double A = matData[3*i+2];

    double tangent = theMat->getInitialTangent();

    double value = tangent * A;
    double vas1 = y*value;
    double vas2 = z*value;
    double vas1as1 = vas1*y;
    double vas1as2 = vas1*z;
    double vas2as2 = vas2*z;

    if( theMat->getInitialTangent() < 25000.0 ){
	kData[0] += value;
	kData[1] -= vas1;
	kData[2] -= vas2;

	kData[12] -= vas1;
	kData[13] += vas1as1;
	kData[14] += vas1as2;

	kData[24] -= vas2;
	kData[25] += vas1as2;
	kData[26] += vas2as2;
    }
    if( theMat->getInitialTangent() > 25000.0 ){
        kData[39] += value;
	kData[40] -= vas1;
	kData[41] -= vas2;

	kData[51] -= vas1;
	kData[52] += vas1as1;
	kData[53] += vas1as2;

	kData[63] -= vas2;
	kData[64] += vas1as2;
	kData[65] += vas2as2;
    }
  }

  kData[1]  = 0.0;
  kData[2]  = 0.0;

  kData[12] = 0.0;
  kData[14] = 0.0;

  kData[24] = 0.0;
  kData[25] = 0.0;

  kData[40] = 0.0;
  kData[41] = 0.0;

  kData[51] = 0.0;
  kData[53] = 0.0;

  kData[63] = 0.0;
  kData[64] = 0.0;
	
  kData[3] = 0.0;
  kData[4] = 0.0;  kData[5] = 0.0;  kData[6] = 0.0;  kData[7] = 0.0;
  kData[8] = 0.0;  kData[9] = 0.0;  kData[10] = 0.0; kData[11] = 0.0;
  kData[15] = 0.0;
  kData[16] = 0.0; kData[17] = 0.0; kData[18] = 0.0; kData[19] = 0.0;
  kData[20] = 0.0; kData[21] = 0.0; kData[22] = 0.0; kData[23] = 0.0;
  kData[27] = 0.0;
  kData[28] = 0.0; kData[29] = 0.0; kData[30] = 0.0; kData[31] = 0.0;
  kData[32] = 0.0; kData[33] = 0.0; kData[34] = 0.0; kData[35] = 0.0;
  kData[36] = 0.0; kData[37] = 0.0; kData[38] = 0.0;
  kData[42] = 0.0; kData[43] = 0.0;
  kData[44] = 0.0; kData[45] = 0.0; kData[46] = 0.0; kData[47] = 0.0;
  kData[48] = 0.0; kData[49] = 0.0; kData[50] = 0.0;
  kData[54] = 0.0; kData[55] = 0.0;
  kData[56] = 0.0; kData[57] = 0.0; kData[58] = 0.0; kData[59] = 0.0;
  kData[60] = 0.0; kData[61] = 0.0; kData[62] = 0.0;
  kData[66] = 0.0; kData[67] = 0.0;
  kData[68] = 0.0; kData[69] = 0.0; kData[70] = 0.0; kData[71] = 0.0;
  kData[72] = 0.0; kData[73] = 0.0; kData[74] = 0.0; kData[75] = 0.0;
  kData[76] = 0.0; kData[77] = 0.0;
  kData[78] = GJ;
  kData[79] = 0.0;
  kData[80] = 0.0; kData[81] = 0.0; kData[82] = 0.0; kData[83] = 0.0;
  kData[84] = 0.0; kData[85] = 0.0; kData[86] = 0.0; kData[87] = 0.0;
  kData[88] = 0.0; kData[89] = 0.0; kData[90] = 0.0;
  kData[91] = D;
  kData[92] = 0.0; kData[93] = 0.0; kData[94] = 0.0; kData[95] = 0.0;
  kData[96] = 0.0; kData[97] = 0.0; kData[98] = 0.0; kData[99] = 0.0;
  kData[100] = 0.0; kData[101] = 0.0; kData[102] = 0.0; kData[103] = 0.0;
  kData[104] = B;
  kData[105] = 0.0; kData[106] = 0.0; kData[107] = 0.0;
  kData[108] = 0.0; kData[109] = 0.0; kData[110] = 0.0; kData[111] = 0.0;
  kData[112] = 0.0; kData[113] = 0.0; kData[114] = 0.0; kData[115] = 0.0;
  kData[116] = 0.0;
  kData[117] = T;
  kData[118] = 0.0; kData[119] = 0.0;
  kData[120] = 0.0; kData[121] = 0.0; kData[122] = 0.0; kData[123] = 0.0;
  kData[124] = 0.0; kData[125] = 0.0; kData[126] = 0.0; kData[127] = 0.0;
  kData[128] = 0.0; kData[129] = 0.0;
  kData[130] = 1.0;
  kData[131] = 0.0;
  kData[132] = 0.0; kData[133] = 0.0; kData[134] = 0.0; kData[135] = 0.0;
  kData[136] = 0.0; kData[137] = 0.0; kData[138] = 0.0; kData[139] = 0.0;
  kData[140] = 0.0; kData[141] = 0.0; kData[142] = 0.0;
  kData[143] = 1.0;  

  (*ks)(0,0) = kData[0]; (*ks)(1,0) = kData[1]; (*ks)(2,0) = kData[2]; (*ks)(3,0) = kData[3];
  (*ks)(4,0) = kData[4]; (*ks)(5,0) = kData[5]; (*ks)(6,0) = kData[6]; (*ks)(7,0) = kData[7];
  (*ks)(8,0) = kData[8]; (*ks)(9,0) = kData[9]; (*ks)(10,0) = kData[10]; (*ks)(11,0) = kData[11];

  (*ks)(0,1) = kData[12]; (*ks)(1,1) = kData[13]; (*ks)(2,1) = kData[14]; (*ks)(3,1) = kData[15];
  (*ks)(4,1) = kData[16]; (*ks)(5,1) = kData[17]; (*ks)(6,1) = kData[18]; (*ks)(7,1) = kData[19];
  (*ks)(8,1) = kData[20]; (*ks)(9,1) = kData[21]; (*ks)(10,1) = kData[22]; (*ks)(11,1) = kData[23];

  (*ks)(0,2) = kData[24]; (*ks)(1,2) = kData[25]; (*ks)(2,2) = kData[26]; (*ks)(3,2) = kData[27];
  (*ks)(4,2) = kData[28]; (*ks)(5,2) = kData[29]; (*ks)(6,2) = kData[30]; (*ks)(7,2) = kData[31];
  (*ks)(8,2) = kData[32]; (*ks)(9,2) = kData[33]; (*ks)(10,2) = kData[34]; (*ks)(11,2) = kData[35];

  (*ks)(0,3) = kData[36]; (*ks)(1,3) = kData[37]; (*ks)(2,3) = kData[38]; (*ks)(3,3) = kData[39];
  (*ks)(4,3) = kData[40]; (*ks)(5,3) = kData[41]; (*ks)(6,3) = kData[42]; (*ks)(7,3) = kData[43];
  (*ks)(8,3) = kData[44]; (*ks)(9,3) = kData[45]; (*ks)(10,3) = kData[46]; (*ks)(11,3) = kData[47];

  (*ks)(0,4) = kData[48]; (*ks)(1,4) = kData[49]; (*ks)(2,4) = kData[50]; (*ks)(3,4) = kData[51];
  (*ks)(4,4) = kData[52]; (*ks)(5,4) = kData[53]; (*ks)(6,4) = kData[54]; (*ks)(7,4) = kData[55];
  (*ks)(8,4) = kData[56]; (*ks)(9,4) = kData[57]; (*ks)(10,4) = kData[58]; (*ks)(11,4) = kData[59];

  (*ks)(0,5) = kData[60]; (*ks)(1,5) = kData[61]; (*ks)(2,5) = kData[62]; (*ks)(3,5) = kData[63];
  (*ks)(4,5) = kData[64]; (*ks)(5,5) = kData[65]; (*ks)(6,5) = kData[66]; (*ks)(7,5) = kData[67];
  (*ks)(8,5) = kData[68]; (*ks)(9,5) = kData[69]; (*ks)(10,5) = kData[70]; (*ks)(11,5) = kData[71];

  (*ks)(0,6) = kData[72]; (*ks)(1,6) = kData[73]; (*ks)(2,6) = kData[74]; (*ks)(3,6) = kData[75];
  (*ks)(4,6) = kData[76]; (*ks)(5,6) = kData[77]; (*ks)(6,6) = kData[78]; (*ks)(7,6) = kData[79];
  (*ks)(8,6) = kData[80]; (*ks)(9,6) = kData[81]; (*ks)(10,6) = kData[82]; (*ks)(11,6) = kData[83];

  (*ks)(0,7) = kData[84]; (*ks)(1,7) = kData[85]; (*ks)(2,7) = kData[86]; (*ks)(3,7) = kData[87];
  (*ks)(4,7) = kData[88]; (*ks)(5,7) = kData[89]; (*ks)(6,7) = kData[90]; (*ks)(7,7) = kData[91];
  (*ks)(8,7) = kData[92]; (*ks)(9,7) = kData[93]; (*ks)(10,7) = kData[94]; (*ks)(11,7) = kData[95];
						
  (*ks)(0,8) = kData[96]; (*ks)(1,8) = kData[97]; (*ks)(2,8) = kData[98]; (*ks)(3,8) = kData[99];
  (*ks)(4,8) = kData[100]; (*ks)(5,8) = kData[101]; (*ks)(6,8) = kData[102]; (*ks)(7,8) = kData[103];
  (*ks)(8,8) = kData[104]; (*ks)(9,8) = kData[105]; (*ks)(10,8) = kData[106]; (*ks)(11,8) = kData[107];

  (*ks)(0,9) = kData[108]; (*ks)(1,9) = kData[109]; (*ks)(2,9) = kData[110]; (*ks)(3,9) = kData[111];
  (*ks)(4,9) = kData[112]; (*ks)(5,9) = kData[113]; (*ks)(6,9) = kData[114]; (*ks)(7,9) = kData[115];
  (*ks)(8,9) = kData[116]; (*ks)(9,9) = kData[117]; (*ks)(10,9) = kData[118]; (*ks)(11,9) = kData[119];

  (*ks)(0,10) = kData[120]; (*ks)(1,10) = kData[121]; (*ks)(2,10) = kData[122]; (*ks)(3,10) = kData[123];
  (*ks)(4,10) = kData[124]; (*ks)(5,10) = kData[125]; (*ks)(6,10) = kData[126]; (*ks)(7,10) = kData[127];
  (*ks)(8,10) = kData[128]; (*ks)(9,10) = kData[129]; (*ks)(10,10) = kData[130]; (*ks)(11,10) = kData[131];

  (*ks)(0,11) = kData[132]; (*ks)(1,11) = kData[133]; (*ks)(2,11) = kData[134]; (*ks)(3,11) = kData[135];
  (*ks)(4,11) = kData[136]; (*ks)(5,11) = kData[137]; (*ks)(6,11) = kData[138]; (*ks)(7,11) = kData[139];
  (*ks)(8,11) = kData[140]; (*ks)(9,11) = kData[141]; (*ks)(10,11) = kData[142]; (*ks)(11,11) = kData[143];
  
  return *ks;
}

const Vector&
RCFTFiberSection3D::getSectionDeformation(void)
{
  return e;
}

const Matrix& 
RCFTFiberSection3D::getSectionCommitTangent(void){
  return ksCommit;
}

const Matrix&
RCFTFiberSection3D::getSectionTangent(void)
{
  (*ks)(0,0) = kData[0]; (*ks)(1,0) = kData[1]; (*ks)(2,0) = kData[2]; (*ks)(3,0) = kData[3]; 
  (*ks)(4,0) = kData[4]; (*ks)(5,0) = kData[5]; (*ks)(6,0) = kData[6]; (*ks)(7,0) = kData[7];
  (*ks)(8,0) = kData[8]; (*ks)(9,0) = kData[9]; (*ks)(10,0) = kData[10]; (*ks)(11,0) = kData[11];

  (*ks)(0,1) = kData[12]; (*ks)(1,1) = kData[13]; (*ks)(2,1) = kData[14]; (*ks)(3,1) = kData[15];
  (*ks)(4,1) = kData[16]; (*ks)(5,1) = kData[17]; (*ks)(6,1) = kData[18]; (*ks)(7,1) = kData[19];
  (*ks)(8,1) = kData[20]; (*ks)(9,1) = kData[21]; (*ks)(10,1) = kData[22]; (*ks)(11,1) = kData[23];

  (*ks)(0,2) = kData[24]; (*ks)(1,2) = kData[25]; (*ks)(2,2) = kData[26]; (*ks)(3,2) = kData[27];
  (*ks)(4,2) = kData[28]; (*ks)(5,2) = kData[29]; (*ks)(6,2) = kData[30]; (*ks)(7,2) = kData[31];
  (*ks)(8,2) = kData[32]; (*ks)(9,2) = kData[33]; (*ks)(10,2) = kData[34]; (*ks)(11,2) = kData[35];

  (*ks)(0,3) = kData[36]; (*ks)(1,3) = kData[37]; (*ks)(2,3) = kData[38]; (*ks)(3,3) = kData[39];
  (*ks)(4,3) = kData[40]; (*ks)(5,3) = kData[41]; (*ks)(6,3) = kData[42]; (*ks)(7,3) = kData[43];
  (*ks)(8,3) = kData[44]; (*ks)(9,3) = kData[45]; (*ks)(10,3) = kData[46]; (*ks)(11,3) = kData[47];

  (*ks)(0,4) = kData[48]; (*ks)(1,4) = kData[49]; (*ks)(2,4) = kData[50]; (*ks)(3,4) = kData[51];
  (*ks)(4,4) = kData[52]; (*ks)(5,4) = kData[53]; (*ks)(6,4) = kData[54]; (*ks)(7,4) = kData[55];
  (*ks)(8,4) = kData[56]; (*ks)(9,4) = kData[57]; (*ks)(10,4) = kData[58]; (*ks)(11,4) = kData[59];

  (*ks)(0,5) = kData[60]; (*ks)(1,5) = kData[61]; (*ks)(2,5) = kData[62]; (*ks)(3,5) = kData[63];
  (*ks)(4,5) = kData[64]; (*ks)(5,5) = kData[65]; (*ks)(6,5) = kData[66]; (*ks)(7,5) = kData[67];
  (*ks)(8,5) = kData[68]; (*ks)(9,5) = kData[69]; (*ks)(10,5) = kData[70]; (*ks)(11,5) = kData[71];

  (*ks)(0,6) = kData[72]; (*ks)(1,6) = kData[73]; (*ks)(2,6) = kData[74]; (*ks)(3,6) = kData[75];
  (*ks)(4,6) = kData[76]; (*ks)(5,6) = kData[77]; (*ks)(6,6) = kData[78]; (*ks)(7,6) = kData[79];
  (*ks)(8,6) = kData[80]; (*ks)(9,6) = kData[81]; (*ks)(10,6) = kData[82]; (*ks)(11,6) = kData[83];

  (*ks)(0,7) = kData[84]; (*ks)(1,7) = kData[85]; (*ks)(2,7) = kData[86]; (*ks)(3,7) = kData[87];
  (*ks)(4,7) = kData[88]; (*ks)(5,7) = kData[89]; (*ks)(6,7) = kData[90]; (*ks)(7,7) = kData[91];
  (*ks)(8,7) = kData[92]; (*ks)(9,7) = kData[93]; (*ks)(10,7) = kData[94]; (*ks)(11,7) = kData[95];

  (*ks)(0,8) = kData[96]; (*ks)(1,8) = kData[97]; (*ks)(2,8) = kData[98]; (*ks)(3,8) = kData[99];
  (*ks)(4,8) = kData[100]; (*ks)(5,8) = kData[101]; (*ks)(6,8) = kData[102]; (*ks)(7,8) = kData[103];
  (*ks)(8,8) = kData[104]; (*ks)(9,8) = kData[105]; (*ks)(10,8) = kData[106]; (*ks)(11,8) = kData[107];

  (*ks)(0,9) = kData[108]; (*ks)(1,9) = kData[109]; (*ks)(2,9) = kData[110]; (*ks)(3,9) = kData[111];
  (*ks)(4,9) = kData[112]; (*ks)(5,9) = kData[113]; (*ks)(6,9) = kData[114]; (*ks)(7,9) = kData[115];
  (*ks)(8,9) = kData[116]; (*ks)(9,9) = kData[117]; (*ks)(10,9) = kData[118]; (*ks)(11,9) = kData[119];

  (*ks)(0,10) = kData[120]; (*ks)(1,10) = kData[121]; (*ks)(2,10) = kData[122]; (*ks)(3,10) = kData[123];
  (*ks)(4,10) = kData[124]; (*ks)(5,10) = kData[125]; (*ks)(6,10) = kData[126]; (*ks)(7,10) = kData[127];
  (*ks)(8,10) = kData[128]; (*ks)(9,10) = kData[129]; (*ks)(10,10) = kData[130]; (*ks)(11,10) = kData[131];

  (*ks)(0,11) = kData[132]; (*ks)(1,11) = kData[133]; (*ks)(2,11) = kData[134]; (*ks)(3,11) = kData[135];
  (*ks)(4,11) = kData[136]; (*ks)(5,11) = kData[137]; (*ks)(6,11) = kData[138]; (*ks)(7,11) = kData[139];
  (*ks)(8,11) = kData[140]; (*ks)(9,11) = kData[141]; (*ks)(10,11) = kData[142]; (*ks)(11,11) = kData[143];
   
  return *ks;
}

const Vector&
RCFTFiberSection3D::getStressResultant(void)
{
  return *s;
}

SectionForceDeformation*
RCFTFiberSection3D::getCopy(void)
{

  RCFTFiberSection3D *theCopy = new RCFTFiberSection3D();
  theCopy->setTag(this->getTag());

  theCopy->numFibers = numFibers;


  if (numFibers != 0) {
    theCopy->theMaterials = new UniaxialMaterial *[numFibers];

    if (theCopy->theMaterials == 0) {
    opserr << "RCFTFiberSection3D::RCFTFiberSection3D -- failed to allocate Material pointers\n";
    exit(-1);
    }
    
    theCopy->matData = new double [numFibers*3];

    if (theCopy->matData == 0) {
    opserr << "RCFTFiberSection3D::RCFTFiberSection3D -- failed to allocate double array for material data\n";
    exit(-1);
    }
       
    for (int i = 0; i < numFibers; i++) {
      theCopy->matData[i*3] = matData[i*3];
      theCopy->matData[i*3+1] = matData[i*3+1];
      theCopy->matData[i*3+2] = matData[i*3+2];
      theCopy->theMaterials[i] = theMaterials[i]->getCopy();

      if (theCopy->theMaterials[i] == 0) {
      opserr << "RCFTFiberSection3D::getCopy -- failed to get copy of a Material\n";
      exit(-1);
      }
    }   	
  }

  theCopy->eCommit = eCommit;
  theCopy->e = e;
  theCopy->yBar = yBar;
  theCopy->zBar = zBar;
  theCopy->ksCommit = ksCommit;

 for (int i=0; i<144; i++)
    theCopy->kData[i] = kData[i];

  theCopy->sData[0] = sData[0];
  theCopy->sData[1] = sData[1];
  theCopy->sData[2] = sData[2];
  theCopy->sData[3] = sData[3];
  theCopy->sData[4] = sData[4];
  theCopy->sData[5] = sData[5];
  theCopy->sData[6] = sData[6];
  theCopy->sData[7] = sData[7];
  theCopy->sData[8] = sData[8];
  theCopy->sData[9] = sData[9];
  theCopy->sData[10] = sData[10];
  theCopy->sData[11] = sData[11];

  theCopy->GJ = GJ;
  theCopy->D = D;
  theCopy->B = B;
  theCopy->T = T;

  return theCopy;
}

const ID&
RCFTFiberSection3D::getType ()
{
  return code;
}

int
RCFTFiberSection3D::getOrder () const
{
  return 12;
}

int
RCFTFiberSection3D::commitState(void)
{
  int err = 0;

  for (int i = 0; i < numFibers; i++){
    err += theMaterials[i]->commitState();
  }

#ifdef COMPOSITE_DEBUG
  ofstream stl36;
  stl36.open("stl36.dat",ios::app);

  ofstream stl40;
  stl40.open("stl40.dat",ios::app);

  ofstream stl100;
  stl100.open("stl100.dat",ios::app);

  ofstream stl108;
  stl108.open("stl108.dat",ios::app);

  ofstream conc0;
  conc0.open("conc0.dat",ios::app);

  ofstream conc2;
  conc2.open("conc2.dat",ios::app);

  for( int i = 0; i < 42; i++ ){
    stringstream number;
    number << i;
    string conc = "strsStrn";
    string type = ".dat";
    string name = conc + number.str() + type;
    ofstream output;
    output.open(name.c_str(),ios::app);
    output<<theMaterials[i]->getStrain()<<" "<<theMaterials[i]->getStress()<<" "<<theMaterials[i]->getTangent()<<endl;
  }

  for( int i = 100; i < 136; i++ ){
    stringstream number;
    number << i;
    string conc = "stl";
    string type = ".dat";
    string name = conc + number.str() + type;
    ofstream output;
    output.open(name.c_str(),ios::app);
    output<<theMaterials[i]->getStrain()<<" "<<theMaterials[i]->getStress()<<" "<<theMaterials[i]->getTangent()<<endl;
  }
 
  stl100<<theMaterials[100]->getStrain()<<" "<<theMaterials[100]->getStress()<<" "<<theMaterials[100]->getTangent()<<endl;
  stl108<<theMaterials[108]->getStrain()<<" "<<theMaterials[108]->getStress()<<" "<<theMaterials[108]->getTangent()<<endl;

  conc0<<theMaterials[0]->getStrain()<<" "<<theMaterials[0]->getStress()<<" "<<theMaterials[0]->getTangent()<<endl;
  conc2<<theMaterials[2]->getStrain()<<" "<<theMaterials[2]->getStress()<<" "<<theMaterials[2]->getTangent()<<endl;
#endif

  eCommit = e;

  ksCommit = *ks;

  return err;
}

int
RCFTFiberSection3D::revertToLastCommit(void)
{
  int err = 0;

  int res = 0;

  // Last committed section deformations
  e = eCommit;

  kData[0] = 0.0;  kData[1] = 0.0;  kData[2] = 0.0;  kData[3] = 0.0;
  kData[4] = 0.0;  kData[5] = 0.0;  kData[6] = 0.0;  kData[7] = 0.0;
  kData[8] = 0.0;  kData[9] = 0.0;  kData[10] = 0.0; kData[11] = 0.0;
  kData[12] = 0.0; kData[13] = 0.0; kData[14] = 0.0; kData[15] = 0.0;
  kData[16] = 0.0; kData[17] = 0.0; kData[18] = 0.0; kData[19] = 0.0;
  kData[20] = 0.0; kData[21] = 0.0; kData[22] = 0.0; kData[23] = 0.0;
  kData[24] = 0.0; kData[25] = 0.0; kData[26] = 0.0; kData[27] = 0.0;
  kData[28] = 0.0; kData[29] = 0.0; kData[30] = 0.0; kData[31] = 0.0;
  kData[32] = 0.0; kData[33] = 0.0; kData[34] = 0.0; kData[35] = 0.0;
  kData[36] = 0.0; kData[37] = 0.0; kData[38] = 0.0; kData[39] = 0.0;
  kData[40] = 0.0; kData[41] = 0.0; kData[42] = 0.0; kData[43] = 0.0;
  kData[44] = 0.0; kData[45] = 0.0; kData[46] = 0.0; kData[47] = 0.0; 
  kData[48] = 0.0; kData[49] = 0.0; kData[50] = 0.0; kData[51] = 0.0;
  kData[52] = 0.0; kData[53] = 0.0; kData[54] = 0.0; kData[55] = 0.0;
  kData[56] = 0.0; kData[57] = 0.0; kData[58] = 0.0; kData[59] = 0.0;
  kData[60] = 0.0; kData[61] = 0.0; kData[62] = 0.0; kData[63] = 0.0;
  kData[64] = 0.0; kData[65] = 0.0; kData[66] = 0.0; kData[67] = 0.0;
  kData[68] = 0.0; kData[69] = 0.0; kData[70] = 0.0; kData[71] = 0.0;
  kData[72] = 0.0; kData[73] = 0.0; kData[74] = 0.0; kData[75] = 0.0;
  kData[76] = 0.0; kData[77] = 0.0; kData[78] = 0.0; kData[79] = 0.0;
  kData[80] = 0.0; kData[81] = 0.0; kData[82] = 0.0; kData[83] = 0.0;
  kData[84] = 0.0; kData[85] = 0.0; kData[86] = 0.0; kData[87] = 0.0;
  kData[88] = 0.0; kData[89] = 0.0; kData[90] = 0.0; kData[91] = 0.0;
  kData[92] = 0.0; kData[93] = 0.0; kData[94] = 0.0; kData[95] = 0.0;
  kData[96] = 0.0; kData[97] = 0.0; kData[98] = 0.0; kData[99] = 0.0;
  kData[100] = 0.0; kData[101] = 0.0; kData[102] = 0.0; kData[103] = 0.0;
  kData[104] = 0.0; kData[105] = 0.0; kData[106] = 0.0; kData[107] = 0.0;
  kData[108] = 0.0; kData[109] = 0.0; kData[110] = 0.0; kData[111] = 0.0;
  kData[112] = 0.0; kData[113] = 0.0; kData[114] = 0.0; kData[115] = 0.0;
  kData[116] = 0.0; kData[117] = 0.0; kData[118] = 0.0; kData[119] = 0.0;
  kData[120] = 0.0; kData[121] = 0.0; kData[122] = 0.0; kData[123] = 0.0;
  kData[124] = 0.0; kData[125] = 0.0; kData[126] = 0.0; kData[127] = 0.0;
  kData[128] = 0.0; kData[129] = 0.0; kData[130] = 0.0; kData[131] = 0.0;
  kData[132] = 0.0; kData[133] = 0.0; kData[134] = 0.0; kData[135] = 0.0;
  kData[136] = 0.0; kData[137] = 0.0; kData[138] = 0.0; kData[139] = 0.0;
  kData[140] = 0.0; kData[141] = 0.0; kData[142] = 0.0; kData[143] = 0.0;

  int loc = 0;

  double d0 = e(0);
  double d1 = e(1);
  double d2 = e(2);
  double d3 = e(3);
  double d4 = e(4);
  double d5 = e(5);

  for (int i = 0; i < numFibers; i++) {
    UniaxialMaterial *theMat = theMaterials[i];
    double y = matData[loc++] - yBar; 
    double z = matData[loc++] - zBar; 
    double A = matData[loc++];

    err += theMat->revertToLastCommit();

    // determine material strain and set it
    double strain; 
    double tangent, stress;
    double value;
    double vas1;
    double vas2;
    double vas1as1;
    double vas1as2;
    double vas2as2;
 
    if( theMat->getInitialTangent() < 25000.0 ){
    
	//strain = d0 - z * d2 - y * d1;
        strain  = theMat->getStrain();
        tangent = theMat->getTangent();
        stress = theMat->getStress();
        value = tangent * A;
	vas1 = y*value;
	vas2 = z*value;
	vas1as1 = vas1*y;
	vas1as2 = vas1*z;
	vas2as2 = vas2*z;

	kData[0] += value;
	kData[1] -= vas1;
	kData[2] -= vas2;

	kData[12] -= vas1;
	kData[13] += vas1as1;
	kData[14] += vas1as2;

	kData[24] -= vas2;
	kData[25] += vas1as2;
	kData[26] += vas2as2;

	double fs0 = stress * A;
	sData[0] += fs0;
	sData[1] -= fs0 * y;
	sData[2] -= fs0 * z;
	sData[3] += fs0 * z * z;
	sData[4] += fs0 * y * y;
	sData[5] += fs0 * y * z;
    }
    if( theMat->getInitialTangent() > 25000.0 ){

        //strain = d3 - z * d5 - y * d4;
        strain = theMat->getStrain();
        tangent = theMat->getTangent();
        stress = theMat->getStress();
        value = tangent * A;
	vas1 = y*value;
	vas2 = z*value;
	vas1as1 = vas1*y;
	vas1as2 = vas1*z;
	vas2as2 = vas2*z;
	
	kData[39] += value;
	kData[40] -= vas1;
	kData[41] -= vas2;

	kData[51] -= vas1;
	kData[52] += vas1as1;
	kData[53] += vas1as2;

	kData[63] -= vas2;
	kData[64] += vas1as2;
	kData[65] += vas2as2;

	double fs0 = stress * A;
	sData[6] += fs0;
	sData[7] -= fs0 * y;
	sData[8] -= fs0 * z;
	sData[9] += fs0 * z * z;
	sData[10] += fs0 * y * y;
	sData[11] += fs0 * y * z;
    }
  }

  kData[3] = 0.0;
  kData[4] = 0.0;  kData[5] = 0.0;  kData[6] = 0.0;  kData[7] = 0.0;
  kData[8] = 0.0;  kData[9] = 0.0;  kData[10] = 0.0; kData[11] = 0.0;
  kData[15] = 0.0;
  kData[16] = 0.0; kData[17] = 0.0; kData[18] = 0.0; kData[19] = 0.0;
  kData[20] = 0.0; kData[21] = 0.0; kData[22] = 0.0; kData[23] = 0.0;
  kData[27] = 0.0;
  kData[28] = 0.0; kData[29] = 0.0; kData[30] = 0.0; kData[31] = 0.0;
  kData[32] = 0.0; kData[33] = 0.0; kData[34] = 0.0; kData[35] = 0.0;
  kData[36] = 0.0; kData[37] = 0.0; kData[38] = 0.0; 
  kData[42] = 0.0; kData[43] = 0.0;
  kData[44] = 0.0; kData[45] = 0.0; kData[46] = 0.0; kData[47] = 0.0;
  kData[48] = 0.0; kData[49] = 0.0; kData[50] = 0.0; 
  kData[54] = 0.0; kData[55] = 0.0;
  kData[56] = 0.0; kData[57] = 0.0; kData[58] = 0.0; kData[59] = 0.0;
  kData[60] = 0.0; kData[61] = 0.0; kData[62] = 0.0; 
  kData[66] = 0.0; kData[67] = 0.0;
  kData[68] = 0.0; kData[69] = 0.0; kData[70] = 0.0; kData[71] = 0.0;
  kData[72] = 0.0; kData[73] = 0.0; kData[74] = 0.0; kData[75] = 0.0;
  kData[76] = 0.0; kData[77] = 0.0; 
  kData[78] = GJ; 
  kData[79] = 0.0;
  kData[80] = 0.0; kData[81] = 0.0; kData[82] = 0.0; kData[83] = 0.0;
  kData[84] = 0.0; kData[85] = 0.0; kData[86] = 0.0; kData[87] = 0.0;
  kData[88] = 0.0; kData[89] = 0.0; kData[90] = 0.0; 
  kData[91] = D;
  kData[92] = 0.0; kData[93] = 0.0; kData[94] = 0.0; kData[95] = 0.0;
  kData[96] = 0.0; kData[97] = 0.0; kData[98] = 0.0; kData[99] = 0.0;
  kData[100] = 0.0; kData[101] = 0.0; kData[102] = 0.0; kData[103] = 0.0;
  kData[104] = B; 
  kData[105] = 0.0; kData[106] = 0.0; kData[107] = 0.0;
  kData[108] = 0.0; kData[109] = 0.0; kData[110] = 0.0; kData[111] = 0.0;
  kData[112] = 0.0; kData[113] = 0.0; kData[114] = 0.0; kData[115] = 0.0;
  kData[116] = 0.0; 
  kData[117] = T; 
  kData[118] = 0.0; kData[119] = 0.0;
  kData[120] = 0.0; kData[121] = 0.0; kData[122] = 0.0; kData[123] = 0.0;
  kData[124] = 0.0; kData[125] = 0.0; kData[126] = 0.0; kData[127] = 0.0;
  kData[128] = 0.0; kData[129] = 0.0; 
  kData[130] = 1.0; 
  kData[131] = 0.0;
  kData[132] = 0.0; kData[133] = 0.0; kData[134] = 0.0; kData[135] = 0.0;
  kData[136] = 0.0; kData[137] = 0.0; kData[138] = 0.0; kData[139] = 0.0;
  kData[140] = 0.0; kData[141] = 0.0; kData[142] = 0.0; 
  kData[143] = 1.0;

  for( int i = 0; i < 144; i++ ){
    if( fabs(kData[i]) < 1.0e-10 ){
        kData[i] = 0.0;   
    }
  }

  //kData[0] = 21600000;
  //kData[13] = 7200000;
  //kData[26] = 7200000;
  //kData[39] = 21600000;
  //kData[52] = 7200000;
  //kData[65] = 7200000;

  //kData[0] = 500.0;
  //kData[13] = 0.5;
  //kData[26] = 0.5;
  //kData[39] = 500.0;
  //kData[52] = 0.5;
  //kData[65] = 0.5;

  //kData[0] = 1875;
  //kData[13] = 2.445;
  //kData[26] = 2.445;
  //kData[39] = 1875;
  //kData[52] = 2.445;
  //kData[65] = 2.445;

  //kData[0] = 500000.0;
  //kData[13] = 500.0;
  //kData[26] = 500.0;
  //kData[39] = 500000.0;
  //kData[52] = 500.0;
  //kData[65] = 500.0;

  //VV>>*ks;

  return err;
}

int
RCFTFiberSection3D::revertToStart(void)
{
  // revert the fibers to start    
  int err = 0;

  kData[0] = 0.0;  kData[1] = 0.0;  kData[2] = 0.0;  kData[3] = 0.0;
  kData[4] = 0.0;  kData[5] = 0.0;  kData[6] = 0.0;  kData[7] = 0.0;
  kData[8] = 0.0;  kData[9] = 0.0;  kData[10] = 0.0; kData[11] = 0.0;
  kData[12] = 0.0; kData[13] = 0.0; kData[14] = 0.0; kData[15] = 0.0;
  kData[16] = 0.0; kData[17] = 0.0; kData[18] = 0.0; kData[19] = 0.0;
  kData[20] = 0.0; kData[21] = 0.0; kData[22] = 0.0; kData[23] = 0.0;
  kData[24] = 0.0; kData[25] = 0.0; kData[26] = 0.0; kData[27] = 0.0;
  kData[28] = 0.0; kData[29] = 0.0; kData[30] = 0.0; kData[31] = 0.0;
  kData[32] = 0.0; kData[33] = 0.0; kData[34] = 0.0; kData[35] = 0.0;
  kData[36] = 0.0; kData[37] = 0.0; kData[38] = 0.0; kData[39] = 0.0;
  kData[40] = 0.0; kData[41] = 0.0; kData[42] = 0.0; kData[43] = 0.0;
  kData[44] = 0.0; kData[45] = 0.0; kData[46] = 0.0; kData[47] = 0.0; kData[48] = 0.0;
	

  sData[0] = 0.0; sData[1] = 0.0;  sData[2] = 0.0; 
  sData[3] = 0.0; sData[4] = 0.0;  sData[5] = 0.0;
  sData[6] = 0.0; sData[7] = 0.0;  sData[8] = 0.0;
  sData[9] = 0.0; sData[10] = 0.0; sData[11] = 0.0; 

  int loc = 0;

  for (int i = 0; i < numFibers; i++) {
    UniaxialMaterial *theMat = theMaterials[i];
    double y = matData[loc++] - yBar;
    double z = matData[loc++] - zBar;
    double A = matData[loc++];

    // invoke revertToStart on the material
    err += theMat->revertToStart();

    double tangent = theMat->getTangent();
    double stress = theMat->getStress();

    double value = tangent * A;
    double vas1 = y*value;
    double vas2 = z*value;
    double vas1as1 = vas1*y;
    double vas1as2 = vas1*z;
    double vas2as2 = vas2*z;

    double fs0 = stress * A;

    if( theMat->getInitialTangent() < 25000.0 ){

	kData[0] += value;
	kData[1] -= vas1;
	kData[2] -= vas2;
    
	kData[7] -= vas1;
	kData[8] += vas1as1;
	kData[9] += vas1as2;

	kData[14] -= vas2;
	kData[15] += vas1as2;
	kData[16] += vas2as2;

	sData[0] += fs0;
	sData[1] -= fs0 * z;
	sData[2] -= fs0 * y;
        sData[3] += fs0 * z * z;
	sData[4] += fs0 * y * y;
	sData[5] += fs0 * y * z;
						
    }
    if( theMat->getInitialTangent() > 25000.0 ){
	kData[24] += value;
	kData[25] -= vas1;
	kData[26] -= vas2;
    
	kData[31] -= vas1;
	kData[32] += vas1as1;
	kData[33] += vas1as2;

	kData[38] -= vas2;
	kData[39] += vas1as2;
	kData[40] += vas2as2;

	sData[6] += fs0;
	sData[7] -= fs0 * z;
	sData[8] -= fs0 * y;
	sData[9] += fs0 * z * z;
	sData[10] += fs0 * y * y;
	sData[11] += fs0 * y * z;
    }
  }
  kData[3] = 0.0;
  kData[4] = 0.0;
  kData[5] = 0.0;
  kData[6] = 0.0;

  kData[10] = 0.0;
  kData[11] = 0.0;
  kData[12] = 0.0;
  kData[13] = 0.0;

  kData[17] = 0.0;
  kData[18] = 0.0;
  kData[19] = 0.0;
  kData[20] = 0.0;

  kData[21] = 0.0;
  kData[22] = 0.0;
  kData[23] = 0.0;
  kData[27] = 0.0;

  kData[28] = 0.0;
  kData[29] = 0.0;
  kData[30] = 0.0;

  kData[34] = 0.0;
  kData[35] = 0.0;
  kData[36] = 0.0;
  kData[37] = 0.0;
  kData[41] = 0.0;
  kData[42] = 0.0;
  kData[43] = 0.0;
  kData[44] = 0.0;
  kData[45] = 0.0;
  kData[46] = 0.0;
  kData[47] = 0.0;
  kData[48] = GJ;
								
  return err;
}

int
RCFTFiberSection3D::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;

  // create an id to send objects tag and numFibers, 
  //     size 3 so no conflict with matData below if just 1 fiber
  static ID data(3);
  data(0) = this->getTag();
  data(1) = numFibers;
  int dbTag = this->getDbTag();
  res += theChannel.sendID(dbTag, commitTag, data);
  if (res < 0) {
    opserr << "RCFTFiberSection3D::sendSelf - failed to send ID data\n";
    return res;
  }    

  if (numFibers != 0) {
    
    // create an id containingg classTag and dbTag for each material & send it
    ID materialData(2*numFibers);
    for (int i=0; i<numFibers; i++) {
      UniaxialMaterial *theMat = theMaterials[i];
      materialData(2*i) = theMat->getClassTag();
      int matDbTag = theMat->getDbTag();
      if (matDbTag == 0) {
	  matDbTag = theChannel.getDbTag();
	  if (matDbTag != 0)
	  theMat->setDbTag(matDbTag);
      }
      materialData(2*i+1) = matDbTag;
    }    
    
    res += theChannel.sendID(dbTag, commitTag, materialData);
    if (res < 0) {
     opserr << "RCFTFiberSection3D::sendSelf - failed to send material data\n";
     return res;
    }    

    // send the fiber data, i.e. area and loc
    Vector fiberData(matData, 3*numFibers);
    res += theChannel.sendVector(dbTag, commitTag, fiberData);
    if (res < 0) {
     opserr << "RCFTFiberSection3D::sendSelf - failed to send material data\n";
     return res;
    }    

    // now invoke send(0 on all the materials
    for (int j=0; j<numFibers; j++)
      theMaterials[j]->sendSelf(commitTag, theChannel);
  }

  return res;
}

int
RCFTFiberSection3D::recvSelf(int commitTag, Channel &theChannel,
			 FEM_ObjectBroker &theBroker)
{
  int res = 0;

  static ID data(3);
  
  int dbTag = this->getDbTag();
  res += theChannel.recvID(dbTag, commitTag, data);

  if (res < 0) {
   opserr << "RCFTFiberSection3D::sendSelf - failed to recv ID data\n";
   return res;
  } 
   
  this->setTag(data(0));

  // recv data about materials objects, classTag and dbTag
  if (data(1) != 0) {
    ID materialData(2*data(1));
    res += theChannel.recvID(dbTag, commitTag, materialData);
    if (res < 0) {
     opserr << "RCFTFiberSection3D::sendSelf - failed to send material data\n";
     return res;
    }    

    // if current arrays not of correct size, release old and resize
    if (theMaterials == 0 || numFibers != data(1)) {
      // delete old stuff if outa date
      if (theMaterials != 0) {
	  for (int i=0; i<numFibers; i++)
	  delete theMaterials[i];
	  delete [] theMaterials;
	  if (matData != 0)
	  delete [] matData;
	  matData = 0;
	  theMaterials = 0;
      }

      // create memory to hold material pointers and fiber data
      numFibers = data(1);
      if (numFibers != 0) {

		theMaterials = new UniaxialMaterial *[numFibers];
	
		if (theMaterials == 0) {
		opserr << "RCFTFiberSection3D::recvSelf -- failed to allocate Material pointers\n";
		exit(-1);
		}

		for (int j=0; j<numFibers; j++)
		theMaterials[j] = 0;
	
		matData = new double [numFibers*3];

		if (matData == 0) {
		opserr << "RCFTFiberSection3D::recvSelf  -- failed to allocate double array for material data\n";
		exit(-1);
		}
      }
    }

    Vector fiberData(matData, 3*numFibers);
    res += theChannel.recvVector(dbTag, commitTag, fiberData);
    if (res < 0) {
     opserr << "RCFTFiberSection3D::sendSelf - failed to send material data\n";
     return res;
    }    
    
    int i;
    for (i=0; i<numFibers; i++) {
      int classTag = materialData(2*i);
      dbTag = materialData(2*i+1);

      // if material pointed to is blank or not of corrcet type, 
      // release old and create a new one
      if (theMaterials[i] == 0)
	  theMaterials[i] = theBroker.getNewUniaxialMaterial(classTag);
	  else if (theMaterials[i]->getClassTag() != classTag) {
		delete theMaterials[i];
		theMaterials[i] = theBroker.getNewUniaxialMaterial(classTag);      
      }

      if (theMaterials[i] == 0) {
		opserr << "RCFTFiberSection3D::recvSelf -- failed to allocate double array for material data\n";
		exit(-1);
      }

      theMaterials[i]->setDbTag(dbTag);
      res += theMaterials[i]->recvSelf(commitTag, theChannel, theBroker);
    }

    double Qz = 0.0;
    double Qy = 0.0;
    double A  = 0.0;
    double yLoc, zLoc, Area;

    // Recompute centroid
    for (i = 0; i < numFibers; i++) {
      yLoc = -matData[3*i];
      zLoc = matData[3*i+1];
      Area = matData[3*i+2];
      A  += Area;
      Qz += yLoc*Area;
      Qy += zLoc*Area;
    }
    
    yBar = -Qz/A;
    zBar = Qy/A;
  }    

  return res;
}

void
RCFTFiberSection3D::Print(OPS_Stream &s, int flag)
{
  if (flag == 2) {
    for (int i = 0; i < numFibers; i++) {
      s << -matData[3*i] << " "  << matData[3*i+1] << " "  << matData[3*i+2] << " " ;
      s << theMaterials[i]->getStress() << " "  << theMaterials[i]->getStrain() << endln;
    } 
  } else {
    s << "\nRCFTFiberSection3D, tag: " << this->getTag() << endln;
    s << "\tSection code: " << code;
    s << "\tNumber of Fibers: " << numFibers << endln;
    s << "\tCentroid: (" << -yBar << ", " << zBar << ')' << endln;
    
    if (flag == 1) {
      for (int i = 0; i < numFibers; i++) {
	  s << "\nLocation (y, z) = (" << -matData[3*i] << ", " << matData[3*i+1] << ")";
	  s << "\nArea = " << matData[3*i+2] << endln;
      theMaterials[i]->Print(s, flag);
      }
    }
  }
}

Response*
RCFTFiberSection3D::setResponse(const char **argv, int argc, OPS_Stream &output)
{
    // See if the response is one of the defaults
    Response *res = SectionForceDeformation::setResponse(argv, argc, output);
    if (res != 0)
      return res;
  
    // Check if fiber response is requested
    else if (strcmp(argv[0],"fiber") == 0) {
    int key = 0;
    int passarg = 2;
    
    if (argc <= 2)          // not enough data input
      return 0;
    else if (argc <= 3)		// RCFTFiber number was input directly
      key = atoi(argv[1]);
    else {                  // RCFTFiber near-to coordinate specified
      double yCoord = atof(argv[1]);
      double zCoord = atof(argv[2]);
      double ySearch = -matData[0];
      double zSearch =  matData[1];
      double closestDist = sqrt( pow(ySearch-yCoord,2) +
                                 pow(zSearch-zCoord,2) );
      double distance;
      for (int j = 1; j < numFibers; j++) {
	  ySearch = -matData[3*j];
  	  zSearch =  matData[3*j+1];
	  distance = sqrt( pow(ySearch-yCoord,2) +
			  pow(zSearch-zCoord,2) );
	  if (distance < closestDist) {
	  closestDist = distance;
	  key = j;
	  }
	  }
      passarg = 3;
	}
    
    if (key < numFibers)
      return theMaterials[key]->setResponse(&argv[passarg],argc-passarg,output);
    else
      return 0;
  }
  
  // otherwise response quantity is unknown for the RCFTFiberSection class
  else
  return 0;
}


int 
RCFTFiberSection3D::getResponse(int responseID, Information &sectInfo)
{
  // Just call the base class method ... don't need to define
  // this function, but keeping it here just for clarity
  return SectionForceDeformation::getResponse(responseID, sectInfo);
}

int
RCFTFiberSection3D::setParameter (const char **argv, int argc, Parameter &param)
{
	// Initial declarations
	int ok = -1;

	// A material parameter
	if (strcmp(argv[0],"material") == 0) {

		// Get the tag of the material
		int paramMatTag = atoi(argv[1]);

		// Loop over fibers to find the right material(s)
		for (int i=0; i<numFibers; i++) {
			if (paramMatTag == theMaterials[i]->getTag()) {
				ok = theMaterials[i]->setParameter(&argv[2], argc-2, param);
			}
		}
		if (ok<0) {
			opserr << "RCFTFiberSection3D::setParameter() - could not set parameter. " << endln;
			return -1;
		}
		else {
			return ok + 100;
		}
	} 
	else
		return -1;
}

int
RCFTFiberSection3D::updateParameter (int parameterID, Information &info)
{
	int ok = -1;

	switch (parameterID) {
	case 1:
		return -1;
	default:
		if (parameterID >= 100) {
			ID *paramIDPtr;
			paramIDPtr = info.theID;
			ID paramID = (*paramIDPtr);
			int paramMatrTag = paramID(1);

			for (int i=0; i<numFibers; i++) {
				if (paramMatrTag == theMaterials[i]->getTag()) {
					ok =theMaterials[i]->updateParameter(parameterID-100, info);
				}
			}
			if (ok < 0) {
				opserr << "RCFTFiberSection3D::updateParameter() - could not update parameter. " << endln;
				return ok;
			}
			else {
				return ok;
			}
		}
		else
			return -1;
	}
}


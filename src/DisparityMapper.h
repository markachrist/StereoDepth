/*
 * DisparityMapper.h
 *
 *  Created on: Oct 25, 2012
 *      Author: mchristopher
 */

#include <math.h>
#include <string.h>
#include <float.h>
#include "itkImageRegionIterator.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkDivideImageFilter.h"
#include "itkDivideOrZeroOutImageFilter.h"
#include "itkInterpolateImagePointsFilter.h"
#include "itkGradientImageFilter.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkCovariantVector.h"
#include "itkNthElementImageAdaptor.h"
#include "itkBinaryMagnitudeImageFilter.h"
#include "itkAtan2ImageFilter.h"
#include "itkAtanImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkOrImageFilter.h"
#include "itkNotImageFilter.h"
#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkThresholdImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkChangeInformationImageFilter.h"
#include "ImageUtils.h"
#include "RGBUtils.h"

#ifndef DISPARITYMAPPER_H_
#define DISPARITYMAPPER_H_

/**
 * Computes 3D information from a stereo image pair by finding the disparity between
 * matching pixels in the left and right images.
 *
 * Uses the method described in the following publication:
 *
 * Tang L, Garvin MK, Lee K, Alward WLM, Kwon YH, Abrˆmoff MD. Robust Multi-Scale Stereo Matching
 * from Fundus Images with Radiometric Differences. 2010 IEEE Transactions on Pattern Analysis
 * and Machine Intelligence.
 *
 */
class DisparityMapper {

protected:

	//Method parameters
	int minSize;
	int maxIt;
	int iHalfWidth;
	int gHalfWidth;
	int noiseParam;
	int intCorrectSize;
	double dispRange[2];
	double ksigma;
	double precision;
	double corrThreshold;

	bool saveIntermediate;

	std::string intermediateDir;

	/** Left image of color stereo pair */
	RGBUtils left;

	/** Left image of color stereo pair */
	RGBUtils right;

	/** Get min of two integers */
	inline int min(int a, int b){
		return (a < b) ? a : b;
	}

	/** Get max of two integers */
	inline int max(int a, int b){
		return (a > b) ? a : b;
	}

	/** Get min of two doubles */
	inline double min(double a, double b){
		return (a > b) ? a : b;
	}

	/** Get max of two doubles */
	inline double max(double a, double b){
		return (a > b) ? a : b;
	}

	std::vector<FloatImageType::Pointer> gradient(FloatImageType::Pointer i);
	FloatImageType::Pointer createFilledImage(int w, int h, double v);
	FloatImageType::Pointer createIndexImage(int w, int h, int dim);
	FloatImageType::Pointer createSubpixelIndexImage(int w, int h, int dim, double inc);
	FloatImageType::Pointer getDownsampleImage(FloatUtils i, int iterNum);
	ByteImageType::Pointer findOutOfRange(FloatImageType::Pointer i, double mx, double mn);
	FloatImageType::Pointer interpolateImagePoints(FloatImageType::Pointer xcoords,
			FloatImageType::Pointer ycoords, FloatImageType::Pointer zimage);
	FloatImageType::Pointer interpolateData(FloatImageType::Pointer x, FloatImageType::Pointer y,
			FloatImageType::Pointer z, FloatImageType::Pointer xq, FloatImageType::Pointer yq);
	double getSubpixelMaximum(FloatImageType::Pointer image, double start, double stop, double inc,
			itk::Point<double, 2> &loc);
	void replaceWithNearest(FloatImageType::Pointer disp, FloatImageType::Pointer score,
			ByteImageType::Pointer mask);
	FloatImageType::Pointer correctIntensity(FloatUtils &disp, FloatUtils &dispVert,
			FloatUtils &left, FloatUtils &right, int nsize, int padW, int padH);


	static double avg(std::vector<double> &v);
	static double correlation(std::vector<double> &v1, std::vector<double> &v2);
	static void abs(std::vector<double> &v);
	static double max(std::vector<double> &v, int *loc);
	static void shift(std::vector<double> &v, double s);
	static void exp(std::vector<double> &v, double s = 1.0);
	static void getPixelValues(itk::NeighborhoodIterator<FloatImageType> &nhood, std::vector<double> &dest);
	static void insertData(std::vector<double> &into, std::vector<double> &from, int idx);

	//Debugging methods
	FloatImageType::Pointer getNhoodImage(itk::NeighborhoodIterator<FloatImageType> &nhood);
	void printVector(std::vector<double> v);
	void fillImage(FloatImageType::Pointer image, double data[]);

public:
	DisparityMapper(RGBImageType::Pointer l, RGBImageType::Pointer r);
	virtual ~DisparityMapper();

//	bool computeDisparity(FloatImageType::Pointer &horz, FloatImageType::Pointer &vert);
	FloatImageType::Pointer computeDisparity();

	//TODO - getters/setters
	inline void setIntermediateDir(std::string &d){
		this->intermediateDir = d;
	}
	inline std::string & getIntermediateDir(){
		return this->intermediateDir;
	}

	inline void setShouldSaveIntermediateResults(bool s){
		this->saveIntermediate = s;
	}
	inline bool shouldSaveIntermediateResults(){
		return this->saveIntermediate;
	}
};

#endif /* DISPARITYMAPPER_H_ */

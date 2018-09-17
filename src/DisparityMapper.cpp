/*
 * DisparityMapper.cpp
 *
 *  Created on: Oct 25, 2012
 *      Author: mchristopher
 */

#include "DisparityMapper.h"

/**
 * Create object to map pixel disparities for the given stereo pair images.
 *
 * @param l
 *   First image of the stereo pair
 * @param r
 *   The other image of the stereo pair
 */
DisparityMapper::DisparityMapper(RGBImageType::Pointer l, RGBImageType::Pointer r)
: left(l), right(r){

	//Default parameter values
	precision = 0.1;
	minSize = 8;
	maxIt = 100;
	iHalfWidth = 3;
	gHalfWidth = 4;
	dispRange[0] = -1.2; dispRange[1] = 1.2;
	ksigma = 1.2;
	corrThreshold = 0.9;
	noiseParam = 40;
	intCorrectSize = 3;

	saveIntermediate = true;
}

/**
 * Applies multiscale algorithm to estimate disparity (depth) from the given stereo image pair. Output is float-valued
 * image with each pixel proportional to the estimated depth at that position.
 *
 * Uses the method described in the following publication:
 *
 * Tang L, Garvin MK, Lee K, Alward WLM, Kwon YH, Abrˆmoff MD. Robust Multi-Scale Stereo Matching
 * from Fundus Images with Radiometric Differences. 2010 IEEE Transactions on Pattern Analysis
 * and Machine Intelligence.
 *
 * @param horz
 *   Output destination for horizontal disparity
 * @param vert
 *   Output destination for vertical disparity
 * @return
 *   True if disparity computed successfully, false otherwise
 */
//bool DisparityMapper::computeDisparity(FloatImageType::Pointer &horz, FloatImageType::Pointer &vert){
FloatImageType::Pointer DisparityMapper::computeDisparity(){

	bool result = true;
	int numIts, minDim;
	int correctNsize = this->intCorrectSize;
	int width = this->left.getWidth();
	int height = this->left.getHeight();

	FloatUtils l(this->left.convertToGrayscale());
	FloatUtils r(this->right.convertToGrayscale());

	//Holds final disparity map
	FloatUtils finalDispMap;
	FloatUtils finalDispMapVert;

	//Holds the estimated disparity map, updated at each scale
	FloatUtils disparityMap;
	FloatUtils disparityMapVert;
	FloatUtils lastDispMap;
	FloatUtils lastDispMapVert;

	//Holds intensity correction factors
	FloatUtils intCorrect;
	FloatUtils intCorrectMir;

	//Find the number of iterations
	numIts = (int)(ceil(log2(this->min(l.getWidth(), l.getHeight())) - log2(minSize)));

	//At each scale
	for(int si = numIts; si >= 0; --si){
		double downScale = pow(2, si);
		double sigma = si*this->ksigma;

		int padW, padH;
		int curW = (int)(width/downScale + 0.5);
		int curH = (int)(height/downScale + 0.5);

		std::cout << "Working on scale " << si << " (" << curW << ", " << curH << ")" << std::endl;

		//Initialize output store on first iteration
		if(si == numIts){
			disparityMap.setImage(this->createFilledImage(curW, curH, 0.0));
			disparityMapVert.setImage(this->createFilledImage(curW, curH, 0.0));

			intCorrect.setImage(this->createFilledImage(curW, curH, 0.0));
			intCorrectMir.setImage(this->createFilledImage(curW, curH, 0.0));
		}

		//Initialize objects to hold images/disp maps for this scale
		FloatUtils curL(l.getImage());
		FloatUtils curR(r.getImage());
		std::vector<FloatUtils> dispMaps;
		std::vector<FloatUtils> dispMapsVert;
		std::vector<FloatUtils> scoreMaps;

		//For each image as the reference
		for(int ii = 0; ii < 2; ++ii){

			FloatUtils tempDisp(this->createFilledImage(curW, curH, 0.0));
			dispMaps.push_back(tempDisp);
			FloatUtils tempDispVert(this->createFilledImage(curW, curH, 0.0));
			dispMapsVert.push_back(tempDispVert);
			FloatUtils tempScore(this->createFilledImage(curW, curH, 0.0));
			scoreMaps.push_back(tempScore);

			//Use right image as reference on second iteration
			if(ii == 1){
				std::cout << "\tUsing right as ref" << std::endl;
				curL.setImage(r.getImage());
				curR.setImage(l.getImage());
			}
			else{
				std::cout << "\tUsing left as ref" << std::endl;
			}

			//Create scaled reference and matching images
			if(si != 0){

				curL.setImage(curL.blur((int)(3.0*sigma), sigma));
				curL.setImage(curL.resize(curW, curH));

				curR.setImage(curR.blur((int)(3.0*sigma), sigma));
				curR.setImage(curR.resize(curW, curH));
			}

			//Get disparity map from previous (coarser) scale, if it exists
			if(si != numIts){

				//Adjust disparity from last scale
				if(ii == 0){

					disparityMap.setImage(lastDispMap.mult(2.0));
					disparityMap.setImage(disparityMap.resize(curW, curH));
					disparityMapVert.setImage(lastDispMapVert.mult(2.0));
					disparityMapVert.setImage(disparityMapVert.resize(curW, curH));

					intCorrect.setImage(intCorrect.resize(curW, curH));

				}
				//When reference/matching images are reversed, need to invert and interpolate
				//data for disparity map going in reverse direction
				else{

					FloatImageType::Pointer xcoords = this->createIndexImage(lastDispMap.getWidth(),
							lastDispMap.getHeight(), 0);
					FloatImageType::Pointer ycoords = this->createIndexImage(lastDispMap.getWidth(),
							lastDispMap.getHeight(), 1);
					FloatUtils xcoordsAdd(xcoords);
					xcoordsAdd.setImage(xcoordsAdd.add(lastDispMap.getImage()));
					disparityMap.setImage(this->interpolateImagePoints(xcoordsAdd.getImage(),
							ycoords, lastDispMap.mult(-1.0)));

					xcoordsAdd.setImage(xcoords);
					xcoordsAdd.setImage(xcoordsAdd.add(lastDispMapVert.getImage()));
					disparityMapVert.setImage(this->interpolateImagePoints(xcoordsAdd.getImage(),
							ycoords, lastDispMapVert.mult(-1.0)));

					disparityMap.setImage(disparityMap.mult(2.0));
					disparityMap.setImage(disparityMap.resize(curW, curH));
					disparityMapVert.setImage(disparityMapVert.mult(2.0));
					disparityMapVert.setImage(disparityMapVert.resize(curW, curH));
				}

				//Find nan/inf values
				double maxPos = this->max(curW, curH);
				ByteImageType::Pointer mask = this->findOutOfRange(disparityMap.getImage(), maxPos, -maxPos);
				this->replaceWithNearest(disparityMap.getImage(), disparityMap.getImage(), mask);
				mask = this->findOutOfRange(disparityMapVert.getImage(), maxPos, -maxPos);
				this->replaceWithNearest(disparityMapVert.getImage(), disparityMapVert.getImage(), mask);
			}

			//Pad boundaries of images using mirroring to avoid edge effects during
			//search for corresponding pixels
			double mx = fabs(disparityMap.max());
			double mn = fabs(disparityMap.min());
			double maxDisp = (mx > mn) ? mx : mn;

			mx = fabs(disparityMapVert.max());
			mn = fabs(disparityMapVert.min());
			double maxDispVert = (mx > mn) ? mx : mn;

			//Find necessary padding
			padW = max(this->iHalfWidth, this->gHalfWidth) +
				   (int)(ceil(max(fabs(this->dispRange[0]), fabs(this->dispRange[0])))) +
				   (int)(ceil(maxDisp));
			padW = (padW >= curW) ? curW - 1: padW;

			padH = max(this->iHalfWidth, this->gHalfWidth) +
				   (int)(ceil(max(fabs(this->dispRange[0]), fabs(this->dispRange[0])))) +
				   (int)(ceil(maxDispVert));
			padH = (padH >= curH) ? curH - 1: padH;

//			std::cout << "\tFuck 10 " << padW << ", " << padH << std::endl;

			std::cout << "padW = " << padW << std::endl;
			std::cout << "padH = " << padH << std::endl;

			//Apply mirror padding
			curL.setImage(curL.mirrorPadImage(padW, padW, padH, padH));
			curR.setImage(curR.mirrorPadImage(padW, padW, padH, padH));
			intCorrectMir.setImage(intCorrect.mirrorPadImage(padW, padW, padH, padH));
			FloatUtils disparityMapMir(disparityMap.mirrorPadImage(padW, padW, padH, padH));
			FloatUtils disparityMapMirVert(disparityMapVert.mirrorPadImage(padW, padW, padH, padH));

//			std::cout << "\tFuck 11" << std::endl;

			//Perform intensity correction on the reference image
			FloatImageType::Pointer toAdd = intCorrectMir.getImage();
			if(ii == 0){
				curL.setImage(curL.add(intCorrectMir.getImage()));
			}
			else{
				curR.setImage(curR.add(intCorrectMir.getImage()));
			}

			std::vector<FloatImageType::Pointer> gradsL = this->gradient(curL.getImage());
			std::vector<FloatImageType::Pointer> gradsR = this->gradient(curR.getImage());

			//Now perform the actual search to determine disparities
			itk::Index<2> start;
			start[0] = padW; start[1] = padH;
			itk::Size<2> size;
			size[0] = curW; size[1] = curH;
			FloatImageType::RegionType region;
			region.SetIndex(start);
			region.SetSize(size);

			//Create output iterators
			itk::ImageRegionIterator<FloatImageType> outIt(dispMaps[ii].getImage(),
					dispMaps[ii].getImage()->GetLargestPossibleRegion());
			itk::ImageRegionIterator<FloatImageType> outItVert(dispMapsVert[ii].getImage(),
					dispMapsVert[ii].getImage()->GetLargestPossibleRegion());
			itk::ImageRegionIterator<FloatImageType> scoreIt(scoreMaps[ii].getImage(),
					scoreMaps[ii].getImage()->GetLargestPossibleRegion());

			//Create a bunch of iterators to access pixels & pixel neighborhoods
			itk::Size<2> intSize;
			intSize[0] = iHalfWidth;
			intSize[1] = iHalfWidth;
			itk::NeighborhoodIterator<FloatImageType> intIt(intSize, curL.getImage(), region);
			itk::NeighborhoodIterator<FloatImageType> dispIntIt(intSize, disparityMapMir.getImage(), region);
			itk::NeighborhoodIterator<FloatImageType> dispIntItVert(intSize, disparityMapMirVert.getImage(), region);

			itk::Size<2> gradSize;
			gradSize[0] = gHalfWidth;
			gradSize[1] = gHalfWidth;
			itk::NeighborhoodIterator<FloatImageType> gradMIt(gradSize, gradsL[0], region);
			itk::NeighborhoodIterator<FloatImageType> gradDIt(gradSize, gradsL[1], region);
			itk::NeighborhoodIterator<FloatImageType> dispGradIt(gradSize, disparityMapMir.getImage(), region);
			itk::NeighborhoodIterator<FloatImageType> dispGradItVert(gradSize, disparityMapMirVert.getImage(), region);

			//Create image to store temp results for each search neighborhood
			itk::Size<2> searchSize;
			searchSize[0] = (int)((ceil(dispRange[1]) - floor(dispRange[0]) + 1));
			searchSize[1] = (int)((ceil(dispRange[1]) - floor(dispRange[0]) + 1));
			FloatUtils nhoodDisp(this->createFilledImage((int)searchSize[0], (int)searchSize[1], 0));
			FloatUtils nhoodDispVert(this->createFilledImage((int)searchSize[0], (int)searchSize[1], 0));

			//Length of the feature vector used to compare pixels
			int length = (2*iHalfWidth + 1)*(2*iHalfWidth + 1) + (2*gHalfWidth + 1)*(2*gHalfWidth + 1);
			std::vector<double> refFeature(length, 0);
			std::vector<double> otherFeature(length, 0);
			std::vector<double> tempI((2*iHalfWidth + 1)*(2*iHalfWidth + 1), 0);
			std::vector<double> weightI((2*iHalfWidth + 1)*(2*iHalfWidth + 1), 0);
			std::vector<double> tempG((2*gHalfWidth + 1)*(2*gHalfWidth + 1), 0);
			std::vector<double> weightG((2*gHalfWidth + 1)*(2*gHalfWidth + 1), 0);
			std::vector<double> gradD((2*gHalfWidth + 1)*(2*gHalfWidth + 1), 0);

			outIt.GoToBegin();
			outItVert.GoToBegin();
			scoreIt.GoToBegin();
			gradMIt.GoToBegin();
			gradDIt.GoToBegin();
			dispIntIt.GoToBegin();
			dispGradIt.GoToBegin();
			dispIntItVert.GoToBegin();
			dispGradItVert.GoToBegin();
			for(intIt.GoToBegin(); !intIt.IsAtEnd(); ++intIt){

				itk::Index<2> refIdx = intIt.GetIndex();
				double curDisp = dispGradIt.GetCenterPixel();
				double curDispVert = dispGradItVert.GetCenterPixel();

				//Extract features from images
				DisparityMapper::getPixelValues(intIt, tempI);
				DisparityMapper::getPixelValues(gradMIt, tempG);
				DisparityMapper::getPixelValues(gradDIt, gradD);

				DisparityMapper::insertData(refFeature, tempI, 0);
				DisparityMapper::insertData(refFeature, tempG, tempI.size());

				//Get distance based weighting
				if(refIdx[1] - 2*padH >= 0 && refIdx[1] - padH + padH < curH &&
				   refIdx[0] - padW - padH >= 0 && refIdx[0] - padH + padH < curW){
					DisparityMapper::getPixelValues(dispIntIt, weightI);
					DisparityMapper::shift(weightI, curDisp);
					DisparityMapper::abs(weightI);
					DisparityMapper::exp(weightI, -1.0);
					DisparityMapper::getPixelValues(dispGradIt, weightG);
					DisparityMapper::shift(weightG, curDisp);
					DisparityMapper::abs(weightG);
					DisparityMapper::exp(weightG, -1.0);
				}
				//Use constant weight around image borders
				else{
					weightI.assign(tempI.size(), 1.0);
					weightG.assign(tempG.size(), 1.0);
				}

				for(int i = 0; i < weightI.size(); ++i){
					refFeature[i] *= weightI[i];
				}
				for(int i = 0; i < weightG.size(); ++i){
					refFeature[i + weightI.size()] *= weightG[i] * gradD[i];
				}

				//Define region to search for correspondence in the target image
				int dispOffsetH = round(curDisp);
				int dispOffsetV = round(curDispVert);

				itk::Index<2> searchIdx;
				searchIdx[0] = refIdx[0] + floor(dispRange[0]) + dispOffsetH;
				searchIdx[1] = refIdx[1] + floor(dispRange[0]) + dispOffsetV;
				FloatImageType::RegionType searchRegion(searchIdx, searchSize);

				//Create search iterators
				itk::NeighborhoodIterator<FloatImageType> intItR(intSize, curR.getImage(), searchRegion);
				itk::NeighborhoodIterator<FloatImageType> gradMItR(gradSize, gradsR[0], searchRegion);
				itk::NeighborhoodIterator<FloatImageType> gradDItR(gradSize, gradsR[1], searchRegion);

				//Create output iterators
				itk::ImageRegionIterator<FloatImageType> nhoodDispIt(nhoodDisp.getImage(),
						nhoodDisp.getImage()->GetLargestPossibleRegion());
				itk::ImageRegionIterator<FloatImageType> nhoodDispItVert(nhoodDispVert.getImage(),
						nhoodDispVert.getImage()->GetLargestPossibleRegion());

				gradMItR.GoToBegin();
				gradDItR.GoToBegin();
				nhoodDispIt.GoToBegin();
				int testi = 0;
				for(intItR.GoToBegin(); !intItR.IsAtEnd(); ++intItR){

					itk::Index<2> curSearchIdx = intItR.GetIndex();

					//Extract features from images
					DisparityMapper::getPixelValues(intItR, tempI);
					DisparityMapper::getPixelValues(gradMItR, tempG);
					DisparityMapper::getPixelValues(gradDItR, gradD);

					DisparityMapper::insertData(otherFeature, tempI, 0);
					DisparityMapper::insertData(otherFeature, tempG, tempI.size());

					for(int i = 0; i < weightI.size(); ++i){
						otherFeature[i] *= weightI[i];
					}
					for(int i = 0; i < weightG.size(); ++i){
						otherFeature[i + weightI.size()] *= weightG[i] * gradD[i];
					}

					double curScore = DisparityMapper::correlation(refFeature, otherFeature);
					nhoodDispIt.Set(curScore);

					++gradMItR;
					++gradDItR;
					++nhoodDispIt;
					++testi;
				}

				//Interpolate values to get subpixel precision
				itk::Point<double, 2> maxLoc;
				double offset = floor(dispRange[0]);
				double maxScore = this->getSubpixelMaximum(nhoodDisp.getImage(), dispRange[0] - offset, dispRange[1] - offset,
						precision, maxLoc);
				maxLoc[0] += offset;
				maxLoc[1] += offset;

				//Save horizontal and vertical disparity to output destinations
				double result = maxLoc[0];
				result = (fabs(result) < 0.00005) ? 0 : result;
				outIt.Set(result + dispOffsetH);
				result = maxLoc[1];
				result = (fabs(result) < 0.00005) ? 0 : result;
				outItVert.Set(result + dispOffsetV);

				scoreIt.Set(maxScore);

				++outIt;
				++outItVert;
				++scoreIt;
				++gradMIt;
				++gradDIt;
				++dispIntIt;
				++dispGradIt;
				++dispIntItVert;
				++dispGradItVert;
			}

			{//Find pixels w/ disparities > 3 st. devs. from mean
				itk::StatisticsImageFilter<FloatImageType>::Pointer stats =
						itk::StatisticsImageFilter<FloatImageType>::New();
				stats->SetInput(dispMaps[ii].getImage());
				stats->Update();

				double mean = stats->GetMean();
				double sd = stats->GetSigma();
				double upper = mean + 3*sd;
				double lower = mean - 3*sd;

				itk::ThresholdImageFilter<FloatImageType>::Pointer thresh1 =
						itk::ThresholdImageFilter<FloatImageType>::New();
				thresh1->SetLower(lower);
				thresh1->SetUpper(upper);
				thresh1->SetOutsideValue(FLT_MAX);
				thresh1->SetInput(dispMaps[ii].getImage());
				thresh1->Update();

				itk::ThresholdImageFilter<FloatImageType>::Pointer thresh2 =
						itk::ThresholdImageFilter<FloatImageType>::New();
				thresh2->SetInput(thresh1->GetOutput());
				thresh2->SetLower(FLT_MAX);
				thresh2->SetOutsideValue(0.0);
				thresh2->Update();
				FloatUtils mask1(thresh2->GetOutput());

				//Find pixels < correlation threshold
				itk::ThresholdImageFilter<FloatImageType>::Pointer thresh =
						itk::ThresholdImageFilter<FloatImageType>::New();
				thresh->SetInput(scoreMaps[ii].getImage());
				thresh->SetLower(this->corrThreshold);
				thresh->SetOutsideValue(FLT_MAX);
				thresh->Update();

				thresh->SetInput(thresh->GetOutput());
				thresh->SetLower(FLT_MAX);
				thresh->SetOutsideValue(0.0);
				thresh->Update();
				FloatUtils mask2(thresh->GetOutput());

				//Combine the masks
				itk::OrImageFilter<ByteImageType>::Pointer combine =
						itk::OrImageFilter<ByteImageType>::New();
				combine->SetInput1(FloatUtils::floatToByteImage(mask1.getImage()));
				combine->SetInput2(FloatUtils::floatToByteImage(mask2.getImage()));
				combine->Update();

				itk::NotImageFilter<ByteImageType, ByteImageType>::Pointer notter =
						itk::NotImageFilter<ByteImageType, ByteImageType>::New();
				notter->SetInput(combine->GetOutput());
				notter->Update();
				ByteUtils mask(notter->GetOutput());

				//Replace masked out pixels in disp and vertical disp maps
				this->replaceWithNearest(dispMaps[ii].getImage(), scoreMaps[ii].getImage(), mask.getImage());
				this->replaceWithNearest(dispMapsVert[ii].getImage(), scoreMaps[ii].getImage(), mask.getImage());
			}

			if(saveIntermediate){
				std::ostringstream outName;
				outName << this->intermediateDir << "/" << si << "-i" << ii << ".vtk";
				std::string out(outName.str());
				dispMaps[ii].writeToPath(out);

				ByteUtils bdisp(FloatUtils::floatToByteImage(dispMaps[ii].getImage()));
				std::ostringstream boutName;
				boutName << this->intermediateDir << "/" << si << "-i" << ii << ".png";
				std::string bout(boutName.str());
				bdisp.writeToPath(bout);
			}
		}

		disparityMap.setImage(dispMaps[0].copyImage());
		disparityMapVert.setImage(dispMapsVert[0].copyImage());

		//Combine the two depth maps
		//Interpolate disparity and score values
		FloatImageType::Pointer xcoords = this->createIndexImage(curW, curH, 0);
		FloatImageType::Pointer ycoords = this->createIndexImage(curW, curH, 1);
		FloatUtils xcoordsAdd(xcoords);
		FloatUtils ycoordsAdd(ycoords);

//		std::cout << "\tFuck 20" << std::endl;

		xcoordsAdd.setImage(xcoordsAdd.add(dispMaps[1].getImage()));
		ycoordsAdd.setImage(ycoordsAdd.add(dispMapsVert[1].getImage()));

//		std::cout << "\tFuck 20 1" << std::endl;

		FloatImageType::Pointer interpDispMap = this->interpolateImagePoints(xcoordsAdd.getImage(),
				ycoordsAdd.getImage(), dispMaps[1].mult(-1.0));

//		std::cout << "\tFuck 20 2" << std::endl;

		FloatImageType::Pointer interpDispMapVert = this->interpolateImagePoints(xcoordsAdd.getImage(),
				ycoordsAdd.getImage(), dispMapsVert[1].mult(-1.0));

//		std::cout << "\tFuck 20 3" << std::endl;

		FloatImageType::Pointer interpScoreMap = this->interpolateImagePoints(xcoordsAdd.getImage(),
				ycoordsAdd.getImage(), scoreMaps[1].getImage());

//		std::cout << "\tFuck 21" << std::endl;

		//Use disparity with best score from each map
		disparityMap.setImage(dispMaps[0].getImage());
		itk::ImageRegionIterator<FloatImageType> outIt(disparityMap.getImage(),
				disparityMap.getImage()->GetLargestPossibleRegion());
		itk::ImageRegionIterator<FloatImageType> outItVert(disparityMapVert.getImage(),
				disparityMapVert.getImage()->GetLargestPossibleRegion());
		itk::ImageRegionIterator<FloatImageType> altIt(interpDispMap,
				interpDispMap->GetLargestPossibleRegion());
		itk::ImageRegionIterator<FloatImageType> altItVert(interpDispMapVert,
				interpDispMapVert->GetLargestPossibleRegion());
		itk::ImageRegionIterator<FloatImageType> leftIt(scoreMaps[0].getImage(),
				scoreMaps[0].getImage()->GetLargestPossibleRegion());
		itk::ImageRegionIterator<FloatImageType> rightIt(interpScoreMap,
				interpScoreMap->GetLargestPossibleRegion());

//		std::cout << "\tFuck 22" << std::endl;

		int altCount = 0;
		outItVert.GoToBegin();
		altIt.GoToBegin();
		altItVert.GoToBegin();
		leftIt.GoToBegin();
		rightIt.GoToBegin();
		for(outIt.GoToBegin(); !outIt.IsAtEnd(); ++outIt){

			double lScore  = leftIt.Get();
			double rScore  = rightIt.Get();

			if(rScore > lScore){
				outIt.Set(altIt.Get());
				outItVert.Set(altItVert.Get());
				++altCount;
			}

			++outItVert;
			++altIt;
			++altItVert;
			++leftIt;
			++rightIt;
		}

//		std::cout << "\tFuck 23" << std::endl;

		double altFraction = ((double)altCount)/(curW*curH);

		if(saveIntermediate){
			std::cout << "% from 2nd map: " << (100*altFraction) << std::endl;

			std::ostringstream outName;
			outName << this->intermediateDir << "/" << si << ".vtk";
			std::string out(outName.str());
			disparityMap.writeToPath(out);

			ByteUtils bdisp(FloatUtils::floatToByteImage(disparityMap.getImage()));
			std::ostringstream boutName;
			boutName << this->intermediateDir << "/" << si << ".png";
			std::string bout(boutName.str());
			bdisp.writeToPath(bout);
		}

		//TODO - weiner filter for noise removal
		//For now, just remove out of bounds values
		itk::StatisticsImageFilter<FloatImageType>::Pointer stats =
				itk::StatisticsImageFilter<FloatImageType>::New();
		stats->SetInput(disparityMap.getImage());
		stats->Update();

		double mean = stats->GetMean();
		double sd = stats->GetSigma();
		double upper = mean + 3*sd;
		double lower = mean - 3*sd;
//		double upper = width;
//		double lower = -upper;

		itk::ThresholdImageFilter<FloatImageType>::Pointer thresh1 =
				itk::ThresholdImageFilter<FloatImageType>::New();
		thresh1->SetLower(lower);
		thresh1->SetUpper(upper);
		thresh1->SetOutsideValue(FLT_MAX);
		thresh1->SetInput(disparityMap.getImage());
		thresh1->Update();

		itk::ThresholdImageFilter<FloatImageType>::Pointer thresh2 =
				itk::ThresholdImageFilter<FloatImageType>::New();
		thresh2->SetInput(thresh1->GetOutput());
		thresh2->SetLower(FLT_MAX);
		thresh2->SetOutsideValue(0.0);
		thresh2->Update();
		FloatUtils mask1(thresh2->GetOutput());

		itk::NotImageFilter<ByteImageType, ByteImageType>::Pointer notter =
				itk::NotImageFilter<ByteImageType, ByteImageType>::New();
		notter->SetInput(FloatUtils::floatToByteImage(mask1.getImage()));
		notter->Update();
		ByteUtils mask(notter->GetOutput());

		//Replace masked out pixels in disp and vertical disp maps
		this->replaceWithNearest(disparityMap.getImage(), scoreMaps[0].getImage(), mask.getImage());
		this->replaceWithNearest(disparityMapVert.getImage(), scoreMaps[0].getImage(), mask.getImage());

		finalDispMap.setImage(disparityMap.copyImage());
		finalDispMapVert.setImage(disparityMapVert.copyImage());

		//Intensity correction (l & r have been switched...)
		if(si != 0){
			intCorrect.setImage(this->correctIntensity(disparityMap, disparityMapVert, curR, curL, correctNsize, padW, padH));
			correctNsize = min(100, (int)(round(2*correctNsize)));
		}

		lastDispMap.setImage(disparityMap.getImage());
		lastDispMapVert.setImage(disparityMapVert.getImage());
	}

//	FloatImageType::Pointer x = finalDispMap.mult(finalDispMap.getImage());
//	FloatImageType::Pointer y = finalDispMapVert.mult(finalDispMapVert.getImage());
//
//	horz = finalDispMap.mult(finalDispMap.getImage());
//	vert = finalDispMapVert.mult(finalDispMapVert.getImage());

//	FloatUtils result(x);
//	result.setImage(result.add(y));
//	result.setImage(result.sqrt());
//
//	result.writeToPath("/Users/mchristopher/Documents/Data/DepthAlignment/testpair/comb-final-disp.vtk");
//	finalDispMap.writeToPath("/Users/mchristopher/Documents/Data/DepthAlignment/testpair/horiz-final-disp.vtk");
//	finalDispMapVert.writeToPath("/Users/mchristopher/Documents/Data/DepthAlignment/testpair/vert-final-disp.vtk");
//
	return finalDispMap.mult(-1.0);

//	return result;
}

/**
 * Creates float image of the given size with every pixel value set to v.
 *
 * @param w
 *   Width of the image to create
 * @param h
 *   Height of the image to create
 * @param v
 *   Value with which to fill image
 * @return
 *   Float image of size w by h, filled with the value v
 */
FloatImageType::Pointer DisparityMapper::createFilledImage(int w, int h, double v){
	itk::Index<2> idx;
	idx.Fill(0);

	itk::Size<2> s;
	s[0] = w; s[1] = h;

	FloatImageType::RegionType region;
	region.SetIndex(idx);
	region.SetSize(s);

	FloatImageType::Pointer img = FloatImageType::New();
	img->SetRegions(region);
	img->Allocate();
	img->FillBuffer(v);

	return img;
}

/**
 * Linearly interpolate the points specified by the values in the xcoords and ycoords images from values
 * in zimage. The images xcoords and ycoords should have the same dimensions and the output will
 * will also have these dimensions.
 *
 * output[x, y] = interpolate( zimage[ xcoords[x,y], ycoords[x,y] ] )
 *
 * Similar to Matlab function griddata.
 *
 * @param xcoords
 *   x-coordinates of the points to interpolate
 * @param ycoords
 *   y-coordinates of the points to interpolate
 * @param zimage
 *   Image data used to perform interpolation
 * @return
 *   Image with values determined by interpolation of points specified by xcoords and ycoords
 */
FloatImageType::Pointer DisparityMapper::interpolateImagePoints(FloatImageType::Pointer xcoords,
		FloatImageType::Pointer ycoords, FloatImageType::Pointer zimage){

	int w = xcoords->GetLargestPossibleRegion().GetSize()[0];
	int h = xcoords->GetLargestPossibleRegion().GetSize()[1];
	FloatUtils out(this->createFilledImage(w, h, 0.0));

//	std::cout << "interp size = " << w << "," << h << std::endl;

	itk::LinearInterpolateImageFunction<FloatImageType>::Pointer interp =
			itk::LinearInterpolateImageFunction<FloatImageType>::New();
	interp->SetInputImage(zimage);

//	std::cout << "interp 2" << std::endl;

	FloatUtils mnmx(zimage);
	double max = mnmx.max();
	double min = mnmx.min();

//	std::cout << "interp minmax = " << min << "," << max << std::endl;
//	/*******************************************/
//	FloatUtils xstat(xcoords);
//	double xmax = xstat.max();
//	double xmin = xstat.min();
//	std::cout << "interp xrange = " << xmin << "," << xmax << std::endl;
//
//	FloatUtils ystat(ycoords);
//	double ymax = ystat.max();
//	double ymin = ystat.min();
//	std::cout << "interp yrange = " << ymin << "," << ymax << std::endl;
//	/*******************************************/

	itk::ImageRegionIterator<FloatImageType> xIt(xcoords, xcoords->GetLargestPossibleRegion());
	itk::ImageRegionIterator<FloatImageType> yIt(ycoords, ycoords->GetLargestPossibleRegion());
	itk::ImageRegionIterator<FloatImageType> outIt(out.getImage(), out.getImage()->GetLargestPossibleRegion());

//	std::cout << "interp 3" << std::endl;

	xIt.GoToBegin();
	yIt.GoToBegin();
	for(outIt.GoToBegin(); !outIt.IsAtEnd(); ++outIt){

		itk::Point<double, 2> p;
		p[0] = xIt.Get();
		p[1] = yIt.Get();

		p[0] = (p[0] < 0) ? 0 : p[0];
		p[0] = (p[0] > w - 1) ? w - 1 : p[0];
		p[1] = (p[1] < 0) ? 0 : p[1];
		p[1] = (p[1] > h - 1) ? h - 1 : p[1];

//		if(!interp->IsInsideBuffer(p)){
//			std::cout << "Problem idx = " << p << std::endl;
//		}

		double zi = interp->Evaluate(p);

		outIt.Set(zi);

		++xIt;
		++yIt;
	}

//	std::cout << "interp 4" << std::endl;

	//Get rid of NaN and out of range values
	ByteImageType::Pointer mask = this->findOutOfRange(out.getImage(), max, min);
	this->replaceWithNearest(out.getImage(), out.getImage(), mask);

//	std::cout << "interp 5" << std::endl;

	return out.getImage();
}

/**
 * Interpolates data for the function z defined at the coords given by x and y. The value
 * of the function at the query points (xq, yq) is determined by (nearest neighbor) interpolation.
 *
 * The images x, y, and z should all be of the same size as should xq and yq.
 *
 * Similar to Matlab function griddata called using the nearest neighbor interpolation method.
 *
 * TODO: Produce (delaunay) triangulation using Voronoi partition to allow better interpolation
 *
 * @param xs
 *   x-coordinates of the values in z
 * @param ys
 *   y-coordinates of the values in z
 * @param zs
 *   Image data used to perform interpolation
 * @param xqs
 *   x-coordinates of the points to interpolate
 * @param yqs
 *   y-coordinates of the points to interpolate
 * @return
 *   Image with values determined by interpolation of points specified by xqcoords and yqcoords
 */
FloatImageType::Pointer DisparityMapper::interpolateData(FloatImageType::Pointer xs, FloatImageType::Pointer ys,
		FloatImageType::Pointer zs, FloatImageType::Pointer xqs, FloatImageType::Pointer yqs){

	int w = xqs->GetLargestPossibleRegion().GetSize()[0];
	int h = yqs->GetLargestPossibleRegion().GetSize()[1];
	FloatUtils zq(w, h, 0.0);

	itk::ImageRegionIterator<FloatImageType> xqIt(xqs, xqs->GetLargestPossibleRegion());
	itk::ImageRegionIterator<FloatImageType> yqIt(yqs, yqs->GetLargestPossibleRegion());
	itk::ImageRegionIterator<FloatImageType> zqIt(zq.getImage(), zq.getImage()->GetLargestPossibleRegion());

	xqIt.GoToBegin();
	yqIt.GoToBegin();
	for(zqIt.GoToBegin(); !zqIt.IsAtEnd(); ++zqIt){

		double xq = xqIt.Get();
		double yq = yqIt.Get();
		double zq = 0.0;
		double dist = DBL_MAX;

		itk::ImageRegionIterator<FloatImageType> xIt(xs, xs->GetLargestPossibleRegion());
		itk::ImageRegionIterator<FloatImageType> yIt(ys, ys->GetLargestPossibleRegion());
		itk::ImageRegionIterator<FloatImageType> zIt(zs, zs->GetLargestPossibleRegion());

		xIt.GoToBegin();
		yIt.GoToBegin();
		for(zIt.GoToBegin(); !zIt.IsAtEnd(); ++zIt){

			double x = xIt.Get() - xq;
			x *= x;
			double y = yIt.Get() - yq;
			y *= y;
			double curDist = sqrt(x + y);

			if(curDist < dist){
				zq = zIt.Get();
				dist = curDist;
			}
			++xIt;
			++yIt;
		}

		zqIt.Set(zq);

		++xqIt;
		++yqIt;
	}

	return zq.getImage();
}


/**
 * Creates an image with each pixel value equal to the value of the its coordinate for the
 * dimension specified by dim.
 *
 * Ex. createIndexImage(4, 3, 0) :
 *
 * 0 1 2 3
 * 0 1 2 3
 * 0 1 2 3
 *
 * @param w
 *   Width of the created image
 * @param h
 *   Height of the created image
 * @param dim
 *   The dimension (0 = x, 1 = y) that should recorded in pixel values
 * @return
 *   Image with the specified dimension coordinate recorded in each pixel value
 */
FloatImageType::Pointer DisparityMapper::createIndexImage(int w, int h, int dim){

	FloatImageType::Pointer r = this->createFilledImage(w, h, 0.0);

	itk::ImageRegionIterator<FloatImageType> it(r, r->GetLargestPossibleRegion());

	for(it.GoToBegin(); !it.IsAtEnd(); ++it){
		it.Set(it.GetIndex()[dim]);
	}

	return r;
}


/**
 * Gets the appropriately downsampled image based on the iteration number.
 *
 * @param i
 *   The image to downsample
 * @param iterNum
 *   The iteration number for whcih the image should be retreived
 */
FloatImageType::Pointer DisparityMapper::getDownsampleImage(FloatUtils i, int iterNum){

	FloatUtils downer(i.copyImage());

	for(int s = 0; s < iterNum; ++s){
		downer.setImage(downer.downSample());
	}

	return downer.getImage();
}

/**
 * Generates mask indicating the values that are (1) out of the indicated range or (2) equal to
 * nan.
 *
 * @param i
 *   The image for which to generate a mask
 * @param mx
 *   The maximum value that should not be masked out
 * @param mn
 *   The minimun value that should not be masked out
 * @return
 *   An with with values equal to 1.0 where
 */
ByteImageType::Pointer DisparityMapper::findOutOfRange(FloatImageType::Pointer i, double mx, double mn){

	itk::Size<2> s = i->GetLargestPossibleRegion().GetSize();
	ByteUtils mask(s[0], s[1], 1);

	itk::ImageRegionIterator<ByteImageType> maskIt(mask.getImage(), i->GetLargestPossibleRegion());
	itk::ImageRegionIterator<FloatImageType> imageIt(i, i->GetLargestPossibleRegion());

	imageIt.GoToBegin();
	for(maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++maskIt){
		float f = imageIt.Get();
		if(std::isnan(f)){
			maskIt.Set(0);
		}
		else if(!(f <= mx && f >= mn)){
			maskIt.Set(0);
		}
		++imageIt;
	}

	return mask.getImage();
}

/**
 * Finds the maximum value (and max location) in the image using subpixel precision.
 * Uses bicubic interpolation to determine subpixel values.
 *
 * @param image
 *   Image for which to find subpixel maximum value/location
 * @param inc
 *   Increment used to determine subpixel precision
 * @param loc
 *   Point to store max location in upon exit, precision determined by inc
 * @return
 *   The maximum subpixel value in image
 */
double DisparityMapper::getSubpixelMaximum(FloatImageType::Pointer image, double start, double stop,
		double inc, itk::Point<double, 2> &loc){

	int r, c, count;
	double cur, rowMax, colMax;
	std::vector<double> interpPoints;
	std::vector<double> maxScores;
	std::vector<double> maxQpoints;
	std::vector<double> maxScoresVert;
	std::vector<double> maxQpointsVert;
	itk::Size<2> s = image->GetLargestPossibleRegion().GetSize();

	//Find interpolation query points
	for(cur = start; !(cur > stop); cur += inc){
		interpPoints.push_back(cur);
	}
	interpPoints.push_back(cur);

	//Create iterator and interpolator
	itk::BSplineInterpolateImageFunction<FloatImageType, double, double>::Pointer interp =
			itk::BSplineInterpolateImageFunction<FloatImageType, double, double>::New();
	interp->SetSplineOrder(3);
	interp->SetInputImage(image);

	//Within row interpolation
	for(r = 0; r < s[0]; ++r){

		double max = -DBL_MAX;
		double maxQ = -DBL_MAX;

		for(int q = 0; q < interpPoints.size(); ++q){

			itk::Point<double, 2> p;
			p[0] = interpPoints[q];
			p[1] = r;

			double curScore = interp->Evaluate(p);

			if(curScore > max){
				max = curScore;
				maxQ = q;
			}

		}
		maxScores.push_back(max);
		maxQpoints.push_back(maxQ);
	}

	rowMax = DisparityMapper::max(maxScores, &r);

	//Within col interpolation
	for(c = 0; c < s[1]; ++c){

		double max = -DBL_MAX;
		double maxQ = -DBL_MAX;

		for(int q = 0; q < interpPoints.size(); ++q){

			itk::Point<double, 2> p;
			p[0] = c;
			p[1] = interpPoints[q];

			double curScore = interp->Evaluate(p);

			if(curScore > max){
				max = curScore;
				maxQ = q;
			}

		}
		maxScoresVert.push_back(max);
		maxQpointsVert.push_back(maxQ);
	}
	colMax = DisparityMapper::max(maxScoresVert, &c);

	//Find the final max score
	if(rowMax >= colMax){
		loc[0] = r;
		loc[1] = interpPoints[(int)maxQpoints[r]];
		cur = rowMax;
	}
	else{
		loc[0] = interpPoints[(int)maxQpointsVert[c]];
		loc[1] = c;
		cur = colMax;
	}

	//Need to fix this...
	double temp;
	temp = loc[0]; loc[0] = loc[1]; loc[1] = temp;

	return cur;
}

/**
 * Replaces values in disp and score based on the mask image. Any pixels with zero mask
 * values are replaced in disp and score. The replacement values are the nearest pixels with
 * non-zero mask values.
 *
 * @param disp
 *   First image in which to replace values, altered on exit
 * @param score
 *   Second image in which to replace values, altered on exit
 * @param mask
 *   Mask image indicating which pixel values to replace
 */
void DisparityMapper::replaceWithNearest(FloatImageType::Pointer disp, FloatImageType::Pointer score,
			ByteImageType::Pointer mask){

	//Find distance transform and mapping in mask
	typedef itk::Image<itk::Offset<2>, 2> OffsetImageType;
	itk::DanielssonDistanceMapImageFilter<ByteImageType, FloatImageType>::Pointer dist =
			itk::DanielssonDistanceMapImageFilter<ByteImageType, FloatImageType>::New();
	dist->SetInput(mask);
	dist->Update();
	OffsetImageType::Pointer offset = dist->GetVectorDistanceMap();

	//Replace values
	itk::ImageRegionConstIterator<ByteImageType> maskIt(mask, mask->GetLargestPossibleRegion());
	itk::ImageRegionIterator<FloatImageType> dispIt(disp, disp->GetLargestPossibleRegion());
	itk::ImageRegionIterator<FloatImageType> scoreIt(score, score->GetLargestPossibleRegion());
	itk::ImageRegionConstIterator<OffsetImageType> offsetIt(offset, offset->GetLargestPossibleRegion());
	dispIt.GoToBegin();
	scoreIt.GoToBegin();
	offsetIt.GoToBegin();
	for(maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++maskIt){
		if(maskIt.Get() == 0){
			itk::Index<2> idx = maskIt.GetIndex();
			OffsetImageType::PixelType off = offsetIt.Get();

			double r = disp->GetPixel(idx + off);

			dispIt.Set(disp->GetPixel(idx + off));
			scoreIt.Set(score->GetPixel(idx + off));
		}

		++dispIt;
		++scoreIt;
		++offsetIt;
	}
}

/**
 * Computes image that can be used to perform additive intensity correction. This correction adjusts
 * intensities in the reference image to move the intensities closer to the estimated matching
 * region in the target image.
 *
 * @param disp
 *   The current horizontal disparity estimate
 * @param dispVert
 *   The current vertical disparity estimate
 * @param left
 *   Reference image for which to compute the correction
 * @param right
 *   Target image from which the correction is computed
 * @param nsize
 *   Size of averaging filter used to smooth the correction image
 * @param padW
 *   Amount of mirror padding that has been applied to the width left and right images
 * @param padH
 *   Amount of mirror padding that has been applied to the height left and right images
 * @return
 *   An image containing intensity correction factors that can be added to the reference image
 */
FloatImageType::Pointer DisparityMapper::correctIntensity(FloatUtils &disp, FloatUtils &dispVert,
			FloatUtils &left, FloatUtils &right, int nsize, int padW, int padH){

	itk::Index<2> idx;
	idx[0] = padW; idx[1] = padH;
	itk::Size<2> size = disp.getImage()->GetLargestPossibleRegion().GetSize();
	FloatImageType::RegionType region(idx, size);

	/************************************************************************************************/
//	disp.printImage();
//	dispVert.printImage();
//	left.printImage();
//	right.printImage();
	/************************************************************************************************/

	FloatUtils correction(size[0], size[1], 0);

	itk::BSplineInterpolateImageFunction<FloatImageType>::Pointer interp =
			itk::BSplineInterpolateImageFunction<FloatImageType>::New();
	interp->SetSplineOrder(3);
	interp->SetInputImage(right.getImage());

	itk::ImageRegionIterator<FloatImageType> dispIt(disp.getImage(),
			disp.getImage()->GetLargestPossibleRegion());
	itk::ImageRegionIterator<FloatImageType> dispItVert(dispVert.getImage(),
			dispVert.getImage()->GetLargestPossibleRegion());
	itk::ImageRegionIterator<FloatImageType> correctIt(correction.getImage(),
			correction.getImage()->GetLargestPossibleRegion());
	itk::ImageRegionIterator<FloatImageType> leftIt(left.getImage(), region);
	itk::ImageRegionIterator<FloatImageType> rightIt(right.getImage(), region);

	itk::Index<2> curIdx;
	itk::Point<double, 2> qp;

	correctIt.GoToBegin();
	dispItVert.GoToBegin();
	leftIt.GoToBegin();
	rightIt.GoToBegin();
	for(dispIt.GoToBegin(); !dispIt.IsAtEnd(); ++dispIt){

		curIdx = rightIt.GetIndex();
		qp[0] = curIdx[0] + dispIt.Get();
		qp[1] = round(curIdx[1] + dispItVert.Get());
		double corrected = interp->Evaluate(qp);
		correctIt.Set(corrected - leftIt.Get());

		++correctIt;
		++dispItVert;
		++leftIt;
		++rightIt;
	}

	/********************************************************************************************************************/
//	correction.printImage();
	/********************************************************************************************************************/

	//Average/blur the correction
	int pad = round(nsize/2.0);
	correction.setImage(correction.mirrorPadImage(pad, pad, pad, pad));
	correction.setImage(correction.average(nsize, nsize));
	correction.setImage(correction.crop(pad, pad, size[0], size[1]));

	return correction.getImage();
}

/**
 * Computes the gradient magnitude and direction of the given image.
 *
 * @param i
 *   The image for which to compute gradient information
 * @return
 *   Vector of float valued images, the first is gradient mag. and the second is direction
 *   (angle in radians from positive x-axis)
 */
std::vector<FloatImageType::Pointer> DisparityMapper::gradient(FloatImageType::Pointer i){

	FloatUtils gradient(i);
	FloatImageType::Pointer xgrad = gradient.simpleXGradient();
	FloatImageType::Pointer ygrad = gradient.simpleYGradient();

	std::vector<FloatImageType::Pointer> grads;

	itk::BinaryMagnitudeImageFilter<FloatImageType, FloatImageType, FloatImageType>::Pointer add =
			itk::BinaryMagnitudeImageFilter<FloatImageType, FloatImageType, FloatImageType>::New();
	add->SetInput1(xgrad);
	add->SetInput2(ygrad);
	add->Update();
	FloatImageType::Pointer mag = add->GetOutput();

	grads.push_back(mag);

//	itk::DivideImageFilter<FloatImageType, FloatImageType, FloatImageType>::Pointer div =
//			itk::DivideImageFilter<FloatImageType, FloatImageType, FloatImageType>::New();
//	div->SetInput1(ygrad);
//	div->SetInput2(xgrad);
//	div->Update();
//
//	itk::AtanImageFilter<FloatImageType, FloatImageType>::Pointer tan =
//			itk::AtanImageFilter<FloatImageType, FloatImageType>::New();
//	tan->SetInput(div->GetOutput());
//	tan->Update();
//	FloatImageType::Pointer dir = tan->GetOutput();

	itk::Atan2ImageFilter<FloatImageType, FloatImageType, FloatImageType>::Pointer tan =
			itk::Atan2ImageFilter<FloatImageType, FloatImageType, FloatImageType>::New();
	tan->SetInput1(ygrad);
	tan->SetInput2(xgrad);
	tan->Update();
	FloatImageType::Pointer dir = tan->GetOutput();

	grads.push_back(dir);

	return grads;
}

/**
 * Computes the correlation coefficient for two n-dimensional vectors.
 *
 * @param v1
 *   The first data vector
 * @param v2
 *   The second vector, should be same length as the first
 * @return
 *   Pearson sample correlation coefficient for the vectors
 */
double DisparityMapper::correlation(std::vector<double> &v1, std::vector<double> &v2){

	int n = v1.size();
	double res12 = 0.0, res1sq = 0.0, res2sq = 0.0;
	double avg1 = DisparityMapper::avg(v1);
	double avg2 = DisparityMapper::avg(v2);

	for(int i = 0; i < n; ++i){
		res12 += (v1[i] - avg1)*(v2[i] - avg2);
		res1sq += (v1[i] - avg1)*(v1[i] - avg1);
		res2sq += (v2[i] - avg2)*(v2[i] - avg2);
	}

	res1sq = sqrt(res1sq);
	res2sq = sqrt(res2sq);

	return (res12)/(res1sq*res2sq);

}

/**
 * Computes the average of a list of values.
 *
 * @param v1
 *   Vector of data
 * @return
 *   Average value of the data
 */
double DisparityMapper::avg(std::vector<double> &v){

	int n = v.size();
	double total = 0.0;

	for(int i = 0; i < n; ++i){
		total += v[i];
	}

	return total/n;
}

/**
 * Shifts each element in the vector by adding the value s.
 *
 * @param v
 *   Data to shift, altered on exit
 * @param s
 *   Value to add to each element of v
 */
void DisparityMapper::shift(std::vector<double> &v, double s){

	int n = v.size();

	for(int i = 0; i < n; ++i){
		v[i] = v[i] + s;
	}
}

/**
 * Takes the absolute value of each element in the vector.
 *
 * @param v
 *   Data to which abs is applied, altered on exit
 */
void DisparityMapper::abs(std::vector<double> &v){

	int n = v.size();

	for(int i = 0; i < n; ++i){
		v[i] = fabs(v[i]);
	}
}

/**
 * Finds value and location of maximum value in the vector.
 *
 * @param v
 *   The vector for which to find the maximum
 * @param loc
 *   On exit, the index of the first occurence of the maximum
 * @return
 *   Value of the max.
 */
double DisparityMapper::max(std::vector<double> &v, int *loc){

	int m = 0;
	double mx = -DBL_MAX;

	for(int i = 0; i < v.size(); ++i){
		if(v[i] > mx){
			mx = v[i];
			*loc = i;
		}
	}

	return mx;
}

/**
 * Applies exponential function to each element in v.
 *
 * An optional scaling factor can be applied to each element before the exponential
 * function is applied.
 *
 * @param v
 *   Data on which to apply exponential function, altered on exit
 * @param s
 *   Scaling factor to apply to each element of v
 */
void DisparityMapper::exp(std::vector<double> &v, double s){

	int n = v.size();

	for(int i = 0; i < n; ++i){
		v[i] = std::exp(s*v[i]);
	}
}

/**
 * Vectorizes pixel data within the current neighborhood of the iterator.
 *
 * @param nhood
 *   Neighborhood iterator of a FloatImageType image
 * @param dest
 *   Destination vector to store pixel data.
 */
void DisparityMapper::getPixelValues(itk::NeighborhoodIterator<FloatImageType> &nhood, std::vector<double> &dest){

	itk::Size<2> s = nhood.GetSize();
	int n = s[0]*s[1];
	dest.resize(n);

	for(int i = 0; i < n; ++i){
		dest[i] = nhood.GetPixel(i);
	}
}

/**
 * Inserts all data in vector from in the vector into starting at index idx.
 *
 * Performs no bounds checking.
 *
 */
void DisparityMapper::insertData(std::vector<double> &into, std::vector<double> &from, int idx){

	for(int i = 0; i < from.size(); ++i){
		into[i + idx] = from[i];
	}
}

/**
 * Creates an image of the current neighborhood. Used for debugging.
 */
FloatImageType::Pointer DisparityMapper::getNhoodImage(itk::NeighborhoodIterator<FloatImageType> &nhood){

	itk::Size<2> s = nhood.GetSize();
	int n = s[0]*s[1];

	FloatImageType::Pointer fi = this->createFilledImage(s[0], s[1], 0.0);
	itk::ImageRegionIterator<FloatImageType> it(fi, fi->GetLargestPossibleRegion());

	for(int i = 0; i < n; ++i){
		it.Set(nhood.GetPixel(i));
		++it;
	}

	return fi;
}

/**
 * Prints double vector to stdout. Used for debugging.
 */
void DisparityMapper::printVector(std::vector<double> v){

	for(int i = 0; i < v.size(); ++i){
		std::cout << v[i] << ", ";
	}

	std::cout << std::endl;
}

/**
 * Fills image with data from vector. Used for debugging.
 */
void DisparityMapper::fillImage(FloatImageType::Pointer image, double data[]){

	int i = -1;
	itk::ImageRegionIterator<FloatImageType> it(image, image->GetLargestPossibleRegion());

	for(it.GoToBegin(); !it.IsAtEnd(); ++it){
		it.Set((float)data[++i]);
	}

}

DisparityMapper::~DisparityMapper() {
	// TODO Auto-generated destructor stub
}

/*
 * main.cpp
 *
 *  Created on: Oct 24, 2012
 *      Author: mchristopher
 */

#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include "RGBUtils.h"
#include "DisparityMapper.h"

/**
 * Starting point for stereo matching program. Takes two input images as a stereo pair and outputs a float-valued depth map image
 * based on disparity computations. Typical usage:
 *
 * %stereo <left image> <right image> <output image> [<intermediate result dir>]
 *
 * The arguments <left image>, <right image>, and <output image> are required input/output paths. The <intermediate result dir> is an
 * optional argument. If provided, intermediate results at each stage of the multi-scale computation will be output to the given
 * directory.
 *
 */
int main(int argc, char **argv){

	if(argc < 4){
		std::cout << "Error: missing arguments!" << std::endl;
		std::cout << "Usage:" << std::endl;
		std::cout << "%" << argv[0] << " <left image> <right image> <output image> [<intermediate result dir>]" << std::endl;
		exit(-1);
	}

	std::string inL = argv[1];
	std::string inR = argv[2];
	std::string outPath = argv[3];

	RGBUtils rgbL(inL);
	RGBUtils rgbR(inR);

	DisparityMapper mapper(rgbL.getImage(), rgbR.getImage());

	bool saveIntermediate = argc > 4;

	if(saveIntermediate){
		struct stat sb;
		if ((stat(argv[4], &sb) == 0 && S_ISDIR(sb.st_mode)) || mkdir(argv[4] , 0777) == 0 ){
			std::string inter = argv[4];
			mapper.setIntermediateDir(inter);
			mapper.setShouldSaveIntermediateResults(true);
		}
	}

	long t0 = time(NULL);

//	FloatImageType::Pointer x, y;
	FloatUtils final(mapper.computeDisparity());
	final.writeToPath(outPath);

	long t1 = time(NULL);

	std::cout << "Finished stereo in " << (t1 - t0) << " seconds." << std::endl;

	return 0;
}

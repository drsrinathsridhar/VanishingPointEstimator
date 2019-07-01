/*
 * sample_usage.cpp
 * This file is part of VanPoints
 *
 * Copyright (C) 2011 - Srinath Sridhar and Yu Xiang
 *
 * VanPoints is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * VanPoints is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with VanPoints; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */


#include <cstdio>
#include "VanPoints.h"

double startt, frmTime;
#define TIMERSTART startt = (double) cvGetTickCount();
#define TIMERSTOP frmTime = ((double)cvGetTickCount() - startt)/(1000.0 * cvGetTickFrequency()), fprintf(stderr, "Time elapsed is %f ms\n\n", frmTime);

int main(int argc, char ** argv)
{
    if(argc != 5)
    {
        fprintf(stderr, "Usage: %s <input_image> <min_length> <jlinkage_model_size> <output_image>\n", argv[0]);
        return -1;
    }

    IplImage * inImage;
    inImage = cvLoadImage(argv[1]);

    VanPoints * VPExtractor = new VanPoints;
    VanPoints::m_VanPointStruct vanpoints;

    TIMERSTART
    vanpoints = VPExtractor->findVanishingPoints(inImage, atoi(argv[2]), atoi(argv[3]));
    TIMERSTOP
    if(VPExtractor->superposeLines(inImage, vanpoints.lines, vanpoints.lineNum, vanpoints.clusters) != 0)
        fprintf(stderr, "Error superposing lines on image. Please check you input.\n");

    cvSaveImage(argv[4], inImage);

    for(int i = 0; i < 3; ++i) // Considering only 3 orthogonal vanishing points
    {
        printf("Vanishing point %d is \n", i);
        printf("%f %f\n", vanpoints.vanPoints[i]->x / vanpoints.vanPoints[i]->z, vanpoints.vanPoints[i]->y / vanpoints.vanPoints[i]->z);
    }
    printf("Vanishing points computed successfully and line membership image written to file.\n");

    return 0;
}

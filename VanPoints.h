// Version 3.0

/*
 * VanPoints.cpp
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


// Version 2.0

// Version 1.0
// This is an implementation of
// Non-Iterative Approach for Fast and Accurate Vanishing m_Point Detection
// by Jean-Philippe Tardif
// for EECS 442, Fall 2010 class project at the
// University of Michigan, Ann Arbor
// Authors: Yu Xiang and Srinath Sridhar

// Last modified 3-June-2011

#ifndef _VANPOINTS_H_
#define _VANPOINTS_H_

#include <string>
#include <cv.h>
#include <highgui.h>
#include <cmath>
#include <ctime>

#ifdef WIN32
#define M_PI 3.141592653589793238462643
#endif

class VanPoints
{
public:
    VanPoints();
    virtual ~VanPoints();

    // Some primitive structs
    // m_Line
    typedef struct
    {
        float x1[2];     /* end point 1 */
        float x2[2];     /* end point 2 */
        float mean[2];   /* middle point */
        float l[3];      /* line representation in homogeneous coordinates */
        float theta;     /* line angle */
        float r;         /* line length */
    } m_Line;
    // m_Point structure
    typedef struct
    {
        float x;
        float y;
        float z;
    } m_Point;

    // Output of the vanishing point detection with all possible results
    typedef struct
    {
        m_Line ** lines; // All lines
        int lineNum; // Number of lines
        int ** clusters; // Individual line clusters
        m_Point ** vanPoints; // Vanishing points
    } m_VanPointStruct;


private:
    m_Line ** extractLongLine(IplImage * img, int minLen, int * lnum);
    m_Point ** findVanPoints(IplImage * img, int modelSize, m_Line ** edges, int lineCount, int * vpNum);

    // JLinkage Stuff
    // m_Model - two lines intersecting at a point is the model
    typedef struct
    {
        m_Line * line1;
        m_Line * line2;
        m_Point * intersectionPt;
    } m_Model;

    typedef struct
    {
        m_Model ** minSet;
        int size;
    } m_MinimalSet;

    float m_clusterThresh;
    int ** m_clusters;

    m_MinimalSet * selectMinimalSets(int modelSize, m_Line ** edges, int lineCount);
    int ** makePSMatrix(m_MinimalSet * RandomMS, m_Line ** edges, int lineCount);
    int ** clusterPSMatrix(int ** PSMatrix, int lineNum, int modelNum, int * clusterNum);
    m_Point * estimateVanPoint(m_Line ** edges, int lineNum, int * cluster);

    // Some handy functions
    m_Point * findIntersection(m_Line * line1, m_Line * line2);
    float findOrthDistance(m_Line * line, m_Point * point);
    float jaccardDist(int * A, int * B, int len);

public:
    int superposeLines(IplImage * img, m_Line ** lines, int lineNum, int ** clusters);
    m_VanPointStruct findVanishingPoints(IplImage * img, int lineLength, int modelSize);
};

#endif

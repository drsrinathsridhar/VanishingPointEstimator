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

#include "VanPoints.h"

VanPoints::VanPoints()
{
    m_clusterThresh = 2; // Clustering threshold
}

VanPoints::~VanPoints()
{

}

VanPoints::m_VanPointStruct VanPoints::findVanishingPoints(IplImage * img, int lineLength, int modelSize)
{
    m_Line ** edges = NULL;
    int lineNum = 0;
    edges = extractLongLine(img, lineLength, &lineNum);

    if(edges == NULL)
        fprintf(stderr, "Error. No edges detected.\n");

    // Tardif's vanishing point detection algorithm starts here
    m_Point ** vps;
    int vpNum;

    vps = findVanPoints(img, modelSize, edges, lineNum, &vpNum);

    m_VanPointStruct outStruct;
    outStruct.vanPoints = vps;
    outStruct.lines = edges;
    outStruct.lineNum = lineNum;
    outStruct.clusters = m_clusters;

    return outStruct;
}

VanPoints::m_Line ** VanPoints::extractLongLine(IplImage * img, int minLen, int * lnum)
{
    if(img == NULL)
    {
        fprintf(stderr, "Error. No input image for long line extraction.\n");
        return NULL;
    }

    int i, j, k;
    m_Line ** lines = NULL;
    int lineNum = 0;

    int width = img->width;
    int height = img->height;
    IplImage * grayImgF;
    IplImage * grayImgU;

    if(img->depth == IPL_DEPTH_32F)
    {
        grayImgF = cvCreateImage(cvSize(width, height), IPL_DEPTH_32F, 1);
        cvCvtColor(img, grayImgF, CV_BGR2GRAY);
    }
    else if(img->depth == IPL_DEPTH_8U)
    {
        grayImgU = cvCreateImage(cvSize(width, height), IPL_DEPTH_8U, 1);
        cvCvtColor(img, grayImgU, CV_BGR2GRAY);
        grayImgF = cvCreateImage(cvSize(width, height), IPL_DEPTH_32F, 1);
        cvConvertScale(grayImgU, grayImgF);
    }
    else
    {
        fprintf(stderr, "VanPoints: Unsupported input image depth. Aborting vanishing point estimation.\n");
        return NULL;
    }

    CvMat * temp = cvCreateMat(height, width, CV_32FC1);
    cvSmooth(grayImgF, temp, CV_GAUSSIAN, 7, 7, 1.5); // Smooth the image with Gaussian kernel
    CvMat * dx = cvCreateMat(height, width, CV_32FC1);
    CvMat * dy = cvCreateMat(height, width, CV_32FC1);
    cvSobel(temp, dx, 1, 0); // Gradient of X
    cvSobel(temp, dy, 0, 1); // Gradient of Y

    CvMat * dst = cvCreateMat(height, width, CV_8UC1);
    cvCanny(grayImgU, dst, 80, 100); // Canny edge extraction

    // Remove boundary lines
    for(i = 0; i < height; ++i)
    {
        dst->data.ptr[i*width] = 0;
        dst->data.ptr[i*width + 1] = 0;
        dst->data.ptr[i*width + width-2] = 0;
        dst->data.ptr[i*width + width-1] = 0;
    }
    for(j = 0; j < width; ++j)
    {
        dst->data.ptr[j] = 0;
        dst->data.ptr[width + j] = 0;
        dst->data.ptr[(height-2)*width + j] = 0;
        dst->data.ptr[(height-1)*width + j] = 0;
    }

    int binNum = 8; // Bin size for edge pixel orientation
    for(i = 0; i < height; ++i)
    {
        for(j = 0; j < width; ++j)
        {
            if(*((unsigned char *)CV_MAT_ELEM_PTR(*dst, i, j)) == 255) // if is edge pixel
            {
                float anger = atan(cvmGet(dy,i,j) / (cvmGet(dx,i,j) + 1e-10)); // gradient orientation
                int index = (int) ((anger+M_PI/2) *binNum / M_PI + 1); // bin index
                cvmSet(temp, i, j, (int) index); // store the bin index to temp
            }
            else
                cvmSet(temp, i, j, 0); // set the bin index for non-edge pixel to zero
        }
    }


    CvMat * used = cvCreateMat(height, width, CV_8UC1);
    CvMat * dirImage = cvCreateMat(height, width, CV_8UC1);
    CvMat * D = cvCreateMat(2, 2, CV_32FC1);
    CvMat * evects = cvCreateMat(2, 2, CV_32FC1);
    CvMat * evals = cvCreateMat(2, 1, CV_32FC1);
    CvMemStorage * mem = cvCreateMemStorage(0);
    CvSeq * contours, * ptr;

    cvZero(dst);
    cvZero(used);
    for(int bin = 1; bin <= binNum; ++bin)
    {
        cvZero(dirImage); // set the direction image to zero
        for(i = 0; i < height; ++i)
        {
            for(j = 0; j < width; ++j)
            {
                int index = (int) cvmGet(temp, i, j); // bin index
                if(*((uchar*)CV_MAT_ELEM_PTR(*used, i, j)) == 0) // if it is not used
                {
                    // cross bin usage
                    if(index == bin || index == (bin == 1 ? binNum : bin-1) || index == (bin == binNum ? 1 : bin+1))
                        *((unsigned char *)CV_MAT_ELEM_PTR(*dirImage, i, j)) = 1;
                }
            }
        }
        cvFindContours(dirImage, mem, &contours); // find contours
        if(contours)
        {
            for(ptr = contours; ptr; ptr = ptr->h_next) // consider each line support region
            {
                if(ptr->total > minLen) // if the region size is larger than threshold
                {
                    float mean_x = 0, mean_y = 0;
                    float max_x = 0, max_y = 0;
                    float min_x = width, min_y = height;
                    for(k = 0; k < ptr->total; k++)
                    {
                        CvPoint *p = (CvPoint*)cvGetSeqElem(ptr, k);
                        mean_x += p->x;
                        mean_y += p->y;
                        if(p->x > max_x)
                            max_x = p->x;
                        if(p->y > max_y)
                            max_y = p->y;
                        if(p->x < min_x)
                            min_x = p->x;
                        if(p->y < min_y)
                            min_y = p->y;
                    }
                    mean_x /= ptr->total;
                    mean_y /= ptr->total;
                    cvZero(D);
                    for(k = 0; k < ptr->total; k++)
                    {
                        CvPoint *p = (CvPoint*)cvGetSeqElem(ptr, k);
                        cvmSet(D, 0, 0, cvmGet(D,0,0)+(p->x - mean_x)*(p->x - mean_x));
                        cvmSet(D, 0, 1, cvmGet(D,0,1)+(p->x - mean_x)*(p->y - mean_y));
                        cvmSet(D, 1, 1, cvmGet(D,1,1)+(p->y - mean_y)*(p->y - mean_y));
                    }
                    cvmSet(D, 1, 0, cvmGet(D, 0, 1));
                    cvEigenVV(D, evects, evals);
                    float theta = atan2(cvmGet(evects, 0, 1), cvmGet(evects, 0, 0));
                    float conf;
                    if(cvmGet(evals, 1, 0) > 0)
                        conf = cvmGet(evals, 0, 0) / cvmGet(evals, 1, 0);
                    else
                        conf = 100000;
                    if(conf >= 400)
                    {
                        float r = sqrt((max_x - min_x)*(max_x - min_x) + (max_y - min_y)*(max_y - min_y));
                        float x1 = mean_x - cos(theta) * r / 2;
                        float x2 = mean_x + cos(theta) * r / 2;
                        float y1 = mean_y - sin(theta) * r / 2;
                        float y2 = mean_y + sin(theta) * r / 2;
                        for(k = 0; k < ptr->total; k++)
                        {
                            CvPoint *p = (CvPoint*)cvGetSeqElem(ptr, k);
                            *((uchar*)CV_MAT_ELEM_PTR(*used, p->y, p->x)) = 1;
                        }

                        if((lines = (m_Line**)realloc(lines, sizeof(m_Line*)*(lineNum+1))) == NULL)
                        {
                            fprintf(stderr, "Realloc error in Graph::addLine.\n");
                            return NULL;
                        }
                        lines[lineNum] = new m_Line();
                        lines[lineNum]->x1[0] = x1;
                        lines[lineNum]->x1[1] = y1;
                        lines[lineNum]->x2[0] = x2;
                        lines[lineNum]->x2[1] = y2;
                        lines[lineNum]->mean[0] = mean_x;
                        lines[lineNum]->mean[1] = mean_y;
                        lines[lineNum]->theta = theta;
                        lines[lineNum]->r = r;
                        lines[lineNum]->l[0] = y1 - y2;
                        lines[lineNum]->l[1] = x2 - x1;
                        lines[lineNum]->l[2] = x1*y2 - x2*y1;
                        lineNum++;
                    }
                } // end if (ptr->toal > minLen)
            } // end for ptr
        } // end if contours
    }

    cvReleaseImage(&grayImgU);
    cvReleaseImage(&grayImgF);
    cvReleaseMat(&dst);
    cvReleaseMat(&temp);
    cvReleaseMat(&dx);
    cvReleaseMat(&dy);
    cvReleaseMat(&used);
    cvReleaseMat(&dirImage);
    cvReleaseMat(&D);
    cvReleaseMat(&evals);
    cvReleaseMat(&evects);
    cvReleaseMemStorage(&mem);

    *lnum = lineNum;
    return lines;
}

VanPoints::m_Point ** VanPoints::findVanPoints(IplImage * img, int modelSize, m_Line ** edges, int lineCount, int * vpNum)
{
#ifdef DEBUG
    printf("Number of edges is %d\n", lineCount);
#endif

    m_MinimalSet * RandomMS;
    // Select minimal sets
    RandomMS = selectMinimalSets(modelSize, edges, lineCount);
    // Construct Preference Set Matrix
    int** PSMatrix;
    PSMatrix = makePSMatrix(RandomMS, edges, lineCount);

    // Next 2 lines just to save the unclustered PSMatrix to an image
    CvMat * scaled = cvCreateMat(lineCount, RandomMS->size, CV_8UC1);
    cvSetZero(scaled);
    for(int i = 0; i < lineCount; ++i)
    {
        for(int j = 0; j < RandomMS->size; ++j)
        {
            if(PSMatrix[i][j] == true)
                *((uchar*)CV_MAT_ELEM_PTR(*scaled, i, j)) = 255;
            else
                *((uchar*)CV_MAT_ELEM_PTR(*scaled, i, j)) = 0;
        }
    }
#ifdef DEBUG
    cvSaveImage("PSMatrix_raw.png", scaled);
#endif
    cvReleaseMat(&scaled);

    //Perform clustering on PSMatrix
    int clusterNum;
    m_clusters = clusterPSMatrix(PSMatrix, lineCount, RandomMS->size, &clusterNum);

#ifdef DEBUG
    superposeLines(img, edges, lineCount, m_clusters);
    fprintf(stderr, "Superimposing lines on input image.\n");
#endif

    m_Point** output = new m_Point*[clusterNum];
    for(int i = 0; i < clusterNum; ++i)
        output[i] = estimateVanPoint(edges, lineCount, m_clusters[i]);

    /*  for(int i = 0; i < clusterNum; ++i)
         free(m_clusters[i]);
      free(m_clusters);*/

    for(int i = 0; i < lineCount; ++i)
        free(PSMatrix[i]);
    free(PSMatrix);

    for(int i = 0; i < RandomMS->size; ++i)
    {
        free(RandomMS->minSet[i]->intersectionPt);
        free(RandomMS->minSet[i]);
    }

    *vpNum = clusterNum;
    return output;
}

VanPoints::m_MinimalSet * VanPoints::selectMinimalSets(int modelSize, m_Line ** edges, int lineCount)
{
    srand((unsigned) time(NULL)); // Seed the rand function
    int modelsM = modelSize;
    m_MinimalSet * countedModels;
    countedModels = new m_MinimalSet;
    m_Model ** allModels;
    allModels = new m_Model * [modelsM];

    for(int i = 0; i < modelsM; ++i)
    {
        allModels[i] = new m_Model;
        // For each modelsM choose, randomly, two lines
        int n1 = rand()%lineCount;
        int n2 = rand()%lineCount;
        while(n2 == n1)
            n2 = rand()%lineCount;
        allModels[i]->line1 = edges[n1];
        allModels[i]->line2 = edges[n2];

        // Find vanishing point (intersection) for chosen model
        allModels[i]->intersectionPt = findIntersection(allModels[i]->line1, allModels[i]->line2);
    }
    countedModels->minSet = allModels;
    countedModels->size = modelsM;
    return countedModels;
}

int ** VanPoints::makePSMatrix(m_MinimalSet * RandomMS, m_Line ** edges, int lineCount)
{
    int** PSMatrix = new int*[lineCount];
    for (int i = 0; i < lineCount; ++i)
        PSMatrix[i] = new int[RandomMS->size];
    float distance;

    for(int i = 0; i < lineCount; ++i)
    {
        for(int j = 0; j < RandomMS->size; ++j)
        {
            // See proximity of ith line with jth model
            // Proximity in this case is perpendicular distance
            distance = findOrthDistance(edges[i], RandomMS->minSet[j]->intersectionPt);
            if(distance <= m_clusterThresh)
                PSMatrix[i][j] = 1;
            else
                PSMatrix[i][j] = 0;
        }
    }
    return PSMatrix;
}
int ** VanPoints::clusterPSMatrix(int ** PSMatrix, int lineNum, int modelNum, int * clusterNum)
{
    /* allocate memory */
    float** distances = new float*[lineNum];
    m_clusters = new int*[lineNum];
    for(int i = 0; i < lineNum; ++i)
    {
        distances[i] = new float[lineNum];
        memset(distances[i], 0, sizeof(float)*lineNum);
        m_clusters[i] = new int[lineNum];
        memset(m_clusters[i], 0, sizeof(int)*lineNum);
    }
    int* indicators = new int[lineNum];

    /* initialization */
    for(int i = 0; i < lineNum; ++i)
        m_clusters[i][i] = 1;

    float minDis = 1;
    int indexA = 0, indexB = 0;
    for(int i = 0; i < lineNum; ++i)
    {
        indicators[i] = 1;
        for(int j = i + 1; j < lineNum; ++j)
        {
            distances[i][j] = jaccardDist(PSMatrix[i], PSMatrix[j], modelNum);
            distances[j][i] = distances[i][j];
            if(distances[i][j] < minDis)
            {
                minDis = distances[i][j];
                indexA = i;
                indexB = j;
            }
        }
    }

    while(minDis != 1)
    {
        /* merge two m_clusters */
        for(int i = 0; i < lineNum; ++i)
        {
            if(m_clusters[indexA][i] == 1 || m_clusters[indexB][i] == 1)
                m_clusters[indexA][i] = m_clusters[indexB][i] = 1;
        }
        indicators[indexB] = 0;
        for(int i = 0; i < modelNum; ++i)
        {
            if(PSMatrix[indexA][i] == 1 && PSMatrix[indexB][i] == 1)
                PSMatrix[indexA][i] = PSMatrix[indexB][i] = 1;
            else
                PSMatrix[indexA][i] = PSMatrix[indexB][i] = 0;
        }

        /* recalculate distance */
        for(int i = 0; i < lineNum; ++i)
        {
            distances[indexA][i] = jaccardDist(PSMatrix[indexA], PSMatrix[i], modelNum);
            distances[i][indexA] = distances[indexA][i];
        }

        /* find minimum distance */
        minDis = 1;
        for(int i = 0; i < lineNum; ++i)
        {
            if(indicators[i] == 0) continue;
            for(int j = i + 1; j < lineNum; ++j)
            {
                if(indicators[j] == 0) continue;
                if(distances[i][j] < minDis)
                {
                    minDis = distances[i][j];
                    indexA = i;
                    indexB = j;
                }
            }
        }
    }

    /* calculate cluster size */
    int* clusterSizes = new int[lineNum];
    for(int i = 0; i < lineNum; ++i)
    {
        clusterSizes[i] = 0;
        if(indicators[i])
        {
            for(int j = 0; j < lineNum; ++j)
            {
                if(m_clusters[i][j])
                    clusterSizes[i]++;
            }
        }
    }

    *clusterNum = 3; 			/* choose the largest three m_clusters */
    int** result = new int*[*clusterNum];
    for(int i = 0; i < *clusterNum; ++i)
        result[i] = new int[lineNum];

    int count = 0;
    while(count < *clusterNum)
    {
        int max_index = 0;
        int max_size = clusterSizes[0];
        for(int i = 1; i < lineNum; ++i)
        {
            if(max_size < clusterSizes[i])
            {
                max_size = clusterSizes[i];
                max_index = i;
            }
        }
        for(int i = 0; i < lineNum; ++i)
            result[count][i] = m_clusters[max_index][i];
        count++;
        clusterSizes[max_index] = 0;
    }

#ifdef DEBUG
    /* print m_clusters */
    for(int i = 0; i < *clusterNum; ++i)
    {
        printf("Cluster %d:\n", i);
        for(int j = 0; j < lineNum; ++j)
        {
            if(result[i][j])
                printf("%d ", j);
        }
        printf("\n");
    }
#endif

    /* free memory */
    for(int i = 0; i < lineNum; ++i)
    {
        free(m_clusters[i]);
        free(distances[i]);
    }
    free(m_clusters);
    free(distances);
    free(indicators);

    return result;
}

VanPoints::m_Point * VanPoints::estimateVanPoint(m_Line ** edges, int lineNum, int * cluster)
{
    int i;
    int num = 0;

    for(i = 0; i < lineNum; ++i)
    {
        if(cluster[i])
            num++;
    }

    CvMat* A = cvCreateMat(num, 3, CV_32FC1);
    int count = 0;
    for(i = 0; i < lineNum; ++i)
    {
        if(cluster[i])
        {
            float l0 = edges[i]->l[0];
            float l1 = edges[i]->l[1];
            float l2 = edges[i]->l[2];
            float nrm = sqrt(l0*l0 + l1*l1 + l2*l2);
            cvmSet(A, count, 0, l0 / nrm);
            cvmSet(A, count, 1, l1 / nrm);
            cvmSet(A, count, 2, l2 / nrm);
            count++;
        }
    }

    CvMat* U = cvCreateMat(num, num, CV_32FC1);
    CvMat* D = cvCreateMat(num, 3, CV_32FC1);
    CvMat* V = cvCreateMat(3, 3, CV_32FC1);
    cvSVD(A, D, U, V);

    m_Point* vp = new m_Point();
    vp->x = cvmGet(V, 0, 2);
    vp->y = cvmGet(V, 1, 2);
    vp->z = cvmGet(V, 2, 2);

    cvReleaseMat(&A);
    cvReleaseMat(&U);
    cvReleaseMat(&D);
    cvReleaseMat(&V);

    return vp;
}

int VanPoints::superposeLines(IplImage * img, m_Line ** lines, int lineNum, int ** m_clusters)
{
    if(lines == NULL)
    {
        printf("superposeLines():: No data in the extracted lines.\n");
        return -1;
    }

    for(int i = 0; i < lineNum; ++i)
    {
        CvScalar color;
        if(m_clusters[0][i] == 1)
            color = CV_RGB(0, 0, 255);
        else if(m_clusters[1][i] == 1)
            color = CV_RGB(0, 255, 0);
        else if(m_clusters[2][i] == 1)
            color = CV_RGB(255, 0, 0);
        else continue;

        float x1 = lines[i]->x1[0];
        float y1 = lines[i]->x1[1];
        float x2 = lines[i]->x2[0];
        float y2 = lines[i]->x2[1];
        cvLine(img, cvPoint((int)x1, (int)y1), cvPoint((int)x2, (int)y2), color, 2);
    }

    return 0;
}

VanPoints::m_Point * VanPoints::findIntersection(m_Line * line1, m_Line * line2)
{
    m_Point *inter = new m_Point();

    inter->x = line1->l[1] * line2->l[2] - line1->l[2] * line2->l[1];
    inter->y = line1->l[2] * line2->l[0] - line1->l[0] * line2->l[2];
    inter->z = line1->l[0] * line2->l[1] - line1->l[1] * line2->l[0];

    return inter;
}

float VanPoints::findOrthDistance(m_Line * line, m_Point * point)
{
    float lhat[3];
    lhat[0] = line->mean[1] * point->z - point->y;
    lhat[1] = point->x - line->mean[0] * point->z;
    lhat[2] = line->mean[0] * point->y - line->mean[1] * point->x;

    float dis = abs(lhat[0]*line->x1[0] + lhat[1]*line->x1[1] + lhat[2]) / sqrt(lhat[0]*lhat[0] + lhat[1]*lhat[1]);
    return dis;
}

float VanPoints::jaccardDist(int * A, int * B, int len)
{
    int n1 = 0;
    int n2 = 0;

    for(int i = 0; i < len; ++i)
    {
        if(A[i] == 1 || B[i] == 1)
            n1++;
        if(A[i] == 1 && B[i] == 1)
            n2++;
    }

    float dis = (float)(n1 - n2) / (float)n1;
    return dis;
}


/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "data-sink.h"
#include "alloc2d.h"
#include "meter.h"

Histogram::Histogram(int n, int ht, int hs, Meter* xm, double xmin, double xmax) : xMeter(xm), nData(n), HistogramSize(hs), HistogramType(ht),  xMin(xmin), xMax(xmax) {
  unset();
  if (HistogramType < 0 || HistogramType > 3) {
    fprintf(stderr, "Unrecognized Histogram type %d\n", HistogramType);
    abort();
  }
  if ((HistogramType == 2 || HistogramType == 3) && HistogramSize % 2 == 1) {
    fprintf(stderr, "Collapsing Histogram requires even numbered Histogram size\n");
    abort();
  }
  if (nData>0) reset();
}

Histogram::~Histogram() {
  dispose();
}

void Histogram::dispose() {
  free2D((void**)data);
  free2D((void**)histogram);
  free(collapseSum);
}

void Histogram::unset() {
  data = nullptr;
  histogram = nullptr;
  collapseSum = nullptr;
}

void Histogram::setHistogramType(int t) {
  HistogramType = t;
  reset();
}

void Histogram::reset() {
  dispose();
  data = (long**)malloc2D(nData, HistogramSize, sizeof(long));
  xData = (double*)malloc(HistogramSize*sizeof(double));
  histogram = (double**)malloc2D(nData, HistogramSize, sizeof(double));
  if (HistogramType == 3) {
    collapseSum = (double*)malloc((nData)*sizeof(double));
    for (int k=0; k<nData; k++) collapseSum[k] = 0;
  }
  for (int j = 0; j<nData;j++) {
    for (int k=0; k<HistogramSize; k++) {
      data[j][k] = 0;
    }
  }

}

void Histogram::collapseDiscard() {
  // discard all odd data, keep the even (starting with 0)
  for (int i=0; i<nData; i++) {
    for (int j=0; j<HistogramSize; j+=2) {
      data[i][j/2] = data[i][j];
    }
  }
}

void Histogram::collapseAverage() {
  // average together each pair (0+1, 2+3, 4+5) of data
  for (int i=0; i<nData; i++) {
    for (int j=0; j<HistogramSize; j+=2) {
      data[i][j/2] = 0.5 * (data[i][j]+data[i][j+1]);
    }
  }
}

void Histogram::addData(double* x) {
  double deltax = (xMax - xMin)/HistogramSize;

  switch (HistogramType) {
    case 0:
      // scrolling
      for (int i=0; i<nData; i++) {
        int idx = (x[i] - xMin) / deltax;
        if (idx >= 0 && idx < HistogramSize) {
          data[i][idx]++;
        }
      }

      break;
    case 1:
      for (int i=0; i<nData; i++) {
        int idx = (x[i] - xMin) / deltax;
        if (x[i] == xMax) {
          idx = HistogramSize-1;
        }
        if (idx < 0 || idx >= HistogramSize) {
          double newxMin = min(xMin, x[i]);
          double newxMax = max(xMax, x[i]);

          int newHistogramSize = ceil((newxMax - newxMin)/deltax);
          printf("realloc %f %d %d\n", x[i], HistogramSize, newHistogramSize);

          long** copy = (long**)copy2D((unsigned char**)data, nData, HistogramSize, sizeof(long));
          data = (long**)realloc2D((void**)data, nData, newHistogramSize, sizeof(long));
          for (int k=0; k<nData; k++) {
            std::copy(copy[k], copy[k]+HistogramSize, data[k]);
          }
          xData = (double*)realloc(xData, HistogramSize*sizeof(double));

          histogram = (double**)realloc2D((void**)data, nData, newHistogramSize, sizeof(double));
          if (x[i] < xMin) {
            for (int j = 0; j<nData;j++) {
              for (int k = HistogramSize-1; k>0; k--) {
                data[j][k+newHistogramSize-HistogramSize]=data[j][k];

              }
              for (int k=0; k<newHistogramSize-HistogramSize; k++) {
                data[j][k] = 0;
              }
            }
            xMin-=(newHistogramSize - HistogramSize)*deltax;

          }
          else {
            xMax+=(newHistogramSize - HistogramSize)*deltax;
            for (int j = 0; j<nData;j++) {
              for (int k=HistogramSize; k<newHistogramSize; k++) {
                data[j][k] = 0;
              }
            }
            for (int k=0; k<newHistogramSize; k++) {
              printf("%d %d\n", k, data[1][k]);
            }

          }
          HistogramSize = newHistogramSize;

          idx = (x[i] - xMin) / deltax;

        }
        data[i][idx]++;
      }
      // complete
      break;
    case 2:
      // collapse average
      break;
    default:
      fprintf(stderr, "Unknown Histogram type %d\n", HistogramType);
      abort();
  }
  count++;
}

int Histogram::getHistogramSize() {
  return HistogramSize;
}

double** Histogram::getHistogram() {
  double deltax = (xMax - xMin)/HistogramSize;
  for (int i=0; i<nData; i++) {
    for (int j=0; j<HistogramSize; j++) {
      histogram[i][j] = data[i][j]/deltax/count;

    }
  }
  return histogram;
}

double* Histogram::getxData() {
  for (int i=0; i<HistogramSize; i++) {
    xData[i] = xMin + i*(xMax - xMin) / HistogramSize;
  }
  return xData;
}

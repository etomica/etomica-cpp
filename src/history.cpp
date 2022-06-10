/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "data-sink.h"
#include "alloc2d.h"
#include "meter.h"

History::History(int n, int ht, int hs, Meter* xm) : xMeter(xm), nData(n), historySize(hs), historyType(ht), collapseSize(1), skipCount(0), collapseSum(0) {
  unset();
  if (historyType < 0 || historyType > 3) {
    fprintf(stderr, "Unrecognized history type %d\n", historyType);
    abort();
  }
  if ((historyType == 2 || historyType == 3) && historySize % 2 == 1) {
    fprintf(stderr, "Collapsing history requires even numbered history size\n");
    abort();
  }
  if (nData>0) reset();
}

History::~History() {
  dispose();
}

void History::dispose() {
  free2D((void**)data);
  free2D((void**)history);
  free(collapseSum);
}

void History::unset() {
  data = nullptr;
  history = nullptr;
  collapseSum = nullptr;
}

void History::reset() {
  data = (double**)malloc2D(1+nData, historySize, sizeof(double));
  history = (double**)malloc2D(1+nData, historySize, sizeof(double));
  if (historyType == 3) collapseSum = (double*)malloc((1+nData)*sizeof(double));
  count = pointer = 0;
}

void History::collapseDiscard() {
  // discard all odd data, keep the even (starting with 0)
  for (int i=0; i<nData; i++) {
    for (int j=0; j<historySize; j+=2) {
      data[i][j/2] = data[i][j];
    }
  }
  pointer = historySize / 2;
}

void History::collapseAverage() {
  // average together each pair (0+1, 2+3, 4+5) of data
  for (int i=0; i<nData; i++) {
    for (int j=0; j<historySize-1; j+=2) {
      data[i][j/2] = 0.5 * (data[i][j]+data[i][j+1]);
    }
  }
  pointer = historySize / 2;
}

void History::addData(double* x) {
  switch (historyType) {
    case 0:
      // scrolling
      pointer = count % historySize;
      data[0][pointer] = xMeter->getData()[0];
      for (int i=0; i<nData; i++) data[1+i][pointer] = x[i];
      break;
    case 1:
      // complete
      if (pointer == historySize) {
        historySize *= 2;
        data = (double**)realloc2D((void**)data, 1+nData, historySize, sizeof(double));
        history = (double**)realloc2D((void**)data, 1+nData, historySize, sizeof(double));
      }
      data[0][pointer] = xMeter->getData()[0];
      for (int i=0; i<nData; i++) data[1+i][pointer] = x[i];
      pointer++;
      break;
    case 2:
      // collapse discard -- keep the first then skip until we want another
      if (skipCount == 0) {
        // first incoming data; take it
        if (pointer >= historySize) {
          // we're full, discard some data
          collapseDiscard();
        }
        data[0][pointer] = xMeter->getData()[0];
        for (int i=0; i<nData; i++) data[1+i][pointer] = x[i];
        pointer++;
      }
      skipCount++;
      if (skipCount == collapseSize) {
        // take the next data
        skipCount = 0;
      }
      break;
    case 3:
      // collapse average
      if (skipCount == 0 && pointer >= historySize) {
        collapseAverage();
      }
      // sum this data for later
      collapseSum[0] = xMeter->getData()[0];
      for (int i=0; i<nData; i++) collapseSum[1+i] += x[i];
      skipCount++;
      if (skipCount < collapseSize) {
        // we still don't have enough data to save an average
        break;
      }
      skipCount = 0;
      for (int i=0; i<nData; i++) {
        data[i][pointer] = collapseSum[i] / collapseSize;
        collapseSum[i] = 0;
      }
      pointer++;
      break;
    default:
      fprintf(stderr, "Unknown history type %d\n", historyType);
      abort();
      break;
  }
  count++;
}

int History::getHistorySize() {
  switch (historyType) {
    case 0:
      return count > historySize ? historySize : pointer;
    case 1:
    case 2:
    case 3:
      return pointer;
  }
  return -1;
}

double** History::getHistory() {
  switch (historyType) {
    case 0:
      {
        int offset = 0;
        if (count > historySize && historySize > pointer) {
          offset = historySize - pointer;
          for (int i=0; i<nData; i++) {
            std::copy(data[i]+pointer, data[i]+historySize, history[i]);
          }
        }
        for (int i=0; i<nData; i++) {
          std::copy(data[i], data[i]+pointer, history[i]+offset);
        }
      }
      break;
    case 1:
    case 2:
    case 3:
      for (int i=0; i<nData; i++) {
        std::copy(data[i], data[i]+pointer, history[i]);
      }
      break;
  }
  return history;
}

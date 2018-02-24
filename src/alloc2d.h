#include <stdlib.h>

inline static void** malloc2D(int rows, int cols, size_t s) {
  void* raw = malloc(s*rows*cols);
  void** array = (void**)malloc(rows*sizeof(void*));
  for (int i=0; i<rows; i++) {
    array[i] = ((char*)raw)+s*cols*i;
  }
  return array;
}

inline static void** malloc2D0(int rows, int totalRows, int cols, size_t s) {
  void* raw = malloc(s*rows*cols);
  void** array = (void**)malloc(totalRows*sizeof(void*));
  for (int i=0; i<rows; i++) {
    array[i] = ((char*)raw)+s*cols*i;
  }
  return array;
}

inline static void** realloc2D(void** array, int rows, int cols, size_t s) {
  if (!array) return malloc2D(rows, cols, s);
  void* newRaw = realloc(*array, s*rows*cols);
  void** newArray = (void**)realloc(array, rows*sizeof(void*));
  for (int i=0; i<rows; i++) {
    newArray[i] = ((char*)newRaw)+s*cols*i;
  }
  return newArray;
}

inline static void** realloc2D0(void** array, int rows, int totalRows, int cols, size_t s) {
  if (!array) return malloc2D0(rows, totalRows, cols, s);
  void* newRaw = realloc(*array, s*rows*cols);
  void** newArray = (void**)realloc(array, totalRows*sizeof(void*));
  for (int i=0; i<rows; i++) {
    newArray[i] = ((char*)newRaw)+s*cols*i;
  }
  return newArray;
}


inline static void free2D(void** array) {
  free(*array);
  free(array);
}

inline static void*** malloc3D(int outer, int rows, int cols, size_t s) {
  //printf("3D alloc %d %d %d %d\n", outer, rows, cols, s);
  void* raw = malloc(s*outer*rows*cols);
  //printf("raw %p total %d\n", raw, s*outer*rows*cols);
  void** raw2d = (void**)malloc(outer*rows*sizeof(void*));
  //printf("raw2d %p total %d\n", raw2d, sizeof(void*)*outer*rows);
  void*** array = (void***)malloc(outer*sizeof(void**));
  //printf("array %p total %d\n", array, sizeof(void**)*outer);
  for (int i=0; i<outer; i++) {
    array[i] = raw2d + rows*i;
    //printf("array[%d] = %p  %d %d %d\n", i, array[i], sizeof(void*), rows, i);
    for (int j=0; j<rows; j++) {
      //printf("%d %d => %d\n", i, j, s*rows*cols*i + s*cols*j);
      array[i][j] = ((char*)raw) + s*rows*cols*i + s*cols*j;
    }
  }
  return array;
}

inline static void free3D(void*** array) {
  free(**array);
  free(*array);
  free(array);
}

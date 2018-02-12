
class ArrayUtil {
  private:
    void* p;
  public:
    ArrayUtil(void* pointer) : p(pointer) {}
    ~ArrayUtil() {}
    double x(int i) { return ((double*)p)[i]; }
    double x2d(int i, int j) { return ((double**)p)[i][j]; }
};

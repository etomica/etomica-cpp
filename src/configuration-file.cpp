
#include "action.h"
#include "box.h"

ConfigurationFile::ConfigurationFile(Box& b, const char* f) : box(b), filename(f) {}

ConfigurationFile::~ConfigurationFile() { }

void ConfigurationFile::go() {
  char buf[512];

  FILE* inp;
  if (!(inp=fopen (filename, "r"))) {
    fprintf (stderr, "unable to open configuration file '%s'\n", filename);
    abort();
  }

  int numAtoms = box.getNumAtoms();
  for (int iAtom=0; iAtom<numAtoms; iAtom++) {
    fgets(buf, 510, inp);
    double *ri = box.getAtomPosition(iAtom);
    for (int j=0; j<3; j++) {
      char* c = buf;
      char* c2;
      ri[j] = strtod(c, &c2);
      if (c == c2) {
        fprintf(stderr, "Could not parse line %d from config file\n", iAtom+1);
        abort();
      }
    }
  }
  fclose(inp);
}

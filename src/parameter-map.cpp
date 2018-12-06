/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <string.h>
#include "parameter-map.h"
#include "util.h"

void ParameterMap::addParameter(string k, string defaultValue) {
  values[k] = defaultValue;
}

void ParameterMap::readFile(const char* filename) {
  FILE* f;
  if (!(f = fopen(filename, "r"))) {
    fprintf(stderr, "Unable to open species file '%s'\n", filename);
    abort();
  }
  char buf[512];
  while (fgets(buf, 510, f)) {
    if (buf[0] == '#') continue; // comment
    char* c = buf;
    char* key = strsep(&c, " ");
    if (!key || key[0] == '#') continue;
    char* value = trim(c);
    string keyS = string(key);
    string valueS = string(value);
    values[keyS] = valueS;
  }
  fclose(f);
}

void ParameterMap::parseArgs(const int argc, char** argv) {
  if (argc%2 != 0) {
    fprintf(stderr, "I need an even number of arguments (I see %d)\n", argc);
    abort();
  }
  for (int i=0; i<argc; i+=2) {
    if (argv[i][0] != '-' || strlen(argv[i]) < 2) {
      fprintf(stderr, "Invalid argument: %s\n", argv[i]);
      abort();
    }

    if (strcmp(argv[i], "-in") == 0) {
      readFile(argv[i+1]);
      continue;
    }

    string keyS = string(argv[i]+1);
    string valueS = string(argv[i+1]);
    values[keyS] = valueS;
  }
}

double ParameterMap::getDouble(string k) {
  return stod(values[k],nullptr);
}

int ParameterMap::getInt(string k) {
  return stoi(values[k],nullptr);
}

int ParameterMap::getLong(string k) {
  return stol(values[k],nullptr);
}

string ParameterMap::getString(string k) {
  return values[k];
}

bool ParameterMap::getBool(string k) {
  return values[k] == "true";
}

bool ParameterMap::hasParameter(string k) {
  return values.find(k) != values.end();
}

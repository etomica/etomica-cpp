/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#pragma once

#include <string>
#include <map>

using namespace std;

class ParameterMap {
  protected:
    map<string,string> values;

  public:
    ParameterMap() {}
    virtual ~ParameterMap() {}
    virtual void addParameter(string k, string defaultValue);
    virtual void readFile(const char* filename);
    virtual void parseArgs(const int argc, char** argv);
    virtual double getDouble(string k);
    virtual int getInt(string k);
    virtual string getString(string k);
    virtual bool hasParameter(string k);
};

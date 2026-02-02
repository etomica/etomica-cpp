//
// Created by raghu on 2/2/2026.
//

#include "diagrams.h"

#include <stdexcept>
#include <vector>


vector<vector <int>> diagrams::getFlipPointsforDiagram(const string &diagram)
{
    if (diagram == "1") {
        return {{0, 1}};
    } else if (diagram =="5c") {
        return {{0, 1}, {0, 2}};

    } else if (diagram =="54c") {
        return {{0, 3}};
    } else if (diagram =="52c") {
        return {{0, 3}, {1, 2}, {0, 1, 2}};
    } else if (diagram =="38c") {
        return {{0, 1}, {0, 3}, {0, 2}};
    }
    else if (diagram =="954c") {
        return {{0, 4}};
    } else if (diagram =="946c") {
        return {{0, 4}};

    } else if (diagram =="952c") {
        return {{0, 4}};

    } else if (diagram =="882c") {
        return {};

    } else if (diagram =="936c") {
        return {{0, 4}, {0, 1, 2, 3}};

    } else if (diagram =="944c") {
        return {{0, 4}, {2, 3}};

    } else if (diagram =="930c") {
        return {{0, 4}};

    } else if (diagram =="818c") {
        return {{0, 3}, {0, 4}};

    } else if (diagram =="928c") {
        return {{0, 4}, {2, 3}, {0, 1, 2, 3}, {1, 2, 3}};

    } else if (diagram =="808c") {
        return {{1, 3}, {1, 2}, {0, 4}, {1, 0, 4}};

    } else if (diagram =="562c") {
        return {{0, 1}, {0, 2}, {0, 3}, {0, 4}};

    } else {
        throw std::runtime_error("unknown diagram " + diagram);
    }
}

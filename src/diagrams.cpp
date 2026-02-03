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
    } if (diagram =="5c") {
        return {{0, 1}, {0, 2}};

    } if (diagram =="54c") {
        return {{0, 3}};
    } if (diagram =="52c") {
        return {{0, 3}, {1, 2}, {0, 1, 2}};
    } if (diagram =="38c") {
        return {{0, 1}, {0, 3}, {0, 2}};
    }
    if (diagram =="954c") {
        return {{0, 4}};
    } if (diagram =="946c") {
        return {{0, 4}};

    } if (diagram =="952c") {
        return {{0, 4}};

    } if (diagram =="882c") {
        return {};

    } if (diagram =="936c") {
        return {{0, 4}, {0, 1, 2, 3}};

    } if (diagram =="944c") {
        return {{0, 4}, {2, 3}};

    } if (diagram =="930c") {
        return {{0, 4}};

    } if (diagram =="818c") {
        return {{0, 3}, {0, 4}};

    } if (diagram =="928c") {
        return {{0, 4}, {2, 3}, {0, 1, 2, 3}, {1, 2, 3}};

    } if (diagram =="808c") {
        return {{1, 3}, {1, 2}, {0, 4}, {1, 0, 4}};

    } if (diagram =="562c") {
        return {{0, 1}, {0, 2}, {0, 3}, {0, 4}};

    } {
        throw std::runtime_error("unknown diagram " + diagram);
    }
}

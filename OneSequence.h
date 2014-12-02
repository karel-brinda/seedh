#ifndef __ONESEQUENCE_H__
#define __ONESEQUENCE_H__

#include <iostream>
#include <string>
#include <cassert>

#include "Seed.h"
#include "utils.h"

using namespace std;

class OneSequence {

private:
    Seed seed;
    string seq;
    int numberOfPositions;
    int hashTableFirstOcc[HASH_TABLE_SIZE + 1];
    int *hashTableLists;

public:
    bool isKmerPresent(const string &kmer) const;

    OneSequence(string seq, Seed seed);

    ~OneSequence();

    void debugPrintHashTable() const;

    void debugPrintPref() const;
};

#endif

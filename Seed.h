#ifndef __SEED_H__
#define __SEED_H__

#include <iostream>
#include <string>
#include <cassert>

using namespace std;

class Seed {
public:
    Seed(string seed);

    Seed(const Seed &seed);

    int weight;
    string seed;

    void apply(const string &text, string &kmer, const int &applPos) const;
};

#endif

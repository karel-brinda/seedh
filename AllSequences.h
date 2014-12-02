#ifndef __ALLSEQUENCES_H__
#define __ALLSEQUENCES_H__

#include "OneSequence.h"

class AllSequences {
private:
    Seed seed;
    vector<OneSequence *> sequences;
    int seqNb;

public:
    AllSequences(fastaDataT &data, Seed seed);

    ~AllSequences();

    bool isKmerPresent(const string &kmer) const;

    void debugPrintHashTables() const;

    void debugPrintPrefs() const;
};

#endif

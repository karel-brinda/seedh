#include "AllSequences.h"

AllSequences::AllSequences(fastaDataT &data, Seed seed) :
        seed(seed) {
    //cout << "test1" << endl;
    for (int i = 0; i < (int) data.size(); i++) {
        //cerr << "testss" << i << endl;

        OneSequence *seq;
        try {
            seq = new OneSequence(data[i].second, seed);

        } catch (std::bad_alloc &ba) {
            cerr << "bad_alloc caught: " << ba.what() << endl;
            cerr << "file " << __FILE__ << ", line " << __LINE__ << endl;
            exit(1);
        }

        //cout << "testst" << i << endl;
        //cerr << "1" << endl;
        sequences.push_back(seq);
        //cerr << "2" << endl;
    }
    //cout << "test2" << endl;
}

AllSequences::~AllSequences() {
    //cerr << "all sequences destructor" << endl;
    for (int i = 0; i < (int) sequences.size(); i++) {
        delete sequences[i];
    }
}

bool AllSequences::isKmerPresent(const string &kmer) const {
    for (int i = 0; i < (int) sequences.size(); i++) {
        if (sequences[i]->isKmerPresent(kmer)) {
            return true;
        }
    }
    return false;
}

void AllSequences::debugPrintHashTables() const {
    cerr << "==========" << endl;
    for (int i = 0; i < (int) sequences.size(); i++) {
        cerr << "hash table " << i << endl;
        sequences[i]->debugPrintHashTable();
    }
    cerr << "==========" << endl;
}

void AllSequences::debugPrintPrefs() const {
    cerr << "==========" << endl;
    for (int i = 0; i < (int) sequences.size(); i++) {
        sequences[i]->debugPrintPref();
    }
    cerr << "==========" << endl;
}

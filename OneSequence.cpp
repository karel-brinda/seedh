#include "OneSequence.h"

OneSequence::OneSequence(string seq, Seed seed) :
        seed(seed), seq(seq) {

    numberOfPositions = (int) seq.length() - (int) seed.seed.length() + 1;

    try {
        hashTableLists = new int[numberOfPositions];
    } catch (std::bad_alloc &ba) {
        cerr << "bad_alloc caught: " << ba.what() << endl;
        cerr << "file " << __FILE__ << ", line " << __LINE__ << endl;
        exit(1);
    }
    int hash;
    string kmer((unsigned int) seed.weight, ' ');

    /*
     * making a spaghetti array
    */
    int *spaghettiArray;
    int lastOcc[HASH_TABLE_SIZE];
    try {
        spaghettiArray = new int[numberOfPositions];
    } catch (std::bad_alloc &ba) {
        cerr << "bad_alloc caught: " << ba.what() << endl;
        cerr << "file " << __FILE__ << ", line " << __LINE__ << endl;
        exit(1);
    }
    for (int hash = 0; hash < HASH_TABLE_SIZE; hash++) {
        //  no lists at beginning
        hashTableFirstOcc[hash] = -1;
    }

    cerr << "Creating spaghetti index" << endl;
    for (int pointer = 0; pointer < numberOfPositions; pointer++) {
        seed.apply(seq, kmer, pointer);
        hash = getHash(kmer);
        assert(0 <= hash && hash < HASH_TABLE_SIZE);

        if (hashTableFirstOcc[hash] == -1) {
            /* first occurrence of the k-mer */
            lastOcc[hash] = pointer;
            hashTableFirstOcc[hash] = pointer;
            spaghettiArray[pointer] = -1;
        } else {
            /* next occurence of the k-mer */
            int lastPos = lastOcc[hash];
            assert((0 <= lastPos && lastPos < pointer) || !(cerr << "lastPos=" << lastPos << ", pointer=" << pointer << ", number of positions=" << numberOfPositions << endl));
            lastOcc[hash] = pointer;
            spaghettiArray[lastPos] = pointer;
            spaghettiArray[pointer] = -1;
        }

    }
    /*
     * making a nice array
     */
    cerr << "Despaghetting index" << endl;
    int pHashTableLists = 0; // current pointer
    for (int hash = 0; hash < HASH_TABLE_SIZE; hash++) {
        int pSpaghetti = hashTableFirstOcc[hash];
        hashTableFirstOcc[hash] = pHashTableLists;
        while (pSpaghetti != -1) {
            assert((0 <= pSpaghetti && pSpaghetti < numberOfPositions) || !(cerr << "pSpaghetti=" << pSpaghetti << endl));
            assert((0 <= pHashTableLists && pHashTableLists < numberOfPositions) || !(cerr << "pSpaghetti=" << pHashTableLists << endl));

            hashTableLists[pHashTableLists++] = pSpaghetti;
            pSpaghetti = spaghettiArray[pSpaghetti];

        }
    }
    assert(pHashTableLists == numberOfPositions);
    // last stopper
    hashTableFirstOcc[HASH_TABLE_SIZE] = numberOfPositions;

    delete[] spaghettiArray;
    cerr << "Index successfully created" << endl;
}

bool OneSequence::isKmerPresent(const string &kmer) const {
    string testKmer(seed.weight, ' ');
    int hash = getHash(kmer);

    if (hashTableFirstOcc[hash] == -1) {
        return false;
    } else {
        int start = hashTableFirstOcc[hash];
        int end = hashTableFirstOcc[hash + 1];

        for (int j = start; j < end; j++) {
            int cand = hashTableLists[j];
            seed.apply(seq, testKmer, cand);
            if (testKmer.compare(kmer) == 0) {
                return true;
            }
        }

    }
    return false;

}

void OneSequence::debugPrintHashTable() const {

    for (int hash = 0; hash < HASH_TABLE_SIZE; hash++) {
        if (hashTableFirstOcc[hash] != -1) {
            int start = hashTableFirstOcc[hash];
            int end = hashTableFirstOcc[hash + 1];

            if (start != end) {
                cerr << getDehash(hash) << "\t" << start << ".." << end - 1 << ":";

                for (int j = start; j < end; j++) {
                    cerr << "\t" << hashTableLists[j];
                }
                cerr << endl;
            }
        }
    }
}

void OneSequence::debugPrintPref() const {
    cout << seq.substr(0, 150) << "..." << endl;
}


OneSequence::~OneSequence() {
    delete[] hashTableLists;
}

#ifndef _SPLATPOOL_H
#define _SPLATPOOL_H

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "types.h"

#define DEFAULT_MAX_SIZE 10000


std::string open_temp (std::string path, std::fstream& f);

class SplatPool {

    public:
        class Reader;

        SplatPool (void);
        SplatPool (int);
        ~SplatPool (void);
        void flush (void);
        void add (BamAlignment&);
        Reader* reader (void);

    private:
        int MAX_SIZE;
        std::vector<BamAlignment> buffer;
        std::vector<std::string> tempFilenames;
};

class SplatPool::Reader {

    // TODO BamAlignment as key makes sense?
    std::map<BamAlignment, std::fstream*> heads;
    void updateHead (std::fstream&);

    public:
        Reader (std::vector<std::string>&);
        ~Reader (void);
        bool hasNext (void);
        BamAlignment read (void);
};

#endif

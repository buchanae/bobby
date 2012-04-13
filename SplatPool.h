#ifndef _SPLATPOOL_H
#define _SPLATPOOL_H

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "api/BamAlignment.h"
#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "types.h"

#define DEFAULT_MAX_SIZE 100000
#define DEFAULT_TMP_DIR "/tmp"

class SplatPool {

    public:
        class Reader;

        SplatPool (const BamTools::RefVector&);
        SplatPool (int, const BamTools::RefVector&);
        SplatPool (string, const BamTools::RefVector&);
        SplatPool (int, string, const BamTools::RefVector&);
        ~SplatPool (void);
        void flush (void);
        void add (BamTools::BamAlignment&);
        Reader* reader(void);

    private:
        std::string TMP_DIR;
        int MAX_SIZE;
        // TODO maybe buffer should be pointers?
        std::vector<BamTools::BamAlignment> buffer;
        std::vector<std::string> filenames;
        BamTools::RefVector references;
};

class SplatPool::Reader : public BamTools::BamMultiReader {
    public:
        bool GetNextAlignment(BamTools::BamAlignment&);
};

#endif

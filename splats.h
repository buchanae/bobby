#ifndef _SPLATS_H
#define _SPLATS_H

#include "api/BamAlignment.h"

#include "types.h"

splat_pos_t getSplatPosition(BamTools::BamAlignment&);

bool isSplat(BamTools::BamAlignment&);

splat_t bam2splat(BamTools::BamAlignment&, BamTools::RefVector&);

#endif

import argparse
from copy import copy
import random
import string


READ_LENGTH = 60
REFS = ['Chr{}'.format(x) for x in xrange(5)]
REF_LENGTH = 3000000
GAP_RANGE = xrange(100, 10000)

parser = argparse.ArgumentParser()
parser.add_argument('out')


def make_seq():
    return ''.join(random.choice('ATCG') for x in xrange(READ_LENGTH))

def make_ID_base():
    return ''.join(random.choice(string.ascii_letters) for x in xrange(20))


class ID(object):
    def __init__(self, base, mateID='1'):
        self.base = base
        self.mateID = mateID

    def sister(self):
        return ID(self.base, '2' if self.mateID == '1' else '1')

    def __str__(self):
        return '{} {}:N:0:TGACCA'.format(self.base, self.mateID)


class SAM(object):

    def __init__(self, ID, ref, pos, rev_strand=False):
        self.ID = ID
        self.ref = ref
        self.pos = pos
        self.rev_strand = rev_strand

        self.seq = make_seq()
        self.cigar = '{}M'.format(len(self.seq))
        self.extra = []

        self.mapq = '255'
        self.rnext = '*'
        self.pnext = '0'
        self.tlen = '0'
        self.qual = '*'

    @property
    def end_pos(self):
        return self.pos + len(self.seq)

    @classmethod
    def random(cls):
        return cls(ID(make_ID_base()), random.choice(REFS), 
                   random.choice(xrange(REF_LENGTH - READ_LENGTH)))

    @classmethod
    def valid_reg_pair(cls):
        a = cls.random()
        b = cls(a.ID.sister(), a.ref, a.end_pos + random.choice(GAP_RANGE),
                not a.rev_strand)
        return a, b

    @classmethod
    def invalid_distance_reg_pair(cls):
        a, b = cls.valid_reg_pair()
        b.pos = a.end_pos + 50
        return a, b

    @classmethod
    def invalid_strand_reg_pair(cls):
        a, b = cls.valid_reg_pair()
        b.rev_strand = a.rev_strand
        return a, b

    @classmethod
    def invalid_orientation_reg_pair(cls):
        a, b = cls.valid_reg_pair()
        a.rev_strand = True
        b.rev_strand = False
        return a, b

    @classmethod
    def valid_splat_pair(cls):
        a, b = cls.valid_reg_pair()
        a.extra.append('XD:Z:GT-AG')
        a.cigar = '30M200N30M'
        b.pos = a.end_pos + 200 + random.choice(GAP_RANGE)
        b.extra.append('XD:Z:GT-AG')
        b.cigar = '30M200N30M'
        return a, b

    @classmethod
    def valid_mergeable_splat_pair(cls):
        a, b = cls.valid_reg_pair()
        a.extra.append('XD:Z:GT-AG')
        a.cigar = '30M200N30M'
        b.pos = a.end_pos + 200 + random.choice(GAP_RANGE)
        b.extra.append('XD:Z:GT-AG')
        b.cigar = '30M200N30M'

        c = copy(a)
        d = copy(b)

        c.ID = ID(make_ID_base())
        d.ID = c.ID.sister()
        return a, b, c, d

    @classmethod
    def valid_unmergeable_strand_splat_pair(cls):
        a, b, c, d = cls.valid_mergeable_splat_pair()
        c.rev_strand = not a.rev_strand
        d.rev_strand = not b.rev_strand
        return a, b, c, d

    @classmethod
    def valid_reg_splat_pair(cls):
        a, b = cls.valid_reg_pair()
        b.pos = a.end_pos + random.choice(GAP_RANGE)
        b.extra.append('XD:Z:GT-AG')
        b.cigar = '30M200N30M'
        return a, b

    @classmethod
    def valid_splat_reg_pair(cls):
        a, b = cls.valid_reg_pair()
        a.extra.append('XD:Z:GT-AG')
        a.cigar = '30M200N30M'
        b.pos = a.end_pos + 200 + random.choice(GAP_RANGE)
        return a, b

    def __str__(self):
        a = '\t'.join([
            str(self.ID), 
            '16' if self.rev_strand else '0', 
            self.ref, 
            str(self.pos), 
            self.mapq,
            self.cigar,
            self.rnext,
            self.pnext,
            self.tlen,
            self.seq,
            self.qual, 
        ])
        if self.extra:
            a = '\t'.join([a] + self.extra)
        return a


if __name__ == '__main__':

    args = parser.parse_args()
    fh = open(args.out, 'w')

    # generate random SAM data
    out = []
    for x in xrange(5):
        #out.extend(SAM.valid_reg_pair())
        #out.extend(SAM.valid_splat_pair())
        #out.extend(SAM.valid_mergeable_splat_pair())
        #out.extend(SAM.valid_reg_splat_pair())
        #out.extend(SAM.valid_splat_reg_pair())
        #out.extend(SAM.invalid_distance_reg_pair())
        #out.extend(SAM.invalid_strand_reg_pair())
        #out.extend(SAM.invalid_orientation_reg_pair())
        out.extend(SAM.valid_unmergeable_strand_splat_pair())

    # print SAM header
    fh.write('@HD\tVN:1.0\n')
    for ref in REFS:
        fh.write('@SQ\tSN:{}\tLN:{}\n'.format(ref, REF_LENGTH))

    # randomize output order and write output
    #random.shuffle(out)
    for o in out:
        fh.write(str(o))
        fh.write('\n')

import argparse
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

    def __init__(self, ID, ref, pos, rev_strand=False, extra=''):
        self.ID = ID
        self.ref = ref
        self.pos = pos
        self.rev_strand = rev_strand

        self.seq = make_seq()
        self.cigar = '{}M'.format(len(self.seq))
        self.extra = extra

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

    def sister(self):
        return SAM(self.ID.sister(), self.ref, self.end_pos + random.choice(GAP_RANGE),
                   not self.rev_strand, self.extra)

    def __str__(self):
        a = '\t'.join([
            str(self.ID), 
            '0' if self.rev_strand else '16', 
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
            a = '\t'.join([a, self.extra])
        return a


if __name__ == '__main__':

    args = parser.parse_args()
    fh = open(args.out, 'w')

    # generate random SAM data
    out = []
    for x in xrange(5):
        s = SAM.random()
        out.append(s)
        out.append(s.sister())


    # print SAM header
    fh.write('@HD\tVN:1.0\n')
    for ref in REFS:
        fh.write('@SQ\tSN:{}\tLN:{}\n'.format(ref, REF_LENGTH))

    # randomize output order and write output
    #random.shuffle(out)
    for o in out:
        fh.write(str(o))
        fh.write('\n')

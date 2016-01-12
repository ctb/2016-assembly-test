#! /usr/bin/env python
from intervaltree import Interval, IntervalTree
import numpy
import screed

def _load_coords(filename, only=None):
    lines = [ x.strip() for x in (open(filename)) ]
    assert lines[1].startswith('NUCMER'), lines[0]

    coords = []
    for line_no in range(5, len(lines)):
        line = lines[line_no].split()
        s1, e1 = int(line[0]), int(line[1])
        s2, e2 = int(line[3]), int(line[4])
        ident = float(line[9])
        name1, name2 = line[11], line[12]
        s1, e1 = min(s1, e1), max(s1, e1)
        s2, e2 = min(s2, e2), max(s2, e2)

        yield (s1, e1, s2, e2, ident, name1, name2)

def load_refsizes(filename):
    d = {}
    for record in screed.open(filename):
        d[record.name.split()[0]] = len(record.sequence)
    return d

class GenomeIntervalsContainer(object):
    def __init__(self, refsizes):
        self.names = refsizes.keys()
        self.refsizes = refsizes
        trees = {}
        for k in self.names:
            trees[k] = numpy.zeros(refsizes[k])
        self.trees = trees

    def load_coords(self, filename, min_length=100, min_ident=99.0):
        trees = self.trees
        for s1, e1, s2, e2, ident, name1, name2 in _load_coords(filename):
            if name1 in trees:
                tree = trees[name1]
                if e1 - s1 + 1 >= min_length and ident >= min_ident:
                    tree[s1 - 1:e1] = numpy.ones(e1 - s1 + 1)

    def calc_uncov(self):
        uncov_d = {}
        for name in self.refsizes:
            length = self.refsizes[name]
            cov = sum(self.trees[name])

            uncov = length - cov
            uncov_d[name] = uncov
        return uncov_d

    def subtract(self, other):
        assert self.refsizes.keys() == other.refsizes.keys()
        new_gic = GenomeIntervalsContainer(self.refsizes)
        otrees = other.trees
        
        for name in self.refsizes:
            ttree = self.trees[name]
            otree = otrees[name]
            newtree = ttree - numpy.logical_and(ttree, otree)
            new_gic.trees[name] = newtree

        return new_gic

    def union(self, other):
        assert self.refsizes.keys() == other.refsizes.keys()
        new_gic = GenomeIntervalsContainer(self.refsizes)
        otrees = other.trees
        
        for name in self.refsizes:
            ttree = self.trees[name]
            otree = otrees[name]
            newtree = numpy.logical_or(ttree, otree)
            new_gic.trees[name] = newtree

        return new_gic
    
    def intersect(self, other):
        assert self.refsizes.keys() == other.refsizes.keys()
        new_gic = GenomeIntervalsContainer(self.refsizes)
        otrees = other.trees
        
        for name in self.refsizes:
            ttree = self.trees[name]
            otree = otrees[name]
            newtree = numpy.logical_and(ttree, otree)
            new_gic.trees[name] = newtree

        return new_gic

def main():
    print 'loading refsizes'
    refsizes = load_refsizes('mircea.fa')

    refsizes2 = {}
    refsizes2['Thermoanaerobacter_pseudethanolicus_ATCC_33223'] = \
       refsizes['Thermoanaerobacter_pseudethanolicus_ATCC_33223']
    refsizes = refsizes2
    
    print 'loading coords'
    gic = GenomeIntervalsContainer(refsizes)
    gic.load_coords('idba.coords')
    uncov_d = gic.calc_uncov()
    
    for name in uncov_d:
        print name, uncov_d[name]

    total_bp = sum(gic.refsizes.values())
    total_uncov = sum(uncov_d.values())
    print total_bp, total_uncov, float(total_uncov) / total_bp

    gic2 = GenomeIntervalsContainer(refsizes)
    gic2.load_coords('megahit.coords')

    sub = gic.subtract(gic2)
    print 'subtract', total_bp, sum(sub.calc_uncov().values())
    
    union = gic.union(gic2)
    print 'union', total_bp, sum(union.calc_uncov().values())
    
    intersect = gic.intersect(gic2)
    print 'intersect', total_bp, sum(intersect.calc_uncov().values())

    
if __name__ == '__main__':
    main()

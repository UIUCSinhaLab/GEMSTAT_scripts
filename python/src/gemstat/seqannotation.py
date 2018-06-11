import numpy as _np

class Site(object):
    def __init__(self, start, end, orientation, motif, energy, wtRatio=None):
        self.start = start
        self.end = end
        self.orientation = orientation
        self.motif = motif
        self.energy = energy
        if None != wtRatio:
            self.wtRatio = wtRatio
        else:
            self.wtRatio = _np.exp(-self.energy)

    def __repr__(self):
        return "{}\t{}\t{}\t{}\t{}".format(self.start,
                                       "+" if self.orientation else "-",
                                       self.motif.name,
                                       self.energy,
                                       self.wtRatio)


def annotate_sequence(seq, motifs, et=0.5, filter=None):

    #HANDLE ET
    if None == et:
        et = 0.99
    if _np.isscalar(et):
        et = [et for i in range(len(motifs))]
    if isinstance(et,list):
        et = _np.array(et)
    if len(motifs) != et.size:
        raise Exception("The energy thresholds did not match the number of motifs. (You can also use a single scalar)")
    #END OF ET SETUP

    annotations = list()

    for one_motif, one_et in zip(motifs,et):
        L = len(seq)
        bestLLR = one_motif.pssm.calculate(one_motif.consensus)

        hits = one_motif.pssm.search(seq)

        hits = [Site(
                        (L + i if i < 0 else i) + 1,
                        (L + i if i < 0 else i) + 1 + one_motif.length - 1,
                        False if i < 0 else True,
                        one_motif,
                        -j+bestLLR
            ) for i,j in hits]


        if None != filter:
            hits = [i for i in hits if filter(i)]
        else:
            hits = [i for i in hits if i.energy <= one_et*bestLLR]

        annotations.extend(hits)

    #sort the annotations first.
    annotations.sort(key=lambda x:x.start)

    return annotations

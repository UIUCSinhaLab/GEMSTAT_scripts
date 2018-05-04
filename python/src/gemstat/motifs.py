import os as _os
import re as _re
import numpy as _np

import Bio.motifs
import Bio.motifs.matrix
import Bio.Alphabet

def to_np_array(inmat,letters="ACGT"):
    return _np.array([[inmat[i][j] for i in letters] for j in range(inmat.length)])

def np_to_biopython(matrix, letters="ACGT"):
    return dict([(letters[l],matrix[:,l]) for l in range(len(letters))])

class GEMSTAT_Motif(Bio.motifs.Motif):

    def __init__(self, name, counts, pseudocount=0.0, comments=None):

        #TODO: Check that the user hasn't provided these in a more BioPython way.
        counts_np = _np.array(counts)#In case the user provided it as a list of lists.

        alphabet = Bio.Alphabet.IUPAC.unambiguous_dna;
        alphabet.letters = "ACGT"

        counts_as_biopython_wants = dict([(alphabet.letters[i],counts_np[:,i]) for i in range(len(alphabet.letters))])

        super(GEMSTAT_Motif,self).__init__(alphabet=alphabet, counts = np_to_biopython(counts_np,alphabet.letters))

        self.name = name
        self.comments = comments if None != comments else ""
        if self.comments.startswith("#"):
            self.comments = self.comments[1:]
        self.pseudocount = pseudocount
        self.pseudocounts = pseudocount



    @property
    def pssm(self):
        np_pwm = to_np_array(self.pwm,self.alphabet.letters)
        np_background = _np.array([[self.background[l] for l in self.alphabet.letters] for i in range(np_pwm.shape[0])])

        return Bio.motifs.matrix.PositionSpecificScoringMatrix(self.alphabet,
                                             np_to_biopython(
                                                 _np.log( np_pwm/np_background

                                                 ) )
                                            )

    def __format_gemstat(self):
        str_list = list()
        str_list.append(">{}\t{}\t{}\t#{}".format(self.name, self.length, self.pseudocount,self.comments))

        for i in range(self.length):
            str_list.append("\t".join([str(self.counts[l][i]) for l in self.alphabet.letters]))

        str_list.append("<")
        return "\n".join(str_list)

    def format(self,*args):
        if 0 == len(args) or args[0] == "gemstat":
            return self.__format_gemstat()
        else:
            return super(GEMSTAT_Motif,self).format(*args)




def read_gemstat_motifs(filelike_or_filename):
    float_re = "(?:[0-9]*[.])?[0-9]+(?:[eE][+-]?[0-9]+)?"

    header_pat = _re.compile("^>(?P<name>\S+)\s+(?P<length>\S+)\s+(?P<pseudocount>" +
                             float_re + ")\s*(?P<comments>.*)$")
    ALPH = "ACGT"
    line_pat_string = "^\s*" + "\s+".join(["(?P<{}>{})".format(i,float_re) for i in ALPH]) + "\s*$"
    line_pat = _re.compile( line_pat_string )

    def make_motif(in_dict):
        #counts_as_biopython_wants = dict( zip(ALPH, map(list, zip(*in_dict["counts"])) ) )
        return GEMSTAT_Motif(in_dict["name"],
                             in_dict["counts"],
                             pseudocount=in_dict["pseudocount"],
                             comments=in_dict["comments"]
                            )


    def gemstat_read_helper(f):
        """This probably doesn't like streams."""

        #
        helper_motif_list = list()
        #read a matrix, but note that Bio.motif needs a dictionary of the counts keyed by the alphabet.
        #like {"A":[1,5],"C":[2,6],"G":[3,7],"T":[4,8]} for a 2 position PWM with 1 A in the first and 5 in the second positions.

        next_motif = None
        state = 0# state 0 is waiting for the first line of a motif, state 1 is reading lines (waiting for last)

        def check_eof(ff):
            #BEGIN EOF SNIPPET FROM : https://stackoverflow.com/a/48679747
            return ff.tell() == _os.fstat(ff.fileno()).st_size
            #END SNIPPET

        try:
            while( not check_eof(f) ):
                aline = f.readline()

                if 0 == state:
                    header_match = header_pat.match(aline)
                    if None == header_match:
                        raise Exception("Error parsing motif header: {}".format(aline))
                    #start creating a new motif.
                    next_motif = {"name":header_match["name"],
                                  "header_length":header_match["length"],
                                  "pseudocount":header_match["pseudocount"],
                                    "comments":header_match["comments"],
                                 "counts":list()}
                    state = 1
                    continue
                elif 1 == state:
                    if aline.startswith("<"):
                           helper_motif_list.append(make_motif(next_motif))
                           next_motif = None
                           state = 0
                           continue
                    #normal line
                    line_match = line_pat.match(aline)
                    if None == line_match:
                        raise Exception("Problem processing line: {}".format(aline))
                    one_line_counts = [float(line_match[i]) for i in ALPH]
                    next_motif["counts"].append(one_line_counts)
                    continue

        except StopIteration:
            pass
        return helper_motif_list



    #Actual function body.

    motiflist = None

    if isinstance(filelike_or_filename, str):
        with open(filelike_or_filename,"r") as infile:
            motiflist = gemstat_read_helper(infile)

    else:
        motiflist = gemstat_read_helper(filelike_or_filename)

    return motiflist

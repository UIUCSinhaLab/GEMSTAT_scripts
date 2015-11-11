import scipy as _S
from Bio.Alphabet.IUPAC import unambiguous_dna as _unamb_dna
from Bio.motifs import Motif as _MOT
import re as _re
import StringIO as _SIO

class GemPWM(_MOT):
	def __init__(self,name,counts,pseudocount="0.1"):
		self.name = name
		countsdict = dict(zip("ACGT",counts.T+pseudocount))
		super(GemPWM, self).__init__(alphabet=_unamb_dna,counts=countsdict)

	@classmethod
	def parse(cls, filename):
		readall = None
		with open(filename,"r") as infile:
			readall = infile.read()


		getter = _re.compile("""^>(\w+)\s*(\d+)\s*([0-9.]+)\s*$([^<]+)<""",_re.MULTILINE | _re.DOTALL)
		for one_motif in getter.findall(readall):
			name, length, pseudocount, counts = one_motif
			counts = _S.loadtxt(_SIO.StringIO(counts))
			yield cls(name, counts, float(pseudocount))

import scipy as _S
from Bio.Alphabet.IUPAC import unambiguous_dna as _unamb_dna
from Bio.motifs import Motif as _MOT
from Bio.Seq import Seq as _Seq
import re as _re
import StringIO as _SIO

class GemPWM(_MOT):
	def __init__(self,name,counts,pseudocount="0.1"):
		self.name = name
		countsdict = dict(zip("ACGT",counts.T+pseudocount))
		super(GemPWM, self).__init__(alphabet=_unamb_dna,counts=countsdict)

	def get_consensus_score(self):
		return self.pssm.calculate(self.consensus) / _S.log2(_S.e)

	def calculate(self, instr):
		vals = None
		if isinstance(instr,_Seq):
			vals = self.pssm.calculate(instr)
		else:
			vals = self.pssm.calculate(_Seq(instr,_unamb_dna))
		energy = self.pssm.calculate(self.consensus) -vals
		return energy /_S.log2(_S.e) # biopython uses log base 2, but GEMSTAT uses log base e #TODO: Make it automatically determine what the base of biopython log is by creating a special pwm.

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

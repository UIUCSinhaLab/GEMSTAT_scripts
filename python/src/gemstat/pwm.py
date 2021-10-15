import numpy as _np
import scipy as _S
#from Bio.Alphabet.IUPAC import unambiguous_dna as _unamb_dna
_unamb_dna = "ACGT"
from Bio.motifs import Motif as _MOT
from Bio.Seq import Seq as _Seq
import re as _re
from io import StringIO

class GemPWM(_MOT):
	def __init__(self,inname,counts,pseudocount=0.1,comment=""):
		
		pseudocount = float(pseudocount)
		countsdict = dict(zip("ACGT",counts.T+pseudocount))
		super(GemPWM, self).__init__(alphabet=_unamb_dna,counts=countsdict)
		
		self.name = inname
		self.comment = comment

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
		"""
		"""
		
		#thanks to : https://www.regular-expressions.info/floatingpoint.html
		float_pattern_str = "[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?"
		float_pattern = _re.compile(float_pattern_str)
		
		#motif_begin  = _re.compile("^>",_re.DOTALL)
		motif_header = _re.compile("^>(?P<name>\w+)\s*(?P<length>\d+)\s*(?P<pseudocount>" + float_pattern_str + ")\s*(?P<comment>#.*)?$",_re.DOTALL)
		motif_line   = _re.compile("^([^<]+)$",_re.DOTALL)
		motif_end    = _re.compile("^<\s*$",_re.DOTALL)
		
		with open(filename,"r") as infile:
			l = infile.readline()
			while l:
				head_match = motif_header.match(l)
				if not head_match:
					raise Exception("Unexpected line in motif file: {} ".format(filename))
				name, length, pseudo, comment = head_match.groups()
				
				#counts_sio = StringIO()
				
				l = infile.readline()
				count_list = list()
				while l:
					if motif_end.match(l):
						break
					#otherwise
					count_list.append(l)
					l = infile.readline()
				
				counts = _np.loadtxt(StringIO("\n".join(count_list)))
				
				yield cls(name,counts,float(pseudo),comment)
				
				l = infile.readline()

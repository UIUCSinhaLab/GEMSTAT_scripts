
import numpy as np

from . import snot

#Inputs
class Factor(object):
	def __init__(self):
		pass

class TranscriptionFactor(Factor):
	def __init__(self, name, conc_coeff=1.0, action=1.0, activator=True, repressor=False):
		self.name = name
		self.conc_coeff = conc_coeff
		self.action = action
		self.activator = activator
		self.repressor = repressor

#Interactions
class FactorInteraction(object):
	def __init__(self):
		pass
	#TODO: More methods, such as "get interaction" which takes the positions of two TBFS, etc.

class CooperativeInteraction(FactorInteraction):
	def __init__(self, factor1, factor2, coopertivity):
		self.factor1 = factor1
		self.factor2 = factor2
		self.coopertivity = coopertivity


#Enhancers and stuff
class Enhancer(object):
	def __init__(self,name,seq=None):
		self.name = name
		self.seq = seq

class Promoter(object):
	def __init__(self,name,basal,seq=None):
		self.name = name
		self.basal = basal
		self.seq = seq

class OutputGene(object):
	def __init__(self,name,enhancers=list(),promoter=None):
		self.name = name
		self.enhancers = enhancers
		self.promoter = promoter

class Model(object):
	def __init__(self):
		self.factors = list()
		self.factor_names_to_positions = dict()
		self.output_genes = list()
		self.pis = list()#NEVER USE PIS
		self.interactions = dict()
		self.betas = dict()#why not roll this into output gene?
		self.factor_thresholds = dict()#not exactly part of the model.

	def add_factor(self, in_fact,threshold=0.5):
		self.factors.append(in_fact)
		self.factor_names_to_positions[in_fact.name] = len(self.factors)-1
		self.factor_thresholds.append(threshold)

	def add_output(self, new_out,beta=1.0,pi=1e-50):
		self.output_genes.append(new_out)
		self.betas.append(beta)
		self.pis.append(pi)

	def add_interaction(self,new_inter):
		self.interactions.append(new_inter)

from collections import OrderedDict

class SetDict(dict):
	def _keyconvert(self,key):
		if isinstance(key,str):
			key = [key]
		return frozenset(key)
	
	def __setitem__(self,key,item):
		return super(SetDict,self).__setitem__(self._keyconvert(key),item)
	
	def __getitem__(self,key):
		return super(SetDict,self).__getitem__(self._keyconvert(key))
	

class ParFile(object):
	"""Reads an intermediate representation of a par file. That can then be transformed into a model.
	Can more easily transform a model into a par file.
	"""
	@classmethod
	def read_1_6a(cls, filename):
		tf_lines = list()
		basal_trans = list()
		pis = list()
		betas = list()
		coop_lines = list()
		annotation_cutoffs = None

		tmp_line = None
		with open(filename, "r") as infile:
			#TF Lines
			header = infile.next()
			if "#GSPAR1.6a" != header.strip().split()[0]:
				raise Exception("File does not contain the correct header!")

			for line in infile:
				if line.startswith("basal_transcription ="):
					tmp_line = line
					break
				tf_lines.append(line)
			#basal_transcription
			basal_trans = tmp_line[tmp_line.find("=")+1:]
			tmp_line = None
			#other promoter stuff
			pis = infile.next()
			betas = infile.next()
			#cooperativity values
			for line in infile:
				coop_lines.append(line)
			#past_eof
			annotation_cutoffs = coop_lines.pop(-1)
		#end
		tf_lines = [i.strip().split() for i in tf_lines]
		tf_lines = [(i[0],map(float,i[1:])) for i in tf_lines]
		basal_trans = map(float,basal_trans.strip().split())
		pis = map(float,pis.strip().split())
		betas = map(float,betas.strip().split())
		coop_lines = [i.strip().split() for i in coop_lines]
		coop_lines = [((a,b),float(c)) for a,b,c in coop_lines]
		annotation_cutoffs = map(float,annotation_cutoffs.strip().split())

		if len(pis) != len(betas):
			raise Exception("number of pis and betas disagree while reading parfile.")

		if len(annotation_cutoffs) != len(tf_lines):
			raise Exception("TF lines and cutoffs disagree in number while reading parfile.")

		#Using an OrderedDict ensures that it prints out in the same order as the parfile.
		retdict = OrderedDict()
		retdict["tfs"] = OrderedDict(tf_lines)
		retdict["qbtms"] = basal_trans
		retdict["pis"] = pis
		retdict["betas"] = betas
		retdict["coops"] = coop_lines
		retdict["thresholds"] = annotation_cutoffs
		return retdict



	@classmethod
	def read_old_par(cls, filename,strict=True):
		retmodel = Model()

		#read lines until basal_transcription
		#read basal_transcrptions into temporary storage.
		#read several floats, one per line.
		#if we see a line as "str	str	float", then it is a coopertivity. start reading coopertivities.
		#if we see a line that has multiple floats, it is the factor thresholds. Confirm that the number of factors
		#if we do not see the previous line, (and mode is strict), Fail.

		def notfloat(in_str):
			try:
				float(in_str)
				return False
			except:
				return True

		class putback_buffer(object):
			def __init__(self, thing):
				self.lines = [l for l in thing]
				self.SIZE = len(self.lines)
				self.position = 0

			def __iter__(self):
				return self

			def __next__(self):
				if self.position == self.SIZE:
					raise StopIteration
				else:
					self.position += 1
					return self.lines[self.position - 1]

			def next(self):
				return self.__next__()
			def putback(self, item):
				self.position-=1

		tf_buffer = list()
		N_TF = -1
		qbtm = list()
		float_buffer = list()
		coop_buffer = list()
		num_thresh = 0
		saw_basal = False
		with open(filename,"r") as infile:
			pb = putback_buffer(infile)
			for line in pb:
				cleanline = line.strip()
				if line.startswith("basal_transcription"):
					qbtm = list(map(float,cleanline.split("=")[1].split()))
					saw_basal = True
					break
				name, other = cleanline.split(None,1)
				tf_buffer.append((name,list(map(float,other.split()))))

			N_TF = len(tf_buffer)
			for line in pb:
				#a single float?
				accepted = False
				try:
					mult_vals = None
					if strict:
						mult_vals = [float(line)]
					else:
						mult_vals = list(map(float,line.strip().split()))
					float_buffer.extend(mult_vals)
					accepted = True
				except:
					pb.putback(line)
					break

			#we are at coopertivities?
			for line in pb:
				#is it a coopertivity?
				brokenline = line.split()
				if len(brokenline) == 3 and notfloat(brokenline[0]) and notfloat(brokenline[1]) and not notfloat(brokenline[2]):
					coop_buffer.append((brokenline[0],brokenline[1],float(brokenline[2])))
				else:
					pb.putback(line)
					break
			#now we are definitely at thresholds
			for line in pb:
				foo = list(map(float,line.strip().split()))
				num_thresh += len(foo)
				float_buffer.extend(foo)

			#CONGRATS, we reached the end of the file!
			#print("TFS: ", len(tf_buffer), N_TF)
			#print("QBTM: ", qbtm)
			#print("FLOATS: ", len(float_buffer), float_buffer)
			#print("saw basal:, ", saw_basal)
			#print("COOPS: ", len(coop_buffer))

			if not saw_basal:
				raise Exception("No basal transcription line while parsing .par file")
			if not (len(float_buffer)-N_TF)%2==0:
				raise Exception("The number of betas and pis is inconsistent")
			if len(qbtm) != 1 and len(qbtm) != (len(float_buffer)-N_TF)%2:
				raise Exception("The number of qbtms, betas, and pis did not match")
			thresh_out = float_buffer[-N_TF:]
			float_buffer = float_buffer[:-N_TF]
			return {"tfs":tf_buffer, "qbtms":qbtm,"pis":float_buffer[:int(len(float_buffer)/2)],"betas":float_buffer[int(len(float_buffer)/2):],"coops":coop_buffer,"thresholds":thresh_out}

def parfile_to_vector(par_dict,factor_roles=None,include_pis=False):
	storage = list() # binding, coop, act, rep, qbtm, PIS(not used), betas, thresholds

	all_bind = list()
	all_act = list()
	all_rep = list()

	all_coops = list()

	tf_names = [i[0] for i in par_dict["tfs"]]

	for name,(bind,act,rep) in par_dict["tfs"]:
		all_bind.append(bind)
		all_act.append(act)
		all_rep.append(rep)

	for name1,name2,coopval in par_dict["coops"]:
		ind1 = tf_names.index(name1)
		ind2 = tf_names.index(name2)

		if ind2 < ind1:
			ind1,ind2 = ind2,ind2
		all_coops.append(((ind1,ind2),coopval))
	all_coops.sort(key=lambda x:x[0][1])
	all_coops.sort(key=lambda x:x[0][0])

	storage.extend(all_bind)
	storage.extend([i[1] for i in all_coops])
	if factor_roles == None:
		storage.extend(all_act)
		storage.extend(all_rep)
	else:
		storage.extend([i for i,j in zip(all_act,tf_names) if factor_roles[j][0]])
		storage.extend([i for i,j in zip(all_rep,tf_names) if factor_roles[j][1]])
	storage.extend(par_dict["qbtms"])
	if include_pis:
		storage.extend(par_dict["pis"])
	storage.extend(par_dict["betas"])
	storage.extend(par_dict["thresholds"])

	return np.array(storage)

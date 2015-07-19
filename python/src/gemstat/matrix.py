from __future__ import print_function

import numpy as NP
import scipy as S
import scipy.io as SIO

def renumber_curve(in_curve,newN,s=0.001):
    import scipy.interpolate as INTERP
    
    orig_ndim = in_curve.ndim
    as_matrix = S.array(in_curve,ndmin=2)
    N,M = as_matrix.shape
    
    old_x = S.linspace(0.0,1.0,M)
    new_x = S.linspace(0.0,1.0,newN)
    
    spline_reps = [INTERP.splrep(old_x,row,s=s) for row in as_matrix]
    reconstituted = S.array([INTERP.splev(new_x,sp_rep) for sp_rep in spline_reps])
    if orig_ndim == 1:
        return reconstituted[0]
    else:
        return reconstituted

class GEMSTAT_Matrix(object):
	def __init__(self):
		self.storage = None
		self.names = S.array([],dtype='|S4')
		self.names_to_rows = dict()
	
	def add_row(self, name, data):
		self.names = S.hstack([self.names, S.array([name],dtype='|S4')])
		self.names_to_rows[name] = len(self.names) - 1

		if self.storage == None:
			self.storage = data.reshape(1,-1)
		else:
			self.storage = S.vstack([self.storage, data])
	
	def __getitem__(self, index):
		if isinstance(index,str):
			#get row(s) by name
			return self.storage[self.names == index]
		else:
			return self.storage[index]

	
	def add_column(self, column):
		self.storage = S.hstack([self.storage, column.reshape(-1,1)])
	
	def hstack_update(self, other):
		desired_order_other = [other.names_to_rows[i] for i in self.names]
		self.storage = S.hstack([self.storage, other[desired_order_other]])
	
	def write(self, filename, format="%.5f"):
		outfile = open(filename,"w") 
		#write the header line
		
		NUMROWS, NUMCOLS = self.storage.shape
		cols = range(1,NUMCOLS+1)
		
		print("ROWS\t" + "\t".join(map(str,cols)),file=outfile)
		
		formatter = lambda x: format % x
		
		for one_name,one_expr in zip(self.names,self.storage):
			one_line_out_str = one_name+"\t"+ "\t".join(map(formatter, one_expr))
			print(one_line_out_str,file=outfile)
		
		outfile.close()
	
	@classmethod
	def load(cls, filename):
		retmat = GEMSTAT_Matrix()
		
		tmp_names = list()
		tmp_namemap = dict()
		def _convert_and_store_name(in_str):
			in_str = str(in_str).decode("ASCII").strip()
			position = len(tmp_names)
			tmp_names.append(in_str)
			tmp_namemap[in_str] = position
			return position
		
		stuff   = S.loadtxt(filename,converters={0:_convert_and_store_name})[:,1:]
		COLNUMS = S.array(stuff[0],dtype=S.int_)
		DATA    = stuff[1:]
		
		#Sanity Check
		import pdb
		pdb.set_trace()
		if "ROWS" != tmp_names[0].upper():
			raise Exception("Currently the GEMSTAT_Matrix parser requires the first row be named 'ROWS' (case insensitive).")
		if any([i != j for i,j in zip(COLNUMS,range(1,len(COLNUMS)+1))]):
			raise Exception("Currently the GEMSTAT_Matirx parser requires that the column numbers be contiguous integers 1-N, sorry.")
		
		retmat.names = tmp_names[1:]
		retmat.names_to_rows = dict([(retmat.names[i], i) for i in range(len(retmat.names))])
		retmat.storage = DATA
		return retmat


if __name__ == "__main__":
	F = FactorExpr()
	F.add("test", S.randn(10))
	F.write("/tmp/factor_test.txt")


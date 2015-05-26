from __future__ import print_function

import scipy as S
import scipy.io as SIO

class GEMSTAT_Matrix(object):
	def __init__(self):
		self.storage = None
		self.names = list()
		self.names_to_rows = dict()
	
	def add_row(self, name, data):
		self.names.append(name)
		self.names_to_rows[name] = len(self.names) - 1

		if self.storage == None:
			self.storage = data.reshape(1,-1)
		else:
			self.storage = S.vstack([self.storage, data])
	
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
		stuff = S.loadtxt(filename,dtype=S.str_ )
		header = stuff[0]
		the_rest = stuff[1:]
			
		retmat.names = list(the_rest[:,0])
		retmat.names_to_rows = dict([(retmat.names[i], i) for i in range(len(retmat.names))])
		retmat.storage = S.array(the_rest[:,1:],dtype=S.float_)
		return retmat


if __name__ == "__main__":
	F = FactorExpr()
	F.add("test", S.randn(10))
	F.write("/tmp/factor_test.txt")


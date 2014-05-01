# junk_aligner.py
import numpy

class JunkAligner:
	"""
	Class for aligning regions of sequences where there
	are no common motifs of a length greater than 3
	"""
	def __init__(self):
		"""
		Constructor for JunkAligner, merely initializes
		the blosum matrix.
		"""
		self.buildBlosum();

	def buildBlosum(self):
		"""
		Loads the Blosum45 substitution matrix
		"""
		self.blosum45 = {};
		self.blosum45["A"]= { "A":5,  "R":-2,  "N":-1,"D": -2,"C": -1,"Q": -1,
						"E": -1,"G":  0,"H": -2,"I": -1,"L": -1,"K": -1,"M": -1,
						"F": -2,"P": -1,"S":  1,"T":  0,"W": -2,"Y": -2,"V":  0,
						"B": -1,"J": -1,"Z": -1,"X": -1 };
		self.blosum45["R"]={ "A":-2, "R":  7, "N": 0,"D": -1,"C": -3,"Q":  1,"E":  0,
						"G": -2,"H":  0,"I": -3,"L": -2,"K":  3,"M": -1,"F": -2,
						"P": -2,"S": -1,"T": -1,"W": -2,"Y": -1,"V": -2,"B": -1,
						"J": -3,"Z":  1,"X": -1 }
		self.blosum45["N"]={ "A":-1 ,"R": 0 , "N": 6,"D":  2,"C": -2,"Q":  0,
						"E":  0,"G":  0,"H":  1,"I": -2,"L": -3,"K":  0,
						"M": -2,"F": -2,"P": -2,"S":  1,"T":  0,"W": -4,
						"Y": -2,"V": -3,"B":  5,"J": -3,"Z":  0,"X": -1 }
		self.blosum45["D"]={ "A":-2 ,"R":-1 , "N": 2,"D":  7,"C": -3,"Q":  0,"E":  2,"G": -1,"H":  0,"I": -4,"L": -3,"K":  0,"M": -3,"F": -4,"P": -1,"S":  0,"T": -1,"W": -4,"Y": -2,"V": -3,"B":  6,"J": -3,"Z":  1,"X": -1 }
		self.blosum45["C"]={ "A":-1 ,"R":-3 , "N":-2,"D": -3,"C": 12,"Q": -3,"E": -3,"G": -3,"H": -3,"I": -3,"L": -2,"K": -3,"M": -2,"F": -2,"P": -4,"S": -1,"T": -1,"W": -5,"Y": -3,"V": -1,"B": -2,"J": -2,"Z": -3,"X": -1 }
		self.blosum45["Q"]={ "A":-1 ,"R": 1 , "N": 0,"D":  0,"C": -3,"Q":  6,"E":  2,"G": -2,"H":  1,"I": -2,"L": -2,"K":  1,"M":  0,"F": -4,"P": -1,"S":  0,"T": -1,"W": -2,"Y": -1,"V": -3,"B":  0,"J": -2,"Z":  4,"X": -1}
		self.blosum45["E"]={ "A":-1 ,"R": 0 , "N": 0,"D":  2,"C": -3,"Q":  2,"E":  6,"G": -2,"H":  0,"I": -3,"L": -2,"K":  1,"M": -2,"F": -3,"P":  0,"S":  0,"T": -1,"W": -3,"Y": -2,"V": -3,"B":  1,"J": -3,"Z":  5,"X": -1}
		self.blosum45["G"]={ "A": 0 ,"R":-2 , "N": 0,"D": -1,"C": -3,"Q": -2,"E": -2,"G":  7,"H": -2,"I": -4,"L": -3,"K": -2,"M": -2,"F": -3,"P": -2,"S":  0,"T": -2,"W": -2,"Y": -3,"V": -3,"B": -1,"J": -4,"Z": -2,"X": -1}
		self.blosum45["H"]={ "A":-2 ,"R": 0 , "N": 1,"D":  0,"C": -3,"Q":  1,"E":  0,"G": -2,"H": 10,"I": -3,"L": -2,"K": -1,"M":  0,"F": -2,"P": -2,"S": -1,"T": -2,"W": -3,"Y":  2,"V": -3,"B":  0,"J": -2,"Z":  0,"X": -1}
		self.blosum45["I"]={ "A":-1 ,"R":-3 , "N":-2,"D": -4,"C": -3,"Q": -2,"E": -3,"G": -4,"H": -3,"I":  5,"L":  2,"K": -3,"M":  2,"F":  0,"P": -2,"S": -2,"T": -1,"W": -2,"Y":  0,"V":  3,"B": -3,"J":  4,"Z": -3,"X": -1}
		self.blosum45["L"]={ "A":-1 ,"R":-2 , "N":-3,"D": -3,"C": -2,"Q": -2,"E": -2,"G": -3,"H": -2,"I":  2,"L":  5,"K": -3,"M":  2,"F":  1,"P": -3,"S": -3,"T": -1,"W": -2,"Y":  0,"V":  1,"B": -3,"J":  4,"Z": -2,"X": -1}
		self.blosum45["K"]={ "A":-1 ,"R": 3 , "N": 0,"D":  0,"C": -3,"Q":  1,"E":  1,"G": -2,"H": -1,"I": -3,"L": -3,"K":  5,"M": -1,"F": -3,"P": -1,"S": -1,"T": -1,"W": -2,"Y": -1,"V": -2,"B":  0,"J": -3,"Z":  1,"X": -1}
		self.blosum45["M"]={ "A":-1 ,"R":-1 , "N":-2,"D": -3,"C": -2,"Q":  0,"E": -2,"G": -2,"H":  0,"I":  2,"L":  2,"K": -1,"M":  6,"F":  0,"P": -2,"S": -2,"T": -1,"W": -2,"Y":  0,"V":  1,"B": -2,"J":  2,"Z": -1,"X": -1}
		self.blosum45["F"]={ "A":-2 ,"R":-2 , "N":-2,"D": -4,"C": -2,"Q": -4,"E": -3,"G": -3,"H": -2,"I":  0,"L":  1,"K": -3,"M":  0,"F":  8,"P": -3,"S": -2,"T": -1,"W":  1,"Y":  3,"V":  0,"B": -3,"J":  1,"Z": -3,"X": -1}
		self.blosum45["P"]={ "A":-1 ,"R":-2 , "N":-2,"D": -1,"C": -4,"Q": -1,"E":  0,"G": -2,"H": -2,"I": -2,"L": -3,"K": -1,"M": -2,"F": -3,"P":  9,"S": -1,"T": -1,"W": -3,"Y": -3,"V": -3,"B": -2,"J": -3,"Z": -1,"X": -1}
		self.blosum45["S"]={ "A": 1 ,"R":-1 , "N": 1,"D":  0,"C": -1,"Q":  0,"E":  0,"G":  0,"H": -1,"I": -2,"L": -3,"K": -1,"M": -2,"F": -2,"P": -1,"S":  4,"T":  2,"W": -4,"Y": -2,"V": -1,"B":  0,"J": -2,"Z":  0,"X": -1}
		self.blosum45["T"]={ "A": 0 ,"R":-1 , "N": 0,"D": -1,"C": -1,"Q": -1,"E": -1,"G": -2,"H": -2,"I": -1,"L": -1,"K": -1,"M": -1,"F": -1,"P": -1,"S":  2,"T":  5,"W": -3,"Y": -1,"V":  0,"B":  0,"J": -1,"Z": -1,"X": -1}
		self.blosum45["W"]={ "A":-2 ,"R":-2 , "N":-4,"D": -4,"C": -5,"Q": -2,"E": -3,"G": -2,"H": -3,"I": -2,"L": -2,"K": -2,"M": -2,"F":  1,"P": -3,"S": -4,"T": -3,"W": 15,"Y":  3,"V": -3,"B": -4,"J": -2,"Z": -2,"X": -1}
		self.blosum45["Y"]={ "A":-2 ,"R":-1 , "N":-2,"D": -2,"C": -3,"Q": -1,"E": -2,"G": -3,"H":  2,"I":  0,"L":  0,"K": -1,"M":  0,"F":  3,"P": -3,"S": -2,"T": -1,"W":  3,"Y":  8,"V": -1,"B": -2,"J":  0,"Z": -2,"X": -1}
		self.blosum45["V"]={ "A": 0 ,"R":-2 , "N":-3,"D": -3,"C": -1,"Q": -3,"E": -3,"G": -3,"H": -3,"I":  3,"L":  1,"K": -2,"M":  1,"F":  0,"P": -3,"S": -1,"T":  0,"W": -3,"Y": -1,"V":  5,"B": -3,"J":  2,"Z": -3,"X": -1}
		self.blosum45["B"]={ "A":-1 ,"R":-1 , "N": 5,"D":  6,"C": -2,"Q":  0,"E":  1,"G": -1,"H":  0,"I": -3,"L": -3,"K":  0,"M": -2,"F": -3,"P": -2,"S":  0,"T":  0,"W": -4,"Y": -2,"V": -3,"B":  5,"J": -3,"Z":  1,"X": -1}
		self.blosum45["J"]={ "A":-1 ,"R":-3 , "N":-3,"D": -3,"C": -2,"Q": -2,"E": -3,"G": -4,"H": -2,"I":  4,"L":  4,"K": -3,"M":  2,"F":  1,"P": -3,"S": -2,"T": -1,"W": -2,"Y":  0,"V":  2,"B": -3,"J":  4,"Z": -2,"X": -1}
		self.blosum45["Z"]={ "A":-1 ,"R": 1 , "N": 0,"D":  1,"C": -3,"Q":  4,"E":  5,"G": -2,"H":  0,"I": -3,"L": -2,"K":  1,"M": -1,"F": -3,"P": -1,"S":  0,"T": -1,"W": -2,"Y": -2,"V": -3,"B":  1,"J": -2,"Z":  5,"X": -1}
		self.blosum45["X"]={ "A":-1 ,"R":-1, "N": -1,"D": -1,"C": -1,"Q": -1,"E": -1,"G": -1,"H": -1,"I": -1,"L": -1,"K": -1,"M": -1,"F": -1,"P": -1,"S": -1,"T": -1,"W": -1,"Y": -1,"V": -1,"B": -1,"J": -1,"Z": -1,"X": -1}

	def costFunction(self, i, j):
		"""
		Returns the cost of substituting j for i,
		or the cost of a gap if either is a gap
		"""
		if i=="-" or j=="-":
			return -5;
		return self.blosum45[i][j];


	def needlemanwunschAlign(self, x,y,s):
		"""
		Aligns the sequences using Needleman-Wunsch
		(first part of code is adapted from lecture
		2013_2_27 Slide 12: "Global alignment: implementation")
		"""
		self._gaps1 = [];
		self._gaps2 = [];
		D = numpy.zeros((len(x) + 1, len(y) + 1), dtype=int)
		for j in xrange(1, len(y) + 1):
			D[0,j] =  j * s('-', y[j-1]);
		for i in xrange(1, len(x) + 1):
			D[i,0] =  i * s('-', x[i-1]);
		for i in xrange(1, len(x) + 1):
			for j in xrange(1, len(y) + 1):
				D[i,j] = max(D[i-1, j-1] + s(x[i-1], y[j-1]),
					D[i-1,j] + s(x[i-1], '-'),
					D[i,j-1] + s('-', y[j-1]))
		alignToReturnX = "";
		alignToReturnY = "";
		i = len(x)-1;
		j = len(y) - 1;
		# print x
		# print y
		alignToReturnX = x[i]
		alignToReturnY = y[j]
		while i>-1 and j>-1:
			if i==0 and j >0:
				j -= 1;
				self._gaps1.append(i)
				alignToReturnX = "-" + alignToReturnX;
				alignToReturnY = y[j] + alignToReturnY;
			elif j==0 and i > 0:
				i -= 1;
				self._gaps2.append(j)
				alignToReturnY = "-" + alignToReturnY;
				alignToReturnX = x[i] + alignToReturnX;
			else:
				topVal = D[i,j+1];
				diaVal = D[i,j];
				lefVal = D[i+1,j];
				if topVal > diaVal and topVal > lefVal:
					i -= 1;
					self._gaps2.append(j)
					alignToReturnX = x[i] + alignToReturnX;
					alignToReturnY = "-" + alignToReturnY;
				elif diaVal > lefVal:
					i -= 1;
					j -= 1;
					if j>-1 and i>-1:
						alignToReturnX = x[i] + alignToReturnX;
						alignToReturnY = y[j] + alignToReturnY;
				else:
					j -= 1;
					self._gaps1.append(i)
					alignToReturnX = "-" + alignToReturnX;
					alignToReturnY = y[j] + alignToReturnY;
		return (alignToReturnX, alignToReturnY, D[len(x), len(y)])

	def getNonZeroScore(self):
		"""
		Finds the first score that is non-zero
		"""
		for i in xrange(len(self._distMatrix)):
			for j in xrange(i + 1, len(self._distMatrix)):
				if self._distMatrix[i,j] != 0:
					return (i,j)
		return None

	def alignSeqs(self, seqs):
		"""
		Aligns the passed sequences and the returns
		them in their aligned format
		"""
		self._seqs = seqs;
		self.buildDistanceMatrix();
		self._matrixVals = {};
		self._aligned = {};
		self._alignedRes = {};
		self._ignore = [];
		self._seqsToAlign = len(seqs)
		for i in xrange(len(self._seqs)):
			self._matrixVals[i] = [i];
		while self._seqsToAlign>1:
			(i,j) = self.findMinScore();
			if i==-1 or j==-1:
				res = self.getNonZeroScore();
				# see if small sequences
				sameLength = True;
				for i in xrange(len(self._seqs) - 1):
					if len(self._seqs[i])!=len(self._seqs[i + 1]):
						sameLength = False;
				if sameLength and len(self._seqs[0]) < 3:
					return seqs;
				if res==None:
					print "error: cannot align sequences due to disimilarity"
					return None;
				(i,j) = res
			self.merge(i,j);
		self._seqs = None;
		return self._alignedRes.values()


	def merge(self, i, j):
		"""
		Merges two sequences together
		"""
		alignedI = "";
		alignedJ = "";
		if i in self._alignedRes and j in self._alignedRes:
			(alignedI, alignedJ, score) = self.needlemanwunschAlign(self._alignedRes[i], self._alignedRes[j], self.costFunction);
			toUpdate = self._aligned[i];
			for v in toUpdate:
				s = self._alignedRes[v]
				for g in self._gaps1:
					s=s[0:g] + "-" + s[g:]
				self._alignedRes[v] = s;
			toUpdate = self._aligned[j]
			for v in toUpdate:
				s = self._alignedRes[v]
				for g in self._gaps2:
					s=s[0:g] + "-" + s[g:]
				self._alignedRes[v] = s;
			self._alignedRes[i] = alignedI
			self._alignedRes[j] = alignedJ
			self._aligned[i] = list(set(self._aligned[i]) | set(self._aligned[j]) | set([j]))
			self._aligned[j] = list(set(self._aligned[i]) | set(self._aligned[j]) | set([i]))

		elif i in self._alignedRes:
			(alignedI, alignedJ, score) = self.needlemanwunschAlign(self._alignedRes[i], self._seqs[j], self.costFunction);
			toUpdate = self._aligned[i];
			for v in toUpdate:
				s = self._alignedRes[v]
				for g in self._gaps1:
					s=s[0:g] + "-" + s[g:]
				self._alignedRes[v] = s;
			self._alignedRes[j] = alignedJ
			self._alignedRes[i] = alignedI
			self._aligned[i] = list(set(self._aligned[i]) | set([j]));
			self._aligned[j] = list(set(self._aligned[i]) | set([i]));
		elif j in self._alignedRes:
			(alignedJ, alignedI, score) = self.needlemanwunschAlign(self._alignedRes[j],self._seqs[i], self.costFunction);
			toUpdate = self._aligned[j];
			for v in toUpdate:
				s = self._alignedRes[v]
				for g in self._gaps1:
					s=s[0:g] + "-" + s[g:]
				self._alignedRes[v] = s;
			self._alignedRes[i] = alignedI
			self._alignedRes[j] = alignedJ
			self._aligned[j] = list(set(self._aligned[j]) | set([i]));
			self._aligned[i] = list(set(self._aligned[j]) | set([j]));
		else:
			(alignedI, alignedJ, score) = self.needlemanwunschAlign(self._seqs[i], self._seqs[j], self.costFunction);
			self._alignedRes[i] = alignedI;
			self._alignedRes[j] = alignedJ;
			self._aligned[j] = [i]
			self._aligned[i] = [j]
		# update distance matrix
		self._ignore.append(i);
		for k in xrange(len(self._distMatrix)):
			for l in xrange(len(self._distMatrix)):
				if k==i or l==i or (k in self._ignore and k!= j) or (l in self._ignore and l!=j):
					self._distMatrix[k,l] = 0;
				elif k==j and k!=l:
					self._distMatrix[k,l] = (self._distMatrixOriginal[k-1,i] + self._distMatrixOriginal[k,l])/2.0
				elif l==j and k!=l:
					self._distMatrix[k,l] = (self._distMatrixOriginal[i,l-1] + self._distMatrixOriginal[k,l])/2.0
		self._seqsToAlign -= 1;


	def findMinScore(self):
		"""
		This actually finds the maximum score in the matrix
		"""
		currBest = 0;
		currBextJ = -1;
		currBestI = -1;
		for i in xrange(len(self._distMatrix)):
			for j in xrange(i + 1, len(self._distMatrix)):
				if self._distMatrix[i,j] > currBest:
					currBest = self._distMatrix[i,j];
					currBextJ = j;
					currBestI = i;
		return (currBestI, currBextJ)

	def buildDistanceMatrix(self):
		"""
		Builds the distance matrix for all of the sequences
		"""
		overall = numpy.zeros((len(self._seqs), len(self._seqs)), dtype=float);
		overall2 = numpy.zeros((len(self._seqs), len(self._seqs)), dtype=float);
		self._aligns = {};
		for i in xrange(len(self._seqs)):
			for j in xrange(len(self._seqs)):
				score = 0;
				if i!=j:
					(alI, alJ, score) = self.needlemanwunschAlign(self._seqs[i], self._seqs[j], self.costFunction)
					self._aligns[(i,j)] = (alI, alJ);
				overall[i,j] = score;
				overall2[i,j] = score;
		self._distMatrix = overall;
		self._distMatrixOriginal = overall2;


if __name__=="__main__":
	a = JunkAligner();
	seqs = ["TACGTCAGC", "TATGTCATGC"]
	a.alignSeqs(seqs)
	seqs = ["GCTATGCGGCTATACGC", "GCGTATGCGGCTAACGC"]
	a.alignSeqs(seqs)
	seqs = ["SCAVPSTDDY", "SSVPSKQTY"]
	a.alignSeqs(seqs)
	seqs = ["SCAVPSTDDY", "SSSVKTY"]
	a.alignSeqs(seqs)
	seqs = ["SCAVPSTDDY", "SSSVKTY", "SSVPSKQTY"]
	a.alignSeqs(seqs)
	seqs = ["TFYFPVLPPM", "TFYPVVLPPM", "LFYFVVPPPM", "TMYFVIPPPT"]
	print seqs
	vals = a.alignSeqs(seqs)
	for vl in vals:
		print vl
	seqs = ["TTFTFLFL", "TTFTLLLLL", "-"]
	vals = a.alignSeqs(seqs)
	for vl in vals:
		print vl







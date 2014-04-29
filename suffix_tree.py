

import fasta as fas
import sys
import subprocess
import time
import junk_aligner


class Delimiter:
	"""
	Class that was entended for more use than jsut generating the
	delimiters but it ended up being used for that one purpose
	"""
	def __init__(self):
		"""
		Loads all delimiters into object
		"""
		self._dels = ["~", "`", "!", "@", "#", "$", "%", "^", "&", "*",
					  "(", ")", "{", "[", "}", "]", ";", ":", "<", ">"];
		self._index = 0;

	def nextDelimiter(self):
		"""
		Returns the next delimeter sequence (currently
		only up to 20 can be aligned)
		"""
		if self._index > 19:
			raise NoMoreDelimitersException;
		toReturn = self._dels[self._index];
		self._index += 1;
		return toReturn;

class SuffixNode:
	"""
	Not even sure if this is used any more?
	"""
	def __init__(self, start, end):
		"""
		Stored the start and endind points of a node
		"""
		self._start = start;
		self._end   = end;

	def getStart(self):
		"""
		Returns the position the node starts at
		"""
		return self._start;
	def getEnd(self): 
		"""
		Returns the position the node ends at
		"""
		return self._end;
	def setEnd(self, newEnd):
		"""
		Sets the position the node ends at
		"""
		self._end = newEnd;

class Edge:
	"""
	The bulk of the data is stored in these Edges which represent
	the branches in the tree
	"""
	def __init__(self, start, end, startNode, endNode, num, is_external):
		"""
		Sets all the edges data members to the passed values
		"""
		self._is_external = is_external;
		self._start = start;
		self._end   = end;
		self._startNode = startNode;
		self._endNode   = endNode;
		self._num = num;

	def getEnd(self):
		"""
		Returns the end position of the edge in the
		suffix tree's sequence
		"""
		if self._is_external:
			return self._end._end;
		return self._end;

	def __repr__(self):
		"""
		Retuns a string representation of the suffix tree
		"""
		return "%-10s%-10s%-10s%-10s%-10s" % (self._start, self.getEnd(), self._startNode, self._endNode, self._num)

class DeepNodes:
	"""
	Used to keep track of internal nodes in 
	the tree so that common motif lookup can 
	be done much faster than actually 
	traversing the tree
	"""
	def __init__(self):
		"""
		Initiliazes memory for the container
		"""
		self._maxCapacity = 15;
		self._nodes=[];
		self._nodeLocation={};
	def insertNode(self, node, height):
		"""
		Inserts a node with the given height into the container
		Nodes are sorted by height so that we can grab the longest
		edges much more quickly
		"""
		# if we're empty
		if len(self._nodes)==0:
			self._nodes.append((node, height));
			self._nodeLocation[node] = 0;
			return;
		# if it would be the smallest possible element
		(lastNode, lastHeight) = self._nodes[len(self._nodes) - 1];
		if lastHeight >= height:
			self._nodes.append((node, height));
			self._nodeLocation[node] = len(self._nodes) - 1;
			return;
		# see if it is as big as biggest
		(firstNode, firstHeight) = self._nodes[0];
		if firstHeight <= height:
			self._nodes[1:len(self._nodes)] = self._nodes;
			self._nodes[0] = (node, height);
			self._nodeLocation[node] = -1
			updateDict = lambda dic : [(k, v+1) for (k,v) in dic.iteritems()];
			revertListtoDict = lambda lis : dict(lis);
			self._nodeLocation = revertListtoDict(updateDict(self._nodeLocation))
			return;
		#try add internally
		for i in xrange(len(self._nodes) - 1, -1, -1):
			(currNode, currHeight) = self._nodes[i];
			if height < currHeight:
				self.insertNodeAt(i + 1, node, height);
				return;
		return;

	def insertNodeAt(self, i, node, height):
		"""
		Inserts a node at a specific location in the container
		"""
		oldLen = len(self._nodes);
		# if oldLen<self._maxCapacity:
		self._nodes.append(self._nodes[oldLen - 1]);
		oldLen += 1;
		for j in xrange(len(self._nodes)-1, i, -1):
			self._nodes[j] = self._nodes[j - 1];
			self._nodeLocation[self._nodes[j]] = j - 1;
		self._nodes[i] = (node, height);
		self._nodeLocation[node] = i;
		return;

	def updateNode(self, node, newHeight):
		"""
		Updates a node to have a new height
		"""
		if not node in self._nodeLocation:
			return; # and throw exception?
		index = self._nodeLocation[node];
		# would be faster to update intelligently (shift somewhere to left, but right
		# now I just want something that works)
		oldVal = self._nodes[index];
		(oldNode, oldHeight) = oldVal;
		if oldHeight > newHeight:
			# print "this means a bug happened"
			return;
		self._nodes.remove(oldVal);
		self._nodeLocation.pop(node)
		self.insertNode(node, newHeight);

	def getShallowest(self):
		"""
		Returns the shallowest node
		(was once used to maintain the container at a specific size)
		"""
		if len(self._nodes) > 0:
			return self._nodes[len(self._nodes) - 1];
		return None;

	def getDeepest(self):
		"""
		Returns the deepest node in the container/tree
		"""
		if len(self._nodes) > 0:
			return self._nodes[0];
		return None;

	def getNodes(self):
		"""
		Returns all internal nodes
		"""
		return self._nodes;

	def __repr__(self):
		"""
		Returns a string representation of the container
		"""
		toReturn = "";
		for (node, height) in self._nodes:
			toReturn += ( "Node #%d ---> height: %d\n" %(node._num, height) );
		return toReturn;

class Node:
	"""
	Used to represent a node in the suffix tree
	to mediate the transition between edges
	"""
	def __init__(self, num, parent):
		"""
		Builds a node based on the given data
		"""
		self._suffixLinks = [];
		self._num = num;
		self._parent = parent;
		self._parentEdge = None;
	def setParentEdge(self, edge):
		"""
		Sets the parents edge to the new value
		"""
		self._parentEdge = edge;
	def getParentEdgeLength(self):
		"""
		Returns a node's parent edge
		"""
		if self._parentEdge == None:
			return 0;
		return self._parentEdge.getEnd() + 1 - self._parentEdge._start;

class ActivePoint:
	"""
	A class for holding onto the current edge and node we are at
	in the tree and for also holding onto where in the sequence we
	are
	"""
	def __init__(self, node, edge, first_index, end_index):
		"""
		Constructs the ActivePoint with the given data
		"""
		self._node = node;
		self._edge = edge;
		self._first_index = first_index
		self._real_first_index = first_index;
		self._end_index = end_index;
	def length(self):
		"""
		Returns the current length of the ActivePoint
		"""
		return self._end_index - self._first_index;

class EndPoint:
	"""
	Used to continually update the end index of all leaf nodes
	"""
	def __init__(self, end):
		"""
		Initializes the EndPoint to the specified index
		"""
		self._end = end;

	def incrementEnd(self):
		"""
		Increments the end index by one
		"""
		self._end +=1;

class SuffixTree:
	"""
	The big kahuna, this class represents the suffix tree
	and holds all the functionality to build itself
	"""
	def __init__(self, seq):
		"""
		Initializes all data members of the tree in order
		for it to build itself
		"""
		self._seq = seq;
		self._nodes=[];
		self._nodes.append(Node(0, -1));
		self._edgeList = [];
		self._edges={};
		self._root = SuffixNode(-1, -1);
		self._active = ActivePoint(0, -1, 0, 0);
		self._remainder = 0;
		self._edge_num = 0;
		self._end_point = EndPoint(-1);
		self._deepNodes = DeepNodes();
		self._nodeChildren = {};
		self._nodeChildren[0] = [];
		self.build();

	def getNodeDepth(self, node):
		"""
		Returns the depth of a given node
		"""
		parentEdgeLength = node.getParentEdgeLength();
		if node._parent == 0:
			return parentEdgeLength;
		parentNode = self._nodes[node._parent];
		return self.getNodeDepth(parentNode) + parentEdgeLength;

	def build(self):
		"""
		Coordinates the main building/extension process
		"""
		for i in xrange(len(self._seq)):
			self._remainder += 1;
			self._end_point.incrementEnd();
			self._active._end_index = i;
			self.addChar(i);

	def addChar(self, i):
		"""
		Extends all suffixes by one letter
		"""
		last_parent_node = -1;
		while True:
			parent_node = self._active._node;
			if self._remainder == 0 :
				return;
			if self._active._end_index < self._active._first_index:
				return;
			if self._active.length()==0 and not (self._active._node, self._seq[i]) in self._edges and last_parent_node==-1:
				edge = Edge(self._active._first_index, self._end_point, self._active._node, -1, self._edge_num, True );
				self._nodeChildren[self._active._node].append(edge._num);
				self._edges[(edge._startNode, self._seq[edge._start])] = edge;
				self._active._first_index += 1;
				self._active._real_first_index += 1;
				self._edge_num += 1;
				self._edgeList.append(edge);
				self._remainder -= 1;
				self._active._node = 0;
				self._active._edge = -1
				self._active._first_index = self._active._end_index + 1 - self._remainder;
				continue;
			# if we have to split a node, then come back and the letter isn't in the tree and active.length()==0
			if self._active.length()==0 and not (self._active._node, self._seq[i]) in self._edges:
				edge = Edge(self._active._first_index, self._end_point, self._active._node, -1, self._edge_num, True );
				self._nodeChildren[self._active._node].append(edge._num);
				self._edges[(edge._startNode, self._seq[edge._start])] = edge;
				self._active._first_index += 1;
				self._active._real_first_index += 1;
				self._edge_num += 1;
				self._edgeList.append(edge);
				self._remainder -= 1;
				if self._active._node != last_parent_node and self._active._node != 0:
					self._nodes[self._active._node]._suffixLinks.append(last_parent_node);
				self._active._node = 0;
				self._active._edge = -1
				self._active._first_index = self._active._end_index + 1 - self._remainder;
				continue;
			# we need to deal with case when we've added an edge when active length > 0, but now active_edge = -1
			if self._active._edge < 0 and self._active.length() > 0:
				if (self._active._node, self._seq[self._active._first_index]) in self._edges:
					edge = self._edges[(self._active._node, self._seq[self._active._first_index])];
					self._nodeChildren[self._active._node].append(edge._num);
					self._active._edge = edge._num;
					# need to try to resync jumping here too
					currEdge = self._edgeList[self._active._edge];
					currEdgeLength = currEdge.getEnd() - currEdge._start;
					if not self._edgeList[self._active._edge]._is_external and \
							self._seq[currEdge._start : currEdge.getEnd() + 1]==self._seq[self._active._first_index  : self._active._first_index + currEdgeLength + 1]:
						currEdge = self._edgeList[self._active._edge];
						currEdgeLength = currEdge.getEnd() - currEdge._start;
						currEndNode = self._edgeList[self._active._edge]._endNode;
						if (currEndNode, self._seq[self._active._first_index + currEdgeLength + 1]) in self._edges:
							self._active._node = currEndNode;
							self._active._edge = self._edges[(currEndNode, self._seq[self._active._first_index + currEdgeLength + 1])]._num
							self._active._first_index = self._active._first_index + currEdgeLength + 1 ;
					continue;
				else:
					edge = Edge(self._active._first_index, self._end_point, self._active._node, -1, self._edge_num, True );
					self._nodeChildren[self._active._node].append(edge._num);
					self._edges[(edge._startNode, self._seq[edge._start])] = edge;
					self._active._first_index += 1;
					self._active._real_first_index += 1;
					self._edge_num += 1;
					self._edgeList.append(edge);
					self._remainder -= 1;
					self._active._edge = -1;
					self._active._node = 0;
					continue;
			# if we get here either our char exists in a single edge, or length() > 0
			if self._active._edge > -1:
				# first we need to see if we're at the end of the current edge
				if not self._edgeList[self._active._edge]._is_external:
					# so if we're in an internal node we need to see if we've jumped from one edge to a child edge
					currEndNode = self._edgeList[self._active._edge]._endNode;
					currEdge = self._edgeList[self._active._edge];
					currEdgeLength = currEdge.getEnd() - currEdge._start;
					if (currEndNode, self._seq[i]) in self._edges and \
						self._seq[currEdge._start : currEdge.getEnd() + 1]==self._seq[self._active._end_index - currEdgeLength - 1 : self._active._end_index ]:
						# we've jumped, but don't we need to make sure last we traversed the entirety of the last edge
						self._active._node = currEndNode;
						self._active._edge = self._edges[(currEndNode, self._seq[i])]._num
						self._active._first_index += currEdgeLength  + 1;
						continue;
				edge = self._edgeList[self._active._edge];
				if self._active.length() >  edge.getEnd() - edge._start:
					self._active._node = edge._endNode;
					self._active._first_index += (edge.getEnd() - edge._start) + 1;
					self._active._edge = -1;
					continue;
				s = self._edgeList[self._active._edge]._start
				e = self._edgeList[self._active._edge].getEnd();
				# we're currently in an edge
				if self._seq[self._edgeList[self._active._edge]._start + self._active.length()]==self._seq[i]:
					#edge exists
					return;
				else:
					# time to split
					if last_parent_node == -1:
						node = Node(len(self._nodes) , self._active._node);
						self._nodes.append(node);
						node_index = len(self._nodes) - 1;
						old_edge = self._edgeList[self._active._edge];
						node.setParentEdge(old_edge);
						old_end = old_edge._endNode;
						e1 = Edge(old_edge._start + self._active.length(), old_edge._end, node_index, old_edge._endNode, self._edge_num, old_edge._is_external);
						self._edge_num += 1;
						e2 = Edge(self._active._end_index, self._end_point, node_index, -1, self._edge_num, True);
						self._edge_num += 1;
						old_edge._endNode = node_index;
						old_edge._end = self._edgeList[self._active._edge]._start + self._active.length()-1
						old_edge._is_external = False;
						self._nodeChildren[node._num] = [e1._num, e2._num];
						if old_end > -1:
							self._nodes[old_end].setParentEdge(e1);
							self._nodes[old_end]._parent = node_index;
							old_height = self.getNodeDepth(self._nodes[old_end]);
							self._deepNodes.updateNode(self._nodes[old_end], old_height);
						height = self.getNodeDepth(node);
						self._deepNodes.insertNode(node, height);
						self._edges[(node_index, self._seq[e1._start])] = e1;
						self._edges[(node_index, self._seq[e2._start])] = e2;
						self._edgeList.append(e1);
						self._edgeList.append(e2);
						last_parent_node = node_index;
						self._active._first_index += 1;
						self._remainder -= 1;
						if len(self._nodes[node_index]._suffixLinks) > 0:
							self._active._edge = -1;
							self._active._node = self._nodes[self._active._node]._suffixLinks[0];
						else:
							self._active._edge = -1;
							self._active._node = 0;
							self._active._first_index = self._active._end_index + 1 - self._remainder;
					else:
						node = Node(len(self._nodes), self._active._node);
						self._nodes.append(node);
						node_index = len(self._nodes) - 1;
						old_edge = self._edgeList[self._active._edge];
						node.setParentEdge(old_edge);
						old_end = old_edge._endNode;
						e1 = Edge(old_edge._start + self._active.length(), old_edge._end, node_index, old_edge._endNode, self._edge_num, old_edge._is_external);
						self._edge_num += 1;
						e2 = Edge(self._active._end_index, self._end_point, node_index, -1, self._edge_num, True);
						self._edge_num += 1;
						old_edge._endNode = node_index;
						old_edge._end = self._edgeList[self._active._edge]._start + self._active.length()-1
						old_edge._is_external = False;
						self._nodeChildren[node._num] = [e1._num, e2._num];
						if old_end > -1:
							self._nodes[old_end].setParentEdge(e1);
							self._nodes[old_end]._parent = node_index;
							old_height = self.getNodeDepth(self._nodes[old_end]);
							self._deepNodes.updateNode(self._nodes[old_end], old_height);
						height = self.getNodeDepth(node);
						self._deepNodes.insertNode(node, height);
						self._edges[(node_index, self._seq[e1._start])] = e1;
						self._edges[(node_index, self._seq[e2._start])] = e2;
						self._edgeList.append(e1);
						self._edgeList.append(e2);
						if node_index != last_parent_node:
							self._nodes[last_parent_node]._suffixLinks.append(node_index);
						last_parent_node = node_index;
						self._active._first_index += 1;
						self._remainder -= 1;
						if len(self._nodes[node_index]._suffixLinks) > 0:
							self._active._edge = -1;
							self._active._node = self._nodes[self._active._node]._suffixLinks[0];
						else:
							self._active._edge = -1;
							self._active._node = 0;
							self._active._first_index = self._active._end_index + 1 - self._remainder;
						continue
			# we can't be in an edge and failed to split or build on edge so move along edge
			elif (self._active._node, self._seq[i]) in self._edges:
				edge = self._edges[(self._active._node, self._seq[i])];
				self._active._edge = edge._num;
				return;
	


	def splitEdge(self, edge_index, i):
		"""
		Splits the specified edge at the specified index
		(no longer used, edge splitting code is embedded
		in addChar)
		"""
		node = Node(self._edges[edge_index]._startNode);
		self._nodes.append(node);
		node_index = len(self._nodes) - 1;
		old_edge = self._edges[edge_index];
		e1 = Edge(old_edge._start + self._active.length(), old_edge._end, node_index, old_edge._endNode, self._edge_num, old_edge._is_external);
		self._edge_num += 1;
		e2 = Edge(i - self._active.length(), self._end_point, node_index, -1, self._edge_num, True);
		self._edge_num += 1;
		old_edge._endNode = node_index;
		self._edges[(node_index, self._seq[e1._start])] = e1;
		self._edges[(node_index, self._seq[e2._start])] = e2;
		self._active._edge = -1;
		self._active._first_index += 1;
		self._remainder -= 1;
		return node_index;

	def toDotFormat(self):
		"""
		Returns the tree in dot format so that it can be visualized
		"""
		toReturn = "digraph X {\n";
		nodes = self._nodes;
		for node in range(len(nodes)):
			toReturn += "node_" + str(node) + " [label=\"" + str(node) + "\"];\n";
		values = self._edges.values();
		values = sorted(values, key=lambda x:x._startNode);
		for value in values:
			toReturn += "\tnode_" + str(value._startNode) + " -> node_";
			if value._endNode == -1:
				toReturn += str(value._startNode) + str(value._num) + "_end";
			else:
				toReturn += str(value._endNode);
			toReturn += " [label=\"" + self._seq[value._start:value.getEnd()+1] + "\"];\n"
		for n in range(len(nodes)):
			node = nodes[n];
			if len(node._suffixLinks) > 0:
				for link in node._suffixLinks:
					toReturn += "\tnode_" + str(n) + " -> node_" + str(link) + " [style=\"dashed\"];\n";
		toReturn += "}\n";
		return toReturn;

	def getSequence(self, node):
		"""
		Returns the sequence which ends at the given node
		"""
		baseString = self._seq[node._parentEdge._start : node._parentEdge._start + node.getParentEdgeLength()];
		if node._parent == 0:
			return baseString;
		parent = self._nodes[node._parent];
		return self.getSequence(parent) + baseString;

	def getNodeChildren(self, node):
		"""
		Returns all the children of a node
		"""
		toReturn = 0;
		for eNum in self._nodeChildren[node._num]:
			edge = self._edgeList[eNum];
			if edge._endNode != -1:
				toReturn += self.getNodeChildren(self._nodes[edge._endNode]);
			else:
				toReturn += 1;
		return toReturn;

	def getMostCommonSequences(self):
		"""
		Returns a string representation of the most common
		sequences (for testing/debugging)
		"""
		toReturn = "";
		deepNodes = self._deepNodes.getNodes();
		for i in range(len(deepNodes)):
			(node, height) = deepNodes[i];
			commSeq = self.getSequence(node);
			toReturn += str(i) + " node #" + str(node._num) + " (" + str(self.getNodeChildren(node))  +" chillins) ---> " + commSeq + " ---> ending at indices [";
			# for eNum in self._nodeChildren[node._num]:
			# 	toReturn += str(self._edgeList[eNum]._start) + ", "
			# indices = self.getIndices(node, self._edgeList[node._parentEdge._num].getEnd() + 1);
			indices = self.getIndices(node, height);
			for seqIndex in indices:
				toReturn += "(" + str(seqIndex) + ", " + str(seqIndex + len(commSeq)) + "), "
			toReturn += "]\n"
		return toReturn;

	def getIndices(self, node, currDepth):
		"""
		Returns all the indices of where the paths start which end at the specified node
		"""
		toReturn = [];
		for eNum in self._nodeChildren[node._num]:
			edge = self._edgeList[eNum];
			if edge._endNode != -1:
				toAdd = self.getIndices(self._nodes[edge._endNode], currDepth + edge.getEnd()  + 1 - edge._start);
				toReturn[len(toReturn):len(toReturn) + len(toAdd)] = toAdd;
			else:
				toReturn.append(edge._start - currDepth);
		return toReturn;

	def bestKSeqSeeds(self, k):
		"""
		Returns all common seeds which occur at least k times
		"""
		toReturn = [];
		deepNodes = self._deepNodes.getNodes();
		for i in range(len(deepNodes)):
			(node, height) = deepNodes[i];
			numChildren = 0
			seqIndices = [];
			indices = self.getIndices(node, height);
			commSeq = self.getSequence(node);
			if len(indices) == k:
				toReturn.append((node, len(commSeq) , indices, self.getNodeChildren(node)));
		return toReturn;



	def __repr__(self):
		"""
		Returns a string representation of the tree for printing out
		"""
		toReturn = ""
		values = self._edges.values();
		values = sorted(values, key=lambda x:x._startNode);
		toReturn += "%-10s%-10s%-10s%-10s%-10s%-10s\n" % ("Start", "End", "S-Node", "E-Node", "ID", "Edge");
		for val in values:
			toReturn+= str(val) + "" + self._seq[val._start:val.getEnd()+1] + "\n"
		return toReturn
			
	def align(self, numSeqs):
		"""
		Aligns the sequences from which the tree was made
		"""
		junker = junk_aligner.JunkAligner();
		seeds = self.bestKSeqSeeds(numSeqs)
		symbols = ["~", "`", "!", "@", "#", "$", "%", "^", "&", "*",
					  "(", ")", "{", "[", "}", "]", ";", ":", "<", ">"];
		if len(seeds)==1:
			# is this case finished???
			forClustal = [];
			toAdd = [];
			numWithLengthZero = 0
			for i in xrange(len(self._seq)):
				if self._seq[i] in symbols:
					numWithLengthZero += 1
					forClustal.append(''.join(toAdd));
					toAdd = [];
				else:
					toAdd.append(self._seq[i])
			if numWithLengthZero > 0:
				print self._seq
				return ['', '', '']
			return junker.alignSeqs(forClustal);
			return seqs
		bestTwo = bestTwoSeeds(seeds);
		if bestTwo==None or len(bestTwo)==0:
			# no seeds so just align with junk aligner
			forClustal = [];
			toAdd = [];
			numWithLengthZero = 0;
			for i in xrange(len(self._seq)):
				if self._seq[i] in symbols:
					if len(toAdd)==0:
						numWithLengthZero += 1;
						# print "got here"
					forClustal.append(''.join(toAdd));
					toAdd = [];
				else:
					toAdd.append(self._seq[i])
			return junker.alignSeqs(forClustal);
		else:
			bestTwo = sorted(bestTwo, key=lambda x:x[0])
			toReturn = [];
			for j in xrange(len(bestTwo[0])):
				toReturn.append("");
			# build alignments from seeds and aligned "junk"/highly variable regions
			for i in xrange(len(bestTwo) - 1):
				seqsToAlign = [];
				validOrder = True;
				for j in xrange(len(bestTwo[i])):
					(s , e ) = bestTwo[i][j];
					(ns, ne) = bestTwo[i+1][j];
					if e>=ns:
						validOrder = False;
				if validOrder:
					for j in xrange(len(bestTwo[i])):
						(s , e ) = bestTwo[i][j];
						(ns, ne) = bestTwo[i+1][j];
						seqsToAlign.append(self._seq[e:ns])
						toReturn[j] = toReturn[j] + self._seq[s:e]
				# add aligned variable region
				alignedJunk = junker.alignSeqs(seqsToAlign);
				for j in xrange(len(alignedJunk)):
					toReturn[j] = toReturn[j] + alignedJunk[j];
			for j in xrange(len(bestTwo[0])):
				(s,e) = bestTwo[-1][j]
				toReturn[j] = toReturn[j] + self._seq[s:e];
			# all alignment between outermost seeds is done
			# now align any regions outside these
			seqStarts=[0];
			seqEnds=[];
			for i in xrange(len(self._seq)):
				if self._seq[i] in symbols:
					seqEnds.append(i);
					if i != len(self._seq)-1:
						seqStarts.append(i+1)
			# get anything to the left of the leftmost seed
			leftSides = [];
			leftExtra = 0;
			for i in xrange(len(bestTwo[0])):
				(s,e) = bestTwo[0][i]
				leftSides.append(self._seq[seqStarts[i]:s])
				if len(leftSides[i]) > leftExtra:
					leftExtra = len(leftSides[i])
			rightSides = [];
			rightExtra = 0;
			for i in xrange(len(bestTwo[0])):
				(s,e) = bestTwo[-1][i]
				rightSides.append(self._seq[e:seqEnds[i]])
				if len(rightSides[i]) > rightExtra:
					rightExtra = len(rightSides[i]);
			if leftExtra==rightExtra==0:
				return toReturn;
			for i in xrange(len(bestTwo[0])):
				if len(leftSides[i])==0:
					leftSides[i] = "-";
			leftalign = junker.alignSeqs(leftSides);
			for i in xrange(len(bestTwo[0])):
				if len(rightSides[i])==0:
					rightSides[i] = "-";
			rightalign = junker.alignSeqs(rightSides);
			for i in xrange(len(bestTwo[0])):
				toReturn[i] = leftalign[i] + toReturn[i] + rightalign[i];
			return toReturn;




def concatSequences(seqs):
	"""
	Concatenates the given sequences into one with delimiters
	between the different sequences
	"""
	if len(seqs) < 1:
		raise EmptySequencesException
	delim = Delimiter();
	currSeq = seqs[0];
	currSeq = [currSeq];
	currSeq = ''.join(currSeq);
	for seq in seqs[1:]:
		currSeq = [currSeq, delim.nextDelimiter(), seq];
		currSeq = ''.join(currSeq);
	currSeq = [currSeq, delim.nextDelimiter()];
	currSeq = ''.join(currSeq);
	return currSeq;



def main():
	"""
	For testing
	"""
	(names, seqs) = fas.parseFastaFile("p150.fa");
	seq = "ABCDEZXCVBPLMOKNFENWAYPARKGHIJKL!ABCDEASDFGPLMOKNVERWSQGORKGHIJKL@ABCDEQWERTYPLMOKNFONDABBAMSGHIJKL#"
	seqs = ["", "", ""]
	names = ["seq1", "seq2", "seq3", "seq4"]
	seqs = ["AABCDETFYFPVLPPMLMKNSSVPSKQTYGHIJKL", "ABCDETFYPVVLPPMLMKNSSVPSKQTYGHIJKL",
			"ABCDELFYFVVPPPMLMKNSSSVKTYGHIJKL", "ABCDETMYFVIPPPTLMKNSCAVPSTDDYGHIJKL"]
	seq = seqs[0] + "~" + seqs[1] + "!" + seqs[2] + "@" + seqs[3] + "#"
	start = time.time()
	tree = SuffixTree(seq)
	aligned = tree.align(len(seqs));
	fas.printToFasta(names, aligned, "g6pd.msa");
	end = time.time() - start
	print end

def bestTwoSeeds(seeds):
	"""
	Incorrectly named (due to change in functionality
	but not from where it was called).
	Returns the best seeds to try to align with
	(will only pick seeds which occur in every sequence
	and are at least of length 3)
	"""
	toReturn = [];
	numAdded = 0;
	# print "====================="
	for (node, seqLen, indices, chillins) in seeds:
		if len(toReturn)==0 and seqLen > 1:
			toAdd = [];
			indices = sorted(indices, key=lambda x:x)
			for index in indices:
				toAdd.append((index, index+seqLen));
			toReturn.append(toAdd);
			numAdded += 1;
		elif seqLen > 1:
			indices = sorted(indices, key=lambda x:x)
			(firstStart, firstEnd) = toReturn[0][0]
			notContained = True;
			for group in toReturn:
				(firstStart, firstEnd) = group[0]
				if indices[0] >= firstStart and indices[0] + seqLen <= firstEnd:
					notContained = False;

			if notContained and seqLen > 2:
				toAdd = [];
				for index in indices:
					toAdd.append((index, index+seqLen));
				toReturn.append(toAdd);
				numAdded += 1;
	return toReturn;
	




if __name__=="__main__":
	main();





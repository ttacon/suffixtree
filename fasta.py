

def parseFastaFile(filename):
    """
    Helper function which sets up for parseFasta
    Returns None if an IOError occurs (i.e. file can't be
    found)
    """
    try:
        f = open(filename, 'r');
        lines = f.readlines();
    except IOError:
        return None;
    return parseFasta(lines, True);

def parseFasta(lines, hasNewLines):
    """
    Parses a FASTA file denoted by the given filename.
    Returns the names and sequences of every entry in the file
    """
    seqs = [];
    names=[];
    if hasNewLines:
        lines = [x[0:len(x)-1] for x in lines];
    currSeq = "";
    lastName = "";
    for i in range(len(lines)):
        line = lines[i];
        if line[0]==">":
            if i > 0:
                seqs.append(currSeq);
                names.append(lastName);
            currSeq = "";
            lastName = line[1:];
        else:
            currSeq += line
    names.append(lastName);
    seqs.append(currSeq);
    return (names, seqs);

def printToFasta(names, seqs, filename):
    """
    Prints the given name and sequence of every entry
    to the file denoted by the given file name.
    Sequences are truncated to the next line at 60
    characters.
    """
    f = open(filename, 'w');
    for i in xrange(len(names)):
        f.write('>');
        f.write(names[i]);
        f.write('\n');
        pos = 0;
        for j in xrange(len(seqs[i])):
            if pos != 0 and pos%60==0:
                f.write('\n');
            f.write(seqs[i][j]);
            pos+=1
        f.write('\n');
    f.close();


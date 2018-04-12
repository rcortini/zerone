import numpy as np

class ZeroneOutput :
    def __init__(self,fname,chromosome_list=None) :
        self.fname = fname
        self.chromosome_list = chromosome_list
        # the first thing we do is parse the output file passed by the user
        self.parse_zerone_output()
        # next, we exclude the values of the array that pertain to chromosomes that are not
        # included in the chromosome list that was passed by the user (if any)
        if chromosome_list is not None :
            self.zerone = np.array([s for s in self.zerone if s['chr'] in chromosome_list])
        # now make the index of the output for quick access to data
        self.make_zerone_output_index()

    def parse_zerone_output(self) :
        """
        Parses a Zerone output and returns a numpy array. The values of the numpy array
        are: chromosome, start, end, enrichment, read_1, read_2, ..., read_n, p.
        The number of `read_i` columns depends on the invocation of Zerone and cannot
        be known beforehand.
        """
        # first, we start by reading the first non-comment line in the Zerone file, to
        # determine the number of `read` columns in the file
        comment = True
        with open(self.fname,'r') as f :
            for line in f :
                if not line.startswith('#') :
                    break
        n_readcols = len(line.split())-6
        zerone_dtype = [('chr','S256'),
                        ('start',np.int64),
                        ('end',np.int64),
                        ('enrichment',np.int32),
                        ('control',np.int64)]
        for i in range(n_readcols) :
            zerone_dtype.append(('read_%d'%(i+1),np.int64))
        zerone_dtype.append(('p',float))
        # now we parse the file using the `genfromtxt` function from numpy
        self.zerone = np.genfromtxt(self.fname,dtype=np.dtype(zerone_dtype))

    def make_zerone_output_index(self) :
        # if the chromosome names are not given, get them
        if self.chromosome_list is None :
            chromosome_list = np.unique(self.zerone['chr'])
        # init the dictionary and init the iteration
        c_idx = {}
        c_start = 0
        prev_c = self.zerone[0]['chr']
        # we don't sort the input array: Zerone already does this by default
        for i,b in enumerate(self.zerone) :
            this_c = b['chr']
            if this_c != prev_c :
                c_end = i-1
                c_idx[prev_c] = (c_start,c_end)
                c_start = i
                prev_c = this_c
        # the last chromosome needs to be manually added
        c_idx[this_c] = (c_start,i)
        self.c_idx = c_idx

    def find_peak(self,chromosome,start,end,bin_size=300) :
        """
        Returns the values of the `zerone` array corresponding to the genomic coordinates
        of the `peak`. Uses the `c_idx` dictionary to rapidly calculate which are the indices
        of the `a` array that correspond to the peak
        """
        c_start,c_end = self.c_idx[chromosome]
        peak_idx_start = start//bin_size
        peak_idx_end = end//bin_size
        if peak_idx_start == peak_idx_end :
            return [self.zerone[c_start+peak_idx_start]]
        else :
            return self.zerone[c_start+peak_idx_start:c_start+peak_idx_end+1]

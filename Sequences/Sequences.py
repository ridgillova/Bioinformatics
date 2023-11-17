class Sequences:
    # common constructor
    def __init__(self, id, seq):
        self._id = id
        self._seq = seq.upper()

    @ property
    def id(self):
        return self._id
    
    @ property
    def seq(self):
        return self._seq
    
    def write_fasta(self):
        fasta = f'>{self.id}\n{self.seq}\n'
        return fasta
    
    # define what to fo when ask len(object)
    def __len__(self):
        return len(self.seq)
    
    # define what to do when ask to print(object)
    def __str__(self):
        return f'ID: {self.id}, Sequence: {self.seq}'
    
    # define what to do when try to iterate over object
#    def __iter__(self):
#        self.length = 0
#        return self
    
#    def __next__(self):
#        if self.length == len(self.seq):
#            raise StopIteration
#        self.length += 1
#        return self.seq
        
    
class DNASequence(Sequences):

    # attributes and methods inherited from superclass Sequence

    # gc content method
    def calc_gc_content(self, dp=2):
        c_count = self.seq.count('C')
        g_count = self.seq.count('G')
        gc_content = (c_count + g_count) / len(self.seq)
        return round(gc_content, dp)

    # translate sequence
    def translate(self):
        bases = "tcag".upper()
        codons = [a + b + c for a in bases for b in bases for c in bases]
        amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
        codon_table = dict(zip(codons, amino_acids))
        translated = ''
        for start in range(0, len(self.seq), 3):
            end = start + 3
            codon = self.seq[start:end]
            if len(codon) == 3:  
                aa = codon_table[codon]
                translated += aa
        return translated
    
    def get_protein_len(self):
        return len(self.seq) // 3
    

# create a calss for protein sequences
class Proteins(Sequences):

    # attributes and methods inherited from superclass Sequence

    # if we would like to add a description for the protein sequence (which is not in the superclass)
    def __init__(self, id, seq, descr = ''): # creating a default for description (of an empty string) - therefore is not necessary
        super().__init__(id, seq)
        self._descr = descr

    @ property
    def descr(self):
        return self._descr
    
    # create method to calculate percentage of hydrophobic residues
    def calc_hydrophobic_perc(self, dp = 2):
        hydrophobic = ['A', 'I', 'L', 'M', 'F', 'W', 'Y', 'V']
        hydrophobic_aas = 0
        for aa in hydrophobic:
            count_aa = self.seq.count(aa)
            hydrophobic_aas += count_aa
        perc = ((hydrophobic_aas)/len(self.seq))*100
        return round(perc, dp)
    
    def get_protein_len(self):
        return len(self.seq)
    
    def __str__(self): # we can choose to overwrite the function to include the description (or leave the superclass method and print just the id and seq)
        return f'ID: {self.id}, Sequence: {self.seq}, Description: {self.descr}'
    
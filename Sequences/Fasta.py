
class Fasta:
    def __init__(self, file):
        self._file = file

    @ property 
    def file(self):
        return self._file
    
    def parse_fasta(self, class_type):
        with open(self.file) as filehandle:
            for line in filehandle:
                id = line.rstrip().lstrip('>')
                seq = next(filehandle).rstrip()
                seq_class = class_type(id, seq)
                yield seq_class
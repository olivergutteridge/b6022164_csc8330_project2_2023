from Bio import SeqIO

class Sequence:

    def __init__(self, file) -> None:
        """Constructore for Sequence class"""
        # specify file type as fasta
        ftype = "fasta"
        self.fname = file
        # parse all records and append to instance of class
        self.records = list(SeqIO.parse(file, ftype))

    def __str__(self) -> str:
        """String representation of Sequence class"""
        # empty string
        output = ""
        for record in self.records:
            # add record to empty string
            output += f"{record}\n"
        return output
    
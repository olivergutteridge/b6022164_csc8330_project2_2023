from Bio import SeqIO

class Sequence:

    def __init__(self, file) -> None:
        """Constructore for Sequence class"""
        ftype = "fasta"
        self.fname = file
        self.records = list(SeqIO.parse(file, ftype))

    def __str__(self) -> str:
        """String representation of Sequence class"""
        output = ""
        for record in self.records:
            output += f"{record}\n"
        return output
    
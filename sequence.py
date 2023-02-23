from Bio import SeqIO

class Sequence:

    def __init__(self, file) -> None:
        ftype = "fasta"
        self.fname = file
        self.records = list(SeqIO.parse(file, ftype))

    def __str__(self) -> str:
        output = ""
        for record in self.records:
            output += f"{record}\n"
        return output
    
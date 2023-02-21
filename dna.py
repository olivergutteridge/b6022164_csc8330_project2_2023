from Bio.Seq import Seq
from Bio import SeqIO
import os

class DNA:

    def __init__(self, file) -> None:
        self.records = list(SeqIO.parse(file, "fasta"))

    def __str__(self) -> str:
        n = 0
        output = ""
        for record in self.records:
            n += 1
            output += f"{n}:\n{record}\n"
        return output

    def statistics(self):
        output = ""
        for record in self.records:
            length = len(record.seq)
            A = record.seq.count("A")
            T = record.seq.count("T")
            C = record.seq.count("C")
            G = record.seq.count("G")
            output += f"{record.id}\nLength: {length}\nA: {A}   ({round(A/length*100)}%)\nT: {T}   ({round(T/length*100)}%)\nC: {C}   ({round(C/length*100)}%)\nG: {G}   ({round(G/length*100)}%)\nGC content: {(G+C)/length}\n\n"
        return output

    def translate(self, prefix):
        path = os.path.join("./", f"{prefix}")
        os.mkdir(path)
        for record in self.records:
            forward = Seq(record.seq)
            reverse = forward.reverse_complement()
            file = open(f"./{prefix}/{record.id}.fa", "w")
            file.write(">"+record.id+"|+1\n")
            file.write(str(forward.translate()))
            file.write("\n>"+record.id+"|+2\n")
            file.write(str(forward[1:].translate()))
            file.write("\n>"+record.id+"|+3\n")
            file.write(str(forward[2:].translate()))
            file.write("\n>"+record.id+"|-1\n")
            file.write(str(reverse.translate()))
            file.write("\n>"+record.id+"|-2\n")
            file.write(str(reverse[1:].translate()))
            file.write("\n>"+record.id+"|-3\n")
            file.write(str(reverse[2:].translate()))
            file.close()
        return "Sequences have been translated."
            
            
if __name__ == "__main__":
    x = DNA("sequences.fa")
    print(x.statistics())
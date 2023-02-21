from Bio.Seq import Seq
from Bio import SeqIO
import os
import matplotlib.pyplot as plt

class DNA:

    def __init__(self, file, f) -> None:
        self.records = list(SeqIO.parse(file, f))

    def __str__(self) -> str:
        n = 0
        output = ""
        for record in self.records:
            n += 1
            output += f"{n}:\n{record}\n"
        return output

    def statistics(self, fname):
        os.mkdir(path = f"./{fname}")
        output = ""
        for record in self.records:
            length = len(record.seq)
            A = record.seq.count("A")
            T = record.seq.count("T")
            C = record.seq.count("C")
            G = record.seq.count("G")
            bases = ["A", "T", "C", "G"]
            values = [A, T, C, G]
            plt.figure(figsize = (15,7), dpi = 600)
            plt.bar(bases, values, width = 0.5, color = ["blue", "red", "yellow", "green"])
            plt.xlabel("Base")
            plt.ylabel("Count")
            plt.title(f"Base counts for {record.id}")
            plt.savefig(f"./{fname}/{record.id}.png")
            output += f"{record.id}\nLength: {length}\nA: {A}   ({round(A/length*100)}%)\nT: {T}   ({round(T/length*100)}%)\nC: {C}   ({round(C/length*100)}%)\nG: {G}   ({round(G/length*100)}%)\nGC: {round((G+C)/length*100)}%\n\n"
        file = open(f"./{fname}/summary_stats.txt", "w")
        file.write(output)
        file.close()
        return "Sattistical analysis is complete"

    def translate(self, fname):
        path = os.path.join("./", f"{fname}")
        os.mkdir(path)
        for record in self.records:
            forward = Seq(record.seq)
            reverse = forward.reverse_complement()
            file = open(f"./{fname}/{record.id}.fa", "w")
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
    x = DNA("sequences.fa", "fasta")
    y = DNA("sequence.gb", "genbank")
    print(x.statistics("test"))
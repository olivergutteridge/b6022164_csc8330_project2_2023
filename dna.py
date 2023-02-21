from Bio.Seq import Seq
from Bio import SeqIO
import os
import matplotlib.pyplot as plt
from geneblocks import CommonBlocks

class DNA:

    def __init__(self, file, ftype) -> None:
        self.records = list(SeqIO.parse(file, ftype))

    def __str__(self) -> str:
        output = ""
        for record in self.records:
            output += f"{record}\n"
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
            plt.bar(bases, values, width = 0.5, color = ["blue"])
            plt.xlabel("Base")
            plt.ylabel("Count")
            plt.title(f"Base counts for {record.id}")
            plt.savefig(f"./{fname}/{record.id}.png")
            output += f"{record.id}\nLength: {length}\nA: {A}   ({round(A/length*100)}%)\nT: {T}   ({round(T/length*100)}%)\nC: {C}   ({round(C/length*100)}%)\nG: {G}   ({round(G/length*100)}%)\nGC: {round((G+C)/length*100)}%\n\n"
        file = open(f"./{fname}/summary_stats.txt", "w")
        file.write(output)
        file.close()
        return "Sattistical analysis complete"
    
    def geneblocks(self, fname):
        sequences = {}
        for record in self.records:
            sequences[f"{record.id}"]=f"{record.seq}"
        common_blocks = CommonBlocks.from_sequences(sequences)
        ax = common_blocks.plot_common_blocks()
        ax[0].figure.savefig(f"{fname}.png", bbox_inches = "tight")
        return "DNA comparison complete"

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
        return "Translation complete"
            
            
if __name__ == "__main__":
    x = DNA("sequences.fa", "fasta")
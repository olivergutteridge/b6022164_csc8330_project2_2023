from Bio.Seq import Seq
from Bio import SeqIO
import os
import csv
import matplotlib.pyplot as plt
import numpy as np
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
        ids = []
        A = []
        T = []
        C = []
        G = []
        f = open(f"./{fname}/stats.csv", "w", encoding="UTF8")
        writer = csv.writer(f)
        headers = ["sequence_id", "length", "A", "T", "C", "G", "GC"]
        writer.writerow(headers)
        for record in self.records:
            ids.append(record.id)
            length = len(record.seq)
            a = record.seq.count("A")
            t = record.seq.count("T")
            c = record.seq.count("C")
            g = record.seq.count("G")
            gc = g + c
            data = [record.id, length, a, t, c, g, gc]
            writer.writerow(data)
            A.append(a)
            T.append(t)
            C.append(c)
            G.append(g)
            output += f"{record.id}\nLength: {length}\nA: {a}   ({round(a/length*100)}%)\nT: {t}   ({round(t/length*100)}%)\nC: {c}   ({round(c/length*100)}%)\nG: {g}   ({round(g/length*100)}%)\nGC: {round((g+c)/length*100)}%\n\n"
        A = np.array(A)
        T = np.array(T)
        C = np.array(C)
        G = np.array(G)
        plt.figure(figsize = (15,7), dpi = 600)
        plt.bar(ids, A, color = "navy")
        plt.bar(ids, T, bottom = A, color = "lime")
        plt.bar(ids, C, bottom = A + T, color = "cyan")
        plt.bar(ids, G, bottom = A + T + C, color = "fuchsia") 
        plt.xlabel("Sequence IDs")
        plt.ylabel("Base Frequency")
        plt.legend(["A", "T", "C", "G"])
        plt.savefig(f"./{fname}/base_graph.png")
        file = open(f"./{fname}/summary.txt", "w")
        file.write(output)
        file.close()
        return "Statistical analysis complete"
    
    def geneblocks(self, fname):
        sequences = {}
        for record in self.records:
            sequences[f"{record.id}"]=f"{record.seq}"
        print(sequences)
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
    x.statistics("test")

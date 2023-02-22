from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import csv
import matplotlib.pyplot as plt
import numpy as np
from geneblocks import CommonBlocks

class Sequence:

    def __init__(self, fname, ftype) -> None:
        self.fname = fname
        self.records = list(SeqIO.parse(fname, ftype))

    def __str__(self) -> str:
        output = ""
        for record in self.records:
            output += f"{record}\n"
        return output

class DNA(Sequence):

    def __init__(self, file, ftype) -> None:
        super().__init__(file, ftype)
        self.lengths = []
        self.A = []
        self.T = []
        self.C = []
        self.G = []
        self.GC = []

    def base_count(self):
        for record in self.records:
            self.lengths.append(len(record.seq))
            self.A.append(record.seq.count("A"))
            self.T.append(record.seq.count("T"))
            self.C.append(record.seq.count("C"))
            self.G.append(record.seq.count("G"))
            self.GC.append(record.seq.count("G")+record.seq.count("C"))
        return "bases counted"
    
    def base_txt(self, fname):
        output = ""
        x = len(self.records)
        for i in range(x):
            output += f"{self.records[i].id}\nLength: {self.lengths[i]}\nA: {self.A[i]}   ({round(self.A[i]/self.lengths[i]*100)}%)\nT: {self.T[i]}   ({round(self.T[i]/self.lengths[i]*100)}%)\nC: {self.C[i]}   ({round(self.C[i]/self.lengths[i]*100)}%)\nG: {self.G[i]}   ({round(self.G[i]/self.lengths[i]*100)}%)\nGC: {round((self.G[i]+self.C[i])/self.lengths[i]*100)}%\n\n"
        file = open(f"./{fname}/base_stats.txt", "w")
        file.write(output)
        file.close()
        return "text file complete"

    def base_csv(self, fname):
        f = open(f"./{fname}/base_stats.csv", "w", encoding="UTF8")
        writer = csv.writer(f)
        headers = ["sequence_id", "length", "A", "T", "C", "G", "GC"]
        writer.writerow(headers)
        x = len(self.records)
        for i in range(x):
            data = [self.records[i].id, self.lengths[i], self.A[i], self.T[i], self.C[i], self.G[i], self.GC[i]]
            writer.writerow(data)
        return "csv file complete"
    
    def base_graph(self, fname):
        ids = []
        x = len(self.records)
        for i in range(x):
            ids.append(self.records[i].id)
        A = np.array(self.A)
        T = np.array(self.T)
        C = np.array(self.C)
        G = np.array(self.G)
        plt.figure(figsize = (15,7), dpi = 600)
        plt.bar(ids, A, color = "navy")
        plt.bar(ids, T, bottom = A, color = "lime")
        plt.bar(ids, C, bottom = A + T, color = "cyan")
        plt.bar(ids, G, bottom = A + T + C, color = "fuchsia") 
        plt.xlabel("Sequence IDs")
        plt.ylabel("Base Frequency")
        plt.legend(["A", "T", "C", "G"])
        plt.xticks(rotation=45, fontsize = 5)
        plt.savefig(f"./{fname}/base_frequency_graph.png")
        return "base frequency graph complete"

    def ORF_graph(self, fname):
        for record in self.records:
            forward = Seq(record.seq)
            reverse = forward.reverse_complement()
            plus1 = forward.translate()
            plus2 = forward[1:].translate()
            plus3 = forward[2:].translate()
            minus1 = reverse.translate()
            minus2 = reverse[1:].translate()
            minus3 = reverse[2:].translate()
                    
    def statistics(self, fname):
        print(f"Starting statistical analysis of {self.fname}")
        os.mkdir(path = f"./{fname}")
        self.base_count()
        self.base_txt(fname=fname)
        self.base_csv(fname=fname)
        self.base_graph(fname=fname)
        return "Statistical analysis complete"
    
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
            plus1 = forward.translate()
            plus2 = forward[1:].translate()
            plus3 = forward[2:].translate()
            minus1 = reverse.translate()
            minus2 = reverse[1:].translate()
            minus3 = reverse[2:].translate()
            translations = [SeqRecord(seq = Seq(plus1), id = f"{record.id}", description = "+1 translation"), SeqRecord(seq = Seq(plus2), id = f"{record.id}", description = "+2 translation"), SeqRecord(seq = Seq(plus3), id = f"{record.id}", description = "+3 translation"), SeqRecord(seq = Seq(minus1), id = f"{record.id}", description = "-1 translation"), SeqRecord(seq = Seq(minus2), id = f"{record.id}", description = "-2 translation"), SeqRecord(seq = Seq(minus3), id = f"{record.id}", description = "-3 translation")]
            file = open(f"./{fname}/{record.id}.fa", "w")
            SeqIO.write(translations, file, "fasta")
            file.close()
        return "Translation complete"
                   
if __name__ == "__main__":
    x = DNA("hoxC_sequences.fa", "fasta")
    print(x.statistics("test"))

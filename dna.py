import statistics
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import csv
import matplotlib.pyplot as plt
import numpy as np

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
        self.plus1 = []
        self.plus2 = []
        self.plus3 = []
        self.minus1 = []
        self.minus2 = []
        self.minus3 = []

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
            output += f"{self.records[i].id}\nLength: {self.lengths[i]}\nA: {self.A[i]}   ({round(self.A[i]/self.lengths[i]*100)}%)\nT: {self.T[i]}   ({round(self.T[i]/self.lengths[i]*100)}%)\nC: {self.C[i]}   ({round(self.C[i]/self.lengths[i]*100)}%)\nG: {self.G[i]}   ({round(self.G[i]/self.lengths[i]*100)}%)\nGC: {round((self.GC[i])/self.lengths[i]*100)}%\n\n"
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
        plt.xticks(rotation = 45, fontsize = 5)
        plt.savefig(f"./{fname}/base_frequency.png")
        plt.close()
        return "base frequency graph complete"
        
    def locate_ORFs(self):
        self.translate()
        x = len(self.records)
        all_orfs = []
        for i in range(x):
            translations = [self.plus1[i], self.plus2[i], self.plus3[i], self.minus1[i], self.minus2[i], self.minus3[i]]
            current_orf = []
            orfs = []
            for tr in translations:
                for aa in tr:
                    if aa == "*":
                        if current_orf:
                            for p in current_orf:
                                if len(p) > 30:
                                    orfs.append(p)
                                else:
                                    pass
                            current_orf = []
                    else:
                        if aa == "M":
                            current_orf.append("")
                        for i in range(len(current_orf)):
                            current_orf[i] += aa
            all_orfs.append(orfs)
        return all_orfs

    def orfGC_graph(self, fname):
        avg_orf_length  =  []
        for orfs in self.locate_ORFs():
            orf_lengths = []
            for orf in orfs:
                orf_lengths.append(len(orf))
            avg_length = statistics.mean(orf_lengths)
            avg_orf_length.append(avg_length)
        gc = []
        for x1, x2 in zip(self.GC, self.lengths):
            gc.append(round(x1/x2*100))
        plt.figure(figsize = (9,9), dpi = 600)
        plt.scatter(gc, avg_orf_length, color = "navy")
        plt.title(f"GC content vs Average ORF length for sequences in f{self.fname}")
        plt.xlabel("GC content (%)")
        plt.ylabel("Average ORF length (aa)")
        plt.savefig(f"./{fname}/GC_ORF.png")
        plt.close()
        return "gc vs avg orf length graph complete"

    def statistics(self, fname):
        print(f"Starting statistical analysis of {self.fname}")
        os.mkdir(path = f"./{fname}")
        self.base_count()
        self.base_txt(fname=fname)
        self.base_csv(fname=fname)
        self.base_graph(fname=fname)
        self.orfGC_graph(fname=fname)
        return "Statistical analysis complete"

    def translate(self):
        for record in self.records:
            forward = Seq(record.seq)
            reverse = forward.reverse_complement()
            self.plus1.append(forward.translate())
            self.plus2.append(forward[1:].translate())
            self.plus3.append(forward[2:].translate())
            self.minus1.append(reverse.translate())
            self.minus2.append(reverse[1:].translate())
            self.minus3.append(reverse[2:].translate())
        return "Sequences translated"

    def translate_tofile(self, fname):
        self.translate()
        path = os.path.join("./", f"{fname}")
        os.mkdir(path)
        x = len(self.records)
        for i in range(x):
            translations = [SeqRecord(seq = Seq(self.plus1[i]), id = f"{self.records[i].id}", description = "+1 translation"), SeqRecord(seq = Seq(self.plus2[i]), id = f"{self.records[i].id}", description = "+2 translation"), SeqRecord(seq = Seq(self.plus3[i]), id = f"{self.records[i].id}", description = "+3 translation"), SeqRecord(seq = Seq(self.minus1[i]), id = f"{self.records[i].id}", description = "-1 translation"), SeqRecord(seq = Seq(self.minus2[i]), id = f"{self.records[i].id}", description = "-2 translation"), SeqRecord(seq = Seq(self.minus3[i]), id = f"{self.records[i].id}", description = "-3 translation")]
            file = open(f"./{fname}/{self.records[i].id}.fa", "w")
            SeqIO.write(translations, file, "fasta")
            file.close()
        return "Translation complete"
                   
if __name__ == "__main__":
    x = DNA("hoxC_sequences.fa", "fasta")
    print(x.statistics("t"))

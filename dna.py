import statistics
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import csv
import matplotlib.pyplot as plt
import numpy as np
from sequence import Sequence
import re

class DNA(Sequence):

    def __init__(self, file, odir, orf_cutoff, basic_stat, complex_stat, translation) -> None:
        super().__init__(file)
        self.file = str(file)
        self.odir = odir
        self.orf_cutoff = orf_cutoff
        self.basic_stat = basic_stat
        self.complex_stat = complex_stat
        self.translation = translation
        self.lengths = []
        self.A = []
        self.T = []
        self.C = []
        self.G = []
        self.GC = []
        self.reverse = []
        self.plus1 = []
        self.plus2 = []
        self.plus3 = []
        self.minus1 = []
        self.minus2 = []
        self.minus3 = []
        self.tr_plus1 = []
        self.tr_plus2 = []
        self.tr_plus3 = []
        self.tr_minus1 = []
        self.tr_minus2 = []
        self.tr_minus3 = []

    def check_dna(self):
        for record in self.records:
            for base in record.seq:
                if base == "A" or "T" or "C" or "G":
                    pass
                else:
                    raise TypeError("Input sequences must be DNA")

    def base_count(self):
        for record in self.records:
            self.lengths.append(len(record.seq))
            self.A.append(record.seq.count("A"))
            self.T.append(record.seq.count("T"))
            self.C.append(record.seq.count("C"))
            self.G.append(record.seq.count("G"))
            self.GC.append(record.seq.count("G")+record.seq.count("C"))
        return "bases counted"

    def reading_frames(self):
        for record in self.records:
            forward = record.seq
            reverse = forward.reverse_complement()
            self.reverse.append(reverse)
            self.plus1.append(forward)
            self.plus2.append(forward[1:])
            self.plus3.append(forward[2:])
            self.minus1.append(reverse)
            self.minus2.append(reverse[1:])
            self.minus3.append(reverse[2:])
        return "sequences shifted into all 6 reading frames"

    def find_orfs(self):
        x = len(self.records)
        all_orfs = []
        pattern = re.compile(r'ATG(?:.{3})*?(?:TAA|TAG|TGA)')
        for i in range(x):
            orfs = []
            reading_frames = [self.plus1[i], self.plus2[i], self.plus3[i], self.minus1[i], self.minus2[i], self.minus3[i]]
            matches = [pattern.finditer(str(frame)) for frame in reading_frames]
            for seq_matches in matches:
                for match in seq_matches:
                    for match in seq_matches:
                        orf = match.group()
                        if len(orf) >= self.orf_cutoff:
                            orfs.append(orf)
            all_orfs.append(orfs)
        return all_orfs
    
    def base_txt(self, odir):
        output = ""
        x = len(self.records)
        for i in range(x):
            output += f"{self.records[i].id}\nLength: {self.lengths[i]}\nA: {self.A[i]}   ({round(self.A[i]/self.lengths[i]*100)}%)\nT: {self.T[i]}   ({round(self.T[i]/self.lengths[i]*100)}%)\nC: {self.C[i]}   ({round(self.C[i]/self.lengths[i]*100)}%)\nG: {self.G[i]}   ({round(self.G[i]/self.lengths[i]*100)}%)\nGC: {round((self.GC[i])/self.lengths[i]*100)}%\n\n"
        file = open(f"./{odir}/stats/base_stats.txt", "w")
        file.write(output)
        file.close()
        return "text file complete"

    def base_csv(self, odir):
        f = open(f"./{odir}/stats/base_stats.csv", "w", encoding="UTF8")
        writer = csv.writer(f)
        headers = ["sequence_id", "length", "A", "T", "C", "G", "GC"]
        writer.writerow(headers)
        x = len(self.records)
        for i in range(x):
            data = [self.records[i].id, self.lengths[i], self.A[i], self.T[i], self.C[i], self.G[i], self.GC[i]]
            writer.writerow(data)
        return "csv file complete"
    
    def base_graph(self, odir):
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
        plt.savefig(f"./{odir}/stats/base_frequency.png")
        plt.close()
        return "base frequency graph complete"

    def orfGC_graph(self, odir):
        avg_orf_length  =  []
        for orfs in self.find_orfs():
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
        plt.title(f"GC content vs Average ORF length for sequences in {self.file}")
        plt.xlabel("GC content (%)")
        plt.ylabel("Average ORF length (base)")
        plt.savefig(f"./{odir}/stats/GC_ORF.png")
        plt.close()
        return "gc vs avg orf length graph complete"

    def r_complement(self, odir):
        x = len(self.records)
        reverses = []
        for i in range(x):
            reverse = SeqRecord(seq = Seq(self.reverse[i]), id = f"{self.records[i].id}", description = "reverse complement")
            reverses.append(reverse)
        file = open(f"./{odir}/translations/rComp_{self.file}.fa", "w")
        SeqIO.write(reverses, file, "fasta")
        file.close()
        return f"reverse complements wrote to {self.odir}_rComp.fa"

    def translate(self):
        for record in self.records:
            forward = Seq(record.seq)
            reverse = forward.reverse_complement()
            self.tr_plus1.append(forward.translate())
            self.tr_plus2.append(forward[1:].translate())
            self.tr_plus3.append(forward[2:].translate())
            self.tr_minus1.append(reverse.translate())
            self.tr_minus2.append(reverse[1:].translate())
            self.tr_minus3.append(reverse[2:].translate())
        return "sequences translated"

    def translate_tofile(self, odir):
        self.translate()
        x = len(self.records)
        for i in range(x):
            translations = [SeqRecord(seq = Seq(self.tr_plus1[i]), id = f"{self.records[i].id}", description = "+1 translation"), SeqRecord(seq = Seq(self.tr_plus2[i]), id = f"{self.records[i].id}", description = "+2 translation"), SeqRecord(seq = Seq(self.tr_plus3[i]), id = f"{self.records[i].id}", description = "+3 translation"), SeqRecord(seq = Seq(self.tr_minus1[i]), id = f"{self.records[i].id}", description = "-1 translation"), SeqRecord(seq = Seq(self.tr_minus2[i]), id = f"{self.records[i].id}", description = "-2 translation"), SeqRecord(seq = Seq(self.tr_minus3[i]), id = f"{self.records[i].id}", description = "-3 translation")]
            file = open(f"./{odir}/translations/{self.records[i].id}_tr.fa", "w")
            SeqIO.write(translations, file, "fasta")
            file.close()
        return "translation complete"

    def b_stat(self):
        print(f"Starting basic statistical analysis of {self.file}")
        os.mkdir(path = f"./{self.odir}/stats")
        self.base_count()
        self.reading_frames()
        self.base_txt(odir=self.odir)
        self.base_csv(odir=self.odir)
        self.base_graph(odir=self.odir)

    def c_stat(self):
        print(f"Starting complex statistical analysis of {self.file}")
        os.mkdir(path = f"./{self.odir}/stats")
        self.base_count()
        self.reading_frames()
        self.base_txt(odir=self.odir)
        self.base_csv(odir=self.odir)
        self.base_graph(odir=self.odir)
        self.orfGC_graph(odir=self.odir)

    def tr(self):
        print(f"Translating sequences in {self.file}")
        os.mkdir(path = f"./{self.odir}/translations")
        self.reading_frames()
        self.r_complement(odir=self.odir)
        self.translate_tofile(odir=self.odir)

    def core(self):
        self.check_dna()
        os.mkdir(path = f"./{self.odir}")
        if self.basic_stat and self.complex_stat:
            print("The arguments '--basic-stat' and '--complex-stat' cannot be selected together")
        elif self.complex_stat and not self.basic_stat and not self.translate:
            self.c_stat()
        elif self.translation and not self.basic_stat and not self.complex_stat:
            self.tr()
        elif self.basic_stat and self.translate and not self.complex_stat:
            self.b_stat()
            self.tr()
        elif self.complex_stat and self.translate and not self.basic_stat:
            self.complex_stat()
            self.tr()
        return "Analysis complete"

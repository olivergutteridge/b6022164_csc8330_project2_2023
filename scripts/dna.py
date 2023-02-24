import statistics
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import shutil
import csv
import matplotlib.pyplot as plt
import numpy as np
from scripts.sequence import Sequence
import re

class DNA(Sequence):

    def __init__(self, file, out_dir, min_orf, basic_stats, complex_stats, translate, save_orfs) -> None:
        """Constructor for DNA class"""
        super().__init__(file)
        self.file = str(file)
        self.out_dir = out_dir
        self.min_orf = min_orf
        self.basic_stats = basic_stats
        self.complex_stats = complex_stats
        self.translate = translate
        self.save_orfs = save_orfs
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
        self.all_orfs = []

    def check_dna(self):
        """Checks if input is DNA"""
        for record in self.records:
            for base in record.seq:
                if base == "A" or "T" or "C" or "G":
                    pass
                else:
                    raise TypeError("Input sequences must be DNA")

    def make_outdir(self):
        """Make out_dir and overwrite if already exists"""
        path = f"./{self.out_dir}"
        if os.path.exists(path=path):
            shutil.rmtree(path=path)
        os.mkdir(path=path)

    def base_count(self):
        """Calculates base counts, lengths and GC content for each sequence in input"""
        for record in self.records:
            self.lengths.append(len(record.seq))
            self.A.append(record.seq.count("A"))
            self.T.append(record.seq.count("T"))
            self.C.append(record.seq.count("C"))
            self.G.append(record.seq.count("G"))
            self.GC.append(record.seq.count("G")+record.seq.count("C"))
        return "bases counted"

    def reading_frames(self):
        """Transforms each sequence in input into all 6 reading frames"""
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
        """Locates ORFs within each sequence of input"""
        pattern = re.compile(r'ATG(?:.{3})*?(?:TAA|TAG|TGA)')
        for i in range(len(self.records)):
            orfs = []
            reading_frames = [self.plus1[i], self.plus2[i], self.plus3[i], self.minus1[i], self.minus2[i], self.minus3[i]]
            matches = [pattern.finditer(str(frame)) for frame in reading_frames]
            for seq_matches in matches:
                for match in seq_matches:
                    for match in seq_matches:
                        orf = match.group()
                        if len(orf) >= self.min_orf:
                            orfs.append(orf)
            self.all_orfs.append(orfs)
        if self.save_orfs:
            os.mkdir(path = f"./{self.out_dir}/orfs")
            for i in range(len(self.records)):
                file = open(f"./{self.out_dir}/orfs/{self.records[i].id}_orfs", "w")
                current_orfs = self.all_orfs[i]
                n = 0
                for orf in current_orfs:
                    orf_tofile = SeqRecord(seq = Seq(orf), id = f"{self.records[i].id}", description = f"orf_{n}")
                    SeqIO.write(orf_tofile, file, "fasta")
                    n += 1
                file.close()
        print(f"ORFs saved to ./{self.out_dir}/orfs")
        return "orfs located"
    
    def base_txt(self, out_dir):
        """Creates a text file containing basic statistics about each sequence in input"""
        output = ""
        for i in range(len(self.records)):
            output += f"{self.records[i].id}\nLength: {self.lengths[i]}\nA: {self.A[i]}   ({round(self.A[i]/self.lengths[i]*100)}%)\nT: {self.T[i]}   ({round(self.T[i]/self.lengths[i]*100)}%)\nC: {self.C[i]}   ({round(self.C[i]/self.lengths[i]*100)}%)\nG: {self.G[i]}   ({round(self.G[i]/self.lengths[i]*100)}%)\nGC: {round((self.GC[i])/self.lengths[i]*100)}%\nORFs: {len(self.all_orfs[i])}\n\n"
        file = open(f"./{out_dir}/stats/base_stats.txt", "w")
        file.write(output)
        file.close()
        return "text file complete"

    def base_csv(self, out_dir):
        """Creates a csv file containing basic statistics about each sequence in input"""
        file = open(f"./{out_dir}/stats/base_stats.csv", "w", encoding="UTF8")
        writer = csv.writer(file)
        headers = ["sequence_id", "length", "A", "T", "C", "G", "GC", "ORFs"]
        writer.writerow(headers)
        for i in range(len(self.records)):
            data = [self.records[i].id, self.lengths[i], self.A[i], self.T[i], self.C[i], self.G[i], self.GC[i], len(self.all_orfs[i])]
            writer.writerow(data)
        file.close()
        return "csv file complete"
    
    def base_freq_graph(self, out_dir):
        """Creates a stacked bar plot of base frequencies for sequence in input"""
        ids = []
        for i in range(len(self.records)):
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
        plt.savefig(f"./{out_dir}/stats/base_frequency.png")
        plt.close()
        return "base frequency graph complete"

    def seq_orfGC_graph(self, out_dir):
        """Creates a scatter plot of GC content vs avg ORF length for sequence in input"""
        avg_orf_length  =  []
        for orfs in self.all_orfs:
            orf_lengths = []
            for orf in orfs:
                orf_lengths.append(len(orf))
            avg_length = statistics.mean(orf_lengths)
            avg_orf_length.append(avg_length)
        gc = []
        for x1, x2 in zip(self.GC, self.lengths):
            gc.append(round(x1/x2*100))
        avg_orf_length = np.array(avg_orf_length)
        gc = np.array(gc)
        plt.figure(figsize = (9,9), dpi = 600)
        plt.scatter(gc, avg_orf_length, color = "navy")
        a, b = np.polyfit(gc, avg_orf_length, 1)
        plt.plot(gc, a*gc+b, color = "lime")
        plt.title(f"Sequence GC content vs Average ORF length for sequences in {self.file}")
        plt.xlabel("Sequence GC content (%)")
        plt.ylabel("Average ORF length (bp)")
        plt.savefig(f"./{out_dir}/stats/GC_ORF.png")
        plt.close()
        return "gc vs avg orf length graph complete"

    def seq_length_orf_graph(self, out_dir):
        """Creates a scatter plot of seq length vs total ORFs for sequence in input"""
        seq_lengths = np.array([length for length in self.lengths])
        orf_count = []
        for orfs in self.all_orfs:
            orf_count.append(len(orfs))
        orf_count = np.array(orf_count)
        plt.figure(figsize = (9,9), dpi = 600)
        plt.scatter(seq_lengths, orf_count, color = "lime")
        a, b = np.polyfit(seq_lengths, orf_count, 1)
        plt.plot(seq_lengths, a*seq_lengths+b, color = "navy")
        plt.title(f"Sequence length vs Number of ORFs for sequences in {self.file}")
        plt.xlabel("Sequence length (bp)")
        plt.ylabel("Number of ORFs")
        plt.savefig(f"./{out_dir}/stats/seqLength_ORF.png")
        plt.close()
        return "seq length vs orf count graph complete"

    def r_complement(self, out_dir):
        """Calculates reverse complement for sequence in input"""
        reverses = []
        for i in range(len(self.records)):
            reverse = SeqRecord(seq = Seq(self.reverse[i]), id = f"{self.records[i].id}", description = "reverse complement")
            reverses.append(reverse)
        file = open(f"./{out_dir}/rComps_{self.file}.fa", "w")
        SeqIO.write(reverses, file, "fasta")
        file.close()
        return f"reverse complements wrote to {self.out_dir}_rComp.fa"

    def translate_seqs(self, out_dir):
        """Translates all sequences in all reading frames for sequence in input"""
        tr_plus1 = []
        tr_plus2 = []
        tr_plus3 = []
        tr_minus1 = []
        tr_minus2 = []
        tr_minus3 = []
        for record in self.records:
            forward = Seq(record.seq)
            reverse = forward.reverse_complement()
            tr_plus1.append(forward.translate())
            tr_plus2.append(forward[1:].translate())
            tr_plus3.append(forward[2:].translate())
            tr_minus1.append(reverse.translate())
            tr_minus2.append(reverse[1:].translate())
            tr_minus3.append(reverse[2:].translate())
        for i in range(len(self.records)):
            translations = [SeqRecord(seq = Seq(tr_plus1[i]), id = f"{self.records[i].id}", description = "+1 translation"), SeqRecord(seq = Seq(tr_plus2[i]), id = f"{self.records[i].id}", description = "+2 translation"), SeqRecord(seq = Seq(tr_plus3[i]), id = f"{self.records[i].id}", description = "+3 translation"), SeqRecord(seq = Seq(tr_minus1[i]), id = f"{self.records[i].id}", description = "-1 translation"), SeqRecord(seq = Seq(tr_minus2[i]), id = f"{self.records[i].id}", description = "-2 translation"), SeqRecord(seq = Seq(tr_minus3[i]), id = f"{self.records[i].id}", description = "-3 translation")]
            file = open(f"./{out_dir}/translations/{self.records[i].id}_tr.fa", "w")
            SeqIO.write(translations, file, "fasta")
            file.close()
        return "translation complete"

    def b_stats(self):
        """Calls each function neccessary for --basic-stats"""
        print(f"Starting basic statistical analysis of {self.file}")
        os.mkdir(path = f"./{self.out_dir}/stats")
        self.base_count()
        self.reading_frames()
        self.find_orfs()
        self.base_txt(out_dir=self.out_dir)
        self.base_csv(out_dir=self.out_dir)
        self.base_freq_graph(out_dir=self.out_dir)
        print(f"Basic statistical analysis saved to ./{self.out_dir}/stats")
        return None

    def c_stats(self):
        """Calls each function neccessary for --complex-stats"""
        print(f"Starting complex statistical analysis of {self.file}")
        os.mkdir(path = f"./{self.out_dir}/stats")
        self.base_count()
        self.reading_frames()
        self.find_orfs()
        self.base_txt(out_dir=self.out_dir)
        self.base_csv(out_dir=self.out_dir)
        self.base_freq_graph(out_dir=self.out_dir)
        self.seq_orfGC_graph(out_dir=self.out_dir)
        self.seq_length_orf_graph(out_dir=self.out_dir)
        print(f"Complex statistical analysis saved to ./{self.out_dir}/stats")
        return None

    def orfs(self):
        """Calls each function neccessary for --save-orfs"""
        print(f"Locating ORFs longer than {self.min_orf} bp in sequences in {self.file}")
        self.reading_frames()
        self.find_orfs()
        return None

    def tr(self):
        """Calls each function neccessary for --translate"""
        print(f"Translating sequences in {self.file}")
        os.mkdir(path = f"./{self.out_dir}/translations")
        self.reading_frames()
        self.r_complement(out_dir=self.out_dir)
        self.translate_seqs(out_dir=self.out_dir)
        print(f"Reverse complements saved to ./{self.out_dir}/rComps_{self.file}.fa")
        print(f"Translations saved to ./{self.out_dir}/translations")
        return None

    def core(self):
        """Logic for arguments in cli.py"""
        self.check_dna()
        self.make_outdir()
        if self.basic_stats and self.complex_stats:
            print("The arguments '--basic-stats' and '--complex-stats' cannot be selected together")
            os.rmdir(path = f"./{self.out_dir}")
        elif not self.basic_stats and not self.complex_stats and not self.translate and not self.save_orfs:
            print("Please select either '--basic-stats'/'--complex-stats' and or '--translate' and or '--save-orfs'")
            os.rmdir(path = f"./{self.out_dir}")
        elif self.basic_stats and not self.complex_stats and not self.translate:
            self.b_stats()
        elif self.complex_stats and not self.basic_stats and not self.translate:
            self.c_stats()
        elif self.save_orfs and not self.basic_stats and not self.complex_stats and not self.translate:
            self.orfs()
        elif self.translate and not self.basic_stats and not self.complex_stats:
            self.tr()
        elif self.basic_stats and self.translate and not self.complex_stats:
            self.b_stats()
            self.tr()
        elif self.complex_stats and self.translate and not self.basic_stats:
            self.c_stats()
            self.tr()
        print("Thanks for using dnaStat")
        return None
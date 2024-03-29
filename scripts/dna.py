import statistics
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW, NCBIXML
import os
import shutil
import csv
import matplotlib.pyplot as plt
import numpy as np
from sequence import Sequence
import re

class DNA(Sequence):

    def __init__(self, file, out_dir, min_orf, basic_stats, complex_stats, translate, save_orfs) -> None:
        """Constructor for DNA class"""
        # inherit constructor from Sequence class
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
        # for loop to access individual records
        for record in self.records:
            # for loop to access bases
            for base in record.seq:
                # if base DNA pass
                if base == "A" or "T" or "C" or "G":
                    pass
                # else raise TypeError
                else:
                    raise TypeError("Input sequences must be DNA")

    def make_outdir(self):
        """Make out_dir and overwrite if already exists"""
        # set out_dir path
        path = f"./{self.out_dir}"
        # if already exists, remove
        if os.path.exists(path=path):
            shutil.rmtree(path=path)
        # create out_dir
        os.mkdir(path=path)
        return f"./{self.out_dir} created"

    def base_count(self):
        """Calculates base counts, lengths and GC content for each sequence in input"""
        # for loop to access individual records
        for record in self.records:
            # append length, base counts and GC content to instance of class
            self.lengths.append(len(record.seq))
            self.A.append(record.seq.count("A"))
            self.T.append(record.seq.count("T"))
            self.C.append(record.seq.count("C"))
            self.G.append(record.seq.count("G"))
            self.GC.append(record.seq.count("G")+record.seq.count("C"))
        return "bases frequencies calculated"

    def reading_frames(self):
        """Transforms each sequence in input into all 6 reading frames"""
        # for loop to access individual records
        for record in self.records:
            # get forwrad and reverse strands
            forward = record.seq
            reverse = forward.reverse_complement()
            # append reading frames to instance of class
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
        # pattern for an ORF using regex
        pattern = re.compile(r'ATG(?:.{3})*?(?:TAA|TAG|TGA)')
        # for loop in length of number of sequences in input
        for i in range(len(self.records)):
            orfs = []
            reading_frames = [self.plus1[i], self.plus2[i], self.plus3[i], self.minus1[i], self.minus2[i], self.minus3[i]]
            # search for matches in all reading frames
            matches = [pattern.finditer(str(frame)) for frame in reading_frames]
            for seq_matches in matches:
                for match in seq_matches:
                    # access individual match (ORF)
                    for match in seq_matches:
                        # return substring of match
                        orf = match.group()
                        # only append ORFs >= min length
                        if len(orf) >= self.min_orf:
                            orfs.append(orf)
            self.all_orfs.append(orfs)
        # if --save-orfs selected
        if self.save_orfs:
            # make directory 
            os.mkdir(path = f"./{self.out_dir}/orfs")
            # for loop in length of number of sequences in input
            for i in range(len(self.records)):
                file = open(f"./{self.out_dir}/orfs/{self.records[i].id}_orfs.fa", "w")
                current_orfs = self.all_orfs[i]
                # counter to label orfs
                n = 0
                # access single ORFs
                for orf in current_orfs:
                    # create SeqRecord object or ORF
                    orf_tofile = SeqRecord(seq = Seq(orf), id = f"{self.records[i].id}", description = f"orf_{n}")
                    # write to file
                    SeqIO.write(orf_tofile, file, "fasta")
                    n += 1
                file.close()
            print(f"ORFs saved to ./{self.out_dir}/orfs")
            return "_orf.fa files created"
        else:   
            return "orfs located"
    
    def base_txt(self, out_dir):
        """Creates a text file containing basic statistics about each sequence in input"""
        # empty output string
        output = ""
        if self.complex_stats:
            avg_orf = self.avg_orf()
            orf_count = self.orf_count()
        # for loop in length of number of sequences in input
        for i in range(len(self.records)):
            # add data to output string
            output += f"{self.records[i].id}\nLength: {self.lengths[i]}\nA: {self.A[i]}   ({round(self.A[i]/self.lengths[i]*100)}%)\nT: {self.T[i]}   ({round(self.T[i]/self.lengths[i]*100)}%)\nC: {self.C[i]}   ({round(self.C[i]/self.lengths[i]*100)}%)\nG: {self.G[i]}   ({round(self.G[i]/self.lengths[i]*100)}%)\nGC: {round((self.GC[i])/self.lengths[i]*100)}%\n"
            if self.complex_stats:
                output += f"ORF(s): {orf_count[i]}\nMean ORF length(bp): {round(avg_orf[i])}\n\n"
            else:
                output += "\n"
        # write output to .txt file
        file = open(f"./{out_dir}/stats/base_stats.txt", "w")
        file.write(output)
        # close file
        file.close()
        return "base_stats.txt created"

    def base_csv(self, out_dir):
        """Creates a csv file containing basic statistics about each sequence in input"""
        # open .csv file
        file = open(f"./{out_dir}/stats/base_stats.csv", "w", encoding="UTF8")
        writer = csv.writer(file)
        # write headers to file
        headers = ["sequence_id", "length", "A", "T", "C", "G", "GC"]
        if self.complex_stats:
            avg_orf = self.avg_orf()
            orf_count = self.orf_count()
            headers.extend(("ORF(s)", "Mean ORF(bp)"))
        writer.writerow(headers)
        # for loop in length of number of sequences in input
        for i in range(len(self.records)):
            # create list of data + write to file
            data = [self.records[i].id, self.lengths[i], self.A[i], self.T[i], self.C[i], self.G[i], self.GC[i]]
            if self.complex_stats:
                data.extend((orf_count[i], round(avg_orf[i])))
            writer.writerow(data)
        # close file
        file.close()
        return "base_stats.csv created"
    
    def base_freq_graph(self, out_dir):
        """Creates a stacked bar plot of base frequencies for sequence in input"""
        # empty ids list
        ids = []
        # for loop in length of number of sequences in input
        for i in range(len(self.records)):
            # append each sequence if to ids list
            ids.append(self.records[i].id)
        # create np array of each base count
        A = np.array(self.A)
        T = np.array(self.T)
        C = np.array(self.C)
        G = np.array(self.G)
        # specify figure size and dpi
        plt.figure(dpi = 600)
        # plot stacked bars for each base
        plt.bar(ids, A, color = "navy")
        plt.bar(ids, T, bottom = A, color = "lime")
        plt.bar(ids, C, bottom = A + T, color = "cyan")
        plt.bar(ids, G, bottom = A + T + C, color = "fuchsia") 
        # add labels and legend
        plt.title(f"Base counts for {self.file}")
        plt.xlabel("Sequence IDs")
        plt.ylabel("Base Frequency")
        plt.legend(["A", "T", "C", "G"])
        plt.xticks(rotation = 90)
        plt.tight_layout()
        # save figure to file and close
        plt.savefig(f"./{out_dir}/stats/base_frequency.png")
        plt.close()
        return "base_frequency.png created"
    
    def avg_orf(self):
        """Calculates avg ORF length for sequence in input"""
        # empty avg_orf_length list
        avg_orf_length  =  []
        # access orfs of each sequence 
        for orfs in self.all_orfs:
            # empty orf lengths list
            orf_lengths = []
            # access individual orfs and append to orf_lengths
            for orf in orfs:
                orf_lengths.append(len(orf))
            # get mean length
            avg_length = statistics.mean(orf_lengths)
            # append to avg_orf_length
            avg_orf_length.append(avg_length)
        return avg_orf_length
    
    def gc_content(self):
        """Calulcates % gc content for sequence in input"""
        # empty gc list
        gc = []
        # for loop over self.GC and self.lengths
        for x1, x2 in zip(self.GC, self.lengths):
            # append calculated gc content
            gc.append(round(x1/x2*100))
        return gc
    
    def orf_count(self):
        """Calculates ORF count for sequence in input"""
        # empty orf_count list
        orf_count = []
        # for loop to access individual orf counts and append orf count orf_count
        for orfs in self.all_orfs:
            orf_count.append(len(orfs))
        return orf_count

    def seq_orfgc_graph(self, out_dir):
        """Creates a scatter plot of GC content vs avg ORF length for sequence in input"""
        # np arrays of avg orf and gc content
        avg_orf_length  =  np.array(self.avg_orf())
        gc = np.array(self.gc_content())
        # specify figure size and dpi
        plt.figure(dpi = 600)
        # plot scatter plot
        plt.scatter(gc, avg_orf_length, marker = "x", color = "navy")
        # caluclate  and add regression line https://www.python-graph-gallery.com/scatterplot-with-regression-fit-in-matplotlib
        a, b = np.polyfit(gc, avg_orf_length, 1)
        plt.plot(gc, a*gc+b, color = "lime")
        # add labels
        plt.title(f"GC content against avg ORF for {self.file}")
        plt.xlabel("Sequence GC content (%)")
        plt.ylabel("Average ORF length (bp)")
        plt.tight_layout()
        # save figure to file and close
        plt.savefig(f"./{out_dir}/stats/GC_ORF.png")
        plt.close()
        return "GC_ORF.png graph complete"

    def seq_length_orf_graph(self, out_dir):
        """Creates a scatter plot of seq length vs total ORFs for sequence in input"""
        # create np array of sequence lengths and orf counts
        seq_lengths = np.array([length for length in self.lengths])
        orf_count = np.array(self.orf_count())
        # specify figure size and dpi
        plt.figure(dpi = 600)
        # plot scatter plot
        plt.scatter(seq_lengths, orf_count, marker = "x", color = "lime")
        # caluclate  and add regression line https://www.python-graph-gallery.com/scatterplot-with-regression-fit-in-matplotlib
        a, b = np.polyfit(seq_lengths, orf_count, 1)
        plt.plot(seq_lengths, a*seq_lengths+b, color = "navy")
        # add labels
        plt.title(f"Seq length against ORF count for {self.file}")
        plt.xlabel("Sequence length (bp)")
        plt.ylabel("Number of ORFs")
        plt.tight_layout()
        # save figure to file and close
        plt.savefig(f"./{out_dir}/stats/seqLength_ORF.png")
        plt.close()
        return f"seqLength_ORF.png created"

    def r_complement(self, out_dir):
        """Calculates reverse complement for sequence in input"""
        # empty revserese list
        reverses = []
        # for loop in length of number of sequences in input
        for i in range(len(self.records)):
            # create SeqRecord object for each reverse complement anmd append to reverses list
            reverse = SeqRecord(seq = Seq(self.reverse[i]), id = f"{self.records[i].id}", description = "reverse complement")
            reverses.append(reverse)
        # open file, write reverses to file in .fa format, close file
        file = open(f"./{out_dir}/rComps_{self.file}", "w")
        SeqIO.write(reverses, file, "fasta")
        file.close()
        return f"reverse complements calculated"

    def translate_seqs(self, out_dir):
        """Translates all sequences in all reading frames for sequence in input"""
        # empty list for translations of all 6 reading frames
        tr_plus1 = []
        tr_plus2 = []
        tr_plus3 = []
        tr_minus1 = []
        tr_minus2 = []
        tr_minus3 = []
        # for loop to access indivdual records
        for record in self.records:
            # create Seq instance of forward and use it to create reverse complement
            forward = Seq(record.seq)
            reverse = forward.reverse_complement()
            # translate in all 6 reading frames and append to lists
            tr_plus1.append(forward.translate())
            tr_plus2.append(forward[1:].translate())
            tr_plus3.append(forward[2:].translate())
            tr_minus1.append(reverse.translate())
            tr_minus2.append(reverse[1:].translate())
            tr_minus3.append(reverse[2:].translate())
        # for loop in length of number of sequences in input
        for i in range(len(self.records)):
            # create SeqRecords for all 6 reading frames
            translations = [SeqRecord(seq = Seq(tr_plus1[i]), id = f"{self.records[i].id}", description = "+1 translation"), SeqRecord(seq = Seq(tr_plus2[i]), id = f"{self.records[i].id}", description = "+2 translation"), SeqRecord(seq = Seq(tr_plus3[i]), id = f"{self.records[i].id}", description = "+3 translation"), SeqRecord(seq = Seq(tr_minus1[i]), id = f"{self.records[i].id}", description = "-1 translation"), SeqRecord(seq = Seq(tr_minus2[i]), id = f"{self.records[i].id}", description = "-2 translation"), SeqRecord(seq = Seq(tr_minus3[i]), id = f"{self.records[i].id}", description = "-3 translation")]
            # open file, write translations to file in .fa format, close file
            file = open(f"./{out_dir}/translations/{self.records[i].id}_tr.fa", "w")
            SeqIO.write(translations, file, "fasta")
            file.close()
        return f"translations for sequences in {self.file} saved to ./{out_dir}/translations/"
    
    def id_seqs(self):
        ids = []
        for record in self.records:
            result_handle = NCBIWWW.qblast("blastn","nt", record.seq)
            blast_result = result_handle.read()
            ids.append(blast_result[0].title)
        return ids

    def b_stats(self):
        """Calls each function neccessary for --basic-stats"""
        print(f"Starting basic statistical analysis of {self.file}")
        # create stats directory in out_dir
        os.mkdir(path = f"./{self.out_dir}/stats")
        # call all functions neccessary for --basic-stats
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
        # create stats directory in out_dir
        os.mkdir(path = f"./{self.out_dir}/stats")
        # call all functions neccessary for --complex-stats
        self.base_count()
        self.reading_frames()
        self.find_orfs()
        self.base_txt(out_dir=self.out_dir)
        self.base_csv(out_dir=self.out_dir)
        self.base_freq_graph(out_dir=self.out_dir)
        self.seq_orfgc_graph(out_dir=self.out_dir)
        self.seq_length_orf_graph(out_dir=self.out_dir)
        print(f"Complex statistical analysis saved to ./{self.out_dir}/stats")
        return None

    def orfs(self):
        """Calls each function neccessary for --save-orfs"""
        print(f"Locating ORFs longer than {self.min_orf} bp in sequences in {self.file}")
        # call all functions neccessary for --save-orfs
        self.reading_frames()
        self.find_orfs()
        return None

    def tr(self):
        """Calls each function neccessary for --translate"""
        print(f"Translating sequences in {self.file}")
        # create translations directory in out_dir
        os.mkdir(path = f"./{self.out_dir}/translations")
        # call all functions neccessary for --translate 
        self.reading_frames()
        self.r_complement(out_dir=self.out_dir)
        self.translate_seqs(out_dir=self.out_dir)
        print(f"Reverse complements saved to ./{self.out_dir}/rComps_{self.file}.fa")
        print(f"Translations saved to ./{self.out_dir}/translations")
        return None

    def core(self):
        """Logic for arguments in cli.py"""
        # check that input is DNA and make out_dir
        self.check_dna()
        self.make_outdir()
        # boolean logic for all arguments, print error message if arguments not compatible
        if self.basic_stats and self.complex_stats:
            print("The arguments '--basic-stats' and '--complex-stats' cannot be selected together")
        elif not self.basic_stats and not self.complex_stats and not self.translate and not self.save_orfs:
            print("Please select one of the following '--basic-stats', '--complex-stats', '--translate', '--save-orfs'")
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
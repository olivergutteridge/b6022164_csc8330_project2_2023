from Bio.Seq import Seq
from Bio import SeqIO

class DNA:

    def __init__(self, file) -> None:
        self.sequences = []
        for sequence in SeqIO.parse(file, "fasta"):
            self.sequences.append(sequence)

    def __str__(self) -> str:
        n = 0
        output = ""
        for sequence in self.sequences:
            n += 1
            output += f"{n}: {sequence}\n"
        return output

    def statistics(self):
        n = 0
        output = ""
        for sequence in self.sequences:
            A = 0
            T = 0
            C = 0
            G = 0
            n += 1
            for base in sequence:
                if base == "A":
                    A += 1
                elif base == "T":
                    T += 1
                elif base == "C":
                    C += 1
                elif base == "G":
                    G += 1
            output += f"{sequence.id}\nLength: {len(sequence.seq)}\nA: {A}   ({round(A/len(sequence.seq)*100)}%)\nT: {T}   ({round(T/len(sequence.seq)*100)}%)\nC: {C}   ({round(C/len(sequence.seq)*100)}%)\nG: {G}   ({round(G/len(sequence.seq)*100)}%)\nGC content: {(G+C)/len(sequence.seq)}\n\n"
        return output
            
            

    

if __name__ == "__main__":
    x = DNA("sequences.fa")
    print(x.statistics())
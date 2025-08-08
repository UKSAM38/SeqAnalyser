import re
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import GC
from Bio.Seq import Seq

class PrimerDesigner:
    def __init__(self):
        self.min_length = 18
        self.max_length = 25
        self.min_gc = 40
        self.max_gc = 60
        self.min_tm = 55
        self.max_tm = 65
        
    def calculate_tm(self, sequence):
        """Calculate melting temperature"""
        return mt.Tm_NN(sequence)
        
    def calculate_gc_content(self, sequence):
        """Calculate GC content"""
        return GC(sequence)
        
    def check_hairpin(self, sequence):
        """Check for potential hairpin structures"""
        # Simple hairpin check - look for reverse complements
        seq = Seq(sequence)
        rev_comp = str(seq.reverse_complement())
        
        # Check for complementarity (simplified)
        for i in range(len(sequence) - 3):
            for j in range(i + 4, len(sequence)):
                if sequence[i:i+4] in rev_comp:
                    return True
        return False
        
    def check_dimer(self, primer1, primer2):
        """Check for primer dimer formation"""
        # Simple dimer check
        for i in range(len(primer1) - 2):
            for j in range(len(primer2) - 2):
                if primer1[i:i+3] == str(Seq(primer2[j:j+3]).reverse_complement()):
                    return True
        return False
        
    def design_primers(self, sequence, target_start=None, target_end=None):
        """Design PCR primers for a given sequence"""
        primers = []
        
        if target_start is None:
            target_start = len(sequence) // 4
        if target_end is None:
            target_end = len(sequence) * 3 // 4
            
        # Design forward primers
        for start in range(max(0, target_start - 200), target_start + 50):
            for length in range(self.min_length, self.max_length + 1):
                if start + length > len(sequence):
                    break
                    
                primer_seq = sequence[start:start + length]
                gc_content = self.calculate_gc_content(primer_seq)
                tm = self.calculate_tm(primer_seq)
                
                if (self.min_gc <= gc_content <= self.max_gc and 
                    self.min_tm <= tm <= self.max_tm and
                    not self.check_hairpin(primer_seq)):
                    
                    primers.append({
                        'type': 'forward',
                        'sequence': primer_seq,
                        'start': start,
                        'end': start + length,
                        'length': length,
                        'gc_content': gc_content,
                        'tm': tm
                    })
                    
        # Design reverse primers
        for end in range(target_end - 50, min(len(sequence), target_end + 200)):
            for length in range(self.min_length, self.max_length + 1):
                if end - length < 0:
                    break
                    
                primer_seq = str(Seq(sequence[end - length:end]).reverse_complement())
                gc_content = self.calculate_gc_content(primer_seq)
                tm = self.calculate_tm(primer_seq)
                
                if (self.min_gc <= gc_content <= self.max_gc and 
                    self.min_tm <= tm <= self.max_tm and
                    not self.check_hairpin(primer_seq)):
                    
                    primers.append({
                        'type': 'reverse',
                        'sequence': primer_seq,
                        'start': end - length,
                        'end': end,
                        'length': length,
                        'gc_content': gc_content,
                        'tm': tm
                    })
                    
        # Sort primers by quality (closest to optimal Tm and GC content)
        optimal_tm = (self.min_tm + self.max_tm) / 2
        optimal_gc = (self.min_gc + self.max_gc) / 2
        
        for primer in primers:
            tm_score = 1 - abs(primer['tm'] - optimal_tm) / optimal_tm
            gc_score = 1 - abs(primer['gc_content'] - optimal_gc) / optimal_gc
            primer['quality_score'] = (tm_score + gc_score) / 2
            
        primers.sort(key=lambda x: x['quality_score'], reverse=True)
        
        return primers[:10]  # Return top 10 primers
        
    def find_primer_pairs(self, primers):
        """Find compatible primer pairs"""
        forward_primers = [p for p in primers if p['type'] == 'forward']
        reverse_primers = [p for p in primers if p['type'] == 'reverse']
        
        pairs = []
        
        for f_primer in forward_primers[:5]:  # Top 5 forward
            for r_primer in reverse_primers[:5]:  # Top 5 reverse
                if (r_primer['start'] > f_primer['end'] + 100 and  # Minimum product size
                    r_primer['start'] < f_primer['end'] + 2000 and  # Maximum product size
                    not self.check_dimer(f_primer['sequence'], r_primer['sequence'])):
                    
                    product_size = r_primer['start'] - f_primer['end']
                    tm_diff = abs(f_primer['tm'] - r_primer['tm'])
                    
                    pairs.append({
                        'forward': f_primer,
                        'reverse': r_primer,
                        'product_size': product_size,
                        'tm_difference': tm_diff,
                        'quality_score': (f_primer['quality_score'] + r_primer['quality_score']) / 2 - tm_diff / 10
                    })
                    
        pairs.sort(key=lambda x: x['quality_score'], reverse=True)
        return pairs[:5]  # Return top 5 pairs
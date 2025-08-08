from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import re

class ORFAnalyzer:
    def __init__(self):
        self.start_codons = ['ATG']
        self.stop_codons = ['TAA', 'TAG', 'TGA']
        self.min_orf_length = 100  # Minimum ORF length in nucleotides
        
    def find_orfs(self, sequence, min_length=None):
        """Find all Open Reading Frames in the sequence"""
        if min_length is None:
            min_length = self.min_orf_length
            
        sequence = sequence.upper()
        orfs = []
        
        # Check all 6 reading frames (3 forward, 3 reverse)
        for strand in [1, -1]:
            if strand == -1:
                seq = str(Seq(sequence).reverse_complement())
            else:
                seq = sequence
                
            for frame in range(3):
                frame_seq = seq[frame:]
                
                # Find ORFs in this frame
                start_pos = 0
                while start_pos < len(frame_seq) - 2:
                    # Look for start codon
                    start_found = False
                    for start_codon in self.start_codons:
                        pos = frame_seq.find(start_codon, start_pos)
                        if pos != -1 and pos % 3 == 0:  # Must be in frame
                            start_pos = pos
                            start_found = True
                            break
                    
                    if not start_found:
                        break
                        
                    # Look for stop codon
                    stop_pos = None
                    for i in range(start_pos + 3, len(frame_seq) - 2, 3):
                        codon = frame_seq[i:i+3]
                        if codon in self.stop_codons:
                            stop_pos = i + 3
                            break
                    
                    if stop_pos is None:
                        stop_pos = len(frame_seq)
                        
                    orf_length = stop_pos - start_pos
                    if orf_length >= min_length:
                        orf_seq = frame_seq[start_pos:stop_pos]
                        
                        # Calculate actual position in original sequence
                        if strand == 1:
                            actual_start = frame + start_pos
                            actual_end = frame + stop_pos
                        else:
                            actual_start = len(sequence) - (frame + stop_pos)
                            actual_end = len(sequence) - (frame + start_pos)
                            
                        orfs.append({
                            'sequence': orf_seq,
                            'start': actual_start,
                            'end': actual_end,
                            'length': orf_length,
                            'strand': strand,
                            'frame': frame + 1,
                            'protein': str(Seq(orf_seq).translate())
                        })
                    
                    start_pos += 3
                    
        # Sort ORFs by length (longest first)
        orfs.sort(key=lambda x: x['length'], reverse=True)
        return orfs
        
    def analyze_protein(self, protein_sequence):
        """Analyze protein properties"""
        if not protein_sequence or '*' in protein_sequence[:-1]:  # Skip if contains stop codons
            return None
            
        try:
            analysis = ProteinAnalysis(protein_sequence.replace('*', ''))  # Remove stop codon
            
            return {
                'length': len(protein_sequence) - 1,  # Exclude stop codon
                'molecular_weight': analysis.molecular_weight(),
                'isoelectric_point': analysis.isoelectric_point(),
                'instability_index': analysis.instability_index(),
                'gravy': analysis.gravy(),  # Grand average of hydropathy
                'amino_acid_percent': analysis.get_amino_acids_percent(),
                'secondary_structure': analysis.secondary_structure_fraction()
            }
        except Exception as e:
            return {'error': str(e)}
            
    def get_orf_summary(self, sequence):
        """Get summary of ORF analysis"""
        orfs = self.find_orfs(sequence)
        
        summary = {
            'total_orfs': len(orfs),
            'longest_orf': orfs[0] if orfs else None,
            'orfs_by_strand': {
                'forward': len([orf for orf in orfs if orf['strand'] == 1]),
                'reverse': len([orf for orf in orfs if orf['strand'] == -1])
            },
            'orfs_by_frame': {}
        }
        
        for frame in range(1, 7):
            if frame <= 3:
                strand = 1
                actual_frame = frame
            else:
                strand = -1
                actual_frame = frame - 3
            summary['orfs_by_frame'][f'Frame {frame}'] = len([
                orf for orf in orfs 
                if orf['strand'] == strand and orf['frame'] == actual_frame
            ])
            
        return summary, orfs
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np
from Bio.SeqUtils import GC
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import seaborn as sns

class SequenceVisualizer:
    def __init__(self):
        self.figure = Figure(figsize=(12, 8))
        self.canvas = FigureCanvas(self.figure)
        
    def plot_quality_scores(self, record):
        """Plot quality scores from AB1 file"""
        self.figure.clear()
        
        if hasattr(record, 'letter_annotations') and 'phred_quality' in record.letter_annotations:
            quality_scores = record.letter_annotations['phred_quality']
            positions = range(1, len(quality_scores) + 1)
            
            ax = self.figure.add_subplot(111)
            ax.plot(positions, quality_scores, 'b-', linewidth=1)
            ax.fill_between(positions, quality_scores, alpha=0.3)
            ax.axhline(y=20, color='r', linestyle='--', label='Quality threshold (Q20)')
            ax.set_xlabel('Position')
            ax.set_ylabel('Quality Score')
            ax.set_title('Sequence Quality Scores')
            ax.legend()
            ax.grid(True, alpha=0.3)
        else:
            ax = self.figure.add_subplot(111)
            ax.text(0.5, 0.5, 'No quality data available', 
                   horizontalalignment='center', verticalalignment='center',
                   transform=ax.transAxes, fontsize=14)
            ax.set_title('Quality Scores')
            
        self.figure.tight_layout()
        self.canvas.draw()
        
    def plot_gc_content(self, sequence, window_size=100):
        """Plot GC content along the sequence"""
        self.figure.clear()
        
        if len(sequence) < window_size:
            window_size = len(sequence) // 4 if len(sequence) > 4 else 1
            
        gc_contents = []
        positions = []
        
        for i in range(0, len(sequence) - window_size + 1, window_size // 2):
            window = sequence[i:i + window_size]
            gc_content = GC(window)
            gc_contents.append(gc_content)
            positions.append(i + window_size // 2)
            
        ax = self.figure.add_subplot(111)
        ax.plot(positions, gc_contents, 'g-', linewidth=2, marker='o', markersize=4)
        ax.axhline(y=50, color='r', linestyle='--', alpha=0.7, label='50% GC')
        ax.set_xlabel('Position')
        ax.set_ylabel('GC Content (%)')
        ax.set_title(f'GC Content (Window size: {window_size})')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_ylim(0, 100)
        
        self.figure.tight_layout()
        self.canvas.draw()
        
    def plot_alignment_overview(self, alignment_data):
        """Plot alignment overview with matches and mismatches"""
        self.figure.clear()
        
        ax = self.figure.add_subplot(111)
        
        # Extract alignment information
        seq1, seq2 = alignment_data['seq1'], alignment_data['seq2']
        matches = [1 if a == b and a != '-' else 0 for a, b in zip(seq1, seq2)]
        gaps = [1 if a == '-' or b == '-' else 0 for a, b in zip(seq1, seq2)]
        
        positions = range(len(matches))
        
        # Create stacked bar chart
        ax.bar(positions, matches, color='green', alpha=0.7, label='Matches')
        ax.bar(positions, gaps, bottom=matches, color='red', alpha=0.7, label='Gaps')
        
        ax.set_xlabel('Position')
        ax.set_ylabel('Match/Gap')
        ax.set_title('Alignment Overview')
        ax.legend()
        
        self.figure.tight_layout()
        self.canvas.draw()
        
    def plot_amino_acid_composition(self, protein_sequence):
        """Plot amino acid composition"""
        self.figure.clear()
        
        if not protein_sequence:
            return
            
        analysis = ProteinAnalysis(protein_sequence)
        aa_percent = analysis.get_amino_acids_percent()
        
        amino_acids = list(aa_percent.keys())
        percentages = list(aa_percent.values())
        
        ax = self.figure.add_subplot(111)
        bars = ax.bar(amino_acids, [p * 100 for p in percentages], 
                     color=plt.cm.Set3(np.linspace(0, 1, len(amino_acids))))
        
        ax.set_xlabel('Amino Acid')
        ax.set_ylabel('Percentage (%)')
        ax.set_title('Amino Acid Composition')
        ax.tick_params(axis='x', rotation=45)
        
        # Add percentage labels on bars
        for bar, pct in zip(bars, percentages):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                   f'{pct*100:.1f}%', ha='center', va='bottom', fontsize=8)
        
        self.figure.tight_layout()
        self.canvas.draw()
        
    def get_canvas(self):
        """Return the matplotlib canvas"""
        return self.canvas
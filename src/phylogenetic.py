from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np

class PhylogeneticAnalyzer:
    def __init__(self):
        self.figure = Figure(figsize=(12, 8))
        self.canvas = FigureCanvas(self.figure)
        self.tree = None
        
    def create_multiple_alignment(self, sequences_dict):
        """Create multiple sequence alignment from dictionary of sequences"""
        records = []
        for name, seq in sequences_dict.items():
            record = SeqRecord(Seq(seq), id=name, description="")
            records.append(record)
            
        # Simple multiple alignment (in practice, you'd use MUSCLE, ClustalW, etc.)
        alignment = MultipleSeqAlignment(records)
        return alignment
        
    def calculate_distance_matrix(self, alignment):
        """Calculate distance matrix from alignment"""
        calculator = DistanceCalculator('identity')
        distance_matrix = calculator.get_distance(alignment)
        return distance_matrix
        
    def construct_tree(self, distance_matrix):
        """Construct phylogenetic tree using UPGMA method"""
        constructor = DistanceTreeConstructor()
        tree = constructor.upgma(distance_matrix)
        self.tree = tree
        return tree
        
    def plot_tree(self, tree=None):
        """Plot phylogenetic tree"""
        if tree is None:
            tree = self.tree
            
        if tree is None:
            return
            
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        
        Phylo.draw(tree, axes=ax, do_show=False)
        ax.set_title('Phylogenetic Tree (UPGMA)')
        
        self.figure.tight_layout()
        self.canvas.draw()
        
    def get_canvas(self):
        """Return the matplotlib canvas"""
        return self.canvas
        
    def analyze_sequences(self, sequences_dict):
        """Complete phylogenetic analysis pipeline"""
        try:
            alignment = self.create_multiple_alignment(sequences_dict)
            distance_matrix = self.calculate_distance_matrix(alignment)
            tree = self.construct_tree(distance_matrix)
            self.plot_tree(tree)
            return tree, distance_matrix
        except Exception as e:
            raise Exception(f"Phylogenetic analysis failed: {str(e)}")